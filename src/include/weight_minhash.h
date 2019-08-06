//
// Created by junior on 2019/8/6.
//

#ifndef LSH_CPP_WEIGHT_MINHASH_H
#define LSH_CPP_WEIGHT_MINHASH_H

#include "lsh_cpp.h"
#include "util.h"

namespace LSH_CPP {
//    using SpraseMatrixXf = Eigen::SparseMatrix<float, Eigen::RowMajor>;
//    using Triplet_f = Eigen::Triplet<float>;

    template<size_t dim, size_t sample_size = 128, size_t seed = 1, typename RandomGenerator = std::mt19937_64>
    struct RandomSample {
        using SampleMatrixType = Eigen::ArrayXXf;// 注意Eigen的二维Array类型是ArrayXXf
        SampleMatrixType r_k, ln_c_k, beta_k;

        // 注意二维Array的维度用 (dim,sample_size), 这样通过 r_k.col(i) 得到的就是一个 ArrayXf (列向量形式),
        // 后面WeightMinHash::update()里就可以直接和其他列向量形式的ArrayXf做element-wise计算了.
        // 如果用 (sample_size,dim),那么r_k.row(i)得到的就是行向量形式的 ArrayXf.transpose(),
        // 直接和其他的ArrayXf计算会变成内积操作,导致结果不对.
        explicit RandomSample() : r_k(SampleMatrixType(dim, sample_size)),
                                  ln_c_k(SampleMatrixType(dim, sample_size)),
                                  beta_k(SampleMatrixType(dim, sample_size)) {
            RandomGenerator generator(seed);
            std::gamma_distribution<float> gamma_dis(2, 1); // Gamma(2,1)
            std::uniform_real_distribution<float> uniform_dis(0, 1); // uniform(0,1)
            for (size_t n_sample = 0; n_sample < sample_size; n_sample++) {
                for (size_t i = 0; i < dim; i++) {
                    r_k(i, n_sample) = gamma_dis(generator);               // r_k ~ Gamma(2,1)
                    ln_c_k(i, n_sample) = std::log(gamma_dis(generator));  // ln_c_k ~ ln(Gamma(2,1))
                    beta_k(i, n_sample) = uniform_dis(generator);
                }
            }
        }
    };

    // WeightMinHash 模板类
    template<
            size_t dim,                               // 全集大小(即权重向量维度)
            typename UpdateInterfaceElementType,      // update接口的参数类型
            size_t sample_size = 128,                 // 采样次数
            size_t seed = 1,                          // 随机数种子
            typename RandomGenerator = std::mt19937_64, // 随机数生成器
            typename SelectByDim = void,              // 根据dim来筛选不同的WeightMinHash实现
            typename SelectByUpdateElementType = void // 对于不同的WeighMinHash实现,需要限制不同的update接口的参数类型
    >
    struct WeightMinHash {
    };

    constexpr size_t dim_gap_for_different_impl = 100000; // dim <= 100000 用dense Array计算; dim > 100000 用稀疏矩阵计算

    // WeightMinHash for dense matrix/vector,仅用于小规模的集合的测试
    // update接口接受权重向量作为参数,权重向量类型必须是数值类型(整数或浮点)
    template<size_t dim, typename UpdateInterfaceElementType, size_t sample_size, size_t seed, typename RandomGenerator>
    struct WeightMinHash<dim, UpdateInterfaceElementType, sample_size, seed, RandomGenerator,
            typename std::enable_if_t<(dim <= dim_gap_for_different_impl)>,
            typename std::enable_if_t<std::is_arithmetic<UpdateInterfaceElementType>::value>> {
    private:
        using WeightType = UpdateInterfaceElementType;
        using MinHashValueFirstType = size_t;        // Type of k*   (k* 是下标,所以用size_t类型)
        using MinHashValueSecondType = int_fast32_t; // Type of t_k* (t_k* 计算过程中最后用了取整操作,所以t_k*是整型)
        using MinHashValueType = std::pair<MinHashValueFirstType, MinHashValueSecondType>; // 记录 k* , t_k*

        static RandomSample<dim, sample_size, seed, RandomGenerator> sample;

    public:
        std::vector<MinHashValueType> hash_values;

    public:
        explicit WeightMinHash() {
            hash_values.resize(sample_size);
        }

        // update by weight vector
        // 注意Eigen计算时的 Aliasing.
        // 见: https://eigen.tuxfamily.org/dox/group__TopicAliasing.html
        // 另外要注意 ArrayXf 是列向量形式的,如果要做element-wise的计算就必须确保其他变量也是列向量形式,否则就会变成内积操作导致结果错误.
        bool update(std::vector<WeightType> &weight_vector) {
            assert(weight_vector.size() == dim);
            // 外部权重类型可以是任意的数值类型比如 int/float/size_t, 但在预处理后统统转换为 ArrayXf w_v, 也就是float类型权重.
            // 所以理论上权重的大小最多不能超过 float::max() 的范围.
            Eigen::Map<Eigen::Array<WeightType, Eigen::Dynamic, 1>> temp(weight_vector.data(), weight_vector.size());
            Eigen::ArrayXf w_v;
            if ((temp == 0).all()) return false; // 如果数据全0,则update失败,返回false
            if ((temp == 0).any()) { // 至少存在一个位置的权重为0,将其改为float的最小值.
                // change zero value in weight_vector to numeric_limits<float>::min() = 1.17549e-38
                // instead of nan (because Eigen process nan will be undefined).
                // select function : C = (expression).select(A,B); => C = expression ?  A : B (element-wise)
                w_v = (temp == 0).select(Eigen::ArrayXf::Constant(temp.size(), std::numeric_limits<float>::min()),
                                         temp);
            } else {
                // 如果一个0都没有,直接把原数据cast为float的ArrayXf即可,即调用 Array.cast<other_type>().
                // 但是因为这里temp的类型本身也是需要经过推导的,所以temp调用cast方法的时候需要加上template前缀.
                w_v = temp.template cast<float>();
            }
            for (size_t i = 0; i < sample_size; i++) {
                Eigen::ArrayXf t_k = (w_v.log() / sample.r_k.col(i) + sample.beta_k.col(i)).floor();
                Eigen::ArrayXf ln_y = (t_k - sample.beta_k.col(i)) * sample.r_k.col(i);
                Eigen::ArrayXf ln_a = sample.ln_c_k.col(i) - ln_y - sample.r_k.col(i);
                Eigen::ArrayXf::Index k_star, garbage; // VectorXf = MatrixXf<float,Dynamic,1>,第二个下标默认是0,随便用某个变量获取即可
                ln_a.minCoeff(&k_star, &garbage);  // 找到ln_a最小值的位置下标 (k*,0)
                hash_values[i] = {static_cast<MinHashValueFirstType >(k_star),                 // k*
                                  static_cast<MinHashValueSecondType >(t_k(k_star, 0))};   // t_k*
            }
            return true;
        }
    };

    template<size_t dim, typename UpdateInterfaceElementType, size_t sample_size, size_t seed, typename RandomGenerator>
    RandomSample<dim, sample_size, seed, RandomGenerator>
            WeightMinHash<dim, UpdateInterfaceElementType, sample_size, seed, RandomGenerator,
                    typename std::enable_if_t<(dim <= dim_gap_for_different_impl)>,
                    typename std::enable_if_t<std::is_arithmetic<UpdateInterfaceElementType>::value>>::sample{};


    // WeightMinHash for sparse matrix/vector,用于大规模的实际文本数据
    // update 接收单一集合作为参数,内部生成权重向量:
    // 其中权重向量的元素位置编码由内部全局哈希表决定(但需要集合元素类型提供接口value()获取集合元素的实际值,集合元素的实际值也可以是任意类型);
    // 权重值需要集合元素类型提供接口获取,权重值类型可以任意(比如整数的重复次数权重,tf-idf浮点权重等).
    template<size_t dim, typename UpdateInterfaceElementType, size_t sample_size, size_t seed, typename RandomGenerator>
    struct WeightMinHash<dim, UpdateInterfaceElementType, sample_size, seed, RandomGenerator,
            typename std::enable_if_t<(dim > dim_gap_for_different_impl)>,
            typename std::enable_if_t<std::is_arithmetic<typename UpdateInterfaceElementType::WeightType>::value>> {
    private:
        using SetElementType = UpdateInterfaceElementType;            // update解释集合作为参数,然后内部生成权重向量(元素位置编码也是内部生成)
        using SetValueType  = typename SetElementType::ValueType;     // 需要SetElementType提供,同时需要提供 value() 接口
        using SetWeightType = typename SetElementType::WeightType;    // 需要SetElementType提供,同时需要提供 weight() 接口
        static HashMap <SetValueType, size_t> global_weight_vector_index_map; //　权重向量里每一个元素的位置编码.
        static size_t global_index;                                          // 初始为0
        static RandomGenerator generator;                                    // 初始为generator(seed)
        static std::uniform_real_distribution<float> uniform_dis;            // uniform dis ~ (0,1)
        static std::gamma_distribution<float> gamma_dis;                     // gamma dis ~ (2,1)

    public:
        explicit WeightMinHash() = default;

        // update by single set
        void update(const HashSet <SetElementType> &) {}
    };

    template<size_t dim, typename T, size_t sample_size, size_t seed, typename RG>
    double weight_minhash_jaccard(const WeightMinHash<dim, T, sample_size, seed, RG, void, void> &A,
                                  const WeightMinHash<dim, T, sample_size, seed, RG, void, void> &B) {
        double count = 0;
        for (size_t i = 0; i < sample_size; i++) {
            if (A.hash_values[i].first == B.hash_values[i].first &&
                A.hash_values[i].second == B.hash_values[i].second) {
                count++;
            }
        }
        return count / (double) (sample_size);
    }

    template<typename WeightType, typename= std::enable_if_t<std::is_arithmetic<WeightType>::value>>
    double generalized_jaccard_similarity(const std::vector<WeightType> &A_weight_set,
                                          const std::vector<WeightType> &B_weight_set) {
        assert(A_weight_set.size() == B_weight_set.size());// dimension should be equal
        double min_ = 0, max_ = 0;
        for (size_t i = 0; i < A_weight_set.size(); i++) {
            min_ += std::min(A_weight_set[i], B_weight_set[i]);
            max_ += std::max(A_weight_set[i], B_weight_set[i]);
        }
        return min_ / max_;
    }
}
#endif //LSH_CPP_WEIGHT_MINHASH_H
