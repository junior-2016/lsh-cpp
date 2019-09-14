//
// Created by junior on 2019/8/6.
//

#ifndef LSH_CPP_WEIGHT_MINHASH_H
#define LSH_CPP_WEIGHT_MINHASH_H

#include "lsh_cpp.h"
#include "util.h"

namespace LSH_CPP {
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
                Eigen::ArrayXf::Index k_star, garbage; // ArrayXf 第二个下标默认是0,随便用某个变量获取即可
                ln_a.minCoeff(&k_star, &garbage);      // 找到ln_a最小值的位置下标 (k*,0)
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


    // 稀疏采样矩阵. 当权重向量位置编号从0逐渐增加时,逐行扩展采样矩阵的分布参数,而不是立刻构建完整的采样矩阵,否则内存开销过大,将无法运行程序.
    template<size_t dim, size_t sample_size, size_t seed, typename RandomGenerator>
    struct SparseSampleMatrix {
        RandomGenerator generator;                           // 初始为generator(seed)
        std::uniform_real_distribution<float> uniform_dis;   // uniform dis ~ (0,1)
        std::gamma_distribution<float> gamma_dis;            // gamma dis ~ (2,1)

        // row: dynamic(最大可能达到dim,不过大多数场合达不到); col: 固定 sample_size 维度.
        // 设计成这种维度的目的是: 非常容易将其转换为 Eigen::ArrayXf 和 Eigen::ArrayXXf, 只要通过 Eigen::Map 转换即可.
        // 比如: 当我们知道当前的集合需要计算 0,2,7,12 这四个位置元素的weight_min_hash, 要将采样矩阵的参数拿出来聚合成 Dense ArrayXXf 以便于后面mkl加速.
        // 只需要拿出 Matrix[0],Matrix[2],Matrix[7],Matrix[12] 这四个位置的std::array的data指针(调用data()方法),
        // 再通过Eigen::Map直接映射为Eigen::ArrayXXf 即可,这个过程只需要耗费遍历集合元素的时间,而且不用多余的拷贝(直接拿原数据指针,快且节省空间).
        using RowType = std::array<float, sample_size>;             // 矩阵每一行的类型
        using Matrix = std::vector<RowType>;                        // 矩阵类型

        Matrix r_k, ln_c_k, beta_k;                                 // 采样矩阵
        static constexpr size_t init_reserve_number = 10000;        // 采样矩阵初始预留的10000行空间,避免push_back频繁alloc

        explicit SparseSampleMatrix() : generator(seed), uniform_dis(0, 1), gamma_dis(2, 1) {
            r_k.reserve(init_reserve_number);
            ln_c_k.reserve(init_reserve_number);
            beta_k.reserve(init_reserve_number);
        }

        /**
         * 单行扩展sample_matrix
         */
        void update_matrix() {
            RowType r_k_temp, ln_c_k_temp, beta_k_temp;
            for (size_t j = 0; j < sample_size; j++) {
                r_k_temp[j] = gamma_dis(generator);                  // r_k ~ Gamma(2,1)
                ln_c_k_temp[j] = std::log(gamma_dis(generator));     // ln_c_k ~ ln(Gamma(2,1))
                beta_k_temp[j] = uniform_dis(generator);             // beta_k ~ uniform(0,1)
            }
            r_k.push_back(r_k_temp);
            ln_c_k.push_back(ln_c_k_temp);
            beta_k.push_back(beta_k_temp);
        }

        /**
         * @param expand_sample_matrix_rows 需要扩展的 sample_matrix 行数
         */
        void update_matrix(size_t sample_matrix_expand_rows) {
            while (sample_matrix_expand_rows--) {
                update_matrix();
            }
        }
    };

    // TODO: WeightMinHash for sparse matrix/vector 的测试结果低于 dense sample matrix 的实现,
    //       可能是稀疏采样矩阵连续取值的问题,即将所有出现元素的编码集中在了一起,
    //       而dense的实现是先对全集元素编码后,再从采样矩阵中取值的,因此集合中出现的元素在全集取到的编码是非集中的,
    //       对应的采样矩阵的值当然也不是集中的,所以可能随机效果更好.
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

        using MinHashValueFirstType = size_t;        // Type of k*   (k* 是下标,所以用size_t类型)
        using MinHashValueSecondType = int_fast32_t; // Type of t_k* (t_k* 计算过程中最后用了取整操作,所以t_k*是整型)
        using MinHashValueType = std::pair<MinHashValueFirstType, MinHashValueSecondType>; // 记录 k* , t_k*

        static HashMap <SetValueType, size_t> global_weight_vector_pos_map; //　记录权重向量中每一个元素的位置编码.
        static SparseSampleMatrix<dim, sample_size, seed, RandomGenerator> sample_matrix; // 静态全局采样矩阵
        static size_t global_pos;                                                         // 全局位置编码,初始为0

    public:
        std::vector<MinHashValueType> hash_values;

    public:
        explicit WeightMinHash() {
            hash_values.resize(sample_size);
        }

        // update by single set
        // 注意 SetElementType 必须有 weight() 和 value() 两个 public 接口,否则编译会失败.
        void update(const HashSet <SetElementType> &set) {
            Eigen::ArrayXf weight_vector(set.size()); // 权重向量中非0部分聚合成Dense Array(注意 ArrayXf 是列向量形式的).
            Eigen::ArrayXXf r_k(set.size(), sample_size);
            Eigen::ArrayXXf ln_c_k(set.size(), sample_size);
            Eigen::ArrayXXf beta_k(set.size(), sample_size);
            int row_index = 0;
            std::pair<typename HashMap<SetValueType, size_t>::const_iterator, bool> ret;
            for (const auto &element:set) {
                // insert and return { const_iterator, bool }.
                ret = global_weight_vector_pos_map.try_emplace(element.value(), global_pos);

                // insert success and update global_pos && update one row of sample_matrix for this new element.
                if (ret.second) {
                    global_pos++;
                    sample_matrix.update_matrix();
                }

                // 无论插入成不成功都可以通过 (*(ret.first)).second 直接获得当前元素element已存在的位置编码pos.
                // 通过这个位置编码pos获取sample_matrix对应位置的参数,并聚合成 Dense ArrayXXf
                // 注意: 这里不需要保证 ArrayXXf 的 r_k,ln_c_k,beta_k 从上到下是按照元素位置编码的顺序摆放的,
                // 它只要保证与 weight_vector(列向量) 每一行对应的元素一致即可.
                auto pos = (*(ret.first)).second;
                r_k.row(row_index) = Eigen::Map<Eigen::ArrayXf>(sample_matrix.r_k[pos].data(), sample_size);
                ln_c_k.row(row_index) = Eigen::Map<Eigen::ArrayXf>(sample_matrix.ln_c_k[pos].data(), sample_size);
                beta_k.row(row_index) = Eigen::Map<Eigen::ArrayXf>(sample_matrix.beta_k[pos].data(), sample_size);

                // 更新 weight_vector (只保留权重非0,即当前集合元素的部分)
                weight_vector(row_index) = element.weight();

                // 更新行号
                row_index++;
            }
            // 剩下的部分就和 dim<=100000 的实现一样了,直接用 Eigen::ArrayXf 之间的 Element-wise 计算即可,从而利用eigen的dense数组优化.
            for (long i = 0; i < static_cast<long>(sample_size); i++) {
                Eigen::ArrayXf t_k = (weight_vector.log() / r_k.col(i) + beta_k.col(i)).floor();
                Eigen::ArrayXf ln_y = (t_k - beta_k.col(i)) * r_k.col(i);
                Eigen::ArrayXf ln_a = ln_c_k.col(i) - ln_y - r_k.col(i);
                Eigen::ArrayXf::Index k_star, garbage; // ArrayXf 第二个下标默认是0,随便用某个变量获取即可
                ln_a.minCoeff(&k_star, &garbage);      // 找到ln_a最小值的位置下标 (k*,0)
                hash_values[i] = {static_cast<MinHashValueFirstType >(k_star),                 // k*
                                  static_cast<MinHashValueSecondType >(t_k(k_star, 0))};   // t_k*
            }
        }
    };

    template<size_t dim, typename T, size_t sample_size, size_t seed, typename RandomGenerator>
    HashMap<typename T::ValueType, size_t> WeightMinHash<dim, T, sample_size, seed, RandomGenerator,
            typename std::enable_if_t<(dim > dim_gap_for_different_impl)>,
            typename std::enable_if_t<std::is_arithmetic<typename T::WeightType>::value>>::global_weight_vector_pos_map{};

    template<size_t dim, typename T, size_t sample_size, size_t seed, typename RandomGenerator>
    SparseSampleMatrix<dim, sample_size, seed, RandomGenerator>
            WeightMinHash<dim, T, sample_size, seed, RandomGenerator,
                    typename std::enable_if_t<(dim > dim_gap_for_different_impl)>,
                    typename std::enable_if_t<std::is_arithmetic<typename T::WeightType>::value>>::sample_matrix{};

    template<size_t dim, typename T, size_t sample_size, size_t seed, typename RandomGenerator>
    size_t WeightMinHash<dim, T, sample_size, seed, RandomGenerator,
            typename std::enable_if_t<(dim > dim_gap_for_different_impl)>,
            typename std::enable_if_t<std::is_arithmetic<typename T::WeightType>::value>>::global_pos = 0;

    // 计算 weight_minhash 的 jaccard_similarity
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

    // 计算带权重集实际的 jaccard_similarity, 需要预先提供权重向量
    // 参考: https://en.wikipedia.org/wiki/Jaccard_index#Generalized_Jaccard_similarity_and_distance
    // 例子: A = { a, a, a, b, b, c } B = { a, a, b, b, b, d }
    // Generalized_Jaccard_Similarity (A,B)
    //   = (min(A_a,B_a)+min(A_b,B_b)+min(A_c,B_c)+min(A_d,B_d)) / (max(A_a,B_a)+max(A_b,B_b)+max(A_c,B_c)+max(A_d,B_d))
    //   = ( 2 + 2 + 0 + 0 ) / ( 3 + 3 + 1 + 1 ) = 4 / 8  = 0.5
    // 当然你也可以这样通过权重向量来计算,预先给 a, b, c, d 一个位置编码 0,1,2,3 代表四个维度,
    // 然后 A_weight_vector: V_A = [ 3 2 1 0 ]; B_weight_vector: V_B = [ 2 3 0 1 ]
    // 然后计算 Generalized_jaccard_similarity = \sigma_{i}^{n=4} min(V_A[i],V_B[i])/ \sigma_{i}^{n=4} max(V_A[i],V_B[i])
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

    // 计算带权重集实际的 jaccard_similarity, 不需要提供权重向量, 直接提供两个集合即可, 计算方法如下:
    // 如果 A 与 B 均不存在的元素位置,必然权重都是0,那么 min(0,0) = max(0,0) = 0, 所以这一部分直接省略;
    // 对于 A存在,B不存在 或者 B存在,A不存在的部分,min必然是0,忽略,而max就是对应存在部分的权重值;
    // 最后交集部门按照基本公式计算即可.
    // 注意 SetElementType 必须有 weight() 接口.
    template<typename SetElementType,
            typename = std::enable_if_t<std::is_arithmetic<typename SetElementType::WeightType>::value>>
    double generalized_jaccard_similarity(const HashSet <SetElementType> &A, const HashSet <SetElementType> &B) {
        double min_ = 0, max_ = 0;
        for (const auto &a:A) {
            if (auto ret = B.find(a);ret != B.end()) { // 交集部分
                min_ += std::min(a.weight(), (*ret).weight());
                max_ += std::max(a.weight(), (*ret).weight());
            } else { // A 独有部分
                max_ += a.weight();
            }
        }
        for (const auto &b:B) {
            if (A.find(b) == A.end()) { // B独有部分
                max_ += b.weight();
            }
        }
        return min_ / max_;
    }
}
#endif //LSH_CPP_WEIGHT_MINHASH_H
