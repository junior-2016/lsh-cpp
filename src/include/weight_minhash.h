//
// Created by junior on 2019/7/30.
//

#ifndef LSH_CPP_WEIGHT_MINHASH_H
#define LSH_CPP_WEIGHT_MINHASH_H

#include "lsh_cpp.h"
#include "util.h"

namespace LSH_CPP {
    // TODO Weight MinHash 算法. 这个算法必须提前知道集合的全集大小,对于基因相似度查询这个场景来说,
    //     Weight MinHash是可行的,比如取 k=9 的 k_mer,单个位置只有A/T/C/G四种可能,所以全集大小是 4^9,是确定的值.
    //     如果不用 WeightMinHash ,而是 MinHash ,那就必须先用 hash_set 过滤并且统计重复个数,这里就存在一个大量的
    //     hash key -> insert -> find 的操作开销; 而且 hash_set 对 集合元素进行哈希后储存, 接着 MinHash.update()
    //     还要对集合元素哈希为整数再做MinHash计算,这就两遍哈希了,开销极大.
    //     使用 Weight MinHash,可以不用 hash_set 过滤储存集合元素,直接用原始的 vector/array 储存即可
    //     (有重复元素也无所谓),然后借助一个公式来计算权重.
    //     paper: [1] Improved Consistent Sampling, Weighted Minhash and L1 Sketching
    //            [2] A Review for Weighted MinHash Algorithms (2018)

    // TODO 当前使用的 Weight MinHash 算法缺点在于无法处理流式数据,必须提前获取全集大小,
    //  并且为了维护全集维度向量/分布采样矩阵耗费内存大(这个问题基本有解决方案,现在内存和运行时间均取决于weight>0的部分,所以已经不算是问题),
    //  所以一个更加需要解决(扩展)的问题是: 提供一个可以处理 stream data 的 weight-minhash 算法.

    template<
            size_t dim,
            size_t sample_size = 128,
            size_t seed = 1,
            typename RandomGenerator = std::mt19937_64,
            typename random_number_type = float>
    struct RandomSample {
        // TODO: 使用 (sample_size,dim) 规模的三个静态浮点矩阵储存, 空间开销过大,
        //        考虑 (128,4^12) 大小的规模, 需要 128*4^12*3*4 / (1024*1024*1024) GB = 20多G,
        //        而且考虑到 gamma(2,1)(范围基本在[0,7]) uniform(0,1)(范围[0,1]) 生成的浮点数都超级小,
        //        如果遵循IEEE754标准存4个字节就太浪费空间了,需要考虑对浮点数做无损压缩.
        //        A方案: 用 zpf Compressed C++ Arrays 储存 r_k,ln_c_k,beta_k;
        //              用 Eigen sparse matrix(vector) 做计算.
        //              weight_min_hash::update(std::vector<pair<index,weight>> & weights){
        //                      Eigen::sparse_vector weight_vector <- weights;
        //                      Eigen::sparse_vector r_k_part <- r_k.get_part_by(weights.index); // ln_ck 和 beta_k同理
        //              }
        //       这样就兼顾了储存(zpf压缩储存完全没有问题,因为r_k,ln_ck,beta_k都是由固定分布生成,具有规律,且数值很小,非常符合zpf压缩算法的要求)
        //               和 运算效率(运算用Eigen稀疏矩阵/向量),
        //       update时只取参数中对应weight>0的部分,其他weight=0的部分的参数不用取出,所以对内存开销没有多大影响.
        //       B方案: r_k,ln_c_k,beta_k 依然用 Eigen::Sparse_matrix 储存,并提前声明维度大小(n_sample,dim),但是此时不存在任何数据,
        //              在 update(std::vector<pair<index,weight>> & weights) 时, 检查 r_k 等矩阵在index位置是否为0,如果是,
        //              就生成对应的参数,如果不是就跳过生成过程. 这样对不同的 weights 向量,使用的 r_k 等分布采样参数矩阵依然是相同的,
        //              只不过是采用了动态生成方法,一点点构造出来而已,但相比原来的算法其实时间和内存开销更少,而且无需考虑浮点数压缩的问题,
        //              除非你运行过程中提供的 weights 向量中的元素能完全覆盖整一个全集,但这种情况出现的概率极小,
        //              因为现实数据基本是大量重复且只能构成全集的某个小的子集而已.
        //              更重要的是, 当 r_k 等分布参数矩阵和 weighs_vector 都是 Eigen::Sparse 类型时,可以直接写矩阵向量乘法,
        //              vectorization 优化更彻底.
        //       从效率角度上选择B方案更好一些.因为A方案依旧需要一次性储存全部的参数(虽然压缩了但空间开销依旧很大)

        std::vector<std::vector<random_number_type>> r_k;
        std::vector<std::vector<random_number_type>> ln_c_k;
        std::vector<std::vector<random_number_type>> beta_k;

        explicit RandomSample() {
            RandomGenerator generator(seed);
            std::gamma_distribution<random_number_type> gamma_dis(2, 1); // Gamma(2,1)
            std::uniform_real_distribution<random_number_type> uniform_dis(0, 1); // uniform(0,1)
            r_k.reserve(sample_size);
            ln_c_k.reserve(sample_size);
            beta_k.reserve(sample_size);
            for (size_t i = 0; i < sample_size; i++) {
                std::vector<random_number_type> r_k_temp, ln_c_k_temp, beta_k_temp;
                r_k_temp.reserve(dim);
                ln_c_k_temp.reserve(dim);
                beta_k_temp.reserve(dim);
                for (size_t j = 0; j < dim; j++) {
                    r_k_temp.push_back(gamma_dis(generator));               // r_k ~ Gamma(2,1)
                    ln_c_k_temp.push_back(std::log(gamma_dis(generator)));  // ln_c_k ~ ln(Gamma(2,1))
                    beta_k_temp.push_back(uniform_dis(generator));
                }
                r_k.push_back(r_k_temp);
                ln_c_k.push_back(ln_c_k_temp);
                beta_k.push_back(beta_k_temp);
            }
        }
    };


    template<
            size_t dim,
            size_t sample_size = 128,
            size_t seed = 1,
            typename RandomGenerator = std::mt19937_64,
            typename random_number_type = float>
    class WeightMinHash {
    private:
        static RandomSample<dim, sample_size, seed, RandomGenerator, random_number_type> random_sample;
    public:
        // WeightSet: W = [ w1 w2 w3 ... wn ]
        // 在 position k 处,[0,wk]范围内做hash, 此时每一个sub_element组成 (k,i) pair (i在[0,wk]范围) ,
        // 令 min_hash(k) = min { hash(k,0), ...., hash(k,wk)}, 并且 yk = arg_min_hash(k), 即 hash(k,yk) = min_hash(k).
        // 在不同的 position 处,有不同的min_hash(k),令 a = min { min_hash(0), min_hash(1), ... min_hash(n) },
        // 即 a 是整个W最终的min_hash, 并且有 k* = arg_min_a, 即 a = min_hash(k*), 也就是说 k* 是最小的min_hash(k)所在的position,
        // 我们最后将 ( k*, yk*) pair 作为 W 的 min_hash_value_vector.
        // 采样不同的参数 rk, ln_ck, beta_k, 得到不同的 min_hash_value_vector.
        using MinHashValueType = std::pair<size_t, int>; // k* , yk* (实际记录的是 k* , tk*)
        std::vector<MinHashValueType> hash_values;

        explicit WeightMinHash() {
            hash_values.resize(sample_size);
        }

        // 输入一个 dim 维度的 weight_set, 计算 weight_min_hash.
        // 每个维度上的值代表当前维度元素的权重,比如权重类型为整型,可以认为这个权重是当前维度元素的重复次数.
        // 得到 weight_set 需要预先知道所有可能元素的个数(即全集的大小),然后给每一种元素编一个号,
        // 接着对某个带重复元素的集合A,用multi-hash-set数据结构统计A中元素的个数(或者通过hash_set统计也行,但是要把权重记录在元素类型字段里),
        // 最后根据元素在全集的编号,往 weight_set_vector 写入权重的值,不存在的元素权重记为0.
        // 当然这样的预处理过程比较耗时(全集大小可能是指数规模,再加上给全集元素编号这个超级耗时的操作),
        // 并且也非常浪费空间(可能存在大量位置元素的权重是0).
        // 所以必须改进得到 weigh_set 的预处理过程,我现在的几个想法:(暂时没有实现,先完成一个naive的)
        // TODO: 1. 给全集元素编号可以通过哈希实现,只要哈希是唯一性哈希,并且全集元素的个数 < 哈希的最大范围(2^61-1),
        //       就可以编码为 pos = hash(string) % (dim). 考虑到 dim 可能不是素数,那么发生碰撞的几率就会变大,
        //       所以需要找一个大于dim的且最接近的梅森素数,这样就能减少冲突.如果发生冲突,可以用线性探测分配新的编号.
        //       当然可以写个测试代码,测试一下在全集元素是(A,T,C,G)长度为12的随机组合下,hash能不能完全不碰撞.
        //      2. 有可能出现大量全0的元素,即weight_set_vector大概率是稀疏向量,那么压缩储存就是必要的,可以用
        //      vector[(pos,weight)]记录下weight>0的部分,等于0的部分完全不用储存.
        //      进一步考虑需要封装一个 sparse_vector. 或者看 Eigen 的 sparse_vector 部分能不能直接用上.
        //      3. naive的实现完全没有numpy快,需要考虑完全向量化的并行计算,Eigen + blas + simd.
        //      Eigen可以用 min_value Array::minCoeff(row*,col*) 来实现 numpy.argmin 的效果;
        //      但Eigen对持有nan的Array做计算是undefined的,所以这里不能将weight=0的部分变成nan,
        //      而是考虑将weight=0的替换为float::max.
        //      另外Eigen好像没有 Array[value==x(condition)]=new_val 的方法,可能需要自己提供 reduce 方法(参考Eigen文档的reduction部分).
        //      现在的想法是直接用 Eigen::sparse_vector,然后提供一个[index,value]集合作为参数构造出稀疏向量,
        //      后面的计算全是向量化的,并且也不用考虑0权重的情况了,因为后面对稀疏向量的遍历计算都是在非0值上计算.
        //
        /**
         * @tparam WeightType 可以是任意的非负数值类型
         * @param weight_set_vector 权重集向量
         */
        template<typename WeightType>
        void update(const std::vector<WeightType> &weight_set_vector) {
            assert(weight_set_vector.size() == dim);
            for (size_t i = 0; i < sample_size; i++) {
                random_number_type ln_a_min = std::numeric_limits<random_number_type>::max();
                size_t k_star = 0; // k*
                random_number_type t_k_star = 0; // t_k*
                for (size_t j = 0; j < dim && weight_set_vector[j] > 0; j++) {
                    random_number_type t_k = std::floor((std::log(weight_set_vector[j]) / random_sample.r_k[i][j]) +
                                                        random_sample.beta_k[i][j]);
                    random_number_type ln_y = (t_k - random_sample.beta_k[i][j]) * random_sample.r_k[i][j];
                    random_number_type ln_a = random_sample.ln_c_k[i][j] - ln_y - random_sample.r_k[i][j];
                    if (ln_a < ln_a_min) {
                        ln_a_min = ln_a;
                        k_star = j;
                        t_k_star = t_k;
                    }
                }
                hash_values[i] = {k_star, static_cast<int>(t_k_star)};// 这里不用push_back,否则多次update同一个WeightMinHash会出错
            }
        }
    };

    template<
            size_t dim,
            size_t sample_size,
            size_t seed,
            typename RandomGenerator,
            typename random_number_type>
    RandomSample<dim, sample_size, seed, RandomGenerator, random_number_type>
            WeightMinHash<dim, sample_size, seed, RandomGenerator, random_number_type>::random_sample{};

    template<
            size_t dim,
            size_t sample_size,
            size_t seed,
            typename RandomGenerator,
            typename random_number_type>
    double weight_minhash_jaccard(const WeightMinHash<dim, sample_size, seed, RandomGenerator, random_number_type> &A,
                                  const WeightMinHash<dim, sample_size, seed, RandomGenerator, random_number_type> &B) {
        double count = 0;
        for (size_t i = 0; i < sample_size; i++) {
            if (A.hash_values[i].first == B.hash_values[i].first &&
                A.hash_values[i].second == B.hash_values[i].second) {
                count++;
            }
        }
        return count / (double) (sample_size);
    }


    // 有权重的jaccard相似度计算(使用Generalized Jaccard Index):
    // https://en.wikipedia.org/wiki/Jaccard_index#Generalized_Jaccard_similarity_and_distance
    // Example: A = { a, a, a, b, b, c } B = { a, a, b, b, b, d }
    // Generalized_Jaccard_Similarity (A,B)
    //   = (min(A_a,B_a)+min(A_b,B_b)+min(A_c,B_c)+min(A_d,B_d)) / (max(A_a,B_a)+max(A_b,B_b)+max(A_c,B_c)+max(A_d,B_d))
    //   = ( 2 + 2 + 0 + 0 ) / ( 3 + 3 + 1 + 1 ) = 4 / 8  = 0.5
    template<typename WeightType>
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
