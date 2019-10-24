//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_MINHASH_H
#define LSH_CPP_MINHASH_H

#include "lsh_cpp.h"
#include "util.h"
#include "io.h"
#include "hash.h"
#include "lru_cache.h"

namespace LSH_CPP {
    constexpr static size_t max_n_permutation = 1024;  // permutation 最大数量.
    constexpr static size_t mersenne_prime = mersenne_prime_for_generate_64_hash;

    /**
     * 生成n个随机hash函数(随机permutation),用于计算 MinHash
     */
    template<size_t Seed,
            // typename Seed = std::random_device,
            typename RandomGenerator,
            size_t n_permutation>
    struct RandomHashPermutation {
        /**
         * permutation-random hash algorithm:
         * a_vector = [ a1 a2 a3 a4 ... an ] (size: n_permutation)
         * b_vector = [ b1 b2 b3 b4 ... bn ] (size: n_permutation)
         * Then compute min_hash of a data-vector D_m = [ d1 d2 d3 ... dm ] (if your data-vector is array
         * of string, please hash each string to integer value before computing min_hash)
         * let min_hash_max_range: (1<<32) - 1 = 0xFFFFFFFF // min_hash is uint32_t
         * let min_hash_vector = [1 1 1 ... 1] * min_hash_max_range (size: n_permutation)
         * let mersenne_prime = (1<<61)-1 // max mersenne prime that approaching 2^64-1
         * for ( di in D_m ){
         *     // Actually it is a linear permutation. Because we want to construct min_hash family :
         *     // \mathcal{H}=\left\{h_{\pi} : h_{\pi}(A)=\max _{a \in A} \pi(a)\right\} (latex),
         *     // so we should construct \pi to hash all word to Number field U = { 0,1, ... ,u } (u large enough to make collisions unlikely.)
         *     // Finally we use linear permutation : \pi(x) = ( a * x + b ) % u, a and b random, u is max(U).
         *     // Here we use U = {0, 1, ... , mersenne_prime }, but we use element_wise_and function to extract low 32bit finally
         *     // (It is equally mod min_hash_max_range).
         *     auto min_hash_vector_now = element_wise_and( (a_vector * di + b_vector) % mersenne_prime ,min_hash_max_range)
         *     // update min_hash_vector
         *     min_hash_vector = min(min_hash_vector_now, min_hash_vector)
         * }
         * 上面的过程中,每个 (a_vector[i],b_vector[i]) 组成一个哈希函数, 哈希过程是:
         * hash (data) = (a_vector[i] * data + b_vector[i]) % mersenne_prime.
         * 然后通过 hash(data) and 0xFFFFFFFF 把后32位拿出来作为最后的哈希值.
         * (当然你也可以直接用64位最小哈希,并且在空间都是储存的64-bit整数)
         */
        std::vector<uint64_t> vector_a;
        std::vector<uint64_t> vector_b;

        explicit RandomHashPermutation() {
            static_assert(n_permutation <= max_n_permutation);
            // Seed seed;
            RandomGenerator generator(Seed);
            std::uniform_int_distribution<uint64_t> dis_a(1, mersenne_prime - 1); // range [1,mersenne_prime-1]
            std::uniform_int_distribution<uint64_t> dis_b(0, mersenne_prime - 1); // range [0,mersenne_prime-1]
            vector_a.reserve(n_permutation);
            vector_b.reserve(n_permutation);
            for (size_t i = 0; i < n_permutation; i++) {
                vector_a.push_back(dis_a(generator));
                vector_b.push_back(dis_b(generator));
            }
        }
    };

    /**
     *
     * @tparam HashFunc 将数据集合的元素(主要字符串类型,但也可以是其他类型)映射到整数的哈希函数,默认用 XXStringViewHash32,
     * 可以在构造函数传递第一个参数来指定其他的HashFunc
     * @tparam MinHashBits 最小哈希值的实际位数(可以是32bit或64bit,但内部储存最小哈希统一用uint64_t)
     * @tparam n_permutation RandomHashPermutation中随机哈希函数的个数
     * @tparam Seed RandomHashPermutation中随机数生成器的种子
     * @tparam RandomGenerator RandomHashPermutation中的随机数生成器,默认用std::mt19937_64
     */
    template<typename HashFunc = XXStringViewHash32,
            size_t MinHashBits = 32,
            size_t n_permutation = 128,
            size_t Seed = 1,
            //typename Seed = std::random_device,
            typename RandomGenerator = std::mt19937_64
    >
    class MinHash {
    public:
        // min_hash value 的储存类型为uint64_t
        using _hash_value_store_type = uint64_t;

        // min_hash 的实际类型取决于MinHashBits,所以这里的哈希值最大值也是由MinHashBits决定
        const uint64_t _max_hash_range = HashValueType<MinHashBits>::max_hash_range;
    private:
        HashFunc hash_func;
        static RandomHashPermutation<Seed, RandomGenerator, n_permutation> permutation;

        // 用 MapArray 将 std::vector 包装为 Eigen::Array,
        // 并且所有对 MapArray 的修改都直接反映在 std::vector 的数据上.
        using Array = Eigen::Array<uint64_t, n_permutation, 1>;
        using MapArray = Eigen::Map<Array>;

        // 对Eigen的对齐问题,见:https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html
        using Cache = lru_cache<uint64_t, Array, phmap::Hash<uint64_t>, phmap::EqualTo<uint64_t>,
                Eigen::aligned_allocator, phmap::flat_hash_map>;
        static Cache cache; // 全局 lru_cache 加速 update 计算.

    public:
        std::vector<_hash_value_store_type> hash_values;

        explicit MinHash(HashFunc &&hash_func = XXStringViewHash32{}) : hash_func(hash_func) {
            static_assert((MinHashBits == 32 || MinHashBits == 64));
            static_assert(n_permutation <= max_n_permutation);
            hash_values.resize(n_permutation, _max_hash_range);
        }

        /**
         * 1.
         * update (T) 兼容任意类型数据,这样就可以轻松扩展到用户自定义数据集合的MinHash计算.
         * 只要特化 LSH_CPP::xxHash<T> 对这种类型进行哈希即可.
         * example:
         * struct K_mer { std::string_view value; Other_member ...  } // 注意这里需要针对 flat_hash_set 给出 K_mer 的哈希函数.
         * phmap::flat_hash_set<K_mer> k_mer_set = { ..... } // init set
         * template <>
         * struct xxHash<K_mer> { // partial specialization for K_mer
         *     uint64_t operator()(const K_mer & k_mer) { return xx_hash<64>(k_mer.value); }
         * }
         *
         * @tparam T 集合数据类型
         * @param val 集合元素
         * update() 方法使用Eigen加速,需要注意一些问题:
         * 1.用 Eigen::Map<Array<Type,Row,Col>> eigen_array(vector.data(),size)与已有的std::vector成员无缝对接,不用破坏原来的代码架构.
         * 2.用上面的方法,如果 eigen_array 的数据被修改,那么vector内的数据也会一样被修改,不需要将eigen_array的数据再次迁移到vector上.
         */
        template<typename T>
        void update(const T &val) {
            MapArray hash_values_array(hash_values.data(), n_permutation);
            MapArray vector_a(permutation.vector_a.data(), n_permutation);
            MapArray vector_b(permutation.vector_b.data(), n_permutation);
            auto value = hash_func(val); // compress element(may be shingling) to integer value.
            if (auto ret = cache.get(value); ret.has_value()) {
                hash_values_array = hash_values_array.min(*ret);
            } else {
                Array temp = (vector_a * value + vector_b).unaryExpr(
                        [&](const auto x) -> uint64_t { return (x % mersenne_prime) & _max_hash_range; });
                cache.put(value, temp);
                hash_values_array = hash_values_array.min(temp);
            }
        }

        /**
         * 兼容任意类型数据集合的MinHash计算,用法同上.不同在于这里是对一个集合对所有元素进行操作,然后更新 min_hash vector.
         */
        template<typename T>
        void update(const HashSet<T> &data_set) {
            MapArray hash_values_array(hash_values.data(), n_permutation);
            MapArray vector_a(permutation.vector_a.data(), n_permutation);
            MapArray vector_b(permutation.vector_b.data(), n_permutation);
            for (const auto &data:data_set) {
                auto value = hash_func(data);
                if (auto ret = cache.get(value); ret.has_value()) {
                    hash_values_array = hash_values_array.min(*ret);
                } else {
                    Array temp = (vector_a * value + vector_b).unaryExpr(
                            [&](const auto x) -> uint64_t { return (x % mersenne_prime) & _max_hash_range; });
                    cache.put(value, temp);
                    hash_values_array = hash_values_array.min(temp);
                }
            }
        }

        [[nodiscard]]constexpr inline size_t length() const {
            return n_permutation;
        }
    };

    constexpr size_t max_cache_size = 10000;
    template<typename HashFunc,
            size_t MinHashBits,
            size_t n_permutation,
            size_t Seed,
            typename RandomGenerator
    >
    typename MinHash<HashFunc, MinHashBits, n_permutation, Seed, RandomGenerator>::Cache
            MinHash<HashFunc, MinHashBits, n_permutation, Seed, RandomGenerator>::cache{max_cache_size};

    // 下面的 jaccard_similarity 计算公式是通过 min_hash_value_vector 估计得到的,
    // 只是从概率上和真实的jaccard_similarity相等,和真实的jaccard_similarity依然存在偏差.
    template<typename H,
            size_t _min_hash_bits,
            size_t _n_permutation,
            size_t _Seed,
            //typename _Seed,
            typename _RandomGenerator>
    static double minhash_jaccard_similarity
            (MinHash<H, _min_hash_bits, _n_permutation, _Seed, _RandomGenerator> &A,
             MinHash<H, _min_hash_bits, _n_permutation, _Seed, _RandomGenerator> &B) {
        using Array = Eigen::Array<uint64_t, _n_permutation, 1>;
        using MapArray = Eigen::Map<Array>;
        // 用 vector.data() 构造 MapArray 时要保证 data() 返回的是 T* 而不是 const T*, 因为 MapArray 可能使用带有副作用的操作.
        // 不过这里的函数不会对 MapArray 做任何修改,所以不用担心 MapArray 影响原数据的值.
        MapArray a_array(A.hash_values.data(), _n_permutation);
        MapArray b_array(B.hash_values.data(), _n_permutation);
        int count = (a_array == b_array)
                .select(one_eigen_array<uint64_t, _n_permutation>, zero_eigen_array<uint64_t, _n_permutation>)
                .sum();
        return (double) count / (double) _n_permutation;
    }

    // jaccard相似度计算(针对无权重集合): jaccard_similarity = (A intersection B) / (A union B). 使用HashSet,复杂度O(N).
    template<typename T>
    double jaccard_similarity(const HashSet<T> &A, const HashSet<T> &B) {
        size_t count = 0;
        for (const auto &a:A) { if (B.contains(a)) count++; }
        return (double) count / (double) (A.size() + B.size() - count);
    }

    template<typename H,
            size_t _min_hash_bits,
            size_t _n_permutation,
            size_t _Seed,
            //typename _Seed,
            typename _RandomGenerator>
    static void print_minhash_table(
            const MinHash<H, _min_hash_bits, _n_permutation, _Seed, _RandomGenerator> &m1,
            const MinHash<H, _min_hash_bits, _n_permutation, _Seed, _RandomGenerator> &m2) {
        printf("\t M_1 \t M_2 \t Eq? \n");
        for (size_t i = 0; i < _n_permutation; i++) {
            printf("\t %ul \t %ul \t %d \n", m1.hash_values[i], m2.hash_values[i],
                   (m1.hash_values[i] == m2.hash_values[i]));
        }
    }

    template<typename HashFunc,
            size_t MinHashBits,
            size_t n_permutation,
            size_t Seed,
            //typename Seed,
            typename RandomGenerator>
    RandomHashPermutation<Seed, RandomGenerator, n_permutation>
            MinHash<HashFunc, MinHashBits, n_permutation, Seed, RandomGenerator>::permutation{};
}
#endif //LSH_CPP_MINHASH_H
