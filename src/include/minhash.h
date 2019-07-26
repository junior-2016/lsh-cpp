//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_MINHASH_H
#define LSH_CPP_MINHASH_H

#include "lsh_cpp.h"
#include "util.h"
#include "io.h"

namespace LSH_CPP {
    constexpr static uint64_t mersenne_prime = (1ull << 61u) - 1u; // 最接近 2^64-1 的梅森素数
    constexpr static size_t max_n_permutation = 256;  // permutation 数量最多256个.

    /**
     * 生成n个随机hash函数(随机permutation),用于计算 MinHash
     */
    template<size_t Seed = 0,
            // typename Seed = std::random_device,
            typename RandomGenerator = std::mt19937_64,
            size_t n_permutation = 128>
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
         *     // use element_wise_and function to extract low 32bit.
         *     // element_wise_and 可以考虑用 sse 指令优化
         *     auto min_hash_vector_now = element_wise_and( (a_vector * di + b_vector) % mersenne_prime ,min_hash_max_range)
         *     // update min_hash_vector
         *     min_hash_vector = min(min_hash_vector_now, min_hash_vector)
         * }
         * 上面的过程中,每个 (a_vector[i],b_vector[i]) 组成一个哈希函数, 哈希过程是:
         * hash (data) = (a_vector[i] * data + b_vector[i]) % mersenne_prime.
         * 然后通过 hash(data) and 0xFFFFFFFF 把后32位拿出来作为最后的哈希值.
         * (当然你也可以直接用64位最小哈希,并且在空间都是储存的64-bit整数)
         */
#ifdef USE_SIMD
        std::vector<uint64_t, XSIMD_DEFAULT_ALLOCATOR(uint64_t) > vector_a;

        std::vector<uint64_t, XSIMD_DEFAULT_ALLOCATOR(uint64_t) > vector_b;
#else
        std::vector<uint64_t> vector_a;
        std::vector<uint64_t> vector_b;
#endif

        explicit RandomHashPermutation() {
            // TODO 考虑 static_assert(n_permutation == 128 || n_permutation == 256) 而不是提供一个范围.
            static_assert(n_permutation <= max_n_permutation);
            // Seed seed;
            RandomGenerator generator(Seed);
            std::uniform_int_distribution<uint64_t> dis_a(1, mersenne_prime);
            std::uniform_int_distribution<uint64_t> dis_b(0, mersenne_prime);
            vector_a.reserve(n_permutation);
            vector_b.reserve(n_permutation);
            for (size_t i = 0; i < n_permutation; i++) {
                vector_a.push_back(dis_a(generator));
                vector_b.push_back(dis_b(generator));
            }
        }
    };

    /**
     * @tparam HashFunc 将数据集的字符串哈希为整数的哈希函数,有三种选择:
     * 1. absl::hash 它是不确定的,也就是说,它只能保证你在某次运行程序的时候,
     * 不同的string有不同的哈希,如果你重新启动程序的时候,所有string的哈希又会重新分配一个新的了,
     * 也就是有salt(随机因子)在Hash函数里.
     * 2. std::hash , phmap::Hash 和 xxHash 都没有不确定性 (phmap::Hash内部用std::hash实现,所以完全没差别).
     * 3. 性能对比:
     * absl::Hash 略优于 xxHash (但是时间上已经非常接近),
     * xxHash 比 phmap::Hash,std::hash 快一倍.
     * 4. 虽然absl::Hash最优,但是由于absl::Hash的不确定性,默认还是使用 xxHash 作为字符串的哈希函数.
     */
    template<typename HashFunc,       // 将数据集合的元素(主要字符串类型,但也可以是其他类型)映射到整数的哈希函数
            size_t MinHashBits = 32,  // 最小哈希值的实际位数(可以是32bit或64bit,但内部储存最小哈希统一用uint64_t)
            size_t n_permutation = 128, // RandomHashPermutation中随机哈希函数的个数
            size_t Seed = 0,            // RandomHashPermutation中随机数生成器的种子
            //typename Seed = std::random_device,
            typename RandomGenerator = std::mt19937_64 // RandomHashPermutation中的随机数生成器,默认用std::mt19937_64.
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
    public:
#ifdef USE_SIMD
        std::vector<_hash_value_store_type, XSIMD_DEFAULT_ALLOCATOR(_hash_value_store_type) > hash_values;
#else
        std::vector<_hash_value_store_type> hash_values;
#endif

        explicit MinHash(HashFunc &&hash_func) : hash_func(hash_func) {
            static_assert((MinHashBits == 32 || MinHashBits == 64));
            // TODO 考虑 static_assert(n_permutation == 128 || n_permutation == 256) 而不是提供一个范围.
            static_assert(n_permutation <= max_n_permutation);
            hash_values.resize(n_permutation, _max_hash_range);
        }

        void update(const std::vector<std::string_view> &data) {
            auto values = hash_func(data);
            for (const auto &value:values) {
#ifdef USE_SIMD
                // 1. 使用组合式的simd函数,这样就不需要额外创建一个result的开销了;
                // 2. 这里假设了你的容器元素数量恰好可以被完全simd化,不会留下多余的部分,也就不需要额外提供处理函数,
                //    如果不能满足这一点,调用__simd_combine_fast__在编译期就会察觉,编译也无法通过.
                // 3. 对齐后的数据进行simd计算会更快些
                __simd_combine_fast__<xsimd::aligned_mode, n_permutation>
                        (permutation.vector_a, permutation.vector_b, hash_values, hash_values,
                         [&](const xsimd::simd_type<uint64_t> &a,
                             const xsimd::simd_type<uint64_t> &b,
                             const xsimd::simd_type<uint64_t> &c) {
                             return xsimd::min(((a * value + b) % mersenne_prime) & _max_hash_range, c);
                         });

#else
                for (size_t i = 0; i < n_permutation; i++) {
                    hash_values[i] = std::min(hash_values[i],
                                              (((value * permutation.vector_a[i] + permutation.vector_b[i]) %
                                                mersenne_prime) & _max_hash_range));
                }
#endif
            }
        }

        void update(std::string_view string) {
            auto value = hash_func(string);
#ifdef USE_SIMD
            __simd_combine_fast__<xsimd::aligned_mode, n_permutation>
                    (permutation.vector_a, permutation.vector_b, hash_values, hash_values,
                     [&](const xsimd::simd_type<uint64_t> &a,
                         const xsimd::simd_type<uint64_t> &b,
                         const xsimd::simd_type<uint64_t> &c) {
                         return xsimd::min(((a * value + b) % mersenne_prime) & _max_hash_range, c);
                     });
#else
            for (size_t i = 0; i < n_permutation; i++) {
                uint64_t ret = ((value * permutation.vector_a[i] + permutation.vector_b[i]) % mersenne_prime) &
                               _max_hash_range;
                hash_values[i] = std::min(hash_values[i], ret);
            }
#endif
        }

        [[nodiscard]] inline size_t length() const {
            return n_permutation;
        }
    };

    template<typename H,
            size_t _min_hash_bits,
            size_t _n_permutation,
            size_t _Seed,
            //typename _Seed,
            typename _RandomGenerator>
    static double jaccard_similarity
            (const MinHash<H, _min_hash_bits, _n_permutation, _Seed, _RandomGenerator> &A,
             const MinHash<H, _min_hash_bits, _n_permutation, _Seed, _RandomGenerator> &B) {
        size_t count = 0;
        // 考虑sse向量优化
        for (size_t i = 0; i < _n_permutation; i++) {
            if (A.hash_values[i] == B.hash_values[i]) count++;
        }
        return (double) count / (double) _n_permutation;
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
