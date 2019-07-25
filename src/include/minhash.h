//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_MINHASH_H
#define LSH_CPP_MINHASH_H

#include "lsh_cpp.h"
#include "util.h"
#include "io.h"

namespace LSH_CPP {
    const uint64_t mersenne_prime = (1ull << 61u) - 1u; // 最接近 2^64-1 的梅森素数
    const size_t max_n_permutation = 256;  // permutation 数量最多256个

    /**
     * 生成n个随机hash函数(随机permutation),用于计算 MinHash
     */
    template<size_t Seed = 1,
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
        std::vector<uint64_t> vector_a;

        std::vector<uint64_t> vector_b;

        explicit RandomHashPermutation() {
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
    template<typename HashFunc,
            size_t MinHashBits = 32,
            size_t n_permutation = 128,
            size_t Seed = 0,
            //typename Seed = std::random_device,
            typename RandomGenerator = std::mt19937_64>
    class MinHash {
    public:
        using _hash_value_store_type = uint64_t; // min_hash value 的储存类型为uint64_t
        // min_hash 的实际类型取决于MinHashBits,所以这里的哈希值最大值也是由MinHashBits决定
        const size_t _max_hash_range = HashValueType<MinHashBits>::max_hash_range;
    private:
        HashFunc hash_func;
        static RandomHashPermutation<Seed, RandomGenerator, n_permutation> permutation;
    public:
        std::vector<_hash_value_store_type> hash_values;

        explicit MinHash(HashFunc &&hash_func) : hash_func(hash_func) {
            static_assert((MinHashBits == 32 || MinHashBits == 64) && n_permutation <= max_n_permutation);
            hash_values.resize(n_permutation, _max_hash_range);
        }

        void update(std::string_view string) {
            auto value = hash_func(string);
            // 考虑sse向量优化
            for (size_t i = 0; i < n_permutation; i++) {
                uint64_t ret = ((value * permutation.vector_a[i] + permutation.vector_b[i]) % mersenne_prime) &
                               _max_hash_range;
                hash_values[i] = std::min(hash_values[i], ret);
            }
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
