//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_MINHASH_H
#define LSH_CPP_MINHASH_H

#include "lsh_cpp.h"

namespace LSH_CPP {
    const uint64_t mersenne_prime = (1ull << 61u) - 1u; // 最接近 2^64-1 的梅森素数
    const size_t max_n_permutation = 256;  // permutation 数量最多256个

    /**
     * 生成n个随机hash函数(随机permutation),用于计算 MinHash
     */
    template<typename Seed = std::random_device, typename RandomGenerator = std::mt19937_64,
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
         * (当然你也可以直接用64位最小哈希,不过64bit内存开销就会比较大)
         */
        std::vector<uint64_t> vector_a;

        std::vector<uint64_t> vector_b;

        explicit RandomHashPermutation() {
            static_assert(n_permutation <= max_n_permutation);
            Seed seed;
            RandomGenerator generator(seed());
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

    template<size_t N>
    struct HashValueType {
    };
    template<>
    struct HashValueType<32> {
        using type = uint32_t;
        const uint64_t max_hash_range = std::numeric_limits<uint32_t>::max(); // 0x00000000FFFFFFFF
    };
    template<>
    struct HashValueType<64> {
        using type = uint64_t;
        const uint64_t max_hash_range = std::numeric_limits<uint64_t>::max(); // 0xFFFFFFFFFFFFFFFF
    };


    template<size_t MinHashBits = 32>
    class MinHash {
    public:
        using HashValueType = typename HashValueType<MinHashBits>::type; // 内部储存的最小哈希类型为uint32
    private:

    public:
        explicit MinHash() {
            static_assert(MinHashBits == 32 || MinHashBits == 64);

        }
    };

    void test() {

    }
}
#endif //LSH_CPP_MINHASH_H
