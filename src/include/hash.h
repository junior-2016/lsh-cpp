//
// Created by junior on 2019/7/27.
//

#ifndef LSH_CPP_HASH_H
#define LSH_CPP_HASH_H

#include "lsh_cpp.h"

namespace xxh {
    // 由于 xxhash 不支持 string_view,所以我仿照它的 string 版本加入了对 string_view 的支持.
    template<size_t N, typename T>
    hash_t <N>
    xxhash(const std::basic_string_view<T> &input, hash_t <N> seed = 0, endianness endian = endianness::unspecified) {
        static_assert(!(N != 32 && N != 64), "You can only call xxhash in 32 or 64 bit mode.");
        return detail::endian_align<N>(static_cast<const void *>(input.data()), input.length() * sizeof(T), seed,
                                       mem_ops::get_endian(endian),
                                       mem_ops::get_alignment<N>(static_cast<const void *>(input.data())));
    }
}

namespace LSH_CPP {
    // mersenne prime https://en.wikipedia.org/wiki/Mersenne_prime
    // mersenne prime format : 2^n - 1 ( n = { 2, 3, 5, 7, 13, 17, 19, 31, 61, ...} )
    // 用于最后生成64bit-hash的mersenne_prime: uint64_t hash_value =  hash_func(data) % mersenne_prime_for_generate_64_hash
    constexpr static uint64_t mersenne_prime_for_generate_64_hash = (1ull << 61u) - 1u;
    // 用于最后生成32bit-hash的mersenne_prime: 用法同上.
    constexpr static uint32_t mersenne_prime_for_generate_32_hash = (1ull << 31u) - 1u;

    template<size_t N>
    struct HashValueType {
    };
    template<>
    struct HashValueType<32> {
        using type = uint32_t;
        static constexpr uint64_t max_hash_range = std::numeric_limits<uint32_t>::max(); // 0x00000000FFFFFFFF
    };
    template<>
    struct HashValueType<64> {
        using type  = uint64_t;
        static constexpr uint64_t max_hash_range = std::numeric_limits<uint64_t>::max(); // 0xFFFFFFFFFFFFFFFF
    };

    /**
     * @tparam Hash : 这里要求 Hash{}(key) 返回的必须是64bit的哈希值!
     * @tparam Bits : 实际返回的哈希值位数.
     * 如果是64bit,就直接返回 Hash{}(key) 的结果;
     * 如果是32bit,就返回 Hash{}(key) % mersenne_prime_for_generate_32_hash
     */
    template<template<typename/* KeyType */ > typename Hash, typename KeyType, size_t Bits = 64>
    struct hash {
        // hash single key
        inline typename HashValueType<Bits>::type operator()(const KeyType &key) const {
            static_assert(Bits == 32 || Bits == 64);
            if constexpr (Bits == 32) {
                return (Hash<KeyType>{}(key)) % mersenne_prime_for_generate_32_hash;
            } else {
                return Hash<KeyType>{}(key);
            }
        }

        // hash whole vector
        inline typename HashValueType<Bits>::type operator()(const std::vector<KeyType> &key) const {
            static_assert(Bits == 32 || Bits == 64);
            if constexpr (Bits == 32) {
                return (Hash<KeyType>{}(key)) % mersenne_prime_for_generate_32_hash;
            } else {
                return Hash<KeyType>{}(key);
            }
        }

        // hash vector in range
        inline typename HashValueType<Bits>::type
        operator()(const std::vector<KeyType> &key, const std::pair<size_t, size_t> &range) const {
            static_assert(Bits == 32 || Bits == 64);
            if constexpr (Bits == 32) {
                return (Hash<KeyType>{}(key, range)) % mersenne_prime_for_generate_32_hash;
            } else {
                return Hash<KeyType>{}(key, range);
            }
        }
    };

    // 注意下面的函数,参数是 const hash<..> & hash_fun,也就是一个const对象,
    // 所以它只能调用 hash 类的const函数, 要把hash类所有的函数都标记为const才能正确调用
    // element_wise_hash for sequence container
    template<
            template<typename /*Element*/, typename /*Alloc*/> typename Container,
            template<typename /* KeyType */ > typename Hash,
            typename KeyType,
            typename Alloc = std::allocator<KeyType>,
            size_t Bits>
    std::vector<typename HashValueType<Bits>::type> element_wise_hash(
            const hash<Hash, KeyType, Bits> &hash_func,
            const Container<KeyType, Alloc> &container) {
        std::vector<typename HashValueType<Bits>::type> ret;
        ret.reserve(container.size()); // 使用 reserve() 提前分配空间但是不初始化
        std::for_each(std::begin(container), std::end(container),
                      [&](const KeyType &key) { ret.push_back(hash_func(key)); });
        return ret;
    }

    // element_wise_hash for ordered set(balanced binary search tree implement set/multiple_set)
    template<
            template<typename /* Key */, typename /* Compare */, typename /* Alloc */> typename Set,
            template<typename /* KeyType */ > typename Hash,
            typename KeyType,
            typename Compare = std::less<KeyType>,
            typename Alloc = std::allocator<KeyType>,
            size_t Bits>
    std::vector<typename HashValueType<Bits>::type>
    element_wise_hash(const hash<Hash, KeyType, Bits> &hash_func,
                      const Set<KeyType, Compare, Alloc> &set) {
        std::vector<typename HashValueType<Bits>::type> ret;
        ret.reserve(set.size()); // 使用 reserve() 提前分配空间但是不初始化
        std::for_each(std::begin(set), std::end(set),
                      [&](const KeyType &key) { ret.push_back(hash_func(key)); });
        return ret;
    }


    // element_wise_hash for unordered_set(hash set/multiple_set)
    template<
            template<typename /* Key */, typename /* Hash */, typename /* Eq */, typename /* Alloc */ > typename HashSet,
            template<typename /* KeyType */ > typename Hash,
            typename KeyType,
            typename _Hash,
            typename _Eq,
            typename _Alloc,
            size_t Bits>
    std::vector<typename HashValueType<Bits>::type> element_wise_hash(
            const hash<Hash, KeyType, Bits> &hash_func,
            const HashSet<KeyType, _Hash, _Eq, _Alloc> &hash_set) {
        std::vector<typename HashValueType<Bits>::type> ret;
        ret.reserve(hash_set.size()); // 使用 reserve() 提前分配空间但是不初始化
        std::for_each(std::begin(hash_set), std::end(hash_set),
                      [&](const KeyType &key) { ret.push_back(hash_func(key)); });
        return ret;
    }

    template<typename T>
    struct xx_Hash {

    };

    template<typename T>
    struct xx_Hash<std::basic_string_view<T>> {
        inline uint64_t operator()(const std::basic_string_view<T> &stringView) {
            return xxh::xxhash<64>(stringView);
        }

        inline uint64_t operator()(const std::vector<std::basic_string_view<T>> &) {
            // TODO: 空实现. 后面再补充,因为这里vector里面储存的是 string_view,并不是基本的数据类型比如int/char/double,
            //  所以不能直接传递整个vector
            return 0;
        }

        inline uint64_t operator()(const std::vector<std::basic_string_view<T>> &,
                                   const std::pair<size_t, size_t> &) {
            return 0;
        }
    };

    template<typename T>
    struct xx_Hash<std::basic_string<T>> {
        inline uint64_t operator()(const std::basic_string<T> &string) { return xxh::xxhash<64>(string); }

        inline uint64_t operator()(const std::vector<std::basic_string<T>> &) { return 0; }

        inline uint64_t operator()(const std::vector<std::basic_string<T>> &,
                                   const std::pair<size_t, size_t> &) { return 0; }
    };

    template<>
    struct xx_Hash<K_mer> {
        inline uint64_t operator()(const K_mer &k_mer) { return xxh::xxhash<64>(k_mer.value); }

        inline uint64_t operator()(const std::vector<K_mer> &) { return 0; }

        inline uint64_t operator()(const std::vector<K_mer> &, const std::pair<size_t, size_t> &) { return 0; }
    };

    template<>
    struct xx_Hash<uint64_t> {
        inline uint64_t operator()(const uint64_t &integer) {
            const void *integer_address = static_cast<const void *>(&integer);
            return xxh::xxhash<64>(integer_address, sizeof(uint64_t));
        }

        inline uint64_t operator()(const std::vector<uint64_t> &integer_vector) {
            return xxh::xxhash<64>(integer_vector); // 具体实现直接看 xxh::xxhash(std::vector<T>& ) 的代码就清楚了.
        }

        inline uint64_t operator()(const std::vector<uint64_t> &integer_vector,
                                   const std::pair<size_t, size_t> &range) {
            auto[start, end] = range;
            // assert will omit on release build mode.
            assert(0 <= start && start < end && end <= integer_vector.size());
            const uint64_t *address = integer_vector.data() + start;
            const size_t interval = end - start;
            // 注意 xxh::xxhash 内部实现将数据解释为uint8_t*, 所以这里必须用: sizeof(uint64_t) * interval
            return xxh::xxhash<64>(static_cast<const void *>(address), sizeof(uint64_t) * interval);
        }
    };

// string hash function
    using XXStringViewHash64 = hash<xx_Hash, std::string_view, 64>; // char_t
    using XXStringViewHash32 = hash<xx_Hash, std::string_view, 32>; // char_t
    using XXWStringViewHash64 = hash<xx_Hash, std::wstring_view, 64>; // wchar_t
    using XXWStringViewHash32 = hash<xx_Hash, std::wstring_view, 32>; // wchar_t

    using XXStringHash64 = hash<xx_Hash, std::string, 64>;
    using XXStringHash32 = hash<xx_Hash, std::string, 32>;
    using XXWStringHash64 = hash<xx_Hash, std::wstring, 64>;
    using XXWStringHash32 = hash<xx_Hash, std::wstring, 32>;

// integer hash function
    using XXUInt64Hash64 = hash<xx_Hash, uint64_t, 64>;
    using XXUInt64Hash32 = hash<xx_Hash, uint64_t, 32>;

//    using StdStringViewHash64 = hash<std::hash, std::string_view, 64>;
//    using StdStringViewHash32 = hash<std::hash, std::string_view, 32>;

}

#endif //LSH_CPP_HASH_H
