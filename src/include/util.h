//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_UTIL_H
#define LSH_CPP_UTIL_H

#include "lsh_cpp.h"

namespace LSH_CPP {
    /**
     * Non-capture lambda can be transfer to function pointer directly.
     * @tparam Func lambda expression type : [](double,void*)->double{...}
     * @param f : integration function for gsl_function
     * @param range : integration interval from range.first to range.second
     * @param params : parameters array for gsl_function
     * gsl integration api:
     * int gsl_integration_qags (const gsl_function * f, double a, double b, double epsabs,
     * double epsrel, size_t limit, gsl_integration_workspace * workspace,
     * double * result, double * abserr)
     */
    template<typename Func>
    double numerical_integration(Func &&f, std::pair<double, double> range, double *params) {
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
        double result, error;
        gsl_function F;
        F.function = static_cast<double (*)(double, void *)> (f);
        F.params = reinterpret_cast<void *>(params);
        gsl_integration_qags(&F, range.first, range.second, 0, 1e-7, 1000, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;
    }

    inline double false_negative_probability(double x, void *param) {
        auto *parameter = reinterpret_cast<double *>(param);
        double b = parameter[0], r = parameter[1];
        return pow(1.0 - pow(x, r), b);      // ∫ (threshold -> 1.0) 1 - (1 - (1 - s^r)^b)
    }

    inline double false_positive_probability(double x, void *param) {
        auto *parameter = reinterpret_cast<double *>(param);
        double b = parameter[0], r = parameter[1];
        return 1.0 - pow(1.0 - pow(x, r), b); // ∫ (0.0 -> threshold) 1 - (1 - s^r)^b
    }

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

    template<template<typename/* KeyType */ > typename Hash, typename KeyType, size_t Bits = 64>
    struct hash {
        static inline typename HashValueType<Bits>::type __hash(const KeyType &key) {
            static_assert(Bits == 32 || Bits == 64);
            return Hash<KeyType>{}(key);
        }

        inline typename HashValueType<Bits>::type operator()(const KeyType &key) {
            static_assert(Bits == 32 || Bits == 64);
            return __hash(key);
        }

        inline std::vector<typename HashValueType<Bits>::type> operator()(const std::vector<KeyType> &array) {
            static_assert(Bits == 32 || Bits == 64);
            std::vector<typename HashValueType<Bits>::type> ret(array.size());
            for (size_t i = 0; i < array.size(); i++) {
                ret[i] = __hash(array[i]);
            }
            return ret;
        }
    };

    template<typename T>
    struct xx_Hash {

    };

    // wrap xx_hash library to LSH_CPP::xx_Hash. only implement for string_view and string.
    // use xxh::hash64_t(actual type is uint64_t) as hash type
    // (because xxh::hash64_t is faster than xxh::hash32_t on x86_64)
    template<>
    struct xx_Hash<std::string_view> {
        inline uint64_t operator()(std::string_view stringView) {
            return xxh::xxhash<64>(stringView.data(), stringView.size());
        }
    };

    template<>
    struct xx_Hash<std::string> {
        inline uint64_t operator()(const std::string &string) {
            return xxh::xxhash<64>(string);
        }
    };

    template<>
    struct xx_Hash<uint64_t> {
        // 下面的接口只适用于对整一个vector做哈希.
        // 如果要哈希vector的某个部分,应该拿: auto data_pointer = vector.data()+bias; auto size = interval;
        // 然后调用下面的 row-pointer with size 的接口.
        inline uint64_t operator()(const std::vector<uint64_t> &int_stream) {
            return xxh::xxhash<64>(int_stream); // this->operator()(int_stream.data(),int_stream.size());
        }

        // row-pointer with size
        inline uint64_t operator()(const uint64_t *int_stream, size_t size) {
            return xxh::xxhash<64>(int_stream, size);
        }
    };

    // string hash function
    using XXStringHash64 = hash<xx_Hash, std::string_view, 64>;
    using XXStringHash32 = hash<xx_Hash, std::string_view, 32>;
    using DefaultStringHash64 = hash<std::hash, std::string_view, 64>;
    using DefaultStringHash32 = hash<std::hash, std::string_view, 32>;
    using AbslStringHash64 = hash<absl::Hash, std::string_view, 64>;
    using AbslStringHash32 = hash<absl::Hash, std::string_view, 32>;

    // integer hash function
    using XXUInt64StreamHash64 = xx_Hash<uint64_t>;

    // integer stream hash vector<uint64_t> between [start,end)
    // usage case: auto ret = int_stream_hash<XXUInt64StreamHash64>( /* int_stream */{1,2,3,4}, /* range */{ 0, 4 })
    template<typename Hash>
    inline uint64_t int_stream_hash(const std::vector<uint64_t> &int_stream, const std::pair<size_t, size_t> &range) {
        auto[start, end] = range;
        // assert will omit on release build mode.
        assert(start >= 0 && start < int_stream.size() && end > start && end <= int_stream.size());
        return Hash{}(int_stream.data() + start, end - start);
    }

    /**
     *
     * split string to k_mer.
     * @param string 源字符串,考虑用 string_view 优化,因为string_view只储存一个data pointer和一个字符串长度,
     * 这样使用substr()将不会有任何的内存开销,速度也相对更快,只需要保证源字符串自始至终都存在(不被析构)就好.
     * 参考: https://stackoverflow.com/questions/42921998/c-efficiently-get-substring-of-string-with-index
     * 经过比较,使用std::string_view比直接用std::string快了一倍!
     * 但是string_view毕竟只是持有原字符串的地址,所以一旦发生源字符串的析构就会导致严重错误,比如我下面的调用代码:
     * int main(){
     *     std::vector<std::string_view> ret;
     *     {
     *        std::string s = "abcdefghijklmnopqrstuvwxyz";
     *        ret = split_k_mer<32>(s, 3);
     *     } // 离开作用域后源字符串s被析构
     *     for (const auto& item:ret){
     *        std::cout<<item<<" ";
     *        // 输出ret的元素,这时string_view持有的地址所在的内容已经析构或者被别的方式处理了,
     *        // 此时输出结果是 undefined behavior (未定义行为).
     *     }
     * }
     * 因此在使用string_view时要格外小心源字符串的生命周期
     */
    // 更快的k_mer_split,但是需要考虑源字符串的生命周期.
    // 使用std::string_view做参数不需要引用,因为开销已经足够小,加不加引用没有区别
    std::vector<std::string_view> split_k_mer_fast(std::string_view string, size_t k) {
        if (k >= string.size()) { return {string}; }
        size_t N = string.size() - k + 1;
        std::vector<std::string_view> result;
        result.reserve(N);
        for (size_t i = 0; i < N; i++) {
            result.push_back(string.substr(i, k));
        }
        return result;
    }

#ifdef USE_SIMD
    template<typename Tag,     // 选择 simd::aligned_mode 还是 std::unaligned_mode
            size_t size,       // 容器数量编译器确定
            typename SimdFunc, // 用于simd处理的函数
            typename Func,     // 用于非simd处理的函数
            template<typename /*T*/, typename /*allocator*/> typename C, // 容器类型
            typename T,          // 容器元素类型
            typename Allocator   // 容器Allocator类型
    >
    void __simd__(const C<T, Allocator> &a,
                  const C<T, Allocator> &b,
                  C<T, Allocator> &res,
                  SimdFunc &&simd_func,
                  Func &&func) {
        // T类型在simd指令的实际类型为simd_type ; 大小为 simd_size(一个simd_size可以存几个基本的T类型元素)
        using simd_type = xsimd::simd_type<T>;
        constexpr std::size_t simd_size = simd_type::size; // 编译期决定值

        // 计算for循环每次加一个simd_size后最大的上界.当继续加size_size超过这个上界时就结束循环
        constexpr std::size_t vec_size = size - size % simd_size; // 编译期决定值

        for (std::size_t i = 0; i < vec_size; i += simd_size) {
            simd_type _a = xsimd::load_simd(&a[i], Tag());
            simd_type _b = xsimd::load_simd(&b[i], Tag());
            simd_type _res = simd_func(_a, _b); // 可simd部分执行simd_func
            xsimd::store_simd(&res[i], _res, Tag());
        }

        if constexpr (vec_size < size) {
            // 用 if constexpr 优化: 如果 vec_size >= size, 下面的代码直接连编译都不用.
            for (std::size_t i = vec_size; i < size; ++i) {
                // 计算剩下不能被simd的部分,使用普通的func
                res[i] = func(a[i], b[i]);
            }
        }
    }

    template<typename Tag,     // simd::aligned_mode or std::unaligned_mode
            size_t size,       // 容器数量编译器确定.
            typename SimdFunc, // 用于simd处理的函数,可以同时处理三个simd元素,即 ret = simd_func(a,b,c)
            template<typename /*T*/, typename /*allocator*/> typename C, // 容器类型
            typename T,
            typename Allocator>
    void __simd_combine_fast__(
            const C<T, Allocator> &a,
            const C<T, Allocator> &b,
            const C<T, Allocator> &c,
            C<T, Allocator> &res,
            SimdFunc &&simd_func) {
        using simd_type = xsimd::simd_type<T>;
        constexpr std::size_t simd_size = simd_type::size; // 编译期决定值
        constexpr std::size_t vec_size = size - size % simd_size; // 编译期决定值
        static_assert(vec_size == size);   // 为了加速计算,这里要求提供的容器数量必须刚刚好可以完全向量化,不留多余部分

        for (std::size_t i = 0; i < vec_size; i += simd_size) {
            simd_type _a = xsimd::load_simd(&a[i], Tag());
            simd_type _b = xsimd::load_simd(&b[i], Tag());
            simd_type _c = xsimd::load_simd(&c[i], Tag());
            simd_type _res = simd_func(_a, _b, _c); // simd部分执行simd_func
            xsimd::store_simd(&res[i], _res, Tag());
        }
    }
#endif

    // 废弃使用 K_mer
//    struct K_mer {
//        const std::string origin_string;
//        std::vector<std::string_view> sub_strings;
//
//        K_mer(std::string origin_string, std::vector<std::string_view> sub_strings) :
//                origin_string(std::move(origin_string)), sub_strings(std::move(sub_strings)) {}
//    };
}
#endif //LSH_CPP_UTIL_H
