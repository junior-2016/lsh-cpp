#include <utility>

#include <utility>

#include <utility>

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

    /**
     * use absl::Hash as default Hash function
     */
    template<template<typename/* KeyType */ > typename Hash, typename KeyType>
    struct hash {
        static inline size_t __hash(const KeyType &key) {
            return Hash<KeyType>{}(key);
        }

        inline size_t operator()(const KeyType &key) {
            return __hash(key);
        }

        inline std::vector<size_t> operator()(const std::vector<KeyType> &array) {
            std::vector<size_t> ret(array.size());
            for (size_t i = 0; i < array.size(); i++) {
                ret[i] = __hash(array[i]);
            }
            return ret;
        }
    };

    using DefaultStringHash = hash<absl::Hash, std::string_view>;
    using StdStringHash = hash<std::hash, std::string_view>;

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
