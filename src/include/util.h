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
 * xxhash 32 bit 和 xxhash 64 bit 与 std::hash 的性能:
 * 64bit 机器下:
 *    xxhash 64 比 std::hash 快一倍;
 *    xxhash 32 和 std::hash 打平.
 * 32bit 机器不清楚,但是根据xxhash的文档,xxhash 32在32bit机器比在64bit机器更快,所以理论上,32bit机器上性能应该是:
 *    xxhash 64 优于 xxhash 32 优于 std::hash.
 * 应该默认采用 xxhash 64 以获得更好的性能.
 */
    inline xxh::hash32_t hash32(const std::string &string) {
        return xxh::xxhash<32>(string);
    }

    inline xxh::hash64_t hash64(const std::string &string) {
        return xxh::xxhash<64>(string);
    }

    inline size_t std_hash(const std::string &string) {
        return std::hash<std::string>{}(string);
    }

    /**
    * split string to k_mer.
    * @tparam ThreadNumber 线程数量
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
    * 为了解决这个问题,需要引入一个数据结构,将源字符串和string_view的子串数组打包起来即可,打包后上面的代码将无法通过编译,
    * 因为打包后的K_mer结构体只能初始化一次且无法进行任何方式的二次赋值.
    *
    * 不过打包后发现性能略微有些下降.如果直接用string_view作为参数,速度可以更快,所以我提供了另一个使用string_view做参数的版本,
    * 需要注意的是这个版本需要时刻考虑源字符串的生命周期.
    * @param k k_mer子串长度
    * @return k_mer集合
    */
    struct K_mer {
        const std::string origin_string;
        std::vector<std::string_view> sub_strings;
    };

    template<size_t ThreadNumber>
    K_mer split_k_mer(const std::string &string, size_t k) {
        std::string_view view = string;
        if (k >= string.size()) { return K_mer{string, {view}}; }
        size_t N = string.size() - k + 1;
        K_mer result{string, std::vector<std::string_view>(N)};
        if (N <= 4 * ThreadNumber) { // k_mer数量不够多的话,依然用单线程.
            for (size_t i = 0; i < N; i++) {
                result.sub_strings[i] = view.substr(i, k);
            }
            return result;
        }
        std::vector<std::thread> threads(ThreadNumber);
        size_t step = N / ThreadNumber;
        for (size_t i = 0; i < ThreadNumber; i++) {
            threads[i] = std::thread([&](size_t id) {
                size_t begin = id * step, end = (id == ThreadNumber - 1) ? N : (id + 1) * step;
                for (size_t j = begin; j < end; j++) {
                    result.sub_strings[j] = view.substr(j, k);
                }
            }, i);
        }
        for (auto &thread:threads) {
            thread.join();
        }
        return result;
    }

    // 更快的k_mer_split,但是需要考虑源字符串的生命周期.
    // 使用std::string_view做参数不需要引用,因为开销已经足够小,加不加引用没有区别
    template<size_t ThreadNumber>
    std::vector<std::string_view> split_k_mer_fast(std::string_view string, size_t k) {
        if (k >= string.size()) { return {string}; }
        size_t N = string.size() - k + 1;
        std::vector<std::string_view> result(N);
        if (N <= 4 * ThreadNumber) { // k_mer数量不够多的话,依然用单线程.
            for (size_t i = 0; i < N; i++) {
                result[i] = string.substr(i, k);
            }
            return result;
        }
        std::vector<std::thread> threads(ThreadNumber);
        size_t step = N / ThreadNumber;
        for (size_t i = 0; i < ThreadNumber; i++) {
            threads[i] = std::thread([&](size_t id) {
                size_t begin = id * step, end = (id == ThreadNumber - 1) ? N : (id + 1) * step;
                for (size_t j = begin; j < end; j++) {
                    result[j] = string.substr(j, k);
                }
            }, i);
        }
        for (auto &thread:threads) {
            thread.join();
        }
        return result;
    }

    // Util内功能函数测试
    namespace Test {
        void test_k_mer_split() {
            std::cout << "============ Test K_mer split. =============\n";
            auto test_function = [](std::string_view string, size_t k) {
                if (k >= string.size()) { return std::vector<std::string_view>{string}; }
                size_t N = string.size() - k + 1;
                std::vector<std::string_view> result(N);
                for (size_t i = 0; i < N; i++) {
                    result[i] = string.substr(i, k);
                }
                return result;
            };
            std::string s;
            size_t k = 6;
            for (int i = 0; i < 10000; i++) {
                s += "ATCGATTCGATTCCCCTTCATCGTTACCGCGATTCACTGCGIJDCJITTCIFJCJIDFJCJIDGCCCAATICJDICJTTCJDIJFCJF";
                s += "JCIIDCJIJIQIOJOICDJIJCIJGINCOIDGUIHICIUGIHECHIODHIHAII*(&*(^%^%^$*&)(&(^%$^$*&^))(_(_(*_";
                s += "4887454548779956123230330356452663262^&%%$%^&%&^(*&*)(&*^%%$%$hJjkjcjdfjhHJFDHHCDHihoHKHK";
            }
            auto ret = split_k_mer_fast<DEFAULT_THREAD_NUMBER>(s, k);
            auto test_ret = test_function(s, k);
            bool test_pass = true;
            if (ret.size() != test_ret.size()) {
                test_pass = false;
            } else {
                for (size_t i = 0; i < ret.size(); i++) {
                    if (test_ret[i] != ret[i]) {
                        test_pass = false;
                        std::cerr << i << " " << test_ret[i] << " " << ret[i] << "\n";
                        break;
                    }
                }
            }
            if (test_pass) std::cout << "Test pass.\n"; else std::cerr << "Test fail.\n";
        }

        void test() {
            test_k_mer_split();
        }
    }
}
#endif //LSH_CPP_UTIL_H
