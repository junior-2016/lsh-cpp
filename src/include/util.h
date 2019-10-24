//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_UTIL_H
#define LSH_CPP_UTIL_H

#include "lsh_cpp.h"

namespace LSH_CPP {
    template<typename T>
    using HashSet = phmap::flat_hash_set<T>;
    template<typename K, typename V>
    using HashMap  = phmap::flat_hash_map<K, V>;

    /**
     * Non-capture lambda can be transfer to function pointer directly,
     * so the first argument can be a non-capture function with <double (*)(double, void *)> signature.
     * @tparam Func lambda expression type : [](double,void*)->double{...}
     * @param f : integration function for gsl_function
     * @param range : integration interval from range.first to range.second
     * @param params : parameters array for gsl_function
     * GSL quad数值积分计算接口:
     * gsl_integration_cquad(const gsl_function * f, double a, double b, double epsabs, double epsrel,
     * gsl_integration_cquad_workspace * workspace, double * result, double * abserr, size_t * nevals)
     * 最后两个参数可以设置成null(不需要计算误差的时候)
     * 这里主要有几个参数可以调整: work_space_alloc_n , epsrel , epsabs,但主要是前面两个.
     * 我初始的设置是 work_space_alloc_n = 1000 && epsrel = 1e-7 && epsabs = 0,虽然精度不错但耗时太长,应该减少一定精度.
     * 最后设置为 work_space_alloc_n = 500 && epsrel = 1e-4, 积分耗时减少.
     * (积分计算耗时其实问题不大,因为LSH只在初始化的时候会优化参数,积分计算最多执行一次)
     */
    template<typename Func>
    double numerical_integration(Func &&f, std::pair<double, double> range, double *params) {
        gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(500);
        double result;
        gsl_function F;
        F.function = static_cast<double (*)(double, void *)> (f);
        F.params = reinterpret_cast<void *>(params);
        gsl_integration_cquad(&F, range.first, range.second, 0, 1e-4, w,
                              &result, nullptr, nullptr);
        gsl_integration_cquad_workspace_free(w);
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

    // 编译期生成固定序列
    // reference: https://stackoverflow.com/questions/45940284/array-initialisation-compile-time-constexpr-sequence
    namespace detail {
        template<size_t ... Is>
        constexpr auto make_constexpr_array_from_sequence_impl(std::index_sequence<Is...>) {
            return std::array{Is...}; // C++ 17 simplify
        }

        template<size_t ...Is, typename generator_rule>
        constexpr auto make_sequence_impl(std::index_sequence<Is...>, generator_rule rule) {
            return std::index_sequence<rule(Is)...>{};
        }
    }

    template<size_t N, typename generator_rule>
    constexpr auto make_sequence(generator_rule rule) {
        return detail::make_sequence_impl(std::make_index_sequence<N>{}, rule);
    }

    // 使用样例见test.h
    template<typename Sequence>
    constexpr auto make_constexpr_array(Sequence sequence) {
        return detail::make_constexpr_array_from_sequence_impl(sequence);
    }

    // 不用 std::array,直接使用 constexpr std::integer_sequence<Type,num1,num2,...> sequence;
    // 然后通过 constexpr Type val = get_sequence_value_by_index(sequence,index); 得到对应位置的值
    // 不过这里使用时 index 必须是 const expression, 不能是变量. 比如下面的:
    // for (size_t index = 0; index < sequence.size(); index++) { constexpr Type val = get_sequence_value_by_index(sequence,index); }
    // 就不能通过编译. 因为index是变量,虽然它的变化范围相当明确,但它语义上不属于constexpr. 要执行这种for-loop的constexpr只能自己实现.
    template<typename T, T ... Numbers>
    constexpr T get_sequence_value_by_index(std::integer_sequence<T, Numbers ...>, size_t index) {
        constexpr T temp[] = {Numbers ...};
        return temp[index];
    }

    // for-loop constexpr
    // reference: https://nilsdeppe.com/posts/for-constexpr
    // test-case 见 test.h
    template<size_t Lower, size_t Upper>
    struct for_bounds {
        static constexpr const size_t lower = Lower;
        static constexpr const size_t upper = Upper;
    };
    namespace detail {
        template<size_t lower, size_t... Is, class F>
        void for_constexpr_impl(F &&f,
                                std::index_sequence<Is...> /*meta*/) {
            (void) std::initializer_list<char>{
                    ((void) f(std::integral_constant<size_t, Is + lower>{}),
                            '0')...};
        }
    }  // namespace for_constexpr_impl
    template<class Bounds0, class F>
    void for_constexpr(F &&f) {
        detail::for_constexpr_impl<Bounds0::lower>(
                std::forward<F>(f),
                std::make_index_sequence<Bounds0::upper - Bounds0::lower>{});
    }

    // static fix-size eigen array declaration:
    template<typename _Scalar, int row>
    static const Eigen::Array<_Scalar, row, 1> one_eigen_array = Eigen::Array<_Scalar, row, 1>::Constant(1);
    template<typename _Scalar, int row>
    static const Eigen::Array<_Scalar, row, 1> zero_eigen_array = Eigen::Array<_Scalar, row, 1>::Constant(0);
    template<typename _Scalar, int row>
    static const Eigen::Array<_Scalar, row, 1> min_eigen_array =
            Eigen::Array<_Scalar, row, 1>::Constant(std::numeric_limits<_Scalar>::min());
    template<typename _Scalar, int row>
    static const Eigen::Array<_Scalar, row, 1> max_eigen_array =
            Eigen::Array<_Scalar, row, 1>::Constant(std::numeric_limits<_Scalar>::max());

    // constexpr power function(only for integer power)
    constexpr int32_t pow(int32_t x, uint32_t y) {
        return y == 0 ? 1 : x * pow(x, y - 1);
    }

    // 统计相关数据量需要的函数
    namespace Statistic {
        using precision_type = double;
        using recall_type = double;
        using precision_recall_t = std::pair<precision_type, recall_type>;

        template<typename T>
        precision_recall_t get_precision_recall(const HashSet<T> &found, const HashSet<T> &truth) {
            double intersection = 0;
            precision_type precision;
            recall_type recall;
            for (const auto &f :found) {
                if (truth.find(f) != truth.end()) {
                    intersection++;
                }
            }
            if (found.size() == 0) {
                precision = 0;
            } else {
                precision = intersection / found.size();
            }
            if (truth.size() == 0) {
                recall = 1.0;
            } else {
                recall = intersection / truth.size();
            }
            if ((found.size() == 0) && (truth.size() == 0)) {
                precision = 1.0;
                recall = 1.0;
            }
            return {precision, recall};
        }

        /**
         * F1 score = (2 * precision * recall) / (precision + recall)
         * @param precision 准确率
         * @param recall    召回率
         */
        double f_score(precision_type precision, recall_type recall) {
            if ((precision == 0) && (recall == 0)) {
                return 0;
            }
            return (2.0 * precision * recall) / (precision + recall);
        }

        double f_score(const precision_recall_t &p_r) {
            auto[p, r] = p_r;
            return f_score(p, r);
        }

        /**
         * 计算一个序列的百分之p分位数.使用 streaming percentiles library.
         * reference:
         * [1] https://en.wikipedia.org/wiki/Percentile
         * [2] https://github.com/sengelha/streaming-percentiles
         */
        double get_percentile(const std::vector<double> &sequence, double p, double epsilon = 0.1) {
            assert(0.0 <= p && p <= 1.0);
            stmpct::gk<double> g(epsilon);
            std::for_each(sequence.begin(), sequence.end(), [&](const auto &d) { g.insert(d); });
            return g.quantile(p);
        }

        // 计算一个序列的算术均值
        double get_mean(const std::vector<double> &sequence) {
#ifdef USE_CXX_PARALLISM_TS
            double sum = std::reduce(std::execution::par, sequence.begin(), sequence.end()); // C++17并行求和
#else
            double sum = std::reduce(sequence.begin(),sequence.end());
#endif
            return sum / (double) sequence.size();
        }
    }
}
#endif //LSH_CPP_UTIL_H
