//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_UTIL_H
#define LSH_CPP_UTIL_H

#include "lsh_cpp.h"

namespace LSH_CPP {
    struct K_mer {
        bool operator==(const K_mer &kMer) const {
            return value == kMer.value;
        }

        std::string_view value;
        mutable std::vector<size_t> pos_list;
        // TODO: pos_list 后面可以考虑用BitSet压缩,不然内存开销会很大.
        //  方法: 先使用相对位置编码,保留pos_list[0],剩下的保存相邻元素的差值,比如原pos_list是{ 0,10,28,35,45 }
        //  相对位置编码为 { 0,10-0,28-10,35-28,45-35 } = { 0,10,18,7,10 } .
        //  然后再用unary_bit_set对上面的序列做位编码(见<<信息检索>>相关压缩方法).

        // 通过weight()方法得到权重(重复元素个数).
        size_t weight() const {
            return pos_list.size();
        }
    };
} // 前置声明
namespace std {
    // inject specialization of std::hash for K_mer
    template<>
    struct hash<LSH_CPP::K_mer> {
        std::size_t operator()(LSH_CPP::K_mer const &k_mer) const {
            return phmap::HashState::combine(0, k_mer.value);
        }
    };
}
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
    HashSet<K_mer> split_k_mer_fast(const std::string_view &string, size_t k) {
        if (k >= string.size()) {
            return {{string, {0}}};
        }
        size_t N = string.size() - k + 1;
        HashSet<K_mer> result;
        for (size_t i = 0; i < N; i++) {
            auto substr = string.substr(i, k);
            // 这里pos_list一开始留空,等后面插入(不管成功还是失败)后,从返回值的first得到迭代器,然后再修改pos_list.
            // 本来 unordered_set.insert(Key) 返回的应该是 { const_iterator, bool },
            // 但因为我们把 K_mer 结构体的 pos_list 成员改为 mutable(可变). 因此,即使是const_iterator,也可以修改它.
            // 这样写可以最大效率实现插入和修改pos_list,完全不需要插入一次,再查找一次,因为插入的时候本身就是在查找.
            // 另外,这种写法完全可以用 map [key(不可变部分)] = value (可变部分) 来代替.
            (*(result.insert(K_mer{substr, {}})).first).pos_list.push_back(i);
        }
        return result;
    }

//    废弃担心std::string生命周期影响std::string_view而引入的K_mer
//    struct K_mer {
//        const std::string origin_string;
//        std::vector<std::string_view> sub_strings;
//
//        K_mer(std::string origin_string, std::vector<std::string_view> sub_strings) :
//                origin_string(std::move(origin_string)), sub_strings(std::move(sub_strings)) {}
//    };

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
}
#endif //LSH_CPP_UTIL_H
