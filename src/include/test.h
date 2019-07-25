//
// Created by junior on 19-7-23.
//

#ifndef LSH_CPP_TEST_H
#define LSH_CPP_TEST_H

#include "lsh_cpp.h"
#include "util.h"
#include "io.h"
#include "minhash.h"

namespace LSH_CPP::Test {
    using double_second = std::chrono::duration<double>;
    using double_millisecond =  std::chrono::duration<double, std::milli>;

    using TimeVar = std::chrono::high_resolution_clock::time_point;

#define duration(a) double_second(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

    using RANDOM_NUMBER_TYPE = uint64_t;
    const RANDOM_NUMBER_TYPE MAX_RANDOM_NUMBER = std::numeric_limits<RANDOM_NUMBER_TYPE>::max();
    const RANDOM_NUMBER_TYPE MIN_RANDOM_NUMBER = 1;
    const size_t MAX_RANDOM_ARRAY_SIZE = 10000000;
    std::vector<uint64_t> random_array;

    template<typename ReturnType, typename F, typename... Args>
    std::pair<ReturnType, double> compute_function_time(F func, Args &&... args) {
        TimeVar t1 = timeNow();
        ReturnType ret = func(std::forward<Args>(args)...);
        return {ret, duration(timeNow() - t1)};
    }

    template<typename F, typename ... Args>
    double compute_function_time(F func, Args &&... args) {
        TimeVar t1 = timeNow();
        func(std::forward<Args>(args)...);
        return duration(timeNow() - t1);
    }

    void init() {
        std::random_device rd;  // 将用于为随机数引擎获得种子
        std::mt19937 gen(rd()); // 以播种标准 mersenne_twister_engine
        std::uniform_int_distribution<RANDOM_NUMBER_TYPE> dis(MIN_RANDOM_NUMBER, MAX_RANDOM_NUMBER);
        random_array.reserve(MAX_RANDOM_ARRAY_SIZE);
        for (size_t i = 0; i < MAX_RANDOM_ARRAY_SIZE; i++) {
            random_array.push_back(dis(gen));
        }
    }

    template<typename Map>
    void hash_map_create_and_insert() {
        Map map;
        for (size_t i = 0; i < MAX_RANDOM_ARRAY_SIZE; i++) {
            map.insert({random_array[i], "a"});
        }
    }

    void test_k_mer_split() {
        std::cout << "============ Test K_mer split. =============\n";
        std::string s;
        size_t k = 9;
        for (int i = 0; i < 1000000; i++) {
            s += "abcdefghijklmnopqrstuvwxyz";
        }
        auto ret = compute_function_time<std::vector<std::string_view >>(split_k_mer_fast, s, k);
        printf("split_k_mer_fast: %.8f seconds\n", ret.second);
        // printf("create k_mer time %.8f seconds\n", compute_function_time([&]() { K_mer kMer(s, ret.first); }));
    }

    void test_hash_map_performance() {
        std::cout << "============ Test HashMap performance. =============\n";
        using std_hash_map = std::unordered_map<uint64_t, std::string>;
        using parallel_hash_map = phmap::flat_hash_map<uint64_t, std::string>;
        printf("std::unordered_map performance: %.8f seconds\n",
               compute_function_time(hash_map_create_and_insert < std_hash_map > ));
        printf("phmap::flat_hash_map performance: %.8f seconds\n",
               compute_function_time(hash_map_create_and_insert < parallel_hash_map > ));
    }

    void test_hash() {
        std::cout << "============ Test Hash Function =============\n";
        std::string s;
        size_t k = 1000;
        for (int i = 0; i < 1000000; i++) {
            s += "abcdefghijklmnopqrstuvwxyz";
        }
        auto string_array = split_k_mer_fast(s, k);
        printf("xx_hash performance %.8f seconds\n", compute_function_time(XXStringHash64{}, string_array));
        printf("absl hash performance %.8f seconds\n", compute_function_time(AbslStringHash64{}, string_array));
        printf("std hash performance %.8f seconds\n", compute_function_time(DefaultStringHash64{}, string_array));
        std::vector<std::string_view> array = {"hello", "hello", "you", "me", "my", "you", "please", "hello"};
        auto ret = XXStringHash64()(array);
        print_sequence_container(array);
        printf("hash result:\n");
        print_sequence_container(ret);
    }

    void test_min_hash() {
        std::vector<std::string_view> data1 = {"minhashs", "are", "a__", "probabilistic", "data", "structure", "for",
                                               "estimating", "the", "similarity", "between", "datasets"};

        std::vector<std::string_view> data2 = {"minhash", "is", "a_", "probability", "data", "structure", "for",
                                               "estimating", "the", "similarity", "between", "documents"};
        std::string s;
        size_t k = 1000;
        for (int i = 0; i < 100000; i++) {
            s += "abcdefghijklmnopqrstuvwxyz";
        }
        auto string_array = split_k_mer_fast(s, k);
        TimeVar start = timeNow();
        MinHash hash1(XXStringHash64{});
        MinHash hash2(XXStringHash64{});
        hash1.update(data1);
        hash2.update(data2);
        double time = duration(timeNow() - start);
        start = timeNow();
        auto ret = jaccard_similarity(hash1, hash2);
        double time_2 = duration(timeNow() - start);
        printf("similarity: %.8f ; update-time : %.8f seconds ; similarity-time %.8f \n", ret, time, time_2);
    }

    namespace xs = xsimd;
    template<typename T>
    using vector_type = std::vector<T, xsimd::aligned_allocator<T, XSIMD_DEFAULT_ALIGNMENT>>;

    template<typename T>
    void mean(const vector_type<T> &a, const vector_type<T> &b, vector_type<T> &res) {
        std::size_t size = a.size();                                  // T类型大小 size
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;  // T类型在simd的大小,一般simd_size是size的几倍
        std::size_t vec_size = size - size % simd_size;  // 计算for循环每次加一个simd_size后,最大的上界.当继续加size_size超过这个上界时就结束循环
        for (std::size_t i = 0; i < vec_size; i += simd_size) {
            // 每一个循环步骤在和原来循环几乎同样的指令周期内计算simd_size大小的数据量.
            // 相当于一次循环展开而且还没有多余的指令开销.
            auto ba = xs::load_aligned(&a[i]);
            auto bb = xs::load_aligned(&b[i]);
            auto bres = (ba + bb) / 2;
            bres.(&res[i]);
        }
        for (std::size_t i = vec_size; i < size; ++i) {  // 将res中剩下的部分计算完(无法放在一个simd_size直接计算),就按照传统的一步步循环计算即可.
            res[i] = (a[i] + b[i]) / 2;
        }
    }


    void mean_(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &res) {
        std::size_t size = res.size();
        for (std::size_t i = 0; i < size; ++i) {
            res[i] = (a[i] + b[i]) / 2;
        }
    }

    void test_simd() {
        // 使用带有侵入性的 xsimd::batch,而且只能使用avx指令(文档好像是这样写)
//        xsimd::batch<double, 4> a(1.5, 2.5, 3.5, 4.5);
//        xsimd::batch<double, 4> b(2.5, 3.5, 4.5, 5.5);
//        auto batch_ret = (a + b) / 2;  // 返回的是 xsimd::batch 类型,可以直接std::cout
//        std::cout << batch_ret << std::endl;

        // 使用不具有侵入性的方法,依然可用stl的容器,至于simd指令由xsimd库选择最优的执行.
        using vector_double = vector_type<double>;
        using origin_vector = std::vector<double>;
        origin_vector o_a, o_b, o_ret;
        vector_double v_a, v_b, v_ret;
        size_t N = 10000000;
        v_ret.resize(N); // 注意必须提前 resize vector_ret
        o_ret.resize(N);
        for (size_t i = 0; i < N; i++) {
            o_a.push_back(i);
            o_b.push_back(i + 1);
            v_a.push_back(i);
            v_b.push_back(i + 1);
        }
        auto time2 = compute_function_time(mean_, o_a, o_b, o_ret);
        auto time1 = compute_function_time(mean<double>, v_a, v_b, v_ret);
        std::cout << time1 << " " << time2 << "\n";
    }

    void test() {
        //init();
        //test_hash_map_performance();
        //test_k_mer_split();
        //test_hash();
        //test_min_hash();
        test_simd();
    }
}
#endif //LSH_CPP_TEST_H
