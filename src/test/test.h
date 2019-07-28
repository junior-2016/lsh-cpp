//
// Created by junior on 19-7-23.
//

#ifndef LSH_CPP_TEST_H
#define LSH_CPP_TEST_H

#include "../include/lsh_cpp.h"
#include "../include/util.h"
#include "../include/io.h"
#include "../include/hash.h"
#include "../include/minhash.h"
#include "../include/lsh.h"

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
        auto ret = compute_function_time<std::vector<std::string_view >>
                (split_k_mer_fast, s, k);
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
//        std::string s;
//        size_t k = 1000;
//        for (int i = 0; i < 1000000; i++) {
//            s += "abcdefghijklmnopqrstuvwxyz";
//        }
//        auto string_array = split_k_mer_fast(s, k);
        std::vector<std::string_view> string_view_array =
                {"hello", "hello", "you", "me", "my", "you", "please", "hello"};
        std::vector<std::string> string_array =
                {"hello", "hello", "you", "me", "my", "you", "please", "hello"};
        std::cout << "hello : " << XXStringViewHash64{}("hello") << "\n";
        print_sequence_container(element_wise_hash(XXStringViewHash64{}, string_view_array));
        print_sequence_container(element_wise_hash(XXStringHash64{}, string_array));
    }

    void test_min_hash() {
        std::vector<std::string_view> data1 = {"minhash", "is", "a", "probabilistic", "data", "structure", "for",
                                               "estimating", "the", "similarity", "between", "datasets"};

        std::vector<std::string_view> data2 = {"minhash", "is", "a", "probability", "data", "structure", "for",
                                               "estimating", "the", "similarity", "between", "documents"};
        std::string s;
        size_t k = 12;
        for (int i = 0; i < 100000; i++) {
            s += "abcdefghijklmnopqrstuvwxyz";
        }
        auto string_array = split_k_mer_fast(s, k);
        MinHash hash1(XXStringViewHash64{});
        MinHash hash2(XXStringViewHash64{});
        TimeVar start = timeNow();
        hash1.update(string_array);
        hash2.update(string_array);
        double time = duration(timeNow() - start);
        start = timeNow();
        auto ret = estimated_jaccard_similarity(hash1, hash2);
        double time_2 = duration(timeNow() - start);
        printf("similarity: %.8f ; update-time : %.8f seconds ; similarity-time %.8f \n", ret, time, time_2);
    }

    void test_lsh_minhash() {
        std::vector<std::string_view> data_1 = {"minhash", "is", "a", "probabilistic", "data", "structure", "for",
                                                "estimating", "the", "similarity", "between", "datasets"};
        std::vector<std::string_view> data_2 = {"minhash", "is", "a", "probability", "data", "structure", "for",
                                                "estimating", "the", "similarity", "between", "documents"};
        std::vector<std::string_view> data_3 = {"minhash", "is", "probability", "data", "structure", "for",
                                                "estimating", "the", "similarity", "between", "documents"};
        MinHash m1(XXStringViewHash32{}), m2(XXStringViewHash32{}), m3(XXStringViewHash32{});
        m1.update(data_1);
        m2.update(data_2);
        m3.update(data_3);
        std::cout << "(m1,m2) jaccard similarity: " << estimated_jaccard_similarity(m1, m2) << "\n";
        // print_minhash_table(m1, m2);
        std::cout << "(m1,m3) jaccard similarity: " << estimated_jaccard_similarity(m1, m3) << "\n";
        // print_minhash_table(m1, m3);
        double threshold = 0.65;
        LSH lsh(threshold, {0.5, 0.5});
        // 尽量以减小false negative为主.(存在false positive也无所谓,因为后面会对candidate set进一步过滤false positive)
        lsh.print_config();
        lsh.insert(m2, "m2");
        lsh.insert(m3, "m3");
        auto ret = lsh.query(m1);
        printf("query result (jaccard similarity > %.2f compare to m1) : ", threshold);
        for (const auto &item:ret) {
            std::cout << item << " ";
        }
    }

    void test_generator_sequence() {
        constexpr auto array = make_constexpr_array(make_sequence<128>(step_rule<1>{}));
        std::cout << array[0] << " " << array[127] << "\n";
    }

    void test() {
        //init();
        //test_hash_map_performance();
        //test_k_mer_split();
        //test_hash();
        //test_min_hash();
        //test_generator_sequence();
        test_lsh_minhash();
    }
}
#endif //LSH_CPP_TEST_H
