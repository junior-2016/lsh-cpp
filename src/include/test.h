//
// Created by junior on 19-7-23.
//

#ifndef LSH_CPP_TEST_H
#define LSH_CPP_TEST_H

#include "lsh_cpp.h"
#include "util.h"


namespace LSH_CPP::Test {
    using double_second = std::chrono::duration<double>;
    using double_millisecond =  std::chrono::duration<double, std::milli>;

    using TimeVar = std::chrono::high_resolution_clock::time_point;

#define duration(a) double_second(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

    using RANDOM_NUMBER_TYPE = uint64_t;
    const RANDOM_NUMBER_TYPE MAX_RANDOM_NUMBER = std::numeric_limits<RANDOM_NUMBER_TYPE>::max();
    const RANDOM_NUMBER_TYPE MIN_RANDOM_NUMBER = 1;
    const size_t MAX_RANDOM_ARRAY_SIZE = 1000000;
    static std::vector<uint64_t> random_array(MAX_RANDOM_ARRAY_SIZE);

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
        for (size_t i = 0; i < MAX_RANDOM_ARRAY_SIZE; i++) {
            random_array[i] = dis(gen);
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
        size_t k = 9;
        for (int i = 0; i < 1000000; i++) {
            s += "abcdefghijklmnopqrstuvwxyz";
        }
        auto ret = compute_function_time<std::vector<std::string_view >>(split_k_mer_fast<DEFAULT_THREAD_NUMBER>, s, k);
        auto test_ret = compute_function_time<std::vector<std::string_view>>(test_function, s, k);
        printf("split_k_mer_fast time:%.8f second\n", ret.second);
        printf("split_k_mer_single_thread:%.8f second\n", test_ret.second);
        bool test_pass = true;
        if (ret.first.size() != test_ret.first.size()) {
            test_pass = false;
        } else {
            for (size_t i = 0; i < ret.first.size(); i++) {
                if (test_ret.first[i] != ret.first[i]) {
                    test_pass = false;
                    std::cerr << i << " " << test_ret.first[i] << " " << ret.first[i] << "\n";
                    break;
                }
            }
        }
        if (test_pass) std::cout << "Test pass.\n"; else std::cerr << "Test fail.\n";
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

    void test() {
        init();
        test_hash_map_performance();
        test_k_mer_split();
    }
}
#endif //LSH_CPP_TEST_H
