//
// Created by junior on 19-7-23.
//

#ifndef LSH_CPP_TEST_H
#define LSH_CPP_TEST_H

#include "../include/lsh_cpp.h"
#include "../include/util.h"
#include "../include/io.h"
#include "../include/time_def.h"
#include "../include/k_shingles.h"
#include "../include/hash.h"
#include "../include/minhash.h"
#include "../include/lsh.h"
#include "../include/weight_minhash.h"
#include "../include/lru_cache.h"

namespace LSH_CPP::Test {
    using RANDOM_NUMBER_TYPE = uint64_t;
    const RANDOM_NUMBER_TYPE MAX_RANDOM_NUMBER = std::numeric_limits<RANDOM_NUMBER_TYPE>::max();
    const RANDOM_NUMBER_TYPE MIN_RANDOM_NUMBER = 1;
    const size_t MAX_RANDOM_ARRAY_SIZE = 10000000;
    std::vector<uint64_t> random_array;

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
        std::string s, temp = "ATCGCCTACTGCTACCCTAATCGGCTAATTCTTTGCTAGCTT";
        size_t k = 9;
        std::uniform_int_distribution<size_t> dis(1, temp.size());
        std::mt19937_64 generator(std::random_device{}());
        for (int i = 0; i < 1000000; i++) {
            std::sample(temp.begin(), temp.end(), std::back_inserter(s), dis(generator),
                        std::mt19937_64{std::random_device{}()});
        }
        std::cout << "string size : " << s.size() << "\n";
        std::cout << "k parameter : " << k << "\n";
        auto ret = compute_function_time<HashSet<K_shingling>>(split_k_shingling_fast, s, k);
        printf("split_k_mer_fast: %.8f seconds\n", ret.second);
        std::cout << "k_mer set size (not repeat element) : " << ret.first.size() << "\n";
        size_t count = 0;
        for (const auto &item:ret.first) { count += item.weight(); }
        std::cout << "k_mer set size (include repeat element) : " << count << "\n";
    }

    void test_hash_map_performance() {
        std::cout << "============ Test HashMap performance. =============\n";
        using std_hash_map = std::unordered_map<uint64_t, std::string>;
        using parallel_hash_map = phmap::flat_hash_map<uint64_t, std::string>;
        printf("std::unordered_map performance: %.8f seconds\n",
               compute_function_time(hash_map_create_and_insert<std_hash_map>));
        printf("phmap::flat_hash_map performance: %.8f seconds\n",
               compute_function_time(hash_map_create_and_insert<parallel_hash_map>));
    }

    void test_hash() {
        std::cout << "============ Test Hash Function =============\n";
        std::vector<std::string_view> string_view_array =
                {"hello", "hello", "you", "me", "my", "you", "please", "hello"};
        std::vector<std::string> string_array =
                {"hello", "hello", "you", "me", "my", "you", "please", "hello"};
        std::cout << "hello : " << XXStringViewHash64{}("hello") << "\n";
        print_sequence_container(element_wise_hash(XXStringViewHash64{}, string_view_array));
        print_sequence_container(element_wise_hash(XXStringHash64{}, string_array));
    }

    void test_min_hash() {
        HashSet<std::string_view> data1 =
                {"minhash", "is", "a", "probabilistic", "data", "structure", "for",
                 "estimating", "the", "similarity", "between", "datasets"};

        HashSet<std::string_view> data2 =
                {"minhash", "is", "a", "probability", "data", "structure", "for",
                 "estimating", "the", "similarity", "between", "documents"};
        MinHash hash1(XXStringViewHash64{});
        MinHash hash2(XXStringViewHash64{});
        TimeVar start = timeNow();
        hash1.update(data1);
        hash2.update(data2);
        double update_time = second_duration(timeNow() - start);
        auto minhash_sm = minhash_jaccard_similarity(hash1, hash2);
        auto sm = jaccard_similarity(data1, data2);
        printf("minhash similarity: %.8f \nactual similarity: %.8f \nupdate-time : %.8f seconds \n",
               minhash_sm, sm, update_time);
    }

    /**
     * 测试得到的相对较好的参数设置:
     * MinHash : Seed = 1; 使用 XXStringViewHash32{} ; n_permutation = 128 (n_permutation可以考虑用 benchmark 测试)
     * LSH : weights = { 0.5 , 0.5 }
     */
    void test_lsh_minhash() {
        HashSet<std::string_view> data_1 = {"minhash", "is", "a", "probabilistic", "data", "structure",
                                            "for", "estimating", "the", "similarity", "between",
                                            "datasets"};
        HashSet<std::string_view> data_2 = {"minhash", "is", "a", "probability", "data", "structure",
                                            "for", "estimating", "the", "similarity", "between",
                                            "documents"};
        HashSet<std::string_view> data_3 = {"minhash", "is", "probability", "data", "structure", "for",
                                            "estimating", "the", "similarity", "between", "documents"};
        MinHash m1, m2, m3; // use default XXStringViewHash32 hash function
        m1.update(data_1);
        m2.update(data_2);
        m3.update(data_3);
        std::cout << "(m1,m2) jaccard similarity: " << minhash_jaccard_similarity(m1, m2) << "\n";
        // print_minhash_table(m1, m2);
        std::cout << "(m1,m3) jaccard similarity: " << minhash_jaccard_similarity(m1, m3) << "\n";
        // print_minhash_table(m1, m3);
        double threshold = 0.7;
        LSH lsh(threshold, {0.5, 0.5});
        lsh.print_config();
        lsh.insert(m2, "m2");
        lsh.insert(m3, "m3");
        auto ret = lsh.query(m1);
        printf("query result (jaccard similarity > %.2f compare to m1) : ", threshold);
        for (const auto &item:ret) {
            std::cout << item << " ";
        }
        std::cout << "\n";
    }

    void test_make_constexpr_array() {
        constexpr auto array = make_constexpr_array(make_sequence<10>([](size_t index) { return index * 2 + 10; }));

        // run-time iteration
        std::for_each(array.begin(), array.end(), [](const auto &item) { std::cout << item << " "; });
        std::cout << std::endl;

        // compile-time iteration (use for constexpr)
        for_constexpr<for_bounds<0, array.size()>>([&](auto index) {
            std::cout << "index: " << index << " value: " << array[index] << "\n";
        });
    }

    // for_constexpr example: create MinHash by constexpr n_permutation array in compile time.
    void test_for_constexpr() {
        std::vector<std::string_view> data;
        for (size_t i = 0; i < 10; i++) { data.push_back(std::to_string(i)); }
        constexpr std::array n_permutations{50, 60, 70, 80, 90, 100, 110, 128};
        for_constexpr<for_bounds<0, n_permutations.size()>>([&](auto index) {
            constexpr size_t n_permutation = n_permutations[index];
            MinHash<XXStringViewHash32, 32, n_permutation> minHash(XXStringViewHash32{});
            minHash.update(data);
            print_sequence_container(minHash.hash_values);
        });
    }

    void test_weight_minhash() {
        std::vector<size_t> data1 = {1, 3, 4, 5, 6, 7, 8, 9, 10, 4};
        std::vector<size_t> data2 = {2, 4, 3, 8, 4, 7, 10, 9, 0, 0};
        WeightMinHash<10, size_t> A, B;
        A.update(data1);
        B.update(data2);
        std::cout << weight_minhash_jaccard(A, B) << "\n";
        std::cout << generalized_jaccard_similarity(data1, data2) << "\n";
    }

    struct A {
        using WeightType = size_t;
        using ValueType = std::string_view;
        ValueType _value;
        WeightType _weight;

        bool operator==(const A &a) const {
            return _value == a._value;
        }

        [[nodiscard]] inline ValueType value() const {
            return _value;
        }

        [[nodiscard]] inline WeightType weight() const {
            return _weight;
        }
    };

    void test_weight_minhash_by_set() {
        HashSet<A> a_set = {
                {"a", 3},
                {"b", 2},
                {"c", 1},
                {"e", 4},
                {"k", 3},
                {"p", 5},
                {"g", 4}
        };
        HashSet<A> b_set = {
                {"a", 2},
                {"b", 3},
                {"d", 1},
                {"k", 5},
                {"f", 3},
                {"p", 2},
                {"m", 9}
        };
        ///////////////////////////a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p
        std::vector<size_t> A_v = {3, 2, 1, 0, 4, 0, 4, 0, 0, 0, 3, 0, 0, 0, 0, 5, 0, 0, 0, 0};
        std::vector<size_t> B_v = {2, 3, 0, 1, 0, 3, 0, 0, 0, 0, 5, 0, 9, 0, 0, 2, 0, 0, 0, 0};
        std::cout << "actual jaccard similarity is : " << generalized_jaccard_similarity(a_set, b_set) << "\n";
        WeightMinHash<200000, A, 128> min_hash_a, min_hash_b;
        min_hash_a.update(a_set);
        min_hash_b.update(b_set);
        std::cout << "weight minhash jaccard similarity is : " << weight_minhash_jaccard(min_hash_a, min_hash_b);
//        WeightMinHash<20, size_t> min_hash_a, min_hash_b;
//        min_hash_a.update(A_v);
//        min_hash_b.update(B_v);
//        std::cout << "weight minhash jaccard similarity is : " << weight_minhash_jaccard(min_hash_a, min_hash_b);
    }

    void test_lru_cache() {
        // 测试 lru_cache 在 Eigen Align 机制下能否工作
        using Cache =  lru_cache<int, Eigen::Vector4i, phmap::Hash<int>, phmap::EqualTo<int>,
                Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4i> >,
                phmap::flat_hash_map>;
        Cache cache(10);
        cache.put(0, Eigen::Vector4i{0, 1, 2, 3});
        cache.put(1, Eigen::Vector4i{2, 4, 6, 8});
        std::cout << std::boolalpha << cache.contains(0) << " " << cache.contains(1) << "\n";
        if (auto ret = cache.get(2); ret.has_value()) {
            std::cout << *ret << "\n";
        } else {
            cache.put(2, Eigen::Vector4i{5, 10, 15, 20});
        }
        if (auto ret = cache.get(2); ret.has_value()) {
            std::cout << *ret << "\n";
        } else {
            cache.put(2, Eigen::Vector4i{5, 10, 15, 20});
        }
    }

    void test_dna_shingling() {
        std::string_view dna_1 = "ATCGTATCGTATCGT", dna_2 = "ATCGTTTACGTATCGTATCG";
        auto data1 = split_dna_shingling<5, no_weight>(dna_1);
        auto data2 = split_dna_shingling<5, no_weight>(dna_2);
        StdDNAShinglingHash64<5> hash64;
        auto func = [&](const auto &item) {
            std::cout << item.value().to_string()
                      << " " << dna_shingling_decode<5, no_weight>(item)
                      // << " "<< item.weight()
                      << " " << hash64(item) << "\n";
        };
        std::for_each(data1.begin(), data1.end(), std::ref(func));
        std::cout << std::endl;
        std::for_each(data2.begin(), data2.end(), std::ref(func));
        std::cout << "sim:" << jaccard_similarity(data1, data2) << "\n";
    }

    void test_parallel_get_mean() {
        std::vector<double> v(10'000'007, 0.5);
        TimeVar start = timeNow();
        auto ret = Statistic::get_mean(v);
        std::cout << "ret:" << ret << " time:" << millisecond_duration((timeNow() - start)) << " ms\n";
    }

    void test() {
        //init();
        //test_hash_map_performance();
        //test_k_mer_split();
        //test_hash();
        //test_min_hash();
        //test_lsh_minhash();
        //test_make_constexpr_array();
        //test_for_constexpr();
        //test_weight_minhash();
        //test_weight_minhash_by_set();
        //test_lru_cache();
        //test_dna_shingling();
        test_parallel_get_mean();
    }
}
namespace std {
    template<>
    struct hash<LSH_CPP::Test::A> {
        std::size_t operator()(LSH_CPP::Test::A const &a) const {
            return phmap::HashState::combine(0, a._value);
        }
    };
}
#endif //LSH_CPP_TEST_H
