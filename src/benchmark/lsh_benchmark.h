//
// Created by junior on 2019/7/28.
//

#ifndef LSH_CPP_LSH_BENCHMARK_H
#define LSH_CPP_LSH_BENCHMARK_H

#include "../include/lsh_cpp.h"
#include "../include/io.h"
#include "../include/util.h"
#include "../include/time.h"
#include "../include/lsh.h"

namespace LSH_CPP::Benchmark {
//  TODO: 当前需要测试的情况有: (指标是 1. 正确率 / 召回率 ; 2. 测试变量包括下面的几个; 3. 用控制变量法; 4. 在 20-newspapers-benchmark 上测试 )
//     修改　n_permutation => 画曲线 => 测试 xx_hash 计算 min_hash 在 n_permutation = ? 时 正确率最高 (注意同时兼顾效率);
//     修改　Hash_function => xx_hash / mur_mur_hash / sha1_hash .. => 测试哪个hash在n_permutation多少下正确率最高 (同时兼顾效率);
//     修改　随机算法 => std::mt19937_64/32 , numpy-cpp-random(查一下numpy-cpp), .. => 测试随机算法对正确性影响 (不过我觉得只要是均匀随机分布应该影响不大才对..)
//     修改 seed 的值, => 测试seed不同下对 min_hash 正确率影响,不同的seed下min_hash_value_vector不同,最后也会影响LSH的结果(
//      导致false_positive和false_negative变化 => 可能出现不该查到的查到了,该查出的反而查不出结果)
//     修改 { false_positive_w , false_negative_w },测试哪种情况下准确率更好
//     => 从逻辑上当然是尽量减少false_negative召回率越高,但是candidate-set准确率会下降(可能引入较多false_positive),
//     不过考虑到后面还要对candidate-set过滤false-positive,所以重点还是减少false_negative,
//     权重测试: {0.0 1.0} {0.1 0.9}  {0.2 0.8} {0.3, ..} ... {1.0,0.0}
//  TODO: 调整测试重点 : ( 测试过程固定 : Seed = 1; weight = { 0.5 0.5 } ; 使用XXStringViewHash32[好像效果更好] )
//   =>　修改 n_permutation [128 -> 1024] 对正确率影响(兼顾效率)
//   =>  加入 sha_1/mur_mur_hash(std::hash) 测试不同hash(32/64)的正确率及效率

    using ReturnResult= std::pair<double, std::vector<HashSet<size_t >>>;

    template<size_t n_sample>
    ReturnResult test_linear_scan_performance(double threshold = 0.9) {
        std::cout << "use n_samples: " << n_sample << "\n";
        auto[train_data, test_data] = DataSpace::global_data;
        auto[train_data_set, train_labels] = train_data;
        auto[test_data_set, test_labels] = test_data;
        using MinHashType = MinHash<XXStringViewHash32, 32, n_sample>;
        using MinHashSets  =  std::vector<MinHashType>;
        MinHashSets train_minhash_sets, test_minhash_sets;
        train_minhash_sets.reserve(train_data_set.size());
        test_minhash_sets.reserve(test_data_set.size());
        for (const auto &item:train_data_set) {
            MinHashType temp;
            temp.update(item);
            train_minhash_sets.push_back(temp);
        }
        for (const auto &item: test_data_set) {
            MinHashType temp;
            temp.update(item);
            test_minhash_sets.push_back(temp);
        }
        double time = 0;
        std::vector<HashSet<size_t>> test_query_result;
        for (size_t test_index = 0; test_index < test_data_set.size(); test_index++) {
            auto test_minhash = test_minhash_sets[test_index];
            HashSet<size_t> query_result;
            TimeVar start = timeNow();
            for (size_t train_index = 0; train_index < train_data_set.size(); train_index++) {
                auto ret = minhash_jaccard_similarity(test_minhash, train_minhash_sets[train_index]);
                if (ret >= threshold) {
                    query_result.insert(train_labels[train_index]);
                }
            }
            time += millisecond_duration(timeNow() - start);
            test_query_result.push_back(query_result);
        }
        time /= test_data_set.size();
        return {time, test_query_result};
    }

    template<size_t n_sample>
    ReturnResult test_lsh_performance(double threshold = 0.9) {
        std::cout << "use n_samples: " << n_sample << "\n";
        auto[train_data, test_data] = DataSpace::global_data;
        auto[train_data_set, train_labels] = train_data;
        auto[test_data_set, test_labels] = test_data;
        using MinHashType = MinHash<XXStringViewHash32, 32, n_sample>;
        using LSH_Type = LSH<XXUInt64Hash64, size_t, 0, 0, n_sample>; // label type is size_t
        using MinHashSets  =  std::vector<MinHashType>;
        MinHashSets test_minhash_sets;
        LSH_Type lsh(threshold, {0.65, 0.35});
        test_minhash_sets.reserve(test_data_set.size());
        for (size_t i = 0; i < train_data_set.size(); i++) {
            MinHashType temp;
            temp.update(train_data_set[i]);
            lsh.insert(temp, train_labels[i]);
        }
        for (const auto &item: test_data_set) {
            MinHashType temp;
            temp.update(item);
            test_minhash_sets.push_back(temp);
        }
        double time = 0;
        std::vector<HashSet<size_t>> test_query_result;
        for (size_t test_index = 0; test_index < test_data_set.size(); test_index++) {
            auto test_minhash = test_minhash_sets[test_index];
            HashSet<size_t> query_result;
            TimeVar start = timeNow();
            query_result = lsh.query(test_minhash);
            time += millisecond_duration(timeNow() - start);
            test_query_result.push_back(query_result);
        }
        time /= test_data_set.size();
        return {time, test_query_result};
    }

    ReturnResult test_ground_truth(double threshold = 0.9) {
        auto[train_data, test_data] = DataSpace::global_data;
        auto[train_data_set, train_labels] = train_data;
        auto[test_data_set, test_labels] = test_data;
        double time = 0;
        std::vector<HashSet<size_t >> ground_truth_result;
        for (const auto &query: test_data_set) {
            HashSet<size_t> temp;
            TimeVar start = timeNow();
            for (size_t i = 0; i < train_data_set.size(); i++) {
                auto ret = jaccard_similarity(train_data_set[i], query);
                if (ret >= threshold) {
                    temp.insert(train_labels[i]);
                }
            }
            time += millisecond_duration((timeNow() - start));
            ground_truth_result.push_back(temp);
        }
        time /= test_data_set.size();
        return {time, ground_truth_result};
    }

    void lsh_benchmark() {
        namespace plt = matplotlibcpp;
        constexpr size_t population_size = 500;
        constexpr size_t train_set_size = 1000;
        constexpr size_t test_set_size = 100;
        auto ret = bootstrap_data<population_size, train_set_size, test_set_size>({10, 500}); // [10,500]
        if (!ret) {
            std::cout << "False bootstrap data!\n";
            std::exit(-1);
        }
        constexpr std::array n_samples{32, 64, 96, 128, 160, 192, 224, 256};
        std::vector<double> x, time_lsh_sets, time_linear_scan_sets,
                lsh_mean_f_scores, linear_scan_mean_f_scores;
        for (int n_sample : n_samples) { x.push_back(n_sample); }
        for_constexpr<for_bounds<0, n_samples.size()>>([&](auto index) {
            constexpr size_t n_sample = n_samples[index];
            std::cout << " Test Linear Scan: \n";
            auto[linear_scan_time, linear_scan_result] = test_linear_scan_performance<n_sample>();
            std::cout << " Test LSH: \n";
            auto[lsh_time, lsh_result] = test_lsh_performance<n_sample>();
            std::cout << " Test Ground truth: \n";
            auto[truth_time, truth_result] = test_ground_truth();
            time_linear_scan_sets.push_back(linear_scan_time);
            time_lsh_sets.push_back(lsh_time);
            double lsh_mean_f_score = 0, linear_scan_mean_f_score = 0;
            for (size_t i = 0; i < test_set_size; i++) {
                linear_scan_mean_f_score += Statistic::f_score(
                        Statistic::get_precision_recall(linear_scan_result[i], truth_result[i]));
                lsh_mean_f_score += Statistic::f_score(
                        Statistic::get_precision_recall(lsh_result[i], truth_result[i]));
            }
            linear_scan_mean_f_score /= test_set_size;
            lsh_mean_f_score /= test_set_size;
            linear_scan_mean_f_scores.push_back(linear_scan_mean_f_score);
            lsh_mean_f_scores.push_back(lsh_mean_f_score);
        });
        plt::figure_size(800, 600);
        plt::subplot(2, 1, 1);
        plt::named_plot("linear_scan time", x, time_linear_scan_sets, "b-");
        plt::named_plot("lsh time", x, time_lsh_sets, "r-");
        plt::xlabel("n_sample");
        plt::ylabel("time (ms)");
        plt::title("lsh and linear_scan performance compare \n "
                   "(Train-Set size = 1000; Query-Set size = 100; threshold = 0.9;\n"
                   "lsh_false_positive_weight:lsh_false_negative_weight = 0.65:0.35)");
        plt::legend();
        plt::subplot(2, 1, 2);
        plt::named_plot("lsh time", x, time_lsh_sets, "r-");
        plt::xlabel("n_sample");
        plt::ylabel("time (ms)");
        plt::legend();
        plt::save("../src/benchmark/lsh_performance.png");

        plt::figure_size(800, 600);
        plt::named_plot("linear_scan", x, linear_scan_mean_f_scores, "b-");
        plt::named_plot("lsh", x, lsh_mean_f_scores, "r-");
        plt::title("lsh and linear_scan f-score \n (Train-Set size = 1000; Query-Set size = 100; threshold = 0.9;\n"
                   "lsh_false_positive_weight:lsh_false_negative_weight = 0.65:0.35)");
        plt::xlabel("n_sample");
        plt::ylabel("f score");
        plt::legend();
        plt::save("../src/benchmark/lsh_f_score.png");
        clear_data();
    }
}
#endif //LSH_CPP_LSH_BENCHMARK_H
