//
// Created by junior on 2019/8/2.
//

#ifndef LSH_CPP_WEIGHT_MINHASH_BENCHMARK_H
#define LSH_CPP_WEIGHT_MINHASH_BENCHMARK_H

#include "../include/lsh_cpp.h"
#include "../include/weight_minhash.h"
#include "../include/time.h"

namespace LSH_CPP::Benchmark {
    namespace detail {
        constexpr auto n_samples = make_constexpr_array(make_sequence<16>([](size_t index) {
            return (index + 1) * 10;
        }));// 10,20,...,160.
        constexpr size_t dim = 5000;
        constexpr size_t repeat_number = 100;

        std::vector<double> performance() {
            std::random_device randomDevice;
            std::mt19937_64 generator(randomDevice());
            std::uniform_real_distribution<float> distribution(0, dim);
            std::vector<std::vector<float>> data;
            data.reserve(repeat_number);
            for (size_t i = 0; i < repeat_number; i++) {
                std::vector<float> temp;
                temp.reserve(dim);
                for (size_t j = 0; j < dim; j++) {
                    temp.push_back(distribution(generator));
                }
                data.push_back(temp);
            }
            std::vector<double> times;
            for_constexpr<for_bounds<0, n_samples.size()>>([&](auto index) {
                constexpr auto sample = n_samples[index];
                double time = 0;
                WeightMinHash<dim, sample> weight_minhash;
                for (size_t i = 0; i < repeat_number; i++) {
                    TimeVar start = timeNow();
                    weight_minhash.update(data[i]);
                    time += millisecond_duration(timeNow() - start);
                }
                time = time / repeat_number;
                times.push_back(time);
                printf("n_sample : %ld , time : %.5f ms\n", sample, time);
            });
            return times;
        }

        std::vector<double> accurate() {
            std::random_device randomDevice;
            std::mt19937_64 generator(randomDevice());
            std::uniform_real_distribution<float> distribution(0, dim);
            std::vector<std::vector<float>> data1, data2;
            data1.reserve(repeat_number);
            data2.reserve(repeat_number);
            for (size_t i = 0; i < repeat_number; i++) {
                std::vector<float> temp1, temp2;
                temp1.reserve(dim);
                temp2.reserve(dim);
                for (size_t j = 0; j < dim; j++) {
                    temp1.push_back(distribution(generator));
                    temp2.push_back(distribution(generator));
                }
                data1.push_back(temp1);
                data2.push_back(temp2);
            }
            std::vector<double> errors;
            for_constexpr<for_bounds<0, n_samples.size()>>([&](auto index) {
                constexpr auto sample = n_samples[index];
                double error = 0;
                WeightMinHash<dim, sample> weight_minhash1, weight_minhash2;
                for (size_t i = 0; i < repeat_number; i++) {
                    weight_minhash1.update(data1[i]);
                    weight_minhash2.update(data2[i]);
                    auto estimated_jaccard = weight_minhash_jaccard(weight_minhash1, weight_minhash2);
                    auto actual_jaccard = generalized_jaccard_similarity(data1[i], data2[i]);
                    error += std::abs(estimated_jaccard - actual_jaccard);
                }
                error = error / repeat_number;
                errors.push_back(error);
                printf("n_sample : %ld , mean_error : %.5f \n", sample, error);
            });
            return errors;
        }
    }

    void weight_minhash_benchmark() {
        namespace plt = matplotlibcpp;
        auto performance_time = detail::performance();
        auto accurate_time = detail::accurate();
        std::vector<double> x;
        for (size_t n_sample : detail::n_samples) {
            x.push_back(n_sample);
        }
        plt::subplot(2, 1, 1);
        plt::plot(x, performance_time, "r-");
        plt::title("weight minhash performance benchmark");
        plt::xlabel("n_sample");
        plt::ylabel("time(ms)");
        plt::grid(true);

        plt::subplot(2, 1, 2);
        plt::plot(x, accurate_time, "k-");
        plt::title("weight minhash accurate benchmark");
        plt::xlabel("n_sample");
        plt::ylabel("abs mean error");
        plt::grid(true);

        plt::tight_layout();
        plt::save("../src/benchmark/weight_minhash_benchmark.png");
    }
}
#endif //LSH_CPP_WEIGHT_MINHASH_BENCHMARK_H
