//
// Created by junior on 2019/7/28.
//

#ifndef LSH_CPP_LSH_BENCHMARK_H
#define LSH_CPP_LSH_BENCHMARK_H

#include "../include/lsh_cpp.h"

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
    void lsh_benchmark() {
        constexpr size_t population_size = 500;
        constexpr size_t train_set_size = 1000;
        constexpr size_t test_set_size = 100;
        auto ret = bootstrap_data<population_size, train_set_size, test_set_size>({10, 500}); // [10,500]
        if (!ret) {
            std::cout << "False bootstrap data!\n";
            std::exit(-1);
        }

        clear_data();
    }
}
#endif //LSH_CPP_LSH_BENCHMARK_H
