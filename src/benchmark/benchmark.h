//
// Created by junior on 2019/7/29.
//

#ifndef LSH_CPP_BENCHMARK_H
#define LSH_CPP_BENCHMARK_H

#include "minhash_benchmark.h"
#include "lsh_benchmark.h"
#include "weight_minhash_benchmark.h"
#include "dna_benchmark.h"

namespace LSH_CPP::Benchmark {
    void run_benchmark() {
        dna_benchmark();
        //lsh_benchmark();
        //weight_minhash_benchmark();
    }
}
#endif //LSH_CPP_BENCHMARK_H
