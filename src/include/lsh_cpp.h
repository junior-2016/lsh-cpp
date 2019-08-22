//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_LSH_CPP_H
#define LSH_CPP_LSH_CPP_H

// Eigen include and mkl macro setting
#define EIGEN_USE_MKL
#define EIGEN_USE_MKL_ALL
#define EIGEN_NO_DEBUG
#define MKL_NUM_THREADS 2
#define EIGEN_VECTORIZE_SSE4_2
#define EIGEN_VECTORIZE_AVX2
#define EIGEN_VECTORIZE_FMA
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "mkl.h"

// IO include
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

// C++ std include
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <bitset>
#include <algorithm>
#include <utility>
#include <memory>
#include <functional>
#include <chrono>
#include <random>
#include <limits>
#include <execution> // C++17 parallelism TS header

// C std include
#include <cstdint>
#include <cstdio>
#include <cstddef>
#include <ctime>
#include <cmath>
#include <cassert>

// Third party include
// phmap默认不使用带随机因子的哈希,所以这里不需要显式设置NON_DETERMINISTIC为0.(当前场景不需要随机性的哈希).
// #define PHMAP_NON_DETERMINISTIC 0
#include "../../third-party/matplotlib-cpp/matplotlibcpp.h"
#include "gsl/gsl_integration.h"
#include "../../third-party/parallel-hashmap/parallel_hashmap/phmap_fwd_decl.h"
#include "../../third-party/parallel-hashmap/parallel_hashmap/phmap.h"
#include "../../third-party/xxhash_cpp/xxhash/xxhash.h"
#include "../../third-party/xxhash_cpp/xxhash/xxhash.hpp"
#include "../../third-party/streaming-percentiles/cpp/src/include/stmpct/gk.hpp"

// Self include
#include "config.h"

namespace LSH_CPP {

}
#endif //LSH_CPP_LSH_CPP_H
