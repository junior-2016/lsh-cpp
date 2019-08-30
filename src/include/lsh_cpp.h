//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_LSH_CPP_H
#define LSH_CPP_LSH_CPP_H

#include "portability.h"

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

#ifdef USE_CXX_PARALLISM_TS

#include <execution> // C++17 parallelism TS header,only be implemented in libstdc++(GCC 9)

#endif

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

// matplotlib-cpp 使用python2.7头文件,但是python2.7使用了register关键字,这个在std c++17已经弃用,
// clang编译器默认开"-Wregister",导致编译python2.7相关的头文件会报错,需要在导入matplotlib-cpp时显式忽略这个编译错误.
// 参考:(gcc也可以用同样的方法忽略某些警告带来的编译错误)
// https://stackoverflow.com/questions/49692794/c-17-compatability-with-python-2-7
// https://stackoverflow.com/questions/22422741/turning-off-the-register-storage-class-specifier-is-deprecated-warning
LSH_CPP_PUSH_WARNING
LSH_CPP_GNU_DISABLE_WARNING("-Wregister")
#include "../../third-party/matplotlib-cpp/matplotlibcpp.h"
LSH_CPP_POP_WARNING

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
