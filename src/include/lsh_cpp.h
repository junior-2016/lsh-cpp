//
// Created by junior on 19-7-22.
//

#ifndef LSH_CPP_LSH_CPP_H
#define LSH_CPP_LSH_CPP_H
// IO include
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// C++ std include
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <memory>
#include <functional>
#include <chrono>
#include <random>
#include <limits>

// simd include
#ifdef USE_SIMD

#include <xsimd/xsimd.hpp>

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
#include "../../third-party/matplotlib-cpp/matplotlibcpp.h"
#include "gsl/gsl_integration.h"
#include "../../third-party/parallel-hashmap/parallel_hashmap/phmap_fwd_decl.h"
#include "../../third-party/parallel-hashmap/parallel_hashmap/phmap.h"
#include "../../third-party/xxhash_cpp/xxhash/xxhash.h"
#include "../../third-party/xxhash_cpp/xxhash/xxhash.hpp"

// Self include
#include "config.h"

namespace LSH_CPP {

}
#endif //LSH_CPP_LSH_CPP_H
