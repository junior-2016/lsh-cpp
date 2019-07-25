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
/**
 * 对于 std::vector 等顺序容器的赋值过程,如果已经提前知道vector需要多少内存,并且不需要做任何初始化赋值
 * (先初始赋值,后面重新分配值会浪费时间,除非你的初始赋值是有作用的.)
 * 最优的策略是使用 reserve(N) 直接分配内存但是不初始化,然后紧接着用 push_back(v) 分配你需要的值.
 * 需要注意的是,如果用了reserve(N),分配值的时候必须用push_back,不能用operator[],否则size不会变化.
 */
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

// C std include
#include <cstdint>
#include <cstdio>
#include <cstddef>
#include <ctime>
#include <cmath>
#include <cassert>

// Third-party include
#define PHMAP_USE_ABSL_HASHEQ      // use absl::Hash as phmap hash framework

// phmap默认不使用带随机因子的哈希,所以这里不需要显式设置NON_DETERMINISTIC为0.
// 当前场景是不需要随机性的哈希的(只有在web服务器为了防止他人哈希攻击才需要引入随机因子)
// #define PHMAP_NON_DETERMINISTIC 0

#include "../../third-party/matplotlib-cpp/matplotlibcpp.h"
#include "gsl/gsl_integration.h"
#include "../../third-party/abseil-cpp/absl/hash/hash.h"
#include "../../third-party/parallel-hashmap/parallel_hashmap/phmap_fwd_decl.h"
#include "../../third-party/parallel-hashmap/parallel_hashmap/phmap.h"
#include "../../third-party/xxhash_cpp/xxhash/xxhash.h"
#include "../../third-party/xxhash_cpp/xxhash/xxhash.hpp"

// Self include
#include "config.h"

namespace LSH_CPP {

}
#endif //LSH_CPP_LSH_CPP_H
