//
// Created by junior on 2019/7/30.
//

#ifndef LSH_CPP_WEIGHT_MINHASH_H
#define LSH_CPP_WEIGHT_MINHASH_H

#include "lsh_cpp.h"
#include "util.h"

namespace LSH_CPP {
    // TODO Weight MinHash 算法. 这个算法必须提前知道集合的全集大小,对于基因相似度查询这个场景来说,
    //     Weight MinHash是可行的,比如取 k=9 的 k_mer,单个位置只有A/T/C/G四种可能,所以全集大小是 4^9,是确定的值.
    //     如果不用 WeightMinHash ,而是 MinHash ,那就必须先用 hash_set 过滤并且统计重复个数,这里就存在一个大量的
    //     hash key -> insert -> find 的操作开销; 而且 hash_set 对 集合元素进行哈希后储存, 接着 MinHash.update()
    //     还要对集合元素哈希为整数再做MinHash计算,这就两遍哈希了,开销极大.
    //     使用 Weight MinHash,可以不用 hash_set 过滤储存集合元素,直接用原始的 vector/array 储存即可
    //     (有重复元素也无所谓),然后借助一个公式来计算权重.
    //     paper: [1] Improved Consistent Sampling, Weighted Minhash and L1 Sketching
    //            [2] A Review for Weighted MinHash Algorithms (2018)
    //     算法思想:


    class WeightMinHash {

    };
}
#endif //LSH_CPP_WEIGHT_MINHASH_H
