//
// Created by junior on 19-7-20.
//

#ifndef LSH_CPP_LSH_H
#define LSH_CPP_LSH_H

#include "lsh_cpp.h"
#include "minhash.h"
#include "util.h"
#include "hash.h"

namespace LSH_CPP {
    /**
     *
     * @tparam BandHashFunc 可以将Band中r个min_hash_value哈希为一个整数(uint64_t),作为当前Band哈希表的key
     * @tparam MinHashLabel  每一个 data pair 都由一个 min_hash_value_vector 和这个data的 label 组成.
     * 即 [ min_hash_value_vector, min_hash_label ] 数据对.
     * MinHashLabel可以是任意类型,比如数据的名称(字符串),或者数据的编号(整数),...,所以这里用模板控制类型.
     * 默认使用字符串作为数据的label,因为扩展性强,整数编号也可以用字符串表示.
     * 默认Label数据类型为 std::string_view. 目的是为了减小内存消耗(std::string_view仅持有数据地址和大小范围)
     * @tparam b number of bands
     * @tparam r number of rows in one band
     * @tparam n_permutation
     */
    template<
            typename BandHashFunc = XXUInt64Hash64,
            typename MinHashLabel = std::string_view,
            size_t b = 0,
            size_t r = 0,
            size_t n_permutation = 128
    >
    class LSH {
    public:
        // 通过BandHashFunc,将某个Band中r行min_hash哈希为一个整数,作为BandHashMap的key
        using BandHashKeyType = uint64_t;

        // BandHashMap以min_hash所在data的label作为value.因为需要通过冲突来寻找相似集合,
        // 所以实际的BandHashValueType是MinHashLabel集合
        using BandHashValueType = std::vector<MinHashLabel>;

        using BandHashMap = HashMap<BandHashKeyType, BandHashValueType>;

        using false_positive_weight = double;
        using false_negative_weight = double;

    private:
        std::pair<size_t, size_t> params = {0, 0}; // { b, r }
        std::vector<BandHashMap> band_hash_maps;                // 每一个Band持有的哈希表集合
        std::vector<std::pair<size_t, size_t>> band_hash_range; // BandHashFunc对min_hash_vector哈希的操作范围集合
        BandHashFunc bandHashFunc;

        // TODO: 引入 std::vector<MinHashLabel> data_set; 来检查有没有重复插入同一个data,但感觉不是特别必要...
        //      std::vector<MinHashLabel> data_set;

        std::pair<size_t, size_t>
        optimal_params(double threshold, const std::pair<false_positive_weight, false_negative_weight> &weights) {
            auto[false_positive_weight, false_negative_weight] = weights; // C++17 特性: auto[a,b] = std::pair<A,B>{a_val,b_val};
            double temp_params[2];
            double min_error = std::numeric_limits<double>::max();
            std::pair<size_t, size_t> ret;
            for (size_t _b = 1; _b <= n_permutation; _b++) {
                size_t max_r = n_permutation / _b;
                temp_params[0] = _b;
                for (size_t _r = 1; _r <= max_r; _r++) {
                    temp_params[1] = _r;
                    auto false_positive = numerical_integration
                            (false_positive_probability, {0.0, threshold}, temp_params);
                    auto false_negative = numerical_integration
                            (false_negative_probability, {threshold, 1.0}, temp_params);
                    auto error = false_positive_weight * false_positive + false_negative_weight * false_negative;
                    if (error < min_error) {
                        min_error = error;
                        ret.first = _b;
                        ret.second = _r;
                    }
                }
            }
            return ret;
        }

    public:
        /**
         * @param params = { b , r }
         * @param threshold: Jaccard similarity threshold. 0.0 <= threshold <= 1.0
         * @param weights: { false_positive_weight, false_negative_weight }. weights.first + weights.second = 1.0.
         * 比重给的越大,代表越希望减少这个方面的误差,比如极限情况下设置 { 0,1 }
         * 此时仅保留false_negative_weight,然后调用optimal_params(..)不断以减少false_negative_error去优化参数.
         */
        explicit LSH(double threshold = 0.9,
                     std::pair<false_positive_weight, false_negative_weight> weights = {0.5, 0.5}) {
            static_assert(n_permutation <= max_n_permutation);
            assert(threshold >= 0 && threshold <= 1.0);
            assert(weights.first >= 0 && weights.second >= 0 && weights.first + weights.second == 1);
            if constexpr (b > 0 && r > 0) { // 编译期判断提供的{b,r}是否可行
                static_assert(b * r <= n_permutation);
                band_hash_maps.resize(b, BandHashMap{});
                for (size_t i = 0; i < b; i++) {
                    band_hash_range.push_back({i * r, (i + 1) * r});
                }
            } else {
                // 根据 threshold 和 weights 得出优化的参数,保存在params.
                // 注意后面都需要先在编译期判断{b,r}是否可用,不可用时再选择params(且没有运行开销).
                params = optimal_params(threshold, weights);
                band_hash_maps.resize(params.first, BandHashMap{});
                for (size_t i = 0; i < params.first; i++) {
                    band_hash_range.push_back({i * params.second, (i + 1) * params.second});
                }
            }
        }

        template<typename HashFunc, size_t MinHashBits, size_t Seed, typename RandomGenerator>
        void insert(const MinHash <HashFunc, MinHashBits, n_permutation, Seed, RandomGenerator> &min_hash,
                    const MinHashLabel &label) {
            // TODO: 检查 label 代表的数据是不是重复插入
            for (size_t i = 0; i < band_hash_maps.size(); i++) {
                auto key = bandHashFunc(min_hash.hash_values, band_hash_range[i]);
                if (auto pos = band_hash_maps[i].find(key); pos == band_hash_maps[i].end()) { // c++17 if (init;cond)
                    band_hash_maps[i].insert({key, {label}});
                } else {
                    (*pos).second.push_back(label);
                }
            }
        }

        template<typename HashFunc, size_t MinHashBits, size_t Seed, typename RandomGenerator>
        HashSet <MinHashLabel>
        query_then_insert(const MinHash <HashFunc, MinHashBits, n_permutation, Seed, RandomGenerator> &min_hash,
                          const MinHashLabel &label) {
            // TODO: 检查 label 代表的数据是不是重复插入
            HashSet<MinHashLabel> candidate_set;
            for (size_t i = 0; i < band_hash_maps.size(); i++) {
                auto key = bandHashFunc(min_hash.hash_values, band_hash_range[i]);
                if (auto pos = band_hash_maps[i].find(key); pos != band_hash_maps[i].end()) {
                    for (const auto &item : (*pos).second) {
                        candidate_set.insert(item);
                    }
                    (*pos).second.push_back(label);
                } else {
                    band_hash_maps[i].insert({key, {label}});
                }
            }
            return candidate_set;
        }

        template<typename HashFunc, size_t MinHashBits, size_t Seed, typename RandomGenerator>
        HashSet <MinHashLabel> // 用hash_set做返回值是为了过滤重复的candidate.
        query(const MinHash <HashFunc, MinHashBits, n_permutation, Seed, RandomGenerator> &min_hash) {
            HashSet<MinHashLabel> candidate_set;
            for (size_t i = 0; i < band_hash_maps.size(); i++) {
                auto key = bandHashFunc(min_hash.hash_values, band_hash_range[i]);
                if (auto pos = band_hash_maps[i].find(key); pos != band_hash_maps[i].end()) {
                    for (const auto &item : (*pos).second) {
                        candidate_set.insert(item);
                    }
                }
            }
            return candidate_set;
        }

        void print_config() {
            std::cout << "===============  LSH config  ===============\n";
            std::cout << "params : b = " << params.first << "  r = " << params.second << "\n";
        }
    };
}
#endif //LSH_CPP_LSH_H
