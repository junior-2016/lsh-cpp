//
// Created by junior on 2019/7/28.
//

#ifndef LSH_CPP_DATA_LOADER_H
#define LSH_CPP_DATA_LOADER_H

#include "../include/lsh_cpp.h"
#include "../include/util.h"
#include "../include/minhash.h"

namespace LSH_CPP::Benchmark {
    namespace DataSpace {
        std::vector<std::string> population; // 采样池 => 变成全局变量,维持内部字符串生命周期.
        using DataLabel = std::vector<size_t>;
        using DataSet = std::vector<HashSet<std::string_view >>;
        using Data = std::pair<DataSet, DataLabel>;
        std::pair<Data, Data> global_data;
        bool is_data_clear = true;
    }

    /**
     * libstdc++ 实现 std::sample 使用了蓄水池抽样算法.
     * Question: 对于未知长度的序列,随机抽出k个数,并且要求每一个数抽出的概率一样.
     * Solution: 1. 建立k个数据的容器C, 用先来的k个数据直接填满.
     *           2. 对于后面来的数据 m ,生成一个随机数 j,如果 j < k, 则让m与C[j]交换.
     */
//    /// Reservoir sampling algorithm.
//    template<typename _InputIterator, typename _RandomAccessIterator,
//            typename _Size, typename _UniformRandomBitGenerator>
//    _RandomAccessIterator
//    __sample(_InputIterator __first, _InputIterator __last, input_iterator_tag,
//             _RandomAccessIterator __out, random_access_iterator_tag,
//             _Size __n, _UniformRandomBitGenerator&& __g)
//    {
//        using __distrib_type = uniform_int_distribution<_Size>;
//        using __param_type = typename __distrib_type::param_type;
//        __distrib_type __d{};
//        _Size __sample_sz = 0;
//        while (__first != __last && __sample_sz != __n)
//        {
//            __out[__sample_sz++] = *__first;
//            ++__first;
//        }
//        for (auto __pop_sz = __sample_sz; __first != __last;
//             ++__first, (void) ++__pop_sz)
//        {
//            const auto __k = __d(__g, __param_type{0, __pop_sz});
//            if (__k < __n)
//                __out[__k] = *__first;
//        }
//        return __out + __sample_sz;
//    }

    /**
     * @tparam population_size  采样池序列长度
     * @tparam data_set_size    数据集大小
     * @tparam data_size_generator_seed 数据大小随机生成器种子
     * @tparam data_sample_generator_seed 采样池采样数据的随机生成器种子
     * @param data_size_random_range 数据大小的随机范围
     * @return
     */
    template<size_t population_size,
            size_t train_set_size,
            size_t test_set_size,
            size_t data_size_generator_seed = 1,
            size_t data_sample_generator_seed = 42>
    bool bootstrap_data(const std::pair<size_t, size_t> &data_size_random_range) {
        if (!DataSpace::is_data_clear) {
            std::cerr << "You call bootstrap data but you don't clear previous data.\n";
            return false; // 返回操作失败
        }
        DataSpace::population.reserve(population_size);
        for_constexpr<for_bounds<0, population_size>>([&](auto index) {
            DataSpace::population.push_back(std::to_string(index));
        });
        DataSpace::DataSet train_set;
        train_set.resize(train_set_size, HashSet<std::string_view>{});
        DataSpace::DataSet test_set;
        test_set.resize(test_set_size, HashSet<std::string_view>{});
        std::mt19937_64 data_size_generator(data_size_generator_seed);
        std::mt19937_64 data_sample_generator(data_sample_generator_seed);
        std::uniform_int_distribution<size_t> data_size_dis(data_size_random_range.first,
                                                            data_size_random_range.second); // [a,b]
        for (size_t i = 0; i < train_set_size; i++) {
            auto data_size = data_size_dis(data_size_generator);
            HashSet<std::string_view> temp;
            std::sample(DataSpace::population.begin(), DataSpace::population.end(),
                        std::inserter(temp, temp.end()), data_size, data_sample_generator); // 注意set容器用std::inserter
            train_set[i] = temp;
        }
        std::vector<size_t> train_labels, test_labels;
        train_labels.reserve(train_set_size);
        for (size_t i = 0; i < train_set_size; i++) { train_labels.push_back(i); }
        std::sample(train_labels.begin(), train_labels.end(),
                    std::back_inserter(test_labels), test_set_size, data_sample_generator);
        for (size_t i = 0; i < test_set_size; i++) {
            test_set[i] = train_set[test_labels[i]];
        }
        DataSpace::Data train_data = {train_set, train_labels};
        DataSpace::Data test_data = {test_set, test_labels};
        DataSpace::global_data = {train_data, test_data};
        DataSpace::is_data_clear = false;
        return true;
    }

    void clear_data() {
        DataSpace::population.clear();
        DataSpace::population = std::vector<std::string>{};
        DataSpace::is_data_clear = true;
    }
}
#endif //LSH_CPP_DATA_LOADER_H
