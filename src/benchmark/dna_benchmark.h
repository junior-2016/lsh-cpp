//
// Created by junior on 2019/8/10.
//

#ifndef LSH_CPP_DNA_BENCHMARK_H
#define LSH_CPP_DNA_BENCHMARK_H

#include "../include/lsh_cpp.h"
#include "../include/util.h"
#include "../include/k_shingles.h"
#include "../include/io.h"
#include "../include/minhash.h"
#include "../include/lsh.h"

namespace LSH_CPP::Benchmark {
    namespace DNA_DATA {
        std::vector<std::string> data;
        constexpr size_t k = 7; // k取太小的话完全没有区分度,会导致所有的都相似,而且相似度阀值0.5可能太小
        constexpr size_t n_sample = 200;
        const double threshold = 0.5;
        using MinHashType = MinHash<XXKShinglingHash64, 32, n_sample>;
        using LSH_Type = LSH<XXUInt64Hash64, size_t, 0, 0, n_sample>;

        std::vector<MinHashType> minhash_set;

        // 权重尽量提高 false negative weight (同时也就减少 false positive weight).
        // 即重点在于减少假阴性,这样可以最大范围得到正确的结果,就算因此引入了更多错误的结果也无所谓,
        // 因为可以进一步计算minhash的jaccard相似度来过滤错误结果.
        // 考虑到lsh得到的candidate set已经很小,因此最后的过滤操作并不会带来多大的开销.
        // 在 n_sample = 200 的情况下:
        //   w = {0.3,0.7}/{0.2,0.8} 时结果一样,因为参数优化下来都是r = 5; 结果文件大小 30.8 M; lsh查询加上过滤结果的时间最多也就一分钟左右;
        //   w = {0.1,0.9} 时, 参数变成 r = 4, 进一步缩小了 band 的宽度,也就是说变成 candidate set 的条件更松了,band的碰撞更多,
        //   这样可以大幅度避免正确的结果被遗漏,但同时也需要花更长的时间来过滤.不过测了一下lsh加上过滤的时间最多两分钟,所以带来的开销并不大.
        //   另外,w = {0.1,0.9} 时, 结果文件 36.4 M,已经非常逼近 minhash linear scan 的结果了.
        const std::pair<double, double> weights = {0.1, 0.9};

        LSH_Type lsh(threshold, weights);

        // std::vector<HashSet < K_shingling>> k_shingling_sets;

        // TODO 碱基压缩为2bit,用std::bit_set<> 储存数据,然后计算.
        //      基于上面的压缩方式,如果k取7,那么一个k_shingling只需要 2*7 = 14 bit,用 uint16_t 储存刚刚好,头两位可以不用.
        //      然后这个 16bit 整数表示的 k_shingling 就可以直接作为 shingling 的 哈希值,不需要用 xxhash 或者其他字符串哈希函数处理,
        //      直接作为哈希值已经可以了,因为不同的shingling压缩的码肯定都不一样(当然也可以直接把bit stream提供给哈希函数进一步哈希);
        //      weight minhash 的位置编码也可以直接用这个压缩的整数值.
        //#define A 0x00
        //#define T 0x01
        //#define C 0x10
        //#define G 0x11
    }

    /**
     * min hash size : 47910
     * exist sim count : 14460
     * max_sim_count : 7043
     */
    void minhash_linear_scan_query() {
        std::cout << "run minhash linear scan ...\n";
        using namespace DNA_DATA;
        std::string output_file_name = "k=" + std::to_string(k) + ",threshold=" + std::to_string(threshold)
                                       + " minhash linear scan result";
        std::ofstream out(output_file_name);
        int exist_sim = 0;
        int max_sim_count = 0;
        for (size_t i = 1; i < minhash_set.size(); i++) {
            out << i << " ";
            int sim_count = 0;
            for (size_t j = 0; j < i; j++) {
                double sim = minhash_jaccard_similarity(minhash_set[i], minhash_set[j]);
                if (sim >= threshold) {
                    sim_count++;
                    out << j << " ";
                }
            }
            max_sim_count = std::max(max_sim_count, sim_count);
            if (sim_count > 0) exist_sim++;
            out << "\n";
        }
        out.close();
        std::cout << "min hash size : " << minhash_set.size() << "\n";
        std::cout << "exist sim count : " << exist_sim << "\n";
        std::cout << "max_sim_count : " << max_sim_count << "\n";
    }

//    void ground_truth_query(double threshold = 0.5) {
//        using namespace DNA_DATA;
//        std::string output_file_name = "k=" + std::to_string(k) + ",threshold=" + std::to_string(threshold)
//                                       + " ground truth result";
//        std::ofstream out(output_file_name);
//        for (size_t i = 1; i < k_shingling_sets.size(); i++) {
//            out << i << " ";
//            for (size_t j = 0; j < i; j++) {
//                double sim = jaccard_similarity(k_shingling_sets[i], k_shingling_sets[j]);
//                if (sim >= threshold) {
//                    out << j << " ";
//                }
//            }
//            out << "\n";
//        }
//    }

    void weight_minhash_query() {

    }

    void lsh_query() {
        std::cout << "run lsh method ... \n";
        using namespace DNA_DATA;
        std::string output_file_name = "k=" + std::to_string(k) + ",threshold=" + std::to_string(threshold)
                                       + ",weights=[" + std::to_string(weights.first) + "," +
                                       std::to_string(weights.second) + "],"
                                       + " minhash lsh result";
        std::ofstream out(output_file_name);
        lsh.insert(minhash_set[0], 0);
        for (size_t i = 1; i < minhash_set.size(); i++) {
            auto ret = lsh.query_then_insert(minhash_set[i], i);
            out << i << " ";
            for (const auto &item:ret) {
                // 可以对 ret 进一步过滤...
                auto sim = minhash_jaccard_similarity(minhash_set[i], minhash_set[item]);
                if (sim >= threshold) {
                    out << item << " ";
                }
            }
            out << "\n";
        }
        out.close();
    }

    void dna_benchmark() {
        using namespace DNA_DATA;
        const char *data_path = "../../dna-data/eco.dna";
        data = get_document_from_file(data_path);
        minhash_set.reserve(data.size());
        lsh.print_config();
        int progress = 0;
        for (const auto &doc : data) {
            std::cout << "process " << ++progress << " doc...\n";
            auto k_shingling_set = split_k_shingling_fast(doc, k);
            MinHashType temp(XXKShinglingHash64{});
            temp.update(k_shingling_set);
            minhash_set.push_back(temp);
        }
        // minhash_linear_scan_query();
        lsh_query();
        // ground_truth_query();
        minhash_set.clear(); // 注意 clear()
    }
}
#endif //LSH_CPP_DNA_BENCHMARK_H
