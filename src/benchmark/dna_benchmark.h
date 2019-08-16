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

        const std::string output_parent_path = "out/";

        /**
         * k=7 sim_threshold 取 0.5 得到一个相对适当的结果,但是取 0.6,0.7 结果急剧减少;
         * k<7 的时候, sim_threshold 不能取 0.5,否则结果过大,相当于每一个之间都存在相似;
         * 这时需要提高 sim_threshold -> 取 0.8 等,既可以抑制结果暴增,也可以得到一个更加正确(sim更大)的结果;
         * k>7 的时候,即使 sim_threshold 取0.5,都无法找到相似的结果,那么只能通过减少 sim threshold 来尝试获得结果,
         * 但是这样得到的结果不能反映太多的相似性,因为 sim threshold 低于0.5,已经不具有参考性.
         * 故而,应该在 k <=7 的范围内进行测试, 然后对于更小的k,应该使用更大的 sim threshold 来抑制结果增长.
         */
        constexpr size_t k = 6; // k取太小的话完全没有区分度,会导致所有的都相似,而且相似度阀值0.5可能太小
        const double threshold = 0.8;

        constexpr size_t n_sample = 256;

        using MinHashType = MinHash<StdDNAShinglingHash64<k>, 32, n_sample>;
        using LSH_Type = LSH<XXUInt64Hash64, size_t, 0, 0, n_sample>;

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

        std::vector<MinHashType> minhash_set;

        LSH_Type lsh(threshold, weights);

        std::vector<HashSet < DNA_Shingling < k>>>
        dna_shingling_sets;
    }

    template<size_t size>
    struct Integer {
        static_assert(size == 8 || size == 16 || size == 32 || size == 64);
    };
    template<>
    struct Integer<8> {
        using type = uint8_t;
    };
    template<>
    struct Integer<16> {
        using type = uint16_t;
    };
    template<>
    struct Integer<32> {
        using type = uint32_t;
    };
    template<>
    struct Integer<64> {
        using type = uint64_t;
    };

    struct FileIO {
        static constexpr size_t Read = 0;
        static constexpr size_t Write = 1;
    };

    template<size_t Flag, size_t IntegerSize>
    struct File {
        // Formally static_assert is a declaration. It is allowed wherever declarations are allowed.
        static_assert(Flag == FileIO::Read || Flag == FileIO::Write,
                      "You should use FileIO::Read or FileIO::Write Flag");
    };


    template<size_t IntegerSize>
    struct File<FileIO::Write, IntegerSize> {
        using IntegerType = typename Integer<IntegerSize>::type;
        std::ofstream out;

        explicit File(const std::string &filename) : out(filename, std::ios::binary) {}

        // 接受任意大小的整数参数,但是实际写入的范围取决于 IntegerType.
        template<typename NumberType>
        void write(NumberType number) {
            out.write(reinterpret_cast<char *>(&number), sizeof(IntegerType));
        }

        void close() { out.close(); }
    };

    template<size_t IntegerSize>
    struct File<FileIO::Read, IntegerSize> {
        using IntegerType = typename Integer<IntegerSize>::type;
        std::ifstream in;

        explicit File(const std::string &filename) : in(filename, std::ios::binary) {}

        std::vector<IntegerType> read_data() {
            std::streampos file_size;
            in.seekg(0, std::ios::end);
            file_size = in.tellg();
            in.seekg(0, std::ios::beg);
            std::vector<uint8_t> binary_data;
            std::vector<IntegerType> data;
            in.read((char *) (&binary_data[0]), file_size);
            int step = sizeof(IntegerType) / sizeof(char);
            int pos = 0;
            while (pos < binary_data.size()) {
                IntegerType temp = 0;
                for (int i = 0; i < step; i++) {
                    temp |= (binary_data[pos + i] << (i * 8));
                }
                data.push_back(temp);
                pos += step;
            }
            in.close();
            return data;
        }
    };

    /**
     * min hash size : 47910
     * exist sim count : 14460
     * max_sim_count : 7043
     */
    void minhash_linear_scan_query() {
        std::cout << "run minhash linear scan ...\n";
        using namespace DNA_DATA;
        std::string output_file_name = output_parent_path
                                       + "k=" + std::to_string(k)
                                       + ",samples=" + std::to_string(n_sample)
                                       + ",threshold=" + std::to_string(threshold)
                                       + " minhash linear scan result";
        File<FileIO::Write, 16> out(output_file_name);
//        int exist_sim = 0;
//        size_t max_sim_count = 0;
        TimeVar start = timeNow();
        for (size_t i = 1; i < minhash_set.size(); i++) {
//            std::cout << i << "...\n";
            std::vector<size_t> temp;
            for (size_t j = 0; j < i; j++) {
                double sim = minhash_jaccard_similarity(minhash_set[i], minhash_set[j]);
                if (sim >= threshold) {
                    temp.push_back(j);
                }
            }
            out.write(i);
            out.write(temp.size());
            for (auto &item:temp) { out.write(item); }
//            max_sim_count = std::max(max_sim_count, temp.size());
//            if (!temp.empty()) exist_sim++;
        }
        std::cout << "time : " << second_duration((timeNow() - start)) << "seconds\n";
        out.close();
//        std::cout << "min hash size : " << minhash_set.size() << "\n";
//        std::cout << "exist sim count : " << exist_sim << "\n";
//        std::cout << "max_sim_count : " << max_sim_count << "\n";
    }

    void ground_truth_query() {
        using namespace DNA_DATA;
        std::string output_file_name = output_parent_path
                                       + "k=" + std::to_string(k)
                                       + ",threshold=" + std::to_string(threshold)
                                       + " ground truth result";
        File<FileIO::Write, 16> out(output_file_name);
        for (size_t i = 1; i < dna_shingling_sets.size(); i++) {
            std::cout << i << "...\n";
            std::vector<size_t> temp;
            for (size_t j = 0; j < i; j++) {
                double sim = jaccard_similarity(dna_shingling_sets[i], dna_shingling_sets[j]);
                if (sim >= threshold) {
                    temp.push_back(j);
                }
            }
            out.write(i);
            out.write(temp.size());
            for (auto &item:temp) { out.write(item); }
        }
        out.close();
    }

    void lsh_query() {
        std::cout << "run lsh method ... \n";
        using namespace DNA_DATA;
        std::string output_file_name = output_parent_path
                                       + "k=" + std::to_string(k)
                                       + ",threshold=" + std::to_string(threshold)
                                       + ",samples=" + std::to_string(n_sample)
                                       + ",weights=[" + std::to_string(weights.first) + ","
                                       + std::to_string(weights.second) + "],"
                                       + " minhash lsh result";
        File<FileIO::Write, 16> out(output_file_name);
        lsh.insert(minhash_set[0], 0);
        TimeVar start = timeNow();
        for (size_t i = 1; i < minhash_set.size(); i++) {
            // std::cout << i << "...\n";
            auto ret = lsh.query_then_insert(minhash_set[i], i);
            std::vector<size_t> temp;
            for (const auto &item:ret) { // 对 ret 进一步过滤...
                auto sim = minhash_jaccard_similarity(minhash_set[i], minhash_set[item]);
                if (sim >= threshold) {
                    temp.push_back(item);
                }
            }
            out.write(i);
            out.write(temp.size());
            for (auto &item:temp) { out.write(item); }
        }
        std::cout << "time : " << second_duration((timeNow() - start)) << "seconds \n";
        out.close();
    }

    void weight_minhash_query() {

    }

    // 将输出文件转换为相似矩阵(用二进制保存)
    void output_file_transfer_matrix(const std::string &input_file_name, const std::string &output_file_name) {
        auto sim_content = get_document_from_file(input_file_name);
        size_t n = sim_content.size() + 1;
        std::vector<std::vector<char>> parse_content;
        parse_content.resize(n);
        for (size_t i = 0; i < n; i++) {
            parse_content[i].resize(n, '0');
        }
        for (const auto &line:sim_content) {
            std::istringstream iss(line);
            std::string s;
            iss >> s;
            size_t start = std::stoull(s);
            for (; iss >> s;) {
                size_t end = std::stoull(s);
                parse_content[start][end] = '1';
                parse_content[end][start] = '1';
            }
        }
        std::ofstream out(output_file_name);
        for (auto &i : parse_content) {
            for (auto &j : i) { out << j; }
            out << "\n";
        }
        out.close();
    }

    void transfer_matrix() {
        output_file_transfer_matrix("out/k=7,samples=256,threshold=0.500000 minhash linear scan result",
                                    "out/k=7,samples=256,threshold=0.500000 minhash linear scan matrix");
    }

    // TODO: 改变测试方法,将原数据集作为一个大数据集,从数据集上面随机采样10个1000大小的组成一个临时数据集,然后在上面测试
    //      minhash(LSH)与 ground truth.
    // TODO: 第二个优化点,储存结果文件,不要用数字字符串储存,那样太耗空间,直接储存一连串的uint16_t(最大支持65536大小)即可.
    //  => PS: 为了可扩展至任意大小编号,应该用bitset储存.

    // [minhash linear scan] 与 [lsh(weight=0.1,0.9) + 过滤] 的比较中, lsh的结果比起minhash只少了一点
    // (因为lsh后面加上了过滤处理,所以lsh结果只会比minhash少,不会比minhash多),时间上lsh减少了一半,加速效率较好.
    // (这里的时间包括了处理和写入文件的所有时间)

#undef TEST_GROUND_TRUTH

    void dna_benchmark() {
        using namespace DNA_DATA;
        const char *data_path = "../../dna-data/eco.dna";
        data = get_document_from_file(data_path);
#ifdef TEST_GROUND_TRUTH
        dna_shingling_sets.reserve(data.size());
#else
        minhash_set.reserve(data.size());
        lsh.print_config();
#endif
        int progress = 0;
        for (const auto &doc : data) {
            std::cout << "process " << ++progress << " doc...\n";
            auto dna_shingling_set = split_dna_shingling<k>(doc);
#ifdef TEST_GROUND_TRUTH
            dna_shingling_sets.push_back(dna_shingling_set);
#else
            MinHashType temp(StdDNAShinglingHash64<k>{});
            temp.update(dna_shingling_set);
            minhash_set.push_back(temp);
#endif
        }
#ifdef TEST_GROUND_TRUTH
        ground_truth_query();
#else
        minhash_linear_scan_query();
        lsh_query();
        minhash_set.clear(); // 注意 clear()
#endif
    }
}
#endif //LSH_CPP_DNA_BENCHMARK_H
