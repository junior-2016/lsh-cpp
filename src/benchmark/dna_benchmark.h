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
    namespace CONFIG {
        const char *data_path = "../../dna-data/eco.dna";

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

        // 权重尽量提高 false negative weight (同时也就减少 false positive weight).
        // 即重点在于减少假阴性,这样可以最大范围得到正确的结果,就算因此引入了更多错误的结果也无所谓,
        // 因为可以进一步计算minhash的jaccard相似度来过滤错误结果.
        // 考虑到lsh得到的candidate set已经很小,因此最后的过滤操作并不会带来多大的开销.
        // 在 n_sample = 200 的情况下:
        //   w = {0.3,0.7}/{0.2,0.8} 时结果一样,因为参数优化下来都是r = 5; 结果文件大小 30.8 M; lsh查询加上过滤结果的时间最多也就一分钟左右;
        //   w = {0.1,0.9} 时, 参数变成 r = 4, 进一步缩小了 band 的宽度,也就是说变成 candidate set 的条件更松了,band的碰撞更多,
        //   这样可以大幅度避免正确的结果被遗漏,但同时也需要花更长的时间来过滤.不过测了一下lsh加上过滤的时间大概两分钟,所以带来的开销并不大.
        //   另外,w = {0.1,0.9} 时, 结果文件 36.4 M,已经非常逼近 minhash linear scan 的结果了.
        const std::pair<double, double> weights = {0.1, 0.9};

        const std::string output_parent_path = "out/";

        const std::string minhash_output_filename = output_parent_path
                                                    + "k=" + std::to_string(k)
                                                    + ",samples=" + std::to_string(n_sample)
                                                    + ",threshold=" + std::to_string(threshold)
                                                    + " minhash linear scan result";
        const std::string lsh_output_filename = output_parent_path
                                                + "k=" + std::to_string(k)
                                                + ",threshold=" + std::to_string(threshold)
                                                + ",samples=" + std::to_string(n_sample)
                                                + ",weights=[" + std::to_string(weights.first) + ","
                                                + std::to_string(weights.second) + "],"
                                                + " minhash lsh result";
        const std::string ground_truth_filename = output_parent_path
                                                  + "k=" + std::to_string(k)
                                                  + ",threshold=" + std::to_string(threshold)
                                                  + " ground truth result";
    }
    namespace DNA_DATA {
        using namespace CONFIG;
        std::vector<std::string> data;
        std::vector<size_t> labels;

        using MinHashType = MinHash<StdDNAShinglingHash64<k>, 32, n_sample>;
        using LSH_Type = LSH<XXUInt64Hash64, size_t, 0, 0, n_sample>;

        LSH_Type lsh(threshold, weights);
        std::vector<MinHashType> minhash_set;
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

        // 接受任意类型的整数参数,但实际写入的类型取决于 IntegerType.
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
            // get binary file size
            std::streampos file_size;
            in.seekg(0, std::ios::end);
            file_size = in.tellg();
            in.seekg(0, std::ios::beg);

            // read binary data
            std::vector<unsigned char> binary_data(static_cast<size_t>(file_size)); // 注意这里需要创建file_size大小的数组
            in.read((char *) (&binary_data[0]), file_size);

            // change binary data to integer array.
            std::vector<IntegerType> data;
            size_t step = sizeof(IntegerType) / sizeof(char);
            size_t pos = 0;
            while (pos < binary_data.size()) {
                IntegerType temp = 0;
                for (size_t i = 0; i < step; i++) temp |= (binary_data[pos + i] << (i * 8));
                data.push_back(temp);
                pos += step;
            }
            in.close();
            return data;
        }
    };

    void minhash_linear_scan_query() {
        std::cout << "run minhash linear scan ...\n";
        using namespace DNA_DATA;
        File<FileIO::Write, 16> out(minhash_output_filename);
        TimeVar start = timeNow();
        for (size_t i = 1; i < labels.size(); i++) {
            std::vector<size_t> temp;
            for (size_t j = 0; j < i; j++) {
                double sim = minhash_jaccard_similarity(minhash_set[labels[i]], minhash_set[labels[j]]);
                if (sim >= threshold) {
                    temp.push_back(labels[j]);
                }
            }
            out.write(labels[i]);
            out.write(temp.size());
            for (auto &item:temp) { out.write(item); }
        }
        std::cout << "time : " << second_duration((timeNow() - start)) << "seconds\n";
        out.close();
    }

    void lsh_query() {
        using namespace DNA_DATA;
        std::cout << "run lsh method ... \n";
        lsh.print_config();
        File<FileIO::Write, 16> out(lsh_output_filename);
        lsh.insert(minhash_set[labels[0]], labels[0]);
        TimeVar start = timeNow();
        for (size_t i = 1; i < labels.size(); i++) {
            auto ret = lsh.query_then_insert(minhash_set[labels[i]], labels[i]);
            std::vector<size_t> temp;
            for (const auto &item:ret) { // 对 ret 进一步过滤,ret储存的已经是label
                auto sim = minhash_jaccard_similarity(minhash_set[labels[i]], minhash_set[item]);
                if (sim >= threshold) {
                    temp.push_back(item);
                }
            }
            out.write(labels[i]);
            out.write(temp.size());
            for (auto &item:temp) { out.write(item); }
        }
        std::cout << "time : " << second_duration((timeNow() - start)) << "seconds \n";
        out.close();
    }

    void weight_minhash_query() {

    }

#define TEST_GROUND_TRUTH

#ifdef TEST_GROUND_TRUTH
    std::vector<HashSet < DNA_Shingling < CONFIG::k>>>
    dna_shingling_sets;

    void ground_truth_query() {
        using namespace DNA_DATA;
        File<FileIO::Write, 16> out(ground_truth_filename);
        for (size_t i = 1; i < dna_shingling_sets.size(); i++) {
            std::vector<size_t> temp;
            for (size_t j = 0; j < i; j++) {
                double sim = jaccard_similarity(dna_shingling_sets[i], dna_shingling_sets[j]);
                if (sim >= threshold) {
                    temp.push_back(labels[j]);
                }
            }
            out.write(labels[i]);
            out.write(temp.size());
            for (auto &item:temp) { out.write(item); }
        }
        out.close();
    }

    template<size_t sample_size, size_t random_sample_seed = 42>
    void sample_data() {
        using namespace DNA_DATA;
        std::mt19937_64 random_generator(random_sample_seed);
        std::vector<size_t> sample_labels;
        for (size_t i = 0; i < data.size(); i++) labels.push_back(i);
        std::sample(labels.begin(), labels.end(), std::back_inserter(sample_labels), sample_size, random_generator);
        minhash_set.reserve(data.size());
        dna_shingling_sets.reserve(sample_size);
        size_t pos = 0;
        for (size_t index = 0; index < sample_labels.size(); index++) {
            while (pos <= sample_labels[index]) {
                std::cout << "process doc " << pos << "...\n";
                auto dna_shingling_set = split_dna_shingling<k>(data[pos]);
                MinHashType temp(StdDNAShinglingHash64<k>{});
                temp.update(dna_shingling_set);
                minhash_set.push_back(temp);
                if (pos == sample_labels[index]) {
                    dna_shingling_sets.push_back(dna_shingling_set);
                }
                pos++;
            }
            if (index == sample_labels.size() - 1 && pos < data.size()) {
                while (pos < data.size()) {
                    std::cout << "process doc " << pos << "...\n";
                    MinHashType temp(StdDNAShinglingHash64<k>{});
                    temp.update(split_dna_shingling<k>(data[pos]));
                    minhash_set.push_back(temp);
                    pos++;
                }
            }
        }
        std::swap(sample_labels, labels); // 后面query函数用到的labels都是采样的labels
    }

    void compute_f_score() {
        std::cout << "compute f score ... \n";
        using namespace CONFIG;
        namespace plt = matplotlibcpp;
        constexpr size_t Integer_size = 16;
        using label_type = Integer<Integer_size>::type;
        using output_t =  std::vector<HashSet<label_type >>;
        output_t minhash_output, lsh_output, ground_truth_output;
        File<FileIO::Read, Integer_size> minhash_output_file(minhash_output_filename);
        File<FileIO::Read, Integer_size> lsh_output_file(lsh_output_filename);
        File<FileIO::Read, Integer_size> ground_truth_file(ground_truth_filename);
        auto minhash_file_data = minhash_output_file.read_data();
        auto lsh_file_data = lsh_output_file.read_data();
        auto ground_truth_data = ground_truth_file.read_data();
        auto process_data = [](const auto &data, output_t &output) -> void {
            for (size_t i = 0; i < data.size();) {
                HashSet<label_type> a;
                size_t a_size = data[++i];
                for (size_t index = 0; index < a_size; index++) a.insert(data[i + 1 + index]);
                output.push_back(a);
                i += (1 + a_size);
            }
        };
        process_data(minhash_file_data, minhash_output);
        process_data(lsh_file_data, lsh_output);
        process_data(ground_truth_data, ground_truth_output);
        size_t size = minhash_output.size();
        std::vector<double> minhash_f_scores, lsh_f_scores;
        for (size_t index = 0; index < size; index++) {
            auto minhash_pr = Statistic::get_precision_recall(minhash_output[index], ground_truth_output[index]);
            auto lsh_pr = Statistic::get_precision_recall(lsh_output[index], ground_truth_output[index]);
            auto minhash_f_score = Statistic::f_score(minhash_pr);
            auto lsh_f_score = Statistic::f_score(lsh_pr);
            minhash_f_scores.push_back(minhash_f_score);
            lsh_f_scores.push_back(lsh_f_score);
        }
        std::cout << "minhash linear scan f1 mean is : " << Statistic::get_mean(minhash_f_scores) << "\n";
        std::cout << "lsh f1 mean is : " << Statistic::get_mean(lsh_f_scores) << "\n\n";
        std::cout << "minhash linear scan f1 90% percentile is : "
                  << Statistic::get_percentile(minhash_f_scores, 0.9) << "\n";
        std::cout << "lsh f1 90% percentile is : " << Statistic::get_percentile(lsh_f_scores, 0.9) << "\n\n";
    }

#endif

    // [minhash linear scan] 与 [lsh(weight=0.1,0.9) + 过滤] 的比较中, lsh的结果比起minhash只少了一点
    // (因为lsh后面加上了过滤处理,所以lsh结果只会比minhash少,不会比minhash多),时间上lsh减少了一半,加速效率较好.
    // (这里的时间包括了处理和写入文件的所有时间)
    void dna_benchmark() {
        using namespace DNA_DATA;
        data = get_document_from_file(data_path);
#ifdef TEST_GROUND_TRUTH
        // TODO: check exist output file then skip sample_data and query process.
        sample_data<1000>();
        minhash_linear_scan_query();
        lsh_query();
        ground_truth_query();
        compute_f_score();
#else
        for (size_t i = 0; i < data.size(); i++) labels.push_back(i);
        minhash_set.reserve(data.size());
        int progress = 0;
        for (const auto &doc : data) {
            std::cout << "process " << ++progress << " doc...\n";
            auto dna_shingling_set = split_dna_shingling<k>(doc);
            MinHashType temp(StdDNAShinglingHash64<k>{});
            temp.update(dna_shingling_set);
            minhash_set.push_back(temp);
        }
        minhash_linear_scan_query();
        lsh_query();
#endif
    }
}
#endif //LSH_CPP_DNA_BENCHMARK_H
