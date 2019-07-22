//
// Created by junior on 19-7-20.
//

#ifndef LSH_CPP_LSH_H
#define LSH_CPP_LSH_H

#include "lsh_cpp.h"
#include "util.h"

namespace LSH_CPP {
    class LSH {
    private:
        size_t b, r, num_of_permutation;

        std::pair<size_t, size_t> optimal_params(double threshold, size_t num_of_permutation) {
            return {0, 0};
        }

    public:
        /**
         *
         * @param params = { b , r }
         * @param threshold: Jaccard similarity threshold. 0.0 <= threshold <= 1.0
         * @param num_of_permutation: number of min_hash permutation.
         * @param weights: { false_positive_weight, false_negative_weight }. weights.first + weights.second = 1.0.
         */
        explicit LSH(double threshold = 0.9,
                     size_t num_of_permutation = 128,
                     std::pair<double, double> weights = {1.0, 1.0},
                     std::pair<size_t, size_t> params = {16, 8}) : num_of_permutation(num_of_permutation) {
            assert(threshold >= 0 && threshold <= 1.0);
            assert(weights.first >= 0 && weights.second >= 0 && weights.first + weights.second == 1);

        }
    };
}
#endif //LSH_CPP_LSH_H
