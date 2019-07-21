#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <gsl/gsl_integration.h>

#include "../third-party/matplotlib-cpp/matplotlibcpp.h"
#include "../third-party/xxhash_cpp/xxhash/xxhash.hpp"

using namespace std;
namespace plt = matplotlibcpp;

/**
 * 测试 xxhash 32 bit 和 xxhash 64 bit 与 std::hash 的性能:
 * 64bit 机器下:
 *    xxhash 64 比 std::hash 快一倍;
 *    xxhash 32 和 std::hash 打平.
 * 32bit 机器不清楚,但是根据xxhash的文档,xxhash 32在32bit机器比在64bit机器更快,所以理论上,32bit机器上性能应该是:
 *    xxhash 64 优于 xxhash 32 优于 std::hash.
 * 应该默认采用 xxhash 64 以获得更好的性能.
 */
template<size_t N>
void test_xxhash() {
    static_assert(N == 32 || N == 64);
    std::cout << "N = " << N << "\n";
    string s;
    for (size_t i = 0; i <= 1000000; i++) {
        s += "a";
        auto ret = xxh::xxhash<N>(s);
        if (i % 100000 == 0) {
            std::cout << ret << "\n";
        }
    }
}

void test_std_hash() {
    string s;
    for (size_t i = 0; i <= 1000000; i++) {
        s += "a";
        auto ret = std::hash<std::string>{}(s);
        if (i % 100000 == 0) {
            std::cout << ret << "\n";
        }
    }
}

void test_matplotlib() {
    std::vector<std::vector<double>> x, y, z;
    for (double i = -5; i <= 5; i += 0.25) {
        std::vector<double> x_row, y_row, z_row;
        for (double j = -5; j <= 5; j += 0.25) {
            x_row.push_back(i);
            y_row.push_back(j);
            //z_row.push_back(sin(hypot(i, j)) + cos(hypot(i, j))); // sin(x^2+y^2)
            z_row.push_back(hypot(i, j));
        }
        x.push_back(x_row);
        y.push_back(y_row);
        z.push_back(z_row);
    }

    plt::plot_surface(x, y, z);
    plt::show();
}

/**
 * Non-capture lambda can be transfer to function pointer directly.
 * @tparam Func lambda expression type
 * @param f function for gsl_function
 */
template<typename Func>
void test_numerical_integration(Func &&f) {
    gsl_integration_workspace *w
            = gsl_integration_workspace_alloc(1000);

    double result, error;
    double expected = -4.0;
    double alpha = 1.0;

    gsl_function F;
    F.function = (double (*)(double, void *)) f;
    F.params = &alpha;

    gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000,
                         w, &result, &error);

    printf("result          = % .18f\n", result);
    printf("exact result    = % .18f\n", expected);
    printf("estimated error = % .18f\n", error);
    printf("actual error    = % .18f\n", result - expected);
    printf("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free(w);
}

int main() {
    test_numerical_integration([](double x, void *params) -> double {
        double alpha = *(double *) params;
        double f = log(alpha * x) / sqrt(x);
        return f;
    });
    return 0;
}