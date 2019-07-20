#include <iostream>
#include <vector>
#include "../third-party/matplotlib-cpp/matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

int main() {
    std::vector<std::vector<double>> x, y, z;
    for (double i = -10; i <= 10; i += 0.25) {
        std::vector<double> x_row, y_row, z_row;
        for (double j = -10; j <= 10; j += 0.25) {
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
    return 0;
}