#include <cmath>
#include <iostream>
#include <vector>

#include "nelder-mead.h"

// The rosenbrock function
// https://en.wikipedia.org/wiki/Rosenbrock_function
template <unsigned int dimension, typename T>
double rosen(const std::vector<T> &m) {
    T res = 0;
    for (unsigned int i = 0; i < dimension - 1; i++) {
        res += 100 * (m[i + 1] - m[i] * m[i]) * (m[i + 1] - m[i] * m[i]) +
               (1 - m[i]) * (1 - m[i]);
    }
    return res;
}

int main() {
    try {
        auto starting_point   = std::vector<double>{1, 2, 3, 4, 5};
        auto starting_simplex = std::vector<std::vector<double>>{
            {2, 3, 4, 5, 1}, {3, 4, 5, 1, 2}, {4, 5, 1, 2, 3},
            {5, 1, 2, 3, 4}, {1, 2, 3, 4, 5}, {0, 0, 0, 0, 0}};
        auto res = nelder_mead::find_min(rosen<5, double>, starting_point, true,
                                         starting_simplex);
        std::cout << "Rosen solution: ";
        for (auto i : res)
            std::cout << i << " ";
        std::cout << std::endl;

    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
