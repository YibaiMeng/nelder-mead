# Nelder Mead in C++

This is a simple C++ implementation of [Nelder Mead method](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method), a numerical method to find the minimum or maximum of an objective function in a multidimensional space. This project implements the Nelder Mead method in 
[Implementing the Nelder-Mead simplex algorithm with adaptive parameters](https://link.springer.com/article/10.1007/s10589-010-9329-3) by Gao and Han, which makes a modification to improve convergence in higher dimensions.

This is a header only library. To use this library, simply include the `nelder_mead.h` file in your program. The library defines a single function `find_min`.
```cpp
#include "nelder_mead.h"

double rosen(const std::vector<double>& x) {
   return std::pow((x[0] - x[1] * x[1]), 2) + std::pow(1 - x[0], 2);
}

std::vector<double> res = nelder_mead::find_min(rosen, std::vector<double>{2,4}); 
std::cout << res[0] << " " << res[1] << std::endl;
```
The output is 
```
1 1
```
See `example.cpp` for a more detailed example.

This is the detail prototype of the functions:
```cpp
template <typename Function, typename T = double>
std::vector<T> find_min(const Function &func, std::vector<T>& initial_point,
                        bool                               adaptive = false,
                        const std::vector<std::vector<T>> &initial_simplex = {},
                        T tol_fun = 1e-8, T tol_x = 1e-8,
                        unsigned int max_iter      = 1000000,
                        unsigned int max_fun_evals = 1000000);
```
Types:
- type `T`: the value type of the objective function.
- type `Function`: a function/function object/function pointer that takes a `std::vector<T>` and returns a `T`.

Arguments:
- `func`:  the objective function to optimize
- `initial_point`: the starting point of the process.
- `adaptive`: whether use the adaptive parameters as described in "Implementing the Nelder-Mead simplex algorithm with adaptive parameters"
- `initial_simplex`: set this when you don't wish to use the default generation method. It must have dimension + 1 elements. When setting this element, `initial_point` will be ignored.
- `tol_fun`, `tol_x`, `max_iter`, `max_fun_evals`: parameters for the termination of the process.

See the paper for a detailed description of the algorithm, including the termination criteria.
