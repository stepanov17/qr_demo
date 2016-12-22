#include <iostream>

#include "SquareMatrix.h"
#include "EigenvalueCalculator.h"

int main(int argc, char** argv) {

    if (argc != 2) {
        std::cout << "usage: eigcalc <input_file>" << std::endl;
        return 1;
    }

    SquareMatrix M;
    try {
        M.read(argv[1]);
    } catch (std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    EigenvalueCalculator calc(M);
    std::pair<SquareMatrix::row_t, SquareMatrix::row_t> ev = calc.getEigenvalues();

    // check real parts
    std::cout << "real parts:" << std::endl;
    for (auto v: ev.first) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    // check imag parts
    std::cout << "imaginary parts:" << std::endl;
    for (auto v: ev.second) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    return 0;
}

