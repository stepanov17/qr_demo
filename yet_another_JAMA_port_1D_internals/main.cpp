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
        std::cout << "reading data...  ";
        M.read(argv[1]);
    } catch (std::exception &ex) {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return 1;
    }
    std::cout << "done. matrix order is " << M.getN() << std::endl;
    std::cout << "calculating eigenvalues..." << std::endl;

    std::pair<SquareMatrix::data_t, SquareMatrix::data_t> ev;

    try {
        EigenvalueCalculator calc(M);
        ev = calc.getEigenvalues();
    } catch (std::exception &ex) {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return 1;
    }

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
