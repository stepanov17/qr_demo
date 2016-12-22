
#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <string>

#include "SquareMatrix.h"
#include "QRIterations.h"


int main(int argc, char** argv) {

    if (argc < 2) {
        std::cerr << "usage: qriterations <matrix_file> [shift_type] [precision]" << std::endl;
        std::cerr << "shift_type = \"none\" | \"rayleygh\" | \"wilkinson\"" << std::endl;
        return 1;
    }

    QRIterations::shift shiftType = QRIterations::WILKINSON;
    std::string sType = "wilkinson";
    if (argc == 3) { // shift type provided
        sType = argv[2];
        std::transform(sType.begin(), sType.end(), sType.begin(), ::tolower);
        if (!sType.compare("none")) {
            shiftType = QRIterations::NONE;
        } else if (!sType.compare("rayleygh")) {
            shiftType = QRIterations::RAYLEIGH;
        } else if (!sType.compare("wilkinson")) {} // ok
        else {
            std::cerr << "invalid shift type: "<< sType << std::endl;
            return 1;
        }
    }

    double precision = 1.e-5;
    if (argc == 4) { // precision provided
        precision = std::stod(argv[3]);
        if (precision < 1.e-16) {
            std::cerr << "please use precsion >= 1.e-16" << std::endl;
            return 1;
        }
    }

    if (argc > 4) {
        std::cerr << "WARNING: extra arguments were found" << std::endl;
    }

    std::cout << "loading matrix..." << std::endl;
    SquareMatrix M;
    try {
        M.read(argv[1]);
        if (M.getN() == 0) {
            std::cerr << "empty matrix" << std::endl;
            return 1;
        }
    } catch (std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown exception" << std::endl;
        return 1;
    }
    std::cout << "done" << std::endl;


    std::cout << "matrix dimension is " <<
         M.getN() << "x" << M.getN() << std::endl;
    std::cout << "shift type = " << sType << std::endl;
    std::cout << "precision = " << precision << std::endl << std::endl;
    std::cout << "starting iterations..." << std::endl;


    try {
        QRIterations::run(M, shiftType, precision);
    } catch (std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown exception" << std::endl;
        return 1;
    }

    SquareMatrix::row_t diag = M.diag();

    std::cout << "done. main diagonal:" << std::endl;
    for (auto v: diag) { std::cout << v << " "; }
    std::cout << std::endl;
    // or
    // M.print();

    return 0;
}

