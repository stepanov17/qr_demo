
#include <cmath>
#include <limits>

#include <iostream>
#include <stdexcept>

#include "QRIterations.h"

const SquareMatrix::element_t EPS = 
    std::numeric_limits<SquareMatrix::element_t>::epsilon();

SquareMatrix::element_t
QRIterations::getShift(const SquareMatrix& M, shift shiftType) {

    if (shift::NONE == shiftType) { return 0.; }

    int n = M.getN();
    if (shift::RAYLEIGH == shiftType) { return M.getData()[n - 1][n - 1]; }
    else if (shift::WILKINSON == shiftType) {

        const SquareMatrix::data_t &data = M.getData();
        SquareMatrix::element_t
            a = data[n - 1][n - 1], b = data[n - 2][n - 2], c = data[n - 1][n - 2];
        SquareMatrix::element_t d = 0.5 * (b - a);
        SquareMatrix::element_t tmp = std::abs(d) + std::hypot(d, c);
        if (tmp <= EPS) { return 0.; }
        double sgn = (d < 0.) ? -1. : 1.;
        return (a - sgn * c * c / tmp);
    }

    return 0; // not reachable
}

void
QRIterations::applyShift(SquareMatrix &M, SquareMatrix::element_t shift) {

    if (std::abs(shift) <= EPS) { return; } // nothing to do

    SquareMatrix::data_t &data = M.accessData();
    for (int i = 0, n = M.getN(); i < n; ++i) { data[i][i] += shift; }
}


void
QRIterations::run(SquareMatrix &M,
                  shift shiftAlgorithm,
                  SquareMatrix::element_t precision,
                  unsigned long maxIter) {

    if (precision <= EPS) {
        throw std::invalid_argument("precision is negative or too low");
    }

    int n = M.getN();
    if (n < 2) { return; } // nothing to do

    std::cout << "transforming to Hessenberg form..." << std::endl;
    M = M.Hessenberg();
    std::cout << "done" << std::endl << "iterating..." << std::endl;

    SquareMatrix::element_t shift = 0.;

    bool stop = false;
    for (unsigned int iter = 0; iter < maxIter; iter++) {

        SquareMatrix::row_t diag0 = M.diag();

        applyShift(M, -shift);
        std::pair<SquareMatrix, SquareMatrix> qr = M.QR();
        M = qr.second * qr.first;
        applyShift(M, shift);

        SquareMatrix::row_t diag = M.diag();

        stop = true;
        for (int i = 0; i < n; ++i) {
            SquareMatrix::element_t v = std::abs(diag0[i] - diag[i]);
            if (v > precision) {
                stop = false;
                break;
            }
        }
        if (stop) {
            std::cout << "stop at step #" << iter << std::endl;
            break;
        }

        shift = getShift(M, shiftAlgorithm);

        //std::cout << iter << "   " << shift << std::endl; // debug
    }
}
