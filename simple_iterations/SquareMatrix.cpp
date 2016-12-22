
#include "SquareMatrix.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include <stdexcept>

#include <cmath>
#include <omp.h>
#include <limits>



SquareMatrix::SquareMatrix(std::size_t n_p) : n(n_p) {

    data.clear(); // just in case
    row_t tmp(n, 0.);
    data.assign(n, tmp);
}

SquareMatrix
SquareMatrix::operator *(const SquareMatrix &M) const {

    if (getN() != M.getN()) { throw std::invalid_argument(
        "matrix multiplication: inconsistent dimensions"); }

    SquareMatrix R(n);
    data_t &res = R.accessData();

    const data_t &other = M.getData();

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                res[i][j] += data[i][k] * other[k][j];
            }
        }
    }

    return R;
}

SquareMatrix&
SquareMatrix::operator =(const SquareMatrix &M) {

    if (n != M.n) { throw std::invalid_argument(
        "matrix assignment: inconsistent dimensions"); }

    if (this != &M) {
        //#pragma omp parallel for
        for (int i = 0; i < n; i++) { data[i] = M.data[i]; }
    }
    return *this;
}

SquareMatrix::row_t
SquareMatrix::diag() const {
    int n = getN();
    row_t d(n, 0.);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) { d[i] = data[i][i]; }
    return d;
}

SquareMatrix
SquareMatrix::Hessenberg() const {

    SquareMatrix H = *this;

    SquareMatrix::data_t &a = H.accessData();

    for (int k = 1; k < n - 1; k++) { // columns
        for (int i = k + 1; i < n; i++) { // rows

            // Givens rotation
            element_t h = std::hypot(a[k][k - 1], a[i][k - 1]);
            if (std::abs(h) <=
                std::numeric_limits<SquareMatrix::element_t>::epsilon()) {
                continue;
            }
            element_t c =  a[k][k - 1] / h;
            element_t s = -a[i][k - 1] / h;

            element_t ai[n], ak[n];

            // it makes little sense to add these pragmas here as the transform
            // is performed only once at the beginning of QR iterations

            #pragma omp parallel for
            for (int j = 0; j < n; j++) {
                ak[j] = c * a[k][j] - s * a[i][j];
                ai[j] = s * a[k][j] + c * a[i][j];
            }

            #pragma omp parallel for
            for (int j = 0; j < n; j++) {
                a[k][j] = ak[j];
                a[i][j] = ai[j];
            }

            #pragma omp parallel for
            for (int j = 0; j < n; j++) {
                ak[j] = c * a[j][k] - s * a[j][i];
                ai[j] = s * a[j][k] + c * a[j][i];
            }

            // copy temp arrays
            #pragma omp parallel for
            for (int j = 0; j < n; j++) {
                a[j][k] = ak[j];
                a[j][i] = ai[j];
            }
        }  // rows
    } // columns

    return H;
}

// see JAMA code
std::pair<SquareMatrix, SquareMatrix>
SquareMatrix::QR() const {

    SquareMatrix QM(n), RM(n);

    SquareMatrix::data_t &Q = QM.accessData();
    SquareMatrix::data_t &R = RM.accessData();

    SquareMatrix VM = *this;
    SquareMatrix::data_t &V = VM.accessData();

    row_t diagR(n, 0.);

    // main loop
    for (int k = 0; k < n; k++) {
        // compute 2-norm of k-th column
        element_t norm = 0.;
        for (int i = k; i < n; i++) { norm = std::hypot(norm, V[i][k]); }

        if (norm > 0.) {
            // form k-th Householder vector
            if (V[k][k] < 0.) { norm = -norm; }

            //slow down
            //#pragma omp parallel for shared(V)
            for (int i = k; i < n; i++) {
                V[i][k] /= norm;
            }
            V[k][k] += 1.;

            // apply transformation to remaining columns
            for (int j = k + 1; j < n; j++) {
                element_t s = 0.;

                //slow down
                //#pragma omp parallel for shared(V) reduction(+:s)
                for (int i = k; i < n; i++) {
                    s += V[i][k] * V[i][j];
                }
                s = -s / V[k][k];

                for (int i = k; i < n; i++) {
                    V[i][j] += s * V[i][k];
                }
            }
        } // norm > 0.
        diagR[k] = -norm;
    } // k

    // form matrix Q
    for (int k = n - 1; k >= 0; k--) {

        #pragma omp simd
        for (int i = 0; i < n; i++) { Q[i][k] = 0.; }

        Q[k][k] = 1.;

        for (int j = k; j < n; j++) {

            if (V[k][k] != 0.) {
                element_t s = 0.;

                //slow down
                //#pragma omp parallel for shared(V, Q) reduction(+:s)
                for (int i = k; i < n; i++) { s += V[i][k] * Q[i][j]; }
                s = -s / V[k][k];

                //slow down
                //#pragma omp parallel for shared(V)
                for (int i = k; i < n; i++) {
                    Q[i][j] += s * V[i][k];
                }
            }
        }
    } // Q

    // form matrix R
    #pragma omp parallel for //?
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i < j) {
                R[i][j] = V[i][j];
            } else if (i == j) {
                R[i][j] = diagR[i];
            } else {
                R[i][j] = 0.;
            }
        }
    }

    return std::make_pair(QM, RM);
}

void
SquareMatrix::print() const {

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << data[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void
SquareMatrix::read(const char *path) {

    data.clear();

    std::ifstream f;
    f.open(path, std::ios_base::in);
    if (!f.is_open()) { throw std::runtime_error("error reading file"); }

    std::string ln;
    bool firstLine = true;
    while (std::getline(f, ln)) {

        std::istringstream buff(ln);
        row_t row((std::istream_iterator<element_t>(buff)),
             std::istream_iterator<element_t>());
        if (firstLine) {
            firstLine = false;
            n = row.size();
            if (n < 1) {
                throw std::range_error("invalid matrix data");
            } else if (n > MAXN) {
                throw std::range_error("max order is exceeded");
            }
        } else if (n != row.size()) {
            throw std::range_error("invalid matrix data");
        }
        data.push_back(row);
    }
    if (data.size() != n) { throw std::range_error("the matrix must be square"); }

    f.close();
}
