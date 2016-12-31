
#include "SquareMatrix.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include <stdexcept>

#include <cmath>
#include <omp.h>
#include <limits>

// check if x is approximately zero
bool isZero(SquareMatrix::element_t x) {

    // hopefully this precision is enough
    // in principle, can use even lower value here, but
    // not sure if that makes sense
    SquareMatrix::element_t eps =
        std::numeric_limits<SquareMatrix::element_t>::epsilon();
    return (std::abs(x) < eps);
}

SquareMatrix::SquareMatrix(std::size_t n_p) : n(n_p) {

    data.clear(); // just in case
    data.assign(n * n, 0.);
}

SquareMatrix
SquareMatrix::operator *(const SquareMatrix &M) const {

    if (getN() != M.getN()) { throw std::invalid_argument(
        "matrix multiplication: inconsistent dimensions"); }

    SquareMatrix R(n);

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                R.at(i, j) += at(i, k) * M.at(k, j);
            }
        }
    }

    return R;
}

SquareMatrix&
SquareMatrix::operator =(const SquareMatrix &M) {

    if (n != M.n) { throw std::invalid_argument(
        "matrix assignment: inconsistent dimensions"); }

    if (this != &M) { data = M.data; }
    return *this;
}

SquareMatrix::data_t
SquareMatrix::diag() const {
    int n = getN();
    data_t d(n, 0.);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {  d[i] = at(i, i); }
    return d;
}

SquareMatrix
SquareMatrix::Hessenberg() const {

    SquareMatrix H = *this;

    for (int k = 1; k < n - 1; k++) { // columns
        for (int i = k + 1; i < n; i++) { // rows

            // Givens rotation
            element_t h = std::hypot(H.at(k, k - 1), H.at(i, k - 1));
            if (isZero(h)) { continue; }
            element_t c =  H.at(k, k - 1) / h;
            element_t s = -H.at(i, k - 1) / h;

            element_t ai[n], ak[n];

            #pragma omp parallel for
            for (int j = 0; j < n; j++) {
                ak[j] = c * H.at(k, j) - s * H.at(i, j);
                ai[j] = s * H.at(k, j) + c * H.at(i, j);
            }

            #pragma omp parallel for
            for (int j = 0; j < n; j++) {
                H.at(k, j) = ak[j];
                H.at(i, j) = ai[j];
            }

            #pragma omp parallel for
            for (int j = 0; j < n; j++) {
                ak[j] = c * H.at(j, k) - s * H.at(j, i);
                ai[j] = s * H.at(j, k) + c * H.at(j, i);
            }

            // copy temp arrays
            #pragma omp parallel for
            for (int j = 0; j < n; j++) {
                H.at(j, k) = ak[j];
                H.at(j, i) = ai[j];
            }
        }  // rows
    } // columns

    return H;
}

std::pair<SquareMatrix, SquareMatrix>
SquareMatrix::QR() const {

    SquareMatrix Q(n), R(n);

    SquareMatrix V = *this;

    data_t diagR(n, 0.);

    // main loop
    for (int k = 0; k < n; k++) {
        // compute 2-norm of k-th column
        element_t norm = 0.;
        for (int i = k; i < n; i++) { norm = std::hypot(norm, V.at(i, k)); }

        if (norm > 0.) {
            // form k-th Householder vector
            if (V.at(k, k) < 0.) { norm = -norm; }

            for (int i = k; i < n; i++) {
                V.at(i, k) /= norm;
            }
            V.at(k, k) += 1.;

            // apply transformation to remaining columns
            for (int j = k + 1; j < n; j++) {
                element_t s = 0.;

                for (int i = k; i < n; i++) {
                    s += V.at(i, k) * V.at(i, j);
                }
                s = -s / V.at(k, k);

                for (int i = k; i < n; i++) {
                    V.at(i, j) += s * V.at(i, k);
                }
            }
        } // norm > 0.
        diagR[k] = -norm;
    } // k

    // form matrix Q
    for (int k = n - 1; k >= 0; k--) {

        #pragma omp simd
        for (int i = 0; i < n; i++) { Q.at(i, k) = 0.; }

        Q.at(k, k) = 1.;

        for (int j = k; j < n; j++) {

            if (!isZero(V.at(k, k)) != 0.) {
                element_t s = 0.;

                for (int i = k; i < n; i++) { s += V.at(i, k) * Q.at(i, j); }
                s = -s / V.at(k, k);

                for (int i = k; i < n; i++) {
                    Q.at(i, j) += s * V.at(i, k);
                }
            }
        }
    } // Q

    // form matrix R
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i < j) {
                R.at(i, j) = V.at(i, j);
            } else if (i == j) {
                R.at(i, j) = diagR[i];
            } else {
                R.at(i, j) = 0.;
            }
        }
    }

    return std::make_pair(Q, R);
}

void
SquareMatrix::read(const char *path) {

    data.clear();

    std::ifstream f;
    f.open(path, std::ios_base::in);
    if (!f.is_open()) { throw std::runtime_error("error reading file"); }

    std::string ln;
    bool firstLine = true;
    int nRows = 0;
    while (std::getline(f, ln)) {

        std::istringstream buff(ln);
        data_t row((std::istream_iterator<element_t>(buff)),
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
        data.insert(data.end(), row.begin(), row.end());
        ++nRows;
    }
    if (nRows != n) { throw std::range_error("the matrix must be square"); }
    else if (nRows == 0) { throw std::range_error("empty matrix"); }

    f.close();
}


std::ostream&
operator <<(std::ostream& os, const SquareMatrix &M) {

    int n = M.getN();
    os << n << "x" << n << " matrix" << std::endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            os << M.at(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}
