
#include "EigenvalueCalculator.h"

#include <cmath>
#include <limits>
#include <iostream>


// utils

bool
EigenvalueCalculator::isEq(element_t x, element_t y) {

    return (std::abs(x - y) < std::numeric_limits<element_t>::epsilon());
}

SquareMatrix::element_t
EigenvalueCalculator::hypot(element_t a, element_t b) {

    element_t r;
    if (std::abs(a) > std::abs(b)) {
        r = b / a;
        r = std::abs(a) * std::sqrt(1 + r * r);
    } else if (!isEq(b, 0.)) {
        r = a / b;
        r = std::abs(b) * std::sqrt(1 + r * r);
    } else {
        r = 0.;
    }
    return r;
}

std::pair<SquareMatrix::element_t, SquareMatrix::element_t>
EigenvalueCalculator::cdiv(element_t xr, element_t xi, element_t yr, element_t yi) {

    element_t r, v, cdivr, cdivi;
    if (std::abs(yr) > std::abs(yi)) {
        r = yi / yr;
        v = yr + r * yi;
        cdivr = (xr + r * xi) / v;
        cdivi = (xi - r * xr) / v;
    } else {
        r = yr / yi;
        v = yi + r * yr;
        cdivr = (r * xr + xi) / v;
        cdivi = (r * xi - xr) / v;
    }
    return std::make_pair(cdivr, cdivi);
}

EigenvalueCalculator::EigenvalueCalculator(const SquareMatrix &M) {

    n = M.getN();
    if (n > MAXN) {
        std::cerr << "the matrix order should not exceed " << MAXN << std::endl;
        std::cerr << "skipping calculations!" << std::endl << std::endl;
        n = 0;
        V = SquareMatrix();
    } else {
        this->M = M;
    }


    V = SquareMatrix(n);
    H = SquareMatrix(n);

    eigIm = eigRe = SquareMatrix::data_t(n, 0.);
}

// utils



std::pair<SquareMatrix::data_t, SquareMatrix::data_t>
EigenvalueCalculator::getEigenvalues() {

    if (n == 0) { // empty matrix
        return std::make_pair(SquareMatrix::data_t(0), SquareMatrix::data_t(0));
    } else if (n == 1) { // nothing to do
        element_t v = M.at(0, 0);
        return std::make_pair(
            SquareMatrix::data_t(1, v), SquareMatrix::data_t(1, 0.));
    }

    if (M.isSymmetric()) {
        V = M;
        Tridiagonalize();
        Diagonalize();

    } else {
        H = M;
        toHessenberg(); // reduce to Hessenberg form
        toSchur(); // reduce to real Schur form.
    }

    return std::make_pair(eigRe, eigIm);
}

void
EigenvalueCalculator::Tridiagonalize() {

    for (int j = 0; j < n; j++) { eigRe[j] = V.at(n - 1, j); }

    // Householder reduction to tridiagonal form
    for (int i = n - 1; i > 0; i--) {

        // scale to avoid under/overflow
        element_t scale = 0.;
        element_t h = 0.;
        for (int k = 0; k < i; k++) {
            scale += std::abs(eigRe[k]);
        }
        if (isEq(scale, 0.)) {
            eigIm[i] = eigRe[i - 1];
            for (int j = 0; j < i; j++) {
                eigRe[j] = V.at(i - 1, j);
                V.at(i, j) = 0.;
                V.at(j, i) = 0.;
            }
        } else {

            // generate Householder vector
            for (int k = 0; k < i; k++) {
                eigRe[k] /= scale;
                h += eigRe[k] * eigRe[k];
            }

            element_t f = eigRe[i - 1];
            element_t g = std::sqrt(h);
            if (f > 0) { g = -g; }
            eigIm[i] = scale * g;
            h = h - f * g;
            eigRe[i - 1] = f - g;
            for (int j = 0; j < i; j++) { eigIm[j] = 0.; }

            // apply similarity transformation to remaining columns
            for (int j = 0; j < i; j++) {
                f = eigRe[j];
                V.at(j, i) = f;
                g = eigIm[j] + V.at(j, j) * f;
                for (int k = j + 1; k <= i - 1; k++) {
                    g += V.at(k, j) * eigRe[k];
                    eigIm[k] += V.at(k, j) * f;
                }
                eigIm[j] = g;
            }
            f = 0.;
            for (int j = 0; j < i; j++) {
                eigIm[j] /= h;
                f += eigIm[j] * eigRe[j];
            }
            element_t hh = f / (h + h);
            for (int j = 0; j < i; j++) {
                eigIm[j] -= hh * eigRe[j];
            }
            for (int j = 0; j < i; j++) {
                f = eigRe[j];
                g = eigIm[j];
                for (int k = j; k <= i - 1; k++) {
                    V.at(k, j) -= (f * eigIm[k] + g * eigRe[k]);
                }
                eigRe[j] = V.at(i - 1, j);
                V.at(i, j) = 0.;
            }
        }
        eigRe[i] = h;
    }

    // accumulate transformations
    //#pragma omp parallel for
    for (int i = 0; i < n - 1; i++) {
        V.at(n - 1, i) = V.at(i, i);
        V.at(i, i) = 1.;
        element_t h = eigRe[i + 1];
        if (!isEq(h, 0.)) {
            #pragma omp parallel for
            for (int k = 0; k <= i; k++) {
                eigRe[k] = V.at(k, i + 1) / h;
            }
            for (int j = 0; j <= i; j++) {
                element_t g = 0.;
                #pragma omp parallel for reduction(+:g)
                for (int k = 0; k <= i; k++) {
                    g += V.at(k, i + 1) * V.at(k, j);
                }
                #pragma omp parallel for
                for (int k = 0; k <= i; k++) {
                    V.at(k, j) -= g * eigRe[k];
                }
            }
        }
        #pragma omp parallel for
        for (int k = 0; k <= i; k++) { V.at(k, i + 1) = 0.; }
    }

    #pragma omp parallel for
    for (int j = 0; j < n; j++) {
        eigRe[j] = V.at(n - 1, j);
        V.at(n - 1, j) = 0.;
    }
    V.at(n - 1, n - 1) = 1.;
    eigIm[0] = 0.;
}

void
EigenvalueCalculator::Diagonalize() {

    for (int i = 1; i < n; i++) { eigIm[i - 1] = eigIm[i]; }
    eigIm[n - 1] = 0.;

    element_t f = 0.;

    for (int l = 0; l < n; l++) {

        // find small subdiagonal element
        element_t tst1 = std::max<element_t>(tst1, std::abs(eigRe[l]) + std::abs(eigIm[l]));
        int m = l;
        while (m < n) {
            if (std::abs(eigIm[m]) <= EPS * tst1) { break; }
            m++;
        }

        // if m == l, eigRe[l] is an eigenvalue
        // otherwise, iterate
        if (m > l) {
            int iter = 0;
            do {
                iter = iter + 1;  // (could check iteration count here)
                // compute implicit shift
                element_t g = eigRe[l];
                element_t p = (eigRe[l + 1] - g) / (2. * eigIm[l]);
                element_t r = hypot(p, 1.);
                if (p < 0) { r = -r; }
                eigRe[l] = eigIm[l] / (p + r);
                eigRe[l + 1] = eigIm[l] * (p + r);
                element_t dl1 = eigRe[l + 1];
                element_t h = g - eigRe[l];
                for (int i = l + 2; i < n; i++) {
                    eigRe[i] -= h;
                }

                f += h;

                // implicit QL transformation
                p = eigRe[m];
                element_t c = 1.;
                element_t c2 = c;
                element_t c3 = c;
                element_t el1 = eigIm[l + 1];
                element_t s = 0.;
                element_t s2 = 0.;
                for (int i = m - 1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * eigIm[i];
                    h = c * p;
                    r = hypot(p, eigIm[i]);
                    eigIm[i + 1] = s * r;
                    s = eigIm[i] / r;
                    c = p / r;
                    p = c * eigRe[i] - s * g;
                    eigRe[i + 1] = h + s * (c * g + s * eigRe[i]);

                    // Accumulate transformation.
                    for (int k = 0; k < n; k++) {
                        h = V.at(k, i + 1);
                        V.at(k, i + 1) = s * V.at(k, i) + c * h;
                        V.at(k, i) = c * V.at(k, i) - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * eigIm[l] / dl1;
                eigIm[l] = s * p;
                eigRe[l] = c * p;

            // check for convergence
            } while (std::abs(eigIm[l]) > EPS * tst1);
        }
        eigRe[l] = eigRe[l] + f;
        eigIm[l] = 0.;
    }

    // sort eigenvalues and corresponding vectors
    #pragma omp parallel for
    for (int i = 0; i < n - 1; i++) {
        int k = i;
        element_t p = eigRe[i];
        for (int j = i + 1; j < n; j++) {
            if (eigRe[j] < p) {
                k = j;
                p = eigRe[j];
            }
        }
        if (k != i) {
            eigRe[k] = eigRe[i];
            eigRe[i] = p;
            for (int j = 0; j < n; j++) {
                p = V.at(j, i);
                V.at(j, i) = V.at(j, k);
                V.at(j, k) = p;
            }
        }
    }
}

void
EigenvalueCalculator::toHessenberg() {

    int low = 0;
    int high = n - 1;

    SquareMatrix::data_t ort(n, 0.);

    for (int m = low + 1; m <= high - 1; m++) {

        // scale column
        element_t scale = 0.;
        #pragma omp parallel for reduction(+:scale)
        for (int i = m; i <= high; i++) {
            scale += std::abs(H.at(i, m - 1));
        }
        if (!isEq(scale, 0.)) {

            // compute Householder transformation

            element_t h = 0.;
            //#pragma omp parallel for
            for (int i = high; i >= m; i--) {
                element_t o = H.at(i, m - 1) / scale;
                ort[i] = o;
                h += o * o;
            }

            element_t g = std::sqrt(h);
            if (ort[m] > 0) {
                g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;

            // apply Householder similarity transformation

            #pragma omp parallel for
            for (int j = m; j < n; j++) {
                element_t f = 0.;

                for (int i = high; i >= m; i--) {
                    f += ort[i] * H.at(i, j);
                }
                f = f / h;
                for (int i = m; i <= high; i++) {
                    H.at(i, j) -= f * ort[i];
                }
            }

            #pragma omp parallel for
            for (int i = 0; i <= high; i++) {
                element_t f = 0.;
                for (int j = high; j >= m; j--) {
                    f += ort[j] * H.at(i, j);
                }
                f = f / h;
                for (int j = m; j <= high; j++) {
                    H.at(i, j) -= f * ort[j];
                }
            }
            ort[m] = scale * ort[m];
            H.at(m, m - 1) = scale * g;
        }
    }

    // accumulate transformations (Algol's ortran).

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            V.at(i, j) = (i == j ? 1. : 0.);
        }
    }

    //#pragma omp parallel for
    for (int m = high - 1; m >= low + 1; m--) {
        if (!isEq(H.at(m, m - 1), 0.)) {
            #pragma omp parallel for
            for (int i = m + 1; i <= high; i++) {
                ort[i] = H.at(i, m - 1);
            }
            for (int j = m; j <= high; j++) {

                element_t g = 0.;
                #pragma omp parallel for reduction(+:g)
                for (int i = m; i <= high; i++) {
                    g += ort[i] * V.at(i, j);
                }
                // double division avoids possible underflow
                g = (g / ort[m]) / H.at(m, m - 1);
                for (int i = m; i <= high; i++) {
                    V.at(i, j) += g * ort[i];
                }
            }
        }
    }
}

void
EigenvalueCalculator::toSchur() {

    // initialize
    int nn = n;
    int n0 = nn - 1;
    int low = 0;
    int high = nn - 1;
    element_t exshift = 0.;
    element_t p = 0., q = 0., r = 0., s = 0., z = 0., t, w, x, y;

    // store roots isolated by balance and compute matrix norm
    element_t norm = 0.;
    for (int i = 0; i < nn; i++) {
        if (i < low || i > high) {
            eigRe[i] = H.at(i, i);
            eigIm[i] = 0.;
        }
        #pragma omp parallel for reduction (+:norm)
        for (int j = std::max(i - 1, 0); j < nn; j++) {
            norm += std::abs(H.at(i, j));
        }
    }

    // outer loop over eigenvalue index
    int iter = 0;
    while (n0 >= low) {

        // look for single small sub-diagonal element
        int l = n0;
        while (l > low) {
            s = std::abs(H.at(l - 1, l - 1)) + std::abs(H.at(l, l));
            if (isEq(s, 0.)) {
                s = norm;
            }
            if (std::abs(H.at(l, l - 1)) < EPS * s) {
                break;
            }
            l--;
        }

        // check for convergence
        // one root found
        if (l == n0) {
            H.at(n0, n0) = H.at(n0, n0) + exshift;
            eigRe[n0] = H.at(n0, n0);
            eigIm[n0] = 0.;
            n0--;
            iter = 0;
        // two roots found
        } else if (l == n0 - 1) {
            w = H.at(n0, n0 - 1) * H.at(n0 - 1, n0);
            p = (H.at(n0 - 1, n0 - 1) - H.at(n0, n0)) / 2.;
            q = p * p + w;
            z = std::sqrt(std::abs(q));
            H.at(n0, n0) = H.at(n0, n0) + exshift;
            H.at(n0 - 1, n0 - 1) = H.at(n0 - 1, n0 - 1) + exshift;
            x = H.at(n0, n0);

            // real pair
            if (q >= 0) {
                if (p >= 0) {
                    z = p + z;
                } else {
                    z = p - z;
                }
                eigRe[n0 - 1] = x + z;
                eigRe[n0] = eigRe[n0 - 1];
                if (!isEq(z, 0.)) {
                    eigRe[n0] = x - w / z;
                }
                eigIm[n0 - 1] = 0.;
                eigIm[n0] = 0.;
                x = H.at(n0, n0 - 1);
                s = std::abs(x) + std::abs(z);
                p = x / s;
                q = z / s;
                r = std::sqrt(p * p + q * q);
                p = p / r;
                q = q / r;

                // row modification
                #pragma omp parallel for
                for (int j = n0 - 1; j < nn; j++) {
                    element_t tmp = H.at(n0 - 1, j);
                    H.at(n0 - 1, j) = q * tmp + p * H.at(n0, j);
                    H.at(n0, j) = q * H.at(n0, j) - p * tmp;
                }

                // column modification
                #pragma omp parallel for
                for (int i = 0; i <= n0; i++) {
                    element_t tmp = H.at(i, n0 - 1);
                    H.at(i, n0 - 1) = q * tmp + p * H.at(i, n0);
                    H.at(i, n0) = q * H.at(i, n0) - p * tmp;
                }

                // accumulate transformations
                #pragma omp parallel for
                for (int i = low; i <= high; i++) {
                    element_t tmp = V.at(i, n0 - 1);
                    V.at(i, n0 - 1) = q * tmp + p * V.at(i, n0);
                    V.at(i, n0) = q * V.at(i, n0) - p * tmp;
                }

            // complex pair
            } else {
                eigRe[n0 - 1] = x + p;
                eigRe[n0] = x + p;
                eigIm[n0 - 1] = z;
                eigIm[n0] = -z;
            }

            n0 = n0 - 2;
            iter = 0;

        // no convergence yet
        } else {

            // form shift
            x = H.at(n0, n0);
            y = 0.;
            w = 0.;
            if (l < n0) {
                y = H.at(n0 - 1, n0 - 1);
                w = H.at(n0, n0 - 1) * H.at(n0 - 1, n0);
            }

            // Wilkinson's original ad hoc shift
            if (iter == 10) {
                exshift += x;
                for (int i = low; i <= n0; i++) {
                    H.at(i, i) -= x;
                }
                s = std::abs(H.at(n0, n0 - 1)) + std::abs(H.at(n0 - 1, n0 - 2));
                x = y = 0.75 * s;
                w = -0.4375 * s * s;
            }

            // MATLAB's ad hoc shift
            if (iter == 30) {
                s = (y - x) / 2.;
                s = s * s + w;
                if (s > 0) {
                    s = std::sqrt(s);
                    if (y < x) {
                        s = -s;
                    }
                    s = x - w / ((y - x) / 2. + s);
                    for (int i = low; i <= n0; i++) {
                        H.at(i, i) -= s;
                    }
                    exshift += s;
                    x = y = w = 0.964;
                }
            }

            iter = iter + 1;

            // look for two consecutive small sub-diagonal elements
            int m = n0 - 2;
            while (m >= l) {
                z = H.at(m, m);
                r = x - z;
                s = y - z;
                p = (r * s - w) / H.at(m + 1, m) + H.at(m, m + 1);
                q = H.at(m + 1, m + 1) - z - r - s;
                r = H.at(m + 2, m + 1);
                s = std::abs(p) + std::abs(q) + std::abs(r);
                p = p / s;
                q = q / s;
                r = r / s;
                if (m == l) {
                    break;
                }
                if (std::abs(H.at(m, m - 1)) * (std::abs(q) + std::abs(r)) <
                    EPS * (std::abs(p) * (std::abs(H.at(m - 1, m - 1)) + std::abs(z) +
                        std::abs(H.at(m + 1, m + 1))))) {
                    break;
                }
                m--;
            }

            for (int i = m + 2; i <= n0; i++) {
                H.at(i, i - 2) = 0.;
                if (i > m + 2) {
                    H.at(i, i - 3) = 0.;
                }
            }

            // double QR step involving rows l:n and columns m:n
            for (int k = m; k <= n0 - 1; k++) {

                bool notlast = (k != n0 - 1);
                if (k != m) {
                    p = H.at(k, k - 1);
                    q = H.at(k + 1, k - 1);
                    r = (notlast ? H.at(k + 2, k - 1) : 0.);
                    x = std::abs(p) + std::abs(q) + std::abs(r);
                    if (isEq(x, 0.)) { continue; }
                    p = p / x;
                    q = q / x;
                    r = r / x;
                }

                s = std::sqrt(p * p + q * q + r * r);
                if (p < 0) {
                    s = -s;
                }
                if (!isEq(s, 0.)) {
                    if (k != m) {
                        H.at(k, k - 1) = -s * x;
                    } else if (l != m) {
                        H.at(k, k - 1) = -H.at(k, k - 1);
                    }
                    p = p + s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q = q / p;
                    r = r / p;

                    // row modification
                    for (int j = k; j < nn; j++) {
                        p = H.at(k, j) + q * H.at(k + 1, j);
                        if (notlast) {
                            p = p + r * H.at(k + 2, j);
                            H.at(k + 2, j) = H.at(k + 2, j) - p * z;
                        }
                        H.at(k, j) = H.at(k, j) - p * x;
                        H.at(k + 1, j) = H.at(k + 1, j) - p * y;
                    }

                    // column modification
                    for (int i = 0; i <= std::min(n0, k + 3); i++) {
                        p = x * H.at(i, k) + y * H.at(i, k + 1);
                        if (notlast) {
                            p += z * H.at(i, k + 2);
                            H.at(i, k + 2) = H.at(i, k + 2) - p * r;
                        }
                        H.at(i, k) = H.at(i, k) - p;
                        H.at(i, k + 1) = H.at(i, k + 1) - p * q;
                    }

                    // accumulate transformations
                    #pragma omp parallel for
                    for (int i = low; i <= high; i++) {
                        element_t tmp = x * V.at(i, k) + y * V.at(i, k + 1);
                        if (notlast) {
                            tmp += z * V.at(i, k + 2);
                            V.at(i, k + 2) = V.at(i, k + 2) - tmp * r;
                        }
                        V.at(i, k) = V.at(i, k) - tmp;
                        V.at(i, k + 1) = V.at(i, k + 1) - tmp * q;
                    }
                }  // (s != 0)
            }  // k loop
        }  // check convergence
    }  // while (n >= low)

    // back substitute to find vectors of upper triangular form
    if (isEq(norm, 0.)) { return; }

    for (n0 = nn - 1; n0 >= 0; n0--) {
        p = eigRe[n0];
        q = eigIm[n0];

        // real vector
        if (isEq(q, 0.)) {
            int l = n0;
            H.at(n0, n0) = 1.;
            for (int i = n0 - 1; i >= 0; i--) {
                w = H.at(i, i) - p;
                r = 0.;

                for (int j = l; j <= n0; j++) {
                    r += H.at(i, j) * H.at(j, n0);
                }
                if (eigIm[i] < 0.) {
                    z = w;
                    s = r;
                } else {
                    l = i;
                    if (isEq(eigIm[i], 0.)) {
                        if (!isEq(w, 0.)) {
                            H.at(i, n0) = -r / w;
                        } else {
                            H.at(i, n0) = -r / (EPS * norm);
                        }

                        // solve real equations
                    } else {
                        x = H.at(i, i + 1);
                        y = H.at(i + 1, i);
                        q = (eigRe[i] - p) * (eigRe[i] - p) + eigIm[i] * eigIm[i];
                        t = (x * s - z * r) / q;
                        H.at(i, n0) = t;
                        if (std::abs(x) > std::abs(z)) {
                            H.at(i + 1, n0) = (-r - w * t) / x;
                        } else {
                            H.at(i + 1, n0) = (-s - y * t) / z;
                        }
                    }

                    // overflow control
                    t = std::abs(H.at(i, n0));
                    if ((EPS * t) * t > 1) {
                        for (int j = i; j <= n0; j++) {
                            H.at(j, n0) = H.at(j, n0) / t;
                        }
                    }
                }
            }

        // complex vector
        } else if (q < 0) {

            int l = n0 - 1;
            // last vector component imaginary so matrix is triangular

            if (std::abs(H.at(n0, n0 - 1)) > std::abs(H.at(n0 - 1, n0))) {
                H.at(n0 - 1, n0 - 1) = q / H.at(n0, n0 - 1);
                H.at(n0 - 1, n0) = -(H.at(n0, n0) - p) / H.at(n0, n0 - 1);
            } else {
                std::pair<element_t, element_t> cdres =
                    cdiv(0., -H.at(n0 - 1, n0), H.at(n0 - 1, n0 - 1) - p, q);
                H.at(n0 - 1, n0 - 1) = cdres.first;
                H.at(n0 - 1, n0) = cdres.second;
            }
            H.at(n0, n0 - 1) = 0.;
            H.at(n0, n0) = 1.;
            for (int i = n0 - 2; i >= 0; i--) {
                element_t ra, sa, vr, vi;
                ra = 0.;
                sa = 0.;
                for (int j = l; j <= n0; j++) {
                    ra += H.at(i, j) * H.at(j, n0 - 1);
                    sa += H.at(i, j) * H.at(j, n0);
                }
                w = H.at(i, i) - p;

                if (eigIm[i] < 0.) {
                    z = w;
                    r = ra;
                    s = sa;
                } else {
                    l = i;
                    if (isEq(eigIm[i], 0.)) {
                        std::pair<element_t, element_t> cdres =
                            cdiv(-ra, -sa, w, q);
                        H.at(i, n0 - 1) = cdres.first;
                        H.at(i, n0) = cdres.second;
                    } else {
                        // Solve complex equations
                        x = H.at(i, i + 1);
                        y = H.at(i + 1, i);
                        vr = (eigRe[i] - p) * (eigRe[i] - p) + eigIm[i] * eigIm[i] - q * q;
                        vi = (eigRe[i] - p) * 2. * q;
                        if (isEq(vr, 0.) && isEq(vi, 0.)) {
                            vr = EPS * norm * (std::abs(w) + std::abs(q) +
                                std::abs(x) + std::abs(y) + std::abs(z));
                        }
                        std::pair<element_t, element_t> cdres =
                            cdiv(x * r - z * ra + q * sa,
                            x * s - z * sa - q * ra, vr, vi);
                        H.at(i, n0 - 1) = cdres.first;
                        H.at(i, n0) = cdres.second;
                        if (std::abs(x) > (std::abs(z) + std::abs(q))) {
                            H.at(i + 1, n0 - 1) = (-ra - w * H.at(i, n0 - 1) + q * H.at(i, n0)) / x;
                            H.at(i + 1, n0) = (-sa - w * H.at(i, n0) - q * H.at(i, n0 - 1)) / x;
                        } else {
                            std::pair<element_t, element_t> cdres =
                                cdiv(-r - y * H.at(i, n0 - 1), -s - y * H.at(i, n0), z, q);
                            H.at(i + 1, n0 - 1) = cdres.first;
                            H.at(i + 1, n0) = cdres.second;
                        }
                    }

                    // overflow control
                    t = std::max(std::abs(H.at(i, n0 - 1)), std::abs(H.at(i, n0)));
                    if ((EPS * t) * t > 1) {
                        for (int j = i; j <= n0; j++) {
                            H.at(j, n0 - 1) = H.at(j, n0 - 1) / t;
                            H.at(j, n0) = H.at(j, n0) / t;
                        }
                    }
                }
            }
        }
    }

    // vectors of isolated roots
    #pragma omp parallel for
    for (int i = 0; i < nn; i++) {
        if (i < low || i > high) {
            for (int j = i; j < nn; j++) {
                V.at(i, j) = H.at(i, j);
            }
        }
    }

    // back transformation to get eigenvectors of original matrix
    for (int j = nn - 1; j >= low; j--) {
        for (int i = low; i <= high; i++) {
            z = 0.;
            #pragma omp parallel for reduction(+:z)
            for (int k = low; k <= std::min(j, high); k++) {
                z += V.at(i, k) * H.at(k, j);
            }
            V.at(i, j) = z;
        }
    }
}
