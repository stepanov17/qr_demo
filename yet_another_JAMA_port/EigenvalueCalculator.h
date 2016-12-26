
#ifndef EIGENVALUECALCULATOR_H
#define EIGENVALUECALCULATOR_H

#include "SquareMatrix.h"

/**
 * @class EigenvalueCalculator
 *
 * Yet another copy-paste of JAMA's EigenvalueDecomposition class.
 *
 */
class EigenvalueCalculator {

public:


    /**
     * A maximum matrix order allowed.
     */
    static const int MAXN = 10000;

    /**
     * Create the eigenvalue calculator for a given square matrix M.
     * The matrix order should not exceed 10000; otherwise
     * no calculations will be performed.
     *
     * @param M input square matrix
     */
    EigenvalueCalculator(const SquareMatrix &M);

    /**
     * Deleted copy constructor.
     *
     * @param orig dummy
     */
    EigenvalueCalculator(const EigenvalueCalculator& orig) = delete;

    /**
     * Get a pair of vectors containing real and imaginary eigenvalues parts.
     *
     * @return real and imaginary eigenvalues parts
     */
    std::pair<SquareMatrix::row_t, SquareMatrix::row_t> getEigenvalues();

    // TODO: return eigenvectors as well


private:

    typedef SquareMatrix::element_t element_t;


    // all the below procedures refer to
    // Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutines in EISPACK.

    /**
     * see JAMA's tred2()
     */
    void Tridiagonalize();
    /**
     * see JAMA's tql2()
     */
    void Diagonalize();
    /**
     * see JAMA's orthes()
     */
    void toHessenberg();
    /**
     * see JAMA's hqr2()
     */
    void toSchur();

    /**
     * matrix order
     */
    int n;

    /**
     * storage for the input matrix
     */
    SquareMatrix mM;

    /**
     * internal storage for eigenvectors
     */
    SquareMatrix mV;

    /**
     * internal storage for a nonsymmetric Hessenberg form
     */
    SquareMatrix mH;


    /**
     * real and imaginary parts of eigenvalues
     */
    SquareMatrix::row_t eigRe, eigIm;

    const element_t EPS = 2.2e-16; // magic JAMA's constant "2^(-52)"

    // some utilities

    /**
     * Check equality of two element_t values.
     *
     * @param x 1st value
     * @param y 2nd value
     * @return  bool (equal or not)
     */
    static inline bool isEq(element_t x, element_t y);

    /**
     * Overflow/underflow safe
     * \code
     * hypot(a, b) = sqrt(a * a + b * b).
     * \endcode
     *
     * @param a 1st argument
     * @param b 2nd argument
     * @return the hypot
     */
    static element_t hypot(element_t a, element_t b);

    /**
     * Overflow/underflow safe complex division (xr + i * xi) / (yr + i * yi).
     *
     * @param xr real part of the 1st argument
     * @param xi imag part of the 1st argument
     * @param yr real part of the 2nd argument
     * @param yi imag part of the 2nd argument
     * @return
     */
    static std::pair<element_t, element_t>
        cdiv(element_t xr, element_t xi, element_t yr, element_t yi);
};

#endif /* EIGENVALUECALCULATOR_H */

