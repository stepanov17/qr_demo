#ifndef QRITERATIONS_H
#define QRITERATIONS_H

#include "SquareMatrix.h"

/**
 * @class QRIterations
 *
 * Perform shifted QR algorithm for a given square matrix M.
 * Assuming a simple case when all eigenvalues of M are real.
 */
class QRIterations {
public:

    /**
     * QR shift type: Rayleigh, Wilkinson or none (zero)
     */
    enum shift { RAYLEIGH, WILKINSON, NONE };

    /**
     * A maximum order of the matrix to be processed.
     */
    static const int MAXN = 10000;

    /**
     * Bring matrix M to a Hessenberg form
     * and then run shifted QR iterations for it:
     * \code
     * (Q, R) := QR(M - shift * I)
     * M := R * Q + shift * I
     * shift := shift(M)
     * \endcode
     * (here I is an identity matrix) until
     * \code
     * max(abs(diag(M_prev) - diag(M))) <= precision.
     * \endcode
     * Initial shift is zero.
     *
     * @param M the matrix to be processed
     * @param shiftAlgorithm shift calculation algorithm;
     *                       "none" means "no shift applied".
     * @param precision precision
     * @param maxIter maximum number of iterations (default value is 1e7)
     * @throw std::invalid_argument if the precision value is non-positive
     * @throw std::invalid_argument if the matrix order exceeds MAXN
     */
    static void run(SquareMatrix            &M,
                    shift                    shiftAlgorithm,
                    SquareMatrix::element_t  precision,
                    unsigned long            maxIter = 10000000);
private:

    /**
     * Get shift for the matrix.
     * @param M  the matrix
     * @param shiftType shift type
     * @return the shift
     */
    static SquareMatrix::element_t getShift(const SquareMatrix &M,
                                            shift shiftType);

    /**
     * Apply shift to matrix:
     * \code
     * M := M + shift * I
     * \endcode
     *
     * @param M the matrix
     * @param shift the shift
     */
    static void applyShift(SquareMatrix &M, SquareMatrix::element_t shift);

    QRIterations();
};

#endif /* QRITERATIONS_H */
