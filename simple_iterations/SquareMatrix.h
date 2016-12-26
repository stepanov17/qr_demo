#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

#include <vector>

/**
 * @class SquareMatrix
 *
 * Square matrix class.
 */
class SquareMatrix {
public:

    /**
     *  Matrix element type.
     */
    typedef /*long*/ double element_t;

    /**
     * Matrix row type.
     */
    typedef std::vector<element_t> row_t;

    /**
     * Matrix data type.
     */
    typedef std::vector<std::vector<element_t> > data_t;

    /**
     * Construct a zero square matrix having a given dimension.
     *
     * @param n matrix order (non-negative)
     */
    SquareMatrix(std::size_t n = 0);

    /**
     * Get a constant reference to the matrix elements.
     *
     * @return the constant reference to the elements
     */
    const data_t &getData() const { return data; }

    /**
     * Get a reference to the matrix elements.
     *
     * @return the reference to the elements
     */
    data_t &accessData() { return data; }

    /**
     * Get the matrix order.
     *
     * @return the matrix order
     */
    int getN() const { return n; }

    /**
     * Matrix multiplication operation.
     *
     * @param M right-side multiplication operand
     * @return the multiplication result
     * @throw std::invalid_argument in case of inconsistent dimensions
     */
    SquareMatrix operator *(const SquareMatrix &M) const;

    /**
     * Assignment operator.
     *
     * @param M a matrix to be copied.
     * @return reference to *this
     * @throw std::invalid_argument in case of inconsistent dimensions
     */
    SquareMatrix& operator =(const SquareMatrix &M);

    /**
     * Return a main diagonal of the matrix.
     *
     * @return the main diagonal of the matrix
     */
    row_t diag() const;

    /**
     * Return a Hessenberg form of the matrix.
     *
     * @return the Hessenberg form
     */
    SquareMatrix Hessenberg() const;

    /**
     * Calculate QR decomposition of the matrix: M = Q * R, where
     * Q is an orthogonal matrix Q and R is an upper triangular matrix.
     *
     * @return the (Q, R) pair
     */
    std::pair<SquareMatrix, SquareMatrix> QR() const;

    /**
     * Print the matrix to a console.
     */
    void print() const;

    /**
     * Read a matrix from a file. It is supposed that the file contains
     * element values separated by spaces (a line per row).
     * No empty lines are allowed.
     * Matrix order must be in range [1, 10000]
     *
     * @param path the file path
     * @throw std::runtime_error in case of file reading issues (e.g., wrong path)
     * @throw std::range_error in case of invalid data ranges
     */
    void read(const char *path);

private:

    const static int MAXN = 10000;

    int n;
    data_t data;
};

#endif /* SQUAREMATRIX_H */
