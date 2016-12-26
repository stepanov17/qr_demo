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
     * Matrix data type.  Use 1D inner matrix data representation
     * to speed up the calculations.
     */
    typedef std::vector<element_t> data_t;

    /**
     * Get a reference to (i, j)-element of the matrix.
     *
     * @return the reference to (i, j)-element
     */
    element_t &at(int i, int j) { return data[i * n + j]; }
    
    /**
     * Get a constant reference to (i, j)-element of the matrix.
     *
     * @return the constant reference to (i, j)-element
     */
    const element_t &at(int i, int j) const { return data[i * n + j]; }

    /**
     * Construct a zero square matrix having a given order.
     *
     * @param n matrix order (non-negative)
     */
    SquareMatrix(std::size_t n = 0);

    /**
     * Get the matrix order.
     *
     * @return the matrix order
     */
    int getN() const { return n; }

    /**
     * Check if the matrix is symmetric.
     *
     * @return true if the matrix is symmetric; otherwise false
     */
    bool isSymmetric() const;

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
