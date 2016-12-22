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
    typedef long double element_t;

    /**
     * Matrix row type.
     */
    typedef std::vector<double> row_t;

    /**
     * Matrix data type.
     */
    typedef std::vector<std::vector<double> > data_t;

    /**
     * Construct a zero square matrix having a given order.
     *
     * @param n matrix order (non-negative)
     */
    SquareMatrix(std::size_t n = 0);

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
