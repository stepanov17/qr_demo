
#include "SquareMatrix.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include <stdexcept>

#include <cmath>
#include <limits>
#include <string>



SquareMatrix::SquareMatrix(std::size_t n_p) : n(n_p) {

    data.clear(); // just in case
    data.assign(n * n, 0.);
}

bool
SquareMatrix::isSymmetric() const {

    element_t eps = std::numeric_limits<SquareMatrix::element_t>::epsilon();

    bool isSymm = true;
    for (int j = 0; (j < n) && isSymm; ++j) {
        for (int i = 0; (i < n) && isSymm; ++i) {
            if (std::abs(at(i, j) - at(j, i)) > eps) { return false; }
        }
    }
    return true;
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
                throw std::range_error(
                    "max order " + std::to_string(MAXN) + " is exceeded");
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
