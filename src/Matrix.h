#ifndef MATRIX_T
#define MATRIX_T

#include "RcppEigenQuiet.h"

// Generic row-major matrix type
template <typename T>
using RowMajorMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// Convert an Rcpp NumericMatrix to row-major Eigen matrix (generic on T)
template <typename T>
RowMajorMatrix<T> as_Matrix(const Rcpp::NumericMatrix &X) {
    // Step 1: get Eigen double matrix from Rcpp::NumericMatrix
    Eigen::MatrixXd tmp = Rcpp::as<Eigen::MatrixXd>(X);  // always double
    // Step 2: convert to desired type T in a row-major layout
    RowMajorMatrix<T> mat = tmp.template cast<T>(); // triggers copy/reorder to row-major
    return mat;
}


// Get a copy of XY in double to avoid losing percision in distance calc
inline RowMajorMatrix<double> get_XY(const Rcpp::NumericMatrix& X) {
    int n = X.nrow();
    if (X.ncol() < 2) Rcpp::stop("Input must have at least 2 columns");

    RowMajorMatrix<double> out(n, 2);
    for (int i = 0; i < n; ++i) {
        out(i, 0) = X(i, 0);
        out(i, 1) = X(i, 1);
    }
    return out;
}


// Filter a matrix view
template <typename Derived>
auto filter_Matrix(const Eigen::MatrixBase<Derived> &matrix,
                   const std::vector<int> &indices) {
    using T = typename Derived::Scalar;
    const int new_rows = indices.size();
    const int new_cols = matrix.cols();
    RowMajorMatrix<T> filtered(new_rows, new_cols);

    for (int i = 0; i < new_rows; ++i) {
        // get the whole row at once (contiguous in row-major!)
        filtered.row(i) = matrix.row(indices[i]);
    }

    return filtered;
}

#endif /* MATRIX_T */
