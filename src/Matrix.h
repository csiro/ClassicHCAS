#ifndef MATRIX_T
#define MATRIX_T
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <RcppEigen.h>

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


template <typename T>
RowMajorMatrix<T> filter_Matrix(const RowMajorMatrix<T> &matrix,
                                const std::vector<int> &indices) {
    const int new_rows = indices.size();
    const int new_cols = matrix.cols();
    RowMajorMatrix<T> filtered(new_rows, new_cols);

    for (int i = 0; i < new_rows; ++i) {
        int orig = indices[i];
        // copy whole row at once (contiguous in row-major!)
        filtered.row(i) = matrix.row(orig);
    }

    return filtered;
}

#endif /* MATRIX_T */