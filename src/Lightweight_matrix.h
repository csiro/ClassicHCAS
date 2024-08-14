#ifndef _LW_MATRIX_H
#define _LW_MATRIX_H

/**
 * Defining a light-weight matrix class and methods. This is mostly to
 * convert the Rcpp::NumericMatrix to std::vector (with the same operators) so it will
 * be more efficient and doesn't cause issues when using OpenMP for parallel processing.
 */

#include <vector>

// a light weight matrix class for faster calculation
template <typename T>
class Lightweight_matrix
{
private:
    using MatrixType = std::vector<T>;

public:
    explicit Lightweight_matrix(int rows, int columns) : rows_(rows), columns_(columns), matrix_(rows * columns) {}

    template <typename Matrix>
    explicit Lightweight_matrix(Matrix matrix) : rows_(matrix.nrow()), columns_(matrix.ncol()), matrix_(matrix.nrow() * matrix.ncol())
    {
        for (int i = 0; i < rows_; ++i)
        {
            for (int j = 0; j < columns_; ++j)
            {
                operator()(i, j) = static_cast<T>(matrix(i, j));
            }
        }
    }

    T operator()(int i, int j) const
    {
        return matrix_[i * columns_ + j];
    }

    T& operator()(int i, int j)
    {
        return matrix_[i * columns_ + j];
    }

    int ncol() const
    {
        return columns_;
    }

    int nrow() const
    {
        return rows_;
    }

    // resize the matrix by number of row and column
    void resize(int new_rows, int new_columns)
    {
        MatrixType new_matrix(new_rows * new_columns);
        int min_rows = std::min(rows_, new_rows);
        int min_columns = std::min(columns_, new_columns);

        for (int i = 0; i < min_rows; ++i)
        {
            for (int j = 0; j < min_columns; ++j)
            {
                new_matrix[i * new_columns + j] = operator()(i, j);
            }
        }

        matrix_ = std::move(new_matrix);
        rows_ = new_rows;
        columns_ = new_columns;
    }

private:
    int rows_;
    int columns_;
    MatrixType matrix_;
};


// convert to NumericMatrix
Rcpp::NumericMatrix as_NumericMatrix(const Lightweight_matrix<double> &matrix)
{
    int nr = matrix.nrow();
    int nc = matrix.ncol();

    Rcpp::NumericMatrix rcpp_matrix(nr, nc);

    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            rcpp_matrix(i, j) = matrix(i, j);
        }
    }

    return rcpp_matrix;
}


// get a subset matrix using a vector of indices
Lightweight_matrix<double> filter_matrix(const Lightweight_matrix<double> &matrix,
                                         const std::vector<int> &indices)
{
    int new_rows = indices.size();
    int new_columns = matrix.ncol();
    Lightweight_matrix<double> filtered_matrix(new_rows, new_columns);

    for (int i = 0; i < new_rows; ++i)
    {
        int original_index = indices[i];
        for (int j = 0; j < new_columns; ++j)
        {
            filtered_matrix(i, j) = matrix(original_index, j);
        }
    }

    return filtered_matrix;
}

#endif // _LW_MATRIX_H
