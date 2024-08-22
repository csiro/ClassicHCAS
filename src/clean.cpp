#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::NumericMatrix clean_cpp(
        const Rcpp::NumericMatrix &x,
        const int trim_size = 400,
        const int offset = 0)
{
    // get input dimensions
    const int nr = x.nrow();
    const int nc = x.ncol();

    // define Gaussian filter weights
    std::vector<double> weights = {
        1, 4,  7,  4,  1,
        4, 16, 26, 16, 4,
        7, 26, 41, 26, 7,
        4, 16, 26, 16, 4,
        1, 4,  7,  4,  1
    };

    // initialize Gaussian weight matrix
    const int size = 5;
    Rcpp::NumericMatrix weights_matrix(size, size);
    int id = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            weights_matrix(i, j) = weights[id++];
        }
    }

    // dimensions of the original and trimmed grid
    const int n_rows = trim_size;
    const int n_cols = trim_size;

    // initialize new_mat
    Rcpp::NumericMatrix new_mat(n_rows, n_cols);

    // reverse column for input grid
    Rcpp::NumericMatrix x_reverse(nr, nc);
    for (int i = 0; i < nr; ++i) {
        Rcpp::NumericVector column = x(_, i);
        std::reverse(column.begin(), column.end());
        x_reverse(_, i) = column;
    }

    // special treatment for x_reverse[0, 0]
    x_reverse(0, 0) = 0;

    // perform Gaussian averaging
    for (int i_row = 0; i_row < n_rows; i_row++) {
        for (int i_col = 0; i_col < n_cols; ++i_col) {
            double d_sum_val = 0.;
            double d_sum_wt = 0.;
            for (int j_row = -2; j_row <= 2; ++j_row) {
                for (int j_col = -2; j_col <= 2; ++j_col) {
                    int row_idx = i_row + j_row + offset;
                    int col_idx = i_col + j_col + offset;
                    if ((row_idx >= 0 && row_idx < nr) && (col_idx >= 0 && col_idx < nc)) {
                        d_sum_val += x_reverse(row_idx, col_idx);
                        d_sum_wt += weights_matrix(j_row + 2, j_col + 2);
                    }
                }
            }
            d_sum_wt = (273 + d_sum_wt) / 2; // 273 is sum of weights
            new_mat(i_row, i_col) = d_sum_val / d_sum_wt;

            // apply noise removal at the bottom
            if (i_row > 395) {
                new_mat(i_row, i_col) = 0;
            }
        }
    }

    double global_max = 0.;
    // normalize columns and calculate median
    for (int j = 0; j < n_cols; ++j) {
        double f_sum = 0.;
        // calculate sum and process values
        for (int i = 0; i < n_rows; i++) {
            if (new_mat(i, j) <= 3) {
                new_mat(i, j) = 0;
            }
            f_sum += new_mat(i, j);
        }

        // normalize values and find median
        for (int i = 0; i < n_rows; ++i) {
            if (f_sum > 0) {
                double new_val = new_mat(i, j) / f_sum;
                new_mat(i, j) = new_val;
                if (new_val > global_max) {
                    global_max = new_val;
                }
            }
        }
    }

    // update the [0,0] point in the histogram
    new_mat(0, 0) = global_max;

    // transpose the matrix to get x as RS and y as ENV
    Rcpp::NumericMatrix trans_mat(n_cols, n_rows);
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            trans_mat(j, i) = new_mat(i, j);
        }
    }

    return trans_mat;
}

