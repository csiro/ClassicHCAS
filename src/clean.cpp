#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix clean(NumericMatrix x, double bin_width, int offset = 0) {
    // define Gaussian filter weights in a vector
    std::vector<double> gaussian_weight = {
        1, 4,  7,  4,  1,
        4, 16, 26, 16, 4,
        7, 26, 41, 26, 7,
        4, 16, 26, 16, 4,
        1, 4,  7,  4,  1
    };

    // initialize d_gaussian_wt as a 5x5 matrix
    int size = 5;
    NumericMatrix d_gaussian_wt(size, size);
    int index = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            d_gaussian_wt(i, j) = gaussian_weight[index++];
        }
    }

    // Dimensions of the original and trimmed grid
    int orig_rows = 650;
    int orig_cols = 650;
    int n_rows = 400;
    int n_cols = 400;

    // Initialize d_new_grid
    NumericMatrix d_new_grid(n_rows, n_cols);

    Rcpp::NumericMatrix d_in_grid(orig_rows, orig_cols);
    for (int i = 0; i < orig_rows; ++i)
    {
        NumericVector row = x(i, _);
        std::reverse(row.begin(), row.end());
        d_in_grid(i, _) = row;
    }

    // Special treatment for d_in_grid[0][0]
    NumericMatrix d_in_grid_mod = clone(d_in_grid);
    d_in_grid_mod(0, 0) = 0;
    for (int col = 0; col < orig_cols; ++col) {
        if (d_in_grid_mod(647, col) > 0) {
            d_in_grid_mod(648, col) = d_in_grid_mod(2, col) / 2;
        } else {
            d_in_grid_mod(648, col) = 0;
        }
    }

    // perform Gaussian averaging
    double d_max = 0;
    double d_grid_max = 0;
    for (int i_row = 0; i_row < n_rows; ++i_row) {
        for (int i_col = 0; i_col < n_cols; ++i_col) {
            double d_sum_val = 0.;
            double d_sum_wt = 0.;
            for (int j_row = -2; j_row <= 2; ++j_row) {
                for (int j_col = -2; j_col <= 2; ++j_col) {
                    int row_idx = i_row + j_row + offset;
                    int col_idx = i_col + j_col + offset;
                    if (row_idx >= 0 && row_idx < orig_rows && col_idx >= 0 && col_idx < orig_cols) {
                        d_sum_val += d_in_grid_mod(row_idx, col_idx);
                        d_sum_wt += d_gaussian_wt(j_row + 2, j_col + 2);
                    }
                }
            }
            d_sum_wt = (273 + d_sum_wt) / 2;
            d_new_grid(i_row, i_col) = d_sum_val / d_sum_wt;

            // Update maximum value
            if (d_new_grid(i_row, i_col) > d_max) {
                d_max = d_new_grid(i_row, i_col);
            }

            // Apply noise removal at the bottom
            if (i_row > 395) {
                d_new_grid(i_row, i_col) = 0;
            }
        }
    }

    // Normalize columns and calculate median or mode
    // NumericVector i_med_row(n_cols);
    for (int i_col = 0; i_col < n_cols; ++i_col) {
        double f_sum = 0;
        for (int i_row = 0; i_row < n_rows; ++i_row) {
            if (d_new_grid(i_row, i_col) <= 3) {
                d_new_grid(i_row, i_col) = 0;
            }
            f_sum += d_new_grid(i_row, i_col);
        }

        // double f_med_sum = 0.;
        double d_grid_max = 0.;
        for (int i_row = 0; i_row < n_rows; ++i_row) {
            if (f_sum > 0) {
                double d_new_val = d_new_grid(i_row, i_col) / f_sum;
                d_new_grid(i_row, i_col) = d_new_val;
                if (d_new_val > d_grid_max) {
                    d_grid_max = d_new_val;
                }
                // if (f_med_sum <= 0.5) {
                //     f_med_sum += d_new_grid(i_row, i_col);
                //     i_med_row[i_col] = i_row;
                // }
            }
            if (d_new_grid(i_row, i_col) > d_max) {
                d_max = d_new_grid(i_row, i_col);
            }
        }
    }

    // Update the [0,0] point in the histogram
    d_new_grid(0, 0) = d_grid_max;


    return d_new_grid;
}
