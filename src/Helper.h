#pragma once

#define PI_Sq 9.869604401089358


// a struct to hold both condition and SU values
struct condition
{
    double hc;
    double su;
};


// scale x and y (columns 0 and 1), and add them to the end of the matrix
void add_xy_scaled(Lightweight_matrix<double> &matrix, const Rcpp::NumericVector &xy_stats, const double penalty)
{
    const int nr = matrix.nrow();
    const int nc = matrix.ncol();

    // resize the matrix to accommodate two new columns
    matrix.resize(nr, nc + 2);

    // scale and add columns; only if penalty is non-zero
    if (penalty != 0.0)
    {
        for (int i = 0; i < nr; i++)
        {
            matrix(i, nc) = (matrix(i, 0) - xy_stats[0]) / xy_stats[2] * penalty;       // x new column
            matrix(i, nc + 1) = (matrix(i, 1) - xy_stats[1]) / xy_stats[3] * penalty;   // y new column
        }
    }
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


// HCAS Cauchy weighted condition using histogram values and env distances
condition get_condition(
    const std::vector<double> &prob_values, // histo probability values
    double prob_max,                        // max probability value of the 20 records
    const std::vector<double> &pred_dists,  // the predicted/modelled distance
    const double confidence,                // the confidence value
    const double lambda)                    // the lambda of the Cauchy weighting
{
    const int n = prob_values.size();
    const double default_hc = -2.0;

    // calculate weights
    std::vector<double> weights(n);
    for (int i = 0; i < n; ++i)
    {
        double p_dist = pred_dists[i];
        if (p_dist > 0)
        {
            double den = PI_Sq * p_dist * lambda * (1.0 + (p_dist * p_dist) / (lambda * lambda));
            weights[i] = 1.0 / den;
        }
        else
        {
            weights[i] = 1.0;
        }
    }

    // calculate weighted mean and sum of weights
    double p_sum = 0.0;
    double w_sum = 0.0;
    for (int j = 0; j < n; ++j)
    {
        p_sum += prob_values[j] * weights[j];
        w_sum += weights[j];
    }

    // calculate p_mean and update prob_max if necessary
    double p_mean = (w_sum > 0) ? p_sum / w_sum : 0.0;
    if (p_mean > prob_max)
    {
        prob_max = p_mean;
    }

    // calculate hc
    double hc = default_hc;
    if (w_sum > 0)
    {
        hc = (prob_max * confidence) + (p_mean * (1.0 - confidence));
    }

    return {hc, w_sum};
}


// get the probability value from histogram
inline double get_prob(const Lightweight_matrix<double>& histo,
                       const double dist_pre,
                       const double dist_obs,
                       const double w_bin,
                       const int n_bin,
                       const int offset)
{
    static const int max_bin = n_bin - 1;
    static const double inverse_w_bin = 1.0 / w_bin;

    // static_cast over std::floor??
    int ii = std::min(static_cast<int>(dist_pre * inverse_w_bin), max_bin);
    int jj = std::min(static_cast<int>(dist_obs * inverse_w_bin), max_bin);
    // make sure there won't be negative values of i and j
    ii = std::max(ii - offset, 0);
    jj = std::max(jj - offset, 0);

    return histo(ii, jj);
}



// find the points within specified geographic distance
// NOTE: the kdtree was several minutes faster than (some version of) this
std::vector<int> nn_search(const Lightweight_matrix<double>& dat, const double x, const double y, const double sq_dist)
{
    std::vector<int> out;
    out.reserve(dat.nrow()); // reserve memory to avoid frequent re-allocations

    int n = dat.nrow();
    for (int i = 0; i < n; i++)
    {
        double dx = dat(i, 0) - x;
        double dy = dat(i, 1) - y;
        double dd = (dx * dx) + (dy * dy);

        if (dd < sq_dist)
        {
            out.push_back(i);
        }
    }

    return out;
}

