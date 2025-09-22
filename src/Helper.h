#pragma once

#define PI_Sq 9.869604401089358

// a struct to hold both condition and SU values
struct condition
{
    double hc;
    double su;
};


Rcpp::NumericMatrix get_XY(const Rcpp::NumericMatrix& X) {
    int n = X.nrow();
    if (X.ncol() < 2) Rcpp::stop("Input must have at least 2 columns");

    Rcpp::NumericMatrix out(n, 2);
    for (int i = 0; i < n; ++i) {
        out(i, 0) = X(i, 0);
        out(i, 1) = X(i, 1);
    }
    return out;
}


// Compute indices of k closest rows in L1 distance
template <typename T, typename QueryType>
std::vector<int> KNN_Search(const RowMajorMatrix<T> &X, const QueryType &q, int k) {
    const int n = X.rows();
    std::vector<std::pair<double,int>> dist_idx;
    dist_idx.reserve(n);

    for (int i = 0; i < n; ++i) {
        auto row_i = X.row(i);
        // Cast both to double for consistent arithmetic and avoid precision issues
        double dist = (row_i.template cast<double>() - q.template cast<double>()).template lpNorm<1>();
        dist_idx.emplace_back(dist, i);
    }

    if (k < n) {
        std::nth_element(dist_idx.begin(), dist_idx.begin() + k, dist_idx.end());
        dist_idx.resize(k);
    }

    std::sort(dist_idx.begin(), dist_idx.end());

    std::vector<int> result(dist_idx.size());
    for (size_t i = 0; i < dist_idx.size(); ++i) {
        result[i] = dist_idx[i].second;
    }
    return result;
}


// HCAS Cauchy weighted condition using histogram values and env distances
inline condition get_condition(
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
inline double get_prob(const RowMajorMatrix<double>& histo,
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

