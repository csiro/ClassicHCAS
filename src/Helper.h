#pragma once
#include <vector>
#include <cstdint>

// a struct to hold both condition and SU values
struct Condition 
{
    double hc;
    double su;
};


// Get a copy of XY in double to avoid losing percision in distance calc
RowMajorMatrix<double> get_XY(const Rcpp::NumericMatrix& X) {
    int n = X.nrow();
    if (X.ncol() < 2) Rcpp::stop("Input must have at least 2 columns");

    RowMajorMatrix<double> out(n, 2);
    for (int i = 0; i < n; ++i) {
        out(i, 0) = X(i, 0);
        out(i, 1) = X(i, 1);
    }
    return out;
}


// finding points in a radius with a fast approximation for lat-long distance
std::vector<int> radius_Search(
    const std::vector<int64_t>& x_data,
    const std::vector<int64_t>& y_data,
    const int64_t query_x,
    const int64_t query_y,
    int64_t r2,
    int64_t cos_scale,
    bool geographic
) {
    const int size = x_data.size();

    std::vector<int> within;
    within.reserve(size);
    
    for (int i = 0; i < size; ++i) {
        const int64_t dx = x_data[i] - query_x;
        const int64_t dy = y_data[i] - query_y;
        
        // Early latitude rejection
        const int64_t dy2 = dy * dy;
        if (dy2 > r2) continue;
        
        int64_t dist2;
        if (geographic) {
            // Scale longitude for geographic coordinates
            const int64_t dx_scaled = (dx * cos_scale) / 1000000;
            dist2 = dy2 + dx_scaled * dx_scaled;
        } else {
            dist2 = dy2 + dx * dx;
        }
        
        if (dist2 <= r2) {
            within.push_back(i);
        }
    }
    
    return within;
}


// Compute indices of k nearest rows in L1 distance
template <typename T, typename Q>
std::vector<int> KNN_Search(const RowMajorMatrix<T> &X, const Q &q, int k) {
    const int n = X.rows();
    std::vector<std::pair<T,int>> dist_idx;
    dist_idx.reserve(n);

    for (int i = 0; i < n; ++i) {
        auto row_i = X.row(i);
        // Cast both to double for consistent arithmetic and avoid precision issues
        T dist = (row_i - q).template lpNorm<1>();
        dist_idx.emplace_back(dist, i);
    }

    if (k < n) {
        std::nth_element(dist_idx.begin(), dist_idx.begin() + k, dist_idx.end());
        dist_idx.resize(k);
    }

    // Sorting the indiceis will help later in filtering
    std::sort(dist_idx.begin(), dist_idx.end());

    std::vector<int> result(dist_idx.size());
    for (size_t i = 0; i < dist_idx.size(); ++i) {
        result[i] = dist_idx[i].second;
    }
    return result;
}


// HCAS Cauchy weighted condition calculation using histogram values and env distances
inline Condition get_Condition(
    const std::vector<double> &prob_values, // histo probability values
    const std::vector<double> &pred_dists,  // the predicted/modelled distance
    double prob_max,                        // max probability value of the 20 records
    const double confidence,                // the confidence value
    const double lambda)                    // the lambda of the Cauchy weighting
{
    const int n = prob_values.size();

    const double PI_SQ = 9.869604401089358;
    const double DEFAULT_HC = -2.0;
    double lambda_sq = lambda * lambda;
    double w_sum = 0.0;
    double p_sum = 0.0;

    // calculate weights
    for (int i = 0; i < n; ++i)
    {
        double p_dist = pred_dists[i];

        double weight = 1.0;
        if (p_dist > 0) {
            weight = 1.0 / (PI_SQ * p_dist * lambda * (1.0 + (p_dist * p_dist) / lambda_sq));
        }

        p_sum += prob_values[i] * weight;
        w_sum += weight;
    }

    // calculate hc    
    double hc = DEFAULT_HC;
    if (w_sum > 0) {
        double p_mean = p_sum / w_sum;
        prob_max = std::max(prob_max, p_mean);
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

