#pragma once
#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility> // for std::pair

// a struct to hold both condition and SU values
struct Condition
{
    double hc;
    double su;
};


// Combine Radius Search (XY) and KNN Search (ENV)
template <typename Q>
std::vector<int> combined_Search(
    const std::vector<int64_t>& sample_x, // Integer and scaled sample coordinates
    const std::vector<int64_t>& sample_y,
    const RowMajorMatrix<float>& samples, // The full original samples matrix
    const Eigen::MatrixBase<Q>& cell_rem, // The target REM row of the raster
    const int64_t query_x,
    const int64_t query_y,
    int64_t r2, // Integer and scaled squared radius
    int k_env, // K for KNN
    int ndim, // Number of REM (Environmental) columns
    int64_t cos_scale,
    bool geographic
) {
    const int size = sample_x.size();
    // Use pair<ENV_Distance, Original_Index>
    std::vector<std::pair<float, int>> dist_idx;
    // Reserve 10000; it's a heuristic, but better than nothing
    dist_idx.reserve(10000);

    // 1. Perform radius check, and calculate ENV distance for KNN later
    for (int i = 0; i < size; ++i) {
        const int64_t dx = sample_x[i] - query_x;
        const int64_t dy = sample_y[i] - query_y;

        const int64_t dy2 = dy * dy;
        // Early latitude rejection
        if (dy2 > r2) continue;

        int64_t dist2;
        if (geographic) {
            const int64_t dx_scaled = (dx * cos_scale) / 1000000;
            dist2 = dy2 + dx_scaled * dx_scaled;
        } else {
            dist2 = dy2 + dx * dx;
        }

        // 2. If within radius, calculate ENV L1 Distance
        if (dist2 <= r2) {
            // Get the REM part of the sample row (first ndim columns)
            // This is for KNN and contain both XY and REM
            const auto sample_rem = samples.row(i).leftCols(ndim);

            // Calculate L1 distance (Minkowski p-norm with p=1)
            float env_dist = (sample_rem - cell_rem).template lpNorm<1>();

            // Store the ENV distance and the *original* sample index (i)
            dist_idx.emplace_back(env_dist, i);
        }
    }

    // 3. Find the K nearest using std::nth_element
    const int n_found = dist_idx.size();
    // If we found more than k_env samples, use nth_element to find the k smallest
    if (k_env < n_found) {
        // Find the k_env'th smallest element (0-indexed)
        std::nth_element(dist_idx.begin(), dist_idx.begin() + k_env, dist_idx.end());
        dist_idx.resize(k_env);
    }

    // 4. Extract final indices
    std::vector<int> result(dist_idx.size());
    // result[i] is the original index 'i' into the full 'samples' matrix
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
inline double get_Prob(const RowMajorMatrix<double>& histo,
                       const float dist_pre,
                       const float dist_obs,
                       const float w_bin,
                       const int n_bin,
                       const int offset)
{
    static const int max_bin = n_bin - 1;
    static const float inverse_w_bin = 1.0 / w_bin;

    // static_cast over std::floor??
    int ii = std::min(static_cast<int>(dist_pre * inverse_w_bin), max_bin);
    int jj = std::min(static_cast<int>(dist_obs * inverse_w_bin), max_bin);
    // make sure there won't be negative values of i and j
    ii = std::max(ii - offset, 0);
    jj = std::max(jj - offset, 0);

    return histo(ii, jj);
}

