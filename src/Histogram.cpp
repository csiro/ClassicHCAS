// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <cstdint>
#include "Matrix.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// calculate squared distance using xy; keep distance as double
inline int64_t distance(
    const int64_t x1, const int64_t y1, const int64_t x2, const int64_t y2,
    int64_t cos_scale, bool geo)
{
    const int64_t dx = x2 - x1;
    const int64_t dy = y2 - y1;

    const int64_t dy2 = dy * dy;

    int64_t dist2;
    if (geo) {
        // Scale longitude for geographic coordinates
        const int64_t dx_scaled = (dx * cos_scale) / 1000000;
        dist2 = dy2 + dx_scaled * dx_scaled;
    } else {
        dist2 = dy2 + dx * dx;
    }

    return dist2;
}


// [[Rcpp::export]]
Rcpp::IntegerMatrix histo_cpp(
    const Rcpp::NumericMatrix &rs_vals,
    const Rcpp::NumericMatrix &pr_vals,
    const Rcpp::NumericMatrix &xy,
    const double within_km = 1000.0, // radius in kilometers to consider ref points
    const double bin_width = 0.05,
    const int bin_num = 650,
    bool geographic = false,
    int num_threads = -1)
{
    // convert all Rcpp matrices to custom C++ matrix [faster computation and avoids OpenMp conflicts]
    RowMajorMatrix<float> rs = as_Matrix<float>(rs_vals);
    RowMajorMatrix<float> pr = as_Matrix<float>(pr_vals);
    // define the output matrix;
    RowMajorMatrix<uint64_t> histo_matrix(bin_num, bin_num);

    const int nr = rs.rows();
    const int ns = xy.nrow();

    // Pre-compute sample coordinates as integers
    std::vector<int64_t> sample_x(ns), sample_y(ns);

    double scale;
    int64_t squared_dist;
    const double radius_m = within_km * 1000.0;
    
    if (geographic) {
        // Geographic coordinates: use micro-degrees
        scale = 1000000.0;
        const double r_deg = radius_m / 111320.0; // METERS_PER_DEG_LAT
        const int64_t r_micro = static_cast<int64_t>(r_deg * scale);
        squared_dist = r_micro * r_micro;
    } else {
        // Projected coordinates: coordinates already in meters, scale for precision
        scale = 100.0;  // 0.01 meter precision
        const int64_t r_scaled = static_cast<int64_t>(radius_m * scale);
        squared_dist = r_scaled * r_scaled;
    }
    
    // Convert sample coordinates to integers
    for (int i = 0; i < ns; ++i) {
        sample_x[i] = static_cast<int64_t>(xy(i, 0) * scale);
        sample_y[i] = static_cast<int64_t>(xy(i, 1) * scale);
    }

    // the inverse of bin-width
    const float bwi = 1.0 / static_cast<float>(bin_width);
    // set the number of threads
    #ifdef _OPENMP
        if (num_threads < 1) num_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);
    #endif

    int64_t cos_scale = 0;
    // iterate over each row and compare it with all other rows 
    // dynamic schedule to balance the load
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nr; i++)
    {
        // get the XY from samples for the first sample
        const auto xi = sample_x[i];
        const auto yi = sample_y[i];
        if (geographic) {
            // Geographic: apply cosine scaling for longitude; only on i sample
            const double ya = static_cast<double>(yi);
            cos_scale = static_cast<int64_t>(std::cos(ya * 0.01745329) * 1000000);
        }
        // second iteration over rows skipping double counting, (i,j) is equal to (j,i)
        for (int j = i + 1; j < nr; j++)
        {
            // get the XY from samples for the second sample
            const auto xj = sample_x[j];
            const auto yj = sample_y[j];

            // // calculate geographic squared distance
            int64_t dist = distance(xi, yi, xj, yj, cos_scale, geographic);
            // only calculate for point within e.g. 1000km distance
            if (dist < squared_dist)
            {
                // Vectorised L1 distance over all columns
                float rsdist = (rs.row(i) - rs.row(j)).template lpNorm<1>();
                float prdist = (pr.row(i) - pr.row(j)).template lpNorm<1>();

                // now update the hist with the distances
                // calculate the row and column; must be floor to correctly get bin number
                int jj = std::floor(rsdist * bwi); // they are both float
                int ii = std::floor(prdist * bwi);
                // for now, ignore the values overshooting; they'll be mostly noise
                if (ii < bin_num && jj < bin_num) {
                    #pragma omp atomic
                    histo_matrix(ii, jj)++;
                }
            }
        }
    }

    // now transpose the matrix; then revert rows to save it in the "correct" position
    Rcpp::IntegerMatrix matrix(bin_num, bin_num);
    for (int i = 0; i < bin_num; ++i)
    {
        int ix = bin_num - 1 - i; // compute the reversed row index directly
        for (int j = 0; j < bin_num; ++j)
        {
            matrix(i, j) = histo_matrix(j, ix); // transpose and reverse in one step
        }
    }

    return matrix;
}

