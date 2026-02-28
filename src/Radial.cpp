// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <cstdint>
#include "Matrix.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static constexpr double DEG_2_RAD = M_PI / 180.0;

// Fast integer-based distance calculation for geographic coordinates
inline int radius_Count(
    const std::vector<int64_t>& x_data,
    const std::vector<int64_t>& y_data,
    int64_t query_x,
    int64_t query_y,
    int64_t r2,
    int64_t cos_scale,
    bool geo
) {
    int count = 0;
    const int size = x_data.size();
    
    for (int i = 0; i < size; ++i) {
        const int64_t dx = x_data[i] - query_x;
        const int64_t dy = y_data[i] - query_y;
        
        // Early latitude rejection
        const int64_t dy2 = dy * dy;
        if (dy2 > r2) continue;
        
        int64_t dist2;
        if (geo) {
            // Scale longitude for geographic coordinates
            const int64_t dx_scaled = (dx * cos_scale) / 1000000;
            dist2 = dy2 + dx_scaled * dx_scaled;
        } else {
            const int64_t dx2 = dx * dx;
            dist2 = dy2 + dx2;
        }
        
        if (dist2 <= r2) count++;
    }
    
    return count;
}


// [[Rcpp::export]]
Rcpp::IntegerVector radial_count_cpp(
    const Rcpp::NumericMatrix &rast,
    const Rcpp::NumericMatrix &xy,
    const double radius_km = 200.0,
    bool geographic = false,
    int num_threads = -1
) {
    RowMajorMatrix<double> raster = as_Matrix<double>(rast);
    RowMajorMatrix<double> samples = as_Matrix<double>(xy);
    
    const int nr = raster.rows();
    const int ns = samples.rows();
    std::vector<int> out(nr);
        
    double scale;
    int64_t r2;
    const double radius_m = radius_km * 1000.0;
    
    if (geographic) {
        // Geographic coordinates: use micro-degrees
        scale = 1000000.0;
        const double r_deg = radius_m / 111320.0; // METERS_PER_DEG_LAT
        const int64_t r_micro = static_cast<int64_t>(r_deg * scale);
        r2 = r_micro * r_micro;
    } else {
        // Projected coordinates: coordinates already in meters, scale for precision
        scale = 100.0;  // 0.01 meter precision
        const int64_t r_scaled = static_cast<int64_t>(radius_m * scale);
        r2 = r_scaled * r_scaled;
    }
    
    // Pre-compute sample coordinates as integers
    std::vector<int64_t> sample_x(ns), sample_y(ns);
    // Convert sample coordinates to integers
    for (int i = 0; i < ns; ++i) {
        sample_x[i] = static_cast<int64_t>(samples(i, 0) * scale);
        sample_y[i] = static_cast<int64_t>(samples(i, 1) * scale);
    }
    
    #ifdef _OPENMP
        if (num_threads < 1) num_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);
    #endif
    
    // Dynamic scheduling improves load balancing when NA-heavy regions
    // make per-cell work uneven across contiguous raster chunks.
    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < nr; i++) {
        const auto row = raster.row(i);
        
        if ((row.array().isNaN()).any()) {
            out[i] = NA_INTEGER;
            continue;
        }
        
        const int64_t query_x = static_cast<int64_t>(row(0) * scale);
        const int64_t query_y = static_cast<int64_t>(row(1) * scale);
        
        int64_t cos_scale = 0;
        bool geo = false;
        if (geographic) {
            cos_scale = static_cast<int64_t>(std::cos(row(1) * DEG_2_RAD) * 1000000);
            geo = true;
        }
        // Find number of samples within the radius
        int count = radius_Count(sample_x, sample_y, query_x, query_y, r2, cos_scale, geo);

        out[i] = count;
    }
    
    return Rcpp::wrap(out);
}
