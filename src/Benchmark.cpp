// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <limits>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif
// Local scripts
#include "Quickisort.h"
#include "Matrix.h"
#include "Helper.h"
#include "Float32_t.h" // importing float32_t

using namespace Rcpp;

static constexpr double DEG_2_RAD = M_PI / 180.0;


// Data type for all distance calcaulations use float for speed imporvement (inputs are normalised)
// [[Rcpp::export]]
Rcpp::NumericMatrix bench_cpp(
    const Rcpp::NumericMatrix &raster_vals, // a raster stack to read everything at once: x,y,rs,env
    const Rcpp::NumericMatrix &sample_vals, // extraction of values of raster using smaples xy: x,y,RS,ENV
    const Rcpp::NumericMatrix &histogram,   // cleaned histogram
    const Rcpp::NumericVector &xy_stats,    // mean(x), mean(y), sd(x), sd(y); ORDER MATTERS
    double xy_penalty = 0.0,          // penalising env nearest neighbour searching for geographic distance
    bool geographic = false,          // correction scale for CRS; 1 if metric, 100,000 otherwise
    double radius_km = 200,           // radius in kilometers to consider ref points
    int k_env = 50,                   // number of ENV nn to select
    int k_rs = 20,                    // number of RS/Prob values to select
    double bin_width = 0.05,          // histogram bin width
    int bin_num = 400,                // number of bins in histogram
    int offset = 0,                   // offset of histogram
    double confidence = 0.5,          // the LDC confidence index; default 0.5
    double lambda = 2.0,              // the lambda of the Cauchy weighting
    bool exclude_slef = true,         // whether to exclude a benchmark sample from assessing itself
    bool make_su = false,             // whether to produce SU map
    int num_threads = -1)                   // -1 or 0 utilises all available threads
{
    // convert all Rcpp matrices to custom C++ matrix [faster computation and avoids OpenMp conflicts]
    RowMajorMatrix<float32_t> raster = as_Matrix<float32_t>(raster_vals);
    RowMajorMatrix<float32_t> samples = as_Matrix<float32_t>(sample_vals);
    // Keep histo as double; not much processing cost with histo query for now...
    RowMajorMatrix<double> histo = as_Matrix<double>(histogram);
    // Get xy in double for calculating cos in degrees
    RowMajorMatrix<double> raster_xy = get_XY(raster_vals);

    // Ensure numeric value coming from R is float32_t
    const float32_t binwidth = static_cast<float32_t>(bin_width);

    const int nr = raster.rows();
    const int ns = samples.rows();
    int nvar = (samples.cols() - 2) / 2; // number of RS vars
    int ndim = nvar + 2; // number of multi-variate space REM + XY

    double scale;
    int64_t r2;
    const double radius_m = radius_km * 1000.0;
    if (geographic) {
        // Geographic coordinates: use micro-degrees; scale for precision 0.11 meter
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

    // Cast the values once; cleaner code
    const float32_t xypenalty = static_cast<float32_t>(xy_penalty);
    std::vector<float32_t> xystats(xy_stats.begin(), xy_stats.end());

    // Normalise XY and apply weight to it
    samples.col(0) = ((samples.col(0).array() - xystats[0]) / xystats[2]) * xypenalty;
    samples.col(1) = ((samples.col(1).array() - xystats[1]) / xystats[3]) * xypenalty;
    raster.col(0) = ((raster.col(0).array() - xystats[0]) / xystats[2]) * xypenalty;
    raster.col(1) = ((raster.col(1).array() - xystats[1]) / xystats[3]) * xypenalty;

    // output condition vector
    std::vector<Condition> condition_vect(nr);

    // set the number of threads
    #ifdef _OPENMP
        if (num_threads < 1) num_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);
    #endif

    // iterate over the raster cells (rows of matrix)
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nr; i++)
    {
        const auto cell_rem = raster.row(i).leftCols(ndim);
        const auto cell_obs = raster.row(i).rightCols(nvar);

        // Take care of NaN cells
        if ((cell_obs.array().isNaN()).any())
        {
            Condition nan_cond {
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()
            };
            // since i is unique for each thread, there's no need for #pragma omp critical
            condition_vect[i] = nan_cond;
        }
        else
        {
            // // get the target points xy from raster stack;
            const int64_t x = static_cast<int64_t>(raster_xy(i, 0) * scale);
            const int64_t y = static_cast<int64_t>(raster_xy(i, 1) * scale);

            int64_t cos_scale = 0;
            if (geographic) {
                // Calcualte only 1 cos() for efficiency
                cos_scale = static_cast<int64_t>(std::cos(raster_xy(i, 1) * DEG_2_RAD) * 1000000);
            }

            // A combined search of raduis and k nearest samples to avoid data copying in memory
            std::vector<int> knn_env = combined_Search(
                sample_x, sample_y, samples, cell_rem,
                x, y, r2, k_env, ndim, cos_scale, geographic
            );

            // rs and env distance of the 50 nearest env neighbour
            std::vector<double> prdist_vect(knn_env.size());
            std::vector<double> histo_vect(knn_env.size());

            int n = 0;
            // 'j' is the original index into the full 'samples' matrix
            for (const auto& j : knn_env)
            {
                // Get the OBS part of the sample row for RS distance calculation
                const auto sub_obs_row = samples.row(j).rightCols(nvar);
                // Vectorised L1 distance over all columns for RS (OBS)
                float32_t rsdist = (cell_obs - sub_obs_row).template lpNorm<1>();

                // The xy coordinates should be ignored for REM dist, so only middleCols(2, nvar)
                const auto sub_rem_row = samples.row(j).middleCols(2, nvar);
                float32_t prdist = (cell_rem.rightCols(nvar) - sub_rem_row).template lpNorm<1>();

                if (exclude_slef && prdist < binwidth)
                {
                    continue;
                }

                prdist_vect[n] =  static_cast<double>(prdist);
                histo_vect[n] = get_Prob(histo, prdist, rsdist, binwidth, bin_num, offset);

                n++;
            }

            // descending sort prob values for selecting the 20 highest values
            // the histo_vect itself will also be sorted
            std::vector<int> sorted_index = qsort_index(histo_vect, true);

            std::vector<double> pr_dist(k_rs);     // ENV distances
            std::vector<double> prob_sorted(k_rs); // probability values
            // get the values of 20 nearest RS neighbors
            for (int k = 0; k < k_rs; k++)
            {
                int id = sorted_index[k];
                pr_dist[k] = prdist_vect[id];
                prob_sorted[k] = histo_vect[k]; // histo_vect is already sorted by qsort_index; just get first 20
            }

            // calculate the Cauchy weighting condition
            Condition wcond = get_Condition(prob_sorted, pr_dist, prob_sorted[0], confidence, lambda);
            condition_vect[i] = wcond;
        }
    }

    // output matrix; doesn't occupy memory until here
    const int dim = make_su ? 2 : 1;
    Rcpp::NumericMatrix out_mat(nr, dim);

    // produce SU map if requested
    int i = 0;
    if (make_su)
    {
        // output matrix - both condition and SU values
        for (const auto& cval : condition_vect)
        {
            out_mat(i, 0) = cval.hc;
            out_mat(i, 1) = std::log(cval.su);
            i++;
        }
    }
    else
    {
        // output matrix - only condition values
        for (const auto& cval : condition_vect)
        {
            out_mat(i, 0) = cval.hc;
            i++;
        }
    }

    return out_mat;
}

