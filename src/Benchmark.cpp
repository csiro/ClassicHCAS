// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#include <Rcpp.h>
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
    const Rcpp::NumericMatrix &ref_density, // cleaned reference density
    const Rcpp::NumericVector &xy_stats,    // mean(x), mean(y), sd(x), sd(y); ORDER MATTERS
    double xy_penalty = 0.0,                // penalising env nearest neighbour searching for geographic distance
    bool geographic = false,                // geographic/unprojected crs?
    double radius_km = 200,                 // radius in kilometers to consider ref points
    int k_env = 50,                         // number of ENV nn to select
    int k_rs = 20,                          // number of RS/Prob values to select
    double bin_width = 0.05,                // reference density bin width
    int bin_num = 400,                      // number of bins in reference density
    int offset = 0,                         // offset of reference density
    double confidence = 0.5,                // the LDC confidence index; default 0.5
    double lambda = 1.0,                    // distance-kernel bandwidth
    bool exclude_slef = true,               // whether to exclude a benchmark sample from assessing itself
    Rcpp::Nullable<Rcpp::NumericVector> temporal_weights = R_NilValue, // Gaussian weights for reference years
    bool make_su = false,                   // whether to produce SU map
    int num_threads = -1,                   // -1 or 0 utilises all available threads
    std::string kernel = "gaussian",        // distance weighting kernel
    Rcpp::Nullable<Rcpp::NumericVector> boost = R_NilValue) // optional maximum-probability weight multiplier
{
    double boost_factor = std::numeric_limits<double>::quiet_NaN();
    if (boost.isNotNull()) {
        Rcpp::NumericVector boost_value(boost);
        if (boost_value.size() != 1) {
            Rcpp::stop("'boost' must be NULL, NA, or one finite number greater than zero.");
        }
        if (!std::isnan(boost_value[0])) {
            if (!std::isfinite(boost_value[0]) || boost_value[0] <= 0.0) {
                Rcpp::stop("'boost' must be NULL, NA, or one finite number greater than zero.");
            }
            boost_factor = boost_value[0];
        }
    }
    if (xy_stats.size() != 4) {
        Rcpp::stop("'xy_stats' must contain exactly four values: mean(x), mean(y), sd(x), sd(y).");
    }
    if (k_env < 1) {
        Rcpp::stop("'k_env' must be >= 1.");
    }
    if (k_rs < 1) {
        Rcpp::stop("'k_rs' must be >= 1.");
    }
    if (bin_num < 2) {
        Rcpp::stop("'bin_num' must be >= 2.");
    }
    if (bin_width <= 0.0) {
        Rcpp::stop("'bin_width' must be > 0.");
    }
    if (!std::isfinite(confidence) || confidence < 0.0 || confidence > 1.0) {
        Rcpp::stop("'confidence' must be a finite number between 0 and 1.");
    }
    if (!std::isfinite(lambda) || lambda <= 0.0) {
        Rcpp::stop("'lambda' must be a finite number greater than zero.");
    }
    kernel = canonical_distance_kernel(kernel);
    if (kernel.empty()) {
        Rcpp::stop("'kernel' must be 'Gaussian'/'gaussian' or 'Cauchy'/'cauchy'.");
    }
    const DistanceKernel kernel_method =
        distance_kernel_from_string(kernel);

    // convert all Rcpp matrices to custom C++ matrix [faster computation and avoids OpenMp conflicts]
    RowMajorMatrix<float32_t> raster = as_Matrix<float32_t>(raster_vals);
    RowMajorMatrix<float32_t> samples = as_Matrix<float32_t>(sample_vals);
    // Keep ref_density as double; not much processing cost with table query for now...
    RowMajorMatrix<double> refdens = as_Matrix<double>(ref_density);
    // Get xy in double for calculating cos in degrees
    RowMajorMatrix<double> raster_xy = get_XY(raster_vals);

    // Ensure numeric value coming from R is float32_t
    const float32_t binwidth = static_cast<float32_t>(bin_width);

    const int nr = raster.rows();
    const int ns = samples.rows();
    if (raster.cols() < 4 || (raster.cols() - 2) % 2 != 0) {
        Rcpp::stop("'raster_vals' must contain x, y, predicted RS, and observed RS columns.");
    }
    const int nvar = (raster.cols() - 2) / 2; // number of RS vars
    int ndim = nvar + 2; // number of multi-variate space REM + XY

    const bool temporal = temporal_weights.isNotNull();
    std::vector<double> year_weights(1, 1.0);
    if (temporal) {
        Rcpp::NumericVector weights(temporal_weights);
        if (weights.size() < 1) {
            Rcpp::stop("'temporal_weights' must contain at least one value.");
        }

        year_weights.resize(weights.size());
        bool any_positive = false;
        for (R_xlen_t i = 0; i < weights.size(); ++i) {
            if (!std::isfinite(weights[i]) || weights[i] < 0.0) {
                Rcpp::stop("'temporal_weights' must contain finite, non-negative values.");
            }
            year_weights[i] = weights[i];
            any_positive = any_positive || weights[i] > 0.0;
        }
        if (!any_positive) {
            Rcpp::stop("'temporal_weights' must contain at least one positive value.");
        }
    }

    const int n_years = static_cast<int>(year_weights.size());
    const int obs_start = 2 + nvar;
    const int expected_sample_cols = 2 + nvar + nvar * n_years;
    if (samples.cols() != expected_sample_cols) {
        Rcpp::stop(
            temporal
                ? "'sample_vals' columns do not match the target features and temporal weights."
                : "'sample_vals' must contain x, y, predicted RS, and observed RS columns matching 'raster_vals'."
        );
    }

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
    if (xystats[2] == 0.0f || xystats[3] == 0.0f) {
        Rcpp::stop("'xy_stats' standard deviations (3rd and 4th elements) must be non-zero.");
    }

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
                // Calculate only 1 cos() for efficiency
                cos_scale = static_cast<int64_t>(std::cos(raster_xy(i, 1) * DEG_2_RAD) * 1000000);
            }

            // A combined search of raduis and k nearest samples to avoid data copying in memory
            std::vector<int> knn_env = combined_Search(
                sample_x, sample_y, samples, cell_rem,
                x, y, r2, k_env, ndim, cos_scale, geographic
            );

            // rs and env distance of nearest env neighbours
            std::vector<double> prdist_vect;
            std::vector<double> prob_vect;
            prdist_vect.reserve(knn_env.size());
            prob_vect.reserve(knn_env.size());

            // 'j' is the original index into the full 'samples' matrix
            for (const auto& j : knn_env)
            {
                // The xy coordinates should be ignored for REM dist, so only middleCols(2, nvar)
                const auto sub_rem_row = samples.row(j).middleCols(2, nvar);
                float32_t prdist = (cell_rem.rightCols(nvar) - sub_rem_row).template lpNorm<1>();

                if (exclude_slef && prdist < binwidth)
                {
                    continue;
                }

                double selected_prob = 0.0;
                double selected_score = -std::numeric_limits<double>::infinity();

                for (int year = 0; year < n_years; ++year)
                {
                    const auto sub_obs_row = samples.row(j).segment(
                        obs_start + year * nvar,
                        nvar
                    );
                    float32_t rsdist = (cell_obs - sub_obs_row).template lpNorm<1>();
                    double raw_prob = get_prob_value(
                        refdens,
                        prdist,
                        rsdist,
                        binwidth,
                        bin_num,
                        offset
                    );
                    double weighted_score = raw_prob * year_weights[year];

                    if (weighted_score > selected_score) {
                        selected_score = weighted_score;
                        selected_prob = weighted_score;
                    }
                }

                prdist_vect.push_back(static_cast<double>(prdist));
                prob_vect.push_back(selected_prob);
            }

            // Not enough candidates for this cell (e.g., sparse/edge tiles).
            if (prob_vect.empty())
            {
                Condition nan_cond {
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN()
                };
                condition_vect[i] = nan_cond;
                continue;
            }

            // descending sort prob values for selecting the 20 highest values
            // the probability vector itself will also be sorted
            std::vector<int> sorted_index = qsort_index(prob_vect, true);

            const int n_keep = std::min(k_rs, static_cast<int>(sorted_index.size()));
            std::vector<double> pr_dist(n_keep);     // ENV distances
            std::vector<double> prob_sorted(n_keep); // probability values

            // get the values of the available nearest RS neighbors
            for (int k = 0; k < n_keep; k++)
            {
                int id = sorted_index[k];
                pr_dist[k] = prdist_vect[id];
                prob_sorted[k] = prob_vect[k]; // prob_vect is already sorted by qsort_index; just get first 20
            }

            // Calculate the selected distance-kernel-weighted condition.
            Condition wcond = get_Condition(
                prob_sorted,
                pr_dist,
                confidence,
                lambda,
                kernel_method,
                boost_factor
            );
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
            out_mat(i, 1) = std::isfinite(cval.log_su) ? cval.log_su : NA_REAL;
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

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
