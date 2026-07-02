// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Float32_t.h"
#include "Matrix.h"
#include "Helper.h"

using namespace Rcpp;

namespace {

static constexpr double DEG_2_RAD = M_PI / 180.0;

int probability_partition(
    std::vector<double>& values,
    std::vector<int>& order,
    int low,
    int high)
{
    const double pivot = values[high];
    int i = low - 1;
    for (int j = low; j < high; ++j) {
        if (values[j] >= pivot) {
            i += 1;
            std::swap(values[i], values[j]);
            std::swap(order[i], order[j]);
        }
    }
    std::swap(values[i + 1], values[high]);
    std::swap(order[i + 1], order[high]);
    return i + 1;
}

void probability_quicksort(
    std::vector<double>& values,
    std::vector<int>& order,
    int low,
    int high)
{
    if (low >= high) {
        return;
    }

    const int pivot = probability_partition(values, order, low, high);
    probability_quicksort(values, order, low, pivot - 1);
    probability_quicksort(values, order, pivot + 1, high);
}

std::vector<int> descending_probability_order(std::vector<double>& values)
{
    std::vector<int> order(values.size());
    for (size_t i = 0; i < order.size(); ++i) {
        order[i] = static_cast<int>(i);
    }
    if (!values.empty()) {
        probability_quicksort(
            values,
            order,
            0,
            static_cast<int>(values.size()) - 1
        );
    }
    return order;
}

} // namespace


// [[Rcpp::export]]
Rcpp::List reference_use_cpp(
    const Rcpp::NumericMatrix &target_vals,
    const Rcpp::NumericMatrix &sample_vals,
    const Rcpp::NumericMatrix &ref_density,
    const Rcpp::NumericVector &xy_stats,
    double xy_penalty = 0.0,
    bool geographic = false,
    double radius_km = 200,
    int k_env = 50,
    int k_rs = 20,
    double bin_width = 0.05,
    int bin_num = 400,
    int offset = 0,
    double confidence = 0.5,
    double lambda = 1.0,
    bool exclude_slef = true,
    int num_threads = -1,
    bool weighted_max = false,
    std::string kernel = "gaussian",
    Rcpp::Nullable<Rcpp::NumericVector> boost = R_NilValue)
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
    if (k_rs < 1 || k_rs > k_env) {
        Rcpp::stop("'k_rs' must be between 1 and 'k_env'.");
    }
    if (bin_num < 2) {
        Rcpp::stop("'bin_num' must be >= 2.");
    }
    if (bin_width <= 0.0) {
        Rcpp::stop("'bin_width' must be > 0.");
    }
    if (confidence < 0.0 || confidence > 1.0) {
        Rcpp::stop("'confidence' must be between 0 and 1.");
    }
    if (!std::isfinite(lambda) || lambda <= 0.0) {
        Rcpp::stop("'lambda' must be finite and > 0.");
    }
    kernel = canonical_distance_kernel(kernel);
    if (kernel.empty()) {
        Rcpp::stop("'kernel' must be 'Gaussian'/'gaussian' or 'Cauchy'/'cauchy'.");
    }
    const DistanceKernel kernel_method =
        distance_kernel_from_string(kernel);
    RowMajorMatrix<float32_t> targets = as_Matrix<float32_t>(target_vals);
    RowMajorMatrix<float32_t> samples = as_Matrix<float32_t>(sample_vals);
    RowMajorMatrix<double> refdens = as_Matrix<double>(ref_density);
    RowMajorMatrix<double> target_xy = get_XY(target_vals);

    const int nr = targets.rows();
    const int ns = samples.rows();
    if (targets.cols() < 4 || (targets.cols() - 2) % 2 != 0) {
        Rcpp::stop("'target_vals' must contain x, y, predicted RS, and observed RS columns.");
    }
    if (samples.cols() != targets.cols()) {
        Rcpp::stop("'sample_vals' must have the same columns as 'target_vals'.");
    }

    const int nvar = (targets.cols() - 2) / 2;
    const int ndim = nvar + 2;
    const float32_t binwidth = static_cast<float32_t>(bin_width);

    double scale;
    int64_t r2;
    const double radius_m = radius_km * 1000.0;
    if (geographic) {
        scale = 1000000.0;
        const double r_deg = radius_m / 111320.0;
        const int64_t r_micro = static_cast<int64_t>(r_deg * scale);
        r2 = r_micro * r_micro;
    } else {
        scale = 100.0;
        const int64_t r_scaled = static_cast<int64_t>(radius_m * scale);
        r2 = r_scaled * r_scaled;
    }

    std::vector<int64_t> sample_x(ns), sample_y(ns);
    for (int i = 0; i < ns; ++i) {
        sample_x[i] = static_cast<int64_t>(samples(i, 0) * scale);
        sample_y[i] = static_cast<int64_t>(samples(i, 1) * scale);
    }

    const float32_t xypenalty = static_cast<float32_t>(xy_penalty);
    std::vector<float32_t> xystats(xy_stats.begin(), xy_stats.end());
    if (xystats[2] == 0.0f || xystats[3] == 0.0f) {
        Rcpp::stop("'xy_stats' standard deviations (3rd and 4th elements) must be non-zero.");
    }

    samples.col(0) = ((samples.col(0).array() - xystats[0]) / xystats[2]) * xypenalty;
    samples.col(1) = ((samples.col(1).array() - xystats[1]) / xystats[3]) * xypenalty;
    targets.col(0) = ((targets.col(0).array() - xystats[0]) / xystats[2]) * xypenalty;
    targets.col(1) = ((targets.col(1).array() - xystats[1]) / xystats[3]) * xypenalty;

    std::vector<double> predicted_use(ns, 0.0);
    std::vector<double> density_use(ns, 0.0);
    std::vector<double> condition_use(ns, 0.0);

    #ifdef _OPENMP
        if (num_threads < 1) num_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);
    #endif

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nr; ++i)
    {
        const auto cell_rem = targets.row(i).leftCols(ndim);
        const auto cell_obs = targets.row(i).rightCols(nvar);

        if ((cell_obs.array().isNaN()).any()) {
            continue;
        }

        const int64_t x = static_cast<int64_t>(target_xy(i, 0) * scale);
        const int64_t y = static_cast<int64_t>(target_xy(i, 1) * scale);
        int64_t cos_scale = 0;
        if (geographic) {
            cos_scale = static_cast<int64_t>(
                std::cos(target_xy(i, 1) * DEG_2_RAD) * 1000000
            );
        }

        std::vector<int> knn_env = combined_Search(
            sample_x, sample_y, samples, cell_rem,
            x, y, r2, k_env, ndim, cos_scale, geographic
        );

        for (const int site : knn_env) {
            #pragma omp atomic update
            predicted_use[site] += 1.0;
        }

        std::vector<int> candidate_sites;
        std::vector<double> predicted_distances;
        std::vector<double> probabilities;
        candidate_sites.reserve(knn_env.size());
        predicted_distances.reserve(knn_env.size());
        probabilities.reserve(knn_env.size());

        for (const int site : knn_env)
        {
            const auto sample_pred = samples.row(site).middleCols(2, nvar);
            const float32_t predicted_distance =
                (cell_rem.rightCols(nvar) - sample_pred).template lpNorm<1>();

            if (exclude_slef && predicted_distance < binwidth) {
                continue;
            }

            const auto sample_obs = samples.row(site).rightCols(nvar);
            const float32_t observed_distance =
                (cell_obs - sample_obs).template lpNorm<1>();
            const double probability = get_prob_value(
                refdens,
                predicted_distance,
                observed_distance,
                binwidth,
                bin_num,
                offset
            );

            candidate_sites.push_back(site);
            predicted_distances.push_back(static_cast<double>(predicted_distance));
            probabilities.push_back(probability);
        }

        if (probabilities.empty()) {
            continue;
        }

        const std::vector<int> order = descending_probability_order(probabilities);
        const int n_keep = std::min(k_rs, static_cast<int>(order.size()));
        std::vector<int> selected_sites(n_keep);
        std::vector<double> selected_distances(n_keep);
        std::vector<double> selected_probabilities(n_keep);

        for (int k = 0; k < n_keep; ++k) {
            const int candidate = order[k];
            selected_sites[k] = candidate_sites[candidate];
            selected_distances[k] = predicted_distances[candidate];
            selected_probabilities[k] = probabilities[k];
        }

        for (const int site : selected_sites) {
            #pragma omp atomic update
            density_use[site] += 1.0;
        }

        double max_log_weight = -std::numeric_limits<double>::infinity();
        for (const double distance : selected_distances) {
            const double log_weight = log_distance_weight(
                distance,
                lambda,
                kernel_method
            );
            if (std::isfinite(log_weight)) {
                max_log_weight = std::max(max_log_weight, log_weight);
            }
        }
        if (!std::isfinite(max_log_weight)) {
            continue;
        }

        std::vector<double> distance_weights(n_keep);
        std::vector<double> maximum_scores(n_keep);
        double weight_sum = 0.0;
        double maximum_score = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < n_keep; ++k) {
            const double log_weight = log_distance_weight(
                selected_distances[k],
                lambda,
                kernel_method
            );
            const double relative_exponent = std::min(
                0.0,
                log_weight - max_log_weight
            );
            distance_weights[k] = std::exp(relative_exponent);
            maximum_scores[k] = weighted_max
                ? selected_probabilities[k] * distance_weights[k]
                : selected_probabilities[k];
            weight_sum += distance_weights[k];
            maximum_score = std::max(maximum_score, maximum_scores[k]);
        }
        if (!(weight_sum > 0.0) || !std::isfinite(weight_sum)) {
            continue;
        }

        if (!std::isnan(boost_factor)) {
            const double log_boost = std::log(boost_factor);
            double max_boosted_log_weight =
                -std::numeric_limits<double>::infinity();
            for (int k = 0; k < n_keep; ++k) {
                const double log_weight = log_distance_weight(
                    selected_distances[k],
                    lambda,
                    kernel_method
                );
                const double boosted_log_weight = log_weight +
                    (k == 0 ? log_boost : 0.0);
                max_boosted_log_weight = std::max(
                    max_boosted_log_weight,
                    boosted_log_weight
                );
            }

            double boosted_weight_sum = 0.0;
            for (int k = 0; k < n_keep; ++k) {
                const double log_weight = log_distance_weight(
                    selected_distances[k],
                    lambda,
                    kernel_method
                );
                const double boosted_log_weight = log_weight +
                    (k == 0 ? log_boost : 0.0);
                distance_weights[k] = std::exp(std::min(
                    0.0,
                    boosted_log_weight - max_boosted_log_weight
                ));
                boosted_weight_sum += distance_weights[k];
            }
            if (!(boosted_weight_sum > 0.0) ||
                !std::isfinite(boosted_weight_sum)) {
                continue;
            }

            for (int k = 0; k < n_keep; ++k) {
                const int site = selected_sites[k];
                const double attribution =
                    distance_weights[k] / boosted_weight_sum;
                #pragma omp atomic update
                condition_use[site] += attribution;
            }
        } else {
            int maximum_ties = 0;
            for (const double score : maximum_scores) {
                if (score == maximum_score) {
                    maximum_ties += 1;
                }
            }

            for (int k = 0; k < n_keep; ++k) {
                double attribution =
                    (1.0 - confidence) * distance_weights[k] / weight_sum;
                if (maximum_scores[k] == maximum_score) {
                    attribution += confidence / maximum_ties;
                }

                const int site = selected_sites[k];
                #pragma omp atomic update
                condition_use[site] += attribution;
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("predicted") = Rcpp::wrap(predicted_use),
        Rcpp::Named("density") = Rcpp::wrap(density_use),
        Rcpp::Named("condition") = Rcpp::wrap(condition_use)
    );
}

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
