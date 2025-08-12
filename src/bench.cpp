// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
// local scripts
#include "KDtree.h"
#include "KDmethods.h"
#include "Quickisort.h"
#include "Lightweight_matrix.h"
#include "Helper.h"

using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::NumericMatrix bench_cpp(
    const Rcpp::NumericMatrix &rast_stack,  // a raster stack to read everything at once: x,y,rs,env
    const Rcpp::NumericMatrix &sample_vals, // extraction of values of rstack using smaples xy: x,y,RS,ENV
    const Rcpp::NumericMatrix &histogram,   // cleaned histogram
    const Rcpp::NumericVector &xy_stats,    // mean(x), mean(y), sd(x), sd(y); ORDER MATTERS
    const double xy_penalty = 0.0,          // penalising env nearest neighbour searching for geographic distance
    unsigned int num_vars = 0,              // can be calculate from the stack
    const double scale = 100000.0,          // correction scale for CRS; 1 if metric, 100,000 otherwise
    const double within_km = 200,           // radius in kilometers to consider ref points
    const int k_env = 50,                   // number of ENV nn to select
    const int k_rs = 20,                    // number of RS/Prob values to select
    const double bin_width = 0.05,          // histogram bin width
    const int bin_num = 400,                // number of bins in histogram
    const int offset = 0,                   // offset of histogram
    const double pnorm = 1,                 // distance fraction; 1 = L1; 2 = L2; <1 fractional distance
    const double confidence = 0.5,          // the LDC confidence index; default 0.5
    const double lambda = 2.0,              // the lambda of the Cauchy weighting
    const bool exclude_slef = true,         // whether to exclude a benchmark sample from assessing itself
    const bool make_su = false,             // whether to produce SU map
    int num_threads = -1)                   // -1 or 0 utilises all available threads
{

    // convert all Rcpp matrices to custom C++ matrix [faster computation and avoids OpenMp conflicts]
    Lightweight_matrix<double> samples(sample_vals);
    Lightweight_matrix<double> rstack(rast_stack);
    Lightweight_matrix<double> histo(histogram);

    // add scaled x and y columns to the end and apply weight to it;
    add_xy_scaled(samples, xy_stats, xy_penalty);
    add_xy_scaled(rstack, xy_stats, xy_penalty);

    // ncol MUST BE after adding xy_scaled, so it has the two extra columns
    const int nc = rstack.ncol();
    const int nr = rstack.nrow();

    // output condition vector
    std::vector<condition> condition_vect(nr);

    // if number of layers is not supplied, define it
    // minus 2 for ij, then divide by 2 for equal layers
    if (num_vars < 3)
        num_vars = (nc - 4) / 2; // 4 for x, y, x_scaled, y_scaled

    // end column for RS varibales; start column for ENV
    const int end_rs = num_vars + 2;
    const int end_pr = nc - 2;

    // define radius in m and divide by scale because the search doesn't correct it
    const double radius = within_km * 1000 / scale;

    // convert points from matrix to vector of XYPoints
    std::vector<XYPoints> points = as_XYPoints(samples);
    // build k-d tree for searching 200km points
    kdt::KDTree<XYPoints> kdtree(points);

    // set the number of threads
    #ifdef _OPENMP
        if (num_threads < 1) num_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);
    #endif

    // iterate over the raster cells (rows of matrix)
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nr; i++)
    {
        // get the target points xy from raster stack;
        const auto x = rstack(i, 0);
        const auto y = rstack(i, 1);
        // define the query point of i; a cell in raster
        XYPoints query(x, y);
        // find all ref samples in a radius e.g. 200km
        std::vector<int> indices = kdtree.radiusSearch(query, radius, 2); // 2 for L2 distance
        // now filter the matrix rows based on indices within 200km radius
        // this sub_sample contain all the points in this radius that can be several thousands
        Lightweight_matrix<double> sub_sample = filter_matrix(samples, indices);

        // find the 50 nearest samples to the target point in ENV dist
        std::vector<int> knn_env = KD_KNN(sub_sample, rstack, i, k_env, 1, end_rs, nc); // 1 for L1 distance

        // rs and env distance of the 50 nearest env neighbour
        std::vector<double> prdist_vect(k_env);
        // a vector to copy histogram values; initialize with size 50
        std::vector<double> histo_vect(k_env);

        // to avoid complexity and overhead of kdtree, use for loop for the 50 records
        int n = 0;
        for (const auto& j : knn_env)
        {
            // calculate fractional distance for RS and ENV of target cell (i) to 50-NN
            double rsdist_sum(0.);
            for (int m = 2; m < end_rs; m++)
            {
                rsdist_sum += std::pow(std::fabs(rstack(i, m) - sub_sample(j, m)), pnorm);
            }
            double prdist_sum(0.);
            for (int m = end_rs; m < end_pr; m++)
            {
                prdist_sum += std::pow(std::fabs(rstack(i, m) - sub_sample(j, m)), pnorm);
            }
            double rsdist = std::pow(rsdist_sum, 1.0 / pnorm);
            double prdist = std::pow(prdist_sum, 1.0 / pnorm);

            // skip the point if exclude-slef is true and
            // prdist is in the first bin i.e < 1 * bw
            if (exclude_slef && prdist < bin_width)
            {
                continue;
            }
            else
            {
                prdist_vect[n] = prdist;
                histo_vect[n] = get_prob(histo, prdist, rsdist, bin_width, bin_num, offset);
                n++;
            }
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
        condition wcond = get_condition(prob_sorted, prob_sorted[0], pr_dist, confidence, lambda);

        #pragma omp critical
        condition_vect[i] = wcond;
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

