// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// local scripts
#include "KDtree.h"
#include "KDmethods.h"
#include "Quickisort.h"
#include "Matrix.h"
#include "Helper.h"

using namespace Rcpp;

using DTYPE = float;


// [[Rcpp::export]]
Rcpp::NumericMatrix bench_cpp(
    const Rcpp::NumericMatrix &raster_vals, // a raster stack to read everything at once: x,y,rs,env
    const Rcpp::NumericMatrix &sample_vals, // extraction of values of rstack using smaples xy: x,y,RS,ENV
    const Rcpp::NumericMatrix &histogram,   // cleaned histogram
    const Rcpp::NumericVector &xy_stats,    // mean(x), mean(y), sd(x), sd(y); ORDER MATTERS
    const double xy_penalty = 0.0,          // penalising env nearest neighbour searching for geographic distance
    const double scale = 100000.0,          // correction scale for CRS; 1 if metric, 100,000 otherwise
    const double within_km = 200,           // radius in kilometers to consider ref points
    const int k_env = 50,                   // number of ENV nn to select
    const int k_rs = 20,                    // number of RS/Prob values to select
    const double bin_width = 0.05,          // histogram bin width
    const int bin_num = 400,                // number of bins in histogram
    const int offset = 0,                   // offset of histogram
    const double confidence = 0.5,          // the LDC confidence index; default 0.5
    const double lambda = 2.0,              // the lambda of the Cauchy weighting
    const bool exclude_slef = true,         // whether to exclude a benchmark sample from assessing itself
    const bool make_su = false,             // whether to produce SU map
    int num_threads = -1)                   // -1 or 0 utilises all available threads
{
    // Copy xy for calculating the correct squared distances;
    // making a copy before reading the rest in the generic format flaot-32
    Rcpp::NumericMatrix xy_r = get_XY(raster_vals);
    Rcpp::NumericMatrix xy_s = get_XY(sample_vals);
    RowMajorMatrix<double> rast_xy = as_Matrix<double>(xy_r);
    RowMajorMatrix<double> sample_xy = as_Matrix<double>(xy_s);

    // convert all Rcpp matrices to custom C++ matrix [faster computation and avoids OpenMp conflicts]
    RowMajorMatrix<DTYPE> rstack = as_Matrix<DTYPE>(raster_vals);
    RowMajorMatrix<DTYPE> samples = as_Matrix<DTYPE>(sample_vals);
    RowMajorMatrix<double> histo = as_Matrix<double>(histogram);

    const int nr = rstack.rows();
    int nvar = (samples.cols() - 2) / 2; // number of RS vars
    int ndim = nvar + 2; // number of multi-variate space REM + XY

    // // Copy xy for calculating the correct squared distances
    // const RowMajorMatrix<double> sample_xy = samples.leftCols(2);
    // const RowMajorMatrix<double> rast_xy = rstack.leftCols(2);
    // Copy OBS and MOD RS from raster cells and samples
    RowMajorMatrix<DTYPE> REM = samples.leftCols(ndim);
    RowMajorMatrix<DTYPE> OBS = samples.rightCols(nvar);

    // Normalise XY and apply weight to it; after filtering XY columns
    // REM.col(0) = ((REM.col(0).array() - xy_stats[0]) / xy_stats[2]) * xy_penalty;
    // REM.col(1) = ((REM.col(1).array() - xy_stats[1]) / xy_stats[3]) * xy_penalty;
    // rstack.col(0) = ((rstack.col(0).array() - xy_stats[0]) / xy_stats[2]) * xy_penalty;
    // rstack.col(1) = ((rstack.col(1).array() - xy_stats[1]) / xy_stats[3]) * xy_penalty;
    REM.col(0) = ((REM.col(0).array() - static_cast<DTYPE>(xy_stats[0])) / static_cast<DTYPE>(xy_stats[2])) * static_cast<DTYPE>(xy_penalty);
    REM.col(1) = ((REM.col(1).array() - static_cast<DTYPE>(xy_stats[1])) / static_cast<DTYPE>(xy_stats[3])) * static_cast<DTYPE>(xy_penalty);
    rstack.col(0) = ((rstack.col(0).array() - static_cast<DTYPE>(xy_stats[0])) / static_cast<DTYPE>(xy_stats[2])) * static_cast<DTYPE>(xy_penalty);
    rstack.col(1) = ((rstack.col(1).array() - static_cast<DTYPE>(xy_stats[1])) / static_cast<DTYPE>(xy_stats[3])) * static_cast<DTYPE>(xy_penalty);


    // define radius in m and divide by scale because the search doesn't correct it
    const double radius = within_km * 1000.0 / scale;

    // convert points from matrix to vector of XYPoints
    std::vector<XYPoints> points = as_XYPoints(sample_xy);
    // build k-d tree for searching 200km points
    kdt::KDTree<XYPoints> kdtree(points);

    // output condition vector
    std::vector<condition> condition_vect(nr);

    // set the number of threads
    #ifdef _OPENMP
        if (num_threads < 1) num_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);
    #endif

    // iterate over the raster cells (rows of matrix)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nr; i++)
    {
        const auto cell_rem = rstack.row(i).leftCols(ndim);
        const auto cell_obs = rstack.row(i).rightCols(nvar);

        // if ((cell_obs.array().isNaN()).any())
        // {
        //     output(i) = std::numeric_limits<double>::quiet_NaN();
        // }
        // else
        // {
        // }

        // get the target points xy from raster stack;
        const double x = rast_xy(i, 0);
        const double y = rast_xy(i, 1);
        // define the query point of i; a cell in raster
        XYPoints query(x, y);
        // find all ref samples in a radius e.g. 200km
        std::vector<int> indices = kdtree.radiusSearch(query, radius, 2); // 2 for L2 distance
        // now filter the matrix rows based on indices within 200km radius
        // this sub_sample contain all the points in this radius that can be several thousands
        RowMajorMatrix<DTYPE> sub_rem = filter_Matrix(REM, indices);
        RowMajorMatrix<DTYPE> sub_obs = filter_Matrix(OBS, indices);

        // find the 50 nearest samples to the target point in ENV dist
        std::vector<int> knn_env = KNN_Search(sub_rem, cell_rem, k_env);

        // rs and env distance of the 50 nearest env neighbour
        std::vector<double> prdist_vect(k_env);
        std::vector<double> histo_vect(k_env);

        // to avoid complexity and overhead of kdtree, use for loop for the 50 records
        int n = 0;
        for (const auto& j : knn_env)
        {
            // vectorised L1 distance over all columns
            double rsdist = (cell_obs.template cast<double>() - sub_obs.row(j).template cast<double>()).template lpNorm<1>();
            // the xy coordinates should be ignored here so only, nvar right columns
            double prdist = (cell_rem.rightCols(nvar).template cast<double>() - sub_rem.row(j).rightCols(nvar).template cast<double>()).template lpNorm<1>();

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

