// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <cstdint>
#include "Lightweight_matrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// calculate squared distance using xy
inline double distance(const double x1, const double y1, const double x2, const double y2)
{
    const double dx = x2 - x1;
    const double dy = y2 - y1;
    double dist = (dx * dx) + (dy * dy);

    return dist;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix histo_cpp(
    const Rcpp::NumericMatrix &rs_vals,
    const Rcpp::NumericMatrix &pr_vals,
    const Rcpp::NumericMatrix &samples_xy,
    const double within_km = 1000.0, // radius in kilometers to consider ref points
    const double scale = 100000.0,   // correction scale for CRS; 1 if metric, 100,000 otherwise
    const double bin_width = 0.05,
    const int bin_num = 650,
    int num_threads = -1)
{
    // convert all Rcpp matrices to custom C++ matrix [faster computation and avoids OpenMp conflicts]
    Lightweight_matrix<double> rs(rs_vals);
    Lightweight_matrix<double> pr(pr_vals);
    Lightweight_matrix<double> samples(samples_xy);
    // define the output matrix;
    Lightweight_matrix<uint64_t> hist(bin_num, bin_num);

    const int nr = rs.nrow();
    const int nc = rs.ncol();

    // squared distance for faster calculation
    // calculating distance squared is faster; one-time rather than many sqrt
    const double squared_dist = std::pow(within_km * 1000 / scale, 2);

    // set the number of threads
    #ifdef _OPENMP
        if (num_threads < 1) num_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);
    #endif

    // iterate over each row and compare it with all other rows (second loop)
    #pragma omp parallel for
    for (int i = 0; i < nr; i++)
    {
        // second iteration over rows
        for (int j = 0; j < nr; j++)
        {
            // get the XY from samples for the first sample
            const auto xi = samples(i, 0);
            const auto yi = samples(i, 1);
            // get the XY from samples for the second sample
            const auto xj = samples(j, 0);
            const auto yj = samples(j, 1);

            // calculate geographic squared distance
            double dist = distance(xi, yi, xj, yj);
            // only calculate for point within e.g. 1000km distance
            if (dist < squared_dist)
            {
                // calculate Manhattan distance
                double rsdist(0.);
                double prdist(0.);
                for (int m = 0; m < nc; m++)
                {
                    rsdist += std::fabs(rs(i, m) - rs(j, m));
                    prdist += std::fabs(pr(i, m) - pr(j, m));
                }

                // now update the hist with the distances
                // calculate the row and column
                int ii = std::floor(prdist / bin_width); // they are both double
                int jj = std::floor(rsdist / bin_width);
                // for now, ignore the values overshooting
                if (ii < bin_num && jj < bin_num)
                {
                    #pragma omp atomic
                    hist(ii, jj)++;
                }
            }
        }
    }

    // now transpose the matrix; then revert rows to save it in the correct position
    Rcpp::IntegerMatrix matrix(bin_num, bin_num);
    for (int i = 0; i < bin_num; ++i)
    {
        int ix = bin_num - 1 - i; // compute the reversed row index directly
        for (int j = 0; j < bin_num; ++j)
        {
            matrix(i, j) = hist(j, ix); // transpose and reverse in one step
        }
    }

    return matrix;
}

