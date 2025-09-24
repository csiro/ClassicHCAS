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
#include "Matrix.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerVector density_cpp(
        const Rcpp::NumericMatrix &rast,
        const Rcpp::NumericMatrix &xy,
        const double radius_km = 200.0,
        const double scale = 100000.0,
        int num_threads = -1
)
{
    RowMajorMatrix<double> raster = as_Matrix<double>(rast);
    RowMajorMatrix<double> samples = as_Matrix<double>(xy);

    const int nr = raster.rows();
    // output vector
    std::vector<int> out(nr);

    // define radius in m and divide by scale because the search doesn't correct it
    const double radius = radius_km * 1000 / scale;

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
        const auto row = raster.row(i);
        // Take care of NaN cells
        if ((row.array().isNaN()).any())
        {
            #pragma omp critical
            out[i] = NA_INTEGER;
        }
        else
        {
            // get the target points xy from raster stack;
            const auto x = raster(i, 0);
            const auto y = raster(i, 1);
            // define the query point of i; a cell in raster
            XYPoints query(x, y);
            // find all ref samples in a radius e.g. 200km
            std::vector<int> indices = kdtree.radiusSearch(query, radius, 2);
            int count = indices.size();

            #pragma omp critical
            out[i] = count;
        }
    }

    return Rcpp::wrap(out);
}

