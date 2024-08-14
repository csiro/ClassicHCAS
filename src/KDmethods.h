#ifndef _KD_POINTS_H
#define _KD_POINTS_H

// dimension of multi-dimensional points (ENV points)
static const int MDP_DIM = 16;

/**
 * A set of point classes and methods for KD-tree algorithm.
 */

#include <vector>
#include <array>
#include <stdexcept>
#include "Lightweight_matrix.h"

// define a 2D point class for KD-tree method
class XYPoints : public std::array<double, 2>
{
public:
    // dimension of the points
    static const int DIM = 2;

    // constructor
    XYPoints(double x, double y)
    {
        (*this)[0] = x;
        (*this)[1] = y;
    }
};

// function to convert Lightweight_matrix to vector of XYPoints; only uses 0-1 columns
std::vector<XYPoints> as_XYPoints(const Lightweight_matrix<double>& mat)
{
    std::vector<XYPoints> points;
    for (int i = 0; i < mat.nrow(); i++)
    {
        double x = mat(i, 0);
        double y = mat(i, 1);
        points.push_back(XYPoints(x, y)); // the push back maight be time consuming! replace by initilized vector?
    }

    return points;
}

/**
 * The dimension of this class must be the same as the dimension of the
 * RS and ENV variables before compiling the code. So, if the number of RS/ENV
 * variables changes this dimension should change accordingly. That would be:
 * the two std::array<double, 7> and DIM = 7; e.g to std::array<double, 14> ...
 */
// define a multi-dimensional point class for KD-tree method
class MDPoints : public std::array<double, MDP_DIM>
{
public:
    // Dimension of the point
    static const int DIM = MDP_DIM;

    // Constructor
    MDPoints() : std::array<double, MDP_DIM>() {}

    // Constructor that takes an array of values
    MDPoints(const double *values)
    {
        for (int i = 0; i < DIM; ++i)
        {
            (*this)[i] = values[i];
        }
    }

    // Constructor that takes a std::vector<double>
    MDPoints(const std::vector<double>& values)
    {
        if (values.size() != DIM)
        {
            throw std::invalid_argument("Vector size doesn't match point dimension!");
        }
        for (int i = 0; i < DIM; ++i)
        {
            (*this)[i] = values[i];
        }
    }
};

// function to create a vector of MDPoints from a subset of columns in Lightweight_matrix
std::vector<MDPoints> as_MDPoints(const Lightweight_matrix<double>& matrix, int start_col, int end_col)
{
    // get the number of rows
    int nrows = matrix.nrow();

    // create a vector to store the MDPoints objects
    std::vector<MDPoints> points(nrows);

    for (int i = 0; i < nrows; ++i)
    {
        MDPoints point;
        // copy data from the subset of columns
        for (int j = start_col; j < end_col; ++j)
        {
            point[j - start_col] = matrix(i, j); // adjust index for data vector
        }
        points[i] = point;
    }

    return points;
}

// create a KD-tree and perform a KNN search
std::vector<int> KD_KNN(const Lightweight_matrix<double>& data,      // main data (sample values)
                        const Lightweight_matrix<double>& query_dat, // data matrix containing query point (raster stack)
                        const int i,                                 // the row of query_dat to be the query point
                        const int k,                                 // number of neighbours
                        const int L,                                 // 1 for L1 (Manhattan) and 2 for L2 (Euclidean) distances
                        const int start,                             // starting column in raster stack matrix
                        const int end)
{

    // define multi-dimensional points only for env variables
    std::vector<MDPoints> points = as_MDPoints(data, start, end);
    // build k-d tree for searching knn
    kdt::KDTree<MDPoints> the_tree(points);
    // define the query point of i; the target cell in raster in multi-domension
    // allocate memory for the data array
    const int dim = MDPoints::DIM; // get the dimension from the point object
    std::vector<double> row_i(dim);
    // copy the target/query row's data (env variables) into the vector
    for (int j = 0; j < dim; ++j)
    {
        row_i[j] = query_dat(i, j + start);
    }
    // define the query point
    MDPoints query_point(row_i);
    // find all 50 nearest-neighbour of target point
    std::vector<int> indices = the_tree.knnSearch(query_point, k, L);

    return indices;
}

#endif // _KD_POINTS_H
