#ifndef KD_POINTS_H
#define KD_POINTS_H

// dimension of multi-dimensional points (ENV points)
static const int MDP_DIM = 16;

/**
 * A set of point classes and methods definitions for KD-tree algorithm.
 */

#include <vector>
#include <array>
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


// functions definitions
std::vector<XYPoints> as_XYPoints(const Lightweight_matrix<double>& mat);

std::vector<MDPoints> as_MDPoints(const Lightweight_matrix<double>& matrix, int start_col, int end_col);

std::vector<int> KD_KNN(const Lightweight_matrix<double>& data,
                        const Lightweight_matrix<double>& query_dat,
                        const int i,
                        const int k,
                        const int L,
                        const int start,
                        const int end);

#endif // KD_POINTS_H
