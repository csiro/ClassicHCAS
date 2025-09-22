#ifndef KD_POINTS_H
#define KD_POINTS_H

#include "Matrix.h"

// dimension of multi-dimensional points (ENV points)
static const int MDP_DIM = 16;

/**
 * A set of point classes and methods definitions for KD-tree algorithm.
 */

#include <vector>
#include <array>

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

// functions definitions
std::vector<XYPoints> as_XYPoints(const RowMajorMatrix<double>& mat);

#endif // KD_POINTS_H
