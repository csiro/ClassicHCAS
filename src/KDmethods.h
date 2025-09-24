#pragma once

#include "Matrix.h"
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


// function to convert RowMajorMatrix to vector of XYPoints; only uses 0-1 columns
inline std::vector<XYPoints> as_XYPoints(const RowMajorMatrix<double> &mat)
{
    std::vector<XYPoints> points;
    points.reserve(mat.rows());

    for (int i = 0; i < mat.rows(); i++)
    {
        double x = mat(i, 0);
        double y = mat(i, 1);
        points.push_back(XYPoints(x, y));
    }

    return points;
}

