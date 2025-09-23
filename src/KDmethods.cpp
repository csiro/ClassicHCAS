#include <stdexcept>
#include "KDmethods.h"
#include "KDtree.h"
#include "Matrix.h"

// function to convert RowMajorMatrix to vector of XYPoints; only uses 0-1 columns
std::vector<XYPoints> as_XYPoints(const RowMajorMatrix<double> &mat)
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

