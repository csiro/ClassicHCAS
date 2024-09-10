// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>


// a function to linearly interpolate between two spline points
double interpolate(
    const std::vector<double> &x,
    const std::vector<double> &y,
    double x_input)
{
    // find the interval x_input falls into
    auto it = std::lower_bound(x.begin(), x.end(), x_input);

    if (it == x.end() || it == x.begin())
    {
        // if x_input is out of range, return the boundary value
        if (x_input <= x.front())
        {
            return y.front();
        }
        else
        {
            return y.back();
        }
    }

    size_t idx = std::distance(x.begin(), it) - 1;

    // linear interpolation
    double x0 = x[idx], x1 = x[idx + 1];
    double y0 = y[idx], y1 = y[idx + 1];

    double y_out = y0 + ((x_input - x0) * (y1 - y0)) / (x1 - x0);

    return y_out;
}

// [[Rcpp::export]]
std::vector<double> spline_rescale(
    const Rcpp::NumericVector &x,
    const std::vector<double> &spline_x,
    const std::vector<double> &spline_y)
{
    // define the output vector
    std::vector<double> out_vect;
    out_vect.reserve(x.size());

    for (const auto &v : x)
    {
        if (std::isnan(v))
        {
            out_vect.push_back(R_NaN);
        }
        else
        {
            double value = interpolate(spline_x, spline_y, v);
            // make sure values doesn't pass 1; the flat line at the end of the graph
            value = (value > 1) ? 1 : value;
            out_vect.push_back(value);
        }
    }

    return out_vect;
}

