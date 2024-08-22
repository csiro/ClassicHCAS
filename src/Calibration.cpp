// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

// [[Rcpp::export]]
std::vector<double> linear_rescale(
    const Rcpp::NumericVector &x,
    const double low_point,
    const double high_point,
    const double max_point,
    const double low_target = 0.1,
    const double high_target = 0.944)
{
    // define the output vector
    std::vector<double> out_vect;
    out_vect.reserve(x.size());

    // calculate the slope of the lines
    const double slope_low = low_target / low_point;
    const double slope_mid = (high_target - low_target) / (high_point - low_point);
    const double slope_hig = (1 - high_target) / (max_point - high_point);

    for (const auto &v : x)
    {
        if (std::isnan(v))
        {
            out_vect.push_back(R_NaN);
        }
        else
        {
            if (v < low_point)
            {
                out_vect.push_back(v * slope_low);
            }
            else
            {
                if (v < high_point)
                {
                    out_vect.push_back(low_target + (v - low_point) * slope_mid);
                }
                else
                {
                    double value = high_target + (v - high_point) * slope_hig;
                    // make sure values doesn't pass 1; the flat line at the end of the graph
                    value = (value > 1) ? 1 : value;
                    out_vect.push_back(value);
                }
            }
        }
    }

    return out_vect;
}


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
