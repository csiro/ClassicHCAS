#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>


enum class DistanceKernel
{
    Gaussian,
    Cauchy
};


inline std::string canonical_distance_kernel(std::string kernel)
{
    std::transform(
        kernel.begin(),
        kernel.end(),
        kernel.begin(),
        [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        }
    );

    return (kernel == "gaussian" || kernel == "cauchy")
        ? kernel
        : "";
}


inline DistanceKernel distance_kernel_from_string(
    const std::string& kernel)
{
    return canonical_distance_kernel(kernel) == "cauchy"
        ? DistanceKernel::Cauchy
        : DistanceKernel::Gaussian;
}


// Log kernel weight for predicted RS L1 distance. Gaussian is
// exp(-(distance / lambda)^2). Cauchy uses the standard Cauchy shape,
// normalised to weight one at zero:
//
//   1 / (1 + (distance / lambda)^2)
inline double log_distance_weight(
    const double predicted_distance,
    const double lambda,
    const DistanceKernel kernel)
{
    const double distance = std::max(predicted_distance, 0.0);

    if (kernel == DistanceKernel::Gaussian) {
        const double distance_scale = distance / lambda;
        return -(distance_scale * distance_scale);
    }

    if (distance == 0.0) {
        return 0.0;
    }

    const double log_distance = std::log(distance);
    const double log_lambda = std::log(lambda);
    const double twice_log_ratio =
        2.0 * (log_distance - log_lambda);
    const double log_one_plus_ratio_sq = twice_log_ratio > 0.0
        ? twice_log_ratio + std::log1p(std::exp(-twice_log_ratio))
        : std::log1p(std::exp(twice_log_ratio));

    return -log_one_plus_ratio_sq;
}


inline double distance_weight(
    const double predicted_distance,
    const double lambda,
    const DistanceKernel kernel)
{
    return std::exp(log_distance_weight(
        predicted_distance,
        lambda,
        kernel
    ));
}
