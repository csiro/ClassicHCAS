// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

struct Node {
    int id;
    int r1;
    int r2;
    int c1;
    int c2;
    double wt;
};

struct Split {
    bool ok;
    bool by_row;
    int point;
    double diff;
    double score;
};

inline bool is_finite(const double x) {
    return std::isfinite(x);
}

double region_sum(
    const NumericMatrix& x,
    const int r1,
    const int r2,
    const int c1,
    const int c2
) {
    double s = 0.0;
    for (int r = r1 - 1; r < r2; ++r) {
        for (int c = c1 - 1; c < c2; ++c) {
            const double v = x(r, c);
            if (is_finite(v)) s += v;
        }
    }
    return s;
}

double aspect_ratio(const int h, const int w) {
    const int mn = std::min(h, w);
    if (mn <= 0) return std::numeric_limits<double>::infinity();
    const int mx = std::max(h, w);
    return static_cast<double>(mx) / static_cast<double>(mn);
}

Split split_point(
    const NumericMatrix& x,
    const Node& node,
    const bool by_row
) {
    const int span = by_row ? (node.r2 - node.r1 + 1) : (node.c2 - node.c1 + 1);
    if (span < 2) {
        return {false, by_row, 0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    }

    std::vector<double> dimsum(static_cast<size_t>(span), 0.0);
    if (by_row) {
        for (int i = 0; i < span; ++i) {
            double s = 0.0;
            const int r = node.r1 + i;
            for (int c = node.c1; c <= node.c2; ++c) {
                const double v = x(r - 1, c - 1);
                if (is_finite(v)) s += v;
            }
            dimsum[static_cast<size_t>(i)] = s;
        }
    } else {
        for (int i = 0; i < span; ++i) {
            double s = 0.0;
            const int c = node.c1 + i;
            for (int r = node.r1; r <= node.r2; ++r) {
                const double v = x(r - 1, c - 1);
                if (is_finite(v)) s += v;
            }
            dimsum[static_cast<size_t>(i)] = s;
        }
    }

    const double total = std::accumulate(dimsum.begin(), dimsum.end(), 0.0);
    const double half = total / 2.0;

    double cs = 0.0;
    double best_abs = std::numeric_limits<double>::infinity();
    int best_pt = 1;
    for (int i = 0; i < (span - 1); ++i) {
        cs += dimsum[static_cast<size_t>(i)];
        const double d = std::fabs(cs - half);
        if (d < best_abs) {
            best_abs = d;
            best_pt = i + 1;
        }
    }

    const double left = std::accumulate(dimsum.begin(), dimsum.begin() + best_pt, 0.0);
    const double right = total - left;
    const double diff = std::fabs(left - right);

    return {true, by_row, best_pt, diff, 0.0};
}

double split_score(const Node& node, const Split& split) {
    const int h = node.r2 - node.r1 + 1;
    const int w = node.c2 - node.c1 + 1;

    int h1 = h;
    int h2 = h;
    int w1 = w;
    int w2 = w;

    if (split.by_row) {
        h1 = split.point;
        h2 = h - split.point;
    } else {
        w1 = split.point;
        w2 = w - split.point;
    }

    const double ratio1 = aspect_ratio(h1, w1);
    const double ratio2 = aspect_ratio(h2, w2);
    const double compactness = ((ratio1 - 1.0) + (ratio2 - 1.0)) / 2.0;
    const double imbalance = (node.wt > 0.0) ? (split.diff / node.wt) : 0.0;

    // Low score means more balanced and less skinny children.
    const double compactness_weight = 0.05;
    return imbalance + (compactness_weight * compactness);
}

Split choose_split(
    const NumericMatrix& x,
    const Node& node,
    const std::string& method
) {
    const bool can_row = node.r2 > node.r1;
    const bool can_col = node.c2 > node.c1;
    if (!can_row && !can_col) {
        return {false, true, 0, 0.0, 0.0};
    }

    Split row = {false, true, 0, 0.0, 0.0};
    Split col = {false, false, 0, 0.0, 0.0};

    if (can_row) {
        row = split_point(x, node, true);
        if (row.ok) row.score = split_score(node, row);
    }
    if (can_col) {
        col = split_point(x, node, false);
        if (col.ok) col.score = split_score(node, col);
    }

    if (method == "row") {
        if (row.ok) return row;
        if (col.ok) return col;
        return {false, true, 0, 0.0, 0.0};
    }

    if (method == "col") {
        if (col.ok) return col;
        if (row.ok) return row;
        return {false, true, 0, 0.0, 0.0};
    }

    if (method == "both") {
        const int h = node.r2 - node.r1 + 1;
        const int w = node.c2 - node.c1 + 1;
        if (h >= w) {
            if (row.ok) return row;
            if (col.ok) return col;
        } else {
            if (col.ok) return col;
            if (row.ok) return row;
        }
        return {false, true, 0, 0.0, 0.0};
    }

    // method == "best"
    if (row.ok && col.ok) {
        if (row.score <= col.score) return row;
        return col;
    }
    if (row.ok) return row;
    if (col.ok) return col;
    return {false, true, 0, 0.0, 0.0};
}

} // namespace

// [[Rcpp::export]]
IntegerMatrix tiling_cpp(
    const NumericMatrix& x,
    const int n_tiles,
    const std::string method = "best",
    const bool exact = true
) {
    if (n_tiles < 1) {
        stop("'n_tiles' must be a positive integer.");
    }
    if (method != "best" && method != "row" && method != "col" && method != "both") {
        stop("'method' must be one of 'best', 'row', 'col', or 'both'.");
    }

    const int n_rows = x.nrow();
    const int n_cols = x.ncol();
    if (n_rows < 1 || n_cols < 1) {
        stop("'x' must be a non-empty matrix.");
    }

    const double total_sum = region_sum(x, 1, n_rows, 1, n_cols);
    const double target_sum = total_sum / static_cast<double>(n_tiles);

    std::vector<Node> nodes;
    nodes.reserve(static_cast<size_t>(n_tiles));
    nodes.push_back({1, 1, n_rows, 1, n_cols, total_sum});

    int next_id = 1;

    while (true) {
        if (exact) {
            if (static_cast<int>(nodes.size()) >= n_tiles) break;
        } else {
            double max_wt = -std::numeric_limits<double>::infinity();
            for (const auto& node : nodes) {
                max_wt = std::max(max_wt, node.wt);
            }
            if (max_wt <= target_sum) break;
        }

        int selected = -1;
        double selected_wt = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < nodes.size(); ++i) {
            const Node& node = nodes[i];
            const bool splittable = (node.r2 > node.r1) || (node.c2 > node.c1);
            if (!splittable) continue;
            if (!exact && node.wt <= target_sum) continue;
            if (node.wt > selected_wt) {
                selected_wt = node.wt;
                selected = static_cast<int>(i);
            }
        }

        if (selected < 0) break;

        const Split split = choose_split(x, nodes[selected], method);
        if (!split.ok) break;

        Node left = nodes[selected];
        Node right = nodes[selected];
        right.id = ++next_id;

        if (split.by_row) {
            left.r2 = left.r1 + split.point - 1;
            right.r1 = left.r2 + 1;
        } else {
            left.c2 = left.c1 + split.point - 1;
            right.c1 = left.c2 + 1;
        }

        left.wt = region_sum(x, left.r1, left.r2, left.c1, left.c2);
        right.wt = region_sum(x, right.r1, right.r2, right.c1, right.c2);

        nodes[selected] = left;
        nodes.push_back(right);
    }

    IntegerMatrix result(n_rows, n_cols);
    std::fill(result.begin(), result.end(), NA_INTEGER);

    for (const auto& node : nodes) {
        for (int r = node.r1 - 1; r < node.r2; ++r) {
            for (int c = node.c1 - 1; c < node.c2; ++c) {
                result(r, c) = node.id;
            }
        }
    }

    return result;
}
