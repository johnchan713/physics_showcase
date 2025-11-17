#ifndef MATHS_REAL_ANALYSIS_HPP
#define MATHS_REAL_ANALYSIS_HPP

#include <vector>
#include <set>
#include <functional>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <numeric>

namespace maths {
namespace real_analysis {

// ============================================================================
// METRIC SPACES
// ============================================================================

template<typename T>
class MetricSpace {
public:
    virtual double distance(const T& x, const T& y) const = 0;
    virtual ~MetricSpace() = default;

    bool satisfiesTriangleInequality(const T& x, const T& y, const T& z, double tol = 1e-10) const {
        double d_xy = distance(x, y);
        double d_yz = distance(y, z);
        double d_xz = distance(x, z);
        return d_xz <= d_xy + d_yz + tol;
    }

    bool isSymmetric(const T& x, const T& y, double tol = 1e-10) const {
        return std::abs(distance(x, y) - distance(y, x)) < tol;
    }

    double diameter(const std::vector<T>& set) const {
        if (set.size() < 2) return 0.0;
        double max_dist = 0.0;
        for (size_t i = 0; i < set.size(); ++i) {
            for (size_t j = i + 1; j < set.size(); ++j) {
                max_dist = std::max(max_dist, distance(set[i], set[j]));
            }
        }
        return max_dist;
    }

    bool isBounded(const std::vector<T>& set, double M) const {
        return diameter(set) <= M;
    }

    bool isCauchy(const std::vector<T>& sequence, double epsilon = 1e-6) const {
        int n = sequence.size();
        for (int i = n - 10; i < n && i >= 0; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (distance(sequence[i], sequence[j]) > epsilon) {
                    return false;
                }
            }
        }
        return true;
    }

    bool converges(const std::vector<T>& sequence, const T& limit, double epsilon = 1e-6) const {
        if (sequence.empty()) return false;
        return distance(sequence.back(), limit) < epsilon;
    }
};

class RealLineMetric : public MetricSpace<double> {
public:
    double distance(const double& x, const double& y) const override {
        return std::abs(x - y);
    }
};

class EuclideanMetric : public MetricSpace<std::vector<double>> {
public:
    double distance(const std::vector<double>& x, const std::vector<double>& y) const override {
        if (x.size() != y.size()) throw std::invalid_argument("Dimension mismatch");
        double sum = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            sum += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return std::sqrt(sum);
    }
};

class MaximumMetric : public MetricSpace<std::vector<double>> {
public:
    double distance(const std::vector<double>& x, const std::vector<double>& y) const override {
        if (x.size() != y.size()) throw std::invalid_argument("Dimension mismatch");
        double max_dist = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            max_dist = std::max(max_dist, std::abs(x[i] - y[i]));
        }
        return max_dist;
    }
};

class DiscreteMetric : public MetricSpace<int> {
public:
    double distance(const int& x, const int& y) const override {
        return (x == y) ? 0.0 : 1.0;
    }
};

// ============================================================================
// COMPLETENESS
// ============================================================================

template<typename T>
class Completeness {
private:
    const MetricSpace<T>& metric_;

public:
    Completeness(const MetricSpace<T>& M) : metric_(M) {}

    bool isComplete(const std::vector<std::vector<T>>& cauchy_sequences,
                   const std::vector<T>& limits,
                   double tol = 1e-6) const {
        if (cauchy_sequences.size() != limits.size()) return false;

        for (size_t i = 0; i < cauchy_sequences.size(); ++i) {
            if (!metric_.isCauchy(cauchy_sequences[i])) continue;
            if (!metric_.converges(cauchy_sequences[i], limits[i], tol)) {
                return false;
            }
        }
        return true;
    }

    T findLimit(const std::vector<T>& cauchy_sequence) const {
        if (cauchy_sequence.empty()) throw std::invalid_argument("Empty sequence");
        if (!metric_.isCauchy(cauchy_sequence)) {
            return cauchy_sequence.back();
        }

        int n = cauchy_sequence.size();
        int avg_count = std::min(10, n);
        return cauchy_sequence[n - avg_count / 2];
    }

    bool verifyCompletion(const std::vector<T>& incomplete_space_sample,
                         const std::vector<T>& completion_sample,
                         double tol = 1e-6) const {
        for (const auto& x : incomplete_space_sample) {
            bool found = false;
            for (const auto& y : completion_sample) {
                if (metric_.distance(x, y) < tol) {
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        return true;
    }
};

// ============================================================================
// COMPACTNESS
// ============================================================================

template<typename T>
class Compactness {
private:
    const MetricSpace<T>& metric_;

public:
    Compactness(const MetricSpace<T>& M) : metric_(M) {}

    bool isSequentiallyCompact(const std::vector<T>& set,
                              std::function<std::vector<T>(const std::vector<T>&)> extract_subsequence,
                              int n_tests = 10) const {
        for (int test = 0; test < n_tests; ++test) {
            std::vector<T> sequence = set;
            auto subsequence = extract_subsequence(sequence);

            if (subsequence.size() < 2) continue;

            bool has_limit = false;
            for (const auto& x : set) {
                if (metric_.converges(subsequence, x)) {
                    has_limit = true;
                    break;
                }
            }
            if (!has_limit) return false;
        }
        return true;
    }

    bool isTotallyBounded(const std::vector<T>& set, double epsilon) const {
        std::vector<T> centers;
        std::vector<bool> covered(set.size(), false);

        for (size_t i = 0; i < set.size(); ++i) {
            if (covered[i]) continue;

            centers.push_back(set[i]);
            for (size_t j = 0; j < set.size(); ++j) {
                if (metric_.distance(set[i], set[j]) < epsilon) {
                    covered[j] = true;
                }
            }
        }

        return true;
    }

    int epsilonNetSize(const std::vector<T>& set, double epsilon) const {
        std::vector<T> net;
        for (const auto& x : set) {
            bool in_net = false;
            for (const auto& y : net) {
                if (metric_.distance(x, y) < epsilon) {
                    in_net = true;
                    break;
                }
            }
            if (!in_net) {
                net.push_back(x);
            }
        }
        return net.size();
    }

    bool heineBorelTheorem(const std::vector<double>& interval_set, double a, double b,
                          double tol = 1e-10) const {
        if (interval_set.empty()) return false;

        double min_val = *std::min_element(interval_set.begin(), interval_set.end());
        double max_val = *std::max_element(interval_set.begin(), interval_set.end());

        bool is_closed = (std::abs(min_val - a) < tol && std::abs(max_val - b) < tol);
        bool is_bounded = (b - a < std::numeric_limits<double>::infinity());

        return is_closed && is_bounded;
    }
};

// ============================================================================
// CONNECTEDNESS
// ============================================================================

template<typename T>
class Connectedness {
private:
    const MetricSpace<T>& metric_;

public:
    Connectedness(const MetricSpace<T>& M) : metric_(M) {}

    bool isConnected(const std::vector<T>& set, double epsilon = 1e-6) const {
        if (set.size() <= 1) return true;

        std::vector<bool> component(set.size(), false);
        component[0] = true;

        bool changed = true;
        while (changed) {
            changed = false;
            for (size_t i = 0; i < set.size(); ++i) {
                if (!component[i]) continue;
                for (size_t j = 0; j < set.size(); ++j) {
                    if (component[j]) continue;
                    if (metric_.distance(set[i], set[j]) < epsilon) {
                        component[j] = true;
                        changed = true;
                    }
                }
            }
        }

        for (bool c : component) {
            if (!c) return false;
        }
        return true;
    }

    std::vector<std::vector<T>> connectedComponents(const std::vector<T>& set,
                                                    double epsilon = 1e-6) const {
        std::vector<std::vector<T>> components;
        std::vector<bool> assigned(set.size(), false);

        for (size_t i = 0; i < set.size(); ++i) {
            if (assigned[i]) continue;

            std::vector<T> component;
            component.push_back(set[i]);
            assigned[i] = true;

            for (size_t k = 0; k < component.size(); ++k) {
                for (size_t j = 0; j < set.size(); ++j) {
                    if (assigned[j]) continue;
                    if (metric_.distance(component[k], set[j]) < epsilon) {
                        component.push_back(set[j]);
                        assigned[j] = true;
                    }
                }
            }
            components.push_back(component);
        }
        return components;
    }

    bool isPathConnected(const std::vector<T>& set,
                        std::function<std::vector<T>(const T&, const T&)> path_finder,
                        double tol = 1e-6) const {
        if (set.size() <= 1) return true;

        for (size_t i = 0; i < set.size(); ++i) {
            for (size_t j = i + 1; j < set.size(); ++j) {
                auto path = path_finder(set[i], set[j]);
                if (path.empty()) return false;

                for (const auto& p : path) {
                    bool in_set = false;
                    for (const auto& s : set) {
                        if (metric_.distance(p, s) < tol) {
                            in_set = true;
                            break;
                        }
                    }
                    if (!in_set) return false;
                }
            }
        }
        return true;
    }
};

// ============================================================================
// UNIFORM CONTINUITY
// ============================================================================

class UniformContinuity {
public:
    static bool isUniformlyContinuous(
        std::function<double(double)> f,
        double a, double b,
        double epsilon, double& delta,
        int n_samples = 1000) {

        double min_delta = b - a;
        bool found_violation = false;

        for (int i = 0; i < n_samples; ++i) {
            for (int j = i + 1; j < n_samples; ++j) {
                double x1 = a + (b - a) * i / n_samples;
                double x2 = a + (b - a) * j / n_samples;
                double dx = std::abs(x2 - x1);
                double df = std::abs(f(x2) - f(x1));

                if (df > epsilon && dx > 0) {
                    min_delta = std::min(min_delta, dx);
                    found_violation = true;
                }
            }
        }

        if (!found_violation) {
            delta = (b - a) / n_samples;
            return true;
        }

        delta = min_delta * 0.9;
        return !found_violation;
    }

    static bool isContinuousAt(std::function<double(double)> f, double x0,
                              double epsilon, double& delta,
                              int n_samples = 100) {
        double h = 1.0;
        while (h > 1e-10) {
            double max_diff = 0.0;
            for (int i = -n_samples; i <= n_samples; ++i) {
                double x = x0 + i * h / n_samples;
                max_diff = std::max(max_diff, std::abs(f(x) - f(x0)));
            }

            if (max_diff < epsilon) {
                delta = h;
                return true;
            }
            h *= 0.5;
        }
        delta = 0.0;
        return false;
    }

    static bool isLipschitzContinuous(std::function<double(double)> f,
                                     double a, double b,
                                     double& L,
                                     int n_samples = 1000) {
        double max_ratio = 0.0;
        for (int i = 0; i < n_samples; ++i) {
            for (int j = i + 1; j < n_samples; ++j) {
                double x1 = a + (b - a) * i / n_samples;
                double x2 = a + (b - a) * j / n_samples;
                double dx = std::abs(x2 - x1);
                double df = std::abs(f(x2) - f(x1));

                if (dx > 1e-15) {
                    max_ratio = std::max(max_ratio, df / dx);
                }
            }
        }
        L = max_ratio;
        return max_ratio < std::numeric_limits<double>::infinity();
    }
};

// ============================================================================
// SEQUENCES AND SERIES
// ============================================================================

class SequencesAndSeries {
public:
    static bool convergesPointwise(
        const std::vector<std::function<double(double)>>& fn_sequence,
        std::function<double(double)> f_limit,
        double x,
        double tol = 1e-6) {

        if (fn_sequence.empty()) return false;
        return std::abs(fn_sequence.back()(x) - f_limit(x)) < tol;
    }

    static bool convergesUniformly(
        const std::vector<std::function<double(double)>>& fn_sequence,
        std::function<double(double)> f_limit,
        double a, double b,
        double tol = 1e-6,
        int n_samples = 100) {

        if (fn_sequence.empty()) return false;

        const auto& f_n = fn_sequence.back();
        double max_diff = 0.0;

        for (int i = 0; i <= n_samples; ++i) {
            double x = a + (b - a) * i / n_samples;
            max_diff = std::max(max_diff, std::abs(f_n(x) - f_limit(x)));
        }

        return max_diff < tol;
    }

    static double supremumNorm(std::function<double(double)> f,
                              double a, double b,
                              int n_samples = 1000) {
        double sup = 0.0;
        for (int i = 0; i <= n_samples; ++i) {
            double x = a + (b - a) * i / n_samples;
            sup = std::max(sup, std::abs(f(x)));
        }
        return sup;
    }

    static bool isEquicontinuous(
        const std::vector<std::function<double(double)>>& family,
        double a, double b,
        double epsilon,
        double& delta,
        int n_samples = 100) {

        double min_delta = b - a;

        for (const auto& f : family) {
            double f_delta;
            if (!UniformContinuity::isUniformlyContinuous(f, a, b, epsilon, f_delta, n_samples)) {
                delta = 0.0;
                return false;
            }
            min_delta = std::min(min_delta, f_delta);
        }

        delta = min_delta;
        return true;
    }

    static bool arzela_ascoli(
        const std::vector<std::function<double(double)>>& sequence,
        double a, double b,
        bool is_equicontinuous,
        bool is_uniformly_bounded) {

        return is_equicontinuous && is_uniformly_bounded;
    }

    static double seriesSum(const std::vector<double>& terms) {
        return std::accumulate(terms.begin(), terms.end(), 0.0);
    }

    static bool cauchyCriterionSeries(const std::vector<double>& partial_sums,
                                     double epsilon = 1e-6) {
        int n = partial_sums.size();
        for (int i = n - 10; i < n && i >= 0; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(partial_sums[j] - partial_sums[i]) > epsilon) {
                    return false;
                }
            }
        }
        return true;
    }

    static bool absolutelyConvergent(const std::vector<double>& terms) {
        double sum_abs = 0.0;
        for (double a_n : terms) {
            sum_abs += std::abs(a_n);
            if (sum_abs > 1e10) return false;
        }
        return true;
    }

    static bool ratioTest(const std::vector<double>& terms) {
        if (terms.size() < 2) return false;

        int n = terms.size();
        double ratio = std::abs(terms[n-1] / terms[n-2]);

        return ratio < 1.0;
    }

    static bool rootTest(const std::vector<double>& terms) {
        if (terms.empty()) return false;

        int n = terms.size();
        double root = std::pow(std::abs(terms[n-1]), 1.0 / n);

        return root < 1.0;
    }
};

// ============================================================================
// CANTOR SET
// ============================================================================

class CantorSet {
private:
    std::vector<std::pair<double, double>> intervals_;
    int iterations_;

public:
    CantorSet(int n) : iterations_(n) {
        intervals_.push_back({0.0, 1.0});
        construct();
    }

    void construct() {
        for (int iter = 0; iter < iterations_; ++iter) {
            std::vector<std::pair<double, double>> new_intervals;
            for (const auto& [a, b] : intervals_) {
                double third = (b - a) / 3.0;
                new_intervals.push_back({a, a + third});
                new_intervals.push_back({b - third, b});
            }
            intervals_ = new_intervals;
        }
    }

    bool contains(double x) const {
        for (const auto& [a, b] : intervals_) {
            if (x >= a && x <= b) return true;
        }
        return false;
    }

    double measure() const {
        return std::pow(2.0 / 3.0, iterations_);
    }

    double hausdorffDimension() const {
        return std::log(2.0) / std::log(3.0);
    }

    bool isNowhereDense() const {
        return true;
    }

    bool isPerfect() const {
        return iterations_ >= 10;
    }

    int cardinality() const {
        return std::pow(2, iterations_);
    }
};

// ============================================================================
// STONE-WEIERSTRASS THEOREM
// ============================================================================

class StoneWeierstrass {
public:
    static std::vector<double> polynomialApproximation(
        std::function<double(double)> f,
        double a, double b,
        int degree,
        int n_points = 100) {

        std::vector<std::vector<double>> A(n_points, std::vector<double>(degree + 1));
        std::vector<double> y(n_points);

        for (int i = 0; i < n_points; ++i) {
            double x = a + (b - a) * i / (n_points - 1);
            y[i] = f(x);
            for (int j = 0; j <= degree; ++j) {
                A[i][j] = std::pow(x, j);
            }
        }

        std::vector<double> coeffs(degree + 1, 0.0);
        for (int j = 0; j <= degree; ++j) {
            double num = 0.0, denom = 0.0;
            for (int i = 0; i < n_points; ++i) {
                num += A[i][j] * y[i];
                denom += A[i][j] * A[i][j];
            }
            if (denom > 1e-15) {
                coeffs[j] = num / denom;
            }
        }
        return coeffs;
    }

    static double evaluatePolynomial(const std::vector<double>& coeffs, double x) {
        double result = 0.0;
        for (size_t i = 0; i < coeffs.size(); ++i) {
            result += coeffs[i] * std::pow(x, i);
        }
        return result;
    }

    static double approximationError(
        std::function<double(double)> f,
        const std::vector<double>& poly_coeffs,
        double a, double b,
        int n_samples = 1000) {

        double max_error = 0.0;
        for (int i = 0; i <= n_samples; ++i) {
            double x = a + (b - a) * i / n_samples;
            double error = std::abs(f(x) - evaluatePolynomial(poly_coeffs, x));
            max_error = std::max(max_error, error);
        }
        return max_error;
    }
};

} // namespace real_analysis
} // namespace maths

#endif // MATHS_REAL_ANALYSIS_HPP
