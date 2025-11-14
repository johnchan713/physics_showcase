#ifndef MATHS_TOPOLOGY_METRIC_SPACES_HPP
#define MATHS_TOPOLOGY_METRIC_SPACES_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <string>

/**
 * @file metric_spaces.hpp
 * @brief Metric spaces, convergence, and topological properties
 *
 * Implements foundational concepts from topology:
 * - Metric spaces and distance functions
 * - Open and closed sets
 * - Convergence of sequences
 * - Completeness and Cauchy sequences
 * - Compactness
 * - Continuity in metric spaces
 */

namespace maths::topology {

/**
 * @class MetricSpace
 * @brief Abstract base class for metric spaces (X, d)
 *
 * A metric space is a set X with a distance function d: X × X → ℝ satisfying:
 * 1. d(x, y) ≥ 0 (non-negativity)
 * 2. d(x, y) = 0 ⟺ x = y (identity of indiscernibles)
 * 3. d(x, y) = d(y, x) (symmetry)
 * 4. d(x, z) ≤ d(x, y) + d(y, z) (triangle inequality)
 */
template<typename T>
class MetricSpace {
public:
    using Point = T;
    using Metric = std::function<double(const T&, const T&)>;

protected:
    Metric metric_;

public:
    /**
     * @brief Construct metric space with given distance function
     */
    explicit MetricSpace(Metric metric) : metric_(metric) {}

    /**
     * @brief Compute distance between two points
     */
    double distance(const T& x, const T& y) const {
        return metric_(x, y);
    }

    /**
     * @brief Check if sequence converges to a point
     *
     * A sequence (xₙ) converges to x if for all ε > 0, there exists N such that
     * for all n ≥ N, d(xₙ, x) < ε
     */
    bool converges(const std::vector<T>& sequence, const T& limit,
                   double tolerance = 1e-10) const {
        if (sequence.empty()) return true;

        // Check last few elements are close to limit
        size_t check_count = std::min(size_t(10), sequence.size());
        for (size_t i = sequence.size() - check_count; i < sequence.size(); ++i) {
            if (distance(sequence[i], limit) > tolerance) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Check if sequence is Cauchy
     *
     * A sequence (xₙ) is Cauchy if for all ε > 0, there exists N such that
     * for all m, n ≥ N, d(xₘ, xₙ) < ε
     */
    bool isCauchy(const std::vector<T>& sequence, double tolerance = 1e-10) const {
        if (sequence.size() < 2) return true;

        // Check last few elements are close to each other
        size_t check_count = std::min(size_t(10), sequence.size());
        size_t start = sequence.size() - check_count;

        for (size_t i = start; i < sequence.size(); ++i) {
            for (size_t j = i + 1; j < sequence.size(); ++j) {
                if (distance(sequence[i], sequence[j]) > tolerance) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Open ball B(x, r) = {y : d(x, y) < r}
     */
    bool inOpenBall(const T& point, const T& center, double radius) const {
        return distance(point, center) < radius;
    }

    /**
     * @brief Closed ball B̄(x, r) = {y : d(x, y) ≤ r}
     */
    bool inClosedBall(const T& point, const T& center, double radius) const {
        return distance(point, center) <= radius;
    }

    /**
     * @brief Sphere S(x, r) = {y : d(x, y) = r}
     */
    bool onSphere(const T& point, const T& center, double radius,
                  double tolerance = 1e-10) const {
        return std::abs(distance(point, center) - radius) < tolerance;
    }
};

/**
 * @class EuclideanSpace
 * @brief ℝⁿ with Euclidean metric
 */
class EuclideanSpace : public MetricSpace<std::vector<double>> {
public:
    EuclideanSpace() : MetricSpace([](const std::vector<double>& x,
                                      const std::vector<double>& y) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }
        double sum = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            double diff = x[i] - y[i];
            sum += diff * diff;
        }
        return std::sqrt(sum);
    }) {}
};

/**
 * @class TaxicabSpace
 * @brief ℝⁿ with Manhattan/taxicab metric: d(x,y) = Σ|xᵢ - yᵢ|
 */
class TaxicabSpace : public MetricSpace<std::vector<double>> {
public:
    TaxicabSpace() : MetricSpace([](const std::vector<double>& x,
                                    const std::vector<double>& y) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }
        double sum = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            sum += std::abs(x[i] - y[i]);
        }
        return sum;
    }) {}
};

/**
 * @class MaximumSpace
 * @brief ℝⁿ with maximum metric: d(x,y) = max|xᵢ - yᵢ|
 */
class MaximumSpace : public MetricSpace<std::vector<double>> {
public:
    MaximumSpace() : MetricSpace([](const std::vector<double>& x,
                                    const std::vector<double>& y) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }
        double max_dist = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            max_dist = std::max(max_dist, std::abs(x[i] - y[i]));
        }
        return max_dist;
    }) {}
};

/**
 * @class DiscreteMetricSpace
 * @brief Discrete metric: d(x,y) = 0 if x=y, 1 otherwise
 */
template<typename T>
class DiscreteMetricSpace : public MetricSpace<T> {
public:
    DiscreteMetricSpace() : MetricSpace<T>([](const T& x, const T& y) {
        return (x == y) ? 0.0 : 1.0;
    }) {}
};

/**
 * @class MetricSpaceProperties
 * @brief Theorems and properties of metric spaces
 */
class MetricSpaceProperties {
public:
    /**
     * @brief Basic metric space axioms
     */
    static std::string axioms() {
        return "Metric Space Axioms:\n"
               "\n"
               "A metric d: X × X → ℝ satisfies:\n"
               "\n"
               "1. Non-negativity: d(x, y) ≥ 0\n"
               "2. Identity: d(x, y) = 0 ⟺ x = y\n"
               "3. Symmetry: d(x, y) = d(y, x)\n"
               "4. Triangle inequality: d(x, z) ≤ d(x, y) + d(y, z)\n"
               "\n"
               "From these, we can derive:\n"
               "- Reverse triangle inequality: |d(x, z) - d(y, z)| ≤ d(x, y)";
    }

    /**
     * @brief Convergence in metric spaces
     */
    static std::string convergence() {
        return "Convergence in Metric Spaces:\n"
               "\n"
               "A sequence (xₙ) converges to x (xₙ → x) if:\n"
               "  ∀ε > 0, ∃N: ∀n ≥ N, d(xₙ, x) < ε\n"
               "\n"
               "Properties:\n"
               "- Limits are unique (if they exist)\n"
               "- Convergent sequences are Cauchy\n"
               "- Convergent sequences are bounded\n"
               "- Subsequences of convergent sequences converge to same limit\n"
               "\n"
               "Cauchy sequence: ∀ε > 0, ∃N: ∀m,n ≥ N, d(xₘ, xₙ) < ε";
    }

    /**
     * @brief Completeness
     */
    static std::string completeness() {
        return "Complete Metric Spaces:\n"
               "\n"
               "A metric space is complete if every Cauchy sequence converges.\n"
               "\n"
               "Examples of complete spaces:\n"
               "- ℝⁿ with Euclidean metric\n"
               "- Closed subsets of complete spaces\n"
               "- ℓᵖ spaces (1 ≤ p ≤ ∞)\n"
               "- C[a,b] with sup norm\n"
               "\n"
               "Incomplete space example:\n"
               "- ℚ (rationals) with Euclidean metric\n"
               "\n"
               "Theorem (Cantor): Nested sequence of closed balls\n"
               "in complete space has non-empty intersection.";
    }

    /**
     * @brief Open and closed sets
     */
    static std::string openClosedSets() {
        return "Open and Closed Sets:\n"
               "\n"
               "Open set U: ∀x ∈ U, ∃r > 0: B(x, r) ⊆ U\n"
               "(contains an open ball around each point)\n"
               "\n"
               "Closed set F: X \\ F is open\n"
               "(or: contains all its limit points)\n"
               "\n"
               "Properties:\n"
               "- ∅ and X are both open and closed\n"
               "- Arbitrary unions of open sets are open\n"
               "- Finite intersections of open sets are open\n"
               "- Arbitrary intersections of closed sets are closed\n"
               "- Finite unions of closed sets are closed\n"
               "\n"
               "Closure: cl(A) = smallest closed set containing A\n"
               "Interior: int(A) = largest open set contained in A\n"
               "Boundary: ∂A = cl(A) \\ int(A)";
    }

    /**
     * @brief Compactness
     */
    static std::string compactness() {
        return "Compactness:\n"
               "\n"
               "A set K is compact if every open cover has a finite subcover.\n"
               "\n"
               "Equivalent definitions (in metric spaces):\n"
               "1. Every sequence has a convergent subsequence\n"
               "2. K is complete and totally bounded\n"
               "3. (Heine-Borel in ℝⁿ): K is closed and bounded\n"
               "\n"
               "Properties:\n"
               "- Continuous image of compact set is compact\n"
               "- Closed subset of compact space is compact\n"
               "- Finite union of compact sets is compact\n"
               "- Compact subsets of Hausdorff spaces are closed\n"
               "- Continuous function on compact set is bounded\n"
               "  and attains max/min (Extreme Value Theorem)";
    }

    /**
     * @brief Continuity
     */
    static std::string continuity() {
        return "Continuity in Metric Spaces:\n"
               "\n"
               "f: (X, d_X) → (Y, d_Y) is continuous at x₀ if:\n"
               "  ∀ε > 0, ∃δ > 0: d_X(x, x₀) < δ ⇒ d_Y(f(x), f(x₀)) < ε\n"
               "\n"
               "Equivalent definitions:\n"
               "1. Sequential: xₙ → x ⇒ f(xₙ) → f(x)\n"
               "2. Topological: f⁻¹(U) open for every open U\n"
               "3. Closed sets: f⁻¹(F) closed for every closed F\n"
               "\n"
               "Uniform continuity: ∀ε > 0, ∃δ > 0: ∀x,y,\n"
               "  d_X(x, y) < δ ⇒ d_Y(f(x), f(y)) < ε\n"
               "\n"
               "Lipschitz continuity: ∃L ≥ 0: ∀x,y,\n"
               "  d_Y(f(x), f(y)) ≤ L·d_X(x, y)\n"
               "\n"
               "Theorem: Continuous function on compact set\n"
               "is uniformly continuous.";
    }
};

/**
 * @class Convergence
 * @brief Different modes of convergence
 */
class Convergence {
public:
    /**
     * @brief Pointwise convergence
     */
    static std::string pointwise() {
        return "Pointwise Convergence:\n"
               "\n"
               "Sequence of functions fₙ: X → Y converges pointwise to f if:\n"
               "  ∀x ∈ X, fₙ(x) → f(x)\n"
               "\n"
               "i.e., ∀x ∈ X, ∀ε > 0, ∃N: ∀n ≥ N, d(fₙ(x), f(x)) < ε\n"
               "\n"
               "Note: N may depend on both ε and x.\n"
               "\n"
               "Warning: Pointwise limit of continuous functions\n"
               "need not be continuous!";
    }

    /**
     * @brief Uniform convergence
     */
    static std::string uniform() {
        return "Uniform Convergence:\n"
               "\n"
               "fₙ converges uniformly to f if:\n"
               "  ∀ε > 0, ∃N: ∀n ≥ N, ∀x ∈ X, d(fₙ(x), f(x)) < ε\n"
               "\n"
               "i.e., sup_{x ∈ X} d(fₙ(x), f(x)) → 0\n"
               "\n"
               "Key difference: N depends only on ε, not on x.\n"
               "\n"
               "Theorem: If fₙ continuous and fₙ → f uniformly,\n"
               "then f is continuous.\n"
               "\n"
               "Metric: d_∞(f, g) = sup_{x ∈ X} |f(x) - g(x)|\n"
               "(uniform convergence = convergence in sup norm)";
    }

    /**
     * @brief Comparison of convergence types
     */
    static std::string comparison() {
        return "Convergence Hierarchy:\n"
               "\n"
               "Uniform ⇒ Pointwise\n"
               "\n"
               "But not conversely!\n"
               "\n"
               "Example (pointwise but not uniform):\n"
               "  fₙ(x) = xⁿ on [0,1]\n"
               "  Converges pointwise to f(x) = 0 (x<1), f(1) = 1\n"
               "  But not uniformly (sup|fₙ(x) - f(x)| = 1)\n"
               "\n"
               "Dini's Theorem: On compact space, if fₙ ↓ f\n"
               "(monotone decreasing) and fₙ, f continuous,\n"
               "then fₙ → f uniformly.";
    }
};

} // namespace maths::topology

#endif // MATHS_TOPOLOGY_METRIC_SPACES_HPP
