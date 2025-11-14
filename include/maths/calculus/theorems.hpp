#ifndef MATHS_CALCULUS_THEOREMS_HPP
#define MATHS_CALCULUS_THEOREMS_HPP

#include <functional>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <limits>

/**
 * @file theorems.hpp
 * @brief Fundamental calculus theorems
 *
 * Implements:
 * - Intermediate Value Theorem (IVT)
 * - Mean Value Theorem (MVT)
 * - Rolle's Theorem
 * - Extreme Value Theorem
 * - Related theorems and applications
 */

namespace maths::calculus {

/**
 * @class IntermediateValueTheorem
 * @brief Intermediate Value Theorem (IVT)
 *
 * If f is continuous on [a, b] and k is between f(a) and f(b),
 * then there exists c in (a, b) such that f(c) = k
 */
class IntermediateValueTheorem {
public:
    /**
     * @brief Statement of IVT
     */
    static std::string statement() {
        return "Intermediate Value Theorem:\n"
               "\n"
               "If f: [a, b] → ℝ is continuous and k is between f(a) and f(b),\n"
               "then ∃c ∈ (a, b) such that f(c) = k\n"
               "\n"
               "Consequence: Continuous functions on closed intervals\n"
               "             take on all intermediate values";
    }

    /**
     * @brief Check IVT conditions
     *
     * @param f Function (must be continuous on [a, b])
     * @param a Left endpoint
     * @param b Right endpoint
     * @param k Target value
     * @return true if k is between f(a) and f(b)
     */
    static bool checkConditions(std::function<double(double)> f,
                                double a, double b, double k) {
        if (a >= b) {
            throw std::invalid_argument("Must have a < b");
        }

        double fa = f(a);
        double fb = f(b);

        // k must be between f(a) and f(b)
        return (fa <= k && k <= fb) || (fb <= k && k <= fa);
    }

    /**
     * @brief Find c such that f(c) = k using bisection method
     *
     * @param f Continuous function on [a, b]
     * @param a Left endpoint
     * @param b Right endpoint
     * @param k Target value
     * @param tolerance Acceptable error
     * @param max_iterations Maximum iterations
     * @return c such that |f(c) - k| < tolerance
     */
    static double findRoot(std::function<double(double)> f,
                          double a, double b, double k = 0.0,
                          double tolerance = 1e-10,
                          int max_iterations = 100) {
        if (!checkConditions(f, a, b, k)) {
            throw std::invalid_argument("IVT conditions not satisfied");
        }

        // Shift function so we're finding f(x) - k = 0
        auto g = [&](double x) { return f(x) - k; };

        double fa = g(a);
        double fb = g(b);

        if (std::abs(fa) < tolerance) return a;
        if (std::abs(fb) < tolerance) return b;

        // Bisection method
        for (int i = 0; i < max_iterations; ++i) {
            double c = (a + b) / 2.0;
            double fc = g(c);

            if (std::abs(fc) < tolerance || (b - a) / 2.0 < tolerance) {
                return c;
            }

            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }
        }

        return (a + b) / 2.0;
    }

    /**
     * @brief Example: Prove √2 exists
     */
    static double proofSqrt2Exists() {
        // f(x) = x² is continuous on [1, 2]
        // f(1) = 1, f(2) = 4, and 2 is between 1 and 4
        // By IVT, ∃c ∈ (1, 2) such that f(c) = 2
        // i.e., c² = 2, so c = √2

        auto f = [](double x) { return x * x; };
        return findRoot(f, 1.0, 2.0, 2.0);  // Returns ≈ 1.414213...
    }

    /**
     * @brief Bolzano's Theorem (special case of IVT)
     *
     * If f is continuous on [a, b] and f(a) and f(b) have opposite signs,
     * then ∃c ∈ (a, b) such that f(c) = 0
     */
    static std::string bolzanoTheorem() {
        return "Bolzano's Theorem (IVT with k = 0):\n"
               "\n"
               "If f is continuous on [a, b] and f(a)·f(b) < 0,\n"
               "then ∃c ∈ (a, b) such that f(c) = 0\n"
               "\n"
               "Application: Root-finding, existence of zeros";
    }
};

/**
 * @class MeanValueTheorem
 * @brief Mean Value Theorem (MVT)
 *
 * If f is continuous on [a, b] and differentiable on (a, b),
 * then ∃c ∈ (a, b) such that f'(c) = [f(b) - f(a)] / (b - a)
 */
class MeanValueTheorem {
public:
    /**
     * @brief Statement of MVT
     */
    static std::string statement() {
        return "Mean Value Theorem:\n"
               "\n"
               "If f: [a, b] → ℝ is continuous on [a, b] and differentiable on (a, b),\n"
               "then ∃c ∈ (a, b) such that:\n"
               "\n"
               "f'(c) = [f(b) - f(a)] / (b - a)\n"
               "\n"
               "Interpretation: ∃ point where instantaneous rate = average rate";
    }

    /**
     * @brief Geometric interpretation
     */
    static std::string geometricInterpretation() {
        return "Geometric interpretation of MVT:\n"
               "\n"
               "Secant line slope: [f(b) - f(a)] / (b - a)\n"
               "Tangent line slope at c: f'(c)\n"
               "\n"
               "MVT: ∃c where tangent line is parallel to secant line\n"
               "\n"
               "Visual: At least one point where tangent || chord";
    }

    /**
     * @brief Calculate average rate of change (secant slope)
     */
    static double averageRateOfChange(std::function<double(double)> f,
                                     double a, double b) {
        if (a == b) {
            throw std::invalid_argument("Must have a ≠ b");
        }
        return (f(b) - f(a)) / (b - a);
    }

    /**
     * @brief Find c satisfying MVT numerically
     *
     * @param f Function
     * @param df Derivative of f
     * @param a Left endpoint
     * @param b Right endpoint
     * @param tolerance Acceptable error
     * @return c such that |f'(c) - avg_rate| < tolerance
     */
    static double findMVTPoint(std::function<double(double)> f,
                               std::function<double(double)> df,
                               double a, double b,
                               double tolerance = 1e-10) {
        double avg_rate = averageRateOfChange(f, a, b);

        // Find c where df(c) = avg_rate
        // i.e., solve df(c) - avg_rate = 0
        auto g = [&](double x) { return df(x) - avg_rate; };

        // Use bisection (assumes g changes sign)
        double ga = g(a + 1e-9);  // Slightly inside interval
        double gb = g(b - 1e-9);

        if (ga * gb > 0) {
            // May need to subdivide interval or use different method
            // For now, sample the interval
            double min_diff = std::abs(ga);
            double best_c = a;

            int n_samples = 1000;
            for (int i = 0; i <= n_samples; ++i) {
                double x = a + (b - a) * i / n_samples;
                double diff = std::abs(g(x));
                if (diff < min_diff) {
                    min_diff = diff;
                    best_c = x;
                }
            }
            return best_c;
        }

        // Bisection method
        double left = a + 1e-9;
        double right = b - 1e-9;

        for (int i = 0; i < 100; ++i) {
            double c = (left + right) / 2.0;
            double gc = g(c);

            if (std::abs(gc) < tolerance) {
                return c;
            }

            if (ga * gc < 0) {
                right = c;
                gb = gc;
            } else {
                left = c;
                ga = gc;
            }
        }

        return (left + right) / 2.0;
    }

    /**
     * @brief Example: f(x) = x² on [0, 2]
     *
     * f'(c) = 2c = [f(2) - f(0)] / 2 = 4/2 = 2
     * Therefore c = 1
     */
    static double exampleQuadratic() {
        auto f = [](double x) { return x * x; };
        auto df = [](double x) { return 2.0 * x; };

        return findMVTPoint(f, df, 0.0, 2.0);  // Returns c = 1.0
    }

    /**
     * @brief Cauchy's Mean Value Theorem (generalized MVT)
     *
     * If f and g are continuous on [a, b] and differentiable on (a, b),
     * and g'(x) ≠ 0 for all x ∈ (a, b),
     * then ∃c ∈ (a, b) such that:
     *
     * f'(c) / g'(c) = [f(b) - f(a)] / [g(b) - g(a)]
     */
    static std::string cauchyMVT() {
        return "Cauchy's Mean Value Theorem:\n"
               "\n"
               "If f, g continuous on [a, b], differentiable on (a, b),\n"
               "and g'(x) ≠ 0 on (a, b), then ∃c ∈ (a, b) such that:\n"
               "\n"
               "f'(c) / g'(c) = [f(b) - f(a)] / [g(b) - g(a)]\n"
               "\n"
               "Reduces to standard MVT when g(x) = x";
    }
};

/**
 * @class RollesTheorem
 * @brief Rolle's Theorem (special case of MVT)
 *
 * If f is continuous on [a, b], differentiable on (a, b),
 * and f(a) = f(b), then ∃c ∈ (a, b) such that f'(c) = 0
 */
class RollesTheorem {
public:
    /**
     * @brief Statement of Rolle's Theorem
     */
    static std::string statement() {
        return "Rolle's Theorem:\n"
               "\n"
               "If f: [a, b] → ℝ is continuous on [a, b],\n"
               "differentiable on (a, b), and f(a) = f(b),\n"
               "then ∃c ∈ (a, b) such that f'(c) = 0\n"
               "\n"
               "Interpretation: ∃ horizontal tangent between equal endpoints";
    }

    /**
     * @brief Geometric interpretation
     */
    static std::string geometricInterpretation() {
        return "Geometric interpretation:\n"
               "\n"
               "If function returns to same height (f(a) = f(b)),\n"
               "then somewhere in between there must be a\n"
               "horizontal tangent (local max or min)\n"
               "\n"
               "Visual: Peak or valley between equal endpoints";
    }

    /**
     * @brief Check Rolle's Theorem conditions
     */
    static bool checkConditions(std::function<double(double)> f,
                                double a, double b,
                                double tolerance = 1e-10) {
        return std::abs(f(a) - f(b)) < tolerance;
    }

    /**
     * @brief Find c where f'(c) = 0
     */
    static double findCriticalPoint(std::function<double(double)> df,
                                   double a, double b,
                                   double tolerance = 1e-10) {
        // Find where df(c) = 0
        // Sample the interval to find sign changes
        int n_samples = 1000;
        double step = (b - a) / n_samples;

        for (int i = 0; i < n_samples; ++i) {
            double x1 = a + i * step;
            double x2 = a + (i + 1) * step;

            double y1 = df(x1);
            double y2 = df(x2);

            if (std::abs(y1) < tolerance) return x1;

            if (y1 * y2 < 0) {
                // Sign change, use bisection
                for (int j = 0; j < 50; ++j) {
                    double mid = (x1 + x2) / 2.0;
                    double ymid = df(mid);

                    if (std::abs(ymid) < tolerance) {
                        return mid;
                    }

                    if (y1 * ymid < 0) {
                        x2 = mid;
                        y2 = ymid;
                    } else {
                        x1 = mid;
                        y1 = ymid;
                    }
                }
                return (x1 + x2) / 2.0;
            }
        }

        // If no sign change found, return point with smallest |df|
        double min_abs = std::abs(df(a));
        double best_c = a;

        for (int i = 0; i <= n_samples; ++i) {
            double x = a + (b - a) * i / n_samples;
            double abs_df = std::abs(df(x));
            if (abs_df < min_abs) {
                min_abs = abs_df;
                best_c = x;
            }
        }

        return best_c;
    }

    /**
     * @brief Example: f(x) = x² - 4x + 3 on [1, 3]
     *
     * f(1) = 0, f(3) = 0
     * f'(x) = 2x - 4 = 0 → c = 2
     */
    static double exampleQuadratic() {
        auto df = [](double x) { return 2.0 * x - 4.0; };
        return findCriticalPoint(df, 1.0, 3.0);  // Returns c = 2.0
    }

    /**
     * @brief Relationship to MVT
     */
    static std::string relationshipToMVT() {
        return "Rolle's Theorem ⊂ Mean Value Theorem:\n"
               "\n"
               "MVT: f'(c) = [f(b) - f(a)] / (b - a)\n"
               "\n"
               "If f(a) = f(b), then [f(b) - f(a)] = 0\n"
               "Therefore f'(c) = 0 (Rolle's Theorem)\n"
               "\n"
               "Rolle's is special case of MVT with f(a) = f(b)";
    }
};

/**
 * @class ExtremeValueTheorem
 * @brief Extreme Value Theorem (EVT)
 *
 * If f is continuous on closed interval [a, b],
 * then f attains both maximum and minimum values
 */
class ExtremeValueTheorem {
public:
    /**
     * @brief Statement of EVT
     */
    static std::string statement() {
        return "Extreme Value Theorem:\n"
               "\n"
               "If f: [a, b] → ℝ is continuous on [a, b],\n"
               "then f attains both absolute maximum and minimum on [a, b]\n"
               "\n"
               "i.e., ∃c₁, c₂ ∈ [a, b] such that:\n"
               "f(c₁) ≤ f(x) ≤ f(c₂) for all x ∈ [a, b]";
    }

    /**
     * @brief Find extreme values numerically
     *
     * @return {min_value, max_value, min_location, max_location}
     */
    static std::vector<double> findExtremes(std::function<double(double)> f,
                                           double a, double b,
                                           int n_samples = 10000) {
        if (a >= b) {
            throw std::invalid_argument("Must have a < b");
        }

        double min_val = f(a);
        double max_val = f(a);
        double min_loc = a;
        double max_loc = a;

        for (int i = 0; i <= n_samples; ++i) {
            double x = a + (b - a) * i / n_samples;
            double fx = f(x);

            if (fx < min_val) {
                min_val = fx;
                min_loc = x;
            }
            if (fx > max_val) {
                max_val = fx;
                max_loc = x;
            }
        }

        return {min_val, max_val, min_loc, max_loc};
    }

    /**
     * @brief Importance of closed and bounded interval
     */
    static std::string importanceOfClosedInterval() {
        return "EVT requires CLOSED and BOUNDED interval:\n"
               "\n"
               "Counterexample (open): f(x) = x on (0, 1)\n"
               "  Continuous but no maximum (sup = 1, not attained)\n"
               "\n"
               "Counterexample (unbounded): f(x) = x on [0, ∞)\n"
               "  Continuous but no maximum (unbounded)\n"
               "\n"
               "Counterexample (discontinuous): f(x) = 1/x on (0, 1]\n"
               "  No minimum (lim_{x→0⁺} f(x) = ∞)";
    }

    /**
     * @brief Applications
     */
    static std::string applications() {
        return "Applications of EVT:\n"
               "\n"
               "1. Optimization: Guarantees optimal solutions exist\n"
               "2. Physics: Equilibrium states (minimum energy)\n"
               "3. Economics: Optimal production levels\n"
               "4. Engineering: Safety margins (maximum stress)\n"
               "\n"
               "EVT ensures extrema exist; calculus finds them!";
    }
};

/**
 * @class FundamentalTheoremOfCalculus
 * @brief Fundamental Theorem of Calculus (FTC)
 *
 * Connects differentiation and integration
 */
class FundamentalTheoremOfCalculus {
public:
    /**
     * @brief First Fundamental Theorem (FTC-1)
     *
     * If f is continuous on [a, b] and F(x) = ∫ₐˣ f(t) dt,
     * then F'(x) = f(x) for all x ∈ (a, b)
     */
    static std::string firstFundamentalTheorem() {
        return "First Fundamental Theorem (FTC-1):\n"
               "\n"
               "If f is continuous on [a, b], define:\n"
               "F(x) = ∫ₐˣ f(t) dt\n"
               "\n"
               "Then F'(x) = f(x) for all x ∈ (a, b)\n"
               "\n"
               "Meaning: Derivative of integral recovers original function";
    }

    /**
     * @brief Second Fundamental Theorem (FTC-2)
     *
     * If f is continuous on [a, b] and F is any antiderivative of f,
     * then ∫ₐᵇ f(x) dx = F(b) - F(a)
     */
    static std::string secondFundamentalTheorem() {
        return "Second Fundamental Theorem (FTC-2):\n"
               "\n"
               "If f is continuous on [a, b] and F' = f, then:\n"
               "∫ₐᵇ f(x) dx = F(b) - F(a)\n"
               "\n"
               "Meaning: Definite integral computed via antiderivatives\n"
               "\n"
               "Notation: F(b) - F(a) = [F(x)]ₐᵇ";
    }

    /**
     * @brief Numerical integration using FTC
     *
     * Approximate ∫ₐᵇ f(x) dx using Riemann sums
     */
    static double integrate(std::function<double(double)> f,
                           double a, double b,
                           int n_intervals = 10000) {
        double h = (b - a) / n_intervals;
        double sum = 0.0;

        // Trapezoidal rule (more accurate than left/right Riemann)
        sum += f(a) / 2.0;
        for (int i = 1; i < n_intervals; ++i) {
            sum += f(a + i * h);
        }
        sum += f(b) / 2.0;

        return h * sum;
    }

    /**
     * @brief Connection to MVT
     */
    static std::string connectionToMVT() {
        return "FTC and MVT connection:\n"
               "\n"
               "MVT for integrals: If f continuous on [a, b], then:\n"
               "∃c ∈ [a, b] such that ∫ₐᵇ f(x) dx = f(c)·(b - a)\n"
               "\n"
               "Interpretation: ∃ point where f(c) = average value of f\n"
               "\n"
               "Average value: (1/(b-a)) ∫ₐᵇ f(x) dx";
    }
};

/**
 * @class LHopitalsRule
 * @brief L'Hôpital's Rule for indeterminate forms
 *
 * Application of Cauchy's MVT
 */
class LHopitalsRule {
public:
    /**
     * @brief Statement of L'Hôpital's Rule
     */
    static std::string statement() {
        return "L'Hôpital's Rule:\n"
               "\n"
               "If lim_{x→a} f(x) = lim_{x→a} g(x) = 0 or ±∞,\n"
               "and lim_{x→a} f'(x)/g'(x) exists, then:\n"
               "\n"
               "lim_{x→a} f(x)/g(x) = lim_{x→a} f'(x)/g'(x)\n"
               "\n"
               "Works for 0/0 and ∞/∞ indeterminate forms";
    }

    /**
     * @brief Indeterminate forms
     */
    static std::string indeterminateForms() {
        return "Indeterminate forms (can use L'Hôpital):\n"
               "\n"
               "Direct: 0/0, ∞/∞\n"
               "\n"
               "Convert to 0/0 or ∞/∞:\n"
               "- 0·∞: Rewrite as (0)/(1/∞) = 0/0\n"
               "- ∞ - ∞: Combine fractions\n"
               "- 1^∞, 0^0, ∞^0: Take logarithm first\n"
               "\n"
               "Not indeterminate: 0/∞ = 0, ∞/0 = ∞, etc.";
    }

    /**
     * @brief Example: lim_{x→0} sin(x)/x = 1
     */
    static std::string exampleSinOverX() {
        return "Example: lim_{x→0} sin(x)/x\n"
               "\n"
               "Form: 0/0 (indeterminate)\n"
               "\n"
               "Apply L'Hôpital:\n"
               "lim_{x→0} sin(x)/x = lim_{x→0} cos(x)/1 = cos(0) = 1\n"
               "\n"
               "Therefore: lim_{x→0} sin(x)/x = 1";
    }

    /**
     * @brief Example: lim_{x→∞} x/eˣ = 0
     */
    static std::string examplePolynomialOverExponential() {
        return "Example: lim_{x→∞} x/eˣ\n"
               "\n"
               "Form: ∞/∞ (indeterminate)\n"
               "\n"
               "Apply L'Hôpital:\n"
               "lim_{x→∞} x/eˣ = lim_{x→∞} 1/eˣ = 0\n"
               "\n"
               "Exponential grows faster than polynomial!";
    }
};

} // namespace maths::calculus

#endif // MATHS_CALCULUS_THEOREMS_HPP
