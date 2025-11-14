#ifndef MATHS_CALCULUS_THEOREMS_HPP
#define MATHS_CALCULUS_THEOREMS_HPP

#include <functional>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <limits>

/**
 * @file theorems.hpp
 * @brief Computational implementations of fundamental calculus theorems
 */

namespace maths::calculus {

/**
 * @class IntermediateValueTheorem
 * @brief Computational tools for Intermediate Value Theorem
 *
 * If f is continuous on [a, b] and k is between f(a) and f(b),
 * then there exists c in (a, b) such that f(c) = k
 */
class IntermediateValueTheorem {
public:
    /**
     * @brief Check if IVT conditions are satisfied
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
        return (fa <= k && k <= fb) || (fb <= k && k <= fa);
    }

    /**
     * @brief Find c such that f(c) = k using bisection method
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

        auto g = [&](double x) { return f(x) - k; };
        double fa = g(a);
        double fb = g(b);

        if (std::abs(fa) < tolerance) return a;
        if (std::abs(fb) < tolerance) return b;

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
     * @brief Compute √2 using IVT
     *
     * f(x) = x² is continuous on [1, 2]
     * f(1) = 1, f(2) = 4, and 2 is between 1 and 4
     * By IVT, ∃c ∈ (1, 2) such that c² = 2
     *
     * @return Approximation of √2 ≈ 1.414213...
     */
    static double computeSqrt2() {
        auto f = [](double x) { return x * x; };
        return findRoot(f, 1.0, 2.0, 2.0);
    }
};

/**
 * @class MeanValueTheorem
 * @brief Computational tools for Mean Value Theorem
 *
 * If f is continuous on [a, b] and differentiable on (a, b),
 * then ∃c ∈ (a, b) such that f'(c) = [f(b) - f(a)] / (b - a)
 */
class MeanValueTheorem {
public:
    /**
     * @brief Calculate average rate of change (secant slope)
     * @param f Function
     * @param a Left endpoint
     * @param b Right endpoint
     * @return [f(b) - f(a)] / (b - a)
     */
    static double averageRateOfChange(std::function<double(double)> f,
                                     double a, double b) {
        if (a == b) {
            throw std::invalid_argument("Must have a ≠ b");
        }
        return (f(b) - f(a)) / (b - a);
    }

    /**
     * @brief Find c satisfying MVT: f'(c) = [f(b)-f(a)]/(b-a)
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
        auto g = [&](double x) { return df(x) - avg_rate; };

        double ga = g(a + 1e-9);
        double gb = g(b - 1e-9);

        if (ga * gb > 0) {
            // Sample the interval to find best approximation
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
     * @brief Example: Find MVT point for f(x) = x² on [0, 2]
     *
     * f'(c) = 2c = [f(2) - f(0)] / 2 = 4/2 = 2
     * Therefore c = 1
     *
     * @return c = 1.0
     */
    static double exampleQuadratic() {
        auto f = [](double x) { return x * x; };
        auto df = [](double x) { return 2.0 * x; };
        return findMVTPoint(f, df, 0.0, 2.0);
    }

    /**
     * @brief Find c satisfying Cauchy MVT
     *
     * f'(c) / g'(c) = [f(b) - f(a)] / [g(b) - g(a)]
     *
     * @param f First function
     * @param df Derivative of f
     * @param g Second function
     * @param dg Derivative of g
     * @param a Left endpoint
     * @param b Right endpoint
     * @param tolerance Acceptable error
     * @return c satisfying Cauchy MVT
     */
    static double findCauchyMVTPoint(std::function<double(double)> f,
                                     std::function<double(double)> df,
                                     std::function<double(double)> g,
                                     std::function<double(double)> dg,
                                     double a, double b,
                                     double tolerance = 1e-10) {
        double ratio = (f(b) - f(a)) / (g(b) - g(a));
        auto h = [&](double x) { return df(x) / dg(x) - ratio; };

        // Sample to find best point
        double min_diff = std::abs(h(a + 1e-9));
        double best_c = a;

        int n_samples = 1000;
        for (int i = 0; i <= n_samples; ++i) {
            double x = a + (b - a) * i / n_samples;
            double diff = std::abs(h(x));
            if (diff < min_diff) {
                min_diff = diff;
                best_c = x;
            }
        }
        return best_c;
    }
};

/**
 * @class RollesTheorem
 * @brief Computational tools for Rolle's Theorem
 *
 * If f is continuous on [a, b], differentiable on (a, b),
 * and f(a) = f(b), then ∃c ∈ (a, b) such that f'(c) = 0
 */
class RollesTheorem {
public:
    /**
     * @brief Check if Rolle's Theorem conditions are satisfied
     * @param f Function
     * @param a Left endpoint
     * @param b Right endpoint
     * @param tolerance Acceptable error
     * @return true if f(a) ≈ f(b)
     */
    static bool checkConditions(std::function<double(double)> f,
                                double a, double b,
                                double tolerance = 1e-10) {
        return std::abs(f(a) - f(b)) < tolerance;
    }

    /**
     * @brief Find c where f'(c) = 0
     * @param df Derivative of function
     * @param a Left endpoint
     * @param b Right endpoint
     * @param tolerance Acceptable error
     * @return c such that |f'(c)| < tolerance
     */
    static double findCriticalPoint(std::function<double(double)> df,
                                   double a, double b,
                                   double tolerance = 1e-10) {
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

        // No sign change found, return point with smallest |df|
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
     *
     * @return c = 2.0
     */
    static double exampleQuadratic() {
        auto df = [](double x) { return 2.0 * x - 4.0; };
        return findCriticalPoint(df, 1.0, 3.0);
    }
};

/**
 * @class ExtremeValueTheorem
 * @brief Computational tools for Extreme Value Theorem
 *
 * If f is continuous on closed interval [a, b],
 * then f attains both maximum and minimum values
 */
class ExtremeValueTheorem {
public:
    /**
     * @brief Find extreme values on [a, b]
     * @param f Continuous function
     * @param a Left endpoint
     * @param b Right endpoint
     * @param n_samples Number of sample points
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
     * @brief Find minimum value and location
     * @return {min_value, min_location}
     */
    static std::vector<double> findMinimum(std::function<double(double)> f,
                                          double a, double b,
                                          int n_samples = 10000) {
        auto extremes = findExtremes(f, a, b, n_samples);
        return {extremes[0], extremes[2]};
    }

    /**
     * @brief Find maximum value and location
     * @return {max_value, max_location}
     */
    static std::vector<double> findMaximum(std::function<double(double)> f,
                                          double a, double b,
                                          int n_samples = 10000) {
        auto extremes = findExtremes(f, a, b, n_samples);
        return {extremes[1], extremes[3]};
    }
};

/**
 * @class FundamentalTheoremOfCalculus
 * @brief Computational integration using FTC
 */
class FundamentalTheoremOfCalculus {
public:
    /**
     * @brief Numerical integration: ∫ₐᵇ f(x) dx
     *
     * Uses trapezoidal rule for accuracy
     *
     * @param f Function to integrate
     * @param a Lower limit
     * @param b Upper limit
     * @param n_intervals Number of intervals
     * @return Approximation of ∫ₐᵇ f(x) dx
     */
    static double integrate(std::function<double(double)> f,
                           double a, double b,
                           int n_intervals = 10000) {
        double h = (b - a) / n_intervals;
        double sum = 0.0;

        // Trapezoidal rule
        sum += f(a) / 2.0;
        for (int i = 1; i < n_intervals; ++i) {
            sum += f(a + i * h);
        }
        sum += f(b) / 2.0;

        return h * sum;
    }

    /**
     * @brief Compute definite integral using antiderivative
     *
     * ∫ₐᵇ f(x) dx = F(b) - F(a) where F' = f
     *
     * @param F Antiderivative of f
     * @param a Lower limit
     * @param b Upper limit
     * @return F(b) - F(a)
     */
    static double evaluateAntiderivative(std::function<double(double)> F,
                                         double a, double b) {
        return F(b) - F(a);
    }

    /**
     * @brief Compute average value of function on [a, b]
     *
     * Average = (1/(b-a)) ∫ₐᵇ f(x) dx
     *
     * @param f Function
     * @param a Lower limit
     * @param b Upper limit
     * @param n_intervals Integration intervals
     * @return Average value of f on [a, b]
     */
    static double averageValue(std::function<double(double)> f,
                              double a, double b,
                              int n_intervals = 10000) {
        return integrate(f, a, b, n_intervals) / (b - a);
    }

    /**
     * @brief Simpson's rule for more accurate integration
     *
     * More accurate than trapezoidal rule
     *
     * @param f Function to integrate
     * @param a Lower limit
     * @param b Upper limit
     * @param n_intervals Number of intervals (must be even)
     * @return Approximation of ∫ₐᵇ f(x) dx
     */
    static double integrateSimpson(std::function<double(double)> f,
                                  double a, double b,
                                  int n_intervals = 10000) {
        if (n_intervals % 2 != 0) n_intervals++;  // Ensure even

        double h = (b - a) / n_intervals;
        double sum = f(a) + f(b);

        for (int i = 1; i < n_intervals; i += 2) {
            sum += 4.0 * f(a + i * h);
        }
        for (int i = 2; i < n_intervals; i += 2) {
            sum += 2.0 * f(a + i * h);
        }

        return (h / 3.0) * sum;
    }
};

/**
 * @class NumericalDerivative
 * @brief Computational derivatives using finite differences
 */
class NumericalDerivative {
public:
    /**
     * @brief Compute derivative using forward difference
     *
     * f'(x) ≈ [f(x+h) - f(x)] / h
     *
     * @param f Function
     * @param x Point
     * @param h Step size (default: 1e-8)
     * @return Approximation of f'(x)
     */
    static double forwardDifference(std::function<double(double)> f,
                                   double x, double h = 1e-8) {
        return (f(x + h) - f(x)) / h;
    }

    /**
     * @brief Compute derivative using central difference
     *
     * f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
     *
     * More accurate than forward/backward difference
     *
     * @param f Function
     * @param x Point
     * @param h Step size (default: 1e-8)
     * @return Approximation of f'(x)
     */
    static double centralDifference(std::function<double(double)> f,
                                   double x, double h = 1e-8) {
        return (f(x + h) - f(x - h)) / (2.0 * h);
    }

    /**
     * @brief Compute second derivative
     *
     * f''(x) ≈ [f(x+h) - 2f(x) + f(x-h)] / h²
     *
     * @param f Function
     * @param x Point
     * @param h Step size (default: 1e-5)
     * @return Approximation of f''(x)
     */
    static double secondDerivative(std::function<double(double)> f,
                                  double x, double h = 1e-5) {
        return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
    }
};

} // namespace maths::calculus

#endif // MATHS_CALCULUS_THEOREMS_HPP
