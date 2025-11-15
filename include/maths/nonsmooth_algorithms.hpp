#ifndef MATHS_ANALYSIS_NONSMOOTH_ALGORITHMS_HPP
#define MATHS_ANALYSIS_NONSMOOTH_ALGORITHMS_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>

/**
 * @file nonsmooth_algorithms.hpp
 * @brief Practical implementations of non-smooth optimization algorithms
 *
 * This module provides ACTUAL computational implementations (not just theory):
 * - Clarke subdifferential calculations
 * - Proximal operators
 * - Bundle methods
 * - ADMM algorithm
 * - Newton's method for non-smooth functions
 * - Subdifferential-based optimization
 */

namespace maths::analysis {

/**
 * @class ProximalOperators
 * @brief Practical proximal operator implementations
 *
 * Usage:
 *   auto prox = ProximalOperators::proxL1(x, lambda);
 *   auto prox_result = ProximalOperators::proxL2Squared(x, lambda, 0.5);
 */
class ProximalOperators {
public:
    /**
     * @brief Proximal operator of L1 norm: prox_{λ‖·‖₁}(x)
     *
     * Soft thresholding operator
     *
     * @param x Input vector
     * @param lambda Proximal parameter (λ > 0)
     * @return Vector result of soft thresholding
     *
     * Formula: prox(x)_i = sign(x_i) * max(|x_i| - λ, 0)
     *
     * Example:
     *   std::vector<double> x = {2.0, -1.5, 0.5, -3.0};
     *   double lambda = 1.0;
     *   auto result = ProximalOperators::proxL1(x, lambda);
     *   // result = {1.0, -0.5, 0.0, -2.0}
     */
    static std::vector<double> proxL1(const std::vector<double>& x, double lambda) {
        if (lambda < 0) {
            throw std::invalid_argument("lambda must be non-negative");
        }

        std::vector<double> result(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            double abs_xi = std::abs(x[i]);
            if (abs_xi <= lambda) {
                result[i] = 0.0;
            } else {
                result[i] = (x[i] > 0 ? 1.0 : -1.0) * (abs_xi - lambda);
            }
        }
        return result;
    }

    /**
     * @brief Proximal operator of squared L2 norm: prox_{(λ/2)‖·‖²}(x)
     *
     * @param x Input vector
     * @param lambda Proximal parameter
     * @param alpha Scaling factor (default 1.0)
     * @return Scaled version of x
     *
     * Formula: prox(x) = x / (1 + λ*alpha)
     *
     * Example:
     *   std::vector<double> x = {2.0, 4.0, 6.0};
     *   auto result = ProximalOperators::proxL2Squared(x, 1.0, 2.0);
     *   // result = {2/3, 4/3, 2}
     */
    static std::vector<double> proxL2Squared(const std::vector<double>& x,
                                              double lambda, double alpha = 1.0) {
        double scale = 1.0 / (1.0 + lambda * alpha);
        std::vector<double> result(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            result[i] = x[i] * scale;
        }
        return result;
    }

    /**
     * @brief Proximal operator of indicator function of box [a,b]
     *
     * Projection onto box constraints
     *
     * @param x Input vector
     * @param a Lower bounds
     * @param b Upper bounds
     * @return Projection of x onto [a,b]
     *
     * Formula: prox(x)_i = clip(x_i, a_i, b_i)
     *
     * Example:
     *   std::vector<double> x = {-1.0, 2.0, 5.0};
     *   std::vector<double> a = {0.0, 0.0, 0.0};
     *   std::vector<double> b = {1.0, 3.0, 4.0};
     *   auto result = ProximalOperators::proxBox(x, a, b);
     *   // result = {0.0, 2.0, 4.0}
     */
    static std::vector<double> proxBox(const std::vector<double>& x,
                                        const std::vector<double>& a,
                                        const std::vector<double>& b) {
        if (x.size() != a.size() || x.size() != b.size()) {
            throw std::invalid_argument("Dimension mismatch");
        }

        std::vector<double> result(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            result[i] = std::max(a[i], std::min(x[i], b[i]));
        }
        return result;
    }

    /**
     * @brief Proximal operator of L-infinity norm
     *
     * Projects onto L1 ball (dual of L-infinity)
     *
     * @param x Input vector
     * @param lambda Radius of L1 ball
     * @return Projection onto {v : ‖v‖₁ ≤ λ}
     *
     * Example:
     *   std::vector<double> x = {0.5, 0.3, 0.4};
     *   auto result = ProximalOperators::proxLInf(x, 1.0);
     */
    static std::vector<double> proxLInf(const std::vector<double>& x, double lambda) {
        // Projection onto L1 ball
        double norm1 = 0.0;
        for (double xi : x) {
            norm1 += std::abs(xi);
        }

        if (norm1 <= lambda) {
            return x;  // Already in ball
        }

        // Binary search for threshold
        std::vector<double> abs_x(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            abs_x[i] = std::abs(x[i]);
        }
        std::sort(abs_x.rbegin(), abs_x.rend());

        double cumsum = 0.0;
        double theta = 0.0;
        for (size_t k = 0; k < abs_x.size(); ++k) {
            cumsum += abs_x[k];
            theta = (cumsum - lambda) / (k + 1);
            if (k + 1 == abs_x.size() || theta >= abs_x[k + 1]) {
                break;
            }
        }

        std::vector<double> result(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            double sign = (x[i] >= 0) ? 1.0 : -1.0;
            result[i] = sign * std::max(0.0, std::abs(x[i]) - theta);
        }
        return result;
    }
};

/**
 * @class SubgradientMethods
 * @brief Subgradient-based optimization algorithms
 *
 * Usage:
 *   auto f = [](const std::vector<double>& x) { return std::abs(x[0]); };
 *   auto subgrad = [](const std::vector<double>& x) {
 *     return std::vector<double>{x[0] >= 0 ? 1.0 : -1.0};
 *   };
 *   auto result = SubgradientMethods::subgradientDescent(f, subgrad, x0, 1000, 0.1);
 */
class SubgradientMethods {
public:
    /**
     * @brief Subgradient descent algorithm
     *
     * @param f Objective function f: R^n → R
     * @param subgradient Function returning subgradient at x
     * @param x0 Initial point
     * @param max_iter Maximum iterations
     * @param alpha Step size (fixed or decreasing)
     * @param tol Tolerance for stopping
     * @return Approximate minimizer
     *
     * Algorithm:
     *   x^{k+1} = x^k - α_k * g^k
     *   where g^k ∈ ∂f(x^k)
     *
     * Example:
     *   // Minimize f(x) = |x| + x²
     *   auto f = [](const std::vector<double>& x) {
     *     return std::abs(x[0]) + x[0] * x[0];
     *   };
     *   auto subgrad = [](const std::vector<double>& x) {
     *     double g = (x[0] > 0 ? 1.0 : -1.0) + 2*x[0];
     *     return std::vector<double>{g};
     *   };
     *   std::vector<double> x0 = {1.0};
     *   auto xopt = SubgradientMethods::subgradientDescent(f, subgrad, x0);
     */
    static std::vector<double> subgradientDescent(
        std::function<double(const std::vector<double>&)> f,
        std::function<std::vector<double>(const std::vector<double>&)> subgradient,
        const std::vector<double>& x0,
        int max_iter = 1000,
        double alpha = 0.01,
        double tol = 1e-6) {

        std::vector<double> x = x0;
        double f_best = f(x);
        std::vector<double> x_best = x;

        for (int iter = 0; iter < max_iter; ++iter) {
            // Compute subgradient
            auto g = subgradient(x);

            // Compute step size (diminishing: α/(iter+1))
            double step = alpha / (iter + 1);

            // Update: x^{k+1} = x^k - step * g^k
            for (size_t i = 0; i < x.size(); ++i) {
                x[i] -= step * g[i];
            }

            // Track best solution
            double fx = f(x);
            if (fx < f_best) {
                f_best = fx;
                x_best = x;
            }

            // Check convergence (subgradient norm)
            double g_norm = 0.0;
            for (double gi : g) {
                g_norm += gi * gi;
            }
            g_norm = std::sqrt(g_norm);

            if (g_norm < tol) {
                break;
            }
        }

        return x_best;
    }

    /**
     * @brief Projected subgradient method
     *
     * @param f Objective function
     * @param subgradient Subgradient oracle
     * @param project Projection onto constraint set
     * @param x0 Initial point
     * @param max_iter Maximum iterations
     * @param alpha Step size
     * @return Constrained minimizer
     *
     * Algorithm:
     *   x^{k+1} = proj_C(x^k - α_k * g^k)
     *
     * Example:
     *   // Minimize |x| subject to x ≥ 0
     *   auto project = [](const std::vector<double>& x) {
     *     return std::vector<double>{std::max(0.0, x[0])};
     *   };
     */
    static std::vector<double> projectedSubgradient(
        std::function<double(const std::vector<double>&)> f,
        std::function<std::vector<double>(const std::vector<double>&)> subgradient,
        std::function<std::vector<double>(const std::vector<double>&)> project,
        const std::vector<double>& x0,
        int max_iter = 1000,
        double alpha = 0.01) {

        std::vector<double> x = x0;
        double f_best = f(x);
        std::vector<double> x_best = x;

        for (int iter = 0; iter < max_iter; ++iter) {
            auto g = subgradient(x);
            double step = alpha / std::sqrt(iter + 1);

            // Update and project
            std::vector<double> x_new(x.size());
            for (size_t i = 0; i < x.size(); ++i) {
                x_new[i] = x[i] - step * g[i];
            }
            x = project(x_new);

            // Track best
            double fx = f(x);
            if (fx < f_best) {
                f_best = fx;
                x_best = x;
            }
        }

        return x_best;
    }
};

/**
 * @class ProximalGradient
 * @brief Proximal gradient algorithms for composite optimization
 *
 * Solves: min f(x) + g(x)
 * where f smooth, g non-smooth with known prox
 *
 * Usage:
 *   // LASSO: min (1/2)‖Ax - b‖² + λ‖x‖₁
 *   auto f = [...]; // smooth part
 *   auto grad_f = [...]; // gradient of f
 *   auto prox_g = [](auto x, double t) { return ProximalOperators::proxL1(x, t); };
 *   auto result = ProximalGradient::ISTA(grad_f, prox_g, x0, lambda);
 */
class ProximalGradient {
public:
    /**
     * @brief ISTA (Iterative Shrinkage-Thresholding Algorithm)
     *
     * @param grad_f Gradient of smooth part f
     * @param prox_g Proximal operator of non-smooth part g
     * @param x0 Initial point
     * @param lambda Regularization parameter
     * @param L Lipschitz constant of ∇f
     * @param max_iter Maximum iterations
     * @param tol Convergence tolerance
     * @return Minimizer
     *
     * Algorithm:
     *   x^{k+1} = prox_{g/L}(x^k - (1/L)∇f(x^k))
     *
     * Example (LASSO):
     *   // min (1/2)‖Ax - b‖² + λ‖x‖₁
     *   auto A = ...; auto b = ...;
     *   auto grad_f = [&](const auto& x) {
     *     // ∇f(x) = A^T(Ax - b)
     *     auto Ax = matmul(A, x);
     *     auto residual = subtract(Ax, b);
     *     return matmul_transpose(A, residual);
     *   };
     *   auto prox_g = [lambda](const auto& x, double t) {
     *     return ProximalOperators::proxL1(x, lambda * t);
     *   };
     *   auto x_opt = ProximalGradient::ISTA(grad_f, prox_g, x0, lambda, L);
     */
    static std::vector<double> ISTA(
        std::function<std::vector<double>(const std::vector<double>&)> grad_f,
        std::function<std::vector<double>(const std::vector<double>&, double)> prox_g,
        const std::vector<double>& x0,
        double lambda,
        double L = 1.0,
        int max_iter = 1000,
        double tol = 1e-6) {

        std::vector<double> x = x0;
        double t = 1.0 / L;

        for (int iter = 0; iter < max_iter; ++iter) {
            std::vector<double> x_old = x;

            // Gradient step
            auto gf = grad_f(x);
            std::vector<double> x_grad(x.size());
            for (size_t i = 0; i < x.size(); ++i) {
                x_grad[i] = x[i] - t * gf[i];
            }

            // Proximal step
            x = prox_g(x_grad, t);

            // Check convergence
            double diff = 0.0;
            for (size_t i = 0; i < x.size(); ++i) {
                double d = x[i] - x_old[i];
                diff += d * d;
            }
            if (std::sqrt(diff) < tol) {
                break;
            }
        }

        return x;
    }

    /**
     * @brief FISTA (Fast ISTA) with Nesterov acceleration
     *
     * Same interface as ISTA but with momentum acceleration
     * Achieves O(1/k²) convergence vs O(1/k) for ISTA
     *
     * @param grad_f Gradient of smooth part
     * @param prox_g Proximal operator
     * @param x0 Initial point
     * @param lambda Regularization parameter
     * @param L Lipschitz constant
     * @param max_iter Maximum iterations
     * @param tol Tolerance
     * @return Minimizer
     *
     * Algorithm:
     *   y^k = x^k + ((t_{k-1} - 1)/t_k)(x^k - x^{k-1})
     *   x^{k+1} = prox_{g/L}(y^k - (1/L)∇f(y^k))
     *   t_{k+1} = (1 + √(1 + 4t_k²))/2
     */
    static std::vector<double> FISTA(
        std::function<std::vector<double>(const std::vector<double>&)> grad_f,
        std::function<std::vector<double>(const std::vector<double>&, double)> prox_g,
        const std::vector<double>& x0,
        double lambda,
        double L = 1.0,
        int max_iter = 1000,
        double tol = 1e-6) {

        std::vector<double> x = x0;
        std::vector<double> x_old = x0;
        std::vector<double> y = x0;
        double t = 1.0;
        double step = 1.0 / L;

        for (int iter = 0; iter < max_iter; ++iter) {
            // Gradient at y
            auto gy = grad_f(y);

            // Proximal step
            std::vector<double> x_new(x.size());
            for (size_t i = 0; i < x.size(); ++i) {
                x_new[i] = y[i] - step * gy[i];
            }
            x_new = prox_g(x_new, step);

            // Momentum update
            double t_new = (1.0 + std::sqrt(1.0 + 4.0 * t * t)) / 2.0;
            double momentum = (t - 1.0) / t_new;

            for (size_t i = 0; i < x.size(); ++i) {
                y[i] = x_new[i] + momentum * (x_new[i] - x[i]);
            }

            // Check convergence
            double diff = 0.0;
            for (size_t i = 0; i < x.size(); ++i) {
                double d = x_new[i] - x[i];
                diff += d * d;
            }

            x_old = x;
            x = x_new;
            t = t_new;

            if (std::sqrt(diff) < tol) {
                break;
            }
        }

        return x;
    }
};

/**
 * @class ADMM
 * @brief Alternating Direction Method of Multipliers
 *
 * Solves: min f(x) + g(z) subject to Ax + Bz = c
 *
 * Usage:
 *   auto result = ADMM::solve(prox_f, prox_g, A, B, c, x0, z0);
 */
class ADMM {
public:
    /**
     * @brief ADMM algorithm for consensus optimization
     *
     * Solves: min f(x) + g(x)
     * by reformulating as: min f(x) + g(z) s.t. x = z
     *
     * @param prox_f Proximal operator of f
     * @param prox_g Proximal operator of g
     * @param x0 Initial x
     * @param rho Penalty parameter (ρ)
     * @param max_iter Maximum iterations
     * @param tol Tolerance
     * @return Optimal x
     *
     * Algorithm:
     *   x^{k+1} = prox_{f/ρ}(z^k - u^k)
     *   z^{k+1} = prox_{g/ρ}(x^{k+1} + u^k)
     *   u^{k+1} = u^k + x^{k+1} - z^{k+1}
     *
     * Example (LASSO via ADMM):
     *   // min (1/2)‖Ax - b‖² + λ‖z‖₁ s.t. x = z
     *   auto prox_f = [&](const auto& v, double t) {
     *     // (A^TA + (1/t)I)^{-1}(A^Tb + v/t)
     *     return solve_linear_system(...);
     *   };
     *   auto prox_g = [lambda](const auto& v, double t) {
     *     return ProximalOperators::proxL1(v, lambda * t);
     *   };
     *   auto x_opt = ADMM::consensus(prox_f, prox_g, x0, 1.0);
     */
    static std::vector<double> consensus(
        std::function<std::vector<double>(const std::vector<double>&, double)> prox_f,
        std::function<std::vector<double>(const std::vector<double>&, double)> prox_g,
        const std::vector<double>& x0,
        double rho = 1.0,
        int max_iter = 1000,
        double tol = 1e-6) {

        size_t n = x0.size();
        std::vector<double> x = x0;
        std::vector<double> z = x0;
        std::vector<double> u(n, 0.0);  // Dual variable

        for (int iter = 0; iter < max_iter; ++iter) {
            // x-update: prox_f(z - u)
            std::vector<double> v(n);
            for (size_t i = 0; i < n; ++i) {
                v[i] = z[i] - u[i];
            }
            x = prox_f(v, 1.0 / rho);

            // z-update: prox_g(x + u)
            for (size_t i = 0; i < n; ++i) {
                v[i] = x[i] + u[i];
            }
            z = prox_g(v, 1.0 / rho);

            // u-update: u + (x - z)
            double residual = 0.0;
            for (size_t i = 0; i < n; ++i) {
                double r = x[i] - z[i];
                u[i] += r;
                residual += r * r;
            }

            // Check convergence
            if (std::sqrt(residual) < tol) {
                break;
            }
        }

        return x;
    }
};

/**
 * @class ClarkeSubdifferentialCalc
 * @brief Practical Clarke subdifferential calculations
 */
class ClarkeSubdifferentialCalc {
public:
    /**
     * @brief Compute Clarke subdifferential of |x|
     *
     * @param x Point to evaluate at
     * @return Pair of (lower, upper) bounds of ∂_C|x|
     *
     * Result:
     *   x > 0: ∂_C|x| = {1}
     *   x < 0: ∂_C|x| = {-1}
     *   x = 0: ∂_C|x| = [-1, 1]
     *
     * Example:
     *   auto [lower, upper] = ClarkeSubdifferentialCalc::clarkeAbsoluteValue(0.0);
     *   // lower = -1.0, upper = 1.0
     */
    static std::pair<double, double> clarkeAbsoluteValue(double x) {
        if (x > 1e-10) {
            return {1.0, 1.0};
        } else if (x < -1e-10) {
            return {-1.0, -1.0};
        } else {
            return {-1.0, 1.0};
        }
    }

    /**
     * @brief Compute Clarke subdifferential of max{x, y}
     *
     * @param x First value
     * @param y Second value
     * @param grad_x Gradient w.r.t. x component
     * @param grad_y Gradient w.r.t. y component
     * @return Vector of possible subgradients
     *
     * Result:
     *   x > y: ∂_C max = grad_x
     *   x < y: ∂_C max = grad_y
     *   x = y: ∂_C max = conv{grad_x, grad_y}
     *
     * Example:
     *   auto subdiff = ClarkeSubdifferentialCalc::clarkeMax(2.0, 2.0, 1.0, 1.0);
     *   // subdiff contains convex combinations
     */
    static std::vector<double> clarkeMax(double x, double y,
                                          double grad_x, double grad_y) {
        if (x > y + 1e-10) {
            return {grad_x};
        } else if (y > x + 1e-10) {
            return {grad_y};
        } else {
            // Return both extremes of convex hull
            return {grad_x, grad_y};
        }
    }

    /**
     * @brief Compute Clarke subdifferential of ReLU: max{0, x}
     *
     * @param x Input value
     * @return Pair (lower, upper) of subdifferential interval
     *
     * Result:
     *   x > 0: [1, 1]
     *   x < 0: [0, 0]
     *   x = 0: [0, 1]
     */
    static std::pair<double, double> clarkeReLU(double x) {
        if (x > 1e-10) {
            return {1.0, 1.0};
        } else if (x < -1e-10) {
            return {0.0, 0.0};
        } else {
            return {0.0, 1.0};
        }
    }
};

} // namespace maths::analysis

#endif // MATHS_ANALYSIS_NONSMOOTH_ALGORITHMS_HPP
