/**
 * @file ode_dynamical_systems.hpp
 * @brief Comprehensive Differential Equations and Dynamical Systems Theory
 *
 * PART 1: CLASSICAL THEORY
 * Chapter 1: Introduction (Newton's equations, qualitative analysis, explicit solutions)
 * Chapter 2: Initial Value Problems (existence, uniqueness, Picard iteration)
 * Chapter 3: Linear Equations (linear systems, Floquet theory, periodic systems)
 * Chapter 4: Differential Equations in Complex Domain (Frobenius method)
 * Chapter 5: Boundary Value Problems (Sturm-Liouville, oscillation theory)
 *
 * PART 2: DYNAMICAL SYSTEMS
 * Chapter 6: Dynamical Systems (flows, orbits, stability, Liapunov)
 * Chapter 7: Local Behavior near Fixed Points (stable/unstable manifolds, Hartman-Grobman)
 * Chapter 8: Planar Dynamical Systems (Poincaré-Bendixson, limit cycles)
 * Chapter 9: Higher Dimensional Systems (attractors, Lorenz, Hamiltonian, KAM)
 *
 * PART 3: CHAOS
 * Chapter 10: Discrete Dynamical Systems (logistic map, period doubling)
 * Chapter 11: Periodic Solutions (Poincaré map, Melnikov method)
 * Chapter 12: 1D Chaos (Sarkovskii, symbolic dynamics, fractals, homoclinic orbits)
 * Chapter 13: Higher Dimensional Chaos (Smale horseshoe, homoclinic theorem)
 */

#ifndef MATHS_DYNAMICAL_SYSTEMS_ODE_HPP
#define MATHS_DYNAMICAL_SYSTEMS_ODE_HPP

#include <vector>
#include <cmath>
#include <functional>
#include <complex>
#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <numeric>

namespace maths::dynamical_systems {

// Type aliases
using VectorField = std::function<std::vector<double>(double, const std::vector<double>&)>;
using ScalarODE = std::function<double(double, double)>;
using ComplexFunction = std::function<std::complex<double>(std::complex<double>)>;

/**
 * ============================================================================
 * PART 1: CLASSICAL THEORY - Chapters 1-5
 * ============================================================================
 */

/**
 * @class NewtonEquations
 * @brief Chapter 1: Newton's second law and first-order systems
 */
class NewtonEquations {
public:
    /**
     * @brief Convert second-order ODE to first-order system
     *
     * x'' = f(t, x, x') becomes:
     * y1' = y2
     * y2' = f(t, y1, y2)
     */
    static VectorField toFirstOrderSystem(
        std::function<double(double, double, double)> f) {

        return [f](double t, const std::vector<double>& y) {
            // y = [x, x']
            return std::vector<double>{y[1], f(t, y[0], y[1])};
        };
    }

    /**
     * @brief Classify ODE: linear, separable, exact, etc.
     */
    enum class ODEType {
        LINEAR,
        SEPARABLE,
        EXACT,
        BERNOULLI,
        RICATTI,
        GENERAL
    };

    /**
     * @brief Check if first-order ODE is separable: dy/dx = g(x)h(y)
     */
    static bool isSeparable(ScalarODE f, double x, double y) {
        // Heuristic check
        double h = 1e-6;
        double fx = f(x + h, y) - f(x, y);
        double fy = f(x, y + h) - f(x, y);

        // If ∂f/∂x * y ≈ ∂f/∂y * x, likely separable
        return std::abs(fx * y - fy * x) < 1e-4;
    }

    /**
     * @brief Autonomous equation: x' = f(x)
     */
    static std::vector<double> findEquilibria(
        std::function<double(double)> f,
        double a, double b, int n_points = 100) {

        std::vector<double> equilibria;
        double dx = (b - a) / n_points;

        for (int i = 0; i < n_points - 1; ++i) {
            double x1 = a + i * dx;
            double x2 = x1 + dx;
            double f1 = f(x1);
            double f2 = f(x2);

            // Sign change indicates equilibrium
            if (f1 * f2 < 0) {
                // Bisection to refine
                double tol = 1e-8;
                while (x2 - x1 > tol) {
                    double xm = 0.5 * (x1 + x2);
                    double fm = f(xm);
                    if (fm * f1 < 0) {
                        x2 = xm;
                        f2 = fm;
                    } else {
                        x1 = xm;
                        f1 = fm;
                    }
                }
                equilibria.push_back(0.5 * (x1 + x2));
            }
        }

        return equilibria;
    }
};

/**
 * @class InitialValueProblem
 * @brief Chapter 2: Existence, uniqueness, Picard iteration, Peano theorem
 */
class InitialValueProblem {
public:
    /**
     * @brief Euler's method: y_{n+1} = y_n + h*f(t_n, y_n)
     */
    static std::vector<std::pair<double, double>> euler(
        ScalarODE f, double t0, double y0, double T, int n_steps) {

        std::vector<std::pair<double, double>> solution;
        solution.reserve(n_steps + 1);

        double h = (T - t0) / n_steps;
        double t = t0, y = y0;

        solution.push_back({t, y});

        for (int i = 0; i < n_steps; ++i) {
            y += h * f(t, y);
            t += h;
            solution.push_back({t, y});
        }

        return solution;
    }

    /**
     * @brief Improved Euler (Heun's method)
     */
    static std::vector<std::pair<double, double>> heun(
        ScalarODE f, double t0, double y0, double T, int n_steps) {

        std::vector<std::pair<double, double>> solution;
        double h = (T - t0) / n_steps;
        double t = t0, y = y0;

        solution.push_back({t, y});

        for (int i = 0; i < n_steps; ++i) {
            double k1 = f(t, y);
            double k2 = f(t + h, y + h * k1);
            y += 0.5 * h * (k1 + k2);
            t += h;
            solution.push_back({t, y});
        }

        return solution;
    }

    /**
     * @brief Runge-Kutta 4th order (RK4)
     */
    static std::vector<std::pair<double, double>> rk4(
        ScalarODE f, double t0, double y0, double T, int n_steps) {

        std::vector<std::pair<double, double>> solution;
        double h = (T - t0) / n_steps;
        double t = t0, y = y0;

        solution.push_back({t, y});

        for (int i = 0; i < n_steps; ++i) {
            double k1 = f(t, y);
            double k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
            double k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
            double k4 = f(t + h, y + h * k3);

            y += (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
            t += h;
            solution.push_back({t, y});
        }

        return solution;
    }

    /**
     * @brief Picard iteration for existence/uniqueness proof
     *
     * φ_{n+1}(t) = y0 + ∫_{t0}^t f(s, φ_n(s)) ds
     */
    static std::vector<double> picard Iteration(
        ScalarODE f, double t0, double y0,
        const std::vector<double>& t_points, int iterations) {

        std::vector<double> phi(t_points.size(), y0);

        for (int iter = 0; iter < iterations; ++iter) {
            std::vector<double> phi_new(t_points.size());
            phi_new[0] = y0;

            for (size_t i = 1; i < t_points.size(); ++i) {
                // Trapezoidal integration
                double integral = 0.0;
                for (size_t j = 0; j < i; ++j) {
                    double dt = t_points[j + 1] - t_points[j];
                    integral += 0.5 * dt * (
                        f(t_points[j], phi[j]) +
                        f(t_points[j + 1], phi[j + 1])
                    );
                }
                phi_new[i] = y0 + integral;
            }

            phi = phi_new;
        }

        return phi;
    }

    /**
     * @brief Check Lipschitz continuity: |f(t,y1) - f(t,y2)| ≤ L|y1-y2|
     */
    static double estimateLipschitzConstant(
        ScalarODE f, double t, double y_min, double y_max, int n_samples = 100) {

        double max_L = 0.0;
        double dy = (y_max - y_min) / n_samples;

        for (int i = 0; i < n_samples; ++i) {
            for (int j = i + 1; j < n_samples; ++j) {
                double y1 = y_min + i * dy;
                double y2 = y_min + j * dy;

                double L = std::abs(f(t, y1) - f(t, y2)) / std::abs(y1 - y2);
                max_L = std::max(max_L, L);
            }
        }

        return max_L;
    }
};

/**
 * @class LinearSystems
 * @brief Chapter 3: Linear first-order systems, fundamental matrices
 */
class LinearSystems {
public:
    /**
     * @brief Solve x' = Ax for constant matrix A
     *
     * Solution: x(t) = e^(tA) x0
     */
    static std::vector<double> solveConstantCoefficients(
        const std::vector<std::vector<double>>& A,
        const std::vector<double>& x0,
        double t) {

        int n = A.size();

        // Compute matrix exponential e^(tA) via series (simplified)
        std::vector<std::vector<double>> exp_tA(n, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> power_tA(n, std::vector<double>(n));

        // Initialize to identity
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                exp_tA[i][j] = (i == j) ? 1.0 : 0.0;
                power_tA[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        // Series: e^(tA) = I + tA + (tA)²/2! + ...
        double factorial = 1.0;
        for (int k = 1; k <= 20; ++k) {  // 20 terms
            factorial *= k;

            // power_tA = power_tA * (tA)
            std::vector<std::vector<double>> temp(n, std::vector<double>(n, 0.0));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    for (int m = 0; m < n; ++m) {
                        temp[i][j] += power_tA[i][m] * t * A[m][j];
                    }
                }
            }
            power_tA = temp;

            // Add to exponential
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    exp_tA[i][j] += power_tA[i][j] / factorial;
                }
            }
        }

        // Multiply by x0
        std::vector<double> result(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result[i] += exp_tA[i][j] * x0[j];
            }
        }

        return result;
    }

    /**
     * @brief Floquet theory for periodic systems: x' = A(t)x, A(t+T) = A(t)
     *
     * Monodromy matrix M and Floquet multipliers
     */
    struct FloquetTheory {
        std::vector<std::vector<double>> monodromy_matrix;
        std::vector<std::complex<double>> multipliers;

        bool isStable() const {
            for (const auto& mu : multipliers) {
                if (std::abs(mu) >= 1.0) return false;
            }
            return true;
        }
    };

    /**
     * @brief Compute fundamental matrix Φ(t) satisfying Φ' = A(t)Φ, Φ(0) = I
     */
    static std::vector<std::vector<double>> fundamentalMatrix(
        std::function<std::vector<std::vector<double>>(double)> A,
        double t, int n) {

        // Use matrix ODE solver (simplified)
        std::vector<std::vector<double>> Phi(n, std::vector<double>(n));
        for (int i = 0; i < n; ++i) {
            Phi[i][i] = 1.0;  // Identity initial condition
        }

        return Phi;  // Placeholder
    }
};

/**
 * @class ComplexODE
 * @brief Chapter 4: Differential equations in complex domain, Frobenius method
 */
class ComplexODE {
public:
    /**
     * @brief Frobenius method for second-order ODEs with regular singular points
     *
     * Solve: z²w'' + p(z)zw' + q(z)w = 0
     * near z = 0 (regular singular point)
     *
     * Frobenius series: w(z) = z^r ∑ aₙ z^n
     */
    struct FrobeniusSolution {
        std::vector<double> indicial_roots;  // Roots of indicial equation
        std::vector<std::vector<double>> coefficients;  // Series coefficients

        std::complex<double> evaluate(std::complex<double> z, int series_idx) const {
            if (series_idx >= static_cast<int>(indicial_roots.size())) {
                throw std::out_of_range("Invalid series index");
            }

            double r = indicial_roots[series_idx];
            const auto& a = coefficients[series_idx];

            std::complex<double> sum = 0.0;
            std::complex<double> z_power = std::pow(z, r);

            for (size_t n = 0; n < a.size(); ++n) {
                sum += a[n] * z_power;
                z_power *= z;
            }

            return sum;
        }
    };

    /**
     * @brief Find indicial equation at regular singular point
     *
     * For z²w'' + p₀w' + q₀w = 0:
     * Indicial equation: r(r-1) + p₀r + q₀ = 0
     */
    static std::vector<double> indicialEquation(double p0, double q0) {
        // r² + (p0-1)r + q0 = 0
        double a = 1.0;
        double b = p0 - 1.0;
        double c = q0;

        double discriminant = b * b - 4 * a * c;

        if (discriminant >= 0) {
            double r1 = (-b + std::sqrt(discriminant)) / (2 * a);
            double r2 = (-b - std::sqrt(discriminant)) / (2 * a);
            return {r1, r2};
        } else {
            // Complex roots
            double real = -b / (2 * a);
            double imag = std::sqrt(-discriminant) / (2 * a);
            return {real, imag};  // Store as [Re, Im] for complex conjugate pair
        }
    }

    /**
     * @brief Bessel's equation: z²w'' + zw' + (z² - ν²)w = 0
     */
    static FrobeniusSolution besselEquation(double nu, int n_terms = 10) {
        FrobeniusSolution sol;
        sol.indicial_roots = {nu, -nu};
        sol.coefficients.resize(2);

        // First solution with r = ν
        sol.coefficients[0].resize(n_terms);
        sol.coefficients[0][0] = 1.0;  // a₀ = 1

        for (int n = 1; n < n_terms; ++n) {
            // Recurrence: aₙ = -aₙ₋₂ / (n(n + 2ν))
            if (n >= 2) {
                sol.coefficients[0][n] = -sol.coefficients[0][n - 2] /
                    (n * (n + 2 * nu));
            }
        }

        return sol;
    }
};

/**
 * ============================================================================
 * PART 2: DYNAMICAL SYSTEMS - Chapters 6-9
 * ============================================================================
 */

/**
 * @class DynamicalSystem
 * @brief Chapter 6: Flows, orbits, invariant sets, stability, Liapunov functions
 */
class DynamicalSystem {
public:
    /**
     * @brief Compute flow φ_t(x0): solution starting at x0 at time t
     */
    static std::vector<double> flow(
        VectorField f, const std::vector<double>& x0, double t, double dt = 0.01) {

        std::vector<double> x = x0;
        int n_steps = static_cast<int>(std::abs(t) / dt);
        double step = (t > 0) ? dt : -dt;

        for (int i = 0; i < n_steps; ++i) {
            // RK4 step
            auto k1 = f(i * step, x);

            std::vector<double> x_temp(x.size());
            for (size_t j = 0; j < x.size(); ++j) {
                x_temp[j] = x[j] + 0.5 * step * k1[j];
            }
            auto k2 = f((i + 0.5) * step, x_temp);

            for (size_t j = 0; j < x.size(); ++j) {
                x_temp[j] = x[j] + 0.5 * step * k2[j];
            }
            auto k3 = f((i + 0.5) * step, x_temp);

            for (size_t j = 0; j < x.size(); ++j) {
                x_temp[j] = x[j] + step * k3[j];
            }
            auto k4 = f((i + 1) * step, x_temp);

            // Update
            for (size_t j = 0; j < x.size(); ++j) {
                x[j] += (step / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
            }
        }

        return x;
    }

    /**
     * @brief Find fixed points: f(x*) = 0
     */
    static std::vector<std::vector<double>> findFixedPoints(
        VectorField f,
        const std::vector<double>& x_min,
        const std::vector<double>& x_max,
        int n_grid = 10) {

        std::vector<std::vector<double>> fixed_points;

        // Grid search followed by Newton refinement
        int dim = x_min.size();
        std::vector<double> dx(dim);
        for (int i = 0; i < dim; ++i) {
            dx[i] = (x_max[i] - x_min[i]) / n_grid;
        }

        // Simplified 2D case
        if (dim == 2) {
            for (int i = 0; i < n_grid; ++i) {
                for (int j = 0; j < n_grid; ++j) {
                    std::vector<double> x = {
                        x_min[0] + i * dx[0],
                        x_min[1] + j * dx[1]
                    };

                    auto fx = f(0, x);
                    double norm = std::sqrt(fx[0] * fx[0] + fx[1] * fx[1]);

                    if (norm < 0.1) {  // Potential fixed point
                        // Newton refinement
                        for (int iter = 0; iter < 10; ++iter) {
                            fx = f(0, x);
                            // Simplified: x := x - f(x) (steepest descent)
                            x[0] -= 0.1 * fx[0];
                            x[1] -= 0.1 * fx[1];
                        }

                        fx = f(0, x);
                        norm = std::sqrt(fx[0] * fx[0] + fx[1] * fx[1]);
                        if (norm < 1e-6) {
                            fixed_points.push_back(x);
                        }
                    }
                }
            }
        }

        return fixed_points;
    }

    /**
     * @brief Liapunov function for stability analysis
     *
     * V(x) > 0 for x ≠ 0, V(0) = 0
     * V'(x) = ∇V · f(x) < 0 => asymptotic stability
     */
    struct LiapunovFunction {
        std::function<double(std::vector<double>)> V;
        std::function<std::vector<double>(std::vector<double>)> grad_V;

        /**
         * @brief Check if V is Liapunov function for system f
         */
        bool verifyStability(VectorField f, const std::vector<double>& x) const {
            auto grad = grad_V(x);
            auto fx = f(0, x);

            double V_dot = 0.0;
            for (size_t i = 0; i < x.size(); ++i) {
                V_dot += grad[i] * fx[i];
            }

            return V_dot < 0;
        }
    };

    /**
     * @brief Quadratic Liapunov function: V(x) = x^T P x
     */
    static LiapunovFunction quadraticLiapunov(
        const std::vector<std::vector<double>>& P) {

        LiapunovFunction lyap;

        lyap.V = [P](const std::vector<double>& x) {
            double result = 0.0;
            for (size_t i = 0; i < x.size(); ++i) {
                for (size_t j = 0; j < x.size(); ++j) {
                    result += x[i] * P[i][j] * x[j];
                }
            }
            return result;
        };

        lyap.grad_V = [P](const std::vector<double>& x) {
            std::vector<double> grad(x.size(), 0.0);
            for (size_t i = 0; i < x.size(); ++i) {
                for (size_t j = 0; j < x.size(); ++j) {
                    grad[i] += 2.0 * P[i][j] * x[j];
                }
            }
            return grad;
        };

        return lyap;
    }
};

/**
 * @class LocalBehavior
 * @brief Chapter 7: Stable/unstable manifolds, Hartman-Grobman theorem, linearization
 */
class LocalBehavior {
public:
    /**
     * @brief Classify fixed point via eigenvalues of Jacobian
     */
    enum class FixedPointType {
        STABLE_NODE,      // All eigenvalues Re(λ) < 0
        UNSTABLE_NODE,    // All eigenvalues Re(λ) > 0
        SADDLE,           // Mixed signs
        STABLE_SPIRAL,    // Complex with Re(λ) < 0
        UNSTABLE_SPIRAL,  // Complex with Re(λ) > 0
        CENTER           // Pure imaginary eigenvalues
    };

    /**
     * @brief Compute Jacobian matrix Df at point x
     */
    static std::vector<std::vector<double>> jacobian(
        VectorField f, const std::vector<double>& x, double t = 0.0, double h = 1e-6) {

        int n = x.size();
        std::vector<std::vector<double>> J(n, std::vector<double>(n));

        auto fx = f(t, x);

        for (int j = 0; j < n; ++j) {
            std::vector<double> x_pert = x;
            x_pert[j] += h;
            auto fx_pert = f(t, x_pert);

            for (int i = 0; i < n; ++i) {
                J[i][j] = (fx_pert[i] - fx[i]) / h;
            }
        }

        return J;
    }

    /**
     * @brief Compute eigenvalues (2x2 case)
     */
    static std::vector<std::complex<double>> eigenvalues2x2(
        const std::vector<std::vector<double>>& A) {

        // Characteristic polynomial: λ² - tr(A)λ + det(A) = 0
        double trace = A[0][0] + A[1][1];
        double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];

        double discriminant = trace * trace - 4 * det;

        if (discriminant >= 0) {
            double lambda1 = 0.5 * (trace + std::sqrt(discriminant));
            double lambda2 = 0.5 * (trace - std::sqrt(discriminant));
            return {lambda1, lambda2};
        } else {
            double real = 0.5 * trace;
            double imag = 0.5 * std::sqrt(-discriminant);
            return {
                std::complex<double>(real, imag),
                std::complex<double>(real, -imag)
            };
        }
    }

    /**
     * @brief Classify 2D fixed point
     */
    static FixedPointType classify2DFixedPoint(
        const std::vector<std::vector<double>>& jacobian) {

        auto eigenvals = eigenvalues2x2(jacobian);

        if (std::abs(eigenvals[0].imag()) < 1e-10) {
            // Real eigenvalues
            double lambda1 = eigenvals[0].real();
            double lambda2 = eigenvals[1].real();

            if (lambda1 < 0 && lambda2 < 0) return FixedPointType::STABLE_NODE;
            if (lambda1 > 0 && lambda2 > 0) return FixedPointType::UNSTABLE_NODE;
            if (lambda1 * lambda2 < 0) return FixedPointType::SADDLE;
        } else {
            // Complex eigenvalues
            double real_part = eigenvals[0].real();

            if (real_part < 0) return FixedPointType::STABLE_SPIRAL;
            if (real_part > 0) return FixedPointType::UNSTABLE_SPIRAL;
            return FixedPointType::CENTER;
        }

        return FixedPointType::CENTER;
    }

    /**
     * @brief Hartman-Grobman theorem: Near hyperbolic fixed point,
     * nonlinear system is topologically conjugate to its linearization
     */
    static bool isHyperbolic(const std::vector<std::complex<double>>& eigenvalues) {
        // Hyperbolic if no eigenvalue has Re(λ) = 0
        for (const auto& lambda : eigenvalues) {
            if (std::abs(lambda.real()) < 1e-10) {
                return false;
            }
        }
        return true;
    }
};

/**
 * ============================================================================
 * PART 3: CHAOS - Chapters 10-13
 * ============================================================================
 */

/**
 * @class DiscreteDynamicalSystems
 * @brief Chapter 10: Logistic map, period doubling, bifurcations
 */
class DiscreteDynamicalSystems {
public:
    /**
     * @brief Logistic map: xₙ₊₁ = rxₙ(1 - xₙ)
     */
    static std::vector<double> logisticMap(double r, double x0, int iterations) {
        std::vector<double> orbit;
        orbit.reserve(iterations);

        double x = x0;
        for (int i = 0; i < iterations; ++i) {
            orbit.push_back(x);
            x = r * x * (1 - x);
        }

        return orbit;
    }

    /**
     * @brief Find fixed points: x* = f(x*)
     */
    static std::vector<double> findFixedPoints(
        std::function<double(double)> f, double a, double b, int n_grid = 100) {

        std::vector<double> fixed_points;
        double dx = (b - a) / n_grid;

        for (int i = 0; i < n_grid; ++i) {
            double x = a + i * dx;
            double fx = f(x);

            // Check if f(x) ≈ x
            if (std::abs(fx - x) < 1e-6) {
                fixed_points.push_back(x);
            }
        }

        return fixed_points;
    }

    /**
     * @brief Period-n orbits: fⁿ(x) = x
     */
    static std::vector<std::vector<double>> findPeriodicOrbits(
        std::function<double(double)> f, int period, double a, double b) {

        std::vector<std::vector<double>> orbits;

        // Compose f with itself period times
        auto f_n = [f, period](double x) {
            double result = x;
            for (int i = 0; i < period; ++i) {
                result = f(result);
            }
            return result;
        };

        // Find fixed points of f^n
        auto points = findFixedPoints(f_n, a, b);

        // Group into orbits
        for (double x : points) {
            std::vector<double> orbit;
            double current = x;
            for (int i = 0; i < period; ++i) {
                orbit.push_back(current);
                current = f(current);
            }
            orbits.push_back(orbit);
        }

        return orbits;
    }

    /**
     * @brief Lyapunov exponent: λ = lim (1/n) ∑ log|f'(xᵢ)|
     *
     * λ > 0: chaos
     * λ = 0: marginal
     * λ < 0: stable
     */
    static double lyapunovExponent(
        std::function<double(double)> f,
        std::function<double(double)> df,
        double x0, int iterations, int transient = 100) {

        double x = x0;

        // Discard transient
        for (int i = 0; i < transient; ++i) {
            x = f(x);
        }

        // Compute Lyapunov exponent
        double sum = 0.0;
        for (int i = 0; i < iterations; ++i) {
            double derivative = df(x);
            if (std::abs(derivative) > 1e-10) {
                sum += std::log(std::abs(derivative));
            }
            x = f(x);
        }

        return sum / iterations;
    }

    /**
     * @brief Bifurcation diagram
     */
    struct BifurcationDiagram {
        std::vector<double> parameters;
        std::vector<std::vector<double>> attractors;
    };

    /**
     * @brief Compute bifurcation diagram for logistic map
     */
    static BifurcationDiagram logisticBifurcation(
        double r_min, double r_max, int n_params = 1000) {

        BifurcationDiagram diagram;
        diagram.parameters.reserve(n_params);
        diagram.attractors.reserve(n_params);

        double dr = (r_max - r_min) / n_params;

        for (int i = 0; i < n_params; ++i) {
            double r = r_min + i * dr;
            diagram.parameters.push_back(r);

            // Iterate to attractor
            double x = 0.5;
            for (int j = 0; j < 500; ++j) {  // Transient
                x = r * x * (1 - x);
            }

            // Collect points on attractor
            std::set<double> attractor_set;
            for (int j = 0; j < 200; ++j) {
                x = r * x * (1 - x);
                attractor_set.insert(x);
            }

            diagram.attractors.push_back(
                std::vector<double>(attractor_set.begin(), attractor_set.end())
            );
        }

        return diagram;
    }
};

/**
 * @class ChaoticSystems
 * @brief Chapters 11-13: Poincaré maps, Melnikov method, symbolic dynamics, fractals
 */
class ChaoticSystems {
public:
    /**
     * @brief Poincaré map for periodic systems: Σ → Σ
     *
     * First return map to surface of section
     */
    static std::vector<double> poincareMap(
        VectorField f,
        const std::vector<double>& x0,
        std::function<bool(std::vector<double>)> surface,
        double period_estimate) {

        // Integrate until surface is crossed again
        std::vector<double> x = x0;
        double t = 0.0;
        double dt = period_estimate / 100.0;

        // Skip first crossing
        bool crossed = false;
        while (t < 10 * period_estimate) {
            auto x_next = DynamicalSystem::flow(f, x, dt);

            bool on_surface = surface(x_next);

            if (on_surface && crossed) {
                return x_next;
            }

            if (on_surface) {
                crossed = true;
            }

            x = x_next;
            t += dt;
        }

        return x;  // Failed to find return
    }

    /**
     * @brief Melnikov method for detecting homoclinic chaos
     *
     * M(t0) = ∫_{-∞}^{∞} f0(γ(t)) ∧ f1(γ(t), t+t0) dt
     *
     * Simple zeros => transverse homoclinic intersection => chaos
     */
    static double melnikovFunction(
        VectorField f0,  // Unperturbed system
        VectorField f1,  // Perturbation
        const std::vector<std::vector<double>>& homoclinic_orbit,
        const std::vector<double>& times,
        double t0) {

        double integral = 0.0;

        for (size_t i = 0; i < homoclinic_orbit.size() - 1; ++i) {
            double t = times[i] + t0;
            double dt = times[i + 1] - times[i];

            auto x = homoclinic_orbit[i];
            auto f0_val = f0(times[i], x);
            auto f1_val = f1(t, x);

            // Wedge product in 2D: f0 ∧ f1 = f0_x * f1_y - f0_y * f1_x
            double wedge = f0_val[0] * f1_val[1] - f0_val[1] * f1_val[0];

            integral += wedge * dt;
        }

        return integral;
    }

    /**
     * @brief Sarkovskii's theorem ordering for periods
     *
     * 3 > 5 > 7 > ... > 2·3 > 2·5 > ... > 2²·3 > ... > 2³ > 2² > 2 > 1
     *
     * If f has period n and n > m in Sarkovskii ordering, then f has period m
     */
    static bool sarkovskiiGreater(int n, int m) {
        // Simplified: 3 implies all periods
        if (n == 3 && m != 3) return true;

        // Powers of 2 ordering
        bool n_power_of_2 = (n & (n - 1)) == 0;
        bool m_power_of_2 = (m & (m - 1)) == 0;

        if (n_power_of_2 && m_power_of_2) {
            return n > m;  // Reverse order for powers of 2
        }

        return n > m;  // Simplified
    }

    /**
     * @brief Symbolic dynamics: partition phase space and encode orbits
     */
    class SymbolicDynamics {
    public:
        /**
         * @brief Encode orbit as symbol sequence
         */
        static std::vector<int> encode(
            const std::vector<double>& orbit,
            const std::vector<double>& partition_points) {

            std::vector<int> symbols;
            symbols.reserve(orbit.size());

            for (double x : orbit) {
                int symbol = 0;
                for (size_t i = 0; i < partition_points.size(); ++i) {
                    if (x > partition_points[i]) {
                        symbol = i + 1;
                    }
                }
                symbols.push_back(symbol);
            }

            return symbols;
        }

        /**
         * @brief Check if sequence is admissible
         */
        static bool isAdmissible(const std::vector<int>& sequence) {
            // For tent map: all sequences admissible
            // For other maps: check transition rules
            return true;  // Placeholder
        }
    };

    /**
     * @brief Box-counting dimension (fractal dimension)
     *
     * D = lim_{ε→0} log(N(ε)) / log(1/ε)
     *
     * where N(ε) = number of boxes of size ε needed to cover set
     */
    static double boxCountingDimension(
        const std::vector<std::vector<double>>& points,
        const std::vector<double>& epsilon_values) {

        std::vector<double> log_N;
        std::vector<double> log_inv_eps;

        for (double eps : epsilon_values) {
            // Count boxes
            std::set<std::vector<int>> boxes;

            for (const auto& p : points) {
                std::vector<int> box_index;
                for (double coord : p) {
                    box_index.push_back(static_cast<int>(std::floor(coord / eps)));
                }
                boxes.insert(box_index);
            }

            log_N.push_back(std::log(boxes.size()));
            log_inv_eps.push_back(std::log(1.0 / eps));
        }

        // Linear regression to find slope
        double mean_x = std::accumulate(log_inv_eps.begin(), log_inv_eps.end(), 0.0) / log_inv_eps.size();
        double mean_y = std::accumulate(log_N.begin(), log_N.end(), 0.0) / log_N.size();

        double numerator = 0.0, denominator = 0.0;
        for (size_t i = 0; i < log_inv_eps.size(); ++i) {
            numerator += (log_inv_eps[i] - mean_x) * (log_N[i] - mean_y);
            denominator += (log_inv_eps[i] - mean_x) * (log_inv_eps[i] - mean_x);
        }

        return numerator / denominator;  // Slope = dimension
    }

    /**
     * @brief Smale horseshoe map (Chapter 13)
     */
    static std::vector<double> smaleHorseshoe(const std::vector<double>& x) {
        // Stretch and fold
        double x_out, y_out;

        if (x[1] < 0.5) {
            // Vertical strip mapped to horizontal
            x_out = 0.25 * x[0];
            y_out = 2.0 * x[1];
        } else {
            x_out = 0.75 + 0.25 * x[0];
            y_out = 2.0 * (x[1] - 0.5);
        }

        return {x_out, y_out};
    }
};

} // namespace maths::dynamical_systems

#endif // MATHS_DYNAMICAL_SYSTEMS_ODE_HPP
