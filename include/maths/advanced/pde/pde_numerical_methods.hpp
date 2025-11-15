/**
 * @file pde_numerical_methods.hpp
 * @brief Numerical Methods and Approximation Techniques for PDEs
 *
 * TAYLOR SERIES EXPANSIONS
 * - Finite difference approximations
 * - Truncation error analysis
 * - Higher order schemes
 *
 * SUCCESSIVE APPROXIMATIONS
 * - Picard iteration
 * - Fixed point methods
 * - Convergence criteria
 *
 * BOUNDARY PERTURBATIONS
 * - Perturbation methods
 * - Regular and singular perturbations
 * - Boundary layer analysis
 *
 * FINITE DIFFERENCE SCHEMES
 * First Order Equations:
 * - Upwind schemes
 * - Lax-Friedrichs
 * - Lax-Wendroff
 * - CFL condition
 *
 * Second Order Equations:
 * - Central differences
 * - Explicit and implicit schemes
 * - ADI (Alternating Direction Implicit)
 * - Stability analysis
 */

#ifndef MATHS_PDE_NUMERICAL_METHODS_HPP
#define MATHS_PDE_NUMERICAL_METHODS_HPP

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdexcept>

namespace maths::pde {

/**
 * ============================================================================
 * TAYLOR SERIES EXPANSIONS
 * ============================================================================
 */

/**
 * @class TaylorExpansions
 * @brief Finite difference approximations via Taylor series
 */
class TaylorExpansions {
public:
    /**
     * @brief Forward difference: f'(x) ≈ (f(x+h) - f(x))/h
     * Truncation error: O(h)
     */
    static double forwardDifference(
        std::function<double(double)> f,
        double x, double h) {
        return (f(x + h) - f(x)) / h;
    }

    /**
     * @brief Backward difference: f'(x) ≈ (f(x) - f(x-h))/h
     * Truncation error: O(h)
     */
    static double backwardDifference(
        std::function<double(double)> f,
        double x, double h) {
        return (f(x) - f(x - h)) / h;
    }

    /**
     * @brief Central difference: f'(x) ≈ (f(x+h) - f(x-h))/(2h)
     * Truncation error: O(h²)
     */
    static double centralDifference(
        std::function<double(double)> f,
        double x, double h) {
        return (f(x + h) - f(x - h)) / (2.0 * h);
    }

    /**
     * @brief Second derivative: f''(x) ≈ (f(x+h) - 2f(x) + f(x-h))/h²
     * Truncation error: O(h²)
     */
    static double secondDerivative(
        std::function<double(double)> f,
        double x, double h) {
        return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
    }

    /**
     * @brief Higher order central difference (4th order accurate)
     * f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h))/(12h)
     */
    static double centralDifference4thOrder(
        std::function<double(double)> f,
        double x, double h) {
        return (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h)) / (12.0 * h);
    }

    /**
     * @brief Truncation error estimate
     */
    static double truncationError(
        std::function<double(double)> f_exact,
        std::function<double(double)> f_approx,
        double x) {
        return std::abs(f_exact(x) - f_approx(x));
    }
};

/**
 * ============================================================================
 * SUCCESSIVE APPROXIMATIONS
 * ============================================================================
 */

/**
 * @class SuccessiveApproximations
 * @brief Iterative methods for solving PDEs
 */
class SuccessiveApproximations {
public:
    /**
     * @brief Picard iteration for IVP: u' = f(t, u), u(t0) = u0
     * u_{n+1}(t) = u0 + ∫_{t0}^t f(s, u_n(s)) ds
     */
    static std::function<double(double)> picardIteration(
        std::function<double(double, double)> f,
        double t0, double u0,
        int n_iterations,
        double t_max) {

        // Start with u_0(t) = u0 (constant)
        auto u_n = [u0](double t) { return u0; };

        for (int iter = 0; iter < n_iterations; ++iter) {
            auto u_next = [f, t0, u0, u_n](double t) {
                // Numerical integration from t0 to t
                int n_points = 100;
                double dt = (t - t0) / n_points;
                double integral = 0.0;

                for (int i = 0; i <= n_points; ++i) {
                    double s = t0 + i * dt;
                    double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
                    integral += weight * f(s, u_n(s)) * dt;
                }

                return u0 + integral;
            };

            u_n = u_next;
        }

        return u_n;
    }

    /**
     * @brief Fixed point iteration for F(u) = 0 written as u = G(u)
     * u_{n+1} = G(u_n)
     */
    static std::vector<double> fixedPointIteration(
        std::function<std::vector<double>(const std::vector<double>&)> G,
        std::vector<double> u0,
        int max_iter = 100,
        double tolerance = 1e-6) {

        std::vector<double> u = u0;

        for (int iter = 0; iter < max_iter; ++iter) {
            std::vector<double> u_new = G(u);

            // Check convergence
            double error = 0.0;
            for (size_t i = 0; i < u.size(); ++i) {
                error += (u_new[i] - u[i]) * (u_new[i] - u[i]);
            }
            error = std::sqrt(error);

            if (error < tolerance) {
                return u_new;
            }

            u = u_new;
        }

        throw std::runtime_error("Fixed point iteration did not converge");
    }

    /**
     * @brief Successive Over-Relaxation (SOR) for linear systems
     * x_{i}^{(k+1)} = (1-ω)x_{i}^{(k)} + (ω/a_{ii})(b_i - ∑_{j<i} a_{ij}x_{j}^{(k+1)} - ∑_{j>i} a_{ij}x_{j}^{(k)})
     */
    static std::vector<double> SOR(
        const std::vector<std::vector<double>>& A,
        const std::vector<double>& b,
        double omega = 1.5,
        int max_iter = 1000,
        double tolerance = 1e-6) {

        int n = b.size();
        std::vector<double> x(n, 0.0);

        for (int iter = 0; iter < max_iter; ++iter) {
            double max_change = 0.0;

            for (int i = 0; i < n; ++i) {
                double sigma = 0.0;
                for (int j = 0; j < n; ++j) {
                    if (j != i) {
                        sigma += A[i][j] * x[j];
                    }
                }

                double x_new = (1.0 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sigma);
                max_change = std::max(max_change, std::abs(x_new - x[i]));
                x[i] = x_new;
            }

            if (max_change < tolerance) {
                return x;
            }
        }

        throw std::runtime_error("SOR did not converge");
    }
};

/**
 * ============================================================================
 * BOUNDARY PERTURBATIONS
 * ============================================================================
 */

/**
 * @class BoundaryPerturbations
 * @brief Perturbation methods for PDEs
 */
class BoundaryPerturbations {
public:
    /**
     * @brief Regular perturbation: u = u₀ + εu₁ + ε²u₂ + ...
     * For εu'' + u = f(x)
     */
    struct RegularPerturbation {
        std::vector<std::function<double(double)>> terms;
        double epsilon;

        /**
         * @brief Evaluate perturbation series
         */
        double evaluate(double x, int n_terms) const {
            double result = 0.0;
            double eps_power = 1.0;

            for (int i = 0; i < std::min(n_terms, (int)terms.size()); ++i) {
                result += eps_power * terms[i](x);
                eps_power *= epsilon;
            }

            return result;
        }
    };

    /**
     * @brief Singular perturbation: εu'' + a(x)u' + b(x)u = f(x)
     * Boundary layer at x = 0 or x = 1
     */
    struct SingularPerturbation {
        double epsilon;

        /**
         * @brief Outer solution (away from boundary layer)
         */
        std::function<double(double)> outerSolution;

        /**
         * @brief Inner solution (boundary layer)
         * Stretched coordinate: ξ = x/√ε
         */
        std::function<double(double)> innerSolution;

        /**
         * @brief Matched asymptotic expansion
         */
        double evaluate(double x) const {
            // Simplified composite solution
            return outerSolution(x) + innerSolution(x / std::sqrt(epsilon));
        }
    };

    /**
     * @brief Boundary layer thickness estimate
     */
    static double boundaryLayerThickness(double epsilon) {
        return std::sqrt(epsilon);  // For second order problems
    }
};

/**
 * ============================================================================
 * FINITE DIFFERENCE SCHEMES FOR FIRST ORDER EQUATIONS
 * ============================================================================
 */

/**
 * @class FirstOrderSchemes
 * @brief Finite difference methods for first order PDEs
 *
 * Linear advection: u_t + c u_x = 0
 */
class FirstOrderSchemes {
public:
    /**
     * @brief Upwind scheme (first order in space and time)
     * u_i^{n+1} = u_i^n - c(Δt/Δx)(u_i^n - u_{i-1}^n) for c > 0
     */
    static std::vector<double> upwind(
        const std::vector<double>& u,
        double c, double dt, double dx) {

        int n = u.size();
        std::vector<double> u_new(n);

        double cfl = c * dt / dx;

        if (c > 0) {
            // Right-moving wave: backward difference in space
            u_new[0] = u[0];  // Boundary condition
            for (int i = 1; i < n; ++i) {
                u_new[i] = u[i] - cfl * (u[i] - u[i-1]);
            }
        } else {
            // Left-moving wave: forward difference in space
            for (int i = 0; i < n - 1; ++i) {
                u_new[i] = u[i] - cfl * (u[i+1] - u[i]);
            }
            u_new[n-1] = u[n-1];  // Boundary condition
        }

        return u_new;
    }

    /**
     * @brief Lax-Friedrichs scheme
     * u_i^{n+1} = ½(u_{i+1}^n + u_{i-1}^n) - (c Δt)/(2Δx)(u_{i+1}^n - u_{i-1}^n)
     */
    static std::vector<double> laxFriedrichs(
        const std::vector<double>& u,
        double c, double dt, double dx) {

        int n = u.size();
        std::vector<double> u_new(n);

        double cfl = c * dt / dx;

        u_new[0] = u[0];  // Boundary
        u_new[n-1] = u[n-1];  // Boundary

        for (int i = 1; i < n - 1; ++i) {
            u_new[i] = 0.5 * (u[i+1] + u[i-1]) - 0.5 * cfl * (u[i+1] - u[i-1]);
        }

        return u_new;
    }

    /**
     * @brief Lax-Wendroff scheme (second order)
     * u_i^{n+1} = u_i^n - (c Δt)/(2Δx)(u_{i+1}^n - u_{i-1}^n)
     *             + (c² Δt²)/(2Δx²)(u_{i+1}^n - 2u_i^n + u_{i-1}^n)
     */
    static std::vector<double> laxWendroff(
        const std::vector<double>& u,
        double c, double dt, double dx) {

        int n = u.size();
        std::vector<double> u_new(n);

        double cfl = c * dt / dx;
        double cfl2 = cfl * cfl;

        u_new[0] = u[0];  // Boundary
        u_new[n-1] = u[n-1];  // Boundary

        for (int i = 1; i < n - 1; ++i) {
            u_new[i] = u[i]
                     - 0.5 * cfl * (u[i+1] - u[i-1])
                     + 0.5 * cfl2 * (u[i+1] - 2.0*u[i] + u[i-1]);
        }

        return u_new;
    }

    /**
     * @brief CFL (Courant-Friedrichs-Lewy) condition
     * For stability: |c| Δt/Δx ≤ 1
     */
    static bool checkCFL(double c, double dt, double dx) {
        return std::abs(c * dt / dx) <= 1.0;
    }

    /**
     * @brief Maximum stable time step
     */
    static double maxTimeStep(double c, double dx) {
        return dx / std::abs(c);
    }
};

/**
 * ============================================================================
 * FINITE DIFFERENCE SCHEMES FOR SECOND ORDER EQUATIONS
 * ============================================================================
 */

/**
 * @class SecondOrderSchemes
 * @brief Finite difference methods for second order PDEs
 */
class SecondOrderSchemes {
public:
    /**
     * @brief Explicit scheme for heat equation: u_t = α u_xx
     * u_i^{n+1} = u_i^n + r(u_{i+1}^n - 2u_i^n + u_{i-1}^n)
     * where r = α Δt/Δx²
     * Stable if r ≤ 1/2
     */
    static std::vector<double> explicitHeat(
        const std::vector<double>& u,
        double alpha, double dt, double dx) {

        int n = u.size();
        std::vector<double> u_new(n);

        double r = alpha * dt / (dx * dx);

        if (r > 0.5) {
            throw std::runtime_error("Explicit heat scheme unstable: r > 0.5");
        }

        u_new[0] = u[0];  // Boundary
        u_new[n-1] = u[n-1];  // Boundary

        for (int i = 1; i < n - 1; ++i) {
            u_new[i] = u[i] + r * (u[i+1] - 2.0*u[i] + u[i-1]);
        }

        return u_new;
    }

    /**
     * @brief Implicit (backward Euler) scheme for heat equation
     * u_i^{n+1} - r(u_{i+1}^{n+1} - 2u_i^{n+1} + u_{i-1}^{n+1}) = u_i^n
     * Unconditionally stable
     */
    static std::vector<double> implicitHeat(
        const std::vector<double>& u,
        double alpha, double dt, double dx) {

        int n = u.size();
        double r = alpha * dt / (dx * dx);

        // Tridiagonal system: A*u_new = u
        std::vector<double> lower(n-1), diag(n), upper(n-1);

        diag[0] = 1.0;  // Boundary
        diag[n-1] = 1.0;  // Boundary

        for (int i = 1; i < n - 1; ++i) {
            lower[i-1] = -r;
            diag[i] = 1.0 + 2.0*r;
            upper[i] = -r;
        }

        return solveTridiagonal(lower, diag, upper, u);
    }

    /**
     * @brief Crank-Nicolson scheme (θ = 1/2)
     * Unconditionally stable, second order in time
     */
    static std::vector<double> crankNicolson(
        const std::vector<double>& u,
        double alpha, double dt, double dx) {

        int n = u.size();
        double r = alpha * dt / (dx * dx);

        // Left side: (I + 0.5r A)u^{n+1}
        // Right side: (I - 0.5r A)u^n

        std::vector<double> lower(n-1), diag(n), upper(n-1);
        std::vector<double> rhs(n);

        diag[0] = 1.0;  rhs[0] = u[0];
        diag[n-1] = 1.0;  rhs[n-1] = u[n-1];

        for (int i = 1; i < n - 1; ++i) {
            lower[i-1] = -0.5 * r;
            diag[i] = 1.0 + r;
            upper[i] = -0.5 * r;

            rhs[i] = 0.5*r*u[i-1] + (1.0-r)*u[i] + 0.5*r*u[i+1];
        }

        return solveTridiagonal(lower, diag, upper, rhs);
    }

    /**
     * @brief ADI (Alternating Direction Implicit) for 2D heat equation
     * Splits 2D problem into sequence of 1D problems
     */
    struct ADI2D {
        /**
         * @brief One time step of ADI
         * Step 1: Implicit in x, explicit in y
         * Step 2: Explicit in x, implicit in y
         */
        static std::vector<std::vector<double>> step(
            const std::vector<std::vector<double>>& u,
            double alpha, double dt, double dx, double dy) {

            int nx = u.size();
            int ny = u[0].size();
            double rx = alpha * dt / (dx * dx);
            double ry = alpha * dt / (dy * dy);

            // Intermediate step: implicit in x
            std::vector<std::vector<double>> u_star(nx, std::vector<double>(ny));

            for (int j = 0; j < ny; ++j) {
                std::vector<double> u_col(nx);
                for (int i = 0; i < nx; ++i) {
                    u_col[i] = u[i][j];
                    if (j > 0 && j < ny-1) {
                        u_col[i] += 0.5 * ry * (u[i][j+1] - 2.0*u[i][j] + u[i][j-1]);
                    }
                }

                // Solve tridiagonal system in x-direction
                std::vector<double> lower(nx-1, -0.5*rx);
                std::vector<double> diag(nx, 1.0 + rx);
                std::vector<double> upper(nx-1, -0.5*rx);

                diag[0] = 1.0; diag[nx-1] = 1.0;
                auto u_new_col = solveTridiagonal(lower, diag, upper, u_col);

                for (int i = 0; i < nx; ++i) {
                    u_star[i][j] = u_new_col[i];
                }
            }

            // Final step: implicit in y
            std::vector<std::vector<double>> u_new(nx, std::vector<double>(ny));

            for (int i = 0; i < nx; ++i) {
                std::vector<double> u_row(ny);
                for (int j = 0; j < ny; ++j) {
                    u_row[j] = u_star[i][j];
                    if (i > 0 && i < nx-1) {
                        u_row[j] += 0.5 * rx * (u_star[i+1][j] - 2.0*u_star[i][j] + u_star[i-1][j]);
                    }
                }

                std::vector<double> lower(ny-1, -0.5*ry);
                std::vector<double> diag(ny, 1.0 + ry);
                std::vector<double> upper(ny-1, -0.5*ry);

                diag[0] = 1.0; diag[ny-1] = 1.0;
                auto u_new_row = solveTridiagonal(lower, diag, upper, u_row);

                for (int j = 0; j < ny; ++j) {
                    u_new[i][j] = u_new_row[j];
                }
            }

            return u_new;
        }
    };

    /**
     * @brief Stability analysis via von Neumann method
     * Amplification factor for explicit heat equation: g = 1 - 4r sin²(kΔx/2)
     */
    static double amplificationFactor(double r, double k, double dx) {
        return 1.0 - 4.0 * r * std::sin(k * dx / 2.0) * std::sin(k * dx / 2.0);
    }

private:
    /**
     * @brief Thomas algorithm for tridiagonal systems
     */
    static std::vector<double> solveTridiagonal(
        const std::vector<double>& lower,
        const std::vector<double>& diag,
        const std::vector<double>& upper,
        const std::vector<double>& rhs) {

        int n = diag.size();
        std::vector<double> c_prime(n-1);
        std::vector<double> d_prime(n);

        // Forward sweep
        c_prime[0] = upper[0] / diag[0];
        d_prime[0] = rhs[0] / diag[0];

        for (int i = 1; i < n - 1; ++i) {
            double denom = diag[i] - lower[i-1] * c_prime[i-1];
            c_prime[i] = upper[i] / denom;
            d_prime[i] = (rhs[i] - lower[i-1] * d_prime[i-1]) / denom;
        }

        d_prime[n-1] = (rhs[n-1] - lower[n-2] * d_prime[n-2]) /
                       (diag[n-1] - lower[n-2] * c_prime[n-2]);

        // Back substitution
        std::vector<double> x(n);
        x[n-1] = d_prime[n-1];
        for (int i = n - 2; i >= 0; --i) {
            x[i] = d_prime[i] - c_prime[i] * x[i+1];
        }

        return x;
    }
};

} // namespace maths::pde

#endif // MATHS_PDE_NUMERICAL_METHODS_HPP
