/**
 * @file pde_classification_solutions.hpp
 * @brief Solutions for Parabolic, Elliptic, and Hyperbolic PDEs
 *
 * PARABOLIC EQUATIONS
 * - Heat equation and diffusion processes
 * - Maximum principles
 * - Fundamental solutions
 * - Separation of variables solutions
 * - Transform methods
 * - Numerical stability and convergence
 *
 * ELLIPTIC EQUATIONS
 * - Laplace and Poisson equations
 * - Harmonic functions
 * - Mean value property
 * - Maximum principles
 * - Green's functions
 * - Boundary value problems
 *
 * HYPERBOLIC EQUATIONS
 * - Wave equation
 * - d'Alembert's solution
 * - Energy methods
 * - Domain of dependence
 * - Characteristics and causality
 * - Finite speed of propagation
 */

#ifndef MATHS_PDE_CLASSIFICATION_SOLUTIONS_HPP
#define MATHS_PDE_CLASSIFICATION_SOLUTIONS_HPP

#include <vector>
#include <cmath>
#include <functional>
#include <complex>
#include <algorithm>
#include <numeric>

namespace maths::pde {

/**
 * ============================================================================
 * PARABOLIC EQUATIONS
 * ============================================================================
 */

/**
 * @class ParabolicEquations
 * @brief Solutions and properties of parabolic PDEs
 *
 * Canonical form: u_t = α Δu (heat equation)
 * Key property: Infinite speed of propagation, smoothing effect
 */
class ParabolicEquations {
public:
    /**
     * @brief One-dimensional heat equation: u_t = α u_xx
     */
    struct HeatEquation1D {
        double alpha;  // Thermal diffusivity

        /**
         * @brief Fundamental solution (heat kernel) on infinite domain
         * G(x, t; ξ) = 1/√(4παt) exp(-(x-ξ)²/4αt)
         */
        static double fundamentalSolution(double x, double t, double xi, double alpha) {
            if (t <= 0) return 0.0;
            double denom = std::sqrt(4.0 * M_PI * alpha * t);
            double exponent = -(x - xi) * (x - xi) / (4.0 * alpha * t);
            return std::exp(exponent) / denom;
        }

        /**
         * @brief Solution via convolution with initial data
         * u(x, t) = ∫ G(x, t; ξ) u₀(ξ) dξ
         */
        static double solutionInfiniteDomain(
            std::function<double(double)> u0,
            double x,
            double t,
            double alpha,
            double xi_max = 50.0) {

            int n_points = 1000;
            double dxi = 2.0 * xi_max / n_points;
            double result = 0.0;

            for (int i = -n_points/2; i <= n_points/2; ++i) {
                double xi = i * dxi;
                double weight = (i == -n_points/2 || i == n_points/2) ? 0.5 : 1.0;
                result += weight * u0(xi) * fundamentalSolution(x, t, xi, alpha) * dxi;
            }

            return result;
        }

        /**
         * @brief Separation of variables on [0, L] with Dirichlet BC
         * u(0,t) = u(L,t) = 0
         * Solution: u(x,t) = ∑ Aₙ exp(-α(nπ/L)²t) sin(nπx/L)
         */
        struct DirichletSolution {
            double alpha, L;
            std::vector<double> A_n;  // Fourier sine coefficients

            static DirichletSolution solve(
                std::function<double(double)> u0,
                double alpha,
                double L,
                int n_terms = 50) {

                DirichletSolution sol;
                sol.alpha = alpha;
                sol.L = L;
                sol.A_n.resize(n_terms);

                // Compute Fourier sine coefficients
                int n_points = 1000;
                double dx = L / n_points;

                for (int n = 1; n <= n_terms; ++n) {
                    double sum = 0.0;
                    for (int i = 0; i <= n_points; ++i) {
                        double x = i * dx;
                        double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
                        sum += weight * u0(x) * std::sin(n * M_PI * x / L) * dx;
                    }
                    sol.A_n[n-1] = (2.0 / L) * sum;
                }

                return sol;
            }

            double evaluate(double x, double t) const {
                double result = 0.0;
                for (size_t n = 1; n <= A_n.size(); ++n) {
                    double lambda_n = n * M_PI / L;
                    result += A_n[n-1] *
                             std::exp(-alpha * lambda_n * lambda_n * t) *
                             std::sin(lambda_n * x);
                }
                return result;
            }
        };

        /**
         * @brief Neumann boundary conditions: u_x(0,t) = u_x(L,t) = 0
         * Solution: u(x,t) = A₀ + ∑ Aₙ exp(-α(nπ/L)²t) cos(nπx/L)
         */
        struct NeumannSolution {
            double alpha, L;
            std::vector<double> A_n;  // Fourier cosine coefficients

            static NeumannSolution solve(
                std::function<double(double)> u0,
                double alpha,
                double L,
                int n_terms = 50) {

                NeumannSolution sol;
                sol.alpha = alpha;
                sol.L = L;
                sol.A_n.resize(n_terms + 1);

                // Compute Fourier cosine coefficients
                int n_points = 1000;
                double dx = L / n_points;

                for (int n = 0; n <= n_terms; ++n) {
                    double sum = 0.0;
                    for (int i = 0; i <= n_points; ++i) {
                        double x = i * dx;
                        double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
                        sum += weight * u0(x) * std::cos(n * M_PI * x / L) * dx;
                    }
                    sol.A_n[n] = (n == 0) ? sum / L : (2.0 / L) * sum;
                }

                return sol;
            }

            double evaluate(double x, double t) const {
                double result = A_n[0];  // Constant term (steady state)
                for (size_t n = 1; n < A_n.size(); ++n) {
                    double lambda_n = n * M_PI / L;
                    result += A_n[n] *
                             std::exp(-alpha * lambda_n * lambda_n * t) *
                             std::cos(lambda_n * x);
                }
                return result;
            }
        };
    };

    /**
     * @brief Maximum principles for parabolic equations
     */
    struct MaximumPrinciples {
        /**
         * @brief Weak maximum principle:
         * max u(x,t) occurs either at t=0 or on the spatial boundary
         */
        static bool weakMaximum(
            std::function<double(double, double)> u,
            double x_min, double x_max,
            double t_max,
            int n_points = 100) {

            double interior_max = -1e100;
            double boundary_max = -1e100;

            // Check interior
            for (int i = 1; i < n_points; ++i) {
                for (int j = 1; j < n_points; ++j) {
                    double x = x_min + (x_max - x_min) * i / n_points;
                    double t = t_max * j / n_points;
                    interior_max = std::max(interior_max, u(x, t));
                }
            }

            // Check boundary (t=0 and spatial boundaries)
            for (int i = 0; i <= n_points; ++i) {
                double x = x_min + (x_max - x_min) * i / n_points;
                boundary_max = std::max(boundary_max, u(x, 0.0));
            }

            for (int j = 0; j <= n_points; ++j) {
                double t = t_max * j / n_points;
                boundary_max = std::max(boundary_max, u(x_min, t));
                boundary_max = std::max(boundary_max, u(x_max, t));
            }

            return interior_max <= boundary_max + 1e-6;
        }

        /**
         * @brief Strong maximum principle:
         * If u attains maximum in interior, u is constant
         */
        static bool isConstant(
            std::function<double(double, double)> u,
            double x_min, double x_max,
            double t_max,
            int n_points = 100) {

            double first_value = u((x_min + x_max) / 2, t_max / 2);
            double tolerance = 1e-6;

            for (int i = 0; i <= n_points; ++i) {
                for (int j = 0; j <= n_points; ++j) {
                    double x = x_min + (x_max - x_min) * i / n_points;
                    double t = t_max * j / n_points;
                    if (std::abs(u(x, t) - first_value) > tolerance) {
                        return false;
                    }
                }
            }

            return true;
        }
    };

    /**
     * @brief Two-dimensional heat equation: u_t = α(u_xx + u_yy)
     */
    struct HeatEquation2D {
        /**
         * @brief Fundamental solution in 2D
         * G(x, y, t; ξ, η) = 1/(4παt) exp(-((x-ξ)² + (y-η)²)/4αt)
         */
        static double fundamentalSolution2D(
            double x, double y, double t,
            double xi, double eta,
            double alpha) {

            if (t <= 0) return 0.0;
            double r_squared = (x - xi) * (x - xi) + (y - eta) * (y - eta);
            return std::exp(-r_squared / (4.0 * alpha * t)) / (4.0 * M_PI * alpha * t);
        }

        /**
         * @brief Separation on rectangle [0,Lx] × [0,Ly]
         * u(x,y,t) = ∑∑ Aₘₙ exp(-α(λₘ² + μₙ²)t) sin(λₘx) sin(μₙy)
         */
        struct RectangleSolution {
            double alpha, Lx, Ly;
            std::vector<std::vector<double>> A_mn;

            double evaluate(double x, double y, double t) const {
                double result = 0.0;
                for (size_t m = 1; m <= A_mn.size(); ++m) {
                    for (size_t n = 1; n <= A_mn[m-1].size(); ++n) {
                        double lambda_m = m * M_PI / Lx;
                        double mu_n = n * M_PI / Ly;
                        double eigenvalue = lambda_m * lambda_m + mu_n * mu_n;
                        result += A_mn[m-1][n-1] *
                                 std::exp(-alpha * eigenvalue * t) *
                                 std::sin(lambda_m * x) * std::sin(mu_n * y);
                    }
                }
                return result;
            }
        };
    };
};

/**
 * ============================================================================
 * ELLIPTIC EQUATIONS
 * ============================================================================
 */

/**
 * @class EllipticEquations
 * @brief Solutions and properties of elliptic PDEs
 *
 * Canonical form: Δu = f (Poisson), Δu = 0 (Laplace)
 * Key property: No time evolution, boundary value problems
 */
class EllipticEquations {
public:
    /**
     * @brief Laplace equation: Δu = u_xx + u_yy = 0
     */
    struct LaplaceEquation {
        /**
         * @brief Mean value property:
         * u(x₀, y₀) = (1/2πr) ∫ u(x₀ + r cos θ, y₀ + r sin θ) dθ
         */
        static double meanValue(
            std::function<double(double, double)> u,
            double x0, double y0,
            double r,
            int n_theta = 100) {

            double dtheta = 2.0 * M_PI / n_theta;
            double sum = 0.0;

            for (int i = 0; i < n_theta; ++i) {
                double theta = i * dtheta;
                double x = x0 + r * std::cos(theta);
                double y = y0 + r * std::sin(theta);
                sum += u(x, y);
            }

            return sum / n_theta;
        }

        /**
         * @brief Verify mean value property holds
         */
        static bool verifyMeanValue(
            std::function<double(double, double)> u,
            double x0, double y0,
            double r,
            double tolerance = 1e-4) {

            double center_value = u(x0, y0);
            double mean = meanValue(u, x0, y0, r);
            return std::abs(center_value - mean) < tolerance;
        }

        /**
         * @brief Solution on rectangle [0, Lx] × [0, Ly] with Dirichlet BC
         */
        struct RectangleDirichlet {
            double Lx, Ly;

            /**
             * @brief Solution with boundary data u(x, 0) = f(x), others = 0
             */
            static std::function<double(double, double)> solveSingleEdge(
                std::function<double(double)> f,
                double Lx, double Ly,
                int n_terms = 50) {

                // Compute Fourier coefficients
                std::vector<double> b_n(n_terms);
                int n_points = 1000;
                double dx = Lx / n_points;

                for (int n = 1; n <= n_terms; ++n) {
                    double sum = 0.0;
                    for (int i = 0; i <= n_points; ++i) {
                        double x = i * dx;
                        double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
                        sum += weight * f(x) * std::sin(n * M_PI * x / Lx) * dx;
                    }
                    b_n[n-1] = (2.0 / Lx) * sum;
                }

                return [b_n, Lx, Ly](double x, double y) {
                    double result = 0.0;
                    for (size_t n = 1; n <= b_n.size(); ++n) {
                        double lambda_n = n * M_PI / Lx;
                        result += b_n[n-1] *
                                 std::sinh(lambda_n * (Ly - y)) / std::sinh(lambda_n * Ly) *
                                 std::sin(lambda_n * x);
                    }
                    return result;
                };
            }
        };

        /**
         * @brief Solution in circular domain r < R
         * Using polar coordinates: u_rr + (1/r)u_r + (1/r²)u_θθ = 0
         */
        struct CircularDomain {
            double R;  // Radius

            /**
             * @brief Poisson integral formula
             * u(r, θ) = (1/2π) ∫₀^(2π) [(R²-r²)/(R²-2Rr cos(θ-φ)+r²)] u(R,φ) dφ
             */
            static double poissonIntegral(
                std::function<double(double)> boundary_data,
                double r, double theta,
                double R,
                int n_phi = 100) {

                double dphi = 2.0 * M_PI / n_phi;
                double sum = 0.0;

                for (int i = 0; i < n_phi; ++i) {
                    double phi = i * dphi;
                    double cos_diff = std::cos(theta - phi);
                    double kernel = (R*R - r*r) / (R*R - 2.0*R*r*cos_diff + r*r);
                    sum += kernel * boundary_data(phi);
                }

                return sum / (2.0 * M_PI) * dphi;
            }
        };
    };

    /**
     * @brief Poisson equation: Δu = f
     */
    struct PoissonEquation {
        /**
         * @brief Green's function for 2D Laplacian on unbounded domain
         * G(x, y; ξ, η) = -(1/2π) ln(r) where r = √((x-ξ)² + (y-η)²)
         */
        static double greensFunctionUnbounded(double x, double y, double xi, double eta) {
            double r_squared = (x - xi) * (x - xi) + (y - eta) * (y - eta);
            if (r_squared < 1e-10) return 0.0;
            return -std::log(std::sqrt(r_squared)) / (2.0 * M_PI);
        }

        /**
         * @brief Solution via Green's function: u(x,y) = ∫∫ G(x,y;ξ,η) f(ξ,η) dξdη
         */
        static double solvePoissonUnbounded(
            std::function<double(double, double)> f,
            double x, double y,
            double domain_size = 10.0,
            int n_points = 50) {

            double dxi = 2.0 * domain_size / n_points;
            double deta = 2.0 * domain_size / n_points;
            double result = 0.0;

            for (int i = -n_points/2; i <= n_points/2; ++i) {
                for (int j = -n_points/2; j <= n_points/2; ++j) {
                    double xi = i * dxi;
                    double eta = j * deta;
                    result += greensFunctionUnbounded(x, y, xi, eta) * f(xi, eta) * dxi * deta;
                }
            }

            return result;
        }
    };

    /**
     * @brief Maximum principles for elliptic equations
     */
    struct MaximumPrinciples {
        /**
         * @brief Maximum principle for Laplace equation:
         * max and min occur on boundary
         */
        static bool boundaryExtremum(
            std::function<double(double, double)> u,
            double x_min, double x_max,
            double y_min, double y_max,
            int n_points = 50) {

            double interior_max = -1e100;
            double interior_min = 1e100;
            double boundary_max = -1e100;
            double boundary_min = 1e100;

            // Check interior
            for (int i = 1; i < n_points; ++i) {
                for (int j = 1; j < n_points; ++j) {
                    double x = x_min + (x_max - x_min) * i / n_points;
                    double y = y_min + (y_max - y_min) * j / n_points;
                    double val = u(x, y);
                    interior_max = std::max(interior_max, val);
                    interior_min = std::min(interior_min, val);
                }
            }

            // Check boundary
            for (int i = 0; i <= n_points; ++i) {
                double x = x_min + (x_max - x_min) * i / n_points;
                double y = y_min + (y_max - y_min) * i / n_points;

                double val1 = u(x, y_min);
                double val2 = u(x, y_max);
                double val3 = u(x_min, y);
                double val4 = u(x_max, y);

                boundary_max = std::max({boundary_max, val1, val2, val3, val4});
                boundary_min = std::min({boundary_min, val1, val2, val3, val4});
            }

            return (interior_max <= boundary_max + 1e-6) &&
                   (interior_min >= boundary_min - 1e-6);
        }
    };
};

/**
 * ============================================================================
 * HYPERBOLIC EQUATIONS
 * ============================================================================
 */

/**
 * @class HyperbolicEquations
 * @brief Solutions and properties of hyperbolic PDEs
 *
 * Canonical form: u_tt = c² Δu (wave equation)
 * Key property: Finite speed of propagation, energy conservation
 */
class HyperbolicEquations {
public:
    /**
     * @brief One-dimensional wave equation: u_tt = c² u_xx
     */
    struct WaveEquation1D {
        double c;  // Wave speed

        /**
         * @brief d'Alembert's solution on infinite domain
         * u(x, t) = ½[f(x+ct) + f(x-ct)] + 1/(2c) ∫_{x-ct}^{x+ct} g(s) ds
         */
        static double dAlembertSolution(
            std::function<double(double)> f,  // Initial displacement
            std::function<double(double)> g,  // Initial velocity
            double x, double t,
            double c) {

            double term1 = 0.5 * (f(x + c*t) + f(x - c*t));

            // Integrate g from x-ct to x+ct
            double x_minus = x - c*t;
            double x_plus = x + c*t;
            int n_points = 100;
            double ds = (x_plus - x_minus) / n_points;
            double integral = 0.0;

            for (int i = 0; i <= n_points; ++i) {
                double s = x_minus + i * ds;
                double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
                integral += weight * g(s) * ds;
            }

            return term1 + integral / (2.0 * c);
        }

        /**
         * @brief Domain of dependence
         * Solution at (x, t) depends only on [x-ct, x+ct]
         */
        struct DomainOfDependence {
            double x, t, c;

            std::pair<double, double> getInterval() const {
                return {x - c*t, x + c*t};
            }

            bool isInDomain(double xi) const {
                return (xi >= x - c*t) && (xi <= x + c*t);
            }
        };

        /**
         * @brief Separation of variables on [0, L] with fixed endpoints
         * u(x,t) = ∑ (Aₙ cos(ωₙt) + Bₙ sin(ωₙt)) sin(nπx/L)
         * where ωₙ = nπc/L
         */
        struct StandingWaveSolution {
            double c, L;
            std::vector<double> A_n, B_n;

            static StandingWaveSolution solve(
                std::function<double(double)> f,  // Initial displacement
                std::function<double(double)> g,  // Initial velocity
                double c, double L,
                int n_terms = 50) {

                StandingWaveSolution sol;
                sol.c = c;
                sol.L = L;
                sol.A_n.resize(n_terms);
                sol.B_n.resize(n_terms);

                int n_points = 1000;
                double dx = L / n_points;

                for (int n = 1; n <= n_terms; ++n) {
                    double omega_n = n * M_PI * c / L;

                    // Aₙ from f(x)
                    double sum_f = 0.0;
                    for (int i = 0; i <= n_points; ++i) {
                        double x = i * dx;
                        double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
                        sum_f += weight * f(x) * std::sin(n * M_PI * x / L) * dx;
                    }
                    sol.A_n[n-1] = (2.0 / L) * sum_f;

                    // Bₙ from g(x)
                    double sum_g = 0.0;
                    for (int i = 0; i <= n_points; ++i) {
                        double x = i * dx;
                        double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
                        sum_g += weight * g(x) * std::sin(n * M_PI * x / L) * dx;
                    }
                    sol.B_n[n-1] = (2.0 / L) * sum_g / omega_n;
                }

                return sol;
            }

            double evaluate(double x, double t) const {
                double result = 0.0;
                for (size_t n = 1; n <= A_n.size(); ++n) {
                    double omega_n = n * M_PI * c / L;
                    result += (A_n[n-1] * std::cos(omega_n * t) +
                              B_n[n-1] * std::sin(omega_n * t)) *
                              std::sin(n * M_PI * x / L);
                }
                return result;
            }
        };
    };

    /**
     * @brief Energy conservation for wave equation
     */
    struct EnergyConservation {
        /**
         * @brief Total energy: E = ½∫[u_t² + c²u_x²]dx
         * For wave equation, energy is conserved
         */
        static double computeEnergy(
            std::function<double(double, double)> u,
            std::function<double(double, double)> u_t,
            std::function<double(double, double)> u_x,
            double c,
            double x_min, double x_max,
            double t,
            int n_points = 1000) {

            double dx = (x_max - x_min) / n_points;
            double energy = 0.0;

            for (int i = 0; i <= n_points; ++i) {
                double x = x_min + i * dx;
                double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;

                double kinetic = u_t(x, t) * u_t(x, t);
                double potential = c*c * u_x(x, t) * u_x(x, t);

                energy += weight * (kinetic + potential) * dx;
            }

            return 0.5 * energy;
        }
    };

    /**
     * @brief Two-dimensional wave equation: u_tt = c²(u_xx + u_yy)
     */
    struct WaveEquation2D {
        double c;  // Wave speed

        /**
         * @brief Separation on rectangle [0,Lx] × [0,Ly]
         * u(x,y,t) = ∑∑ (Aₘₙ cos(ωₘₙt) + Bₘₙ sin(ωₘₙt)) sin(λₘx) sin(μₙy)
         * where ωₘₙ = c√(λₘ² + μₙ²)
         */
        struct RectangleSolution {
            double c, Lx, Ly;
            std::vector<std::vector<double>> A_mn, B_mn;

            double evaluate(double x, double y, double t) const {
                double result = 0.0;
                for (size_t m = 1; m <= A_mn.size(); ++m) {
                    for (size_t n = 1; n <= A_mn[m-1].size(); ++n) {
                        double lambda_m = m * M_PI / Lx;
                        double mu_n = n * M_PI / Ly;
                        double omega_mn = c * std::sqrt(lambda_m*lambda_m + mu_n*mu_n);

                        result += (A_mn[m-1][n-1] * std::cos(omega_mn * t) +
                                  B_mn[m-1][n-1] * std::sin(omega_mn * t)) *
                                  std::sin(lambda_m * x) * std::sin(mu_n * y);
                    }
                }
                return result;
            }
        };

        /**
         * @brief Characteristics and causality
         * Disturbances propagate along characteristic cones
         */
        struct CharacteristicCones {
            double x0, y0, t0;  // Source point
            double c;           // Wave speed

            /**
             * @brief Check if point (x, y, t) is in future light cone
             */
            bool isInFutureCone(double x, double y, double t) const {
                if (t < t0) return false;
                double r = std::sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0));
                return r <= c * (t - t0);
            }

            /**
             * @brief Check if point is in past light cone
             */
            bool isInPastCone(double x, double y, double t) const {
                if (t > t0) return false;
                double r = std::sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0));
                return r <= c * (t0 - t);
            }
        };
    };

    /**
     * @brief Finite speed of propagation
     * Disturbances travel at speed c, not instantaneously
     */
    struct FiniteSpeedPropagation {
        /**
         * @brief Maximum distance traveled by wave in time t
         */
        static double maxDistance(double c, double t) {
            return c * t;
        }

        /**
         * @brief Check if point x is affected by source at x0 at time t
         */
        static bool isAffected(double x, double x0, double t, double c) {
            return std::abs(x - x0) <= c * t;
        }
    };
};

} // namespace maths::pde

#endif // MATHS_PDE_CLASSIFICATION_SOLUTIONS_HPP
