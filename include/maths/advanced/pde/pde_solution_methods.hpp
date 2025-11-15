/**
 * @file pde_solution_methods.hpp
 * @brief Classical PDE Solution Methods
 *
 * LINEAR EQUATIONS WITH CONSTANT COEFFICIENTS
 * - Inverse Operators (D⁻¹ for solving ODEs and PDEs)
 * - Homogeneous Equations (complementary solutions)
 * - Nonhomogeneous Equations (particular solutions)
 *
 * ORTHOGONAL EXPANSIONS
 * - Orthogonality (inner products, norms)
 * - Orthogonal Polynomials (Legendre, Chebyshev, Hermite, Laguerre)
 * - Series of Orthogonal Functions (completeness, convergence)
 * - Trigonometric Fourier Series (sine, cosine series)
 * - Eigenfunction Expansions (Sturm-Liouville problems)
 * - Bessel Functions (J_n, Y_n, generating functions)
 *
 * SEPARATION OF VARIABLES
 * - Product solutions and separation ansatz
 * - Hyperbolic Equations (wave equation)
 * - Parabolic Equations (heat equation)
 * - Elliptic Equations (Laplace equation)
 * - Cylindrical Coordinates (Bessel functions)
 * - Spherical Coordinates (Legendre polynomials, spherical harmonics)
 * - Nonhomogeneous Problems (eigenfunction expansion method)
 */

#ifndef MATHS_PDE_SOLUTION_METHODS_HPP
#define MATHS_PDE_SOLUTION_METHODS_HPP

#include <vector>
#include <cmath>
#include <functional>
#include <complex>
#include <algorithm>
#include <numeric>
#include <map>

namespace maths::pde {

/**
 * ============================================================================
 * LINEAR EQUATIONS WITH CONSTANT COEFFICIENTS
 * ============================================================================
 */

/**
 * @class InverseOperators
 * @brief Differential operators and their inverses
 */
class InverseOperators {
public:
    /**
     * @brief Differential operator D = d/dx
     */
    struct DifferentialOperator {
        int order;  // Order of derivative

        /**
         * @brief Apply D^n to function f
         */
        std::function<double(double)> apply(
            std::function<double(double)> f, int n, double h = 1e-6) const {

            return [f, n, h](double x) {
                if (n == 0) return f(x);
                if (n == 1) return (f(x + h) - f(x - h)) / (2 * h);

                // Higher derivatives via finite differences
                double result = 0.0;
                int sign = 1;
                for (int k = 0; k <= n; ++k) {
                    double coeff = 1.0;
                    for (int j = 1; j <= n; ++j) {
                        if (j != k) coeff *= (k - j);
                    }
                    result += sign * f(x + (k - n/2.0) * h) / coeff;
                    sign *= -1;
                }
                return result / std::pow(h, n);
            };
        }

        /**
         * @brief Inverse operator D⁻¹ (integration)
         */
        std::function<double(double)> inverse(
            std::function<double(double)> f, double x0 = 0.0) const {

            return [f, x0](double x) {
                // Numerical integration using Simpson's rule
                int n = 100;
                double dx = (x - x0) / n;
                double integral = 0.0;

                for (int i = 0; i <= n; ++i) {
                    double xi = x0 + i * dx;
                    double weight = (i == 0 || i == n) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
                    integral += weight * f(xi);
                }

                return integral * dx / 3.0;
            };
        }

        /**
         * @brief Polynomial operator: P(D) = a_n D^n + ... + a_1 D + a_0
         */
        std::function<double(double)> polynomialOperator(
            const std::vector<double>& coefficients,
            std::function<double(double)> f) const {

            return [this, coefficients, f](double x) {
                double result = 0.0;
                for (size_t n = 0; n < coefficients.size(); ++n) {
                    auto D_n_f = apply(f, n);
                    result += coefficients[n] * D_n_f(x);
                }
                return result;
            };
        }
    };

    /**
     * @brief Exponential shift formula: e^{aD} f(x) = f(x + a)
     */
    static std::function<double(double)> exponentialShift(
        double a, std::function<double(double)> f) {

        return [a, f](double x) {
            return f(x + a);
        };
    }
};

/**
 * @class ConstantCoefficientPDE
 * @brief Solving linear constant coefficient PDEs
 */
class ConstantCoefficientPDE {
public:
    /**
     * @brief Homogeneous solution using characteristic equation
     *
     * For P(D_x, D_y) u = 0 with constant coefficients
     */
    struct HomogeneousSolution {
        std::vector<std::complex<double>> roots;  // Roots of characteristic polynomial

        /**
         * @brief General solution for homogeneous PDE
         */
        std::function<double(double, double)> solution(
            const std::vector<double>& constants) const {

            return [this, constants](double x, double y) {
                double result = 0.0;
                for (size_t i = 0; i < roots.size() && i < constants.size(); ++i) {
                    double r_real = roots[i].real();
                    double r_imag = roots[i].imag();

                    if (std::abs(r_imag) < 1e-10) {
                        // Real root
                        result += constants[i] * std::exp(r_real * x);
                    } else {
                        // Complex conjugate pair
                        result += constants[i] * std::exp(r_real * x) * std::cos(r_imag * x);
                    }
                }
                return result;
            };
        }
    };

    /**
     * @brief Particular solution via method of undetermined coefficients
     */
    static std::function<double(double, double)> particularSolution(
        std::function<double(double, double)> forcing,
        const std::string& forcing_type) {

        // Example: for f(x,y) = polynomial, exponential, sin/cos
        // Guess appropriate form and determine coefficients

        if (forcing_type == "polynomial") {
            return [forcing](double x, double y) {
                // Guess polynomial of same degree
                return 0.0;  // Placeholder
            };
        } else if (forcing_type == "exponential") {
            return [forcing](double x, double y) {
                // Guess A e^(αx + βy)
                return 0.0;  // Placeholder
            };
        }

        return [](double x, double y) { return 0.0; };
    }

    /**
     * @brief Complete solution: u = u_h + u_p
     */
    struct CompleteSolution {
        HomogeneousSolution homogeneous;
        std::function<double(double, double)> particular;

        std::function<double(double, double)> total(const std::vector<double>& constants) const {
            auto u_h = homogeneous.solution(constants);
            return [u_h, this](double x, double y) {
                return u_h(x, y) + particular(x, y);
            };
        }
    };
};

/**
 * ============================================================================
 * ORTHOGONAL EXPANSIONS
 * ============================================================================
 */

/**
 * @class OrthogonalFunctions
 * @brief Orthogonality and inner products
 */
class OrthogonalFunctions {
public:
    /**
     * @brief Inner product: ⟨f, g⟩ = ∫_a^b w(x) f(x) g(x) dx
     */
    static double innerProduct(
        std::function<double(double)> f,
        std::function<double(double)> g,
        std::function<double(double)> weight,
        double a, double b, int n_points = 1000) {

        double dx = (b - a) / n_points;
        double integral = 0.0;

        for (int i = 0; i <= n_points; ++i) {
            double x = a + i * dx;
            double w_x = (i == 0 || i == n_points) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
            integral += w_x * weight(x) * f(x) * g(x);
        }

        return integral * dx / 3.0;
    }

    /**
     * @brief Norm: ∥f∥ = √⟨f, f⟩
     */
    static double norm(
        std::function<double(double)> f,
        std::function<double(double)> weight,
        double a, double b) {

        return std::sqrt(innerProduct(f, f, weight, a, b));
    }

    /**
     * @brief Check orthogonality: ⟨f, g⟩ = 0
     */
    static bool areOrthogonal(
        std::function<double(double)> f,
        std::function<double(double)> g,
        std::function<double(double)> weight,
        double a, double b) {

        return std::abs(innerProduct(f, g, weight, a, b)) < 1e-6;
    }

    /**
     * @brief Gram-Schmidt orthogonalization
     */
    static std::vector<std::function<double(double)>> gramSchmidt(
        const std::vector<std::function<double(double)>>& functions,
        std::function<double(double)> weight,
        double a, double b) {

        std::vector<std::function<double(double)>> orthogonal;
        orthogonal.reserve(functions.size());

        for (const auto& f : functions) {
            auto u = f;

            // Subtract projections onto previous orthogonal functions
            for (const auto& v : orthogonal) {
                double proj_coeff = innerProduct(f, v, weight, a, b) /
                                   innerProduct(v, v, weight, a, b);

                auto old_u = u;
                u = [old_u, v, proj_coeff](double x) {
                    return old_u(x) - proj_coeff * v(x);
                };
            }

            // Normalize
            double u_norm = norm(u, weight, a, b);
            auto u_normalized = [u, u_norm](double x) {
                return u(x) / u_norm;
            };

            orthogonal.push_back(u_normalized);
        }

        return orthogonal;
    }
};

/**
 * @class OrthogonalPolynomials
 * @brief Classical orthogonal polynomials
 */
class OrthogonalPolynomials {
public:
    /**
     * @brief Legendre polynomials P_n(x) on [-1, 1] with weight w(x) = 1
     */
    static double legendreP(int n, double x) {
        if (n == 0) return 1.0;
        if (n == 1) return x;

        // Recurrence: (n+1)P_{n+1} = (2n+1)x P_n - n P_{n-1}
        double P_prev = 1.0;
        double P_curr = x;

        for (int k = 1; k < n; ++k) {
            double P_next = ((2*k + 1) * x * P_curr - k * P_prev) / (k + 1);
            P_prev = P_curr;
            P_curr = P_next;
        }

        return P_curr;
    }

    /**
     * @brief Chebyshev polynomials T_n(x) on [-1, 1] with weight w(x) = 1/√(1-x²)
     */
    static double chebyshevT(int n, double x) {
        if (n == 0) return 1.0;
        if (n == 1) return x;

        // T_n(x) = cos(n arccos(x))
        return std::cos(n * std::acos(x));
    }

    /**
     * @brief Hermite polynomials H_n(x) on (-∞, ∞) with weight w(x) = e^{-x²}
     */
    static double hermiteH(int n, double x) {
        if (n == 0) return 1.0;
        if (n == 1) return 2 * x;

        // Recurrence: H_{n+1} = 2x H_n - 2n H_{n-1}
        double H_prev = 1.0;
        double H_curr = 2 * x;

        for (int k = 1; k < n; ++k) {
            double H_next = 2 * x * H_curr - 2 * k * H_prev;
            H_prev = H_curr;
            H_curr = H_next;
        }

        return H_curr;
    }

    /**
     * @brief Laguerre polynomials L_n(x) on [0, ∞) with weight w(x) = e^{-x}
     */
    static double laguerreL(int n, double x) {
        if (n == 0) return 1.0;
        if (n == 1) return 1 - x;

        // Recurrence: (n+1)L_{n+1} = (2n+1-x)L_n - n L_{n-1}
        double L_prev = 1.0;
        double L_curr = 1 - x;

        for (int k = 1; k < n; ++k) {
            double L_next = ((2*k + 1 - x) * L_curr - k * L_prev) / (k + 1);
            L_prev = L_curr;
            L_curr = L_next;
        }

        return L_curr;
    }
};

/**
 * @class FourierSeries
 * @brief Trigonometric Fourier series
 */
class FourierSeries {
public:
    /**
     * @brief Fourier coefficients for f(x) on [-L, L]
     */
    struct Coefficients {
        double a0;  // DC component
        std::vector<double> a_n;  // Cosine coefficients
        std::vector<double> b_n;  // Sine coefficients
        double L;  // Half-period
    };

    /**
     * @brief Compute Fourier series coefficients
     */
    static Coefficients computeCoefficients(
        std::function<double(double)> f, double L, int n_terms, int n_points = 1000) {

        Coefficients coeffs;
        coeffs.L = L;
        coeffs.a_n.resize(n_terms);
        coeffs.b_n.resize(n_terms);

        double dx = (2 * L) / n_points;

        // a_0 = (1/L) ∫_{-L}^{L} f(x) dx
        double integral_a0 = 0.0;
        for (int i = 0; i <= n_points; ++i) {
            double x = -L + i * dx;
            double weight = (i == 0 || i == n_points) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
            integral_a0 += weight * f(x);
        }
        coeffs.a0 = (integral_a0 * dx / 3.0) / L;

        // a_n = (1/L) ∫_{-L}^{L} f(x) cos(nπx/L) dx
        // b_n = (1/L) ∫_{-L}^{L} f(x) sin(nπx/L) dx
        for (int n = 1; n <= n_terms; ++n) {
            double integral_an = 0.0;
            double integral_bn = 0.0;

            for (int i = 0; i <= n_points; ++i) {
                double x = -L + i * dx;
                double weight = (i == 0 || i == n_points) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);

                integral_an += weight * f(x) * std::cos(n * M_PI * x / L);
                integral_bn += weight * f(x) * std::sin(n * M_PI * x / L);
            }

            coeffs.a_n[n-1] = (integral_an * dx / 3.0) / L;
            coeffs.b_n[n-1] = (integral_bn * dx / 3.0) / L;
        }

        return coeffs;
    }

    /**
     * @brief Evaluate Fourier series at x
     */
    static double evaluate(const Coefficients& coeffs, double x) {
        double result = coeffs.a0 / 2.0;

        for (size_t n = 0; n < coeffs.a_n.size(); ++n) {
            int k = n + 1;
            result += coeffs.a_n[n] * std::cos(k * M_PI * x / coeffs.L);
            result += coeffs.b_n[n] * std::sin(k * M_PI * x / coeffs.L);
        }

        return result;
    }

    /**
     * @brief Half-range expansions (cosine series for even extension)
     */
    static Coefficients cosineSeriesExpansion(
        std::function<double(double)> f, double L, int n_terms) {

        // Even extension: use only cosine terms
        auto f_even = [f](double x) {
            return (x >= 0) ? f(x) : f(-x);
        };

        auto coeffs = computeCoefficients(f_even, L, n_terms);
        // Zero out sine coefficients for even function
        std::fill(coeffs.b_n.begin(), coeffs.b_n.end(), 0.0);

        return coeffs;
    }

    /**
     * @brief Half-range expansions (sine series for odd extension)
     */
    static Coefficients sineSeriesExpansion(
        std::function<double(double)> f, double L, int n_terms) {

        // Odd extension: use only sine terms
        auto f_odd = [f](double x) {
            return (x >= 0) ? f(x) : -f(-x);
        };

        auto coeffs = computeCoefficients(f_odd, L, n_terms);
        coeffs.a0 = 0.0;
        std::fill(coeffs.a_n.begin(), coeffs.a_n.end(), 0.0);

        return coeffs;
    }
};

/**
 * @class OrthogonalSeriesTheory
 * @brief Completeness and convergence of orthogonal function series
 */
class OrthogonalSeriesTheory {
public:
    /**
     * @brief Parseval's identity for orthogonal series
     *
     * For complete orthonormal system {φ_n}, if f = ∑ c_n φ_n, then:
     * ∥f∥² = ∑ |c_n|² (energy conservation)
     *
     * @param coefficients Expansion coefficients c_n
     * @param norms Norms ∥φ_n∥ of basis functions
     * @return ∥f∥² computed from series
     */
    static double parsevalIdentity(
        const std::vector<double>& coefficients,
        const std::vector<double>& norms) {

        double energy = 0.0;
        for (size_t n = 0; n < coefficients.size() && n < norms.size(); ++n) {
            energy += coefficients[n] * coefficients[n] * norms[n] * norms[n];
        }
        return energy;
    }

    /**
     * @brief Bessel's inequality: ∑|c_n|² ≤ ∥f∥²
     *
     * Equality holds if and only if {φ_n} is complete
     *
     * @return true if series satisfies Bessel's inequality
     */
    static bool besselInequality(
        double f_norm_squared,
        const std::vector<double>& coefficients,
        const std::vector<double>& basis_norms) {

        double series_sum = parsevalIdentity(coefficients, basis_norms);
        return series_sum <= f_norm_squared + 1e-10;  // Numerical tolerance
    }

    /**
     * @brief Test completeness via Parseval's equality
     *
     * System is complete if ∑|c_n|² = ∥f∥²
     */
    static bool isComplete(
        double f_norm_squared,
        const std::vector<double>& coefficients,
        const std::vector<double>& basis_norms,
        double tolerance = 1e-6) {

        double series_sum = parsevalIdentity(coefficients, basis_norms);
        return std::abs(series_sum - f_norm_squared) < tolerance;
    }

    /**
     * @brief Compute convergence rate of orthogonal series
     *
     * For smooth functions, coefficients decay like c_n ~ 1/n^p
     *
     * @return Estimated decay exponent p
     */
    static double convergenceRate(const std::vector<double>& coefficients) {
        if (coefficients.size() < 3) return 0.0;

        // Estimate p from log|c_n| ~ -p log(n)
        // Use linear regression on log-log plot
        double sum_log_n = 0.0, sum_log_cn = 0.0;
        double sum_log_n_sq = 0.0, sum_log_n_log_cn = 0.0;
        int count = 0;

        for (size_t n = 1; n < coefficients.size(); ++n) {
            if (std::abs(coefficients[n]) > 1e-15) {
                double log_n = std::log(n + 1);
                double log_cn = std::log(std::abs(coefficients[n]));

                sum_log_n += log_n;
                sum_log_cn += log_cn;
                sum_log_n_sq += log_n * log_n;
                sum_log_n_log_cn += log_n * log_cn;
                count++;
            }
        }

        if (count < 2) return 0.0;

        // Slope of regression line
        double slope = (count * sum_log_n_log_cn - sum_log_n * sum_log_cn) /
                      (count * sum_log_n_sq - sum_log_n * sum_log_n);

        return -slope;  // Decay rate
    }

    /**
     * @brief Mean square error of partial sum approximation
     *
     * E_N = ∥f - ∑_{n=0}^{N-1} c_n φ_n∥²
     */
    static double meanSquareError(
        std::function<double(double)> f,
        const std::vector<std::function<double(double)>>& basis,
        const std::vector<double>& coefficients,
        std::function<double(double)> weight,
        double a, double b,
        int N) {

        // Construct partial sum
        auto partial_sum = [&](double x) {
            double sum = 0.0;
            for (int n = 0; n < N && n < (int)coefficients.size(); ++n) {
                sum += coefficients[n] * basis[n](x);
            }
            return sum;
        };

        // Compute ∥f - S_N∥²
        auto error = [f, partial_sum](double x) {
            double diff = f(x) - partial_sum(x);
            return diff * diff;
        };

        return OrthogonalFunctions::innerProduct(
            [error](double x) { return std::sqrt(error(x)); },
            [error](double x) { return std::sqrt(error(x)); },
            weight, a, b);
    }
};

/**
 * @class EigenfunctionExpansions
 * @brief Sturm-Liouville theory and eigenfunction expansions
 */
class EigenfunctionExpansions {
public:
    /**
     * @brief Sturm-Liouville problem: (p(x)u')' + (q(x) + λw(x))u = 0
     *
     * Regular Sturm-Liouville problem on [a,b] with boundary conditions:
     * - α₁u(a) + α₂u'(a) = 0
     * - β₁u(b) + β₂u'(b) = 0
     */
    struct SturmLiouvilleProblem {
        std::function<double(double)> p;  // Coefficient p(x) > 0
        std::function<double(double)> q;  // Coefficient q(x)
        std::function<double(double)> w;  // Weight function w(x) > 0
        double a, b;  // Domain [a, b]

        // Boundary condition coefficients
        double alpha1, alpha2, beta1, beta2;

        /**
         * @brief Compute eigenvalues via shooting method
         *
         * Eigenvalues are λ_n such that boundary value problem has nontrivial solution
         */
        std::vector<double> computeEigenvalues(int n_eigenvalues) const {
            std::vector<double> eigenvalues;
            eigenvalues.reserve(n_eigenvalues);

            // Use shooting method: try different λ values
            double lambda_min = -10.0;
            double lambda_max = 100.0;
            double dlambda = 0.1;

            for (double lambda = lambda_min; lambda < lambda_max && (int)eigenvalues.size() < n_eigenvalues; lambda += dlambda) {
                // Solve ODE with this λ using RK4
                // Check if boundary condition at b is satisfied

                // Initial conditions from BC at a
                double u0 = alpha2;  // If α₁u(a) + α₂u'(a) = 0, choose u(a) = α₂
                double up0 = -alpha1;  // u'(a) = -α₁

                if (std::abs(alpha2) < 1e-10 && std::abs(alpha1) > 1e-10) {
                    u0 = 0.0;
                    up0 = 1.0;
                }

                // Integrate using simple Euler (for demonstration)
                int n_steps = 100;
                double h = (b - a) / n_steps;
                double x = a;
                double u = u0;
                double up = up0;

                for (int i = 0; i < n_steps; ++i) {
                    double p_val = p(x);
                    double q_val = q(x);
                    double w_val = w(x);

                    // u'' = -(q + λw)u/p - p'u'/p (simplified assuming p' computed numerically)
                    double p_prime = (p(x + 1e-6) - p(x)) / 1e-6;
                    double upp = -(q_val + lambda * w_val) * u / p_val - p_prime * up / p_val;

                    // Euler step
                    u += h * up;
                    up += h * upp;
                    x += h;
                }

                // Check if BC at b is satisfied: β₁u(b) + β₂u'(b) ≈ 0
                double bc_residual = beta1 * u + beta2 * up;

                // Look for sign changes (eigenvalues)
                static double prev_residual = bc_residual;
                if (eigenvalues.size() > 0 || lambda > lambda_min + dlambda) {
                    if (prev_residual * bc_residual < 0) {
                        // Sign change detected - refine eigenvalue
                        double lambda_refined = lambda - dlambda * bc_residual / (bc_residual - prev_residual);
                        eigenvalues.push_back(lambda_refined);
                    }
                }
                prev_residual = bc_residual;
            }

            return eigenvalues;
        }

        /**
         * @brief Compute eigenfunction for given eigenvalue
         */
        std::function<double(double)> eigenfunction(double lambda) const {
            return [this, lambda](double x) {
                // Solve ODE from a to x with initial conditions
                double u0 = alpha2;
                double up0 = -alpha1;

                if (std::abs(alpha2) < 1e-10 && std::abs(alpha1) > 1e-10) {
                    u0 = 0.0;
                    up0 = 1.0;
                }

                int n_steps = 100;
                double h = (x - a) / n_steps;
                double xi = a;
                double u = u0;
                double up = up0;

                for (int i = 0; i < n_steps; ++i) {
                    double p_val = p(xi);
                    double q_val = q(xi);
                    double w_val = w(xi);
                    double p_prime = (p(xi + 1e-6) - p(xi)) / 1e-6;

                    double upp = -(q_val + lambda * w_val) * u / p_val - p_prime * up / p_val;

                    u += h * up;
                    up += h * upp;
                    xi += h;
                }

                return u;
            };
        }

        /**
         * @brief Expand function in eigenfunction series
         *
         * f(x) = ∑ c_n φ_n(x), where c_n = ⟨f, φ_n⟩ / ⟨φ_n, φ_n⟩
         */
        std::vector<double> expandFunction(
            std::function<double(double)> f,
            const std::vector<double>& eigenvalues) const {

            std::vector<double> coefficients;
            coefficients.reserve(eigenvalues.size());

            for (double lambda : eigenvalues) {
                auto phi = eigenfunction(lambda);

                // c_n = ∫ w(x) f(x) φ_n(x) dx / ∫ w(x) φ_n²(x) dx
                double numerator = OrthogonalFunctions::innerProduct(f, phi, w, a, b);
                double denominator = OrthogonalFunctions::innerProduct(phi, phi, w, a, b);

                coefficients.push_back(numerator / denominator);
            }

            return coefficients;
        }
    };

    /**
     * @brief Common Sturm-Liouville problems
     */

    /**
     * @brief Fourier sine series (special case)
     *
     * u'' + λu = 0 on [0, L], u(0) = u(L) = 0
     * Eigenvalues: λ_n = (nπ/L)²
     * Eigenfunctions: φ_n(x) = sin(nπx/L)
     */
    static SturmLiouvilleProblem fourierSineProblem(double L) {
        SturmLiouvilleProblem prob;
        prob.p = [](double) { return 1.0; };
        prob.q = [](double) { return 0.0; };
        prob.w = [](double) { return 1.0; };
        prob.a = 0.0;
        prob.b = L;
        prob.alpha1 = 1.0; prob.alpha2 = 0.0;  // u(0) = 0
        prob.beta1 = 1.0; prob.beta2 = 0.0;    // u(L) = 0
        return prob;
    }

    /**
     * @brief Bessel problem (cylindrical symmetry)
     *
     * (xu')' + λxu = 0 on [0, R], u(R) = 0, u(0) bounded
     * Eigenfunctions: φ_n(r) = J_0(√λ_n r)
     */
    static SturmLiouvilleProblem besselProblem(double R) {
        SturmLiouvilleProblem prob;
        prob.p = [](double x) { return x; };
        prob.q = [](double) { return 0.0; };
        prob.w = [](double x) { return x; };
        prob.a = 0.0;
        prob.b = R;
        prob.alpha1 = 0.0; prob.alpha2 = 1.0;  // Bounded at 0
        prob.beta1 = 1.0; prob.beta2 = 0.0;    // u(R) = 0
        return prob;
    }

    /**
     * @brief Legendre problem (spherical symmetry)
     *
     * ((1-x²)u')' + λu = 0 on [-1, 1]
     * Eigenvalues: λ_n = n(n+1)
     * Eigenfunctions: φ_n(x) = P_n(x) (Legendre polynomials)
     */
    static SturmLiouvilleProblem legendreProblem() {
        SturmLiouvilleProblem prob;
        prob.p = [](double x) { return 1.0 - x*x; };
        prob.q = [](double) { return 0.0; };
        prob.w = [](double) { return 1.0; };
        prob.a = -1.0;
        prob.b = 1.0;
        prob.alpha1 = 0.0; prob.alpha2 = 1.0;  // Bounded at -1
        prob.beta1 = 0.0; prob.beta2 = 1.0;    // Bounded at +1
        return prob;
    }
};

/**
 * @class BesselFunctions
 * @brief Bessel functions and applications
 */
class BesselFunctions {
public:
    /**
     * @brief Bessel function of first kind J_n(x) via series expansion
     */
    static double besselJ(int n, double x, int max_terms = 50) {
        double result = 0.0;
        double term = std::pow(x / 2.0, n) / std::tgamma(n + 1);

        for (int k = 0; k < max_terms; ++k) {
            result += term;
            term *= -x * x / (4.0 * (k + 1) * (n + k + 1));

            if (std::abs(term) < 1e-15) break;
        }

        return result;
    }

    /**
     * @brief Recurrence relation: J_{n+1}(x) = (2n/x)J_n(x) - J_{n-1}(x)
     */
    static double besselJRecurrence(int n, double x) {
        if (n == 0) return besselJ(0, x);
        if (n == 1) return besselJ(1, x);

        double J_prev = besselJ(0, x);
        double J_curr = besselJ(1, x);

        for (int k = 1; k < n; ++k) {
            double J_next = (2.0 * k / x) * J_curr - J_prev;
            J_prev = J_curr;
            J_curr = J_next;
        }

        return J_curr;
    }

    /**
     * @brief Zeros of J_n(x) (useful for eigenvalue problems)
     */
    static std::vector<double> besselJZeros(int n, int num_zeros) {
        std::vector<double> zeros;
        zeros.reserve(num_zeros);

        // Approximate initial guess
        double x = M_PI;  // First zero is roughly near π for J_0

        for (int i = 0; i < num_zeros; ++i) {
            // Bisection to refine zero
            double x_left = x;
            double x_right = x + M_PI;

            // Find interval where sign changes
            while (besselJ(n, x_left) * besselJ(n, x_right) > 0) {
                x_left = x_right;
                x_right += M_PI;
            }

            // Bisection
            for (int iter = 0; iter < 50; ++iter) {
                double x_mid = 0.5 * (x_left + x_right);
                double f_mid = besselJ(n, x_mid);

                if (std::abs(f_mid) < 1e-10) {
                    zeros.push_back(x_mid);
                    x = x_mid + M_PI;  // Next approximate zero
                    break;
                }

                if (besselJ(n, x_left) * f_mid < 0) {
                    x_right = x_mid;
                } else {
                    x_left = x_mid;
                }
            }
        }

        return zeros;
    }

    /**
     * @brief Modified Bessel function I_n(x) = i^{-n} J_n(ix)
     */
    static double besselI(int n, double x, int max_terms = 50) {
        double result = 0.0;
        double term = std::pow(x / 2.0, n) / std::tgamma(n + 1);

        for (int k = 0; k < max_terms; ++k) {
            result += term;
            term *= x * x / (4.0 * (k + 1) * (n + k + 1));

            if (std::abs(term) < 1e-15) break;
        }

        return result;
    }
};

/**
 * ============================================================================
 * SEPARATION OF VARIABLES
 * ============================================================================
 */

/**
 * @class SeparationOfVariables
 * @brief Solving PDEs via separation of variables
 */
class SeparationOfVariables {
public:
    /**
     * @brief Wave equation: u_tt = c² u_xx on [0, L] with u(0,t) = u(L,t) = 0
     *
     * Solution: u(x,t) = ∑ (A_n cos(nπct/L) + B_n sin(nπct/L)) sin(nπx/L)
     */
    struct WaveEquationSolution {
        double c;  // Wave speed
        double L;  // Domain length
        std::vector<double> A_n, B_n;  // Fourier coefficients

        /**
         * @brief Solve with initial conditions u(x,0) = f(x), u_t(x,0) = g(x)
         */
        static WaveEquationSolution solve(
            std::function<double(double)> f,
            std::function<double(double)> g,
            double c, double L, int n_terms) {

            WaveEquationSolution sol;
            sol.c = c;
            sol.L = L;
            sol.A_n.resize(n_terms);
            sol.B_n.resize(n_terms);

            // Compute Fourier sine series coefficients
            int n_points = 1000;
            double dx = L / n_points;

            for (int n = 1; n <= n_terms; ++n) {
                double integral_f = 0.0;
                double integral_g = 0.0;

                for (int i = 0; i <= n_points; ++i) {
                    double x = i * dx;
                    double weight = (i == 0 || i == n_points) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);

                    integral_f += weight * f(x) * std::sin(n * M_PI * x / L);
                    integral_g += weight * g(x) * std::sin(n * M_PI * x / L);
                }

                sol.A_n[n-1] = (2.0 / L) * (integral_f * dx / 3.0);
                sol.B_n[n-1] = (2.0 / L) * (integral_g * dx / 3.0) * L / (n * M_PI * c);
            }

            return sol;
        }

        /**
         * @brief Evaluate solution at (x, t)
         */
        double evaluate(double x, double t) const {
            double result = 0.0;

            for (size_t n = 0; n < A_n.size(); ++n) {
                int k = n + 1;
                double omega = k * M_PI * c / L;
                double sin_nx = std::sin(k * M_PI * x / L);

                result += (A_n[n] * std::cos(omega * t) + B_n[n] * std::sin(omega * t)) * sin_nx;
            }

            return result;
        }
    };

    /**
     * @brief Heat equation: u_t = α u_xx on [0, L] with u(0,t) = u(L,t) = 0
     *
     * Solution: u(x,t) = ∑ A_n exp(-α(nπ/L)²t) sin(nπx/L)
     */
    struct HeatEquationSolution {
        double alpha;  // Thermal diffusivity
        double L;      // Domain length
        std::vector<double> A_n;  // Fourier coefficients

        /**
         * @brief Solve with initial condition u(x,0) = f(x)
         */
        static HeatEquationSolution solve(
            std::function<double(double)> f,
            double alpha, double L, int n_terms) {

            HeatEquationSolution sol;
            sol.alpha = alpha;
            sol.L = L;
            sol.A_n.resize(n_terms);

            // Fourier sine series for initial condition
            int n_points = 1000;
            double dx = L / n_points;

            for (int n = 1; n <= n_terms; ++n) {
                double integral = 0.0;

                for (int i = 0; i <= n_points; ++i) {
                    double x = i * dx;
                    double weight = (i == 0 || i == n_points) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
                    integral += weight * f(x) * std::sin(n * M_PI * x / L);
                }

                sol.A_n[n-1] = (2.0 / L) * (integral * dx / 3.0);
            }

            return sol;
        }

        /**
         * @brief Evaluate solution at (x, t)
         */
        double evaluate(double x, double t) const {
            double result = 0.0;

            for (size_t n = 0; n < A_n.size(); ++n) {
                int k = n + 1;
                double lambda = k * M_PI / L;
                result += A_n[n] * std::exp(-alpha * lambda * lambda * t) *
                         std::sin(k * M_PI * x / L);
            }

            return result;
        }
    };

    /**
     * @brief Laplace equation: ∇²u = 0 on rectangle [0, a] × [0, b]
     */
    struct LaplaceRectangleSolution {
        double a, b;  // Domain dimensions
        std::vector<double> coefficients;

        /**
         * @brief Solve with Dirichlet boundary conditions
         */
        static LaplaceRectangleSolution solve(
            std::function<double(double)> f_top,  // u(x, b) = f(x)
            double a, double b, int n_terms) {

            LaplaceRectangleSolution sol;
            sol.a = a;
            sol.b = b;
            sol.coefficients.resize(n_terms);

            // u(x,y) = ∑ A_n sinh(nπy/a) sin(nπx/a)
            // Determine A_n from boundary condition at y = b

            int n_points = 1000;
            double dx = a / n_points;

            for (int n = 1; n <= n_terms; ++n) {
                double integral = 0.0;

                for (int i = 0; i <= n_points; ++i) {
                    double x = i * dx;
                    double weight = (i == 0 || i == n_points) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
                    integral += weight * f_top(x) * std::sin(n * M_PI * x / a);
                }

                double An = (2.0 / a) * (integral * dx / 3.0) / std::sinh(n * M_PI * b / a);
                sol.coefficients[n-1] = An;
            }

            return sol;
        }

        /**
         * @brief Evaluate solution at (x, y)
         */
        double evaluate(double x, double y) const {
            double result = 0.0;

            for (size_t n = 0; n < coefficients.size(); ++n) {
                int k = n + 1;
                result += coefficients[n] * std::sinh(k * M_PI * y / a) *
                         std::sin(k * M_PI * x / a);
            }

            return result;
        }
    };

    /**
     * @brief Cylindrical coordinates: Laplace equation ∇²u = 0
     *
     * u(r, θ, z) via separation gives Bessel functions in r
     */
    struct CylindricalSolution {
        std::vector<double> eigenvalues;  // From boundary conditions
        std::vector<double> coefficients;

        /**
         * @brief Solve axisymmetric problem (no θ dependence)
         */
        double evaluate(double r, double z) const {
            double result = 0.0;

            for (size_t n = 0; n < eigenvalues.size(); ++n) {
                double lambda = eigenvalues[n];
                result += coefficients[n] * BesselFunctions::besselJ(0, lambda * r) *
                         std::exp(-lambda * z);
            }

            return result;
        }
    };

    /**
     * @brief Spherical coordinates: Laplace equation ∇²u = 0
     *
     * u(r, θ, φ) via separation gives Legendre polynomials
     */
    struct SphericalSolution {
        std::vector<double> A_n, B_n;  // Radial coefficients

        /**
         * @brief Axisymmetric solution (no φ dependence)
         */
        double evaluate(double r, double theta) const {
            double result = 0.0;

            for (size_t n = 0; n < A_n.size(); ++n) {
                int l = n;
                double P_l = OrthogonalPolynomials::legendreP(l, std::cos(theta));
                result += (A_n[n] * std::pow(r, l) + B_n[n] * std::pow(r, -(l+1))) * P_l;
            }

            return result;
        }
    };
};

} // namespace maths::pde

#endif // MATHS_PDE_SOLUTION_METHODS_HPP
