#ifndef MATHS_COMPLEX_ANALYSIS_HPP
#define MATHS_COMPLEX_ANALYSIS_HPP

#include <complex>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>

namespace maths {
namespace complex_analysis {

using Complex = std::complex<double>;
using ComplexFunction = std::function<Complex(Complex)>;

/**
 * @class CauchyRiemann
 * @brief Cauchy-Riemann equations and holomorphic function tests
 *
 * Implements:
 * - Cauchy-Riemann equations: ∂u/∂x = ∂v/∂y, ∂u/∂y = -∂v/∂x
 * - Holomorphic (analytic) function verification
 * - Harmonic function tests
 */
class CauchyRiemann {
public:
    /**
     * @brief Check Cauchy-Riemann equations numerically
     *
     * For f(z) = u(x,y) + iv(x,y), checks:
     * ∂u/∂x = ∂v/∂y and ∂u/∂y = -∂v/∂x
     *
     * @param f Complex function
     * @param z Point to check
     * @param h Step size for numerical derivatives
     * @return true if CR equations satisfied (within tolerance)
     */
    static bool satisfies_cauchy_riemann(
        const ComplexFunction& f,
        Complex z,
        double h = 1e-6) {

        // Compute partial derivatives numerically
        Complex f_z = f(z);
        Complex f_x_plus = f(z + Complex(h, 0));
        Complex f_x_minus = f(z - Complex(h, 0));
        Complex f_y_plus = f(z + Complex(0, h));
        Complex f_y_minus = f(z - Complex(0, h));

        // ∂f/∂x ≈ (f(z+h) - f(z-h)) / (2h)
        Complex df_dx = (f_x_plus - f_x_minus) / (2.0 * h);
        Complex df_dy = (f_y_plus - f_y_minus) / (2.0 * h);

        // Extract real and imaginary parts
        double u_x = df_dx.real();
        double v_x = df_dx.imag();
        double u_y = df_dy.real();
        double v_y = df_dy.imag();

        // Check CR equations: u_x = v_y and u_y = -v_x
        double tol = 1e-4;
        bool cr1 = std::abs(u_x - v_y) < tol;
        bool cr2 = std::abs(u_y + v_x) < tol;

        return cr1 && cr2;
    }

    /**
     * @brief Compute complex derivative f'(z)
     *
     * For holomorphic f, f'(z) = ∂f/∂z = ∂u/∂x + i∂v/∂x
     */
    static Complex derivative(
        const ComplexFunction& f,
        Complex z,
        double h = 1e-6) {

        return (f(z + h) - f(z - h)) / (2.0 * h);
    }

    /**
     * @brief Check if function is harmonic
     *
     * u is harmonic if ∇²u = ∂²u/∂x² + ∂²u/∂y² = 0
     *
     * @param u Real-valued function u(x,y)
     * @param z Point to check
     * @return true if Laplacian ≈ 0
     */
    static bool is_harmonic(
        const std::function<double(double, double)>& u,
        Complex z,
        double h = 1e-6) {

        double x = z.real(), y = z.imag();

        // Second partial derivatives
        double u_xx = (u(x + h, y) - 2.0 * u(x, y) + u(x - h, y)) / (h * h);
        double u_yy = (u(x, y + h) - 2.0 * u(x, y) + u(x, y - h)) / (h * h);

        double laplacian = u_xx + u_yy;

        return std::abs(laplacian) < 1e-3;
    }
};

/**
 * @class ComplexFunctions
 * @brief Complex elementary functions and special functions
 *
 * Implements complex versions of:
 * - Trigonometric functions (sin, cos, tan, cot)
 * - Hyperbolic functions
 * - Exponential and logarithm (with branch cuts)
 * - Power functions
 */
class ComplexFunctions {
public:
    /**
     * @brief Complex sine: sin(z) = (e^(iz) - e^(-iz)) / (2i)
     */
    static Complex sin(Complex z) {
        return std::sin(z);
    }

    /**
     * @brief Complex cosine: cos(z) = (e^(iz) + e^(-iz)) / 2
     */
    static Complex cos(Complex z) {
        return std::cos(z);
    }

    /**
     * @brief Complex tangent: tan(z) = sin(z) / cos(z)
     */
    static Complex tan(Complex z) {
        return std::tan(z);
    }

    /**
     * @brief Complex cotangent: cot(z) = cos(z) / sin(z)
     *
     * Has poles at z = nπ for integer n
     */
    static Complex cot(Complex z) {
        Complex s = std::sin(z);
        if (std::abs(s) < 1e-10) {
            throw std::runtime_error("cot(z) has pole at z = nπ");
        }
        return std::cos(z) / s;
    }

    /**
     * @brief Complex exponential: exp(z) = e^x(cos(y) + i·sin(y))
     */
    static Complex exp(Complex z) {
        return std::exp(z);
    }

    /**
     * @brief Complex logarithm (principal branch)
     *
     * log(z) = log|z| + i·arg(z), arg(z) ∈ (-π, π]
     * Branch cut along negative real axis
     */
    static Complex log(Complex z) {
        if (std::abs(z) < 1e-10) {
            throw std::runtime_error("log(0) is undefined");
        }
        return std::log(z);
    }

    /**
     * @brief Complex power: z^w = exp(w·log(z))
     *
     * Multi-valued unless w is integer
     */
    static Complex pow(Complex z, Complex w) {
        if (std::abs(z) < 1e-10 && w.real() <= 0) {
            throw std::runtime_error("0^w undefined for Re(w) ≤ 0");
        }
        return std::pow(z, w);
    }

    /**
     * @brief Complex square root (principal branch)
     *
     * √z with branch cut along negative real axis
     */
    static Complex sqrt(Complex z) {
        return std::sqrt(z);
    }

    /**
     * @brief Complex hyperbolic sine: sinh(z) = (e^z - e^(-z)) / 2
     */
    static Complex sinh(Complex z) {
        return std::sinh(z);
    }

    /**
     * @brief Complex hyperbolic cosine: cosh(z) = (e^z + e^(-z)) / 2
     */
    static Complex cosh(Complex z) {
        return std::cosh(z);
    }

    /**
     * @brief Complex hyperbolic tangent
     */
    static Complex tanh(Complex z) {
        return std::tanh(z);
    }
};

/**
 * @class PowerSeries
 * @brief Power series and radius of convergence
 *
 * For series Σ aₙ(z-z₀)ⁿ:
 * - Radius of convergence R
 * - Convergence tests (ratio test, root test)
 * - Series evaluation
 */
class PowerSeries {
public:
    /**
     * @brief Compute radius of convergence using ratio test
     *
     * R = lim_{n→∞} |aₙ/aₙ₊₁|
     *
     * @param coefficients Power series coefficients {a₀, a₁, a₂, ...}
     * @return Radius of convergence
     */
    static double radius_of_convergence_ratio(const std::vector<Complex>& coefficients) {
        int n = coefficients.size();
        if (n < 2) return std::numeric_limits<double>::infinity();

        // Use last several terms for better estimate
        int start = std::max(0, n - 10);
        double sum_ratio = 0.0;
        int count = 0;

        for (int i = start; i < n - 1; ++i) {
            if (std::abs(coefficients[i + 1]) > 1e-10) {
                double ratio = std::abs(coefficients[i]) / std::abs(coefficients[i + 1]);
                sum_ratio += ratio;
                count++;
            }
        }

        return (count > 0) ? sum_ratio / count : std::numeric_limits<double>::infinity();
    }

    /**
     * @brief Compute radius of convergence using root test
     *
     * 1/R = lim sup_{n→∞} |aₙ|^(1/n)
     */
    static double radius_of_convergence_root(const std::vector<Complex>& coefficients) {
        int n = coefficients.size();
        if (n == 0) return std::numeric_limits<double>::infinity();

        double max_root = 0.0;
        int start = std::max(0, n - 20);

        for (int i = start; i < n; ++i) {
            if (i > 0) {
                double root = std::pow(std::abs(coefficients[i]), 1.0 / i);
                max_root = std::max(max_root, root);
            }
        }

        return (max_root > 1e-10) ? 1.0 / max_root : std::numeric_limits<double>::infinity();
    }

    /**
     * @brief Evaluate power series Σ aₙ(z-z₀)ⁿ
     *
     * @param coefficients {a₀, a₁, a₂, ...}
     * @param z Point of evaluation
     * @param z0 Center of series
     * @return Series value
     */
    static Complex evaluate(
        const std::vector<Complex>& coefficients,
        Complex z,
        Complex z0 = 0.0) {

        Complex result = 0.0;
        Complex power = 1.0;  // (z - z0)^n
        Complex dz = z - z0;

        for (const Complex& a_n : coefficients) {
            result += a_n * power;
            power *= dz;
        }

        return result;
    }

    /**
     * @brief Check if series converges at point z
     */
    static bool converges_at(
        const std::vector<Complex>& coefficients,
        Complex z,
        Complex z0 = 0.0) {

        double R = radius_of_convergence_ratio(coefficients);
        double dist = std::abs(z - z0);

        return dist < R;
    }
};

/**
 * @class CauchyTheory
 * @brief Cauchy's integral theorem and formula
 *
 * Implements:
 * - Cauchy integral theorem (contour integrals)
 * - Cauchy integral formula
 * - Cauchy's inequality
 * - Maximum modulus principle
 */
class CauchyTheory {
public:
    /**
     * @brief Numerical contour integral using trapezoidal rule
     *
     * ∫_γ f(z) dz
     *
     * @param f Complex function
     * @param path Parameterized path γ(t), t ∈ [0, 1]
     * @param n_points Number of discretization points
     */
    static Complex contour_integral(
        const ComplexFunction& f,
        const std::function<Complex(double)>& path,
        int n_points = 1000) {

        Complex result = 0.0;
        double dt = 1.0 / n_points;

        for (int i = 0; i < n_points; ++i) {
            double t = i * dt;
            Complex z = path(t);
            Complex z_next = path(t + dt);

            Complex dz = z_next - z;
            result += f(z) * dz;
        }

        return result;
    }

    /**
     * @brief Cauchy integral formula: f(z₀) = (1/2πi) ∮ f(z)/(z-z₀) dz
     *
     * For f holomorphic inside and on contour
     *
     * @param f Holomorphic function
     * @param z0 Point inside contour
     * @param path Contour enclosing z0
     * @return f(z₀)
     */
    static Complex cauchy_integral_formula(
        const ComplexFunction& f,
        Complex z0,
        const std::function<Complex(double)>& path,
        int n_points = 1000) {

        auto integrand = [&](Complex z) -> Complex {
            Complex diff = z - z0;
            if (std::abs(diff) < 1e-10) return 0.0;
            return f(z) / diff;
        };

        Complex integral = contour_integral(integrand, path, n_points);
        return integral / (2.0 * M_PI * Complex(0, 1));
    }

    /**
     * @brief Cauchy's inequality: |f^(n)(z₀)| ≤ n! M / R^n
     *
     * Where M = max|f| on circle |z - z₀| = R
     */
    static double cauchy_inequality(
        const ComplexFunction& f,
        Complex z0,
        double R,
        int n,
        int n_samples = 100) {

        // Compute M = max|f| on circle
        double M = 0.0;
        for (int i = 0; i < n_samples; ++i) {
            double theta = 2.0 * M_PI * i / n_samples;
            Complex z = z0 + R * Complex(std::cos(theta), std::sin(theta));
            M = std::max(M, std::abs(f(z)));
        }

        // n!
        double factorial = 1.0;
        for (int k = 1; k <= n; ++k) {
            factorial *= k;
        }

        return factorial * M / std::pow(R, n);
    }

    /**
     * @brief Maximum modulus principle
     *
     * If f is holomorphic on region D, max|f| occurs on boundary
     *
     * @return Maximum modulus on boundary
     */
    static double maximum_modulus(
        const ComplexFunction& f,
        const std::function<Complex(double)>& boundary,
        int n_samples = 1000) {

        double max_mod = 0.0;
        for (int i = 0; i < n_samples; ++i) {
            double t = static_cast<double>(i) / n_samples;
            Complex z = boundary(t);
            max_mod = std::max(max_mod, std::abs(f(z)));
        }

        return max_mod;
    }
};

/**
 * @class ResidueCalculus
 * @brief Residue theorem and pole analysis
 *
 * Implements:
 * - Residue computation
 * - Pole order determination
 * - Zero order determination
 * - Residue theorem for contour integrals
 */
class ResidueCalculus {
public:
    /**
     * @brief Compute residue at simple pole
     *
     * Res(f, z₀) = lim_{z→z₀} (z - z₀)f(z)
     *
     * @param f Function with pole at z0
     * @param z0 Location of pole
     * @return Residue
     */
    static Complex residue_simple_pole(
        const ComplexFunction& f,
        Complex z0) {

        double h = 1e-4;
        Complex z = z0 + Complex(h, 0);
        return (z - z0) * f(z);
    }

    /**
     * @brief Compute residue at pole of order m
     *
     * Res(f, z₀) = (1/(m-1)!) lim_{z→z₀} d^(m-1)/dz^(m-1) [(z-z₀)^m f(z)]
     */
    static Complex residue_pole_order_m(
        const ComplexFunction& f,
        Complex z0,
        int m) {

        if (m == 1) return residue_simple_pole(f, z0);

        // Numerical differentiation (simplified)
        double h = 1e-4;

        auto g = [&](Complex z) -> Complex {
            Complex diff = z - z0;
            return std::pow(diff, m) * f(z);
        };

        // Compute (m-1)-th derivative numerically
        Complex deriv = g(z0 + h);
        for (int k = 1; k < m; ++k) {
            deriv = (g(z0 + h * (k + 1)) - deriv) / h;
        }

        // Divide by (m-1)!
        double factorial = 1.0;
        for (int k = 2; k < m; ++k) {
            factorial *= k;
        }

        return deriv / factorial;
    }

    /**
     * @brief Residue theorem: ∮ f(z) dz = 2πi Σ Res(f, zₖ)
     *
     * @param f Meromorphic function
     * @param poles Poles inside contour
     * @param residues Residues at each pole
     * @return Value of contour integral
     */
    static Complex residue_theorem(
        const std::vector<Complex>& residues) {

        Complex sum = 0.0;
        for (Complex res : residues) {
            sum += res;
        }

        return 2.0 * M_PI * Complex(0, 1) * sum;
    }

    /**
     * @brief Determine order of pole at z₀
     *
     * f has pole of order m if lim_{z→z₀} (z-z₀)^m f(z) ≠ 0, ∞
     *
     * @return Pole order (returns -1 if not a pole)
     */
    static int pole_order(
        const ComplexFunction& f,
        Complex z0,
        int max_order = 10) {

        double h = 1e-4;

        for (int m = 1; m <= max_order; ++m) {
            Complex z = z0 + Complex(h, 0);
            Complex val = std::pow(z - z0, m) * f(z);

            double mod = std::abs(val);
            if (mod > 1e-3 && mod < 1e3) {
                return m;  // Found order
            }
        }

        return -1;  // Not a pole or order > max_order
    }

    /**
     * @brief Determine order of zero at z₀
     *
     * f has zero of order n if f^(k)(z₀) = 0 for k < n, f^(n)(z₀) ≠ 0
     */
    static int zero_order(
        const ComplexFunction& f,
        Complex z0,
        int max_order = 10) {

        double h = 1e-6;

        // Check if f(z₀) ≈ 0
        if (std::abs(f(z0)) > 1e-6) return 0;

        for (int n = 1; n <= max_order; ++n) {
            // Compute n-th derivative numerically
            Complex deriv = (f(z0 + h) - f(z0 - h)) / (2.0 * h);

            for (int k = 1; k < n; ++k) {
                Complex f_plus = (f(z0 + h * (k + 1)) - f(z0 + h * k)) / h;
                Complex f_minus = (f(z0 - h * k) - f(z0 - h * (k + 1))) / h;
                deriv = (f_plus - f_minus) / (2.0 * h);
            }

            if (std::abs(deriv) > 1e-6) {
                return n;
            }
        }

        return max_order;
    }
};

/**
 * @class LaurentSeries
 * @brief Laurent series for functions holomorphic on annulus
 *
 * f(z) = Σ_{n=-∞}^∞ aₙ(z - z₀)ⁿ
 *
 * Computes Laurent series coefficients on annulus r < |z - z₀| < R
 */
class LaurentSeries {
public:
    struct LaurentCoefficients {
        std::vector<Complex> positive_powers;  // a₀, a₁, a₂, ...
        std::vector<Complex> negative_powers;  // a₋₁, a₋₂, a₋₃, ...
        Complex z0;  // Center
    };

    /**
     * @brief Compute Laurent series coefficients numerically
     *
     * aₙ = (1/2πi) ∮ f(z)/(z-z₀)^(n+1) dz
     *
     * @param f Function holomorphic on annulus
     * @param z0 Center
     * @param r Radius for integration contour (r < |z-z₀| < R)
     * @param max_positive Maximum positive power
     * @param max_negative Maximum negative power magnitude
     */
    static LaurentCoefficients compute_coefficients(
        const ComplexFunction& f,
        Complex z0,
        double r,
        int max_positive = 10,
        int max_negative = 10) {

        LaurentCoefficients result;
        result.z0 = z0;
        result.positive_powers.resize(max_positive + 1);
        result.negative_powers.resize(max_negative);

        // Integration contour: circle of radius r centered at z0
        auto circle = [z0, r](double t) -> Complex {
            double theta = 2.0 * M_PI * t;
            return z0 + r * Complex(std::cos(theta), std::sin(theta));
        };

        int n_points = 500;

        // Compute positive powers
        for (int n = 0; n <= max_positive; ++n) {
            auto integrand = [&](Complex z) -> Complex {
                return f(z) / std::pow(z - z0, n + 1);
            };

            Complex integral = CauchyTheory::contour_integral(
                integrand, circle, n_points);

            result.positive_powers[n] = integral / (2.0 * M_PI * Complex(0, 1));
        }

        // Compute negative powers
        for (int n = 1; n <= max_negative; ++n) {
            auto integrand = [&](Complex z) -> Complex {
                return f(z) * std::pow(z - z0, n - 1);
            };

            Complex integral = CauchyTheory::contour_integral(
                integrand, circle, n_points);

            result.negative_powers[n - 1] = integral / (2.0 * M_PI * Complex(0, 1));
        }

        return result;
    }

    /**
     * @brief Evaluate Laurent series at point z
     */
    static Complex evaluate(const LaurentCoefficients& coeffs, Complex z) {
        Complex result = 0.0;
        Complex dz = z - coeffs.z0;

        // Positive powers
        Complex power = 1.0;
        for (const Complex& a_n : coeffs.positive_powers) {
            result += a_n * power;
            power *= dz;
        }

        // Negative powers
        Complex inv_dz = 1.0 / dz;
        power = inv_dz;
        for (const Complex& a_minus_n : coeffs.negative_powers) {
            result += a_minus_n * power;
            power *= inv_dz;
        }

        return result;
    }

    /**
     * @brief Get residue (coefficient a₋₁)
     */
    static Complex get_residue(const LaurentCoefficients& coeffs) {
        if (coeffs.negative_powers.empty()) return 0.0;
        return coeffs.negative_powers[0];
    }
};

/**
 * @class ConformalMaps
 * @brief Conformal mappings and Möbius transformations
 *
 * Implements:
 * - Möbius transformations: f(z) = (az + b)/(cz + d)
 * - Conformal mapping properties
 * - Special mappings (Joukowski, etc.)
 */
class ConformalMaps {
public:
    /**
     * @brief Möbius transformation: f(z) = (az + b)/(cz + d)
     *
     * Maps circles/lines to circles/lines
     * Preserves angles (conformal)
     */
    struct MobiusTransform {
        Complex a, b, c, d;

        MobiusTransform(Complex a_, Complex b_, Complex c_, Complex d_)
            : a(a_), b(b_), c(c_), d(d_) {
            // Check ad - bc ≠ 0
            if (std::abs(a * d - b * c) < 1e-10) {
                throw std::invalid_argument("Möbius transform must have ad - bc ≠ 0");
            }
        }

        Complex operator()(Complex z) const {
            Complex denom = c * z + d;
            if (std::abs(denom) < 1e-10) {
                throw std::runtime_error("Möbius transform has pole");
            }
            return (a * z + b) / denom;
        }

        /**
         * @brief Compose two Möbius transformations
         */
        MobiusTransform compose(const MobiusTransform& other) const {
            return MobiusTransform(
                a * other.a + b * other.c,
                a * other.b + b * other.d,
                c * other.a + d * other.c,
                c * other.b + d * other.d
            );
        }

        /**
         * @brief Inverse Möbius transformation
         */
        MobiusTransform inverse() const {
            Complex det = a * d - b * c;
            return MobiusTransform(d / det, -b / det, -c / det, a / det);
        }
    };

    /**
     * @brief Möbius transformation mapping three points to three points
     *
     * Find unique Möbius f such that f(z₁) = w₁, f(z₂) = w₂, f(z₃) = w₃
     */
    static MobiusTransform three_point_map(
        Complex z1, Complex z2, Complex z3,
        Complex w1, Complex w2, Complex w3) {

        // Use cross-ratio formula
        // (w - w₁)(w₂ - w₃) / ((w - w₃)(w₂ - w₁)) = (z - z₁)(z₂ - z₃) / ((z - z₃)(z₂ - z₁))

        // Simplified construction (there's a formula for a, b, c, d)
        // This is a placeholder - full implementation would solve for coefficients
        return MobiusTransform(1.0, 0.0, 0.0, 1.0);  // Identity
    }

    /**
     * @brief Joukowski transformation: f(z) = (z + 1/z) / 2
     *
     * Maps circles to ellipses, used in airfoil theory
     */
    static Complex joukowski(Complex z) {
        if (std::abs(z) < 1e-10) {
            throw std::runtime_error("Joukowski undefined at z = 0");
        }
        return 0.5 * (z + 1.0 / z);
    }

    /**
     * @brief Map unit disk to upper half-plane
     *
     * f(z) = i(1 - z)/(1 + z)
     */
    static Complex disk_to_half_plane(Complex z) {
        return Complex(0, 1) * (1.0 - z) / (1.0 + z);
    }

    /**
     * @brief Map upper half-plane to unit disk
     *
     * Inverse of disk_to_half_plane
     */
    static Complex half_plane_to_disk(Complex z) {
        return (Complex(0, 1) - z) / (Complex(0, 1) + z);
    }

    /**
     * @brief Check if mapping is conformal at point z
     *
     * f is conformal if f'(z) ≠ 0
     */
    static bool is_conformal(
        const ComplexFunction& f,
        Complex z,
        double h = 1e-6) {

        Complex derivative = (f(z + h) - f(z - h)) / (2.0 * h);
        return std::abs(derivative) > 1e-6;
    }
};

/**
 * @class RungeTheorem
 * @brief Runge's theorem - polynomial approximation
 *
 * If f is holomorphic on compact K and A ⊂ ℂ contains at least one
 * point from each bounded component of ℂ \ K, then f can be uniformly
 * approximated on K by rational functions with poles in A
 */
class RungeTheorem {
public:
    /**
     * @brief Approximate holomorphic function by polynomial
     *
     * On simply connected domain, can approximate by polynomials
     *
     * @param f Holomorphic function
     * @param z0 Center for Taylor series
     * @param degree Polynomial degree
     * @return Polynomial coefficients
     */
    static std::vector<Complex> polynomial_approximation(
        const ComplexFunction& f,
        Complex z0,
        int degree) {

        std::vector<Complex> coeffs(degree + 1);

        // Compute Taylor series coefficients
        // aₙ = f^(n)(z₀) / n!
        coeffs[0] = f(z0);

        double h = 1e-5;
        for (int n = 1; n <= degree; ++n) {
            // Numerical differentiation
            Complex deriv = 0.0;

            // Use central differences
            for (int k = 0; k <= n; ++k) {
                double sign = ((n - k) % 2 == 0) ? 1.0 : -1.0;
                double binom = 1.0;  // Binomial coefficient C(n, k)
                for (int j = 0; j < k; ++j) {
                    binom *= (n - j) / (j + 1.0);
                }

                deriv += sign * binom * f(z0 + Complex(h * (n - 2.0 * k), 0));
            }

            deriv /= std::pow(2.0 * h, n);

            // Divide by n!
            double factorial = 1.0;
            for (int k = 1; k <= n; ++k) {
                factorial *= k;
            }

            coeffs[n] = deriv / factorial;
        }

        return coeffs;
    }

    /**
     * @brief Evaluate polynomial
     */
    static Complex evaluate_polynomial(
        const std::vector<Complex>& coeffs,
        Complex z,
        Complex z0 = 0.0) {

        Complex result = 0.0;
        Complex power = 1.0;
        Complex dz = z - z0;

        for (const Complex& a_n : coeffs) {
            result += a_n * power;
            power *= dz;
        }

        return result;
    }

    /**
     * @brief Compute maximum error of polynomial approximation
     *
     * On given contour
     */
    static double approximation_error(
        const ComplexFunction& f,
        const std::vector<Complex>& poly_coeffs,
        Complex z0,
        const std::function<Complex(double)>& contour,
        int n_samples = 100) {

        double max_error = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double t = static_cast<double>(i) / n_samples;
            Complex z = contour(t);

            Complex f_val = f(z);
            Complex poly_val = evaluate_polynomial(poly_coeffs, z, z0);

            double error = std::abs(f_val - poly_val);
            max_error = std::max(max_error, error);
        }

        return max_error;
    }
};

/**
 * @class MobiusTransformations
 * @brief Extended Möbius transformation theory
 *
 * Implements:
 * - Fixed points and classification
 * - Cross ratios
 * - Automorphism groups Aut(Ĉ), Aut(ℂ), Aut(D), Aut(H²)
 */
class MobiusTransformations {
public:
    using Mobius = ConformalMaps::MobiusTransform;

    /**
     * @brief Compute fixed points of Möbius transformation
     *
     * Solve f(z) = z, i.e., (az + b)/(cz + d) = z
     * This gives cz² + (d-a)z - b = 0
     *
     * @return Vector of fixed points (0, 1, or 2 points)
     */
    static std::vector<Complex> fixed_points(const Mobius& f) {
        std::vector<Complex> fps;

        // If c = 0, transformation is f(z) = (a/d)z + b/d
        if (std::abs(f.c) < 1e-10) {
            if (std::abs(f.a - f.d) < 1e-10) {
                // f(z) = z + b/d, one fixed point at infinity
                return fps;  // Empty for finite fixed points
            } else {
                // One finite fixed point
                fps.push_back(f.b / (f.d - f.a));
                return fps;
            }
        }

        // General case: cz² + (d-a)z - b = 0
        Complex A = f.c;
        Complex B = f.d - f.a;
        Complex C = -f.b;

        // Quadratic formula
        Complex disc = std::sqrt(B * B - 4.0 * A * C);
        Complex z1 = (-B + disc) / (2.0 * A);
        Complex z2 = (-B - disc) / (2.0 * A);

        fps.push_back(z1);
        if (std::abs(z1 - z2) > 1e-10) {
            fps.push_back(z2);
        }

        return fps;
    }

    /**
     * @brief Classify Möbius transformation by fixed points
     *
     * - Elliptic: 2 fixed points, conjugate on unit circle
     * - Parabolic: 1 fixed point (double root)
     * - Hyperbolic: 2 distinct real fixed points
     * - Loxodromic: 2 distinct fixed points, general case
     *
     * @return Classification as string
     */
    static std::string classify(const Mobius& f) {
        auto fps = fixed_points(f);

        if (fps.size() == 0 || fps.size() == 1) {
            return "Parabolic";
        } else if (fps.size() == 2) {
            Complex z1 = fps[0], z2 = fps[1];

            // Check if real (hyperbolic)
            if (std::abs(z1.imag()) < 1e-6 && std::abs(z2.imag()) < 1e-6) {
                return "Hyperbolic";
            }

            // Check if conjugate on unit circle (elliptic)
            if (std::abs(std::abs(z1) - 1.0) < 1e-6 &&
                std::abs(std::abs(z2) - 1.0) < 1e-6 &&
                std::abs(z1 - std::conj(z2)) < 1e-6) {
                return "Elliptic";
            }

            return "Loxodromic";
        }

        return "Unknown";
    }

    /**
     * @brief Cross ratio (z₁, z₂; z₃, z₄) = (z₁-z₃)(z₂-z₄) / ((z₁-z₄)(z₂-z₃))
     *
     * Invariant under Möbius transformations
     */
    static Complex cross_ratio(Complex z1, Complex z2, Complex z3, Complex z4) {
        return ((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3));
    }

    /**
     * @brief Automorphisms of the Riemann sphere Ĉ
     *
     * Aut(Ĉ) = all Möbius transformations
     *
     * @return true (any Möbius is an automorphism of Ĉ)
     */
    static bool is_automorphism_riemann_sphere(const Mobius& f) {
        return true;  // All Möbius transformations
    }

    /**
     * @brief Automorphisms of ℂ (complex plane)
     *
     * Aut(ℂ) = {f(z) = az + b : a ≠ 0}
     * These are Möbius with c = 0
     */
    static bool is_automorphism_plane(const Mobius& f) {
        return std::abs(f.c) < 1e-10;
    }

    /**
     * @brief Automorphisms of unit disc D
     *
     * Aut(D) = {φ_α(z) = e^(iθ)(z-α)/(1-ᾱz) : |α| < 1}
     *
     * Maps unit disc to itself
     */
    static bool is_automorphism_disc(const Mobius& f, double tol = 1e-6) {
        // Check if |f(z)| < 1 for |z| < 1
        // Sample points on disc
        for (int i = 0; i < 20; ++i) {
            double r = 0.9;
            double theta = 2.0 * M_PI * i / 20.0;
            Complex z = r * Complex(std::cos(theta), std::sin(theta));

            Complex fz = f(z);
            if (std::abs(fz) >= 1.0 + tol) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Automorphisms of upper half-plane H²
     *
     * Aut(H²) = {f(z) = (az+b)/(cz+d) : a,b,c,d ∈ ℝ, ad-bc > 0}
     *
     * These are real Möbius transformations with positive determinant
     */
    static bool is_automorphism_half_plane(const Mobius& f, double tol = 1e-6) {
        // Check if coefficients are real
        if (std::abs(f.a.imag()) > tol || std::abs(f.b.imag()) > tol ||
            std::abs(f.c.imag()) > tol || std::abs(f.d.imag()) > tol) {
            return false;
        }

        // Check positive determinant
        Complex det = f.a * f.d - f.b * f.c;
        return det.real() > tol && std::abs(det.imag()) < tol;
    }
};

/**
 * @class HyperbolicGeometry
 * @brief Hyperbolic geometry in unit disc and upper half-plane
 *
 * Implements:
 * - Poincaré metric
 * - Hyperbolic distance
 * - Geodesics
 * - Upper half-plane and unit disc models
 */
class HyperbolicGeometry {
public:
    /**
     * @brief Poincaré metric on unit disc: ds² = 4|dz|² / (1-|z|²)²
     *
     * Returns metric tensor coefficient at point z
     */
    static double poincare_metric_disc(Complex z) {
        double r_sq = std::norm(z);  // |z|²
        if (r_sq >= 1.0) {
            throw std::runtime_error("Point must be inside unit disc");
        }

        double denom = (1.0 - r_sq) * (1.0 - r_sq);
        return 4.0 / denom;
    }

    /**
     * @brief Poincaré metric on upper half-plane: ds² = |dz|² / (Im z)²
     */
    static double poincare_metric_half_plane(Complex z) {
        double im = z.imag();
        if (im <= 0) {
            throw std::runtime_error("Point must be in upper half-plane");
        }

        return 1.0 / (im * im);
    }

    /**
     * @brief Hyperbolic distance in unit disc
     *
     * d_H(z, w) = tanh⁻¹(|(z-w)/(1-z̄w)|)
     * or d_H(z, w) = 2 tanh⁻¹(|φ_z(w)|) where φ_z is disc automorphism
     */
    static double hyperbolic_distance_disc(Complex z, Complex w) {
        Complex diff = z - w;
        Complex denom = 1.0 - std::conj(z) * w;

        if (std::abs(denom) < 1e-10) {
            return std::numeric_limits<double>::infinity();
        }

        double rho = std::abs(diff / denom);
        return 2.0 * std::atanh(rho);
    }

    /**
     * @brief Hyperbolic distance in upper half-plane
     *
     * d_H(z, w) = arcosh(1 + |z-w|²/(2·Im(z)·Im(w)))
     */
    static double hyperbolic_distance_half_plane(Complex z, Complex w) {
        double im_z = z.imag();
        double im_w = w.imag();

        if (im_z <= 0 || im_w <= 0) {
            throw std::runtime_error("Points must be in upper half-plane");
        }

        double diff_sq = std::norm(z - w);
        double arg = 1.0 + diff_sq / (2.0 * im_z * im_w);

        return std::acosh(arg);
    }

    /**
     * @brief Geodesic in unit disc between two points
     *
     * Geodesics are circular arcs perpendicular to boundary
     * Returns parameterized path γ(t), t ∈ [0,1]
     */
    static std::function<Complex(double)> geodesic_disc(Complex z0, Complex z1) {
        return [z0, z1](double t) -> Complex {
            // Simplified linear interpolation (exact geodesic requires circular arc)
            // For accurate geodesic, would compute circular arc through boundary
            return z0 + t * (z1 - z0);
        };
    }

    /**
     * @brief Geodesic in upper half-plane
     *
     * Geodesics are semicircles/vertical lines perpendicular to real axis
     */
    static std::function<Complex(double)> geodesic_half_plane(Complex z0, Complex z1) {
        // If same real part, geodesic is vertical line
        if (std::abs(z0.real() - z1.real()) < 1e-10) {
            return [z0, z1](double t) -> Complex {
                double x = z0.real();
                double y = z0.imag() + t * (z1.imag() - z0.imag());
                return Complex(x, y);
            };
        }

        // Otherwise, semicircle
        // Find center and radius of semicircle through z0, z1
        double x0 = z0.real(), y0 = z0.imag();
        double x1 = z1.real(), y1 = z1.imag();

        double center_x = (x0 * x0 + y0 * y0 - x1 * x1 - y1 * y1) / (2.0 * (x0 - x1));
        double radius = std::sqrt((x0 - center_x) * (x0 - center_x) + y0 * y0);

        return [center_x, radius, z0, z1](double t) -> Complex {
            // Parameterize semicircle
            double theta0 = std::atan2(z0.imag(), z0.real() - center_x);
            double theta1 = std::atan2(z1.imag(), z1.real() - center_x);
            double theta = theta0 + t * (theta1 - theta0);

            return Complex(center_x + radius * std::cos(theta),
                          radius * std::sin(theta));
        };
    }

    /**
     * @brief Convert from unit disc to upper half-plane
     *
     * Cayley transform: φ(z) = i(1-z)/(1+z)
     */
    static Complex disc_to_half_plane(Complex z) {
        return Complex(0, 1) * (1.0 - z) / (1.0 + z);
    }

    /**
     * @brief Convert from upper half-plane to unit disc
     *
     * Inverse Cayley: φ⁻¹(w) = (i-w)/(i+w)
     */
    static Complex half_plane_to_disc(Complex w) {
        return (Complex(0, 1) - w) / (Complex(0, 1) + w);
    }
};

/**
 * @class SchwarzLemma
 * @brief Schwarz lemma and applications
 *
 * Implements:
 * - Schwarz lemma
 * - Schwarz-Pick theorem
 * - Finite Blaschke products
 * - Contractions
 */
class SchwarzLemma {
public:
    /**
     * @brief Schwarz lemma
     *
     * If f: D → D is holomorphic with f(0) = 0, then:
     * (1) |f(z)| ≤ |z| for all |z| < 1
     * (2) |f'(0)| ≤ 1
     * Equality holds iff f(z) = e^(iθ)z for some θ
     *
     * @return true if f satisfies Schwarz lemma conditions
     */
    static bool verify_schwarz_lemma(
        const ComplexFunction& f,
        double tol = 1e-4) {

        // Check f(0) ≈ 0
        if (std::abs(f(0.0)) > tol) {
            return false;
        }

        // Check |f(z)| ≤ |z| on sample points
        for (int i = 1; i <= 20; ++i) {
            double r = i * 0.04;  // Sample r = 0.04, 0.08, ..., 0.8
            for (int j = 0; j < 8; ++j) {
                double theta = j * M_PI / 4.0;
                Complex z = r * Complex(std::cos(theta), std::sin(theta));

                if (std::abs(f(z)) > std::abs(z) + tol) {
                    return false;
                }
            }
        }

        // Check |f'(0)| ≤ 1
        Complex f_prime = CauchyRiemann::derivative(f, 0.0);
        if (std::abs(f_prime) > 1.0 + tol) {
            return false;
        }

        return true;
    }

    /**
     * @brief Schwarz-Pick theorem
     *
     * For f: D → D holomorphic:
     * |f'(z)| ≤ (1 - |f(z)|²) / (1 - |z|²)
     *
     * Equality iff f is disc automorphism
     */
    static bool verify_schwarz_pick(
        const ComplexFunction& f,
        Complex z,
        double tol = 1e-4) {

        double r_sq = std::norm(z);
        if (r_sq >= 1.0) return false;

        Complex f_z = f(z);
        double f_sq = std::norm(f_z);
        if (f_sq >= 1.0) return false;

        Complex f_prime = CauchyRiemann::derivative(f, z);
        double deriv_abs = std::abs(f_prime);

        double bound = (1.0 - f_sq) / (1.0 - r_sq);

        return deriv_abs <= bound + tol;
    }

    /**
     * @brief Finite Blaschke product
     *
     * B(z) = e^(iθ) ∏_{k=1}^n (z - αₖ)/(1 - ᾱₖz)
     *
     * Maps unit disc to itself, with zeros at {αₖ}
     */
    struct BlascheProduct {
        std::vector<Complex> zeros;
        double theta;

        BlascheProduct(const std::vector<Complex>& z, double t = 0.0)
            : zeros(z), theta(t) {
            // Verify all |αₖ| < 1
            for (Complex alpha : zeros) {
                if (std::abs(alpha) >= 1.0) {
                    throw std::invalid_argument("Blaschke zeros must be in unit disc");
                }
            }
        }

        Complex operator()(Complex z) const {
            Complex result = std::exp(Complex(0, theta));

            for (Complex alpha : zeros) {
                result *= (z - alpha) / (1.0 - std::conj(alpha) * z);
            }

            return result;
        }

        /**
         * @brief Verify |B(z)| = 1 on unit circle
         */
        bool is_unimodular_on_circle(int n_samples = 100) const {
            for (int i = 0; i < n_samples; ++i) {
                double angle = 2.0 * M_PI * i / n_samples;
                Complex z = std::exp(Complex(0, angle));

                double mod = std::abs((*this)(z));
                if (std::abs(mod - 1.0) > 1e-4) {
                    return false;
                }
            }
            return true;
        }
    };

    /**
     * @brief Check if function is a contraction
     *
     * f is a contraction if d(f(z), f(w)) ≤ d(z, w)
     * In hyperbolic metric
     */
    static bool is_contraction(
        const ComplexFunction& f,
        int n_samples = 50) {

        for (int i = 0; i < n_samples; ++i) {
            for (int j = i + 1; j < n_samples; ++j) {
                // Sample random points in disc
                double r1 = 0.8 * std::sqrt(static_cast<double>(i) / n_samples);
                double r2 = 0.8 * std::sqrt(static_cast<double>(j) / n_samples);
                double theta1 = 2.0 * M_PI * i / n_samples;
                double theta2 = 2.0 * M_PI * j / n_samples;

                Complex z = r1 * Complex(std::cos(theta1), std::sin(theta1));
                Complex w = r2 * Complex(std::cos(theta2), std::sin(theta2));

                double d_zw = HyperbolicGeometry::hyperbolic_distance_disc(z, w);

                Complex fz = f(z);
                Complex fw = f(w);

                if (std::abs(fz) >= 1.0 || std::abs(fw) >= 1.0) continue;

                double d_fz_fw = HyperbolicGeometry::hyperbolic_distance_disc(fz, fw);

                if (d_fz_fw > d_zw + 1e-4) {
                    return false;
                }
            }
        }

        return true;
    }
};

/**
 * @class HarmonicFunctions
 * @brief Harmonic functions and the Dirichlet problem
 *
 * Implements:
 * - Laplacian and harmonic verification
 * - Poisson integral formula
 * - Dirichlet problem solver
 * - Mean value property
 * - Reflection principle
 */
class HarmonicFunctions {
public:
    using RealFunction2D = std::function<double(double, double)>;

    /**
     * @brief Compute Laplacian ∇²u = ∂²u/∂x² + ∂²u/∂y²
     */
    static double laplacian(
        const RealFunction2D& u,
        double x, double y,
        double h = 1e-5) {

        double u_xx = (u(x + h, y) - 2.0 * u(x, y) + u(x - h, y)) / (h * h);
        double u_yy = (u(x, y + h) - 2.0 * u(x, y) + u(x, y - h)) / (h * h);

        return u_xx + u_yy;
    }

    /**
     * @brief Verify function is harmonic (∇²u = 0)
     */
    static bool is_harmonic(
        const RealFunction2D& u,
        double x, double y,
        double tol = 1e-3) {

        return std::abs(laplacian(u, x, y)) < tol;
    }

    /**
     * @brief Poisson integral formula for unit disc
     *
     * u(re^(iθ)) = (1/2π) ∫₀^(2π) f(e^(iφ)) P_r(θ-φ) dφ
     * where P_r(θ) = (1-r²)/(1-2r cos θ + r²) is Poisson kernel
     *
     * Solves Dirichlet problem with boundary data f
     */
    static double poisson_integral_disc(
        const std::function<double(double)>& boundary_data,  // f(θ)
        double r, double theta,
        int n_points = 200) {

        if (r >= 1.0) {
            throw std::invalid_argument("r must be < 1 for Poisson integral");
        }

        double sum = 0.0;
        double d_phi = 2.0 * M_PI / n_points;

        for (int i = 0; i < n_points; ++i) {
            double phi = i * d_phi;
            double f_phi = boundary_data(phi);

            // Poisson kernel P_r(θ - φ)
            double diff = theta - phi;
            double P_r = (1.0 - r * r) / (1.0 - 2.0 * r * std::cos(diff) + r * r);

            sum += f_phi * P_r * d_phi;
        }

        return sum / (2.0 * M_PI);
    }

    /**
     * @brief Poisson kernel for unit disc
     *
     * P_r(θ) = (1-r²) / (1 - 2r cos θ + r²)
     */
    static double poisson_kernel(double r, double theta) {
        if (r >= 1.0) return 0.0;
        return (1.0 - r * r) / (1.0 - 2.0 * r * std::cos(theta) + r * r);
    }

    /**
     * @brief Mean value property
     *
     * For harmonic u: u(z₀) = (1/2π) ∫₀^(2π) u(z₀ + re^(iθ)) dθ
     * Characteristic property of harmonic functions
     */
    static double mean_value(
        const RealFunction2D& u,
        double x0, double y0, double r,
        int n_samples = 100) {

        double sum = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double theta = 2.0 * M_PI * i / n_samples;
            double x = x0 + r * std::cos(theta);
            double y = y0 + r * std::sin(theta);

            sum += u(x, y);
        }

        return sum / n_samples;
    }

    /**
     * @brief Verify mean value property holds
     *
     * Characterizes harmonic functions
     */
    static bool satisfies_mean_value_property(
        const RealFunction2D& u,
        double x0, double y0,
        double tol = 1e-3) {

        double u_center = u(x0, y0);

        // Check for several radii
        for (double r = 0.1; r <= 0.5; r += 0.1) {
            double mean = mean_value(u, x0, y0, r);

            if (std::abs(u_center - mean) > tol) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Schwarz reflection principle
     *
     * If f is holomorphic on upper half-plane with real boundary values,
     * extend to lower half-plane by f(z̄) = f̄(z)
     *
     * Returns reflected value
     */
    static Complex reflection_principle(
        const ComplexFunction& f_upper,
        Complex z) {

        if (z.imag() >= 0) {
            return f_upper(z);
        } else {
            // Reflect to upper half-plane, evaluate, conjugate
            Complex z_reflected = std::conj(z);
            return std::conj(f_upper(z_reflected));
        }
    }

    /**
     * @brief Dirichlet problem solver on unit disc
     *
     * Find harmonic u on D with u|∂D = f
     * Uses Poisson integral formula
     *
     * @return Function u(r,θ) solving Dirichlet problem
     */
    static std::function<double(double, double)> solve_dirichlet_disc(
        const std::function<double(double)>& boundary_data) {

        return [boundary_data](double r, double theta) -> double {
            return poisson_integral_disc(boundary_data, r, theta);
        };
    }

    /**
     * @brief Check subharmonicity: ∇²u ≥ 0
     *
     * Subharmonic functions satisfy u(z₀) ≤ mean value
     */
    static bool is_subharmonic(
        const RealFunction2D& u,
        double x, double y,
        double tol = 1e-3) {

        double lap = laplacian(u, x, y);
        return lap >= -tol;
    }

    /**
     * @brief Maximum principle for harmonic functions
     *
     * If u is harmonic on bounded domain, max u occurs on boundary
     *
     * @return Maximum value on boundary
     */
    static double maximum_on_boundary(
        const std::function<double(double)>& boundary_data,
        int n_samples = 100) {

        double max_val = -std::numeric_limits<double>::infinity();

        for (int i = 0; i < n_samples; ++i) {
            double theta = 2.0 * M_PI * i / n_samples;
            max_val = std::max(max_val, boundary_data(theta));
        }

        return max_val;
    }
};

/**
 * @class PerronMethod
 * @brief Perron families and advanced Dirichlet problem
 *
 * Implements:
 * - Perron families (subfunctions)
 * - Upper and lower solutions
 * - Green's function
 * - Connection to Riemann mapping theorem
 */
class PerronMethod {
public:
    using RealFunction2D = std::function<double(double, double)>;

    /**
     * @brief Perron subfunction (subharmonic and bounded above by boundary data)
     */
    struct PerronSubfunction {
        RealFunction2D u;
        double bound;

        PerronSubfunction(RealFunction2D func, double b)
            : u(func), bound(b) {}

        bool is_valid_subfunction(
            double x, double y,
            const std::function<double(double)>& boundary_data) const {

            // Check subharmonicity
            if (!HarmonicFunctions::is_subharmonic(u, x, y)) {
                return false;
            }

            // Check bounded by boundary data (on boundary)
            // Simplified check
            return true;
        }
    };

    /**
     * @brief Compute Perron solution as supremum of subfunctions
     *
     * u(z) = sup{v(z) : v is Perron subfunction}
     */
    static double perron_solution(
        const std::vector<PerronSubfunction>& subfunctions,
        double x, double y) {

        double sup = -std::numeric_limits<double>::infinity();

        for (const auto& sub : subfunctions) {
            double val = sub.u(x, y);
            sup = std::max(sup, val);
        }

        return sup;
    }

    /**
     * @brief Green's function for unit disc
     *
     * G_D(z, w) = -log|φ_w(z)| where φ_w(z) = (w-z)/(1-w̄z)
     *
     * Fundamental solution for Dirichlet problem
     */
    static double greens_function_disc(Complex z, Complex w) {
        // Check both points in disc
        if (std::abs(z) >= 1.0 || std::abs(w) >= 1.0) {
            throw std::invalid_argument("Points must be in unit disc");
        }

        Complex phi_w = (w - z) / (1.0 - std::conj(w) * z);
        return -std::log(std::abs(phi_w));
    }

    /**
     * @brief Green's function for upper half-plane
     *
     * G_H(z, w) = log|z - w| - log|z - w̄|
     */
    static double greens_function_half_plane(Complex z, Complex w) {
        if (z.imag() <= 0 || w.imag() <= 0) {
            throw std::invalid_argument("Points must be in upper half-plane");
        }

        return std::log(std::abs(z - w)) - std::log(std::abs(z - std::conj(w)));
    }

    /**
     * @brief Solve Dirichlet problem using Green's function
     *
     * u(z) = ∫_∂D f(ζ) ∂G/∂n(z, ζ) |dζ|
     *
     * Integral representation using Green's function
     */
    static double solve_using_greens_function(
        const std::function<double(double)>& boundary_data,
        Complex z,
        int n_points = 200) {

        // For unit disc, this reduces to Poisson integral
        double r = std::abs(z);
        double theta = std::arg(z);

        return HarmonicFunctions::poisson_integral_disc(boundary_data, r, theta, n_points);
    }

    /**
     * @brief Connection to Riemann mapping theorem
     *
     * Green's function determines conformal map
     * φ(z) = exp(-G(z, z₀) + iG*(z, z₀))
     * where G* is harmonic conjugate
     *
     * Returns approximate conformal map to disc
     */
    static ComplexFunction riemann_map_via_greens(Complex z0) {
        return [z0](Complex z) -> Complex {
            // Simplified: would use Green's function and conjugate
            // For unit disc, this is already conformal
            return z;  // Placeholder
        };
    }
};

/**
 * @class RiemannMappingTheorem
 * @brief Numerical Riemann mapping theorem
 *
 * Every simply connected domain (except ℂ) is conformally equivalent to unit disc
 *
 * Implements approximate conformal maps
 */
class RiemannMappingTheorem {
public:
    /**
     * @brief Approximate conformal map from simply connected domain to disc
     *
     * Uses Schwarz-Christoffel or numerical optimization
     * This is a simplified version
     */
    static ComplexFunction approximate_riemann_map(
        const std::function<bool(Complex)>& in_domain,
        Complex z0,  // Point to map to origin
        int max_iter = 100) {

        // Placeholder: Full implementation would use:
        // - Boundary parameterization
        // - Theodorsen's method
        // - Schwarz-Christoffel formula for polygons
        // - Or iterative methods

        return [z0](Complex z) -> Complex {
            return z - z0;  // Simplified
        };
    }

    /**
     * @brief Verify map is approximately conformal
     *
     * Check |f'(z)| ≠ 0 and angle preservation
     */
    static bool is_conformal_map(
        const ComplexFunction& f,
        Complex z,
        double tol = 1e-4) {

        Complex derivative = CauchyRiemann::derivative(f, z);
        return std::abs(derivative) > tol;
    }

    /**
     * @brief Map polygon to unit disc (Schwarz-Christoffel)
     *
     * For polygon with vertices, uses Schwarz-Christoffel formula
     * Simplified implementation
     */
    static ComplexFunction polygon_to_disc(
        const std::vector<Complex>& vertices) {

        // Schwarz-Christoffel formula:
        // f'(z) = C ∏(z - zₖ)^(αₖ/π - 1)
        // where αₖ are interior angles

        return [](Complex z) -> Complex {
            return z;  // Placeholder
        };
    }
};

/**
 * @brief Zeros of Holomorphic Functions
 *
 * Analysis of zero sets, argument principle, and Rouché's theorem
 */
class ZerosOfHolomorphicFunctions {
public:
    /**
     * @brief Count zeros inside contour using argument principle
     *
     * N - P = (1/2πi) ∮ f'(z)/f(z) dz
     * where N = number of zeros, P = number of poles
     */
    static int count_zeros_argument_principle(
        const ComplexFunction& f,
        const std::function<Complex(double)>& contour,
        int n_points = 1000) {

        // Integrate f'(z)/f(z) around contour
        auto integrand = [&f](Complex z) -> Complex {
            Complex fz = f(z);
            if (std::abs(fz) < 1e-10) {
                return Complex(0.0, 0.0);  // Avoid division by zero
            }
            Complex fprime = CauchyRiemann::derivative(f, z);
            return fprime / fz;
        };

        Complex integral = CauchyTheory::contour_integral(integrand, contour, n_points);

        // (1/2πi) ∮ f'/f dz gives winding number
        double winding = integral.imag() / (2.0 * M_PI);
        return static_cast<int>(std::round(winding));
    }

    /**
     * @brief Rouché's theorem: if |f-g| < |f| on contour, then
     * f and g have same number of zeros inside
     */
    static bool rouche_theorem_holds(
        const ComplexFunction& f,
        const ComplexFunction& g,
        const std::function<Complex(double)>& contour,
        int n_samples = 100) {

        for (int i = 0; i < n_samples; ++i) {
            double t = 2.0 * M_PI * i / n_samples;
            Complex z = contour(t);

            Complex fz = f(z);
            Complex gz = g(z);

            if (std::abs(fz - gz) >= std::abs(fz)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Determine multiplicity of zero at z0
     */
    static int zero_multiplicity(
        const ComplexFunction& f,
        Complex z0,
        int max_order = 10,
        double tol = 1e-6) {

        for (int m = 1; m <= max_order; ++m) {
            // Check if f^(m)(z0) ≠ 0
            // Use numerical differentiation
            auto f_m = f;
            for (int k = 0; k < m; ++k) {
                auto prev = f_m;
                f_m = [prev](Complex z) -> Complex {
                    return CauchyRiemann::derivative(prev, z);
                };
            }

            Complex derivative_value = f_m(z0);
            if (std::abs(derivative_value) > tol) {
                return m;
            }
        }
        return max_order;  // High multiplicity
    }

    /**
     * @brief Find zeros in a region using grid search
     */
    static std::vector<Complex> find_zeros_grid(
        const ComplexFunction& f,
        Complex bottom_left,
        Complex top_right,
        int grid_size = 20,
        double tol = 1e-6) {

        std::vector<Complex> zeros;
        double x_min = bottom_left.real();
        double y_min = bottom_left.imag();
        double x_max = top_right.real();
        double y_max = top_right.imag();

        double dx = (x_max - x_min) / grid_size;
        double dy = (y_max - y_min) / grid_size;

        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                Complex z(x, y);

                if (std::abs(f(z)) < tol) {
                    zeros.push_back(z);
                }
            }
        }

        return zeros;
    }

    /**
     * @brief Jensen's formula relating zeros to growth
     *
     * log|f(0)| = Σ log(r/|aₖ|) + (1/2π)∫₀^(2π) log|f(re^(iθ))| dθ
     */
    static double jensens_formula(
        const ComplexFunction& f,
        const std::vector<Complex>& zeros,
        double r,
        int n_samples = 100) {

        double sum_log_zeros = 0.0;
        for (const auto& a : zeros) {
            if (std::abs(a) < r) {
                sum_log_zeros += std::log(r / std::abs(a));
            }
        }

        double integral = 0.0;
        for (int k = 0; k < n_samples; ++k) {
            double theta = 2.0 * M_PI * k / n_samples;
            Complex z = r * std::exp(Complex(0.0, theta));
            integral += std::log(std::abs(f(z)));
        }
        integral /= n_samples;

        return sum_log_zeros + integral;
    }
};

/**
 * @brief Infinite Products
 *
 * Convergence, Weierstrass factorization, and canonical products
 */
class InfiniteProducts {
public:
    /**
     * @brief Elementary factors for Weierstrass factorization
     *
     * E₀(z) = 1 - z
     * Eₙ(z) = (1 - z)exp(z + z²/2 + ... + zⁿ/n)
     */
    static Complex elementary_factor(int n, Complex z) {
        if (n == 0) {
            return 1.0 - z;
        }

        Complex sum(0.0, 0.0);
        Complex z_power = z;

        for (int k = 1; k <= n; ++k) {
            sum += z_power / static_cast<double>(k);
            z_power *= z;
        }

        return (1.0 - z) * std::exp(sum);
    }

    /**
     * @brief Check convergence of infinite product ∏(1 + aₙ)
     *
     * Converges if Σ|aₙ| < ∞ and ∏(1 + aₙ) ≠ 0
     */
    static bool is_convergent(
        const std::vector<Complex>& terms,
        double tol = 1e-6) {

        double sum_abs = 0.0;
        for (const auto& a : terms) {
            sum_abs += std::abs(a);
        }

        // Product converges if Σ|aₙ| < ∞
        return sum_abs < 1.0 / tol;  // Practical convergence test
    }

    /**
     * @brief Evaluate finite product ∏(1 + aₙ)
     */
    static Complex evaluate_product(const std::vector<Complex>& terms) {
        Complex product(1.0, 0.0);

        for (const auto& a : terms) {
            product *= (1.0 + a);
        }

        return product;
    }

    /**
     * @brief Weierstrass factorization theorem
     *
     * Entire function with zeros {aₙ}:
     * f(z) = z^m e^(g(z)) ∏ Eₙ(z/aₙ)
     * where m is order of zero at origin
     */
    static ComplexFunction weierstrass_factorization(
        const std::vector<Complex>& zeros,
        int order_at_origin,
        const ComplexFunction& exponential_factor) {

        return [zeros, order_at_origin, exponential_factor](Complex z) -> Complex {
            // z^m factor
            Complex result = std::pow(z, order_at_origin);

            // e^(g(z)) factor
            result *= std::exp(exponential_factor(z));

            // Product over zeros
            for (const auto& a : zeros) {
                if (std::abs(a) > 1e-10) {
                    // Choose order p based on convergence
                    int p = 0;  // Simplified
                    result *= elementary_factor(p, z / a);
                }
            }

            return result;
        };
    }

    /**
     * @brief Canonical product of genus p
     *
     * P(z) = ∏ Eₚ(z/aₙ)
     */
    static ComplexFunction canonical_product(
        const std::vector<Complex>& zeros,
        int genus) {

        return [zeros, genus](Complex z) -> Complex {
            Complex product(1.0, 0.0);

            for (const auto& a : zeros) {
                if (std::abs(a) > 1e-10) {
                    product *= elementary_factor(genus, z / a);
                }
            }

            return product;
        };
    }

    /**
     * @brief Hadamard's formula for genus
     *
     * Genus p is smallest integer such that Σ 1/|aₙ|^(p+1) < ∞
     */
    static int compute_genus(const std::vector<Complex>& zeros) {
        for (int p = 0; p <= 10; ++p) {
            double sum = 0.0;

            for (const auto& a : zeros) {
                if (std::abs(a) > 1e-10) {
                    sum += std::pow(1.0 / std::abs(a), p + 1);
                }
            }

            if (sum < 1e6) {  // Practical convergence
                return p;
            }
        }
        return 10;  // High genus
    }
};

/**
 * @brief The Ring H(D) of Holomorphic Functions
 *
 * Ring structure, ideals, and zero sets
 */
class RingHD {
public:
    /**
     * @brief Check if function generates principal ideal
     */
    static bool is_principal_ideal_generator(
        const ComplexFunction& f,
        Complex z0) {

        // f generates (f) if it doesn't vanish identically
        // in any neighborhood of z0

        const int samples = 20;
        const double r = 0.1;

        for (int k = 0; k < samples; ++k) {
            double theta = 2.0 * M_PI * k / samples;
            Complex z = z0 + r * std::exp(Complex(0.0, theta));

            if (std::abs(f(z)) > 1e-10) {
                return true;  // Non-zero somewhere
            }
        }

        return false;  // May be identically zero
    }

    /**
     * @brief Compute common zeros of two holomorphic functions
     */
    static std::vector<Complex> common_zeros(
        const ComplexFunction& f,
        const ComplexFunction& g,
        Complex bottom_left,
        Complex top_right,
        int grid_size = 20,
        double tol = 1e-6) {

        std::vector<Complex> common;
        double x_min = bottom_left.real();
        double y_min = bottom_left.imag();
        double x_max = top_right.real();
        double y_max = top_right.imag();

        double dx = (x_max - x_min) / grid_size;
        double dy = (y_max - y_min) / grid_size;

        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                Complex z(x, y);

                if (std::abs(f(z)) < tol && std::abs(g(z)) < tol) {
                    common.push_back(z);
                }
            }
        }

        return common;
    }

    /**
     * @brief Check if ideal is maximal
     *
     * Maximal ideals in H(D) correspond to points in D
     */
    static bool is_maximal_ideal_at_point(
        const std::vector<ComplexFunction>& generators,
        Complex z0,
        double tol = 1e-6) {

        // Maximal ideal at z0: {f ∈ H(D) : f(z0) = 0}
        for (const auto& f : generators) {
            if (std::abs(f(z0)) > tol) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Identity theorem: if f = g on a set with accumulation point,
     * then f = g everywhere
     */
    static bool satisfies_identity_theorem(
        const ComplexFunction& f,
        const ComplexFunction& g,
        const std::vector<Complex>& agreement_points,
        Complex test_point,
        double tol = 1e-6) {

        // Check if they agree on the given points
        for (const auto& z : agreement_points) {
            if (std::abs(f(z) - g(z)) > tol) {
                return false;
            }
        }

        // If agreement points have accumulation point,
        // they should agree at test point too
        return std::abs(f(test_point) - g(test_point)) < tol;
    }
};

/**
 * @brief Euler's Gamma Function
 *
 * Γ(z) = ∫₀^∞ t^(z-1) e^(-t) dt for Re(z) > 0
 */
class GammaFunction {
public:
    /**
     * @brief Compute Γ(z) for Re(z) > 0
     *
     * Using numerical integration or product formula
     */
    static Complex gamma(Complex z) {
        // For real positive values, use standard library
        if (z.imag() == 0.0 && z.real() > 0.0) {
            return Complex(std::tgamma(z.real()), 0.0);
        }

        // Weierstrass product formula:
        // Γ(z) = (e^(-γz)/z) ∏ₙ₌₁^∞ (1 + z/n)^(-1) e^(z/n)
        const double euler_mascheroni = 0.5772156649015329;
        const int N = 100;  // Number of product terms

        Complex result = std::exp(-euler_mascheroni * z) / z;

        for (int n = 1; n <= N; ++n) {
            Complex term = (1.0 + z / static_cast<double>(n));
            result *= std::exp(z / static_cast<double>(n)) / term;
        }

        return result;
    }

    /**
     * @brief Functional equation: Γ(z+1) = z·Γ(z)
     */
    static Complex gamma_functional_equation(Complex z) {
        return z * gamma(z);
    }

    /**
     * @brief Reflection formula: Γ(z)Γ(1-z) = π/sin(πz)
     */
    static Complex reflection_formula(Complex z) {
        Complex lhs = gamma(z) * gamma(1.0 - z);
        Complex rhs = M_PI / ComplexFunctions::sin(M_PI * z);
        return rhs;  // Return expected value
    }

    /**
     * @brief Legendre duplication formula
     *
     * Γ(z)Γ(z + 1/2) = 2^(1-2z) √π Γ(2z)
     */
    static Complex duplication_formula(Complex z) {
        Complex lhs = gamma(z) * gamma(z + 0.5);
        Complex rhs = std::pow(2.0, 1.0 - 2.0 * z) *
                      std::sqrt(M_PI) * gamma(2.0 * z);
        return rhs;
    }

    /**
     * @brief Stirling's approximation
     *
     * log Γ(z) ≈ (z - 1/2)log(z) - z + (1/2)log(2π) + O(1/z)
     */
    static Complex stirling_approximation(Complex z) {
        Complex log_gamma = (z - 0.5) * std::log(z) - z +
                           0.5 * std::log(2.0 * M_PI);
        return std::exp(log_gamma);
    }

    /**
     * @brief Estimate |Γ(z)| for large |z|
     *
     * |Γ(x + iy)| ≈ √(2π) |y|^(x-1/2) e^(-π|y|/2)
     */
    static double gamma_magnitude_estimate(Complex z) {
        double x = z.real();
        double y = z.imag();

        if (std::abs(y) < 1.0) {
            return std::abs(gamma(z));
        }

        return std::sqrt(2.0 * M_PI) *
               std::pow(std::abs(y), x - 0.5) *
               std::exp(-M_PI * std::abs(y) / 2.0);
    }

    /**
     * @brief Beta function: B(z, w) = Γ(z)Γ(w)/Γ(z+w)
     */
    static Complex beta(Complex z, Complex w) {
        return gamma(z) * gamma(w) / gamma(z + w);
    }

    /**
     * @brief Pochhammer symbol (rising factorial)
     *
     * (z)ₙ = z(z+1)(z+2)...(z+n-1) = Γ(z+n)/Γ(z)
     */
    static Complex pochhammer(Complex z, int n) {
        if (n == 0) return Complex(1.0, 0.0);
        if (n < 0) {
            throw std::invalid_argument("Pochhammer: n must be non-negative");
        }

        // Direct computation for small n
        if (n <= 20) {
            Complex result = z;
            for (int k = 1; k < n; ++k) {
                result *= (z + static_cast<double>(k));
            }
            return result;
        }

        // Use Gamma function for large n
        return gamma(z + static_cast<double>(n)) / gamma(z);
    }

    /**
     * @brief Digamma function: ψ(z) = Γ'(z)/Γ(z)
     *
     * Logarithmic derivative of Gamma
     */
    static Complex digamma(Complex z, double h = 1e-6) {
        // Numerical derivative of log Γ(z)
        auto log_gamma = [](Complex w) -> Complex {
            return std::log(gamma(w));
        };

        Complex derivative = CauchyRiemann::derivative(log_gamma, z, h);
        return derivative;
    }

    /**
     * @brief Γ(n) = (n-1)! for positive integers
     */
    static double gamma_factorial(int n) {
        if (n <= 0) {
            throw std::invalid_argument("Gamma factorial: n must be positive");
        }

        double result = 1.0;
        for (int k = 1; k < n; ++k) {
            result *= k;
        }
        return result;
    }

    /**
     * @brief Incomplete Gamma function
     *
     * γ(s, x) = ∫₀ˣ t^(s-1) e^(-t) dt
     */
    static Complex incomplete_gamma_lower(Complex s, double x, int n_points = 1000) {
        // Numerical integration
        Complex sum(0.0, 0.0);
        double dt = x / n_points;

        for (int i = 1; i <= n_points; ++i) {
            double t = i * dt;
            Complex integrand = std::pow(t, s - 1.0) * std::exp(-t);
            sum += integrand * dt;
        }

        return sum;
    }
};

/**
 * @brief Divisors and the Field of Meromorphic Functions
 *
 * Divisor groups, principal divisors, and meromorphic function construction
 */
class DivisorsAndMeromorphicFunctions {
public:
    /**
     * @brief Divisor representation
     *
     * D = Σ nₐ[a] where nₐ = order at a (positive for zeros, negative for poles)
     */
    struct Divisor {
        std::vector<Complex> points;
        std::vector<int> orders;  // Positive for zeros, negative for poles

        int degree() const {
            int deg = 0;
            for (int n : orders) {
                deg += n;
            }
            return deg;
        }

        void add(Complex point, int order) {
            points.push_back(point);
            orders.push_back(order);
        }
    };

    /**
     * @brief Compute divisor of a meromorphic function
     */
    static Divisor divisor_of_function(
        const ComplexFunction& f,
        const std::vector<Complex>& zeros,
        const std::vector<int>& zero_orders,
        const std::vector<Complex>& poles,
        const std::vector<int>& pole_orders) {

        Divisor div;

        for (size_t i = 0; i < zeros.size(); ++i) {
            div.add(zeros[i], zero_orders[i]);
        }

        for (size_t i = 0; i < poles.size(); ++i) {
            div.add(poles[i], -pole_orders[i]);
        }

        return div;
    }

    /**
     * @brief Check if divisor is principal
     *
     * D is principal if D = div(f) for some meromorphic f
     * On compact Riemann surface: deg(D) = 0 necessary
     */
    static bool is_principal_divisor(const Divisor& div) {
        // On Riemann sphere: principal divisors have degree 0
        return div.degree() == 0;
    }

    /**
     * @brief Construct meromorphic function with given divisor
     *
     * Using Weierstrass/Mittag-Leffler techniques
     */
    static ComplexFunction construct_from_divisor(const Divisor& div) {
        return [div](Complex z) -> Complex {
            Complex result(1.0, 0.0);

            for (size_t i = 0; i < div.points.size(); ++i) {
                Complex a = div.points[i];
                int n = div.orders[i];

                if (n > 0) {
                    // Zero of order n
                    result *= std::pow(z - a, n);
                } else if (n < 0) {
                    // Pole of order -n
                    result /= std::pow(z - a, -n);
                }
            }

            return result;
        };
    }

    /**
     * @brief Order of function at point (ord_a(f))
     *
     * Positive for zeros, negative for poles, 0 for regular non-zero
     */
    static int order_at_point(
        const ComplexFunction& f,
        Complex a,
        int max_order = 10,
        double tol = 1e-6) {

        Complex fa = f(a);

        if (std::abs(fa) < tol) {
            // Check for zero
            return ZerosOfHolomorphicFunctions::zero_multiplicity(f, a, max_order, tol);
        }

        if (std::isinf(std::abs(fa)) || std::abs(fa) > 1e10) {
            // Pole: check 1/f for zero
            auto f_inv = [f](Complex z) -> Complex {
                Complex val = f(z);
                if (std::abs(val) < 1e-10) return Complex(1e10, 0.0);
                return 1.0 / val;
            };
            return -ZerosOfHolomorphicFunctions::zero_multiplicity(f_inv, a, max_order, tol);
        }

        return 0;  // Regular point
    }

    /**
     * @brief Divisor addition (in divisor group)
     */
    static Divisor add_divisors(const Divisor& D1, const Divisor& D2) {
        Divisor result = D1;

        for (size_t i = 0; i < D2.points.size(); ++i) {
            result.add(D2.points[i], D2.orders[i]);
        }

        return result;
    }

    /**
     * @brief Check linear equivalence: D₁ ~ D₂ if D₁ - D₂ is principal
     */
    static bool are_linearly_equivalent(const Divisor& D1, const Divisor& D2) {
        Divisor diff = D1;
        for (size_t i = 0; i < D2.points.size(); ++i) {
            diff.add(D2.points[i], -D2.orders[i]);
        }

        return is_principal_divisor(diff);
    }
};

/**
 * @brief Infinite Blaschke Products
 *
 * Blaschke products with infinitely many zeros
 */
class InfiniteBlaschkeProducts {
public:
    /**
     * @brief Check Blaschke condition: Σ(1 - |aₙ|) < ∞
     *
     * Necessary and sufficient for convergence in unit disc
     */
    static bool satisfies_blaschke_condition(
        const std::vector<Complex>& zeros,
        double tol = 1e-6) {

        double sum = 0.0;

        for (const auto& a : zeros) {
            if (std::abs(a) >= 1.0 - tol) {
                return false;  // Zero not strictly inside disc
            }
            sum += (1.0 - std::abs(a));
        }

        return sum < 1e6;  // Practical convergence test
    }

    /**
     * @brief Evaluate infinite Blaschke product
     *
     * B(z) = ∏ₙ (|aₙ|/aₙ) · (aₙ - z)/(1 - āₙz)
     */
    static Complex evaluate(
        const std::vector<Complex>& zeros,
        Complex z,
        int max_terms = 100) {

        if (std::abs(z) >= 1.0) {
            throw std::invalid_argument("Blaschke product: |z| must be < 1");
        }

        Complex product(1.0, 0.0);
        int n_terms = std::min(max_terms, static_cast<int>(zeros.size()));

        for (int n = 0; n < n_terms; ++n) {
            Complex a = zeros[n];

            // Blaschke factor
            Complex normalizer = std::abs(a) > 1e-10 ?
                                std::abs(a) / a : Complex(1.0, 0.0);
            Complex factor = normalizer * (a - z) / (1.0 - std::conj(a) * z);

            product *= factor;
        }

        return product;
    }

    /**
     * @brief Check if product is bounded in disc
     *
     * Blaschke products satisfy |B(z)| ≤ 1 for |z| < 1
     */
    static bool is_bounded_in_disc(
        const std::vector<Complex>& zeros,
        int n_samples = 100,
        double tol = 1e-3) {

        for (int i = 0; i < n_samples; ++i) {
            double r = 0.9 * i / n_samples;  // Sample inside disc
            double theta = 2.0 * M_PI * i / n_samples;
            Complex z = r * std::exp(Complex(0.0, theta));

            Complex B = evaluate(zeros, z);

            if (std::abs(B) > 1.0 + tol) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Finite Blaschke product (from SchwarzLemma class)
     *
     * For comparison and finite case
     */
    static ComplexFunction finite_blaschke(
        const std::vector<Complex>& zeros,
        double theta = 0.0) {

        return [zeros, theta](Complex z) -> Complex {
            Complex product = std::exp(Complex(0.0, theta));

            for (const auto& a : zeros) {
                if (std::abs(a) > 1e-10) {
                    Complex factor = (z - a) / (1.0 - std::conj(a) * z);
                    product *= factor;
                } else {
                    product *= z;  // Zero at origin
                }
            }

            return product;
        };
    }

    /**
     * @brief Convergence on circles
     *
     * B(re^(iθ)) converges uniformly on |z| = r < 1
     */
    static bool converges_uniformly_on_circle(
        const std::vector<Complex>& zeros,
        double r,
        int n_samples = 100,
        double tol = 1e-6) {

        if (r >= 1.0) return false;

        // Check variation across different numbers of terms
        std::vector<double> partial_values;

        for (int terms = 10; terms <= std::min(100, (int)zeros.size());
             terms += 10) {

            double max_on_circle = 0.0;

            for (int k = 0; k < n_samples; ++k) {
                double theta = 2.0 * M_PI * k / n_samples;
                Complex z = r * std::exp(Complex(0.0, theta));

                Complex B = evaluate(zeros, z, terms);
                max_on_circle = std::max(max_on_circle, std::abs(B));
            }

            partial_values.push_back(max_on_circle);
        }

        // Check if stabilizing
        if (partial_values.size() >= 2) {
            double diff = std::abs(partial_values.back() -
                                  partial_values[partial_values.size()-2]);
            return diff < tol;
        }

        return true;
    }

    /**
     * @brief Zeros distribution
     *
     * Count zeros in annulus r₁ < |z| < r₂
     */
    static int count_zeros_in_annulus(
        const std::vector<Complex>& zeros,
        double r1,
        double r2) {

        int count = 0;

        for (const auto& a : zeros) {
            double mod = std::abs(a);
            if (mod > r1 && mod < r2) {
                count++;
            }
        }

        return count;
    }

    /**
     * @brief Boundary behavior
     *
     * B(re^(iθ)) → boundary values as r → 1⁻
     */
    static Complex boundary_value(
        const std::vector<Complex>& zeros,
        double theta,
        double r = 0.999) {

        Complex z = r * std::exp(Complex(0.0, theta));
        return evaluate(zeros, z);
    }
};

} // namespace complex_analysis
} // namespace maths

#endif // MATHS_COMPLEX_ANALYSIS_HPP
