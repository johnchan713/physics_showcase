/**
 * @file pde_transform_methods.hpp
 * @brief Transform Methods for Solving PDEs
 *
 * LAPLACE TRANSFORMS
 * - Definition and notation
 * - Basic Laplace transform pairs
 * - Properties (linearity, shifting, derivatives)
 * - Inverse Laplace transform
 * - Convolution theorem
 * - Application to ODEs and PDEs
 *
 * FOURIER TRANSFORMS
 * - Fourier integral representation
 * - Fourier transform pairs
 * - Properties (linearity, shifting, scaling, differentiation)
 * - Fourier sine and cosine transforms
 * - Finite Fourier transforms
 * - Parseval's theorem
 * - Application to boundary value problems
 */

#ifndef MATHS_PDE_TRANSFORM_METHODS_HPP
#define MATHS_PDE_TRANSFORM_METHODS_HPP

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
 * LAPLACE TRANSFORMS
 * ============================================================================
 */

/**
 * @class LaplaceTransform
 * @brief Laplace transform and its properties
 *
 * Definition: L{f(t)} = F(s) = ∫₀^∞ e^(-st) f(t) dt
 */
class LaplaceTransform {
public:
    /**
     * @brief Numerical Laplace transform via quadrature
     * @param f Function f(t)
     * @param s Complex frequency
     * @param T Upper limit for integration
     * @return F(s) = L{f(t)}
     */
    static std::complex<double> transform(
        std::function<double(double)> f,
        std::complex<double> s,
        double T = 20.0,
        int n_points = 1000) {

        double dt = T / n_points;
        std::complex<double> result(0.0, 0.0);

        // Trapezoidal rule
        for (int i = 0; i <= n_points; ++i) {
            double t = i * dt;
            double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
            result += weight * f(t) * std::exp(-s * t) * dt;
        }

        return result;
    }

    /**
     * @brief Standard Laplace transform pairs
     */
    struct TransformPairs {
        /**
         * @brief L{1} = 1/s, s > 0
         */
        static std::complex<double> constant(std::complex<double> s) {
            return 1.0 / s;
        }

        /**
         * @brief L{t^n} = n!/s^(n+1), s > 0
         */
        static std::complex<double> power(int n, std::complex<double> s) {
            double factorial = 1.0;
            for (int i = 1; i <= n; ++i) factorial *= i;
            return factorial / std::pow(s, n + 1);
        }

        /**
         * @brief L{e^(at)} = 1/(s-a), s > a
         */
        static std::complex<double> exponential(double a, std::complex<double> s) {
            return 1.0 / (s - a);
        }

        /**
         * @brief L{sin(ωt)} = ω/(s² + ω²), s > 0
         */
        static std::complex<double> sine(double omega, std::complex<double> s) {
            return omega / (s*s + omega*omega);
        }

        /**
         * @brief L{cos(ωt)} = s/(s² + ω²), s > 0
         */
        static std::complex<double> cosine(double omega, std::complex<double> s) {
            return s / (s*s + omega*omega);
        }

        /**
         * @brief L{sinh(at)} = a/(s² - a²), s > |a|
         */
        static std::complex<double> sinh_func(double a, std::complex<double> s) {
            return a / (s*s - a*a);
        }

        /**
         * @brief L{cosh(at)} = s/(s² - a²), s > |a|
         */
        static std::complex<double> cosh_func(double a, std::complex<double> s) {
            return s / (s*s - a*a);
        }

        /**
         * @brief L{t·e^(at)} = 1/(s-a)², s > a
         */
        static std::complex<double> t_exponential(double a, std::complex<double> s) {
            std::complex<double> denom = s - a;
            return 1.0 / (denom * denom);
        }
    };

    /**
     * @brief Properties of Laplace transform
     */
    struct Properties {
        /**
         * @brief Linearity: L{af + bg} = aL{f} + bL{g}
         */
        static std::complex<double> linearity(
            double a, std::complex<double> F,
            double b, std::complex<double> G) {
            return a * F + b * G;
        }

        /**
         * @brief First shifting theorem: L{e^(at)f(t)} = F(s-a)
         */
        static std::complex<double> firstShift(std::complex<double> F, double a, std::complex<double> s) {
            // If F is F(s), return F(s-a)
            return F;  // Symbolic representation
        }

        /**
         * @brief Second shifting theorem (time delay): L{f(t-a)u(t-a)} = e^(-as)F(s)
         */
        static std::complex<double> secondShift(std::complex<double> F, double a, std::complex<double> s) {
            return std::exp(-a * s) * F;
        }

        /**
         * @brief Transform of derivative: L{f'(t)} = sF(s) - f(0)
         */
        static std::complex<double> derivative(std::complex<double> F, std::complex<double> s, double f0) {
            return s * F - f0;
        }

        /**
         * @brief Transform of second derivative: L{f''(t)} = s²F(s) - sf(0) - f'(0)
         */
        static std::complex<double> secondDerivative(
            std::complex<double> F, std::complex<double> s,
            double f0, double f_prime_0) {
            return s*s * F - s * f0 - f_prime_0;
        }

        /**
         * @brief Transform of n-th derivative
         */
        static std::complex<double> nthDerivative(
            std::complex<double> F, std::complex<double> s,
            const std::vector<double>& initial_conditions) {

            std::complex<double> result = std::pow(s, initial_conditions.size()) * F;

            for (size_t k = 0; k < initial_conditions.size(); ++k) {
                result -= std::pow(s, initial_conditions.size() - 1 - k) * initial_conditions[k];
            }

            return result;
        }

        /**
         * @brief Transform of integral: L{∫₀ᵗ f(τ)dτ} = F(s)/s
         */
        static std::complex<double> integral(std::complex<double> F, std::complex<double> s) {
            return F / s;
        }
    };

    /**
     * @brief Convolution theorem: L{f * g} = F(s)G(s)
     */
    struct Convolution {
        /**
         * @brief Convolution of two functions: (f * g)(t) = ∫₀ᵗ f(τ)g(t-τ)dτ
         */
        static std::function<double(double)> convolve(
            std::function<double(double)> f,
            std::function<double(double)> g) {

            return [f, g](double t) {
                int n = 100;
                double dt = t / n;
                double result = 0.0;

                for (int i = 0; i <= n; ++i) {
                    double tau = i * dt;
                    double weight = (i == 0 || i == n) ? 0.5 : 1.0;
                    result += weight * f(tau) * g(t - tau) * dt;
                }

                return result;
            };
        }

        /**
         * @brief L{f * g} = L{f} · L{g}
         */
        static std::complex<double> transformOfConvolution(
            std::complex<double> F, std::complex<double> G) {
            return F * G;
        }
    };

    /**
     * @brief Inverse Laplace transform (symbolic)
     */
    struct InverseTransform {
        /**
         * @brief Partial fraction decomposition helper
         * For rational functions F(s) = P(s)/Q(s)
         */
        struct PartialFraction {
            std::vector<double> poles;
            std::vector<double> residues;
        };

        /**
         * @brief Inverse of 1/(s-a): e^(at)
         */
        static std::function<double(double)> simpleExponential(double a) {
            return [a](double t) { return std::exp(a * t); };
        }

        /**
         * @brief Inverse of 1/(s-a)^n: t^(n-1)e^(at)/(n-1)!
         */
        static std::function<double(double)> repeatedPole(double a, int n) {
            return [a, n](double t) {
                double factorial = 1.0;
                for (int i = 1; i < n; ++i) factorial *= i;
                return std::pow(t, n-1) * std::exp(a * t) / factorial;
            };
        }

        /**
         * @brief Inverse of ω/(s² + ω²): sin(ωt)
         */
        static std::function<double(double)> inverseSine(double omega) {
            return [omega](double t) { return std::sin(omega * t); };
        }

        /**
         * @brief Inverse of s/(s² + ω²): cos(ωt)
         */
        static std::function<double(double)> inverseCosine(double omega) {
            return [omega](double t) { return std::cos(omega * t); };
        }
    };
};

/**
 * ============================================================================
 * FOURIER TRANSFORMS
 * ============================================================================
 */

/**
 * @class FourierTransform
 * @brief Fourier transform and integral theorems
 *
 * Definition: F{f(x)} = F(k) = ∫₋∞^∞ f(x) e^(-ikx) dx
 * Inverse: f(x) = (1/2π) ∫₋∞^∞ F(k) e^(ikx) dk
 */
class FourierTransform {
public:
    using ComplexFunction = std::function<std::complex<double>(double)>;

    /**
     * @brief Numerical Fourier transform
     */
    static std::complex<double> transform(
        std::function<double(double)> f,
        double k,
        double x_max = 20.0,
        int n_points = 1000) {

        double dx = 2.0 * x_max / n_points;
        std::complex<double> result(0.0, 0.0);
        std::complex<double> i_unit(0.0, 1.0);

        for (int j = -n_points/2; j <= n_points/2; ++j) {
            double x = j * dx;
            double weight = (j == -n_points/2 || j == n_points/2) ? 0.5 : 1.0;
            result += weight * f(x) * std::exp(-i_unit * k * x) * dx;
        }

        return result;
    }

    /**
     * @brief Inverse Fourier transform
     */
    static double inverseTransform(
        ComplexFunction F,
        double x,
        double k_max = 20.0,
        int n_points = 1000) {

        double dk = 2.0 * k_max / n_points;
        std::complex<double> result(0.0, 0.0);
        std::complex<double> i_unit(0.0, 1.0);

        for (int j = -n_points/2; j <= n_points/2; ++j) {
            double k = j * dk;
            double weight = (j == -n_points/2 || j == n_points/2) ? 0.5 : 1.0;
            result += weight * F(k) * std::exp(i_unit * k * x) * dk;
        }

        return result.real() / (2.0 * M_PI);
    }

    /**
     * @brief Standard Fourier transform pairs
     */
    struct TransformPairs {
        /**
         * @brief F{e^(-a|x|)} = 2a/(a² + k²), a > 0
         */
        static double doubleExponential(double a, double k) {
            return 2.0 * a / (a*a + k*k);
        }

        /**
         * @brief F{e^(-ax²)} = √(π/a) e^(-k²/4a), a > 0
         */
        static std::complex<double> gaussian(double a, double k) {
            return std::sqrt(M_PI / a) * std::exp(-k*k / (4.0*a));
        }

        /**
         * @brief F{rect(x)} = sinc(k) = sin(k)/k for rect(x) = 1 if |x|<1/2, 0 otherwise
         */
        static double rectangularPulse(double k) {
            if (std::abs(k) < 1e-10) return 1.0;
            return std::sin(k) / k;
        }

        /**
         * @brief F{δ(x)} = 1 (Dirac delta)
         */
        static std::complex<double> diracDelta(double k) {
            return std::complex<double>(1.0, 0.0);
        }

        /**
         * @brief F{1} = 2πδ(k) (constant function)
         */
        static double constant() {
            return 2.0 * M_PI;
        }
    };

    /**
     * @brief Properties of Fourier transform
     */
    struct Properties {
        /**
         * @brief Linearity: F{af + bg} = aF{f} + bF{g}
         */
        static std::complex<double> linearity(
            std::complex<double> a, std::complex<double> F_f,
            std::complex<double> b, std::complex<double> F_g) {
            return a * F_f + b * F_g;
        }

        /**
         * @brief Time shifting: F{f(x-a)} = e^(-ika)F(k)
         */
        static std::complex<double> timeShift(std::complex<double> F_k, double k, double a) {
            std::complex<double> i_unit(0.0, 1.0);
            return std::exp(-i_unit * k * a) * F_k;
        }

        /**
         * @brief Frequency shifting: F{e^(iax)f(x)} = F(k-a)
         */
        static std::complex<double> frequencyShift(std::complex<double> F_k, double a) {
            // Symbolic: F(k-a)
            return F_k;
        }

        /**
         * @brief Scaling: F{f(ax)} = (1/|a|)F(k/a)
         */
        static std::complex<double> scaling(std::complex<double> F_k, double a) {
            return F_k / std::abs(a);
        }

        /**
         * @brief Differentiation in time: F{f'(x)} = ikF(k)
         */
        static std::complex<double> differentiation(std::complex<double> F_k, double k) {
            std::complex<double> i_unit(0.0, 1.0);
            return i_unit * k * F_k;
        }

        /**
         * @brief n-th derivative: F{f^(n)(x)} = (ik)^n F(k)
         */
        static std::complex<double> nthDerivative(std::complex<double> F_k, double k, int n) {
            std::complex<double> i_unit(0.0, 1.0);
            return std::pow(i_unit * k, n) * F_k;
        }

        /**
         * @brief Differentiation in frequency: F{xf(x)} = i F'(k)
         */
        static std::complex<double> multiplicationByX(std::complex<double> F_prime_k) {
            std::complex<double> i_unit(0.0, 1.0);
            return i_unit * F_prime_k;
        }
    };

    /**
     * @brief Parseval's theorem: ∫ |f(x)|² dx = (1/2π) ∫ |F(k)|² dk
     */
    static double parsevalEnergy(std::function<double(double)> f, double x_max = 20.0) {
        int n = 1000;
        double dx = 2.0 * x_max / n;
        double energy = 0.0;

        for (int i = -n/2; i <= n/2; ++i) {
            double x = i * dx;
            energy += f(x) * f(x) * dx;
        }

        return energy;
    }

    /**
     * @brief Convolution theorem: F{f * g} = F{f} · F{g}
     */
    static std::complex<double> convolutionTheorem(
        std::complex<double> F_f, std::complex<double> F_g) {
        return F_f * F_g;
    }
};

/**
 * @class FourierSineCosineTransforms
 * @brief Fourier sine and cosine transforms for semi-infinite domains
 */
class FourierSineCosineTransforms {
public:
    /**
     * @brief Fourier sine transform: Fs{f(x)} = ∫₀^∞ f(x) sin(kx) dx
     * Used for odd extensions of f(x) on [0,∞)
     */
    static double sineTransform(
        std::function<double(double)> f,
        double k,
        double x_max = 20.0,
        int n_points = 1000) {

        double dx = x_max / n_points;
        double result = 0.0;

        for (int i = 0; i <= n_points; ++i) {
            double x = i * dx;
            double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
            result += weight * f(x) * std::sin(k * x) * dx;
        }

        return result;
    }

    /**
     * @brief Fourier cosine transform: Fc{f(x)} = ∫₀^∞ f(x) cos(kx) dx
     * Used for even extensions of f(x) on [0,∞)
     */
    static double cosineTransform(
        std::function<double(double)> f,
        double k,
        double x_max = 20.0,
        int n_points = 1000) {

        double dx = x_max / n_points;
        double result = 0.0;

        for (int i = 0; i <= n_points; ++i) {
            double x = i * dx;
            double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
            result += weight * f(x) * std::cos(k * x) * dx;
        }

        return result;
    }

    /**
     * @brief Inverse sine transform: f(x) = (2/π) ∫₀^∞ Fs(k) sin(kx) dk
     */
    static double inverseSineTransform(
        std::function<double(double)> Fs,
        double x,
        double k_max = 20.0,
        int n_points = 1000) {

        double dk = k_max / n_points;
        double result = 0.0;

        for (int i = 0; i <= n_points; ++i) {
            double k = i * dk;
            double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
            result += weight * Fs(k) * std::sin(k * x) * dk;
        }

        return (2.0 / M_PI) * result;
    }

    /**
     * @brief Inverse cosine transform: f(x) = (2/π) ∫₀^∞ Fc(k) cos(kx) dk
     */
    static double inverseCosineTransform(
        std::function<double(double)> Fc,
        double x,
        double k_max = 20.0,
        int n_points = 1000) {

        double dk = k_max / n_points;
        double result = 0.0;

        for (int i = 0; i <= n_points; ++i) {
            double k = i * dk;
            double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
            result += weight * Fc(k) * std::cos(k * x) * dk;
        }

        return (2.0 / M_PI) * result;
    }

    /**
     * @brief Properties specific to sine transform
     */
    struct SineProperties {
        /**
         * @brief Fs{f'(x)} = -kFc{f(x)} + f(0)
         */
        static double derivative(double k, double Fc_f, double f0) {
            return -k * Fc_f + f0;
        }

        /**
         * @brief Fs{f''(x)} = -k²Fs{f(x)} - kf(0)
         */
        static double secondDerivative(double k, double Fs_f, double f0) {
            return -k*k * Fs_f - k * f0;
        }
    };

    /**
     * @brief Properties specific to cosine transform
     */
    struct CosineProperties {
        /**
         * @brief Fc{f'(x)} = kFs{f(x)} - f(0)
         */
        static double derivative(double k, double Fs_f, double f0) {
            return k * Fs_f - f0;
        }

        /**
         * @brief Fc{f''(x)} = -k²Fc{f(x)} + f'(0)
         */
        static double secondDerivative(double k, double Fc_f, double f_prime_0) {
            return -k*k * Fc_f + f_prime_0;
        }
    };
};

/**
 * @class FiniteFourierTransforms
 * @brief Fourier transforms on finite intervals [0, L]
 */
class FiniteFourierTransforms {
public:
    /**
     * @brief Finite sine transform: Fsn = ∫₀^L f(x) sin(nπx/L) dx
     */
    static std::vector<double> finiteSineTransform(
        std::function<double(double)> f,
        double L,
        int n_terms) {

        std::vector<double> coefficients(n_terms);
        int n_points = 1000;
        double dx = L / n_points;

        for (int n = 1; n <= n_terms; ++n) {
            double sum = 0.0;
            for (int i = 0; i <= n_points; ++i) {
                double x = i * dx;
                double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
                sum += weight * f(x) * std::sin(n * M_PI * x / L) * dx;
            }
            coefficients[n-1] = sum;
        }

        return coefficients;
    }

    /**
     * @brief Finite cosine transform: Fcn = ∫₀^L f(x) cos(nπx/L) dx
     */
    static std::vector<double> finiteCosineTransform(
        std::function<double(double)> f,
        double L,
        int n_terms) {

        std::vector<double> coefficients(n_terms + 1);
        int n_points = 1000;
        double dx = L / n_points;

        for (int n = 0; n <= n_terms; ++n) {
            double sum = 0.0;
            for (int i = 0; i <= n_points; ++i) {
                double x = i * dx;
                double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
                sum += weight * f(x) * std::cos(n * M_PI * x / L) * dx;
            }
            coefficients[n] = sum;
        }

        return coefficients;
    }

    /**
     * @brief Inverse finite sine transform: f(x) = (2/L) ∑ Fsn sin(nπx/L)
     */
    static std::function<double(double)> inverseFiniteSine(
        const std::vector<double>& coefficients,
        double L) {

        return [coefficients, L](double x) {
            double result = 0.0;
            for (size_t n = 1; n <= coefficients.size(); ++n) {
                result += coefficients[n-1] * std::sin(n * M_PI * x / L);
            }
            return (2.0 / L) * result;
        };
    }

    /**
     * @brief Inverse finite cosine transform: f(x) = (1/L)Fc0 + (2/L) ∑ Fcn cos(nπx/L)
     */
    static std::function<double(double)> inverseFiniteCosine(
        const std::vector<double>& coefficients,
        double L) {

        return [coefficients, L](double x) {
            double result = coefficients[0] / 2.0;
            for (size_t n = 1; n < coefficients.size(); ++n) {
                result += coefficients[n] * std::cos(n * M_PI * x / L);
            }
            return (2.0 / L) * result;
        };
    }

    /**
     * @brief Application to heat equation: u_t = α u_xx on [0,L]
     *
     * Solution using finite sine transform
     */
    struct HeatEquationSolution {
        double alpha;  // Thermal diffusivity
        double L;      // Domain length
        std::vector<double> initial_coeffs;  // Sine coefficients of initial condition

        /**
         * @brief Solve heat equation with Dirichlet BCs: u(0,t) = u(L,t) = 0
         */
        static HeatEquationSolution solve(
            std::function<double(double)> u0,
            double alpha,
            double L,
            int n_terms = 50) {

            HeatEquationSolution sol;
            sol.alpha = alpha;
            sol.L = L;
            sol.initial_coeffs = finiteSineTransform(u0, L, n_terms);
            return sol;
        }

        /**
         * @brief Evaluate solution at (x, t)
         */
        double evaluate(double x, double t) const {
            double result = 0.0;
            for (size_t n = 1; n <= initial_coeffs.size(); ++n) {
                double lambda_n = n * M_PI / L;
                result += initial_coeffs[n-1] *
                         std::exp(-alpha * lambda_n * lambda_n * t) *
                         std::sin(lambda_n * x);
            }
            return (2.0 / L) * result;
        }
    };
};

} // namespace maths::pde

#endif // MATHS_PDE_TRANSFORM_METHODS_HPP
