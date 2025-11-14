#ifndef MATHS_TRANSFORMS_FOURIER_HPP
#define MATHS_TRANSFORMS_FOURIER_HPP

#include <complex>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>

/**
 * @file fourier.hpp
 * @brief Fourier transforms and series
 *
 * Implements:
 * - Fourier series (periodic functions)
 * - Fourier transform (non-periodic)
 * - Discrete Fourier Transform (DFT)
 * - Fast Fourier Transform (FFT)
 * - Properties and applications
 */

namespace maths::transforms {

/**
 * @class FourierSeries
 * @brief Fourier series decomposition of periodic functions
 *
 * f(x) = a₀/2 + Σ[aₙ cos(nx) + bₙ sin(nx)]
 */
class FourierSeries {
public:
    /**
     * @brief Fourier series definition
     */
    static std::string definition() {
        return "Fourier Series:\n"
               "\n"
               "For periodic function f(x) with period 2π:\n"
               "\n"
               "f(x) = a₀/2 + Σₙ₌₁^∞ [aₙ cos(nx) + bₙ sin(nx)]\n"
               "\n"
               "Coefficients:\n"
               "a₀ = (1/π) ∫₋π^π f(x) dx\n"
               "aₙ = (1/π) ∫₋π^π f(x) cos(nx) dx\n"
               "bₙ = (1/π) ∫₋π^π f(x) sin(nx) dx";
    }

    /**
     * @brief Complex Fourier series
     */
    static std::string complexForm() {
        return "Complex Fourier Series:\n"
               "\n"
               "f(x) = Σₙ₌₋∞^∞ cₙ e^(inx)\n"
               "\n"
               "where cₙ = (1/2π) ∫₋π^π f(x) e^(-inx) dx\n"
               "\n"
               "Relation to real form:\n"
               "c₀ = a₀/2\n"
               "cₙ = (aₙ - ibₙ)/2  (n > 0)\n"
               "c₋ₙ = (aₙ + ibₙ)/2  (n > 0)";
    }

    /**
     * @brief Compute a₀ coefficient numerically
     */
    static double computeA0(std::function<double(double)> f,
                           double period = 2.0 * M_PI,
                           int n_samples = 1000) {
        double L = period / 2.0;
        double dx = period / n_samples;
        double sum = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double x = -L + i * dx;
            sum += f(x);
        }

        return sum * dx / L;
    }

    /**
     * @brief Compute aₙ coefficient numerically
     */
    static double computeAn(std::function<double(double)> f, int n,
                           double period = 2.0 * M_PI,
                           int n_samples = 1000) {
        double L = period / 2.0;
        double dx = period / n_samples;
        double sum = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double x = -L + i * dx;
            sum += f(x) * std::cos(n * M_PI * x / L);
        }

        return sum * dx / L;
    }

    /**
     * @brief Compute bₙ coefficient numerically
     */
    static double computeBn(std::function<double(double)> f, int n,
                           double period = 2.0 * M_PI,
                           int n_samples = 1000) {
        double L = period / 2.0;
        double dx = period / n_samples;
        double sum = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double x = -L + i * dx;
            sum += f(x) * std::sin(n * M_PI * x / L);
        }

        return sum * dx / L;
    }

    /**
     * @brief Dirichlet conditions (sufficient for convergence)
     */
    static std::string dirichletConditions() {
        return "Dirichlet Conditions (sufficient for convergence):\n"
               "\n"
               "1. f(x) is periodic\n"
               "2. f(x) is piecewise continuous\n"
               "3. f(x) has finite number of maxima and minima\n"
               "4. f(x) has finite number of discontinuities\n"
               "\n"
               "If satisfied: Fourier series converges to:\n"
               "- f(x) at continuity points\n"
               "- [f(x⁺) + f(x⁻)]/2 at discontinuities";
    }

    /**
     * @brief Parseval's theorem
     */
    static std::string parsevalTheorem() {
        return "Parseval's Theorem (energy conservation):\n"
               "\n"
               "(1/2π) ∫₋π^π |f(x)|² dx = Σₙ₌₋∞^∞ |cₙ|²\n"
               "\n"
               "Or in real form:\n"
               "(1/π) ∫₋π^π |f(x)|² dx = a₀²/2 + Σₙ₌₁^∞ (aₙ² + bₙ²)\n"
               "\n"
               "Energy in time domain = Energy in frequency domain";
    }
};

/**
 * @class FourierTransform
 * @brief Continuous Fourier transform
 *
 * F(ω) = ∫₋∞^∞ f(t) e^(-iωt) dt
 */
class FourierTransform {
public:
    /**
     * @brief Fourier transform definition
     */
    static std::string definition() {
        return "Fourier Transform:\n"
               "\n"
               "Forward transform:\n"
               "F(ω) = ∫₋∞^∞ f(t) e^(-iωt) dt\n"
               "\n"
               "Inverse transform:\n"
               "f(t) = (1/2π) ∫₋∞^∞ F(ω) e^(iωt) dω\n"
               "\n"
               "Decomposes signal into frequency components";
    }

    /**
     * @brief Properties of Fourier transform
     */
    static std::string properties() {
        return "Fourier Transform Properties:\n"
               "\n"
               "1. Linearity: ℱ{af + bg} = aℱ{f} + bℱ{g}\n"
               "2. Time shift: ℱ{f(t - t₀)} = e^(-iωt₀) F(ω)\n"
               "3. Frequency shift: ℱ{e^(iω₀t) f(t)} = F(ω - ω₀)\n"
               "4. Scaling: ℱ{f(at)} = (1/|a|) F(ω/a)\n"
               "5. Derivative: ℱ{f'(t)} = iω F(ω)\n"
               "6. Convolution: ℱ{f * g} = F(ω) G(ω)\n"
               "7. Parseval: ∫|f(t)|² dt = (1/2π) ∫|F(ω)|² dω";
    }

    /**
     * @brief Common Fourier transform pairs
     */
    static std::string commonPairs() {
        return "Common Fourier Transform Pairs:\n"
               "\n"
               "f(t)                     F(ω)\n"
               "--------------------------------------------\n"
               "δ(t)                     1\n"
               "1                        2πδ(ω)\n"
               "e^(-at) u(t)  (a>0)      1/(a + iω)\n"
               "e^(-a|t|)     (a>0)      2a/(a² + ω²)\n"
               "rect(t/a)                a sinc(ωa/2)\n"
               "e^(-t²/2)                √(2π) e^(-ω²/2)\n"
               "cos(ω₀t)                 π[δ(ω-ω₀) + δ(ω+ω₀)]\n"
               "sin(ω₀t)                 (π/i)[δ(ω-ω₀) - δ(ω+ω₀)]";
    }

    /**
     * @brief Uncertainty principle
     */
    static std::string uncertaintyPrinciple() {
        return "Uncertainty Principle:\n"
               "\n"
               "Δt · Δω ≥ 1/2\n"
               "\n"
               "where Δt is time spread, Δω is frequency spread\n"
               "\n"
               "Interpretation:\n"
               "- Localized in time → spread in frequency\n"
               "- Localized in frequency → spread in time\n"
               "- Cannot be simultaneously narrow in both!\n"
               "\n"
               "Analogous to Heisenberg uncertainty principle";
    }

    /**
     * @brief Convolution theorem
     */
    static std::string convolutionTheorem() {
        return "Convolution Theorem:\n"
               "\n"
               "Time domain: (f * g)(t) = ∫ f(τ) g(t - τ) dτ\n"
               "Frequency domain: ℱ{f * g} = F(ω) · G(ω)\n"
               "\n"
               "Dual:\n"
               "ℱ{f(t) · g(t)} = (1/2π) F(ω) * G(ω)\n"
               "\n"
               "Multiplication ↔ Convolution\n"
               "\n"
               "Application: Filtering, signal processing";
    }
};

/**
 * @class DiscreteFourierTransform
 * @brief Discrete Fourier Transform (DFT)
 *
 * For discrete signals (digital signal processing)
 */
class DiscreteFourierTransform {
public:
    /**
     * @brief DFT definition
     */
    static std::string definition() {
        return "Discrete Fourier Transform:\n"
               "\n"
               "Forward DFT:\n"
               "X[k] = Σₙ₌₀^(N-1) x[n] e^(-i2πkn/N)\n"
               "\n"
               "Inverse DFT:\n"
               "x[n] = (1/N) Σₖ₌₀^(N-1) X[k] e^(i2πkn/N)\n"
               "\n"
               "where N = number of samples\n"
               "      k = frequency index\n"
               "      n = time index";
    }

    /**
     * @brief Compute DFT (naive O(N²) algorithm)
     */
    static std::vector<std::complex<double>> compute(
        const std::vector<std::complex<double>>& x) {

        size_t N = x.size();
        std::vector<std::complex<double>> X(N);

        for (size_t k = 0; k < N; ++k) {
            std::complex<double> sum(0.0, 0.0);
            for (size_t n = 0; n < N; ++n) {
                double angle = -2.0 * M_PI * k * n / N;
                std::complex<double> w(std::cos(angle), std::sin(angle));
                sum += x[n] * w;
            }
            X[k] = sum;
        }

        return X;
    }

    /**
     * @brief Compute inverse DFT
     */
    static std::vector<std::complex<double>> inverse(
        const std::vector<std::complex<double>>& X) {

        size_t N = X.size();
        std::vector<std::complex<double>> x(N);

        for (size_t n = 0; n < N; ++n) {
            std::complex<double> sum(0.0, 0.0);
            for (size_t k = 0; k < N; ++k) {
                double angle = 2.0 * M_PI * k * n / N;
                std::complex<double> w(std::cos(angle), std::sin(angle));
                sum += X[k] * w;
            }
            x[n] = sum / static_cast<double>(N);
        }

        return x;
    }

    /**
     * @brief Magnitude spectrum
     */
    static std::vector<double> magnitudeSpectrum(
        const std::vector<std::complex<double>>& X) {

        std::vector<double> magnitude(X.size());
        for (size_t i = 0; i < X.size(); ++i) {
            magnitude[i] = std::abs(X[i]);
        }
        return magnitude;
    }

    /**
     * @brief Phase spectrum
     */
    static std::vector<double> phaseSpectrum(
        const std::vector<std::complex<double>>& X) {

        std::vector<double> phase(X.size());
        for (size_t i = 0; i < X.size(); ++i) {
            phase[i] = std::arg(X[i]);
        }
        return phase;
    }

    /**
     * @brief Sampling theorem (Nyquist)
     */
    static std::string nyquistTheorem() {
        return "Nyquist-Shannon Sampling Theorem:\n"
               "\n"
               "To avoid aliasing:\n"
               "f_sample ≥ 2 × f_max\n"
               "\n"
               "where f_max is highest frequency in signal\n"
               "\n"
               "Nyquist frequency: f_N = f_sample / 2\n"
               "\n"
               "Frequencies above f_N will alias!\n"
               "\n"
               "Example: CD audio samples at 44.1 kHz\n"
               "         → can represent up to 22.05 kHz";
    }
};

/**
 * @class FastFourierTransform
 * @brief FFT algorithm (Cooley-Tukey)
 *
 * O(N log N) instead of O(N²)
 */
class FastFourierTransform {
public:
    /**
     * @brief FFT algorithm description
     */
    static std::string algorithm() {
        return "Fast Fourier Transform (FFT):\n"
               "\n"
               "Cooley-Tukey algorithm:\n"
               "- Divide-and-conquer approach\n"
               "- Splits DFT into even/odd indices\n"
               "- Complexity: O(N log N)\n"
               "- Requires N = 2^k (power of 2)\n"
               "\n"
               "Speedup over naive DFT:\n"
               "N=1024: 100x faster\n"
               "N=1,000,000: 50,000x faster!\n"
               "\n"
               "Revolutionary algorithm (Cooley & Tukey, 1965)";
    }

    /**
     * @brief Compute FFT (recursive Cooley-Tukey)
     */
    static std::vector<std::complex<double>> compute(
        const std::vector<std::complex<double>>& x) {

        size_t N = x.size();

        // Base case
        if (N <= 1) return x;

        // Check if N is power of 2
        if ((N & (N - 1)) != 0) {
            throw std::invalid_argument("FFT requires size to be power of 2");
        }

        // Divide into even and odd
        std::vector<std::complex<double>> even(N / 2);
        std::vector<std::complex<double>> odd(N / 2);

        for (size_t i = 0; i < N / 2; ++i) {
            even[i] = x[2 * i];
            odd[i] = x[2 * i + 1];
        }

        // Recursive FFT
        std::vector<std::complex<double>> even_fft = compute(even);
        std::vector<std::complex<double>> odd_fft = compute(odd);

        // Combine
        std::vector<std::complex<double>> result(N);
        for (size_t k = 0; k < N / 2; ++k) {
            double angle = -2.0 * M_PI * k / N;
            std::complex<double> w(std::cos(angle), std::sin(angle));
            std::complex<double> t = w * odd_fft[k];

            result[k] = even_fft[k] + t;
            result[k + N/2] = even_fft[k] - t;
        }

        return result;
    }

    /**
     * @brief Applications of FFT
     */
    static std::string applications() {
        return "FFT Applications:\n"
               "\n"
               "1. Audio processing:\n"
               "   - Spectral analysis\n"
               "   - Equalization\n"
               "   - Compression (MP3)\n"
               "\n"
               "2. Image processing:\n"
               "   - JPEG compression\n"
               "   - Filtering\n"
               "   - Edge detection\n"
               "\n"
               "3. Communications:\n"
               "   - OFDM (WiFi, 4G/5G)\n"
               "   - Channel equalization\n"
               "\n"
               "4. Scientific computing:\n"
               "   - Solving PDEs\n"
               "   - Convolution\n"
               "   - Correlation analysis";
    }
};

/**
 * @class WindowFunctions
 * @brief Window functions for spectral analysis
 */
class WindowFunctions {
public:
    /**
     * @brief Rectangular window
     */
    static std::vector<double> rectangular(size_t N) {
        return std::vector<double>(N, 1.0);
    }

    /**
     * @brief Hann window (raised cosine)
     */
    static std::vector<double> hann(size_t N) {
        std::vector<double> window(N);
        for (size_t n = 0; n < N; ++n) {
            window[n] = 0.5 * (1.0 - std::cos(2.0 * M_PI * n / (N - 1)));
        }
        return window;
    }

    /**
     * @brief Hamming window
     */
    static std::vector<double> hamming(size_t N) {
        std::vector<double> window(N);
        for (size_t n = 0; n < N; ++n) {
            window[n] = 0.54 - 0.46 * std::cos(2.0 * M_PI * n / (N - 1));
        }
        return window;
    }

    /**
     * @brief Purpose of windowing
     */
    static std::string purpose() {
        return "Window Functions:\n"
               "\n"
               "Purpose: Reduce spectral leakage in FFT\n"
               "\n"
               "Problem: Finite-length signal creates discontinuities\n"
               "Solution: Multiply by smooth window\n"
               "\n"
               "Trade-off:\n"
               "- Rectangular: narrow main lobe, high sidelobes\n"
               "- Hann/Hamming: wider main lobe, low sidelobes\n"
               "\n"
               "Choice depends on application!";
    }
};

} // namespace maths::transforms

#endif // MATHS_TRANSFORMS_FOURIER_HPP
