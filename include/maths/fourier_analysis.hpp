#ifndef MATHS_ANALYSIS_FOURIER_ANALYSIS_HPP
#define MATHS_ANALYSIS_FOURIER_ANALYSIS_HPP

#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdexcept>

/**
 * @file fourier_analysis.hpp
 * @brief Computational Fourier analysis: DFT, FFT, wavelets, operators
 *
 * Implements discrete Fourier transforms, time-frequency analysis,
 * wavelet transforms, and Hilbert space operators.
 */

namespace maths::analysis {

using Complex = std::complex<double>;
using Signal = std::vector<double>;
using ComplexSignal = std::vector<Complex>;
using Matrix = std::vector<std::vector<double>>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

/**
 * @class DiscreteFourierTransform
 * @brief Discrete Fourier Transform (DFT) and inverse
 *
 * DFT: X[k] = Σ_{n=0}^{N-1} x[n] e^{-2πikn/N}
 */
class DiscreteFourierTransform {
public:
    /**
     * @brief Compute DFT using direct formula
     *
     * Time complexity: O(N²)
     *
     * @param signal Input signal x[n]
     * @return Frequency domain X[k]
     */
    static ComplexSignal dft(const Signal& signal) {
        size_t N = signal.size();
        ComplexSignal result(N);
        const double PI = M_PI;

        for (size_t k = 0; k < N; ++k) {
            Complex sum(0.0, 0.0);
            for (size_t n = 0; n < N; ++n) {
                double angle = -2.0 * PI * k * n / N;
                sum += signal[n] * std::exp(Complex(0.0, angle));
            }
            result[k] = sum;
        }

        return result;
    }

    /**
     * @brief Compute DFT for complex input
     *
     * @param signal Complex input signal
     * @return Complex frequency domain
     */
    static ComplexSignal dft(const ComplexSignal& signal) {
        size_t N = signal.size();
        ComplexSignal result(N);
        const double PI = M_PI;

        for (size_t k = 0; k < N; ++k) {
            Complex sum(0.0, 0.0);
            for (size_t n = 0; n < N; ++n) {
                double angle = -2.0 * PI * k * n / N;
                sum += signal[n] * std::exp(Complex(0.0, angle));
            }
            result[k] = sum;
        }

        return result;
    }

    /**
     * @brief Compute inverse DFT (IDFT)
     *
     * IDFT: x[n] = (1/N) Σ_{k=0}^{N-1} X[k] e^{2πikn/N}
     *
     * @param spectrum Frequency domain X[k]
     * @return Time domain signal x[n]
     */
    static ComplexSignal idft(const ComplexSignal& spectrum) {
        size_t N = spectrum.size();
        ComplexSignal result(N);
        const double PI = M_PI;

        for (size_t n = 0; n < N; ++n) {
            Complex sum(0.0, 0.0);
            for (size_t k = 0; k < N; ++k) {
                double angle = 2.0 * PI * k * n / N;
                sum += spectrum[k] * std::exp(Complex(0.0, angle));
            }
            result[n] = sum / double(N);
        }

        return result;
    }

    /**
     * @brief Compute power spectrum |X[k]|²
     *
     * @param signal Input signal
     * @return Power spectrum
     */
    static Signal powerSpectrum(const Signal& signal) {
        ComplexSignal spectrum = dft(signal);
        Signal power(spectrum.size());

        for (size_t k = 0; k < spectrum.size(); ++k) {
            power[k] = std::norm(spectrum[k]);  // |X[k]|²
        }

        return power;
    }

    /**
     * @brief Compute magnitude spectrum |X[k]|
     *
     * @param signal Input signal
     * @return Magnitude spectrum
     */
    static Signal magnitudeSpectrum(const Signal& signal) {
        ComplexSignal spectrum = dft(signal);
        Signal magnitude(spectrum.size());

        for (size_t k = 0; k < spectrum.size(); ++k) {
            magnitude[k] = std::abs(spectrum[k]);
        }

        return magnitude;
    }

    /**
     * @brief Compute phase spectrum arg(X[k])
     *
     * @param signal Input signal
     * @return Phase spectrum in radians
     */
    static Signal phaseSpectrum(const Signal& signal) {
        ComplexSignal spectrum = dft(signal);
        Signal phase(spectrum.size());

        for (size_t k = 0; k < spectrum.size(); ++k) {
            phase[k] = std::arg(spectrum[k]);
        }

        return phase;
    }
};

/**
 * @class FastFourierTransform
 * @brief Fast Fourier Transform (FFT) using Cooley-Tukey algorithm
 *
 * Time complexity: O(N log N)
 * Requires N to be a power of 2
 */
class FastFourierTransform {
public:
    /**
     * @brief Compute FFT using Cooley-Tukey radix-2 algorithm
     *
     * @param signal Input signal (size must be power of 2)
     * @return Frequency domain
     */
    static ComplexSignal fft(const Signal& signal) {
        ComplexSignal complex_signal(signal.size());
        for (size_t i = 0; i < signal.size(); ++i) {
            complex_signal[i] = Complex(signal[i], 0.0);
        }
        return fft(complex_signal);
    }

    /**
     * @brief FFT for complex input
     *
     * @param signal Complex input signal
     * @return Complex frequency domain
     */
    static ComplexSignal fft(ComplexSignal signal) {
        size_t N = signal.size();

        // Check if power of 2
        if (N == 0 || (N & (N - 1)) != 0) {
            throw std::invalid_argument("FFT size must be power of 2");
        }

        if (N == 1) {
            return signal;
        }

        // Divide: separate even and odd indices
        ComplexSignal even(N / 2);
        ComplexSignal odd(N / 2);

        for (size_t i = 0; i < N / 2; ++i) {
            even[i] = signal[2 * i];
            odd[i] = signal[2 * i + 1];
        }

        // Conquer: recursive FFT
        ComplexSignal fft_even = fft(even);
        ComplexSignal fft_odd = fft(odd);

        // Combine
        ComplexSignal result(N);
        const double PI = M_PI;

        for (size_t k = 0; k < N / 2; ++k) {
            double angle = -2.0 * PI * k / N;
            Complex twiddle = std::exp(Complex(0.0, angle));
            Complex t = twiddle * fft_odd[k];

            result[k] = fft_even[k] + t;
            result[k + N / 2] = fft_even[k] - t;
        }

        return result;
    }

    /**
     * @brief Compute inverse FFT
     *
     * @param spectrum Frequency domain
     * @return Time domain signal
     */
    static ComplexSignal ifft(ComplexSignal spectrum) {
        size_t N = spectrum.size();

        // Conjugate input
        for (auto& s : spectrum) {
            s = std::conj(s);
        }

        // Apply FFT
        ComplexSignal result = fft(spectrum);

        // Conjugate output and scale
        for (auto& r : result) {
            r = std::conj(r) / double(N);
        }

        return result;
    }

    /**
     * @brief Zero-pad signal to next power of 2
     *
     * @param signal Input signal
     * @return Zero-padded signal
     */
    static Signal zeroPad(const Signal& signal) {
        size_t N = signal.size();
        size_t N_padded = 1;
        while (N_padded < N) {
            N_padded *= 2;
        }

        Signal padded(N_padded, 0.0);
        std::copy(signal.begin(), signal.end(), padded.begin());
        return padded;
    }

    /**
     * @brief Compute 2D FFT
     *
     * @param image 2D signal (matrix)
     * @return 2D frequency domain
     */
    static ComplexMatrix fft2D(const Matrix& image) {
        size_t rows = image.size();
        size_t cols = image[0].size();

        // Convert to complex
        ComplexMatrix complex_image(rows, ComplexSignal(cols));
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                complex_image[i][j] = Complex(image[i][j], 0.0);
            }
        }

        // FFT along rows
        for (size_t i = 0; i < rows; ++i) {
            complex_image[i] = fft(complex_image[i]);
        }

        // FFT along columns
        for (size_t j = 0; j < cols; ++j) {
            ComplexSignal column(rows);
            for (size_t i = 0; i < rows; ++i) {
                column[i] = complex_image[i][j];
            }
            column = fft(column);
            for (size_t i = 0; i < rows; ++i) {
                complex_image[i][j] = column[i];
            }
        }

        return complex_image;
    }
};

/**
 * @class CirculantMatrix
 * @brief Circulant matrices and convolution operators
 *
 * A circulant matrix is diagonalized by the DFT matrix
 */
class CirculantMatrix {
public:
    /**
     * @brief Create circulant matrix from first column
     *
     * C[i,j] = c[(i-j) mod N]
     *
     * @param first_column First column of circulant matrix
     * @return Circulant matrix
     */
    static Matrix createCirculant(const Signal& first_column) {
        size_t N = first_column.size();
        Matrix C(N, Signal(N));

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                C[i][j] = first_column[(i - j + N) % N];
            }
        }

        return C;
    }

    /**
     * @brief Compute eigenvalues of circulant matrix
     *
     * Eigenvalues are the DFT of the first column
     *
     * @param first_column First column
     * @return Eigenvalues (DFT of first column)
     */
    static ComplexSignal eigenvalues(const Signal& first_column) {
        return DiscreteFourierTransform::dft(first_column);
    }

    /**
     * @brief Multiply circulant matrix by vector
     *
     * Efficient using FFT: C*x = IFFT(FFT(c) .* FFT(x))
     *
     * @param first_column First column of C
     * @param x Vector to multiply
     * @return C*x
     */
    static Signal multiply(const Signal& first_column, const Signal& x) {
        if (first_column.size() != x.size()) {
            throw std::invalid_argument("Size mismatch");
        }

        // Pad to power of 2
        Signal c_padded = FastFourierTransform::zeroPad(first_column);
        Signal x_padded = FastFourierTransform::zeroPad(x);

        // FFT
        ComplexSignal fft_c = FastFourierTransform::fft(c_padded);
        ComplexSignal fft_x = FastFourierTransform::fft(x_padded);

        // Element-wise multiplication
        ComplexSignal fft_result(fft_c.size());
        for (size_t i = 0; i < fft_c.size(); ++i) {
            fft_result[i] = fft_c[i] * fft_x[i];
        }

        // IFFT
        ComplexSignal result_complex = FastFourierTransform::ifft(fft_result);

        // Extract real part and truncate
        Signal result(first_column.size());
        for (size_t i = 0; i < result.size(); ++i) {
            result[i] = result_complex[i].real();
        }

        return result;
    }
};

/**
 * @class Convolution
 * @brief Convolution operators and fast convolution
 */
class Convolution {
public:
    /**
     * @brief Compute circular convolution (direct)
     *
     * (f * g)[n] = Σ_m f[m] g[(n-m) mod N]
     *
     * @param f First signal
     * @param g Second signal
     * @return Circular convolution f * g
     */
    static Signal circular(const Signal& f, const Signal& g) {
        if (f.size() != g.size()) {
            throw std::invalid_argument("Signals must have same length");
        }

        size_t N = f.size();
        Signal result(N, 0.0);

        for (size_t n = 0; n < N; ++n) {
            for (size_t m = 0; m < N; ++m) {
                result[n] += f[m] * g[(n - m + N) % N];
            }
        }

        return result;
    }

    /**
     * @brief Fast circular convolution using FFT
     *
     * Convolution theorem: f * g = IFFT(FFT(f) · FFT(g))
     *
     * @param f First signal
     * @param g Second signal
     * @return Circular convolution f * g
     */
    static Signal fastCircular(const Signal& f, const Signal& g) {
        if (f.size() != g.size()) {
            throw std::invalid_argument("Signals must have same length");
        }

        // Pad to power of 2
        Signal f_padded = FastFourierTransform::zeroPad(f);
        Signal g_padded = FastFourierTransform::zeroPad(g);

        // FFT
        ComplexSignal fft_f = FastFourierTransform::fft(f_padded);
        ComplexSignal fft_g = FastFourierTransform::fft(g_padded);

        // Pointwise multiplication
        ComplexSignal fft_result(fft_f.size());
        for (size_t i = 0; i < fft_f.size(); ++i) {
            fft_result[i] = fft_f[i] * fft_g[i];
        }

        // IFFT
        ComplexSignal result_complex = FastFourierTransform::ifft(fft_result);

        // Extract real part
        Signal result(f.size());
        for (size_t i = 0; i < result.size(); ++i) {
            result[i] = result_complex[i].real();
        }

        return result;
    }

    /**
     * @brief Linear convolution
     *
     * @param f First signal
     * @param g Second signal (filter)
     * @return Linear convolution
     */
    static Signal linear(const Signal& f, const Signal& g) {
        size_t N = f.size();
        size_t M = g.size();
        Signal result(N + M - 1, 0.0);

        for (size_t n = 0; n < result.size(); ++n) {
            for (size_t m = 0; m < M; ++m) {
                if (n >= m && n - m < N) {
                    result[n] += f[n - m] * g[m];
                }
            }
        }

        return result;
    }

    /**
     * @brief Cross-correlation
     *
     * (f ⋆ g)[n] = Σ_m f[m] g[m+n]
     *
     * @param f First signal
     * @param g Second signal
     * @return Cross-correlation
     */
    static Signal crossCorrelation(const Signal& f, const Signal& g) {
        // Correlation = convolution with time-reversed signal
        Signal g_reversed = g;
        std::reverse(g_reversed.begin(), g_reversed.end());
        return fastCircular(f, g_reversed);
    }

    /**
     * @brief Auto-correlation
     *
     * @param f Signal
     * @return Auto-correlation f ⋆ f
     */
    static Signal autoCorrelation(const Signal& f) {
        return crossCorrelation(f, f);
    }
};

/**
 * @class WaveletTransform
 * @brief Wavelet transforms: Haar and Daubechies
 */
class WaveletTransform {
public:
    /**
     * @brief Haar wavelet transform (simplest wavelet)
     *
     * @param signal Input signal (length must be power of 2)
     * @return Wavelet coefficients
     */
    static Signal haarTransform(const Signal& signal) {
        size_t N = signal.size();
        if (N == 0 || (N & (N - 1)) != 0) {
            throw std::invalid_argument("Size must be power of 2");
        }

        Signal coeffs = signal;
        Signal temp(N);

        for (size_t len = N; len >= 2; len /= 2) {
            for (size_t i = 0; i < len / 2; ++i) {
                // Averaging (low-pass)
                temp[i] = (coeffs[2 * i] + coeffs[2 * i + 1]) / std::sqrt(2.0);
                // Differencing (high-pass)
                temp[len / 2 + i] = (coeffs[2 * i] - coeffs[2 * i + 1]) / std::sqrt(2.0);
            }
            std::copy(temp.begin(), temp.begin() + len, coeffs.begin());
        }

        return coeffs;
    }

    /**
     * @brief Inverse Haar wavelet transform
     *
     * @param coeffs Wavelet coefficients
     * @return Reconstructed signal
     */
    static Signal haarInverse(const Signal& coeffs) {
        size_t N = coeffs.size();
        Signal signal = coeffs;
        Signal temp(N);

        for (size_t len = 2; len <= N; len *= 2) {
            for (size_t i = 0; i < len / 2; ++i) {
                // Reconstruct from averages and differences
                temp[2 * i] = (signal[i] + signal[len / 2 + i]) / std::sqrt(2.0);
                temp[2 * i + 1] = (signal[i] - signal[len / 2 + i]) / std::sqrt(2.0);
            }
            std::copy(temp.begin(), temp.begin() + len, signal.begin());
        }

        return signal;
    }

    /**
     * @brief Daubechies-4 wavelet transform
     *
     * Uses Daubechies D4 scaling and wavelet filters
     *
     * @param signal Input signal (length must be power of 2)
     * @return Wavelet coefficients
     */
    static Signal daubechies4Transform(const Signal& signal) {
        size_t N = signal.size();
        if (N < 4 || (N & (N - 1)) != 0) {
            throw std::invalid_argument("Size must be power of 2 and >= 4");
        }

        // Daubechies D4 coefficients
        const double sqrt3 = std::sqrt(3.0);
        const double c0 = (1.0 + sqrt3) / (4.0 * std::sqrt(2.0));
        const double c1 = (3.0 + sqrt3) / (4.0 * std::sqrt(2.0));
        const double c2 = (3.0 - sqrt3) / (4.0 * std::sqrt(2.0));
        const double c3 = (1.0 - sqrt3) / (4.0 * std::sqrt(2.0));

        Signal coeffs = signal;
        Signal temp(N);

        for (size_t len = N; len >= 4; len /= 2) {
            for (size_t i = 0; i < len / 2; ++i) {
                size_t i0 = 2 * i;
                size_t i1 = (2 * i + 1) % len;
                size_t i2 = (2 * i + 2) % len;
                size_t i3 = (2 * i + 3) % len;

                // Low-pass (scaling)
                temp[i] = c0 * coeffs[i0] + c1 * coeffs[i1] +
                         c2 * coeffs[i2] + c3 * coeffs[i3];

                // High-pass (wavelet)
                temp[len / 2 + i] = c3 * coeffs[i0] - c2 * coeffs[i1] +
                                   c1 * coeffs[i2] - c0 * coeffs[i3];
            }
            std::copy(temp.begin(), temp.begin() + len, coeffs.begin());
        }

        return coeffs;
    }

    /**
     * @brief Compute wavelet packet decomposition
     *
     * @param signal Input signal
     * @param level Decomposition level
     * @return Vector of coefficient vectors at each level
     */
    static std::vector<Signal> packetDecomposition(const Signal& signal, int level) {
        std::vector<Signal> packets;
        Signal current = signal;

        for (int l = 0; l < level; ++l) {
            current = haarTransform(current);
            packets.push_back(current);
        }

        return packets;
    }
};

/**
 * @class FourierSeries
 * @brief Fourier series for periodic functions on S¹
 */
class FourierSeries {
public:
    /**
     * @brief Compute Fourier series coefficients
     *
     * c_k = (1/N) Σ_n f(2πn/N) e^{-2πikn/N}
     *
     * @param samples Function samples on [0, 2π)
     * @return Fourier coefficients
     */
    static ComplexSignal coefficients(const Signal& samples) {
        return DiscreteFourierTransform::dft(samples);
    }

    /**
     * @brief Reconstruct function from Fourier coefficients
     *
     * f(x) ≈ Σ_k c_k e^{ikx}
     *
     * @param coeffs Fourier coefficients
     * @param n_points Number of evaluation points
     * @return Reconstructed function values
     */
    static Signal reconstruct(const ComplexSignal& coeffs, size_t n_points) {
        ComplexSignal padded_coeffs(n_points);
        size_t N = std::min(coeffs.size(), n_points);

        for (size_t i = 0; i < N; ++i) {
            padded_coeffs[i] = coeffs[i];
        }

        ComplexSignal reconstructed = DiscreteFourierTransform::idft(padded_coeffs);

        Signal result(n_points);
        for (size_t i = 0; i < n_points; ++i) {
            result[i] = reconstructed[i].real();
        }

        return result;
    }

    /**
     * @brief Apply Fourier multiplier
     *
     * T_m f(x) = Σ_k m(k) ĉ_k e^{ikx}
     *
     * @param samples Function samples
     * @param multiplier Multiplier function m(k)
     * @return Transformed function
     */
    static Signal applyMultiplier(const Signal& samples,
                                  std::function<Complex(int)> multiplier) {
        ComplexSignal coeffs = DiscreteFourierTransform::dft(samples);

        // Apply multiplier
        for (size_t k = 0; k < coeffs.size(); ++k) {
            int k_signed = (k <= coeffs.size() / 2) ? k : k - coeffs.size();
            coeffs[k] *= multiplier(k_signed);
        }

        // Reconstruct
        ComplexSignal result_complex = DiscreteFourierTransform::idft(coeffs);

        Signal result(samples.size());
        for (size_t i = 0; i < samples.size(); ++i) {
            result[i] = result_complex[i].real();
        }

        return result;
    }

    /**
     * @brief Compute derivative using Fourier multiplier
     *
     * D_x f = Σ_k (ik) ĉ_k e^{ikx}
     *
     * @param samples Function samples
     * @return Derivative samples
     */
    static Signal derivative(const Signal& samples) {
        return applyMultiplier(samples, [](int k) { return Complex(0.0, k); });
    }

    /**
     * @brief Apply fractional Laplacian (-Δ)^s
     *
     * @param samples Function samples
     * @param s Fractional power
     * @return Result of (-Δ)^s f
     */
    static Signal fractionalLaplacian(const Signal& samples, double s) {
        return applyMultiplier(samples, [s](int k) {
            return std::pow(std::abs(k), 2.0 * s);
        });
    }
};

/**
 * @class HilbertSpaceOperators
 * @brief Operators on discrete L² spaces
 */
class HilbertSpaceOperators {
public:
    /**
     * @brief Compute inner product ⟨f, g⟩
     *
     * @param f First signal
     * @param g Second signal
     * @return Inner product
     */
    static double innerProduct(const Signal& f, const Signal& g) {
        if (f.size() != g.size()) {
            throw std::invalid_argument("Size mismatch");
        }

        double sum = 0.0;
        for (size_t i = 0; i < f.size(); ++i) {
            sum += f[i] * g[i];
        }
        return sum;
    }

    /**
     * @brief Compute L² norm ‖f‖₂
     *
     * @param f Signal
     * @return L² norm
     */
    static double l2Norm(const Signal& f) {
        return std::sqrt(innerProduct(f, f));
    }

    /**
     * @brief Apply self-adjoint operator (symmetric matrix)
     *
     * @param A Symmetric matrix
     * @param x Vector
     * @return A*x
     */
    static Signal applySelfAdjoint(const Matrix& A, const Signal& x) {
        size_t N = A.size();
        Signal result(N, 0.0);

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                result[i] += A[i][j] * x[j];
            }
        }

        return result;
    }

    /**
     * @brief Compute trace of operator
     *
     * tr(A) = Σ_i A_ii
     *
     * @param A Matrix
     * @return Trace
     */
    static double trace(const Matrix& A) {
        double tr = 0.0;
        for (size_t i = 0; i < A.size(); ++i) {
            tr += A[i][i];
        }
        return tr;
    }

    /**
     * @brief Compute Frobenius norm (Schatten 2-norm)
     *
     * ‖A‖_F = √(tr(A*A))
     *
     * @param A Matrix
     * @return Frobenius norm
     */
    static double frobeniusNorm(const Matrix& A) {
        double sum = 0.0;
        for (const auto& row : A) {
            for (double val : row) {
                sum += val * val;
            }
        }
        return std::sqrt(sum);
    }

    /**
     * @brief Compute Schatten p-norm
     *
     * ‖A‖_p = (Σ σ_i^p)^{1/p} where σ_i are singular values
     * Approximate using eigenvalues for symmetric matrices
     *
     * @param A Matrix
     * @param p Schatten parameter
     * @return Schatten p-norm (approximation)
     */
    static double schattenNorm(const Matrix& A, double p) {
        // For diagonal approximation
        double sum = 0.0;
        for (size_t i = 0; i < A.size(); ++i) {
            sum += std::pow(std::abs(A[i][i]), p);
        }
        return std::pow(sum, 1.0 / p);
    }

    /**
     * @brief Check if operator is compact (finite rank approximation)
     *
     * @param A Matrix
     * @param rank Rank threshold
     * @param tolerance Singular value tolerance
     * @return true if effectively finite rank
     */
    static bool isCompact(const Matrix& A, size_t rank, double tolerance = 1e-10) {
        // Check if small eigenvalues below tolerance
        size_t small_values = 0;
        for (size_t i = 0; i < A.size(); ++i) {
            if (std::abs(A[i][i]) < tolerance) {
                small_values++;
            }
        }
        return (A.size() - small_values) <= rank;
    }
};

/**
 * @class TimeFrequencyAnalysis
 * @brief Short-time Fourier transform and spectrograms
 */
class TimeFrequencyAnalysis {
public:
    /**
     * @brief Compute short-time Fourier transform (STFT)
     *
     * @param signal Input signal
     * @param window_size Window size
     * @param hop_size Hop size between windows
     * @return STFT matrix (time × frequency)
     */
    static ComplexMatrix stft(const Signal& signal,
                             size_t window_size,
                             size_t hop_size) {
        size_t N = signal.size();
        size_t n_frames = (N - window_size) / hop_size + 1;
        ComplexMatrix stft_matrix(n_frames);

        // Hann window
        Signal window = hannWindow(window_size);

        for (size_t frame = 0; frame < n_frames; ++frame) {
            size_t start = frame * hop_size;
            Signal windowed(window_size);

            for (size_t i = 0; i < window_size; ++i) {
                windowed[i] = signal[start + i] * window[i];
            }

            stft_matrix[frame] = DiscreteFourierTransform::dft(windowed);
        }

        return stft_matrix;
    }

    /**
     * @brief Compute spectrogram |STFT|²
     *
     * @param signal Input signal
     * @param window_size Window size
     * @param hop_size Hop size
     * @return Spectrogram (time × frequency)
     */
    static Matrix spectrogram(const Signal& signal,
                             size_t window_size,
                             size_t hop_size) {
        ComplexMatrix stft_result = stft(signal, window_size, hop_size);
        Matrix spec(stft_result.size());

        for (size_t t = 0; t < stft_result.size(); ++t) {
            spec[t].resize(stft_result[t].size());
            for (size_t f = 0; f < stft_result[t].size(); ++f) {
                spec[t][f] = std::norm(stft_result[t][f]);
            }
        }

        return spec;
    }

    /**
     * @brief Hann window
     *
     * @param size Window size
     * @return Hann window coefficients
     */
    static Signal hannWindow(size_t size) {
        Signal window(size);
        for (size_t i = 0; i < size; ++i) {
            window[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (size - 1)));
        }
        return window;
    }

    /**
     * @brief Hamming window
     *
     * @param size Window size
     * @return Hamming window coefficients
     */
    static Signal hammingWindow(size_t size) {
        Signal window(size);
        for (size_t i = 0; i < size; ++i) {
            window[i] = 0.54 - 0.46 * std::cos(2.0 * M_PI * i / (size - 1));
        }
        return window;
    }

    /**
     * @brief Gabor transform (STFT with Gaussian window)
     *
     * @param signal Input signal
     * @param window_size Window size
     * @param hop_size Hop size
     * @param sigma Gaussian width
     * @return Gabor transform
     */
    static ComplexMatrix gaborTransform(const Signal& signal,
                                       size_t window_size,
                                       size_t hop_size,
                                       double sigma) {
        // Create Gaussian window
        Signal window(window_size);
        double center = window_size / 2.0;

        for (size_t i = 0; i < window_size; ++i) {
            double x = i - center;
            window[i] = std::exp(-x * x / (2.0 * sigma * sigma));
        }

        // Normalize
        double sum = 0.0;
        for (double w : window) sum += w;
        for (double& w : window) w /= sum;

        // Apply STFT with Gaussian window
        size_t N = signal.size();
        size_t n_frames = (N - window_size) / hop_size + 1;
        ComplexMatrix result(n_frames);

        for (size_t frame = 0; frame < n_frames; ++frame) {
            size_t start = frame * hop_size;
            Signal windowed(window_size);

            for (size_t i = 0; i < window_size; ++i) {
                windowed[i] = signal[start + i] * window[i];
            }

            result[frame] = DiscreteFourierTransform::dft(windowed);
        }

        return result;
    }
};

} // namespace maths::analysis

#endif // MATHS_ANALYSIS_FOURIER_ANALYSIS_HPP
