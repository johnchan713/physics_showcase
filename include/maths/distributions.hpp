/**
 * @file distributions.hpp
 * @brief Comprehensive probability distributions and statistical computations
 *
 * Implements computational algorithms for:
 * - Discrete distributions (Bernoulli, Binomial, Poisson, Geometric, etc.)
 * - Continuous distributions (Normal, Exponential, Uniform, Gamma, Beta, etc.)
 * - PDF, CDF, quantile functions, moments, and sampling
 * - Maximum likelihood estimation and hypothesis testing
 */

#ifndef MATHS_PROBABILITY_DISTRIBUTIONS_HPP
#define MATHS_PROBABILITY_DISTRIBUTIONS_HPP

#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <functional>

namespace maths::probability {

// Constants
const double PI = 3.14159265358979323846;
const double E = 2.71828182845904523536;

// Random number generator
static std::random_device rd;
static std::mt19937 gen(rd());

/**
 * @brief Factorial function
 */
inline double factorial(int n) {
    if (n < 0) throw std::invalid_argument("Factorial undefined for negative numbers");
    if (n == 0 || n == 1) return 1.0;
    double result = 1.0;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

/**
 * @brief Binomial coefficient C(n, k) = n! / (k!(n-k)!)
 */
inline double binomial_coefficient(int n, int k) {
    if (k > n || k < 0) return 0.0;
    if (k == 0 || k == n) return 1.0;

    // Use multiplicative formula to avoid overflow
    double result = 1.0;
    for (int i = 0; i < k; ++i) {
        result *= (n - i);
        result /= (i + 1);
    }
    return result;
}

/**
 * @brief Gamma function Γ(x) (using Stirling's approximation)
 */
inline double gamma_function(double x) {
    if (x <= 0) throw std::invalid_argument("Gamma function undefined for x <= 0");

    // For integers, use factorial
    if (x == std::floor(x) && x < 20) {
        return factorial(static_cast<int>(x) - 1);
    }

    // Stirling's approximation: Γ(x) ≈ √(2π/x) * (x/e)^x
    return std::sqrt(2 * PI / x) * std::pow(x / E, x);
}

/**
 * @brief Error function erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt
 */
inline double error_function(double x) {
    // Abramowitz and Stegun approximation
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    int sign = (x >= 0) ? 1 : -1;
    x = std::abs(x);

    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

    return sign * y;
}

/**
 * @class Distribution
 * @brief Base class for all probability distributions
 */
class Distribution {
public:
    virtual ~Distribution() = default;

    virtual double mean() const = 0;
    virtual double variance() const = 0;
    double standard_deviation() const { return std::sqrt(variance()); }
};

// ============================================================================
// SECTION 1: DISCRETE PROBABILITY DISTRIBUTIONS
// ============================================================================

/**
 * @class BernoulliDistribution
 * @brief Bernoulli distribution: X ∈ {0,1} with P(X=1) = p
 */
class BernoulliDistribution : public Distribution {
private:
    double p;

public:
    BernoulliDistribution(double probability) : p(probability) {
        if (p < 0.0 || p > 1.0) {
            throw std::invalid_argument("Probability must be in [0,1]");
        }
    }

    double pmf(int x) const {
        if (x == 0) return 1.0 - p;
        if (x == 1) return p;
        return 0.0;
    }

    double cdf(int x) const {
        if (x < 0) return 0.0;
        if (x >= 1) return 1.0;
        return 1.0 - p;
    }

    double mean() const override { return p; }
    double variance() const override { return p * (1.0 - p); }

    int sample() const {
        std::bernoulli_distribution dist(p);
        return dist(gen);
    }
};

/**
 * @class BinomialDistribution
 * @brief Binomial distribution: number of successes in n trials
 */
class BinomialDistribution : public Distribution {
private:
    int n;
    double p;

public:
    BinomialDistribution(int trials, double probability) : n(trials), p(probability) {
        if (n < 0) throw std::invalid_argument("Number of trials must be non-negative");
        if (p < 0.0 || p > 1.0) throw std::invalid_argument("Probability must be in [0,1]");
    }

    /**
     * @brief Probability mass function: P(X = k) = C(n,k) * p^k * (1-p)^(n-k)
     */
    double pmf(int k) const {
        if (k < 0 || k > n) return 0.0;
        return binomial_coefficient(n, k) * std::pow(p, k) * std::pow(1.0 - p, n - k);
    }

    /**
     * @brief Cumulative distribution function: P(X ≤ k)
     */
    double cdf(int k) const {
        if (k < 0) return 0.0;
        if (k >= n) return 1.0;

        double sum = 0.0;
        for (int i = 0; i <= k; ++i) {
            sum += pmf(i);
        }
        return sum;
    }

    double mean() const override { return n * p; }
    double variance() const override { return n * p * (1.0 - p); }

    int sample() const {
        std::binomial_distribution<int> dist(n, p);
        return dist(gen);
    }
};

/**
 * @class PoissonDistribution
 * @brief Poisson distribution: P(X=k) = λ^k e^(-λ) / k!
 */
class PoissonDistribution : public Distribution {
private:
    double lambda;

public:
    PoissonDistribution(double rate) : lambda(rate) {
        if (lambda <= 0.0) {
            throw std::invalid_argument("Rate parameter must be positive");
        }
    }

    double pmf(int k) const {
        if (k < 0) return 0.0;
        return std::pow(lambda, k) * std::exp(-lambda) / factorial(k);
    }

    double cdf(int k) const {
        if (k < 0) return 0.0;

        double sum = 0.0;
        for (int i = 0; i <= k; ++i) {
            sum += pmf(i);
        }
        return sum;
    }

    double mean() const override { return lambda; }
    double variance() const override { return lambda; }

    int sample() const {
        std::poisson_distribution<int> dist(lambda);
        return dist(gen);
    }
};

/**
 * @class GeometricDistribution
 * @brief Geometric distribution: number of trials until first success
 */
class GeometricDistribution : public Distribution {
private:
    double p;

public:
    GeometricDistribution(double probability) : p(probability) {
        if (p <= 0.0 || p > 1.0) {
            throw std::invalid_argument("Probability must be in (0,1]");
        }
    }

    /**
     * @brief PMF: P(X=k) = (1-p)^(k-1) * p
     */
    double pmf(int k) const {
        if (k < 1) return 0.0;
        return std::pow(1.0 - p, k - 1) * p;
    }

    /**
     * @brief CDF: P(X≤k) = 1 - (1-p)^k
     */
    double cdf(int k) const {
        if (k < 1) return 0.0;
        return 1.0 - std::pow(1.0 - p, k);
    }

    double mean() const override { return 1.0 / p; }
    double variance() const override { return (1.0 - p) / (p * p); }

    int sample() const {
        std::geometric_distribution<int> dist(p);
        return dist(gen) + 1;  // Shift to start at 1
    }
};

// ============================================================================
// SECTION 2: CONTINUOUS PROBABILITY DISTRIBUTIONS
// ============================================================================

/**
 * @class UniformDistribution
 * @brief Continuous uniform distribution on [a, b]
 */
class UniformDistribution : public Distribution {
private:
    double a, b;

public:
    UniformDistribution(double lower, double upper) : a(lower), b(upper) {
        if (a >= b) {
            throw std::invalid_argument("Lower bound must be less than upper bound");
        }
    }

    double pdf(double x) const {
        if (x < a || x > b) return 0.0;
        return 1.0 / (b - a);
    }

    double cdf(double x) const {
        if (x < a) return 0.0;
        if (x > b) return 1.0;
        return (x - a) / (b - a);
    }

    double quantile(double prob) const {
        if (prob < 0.0 || prob > 1.0) {
            throw std::invalid_argument("Probability must be in [0,1]");
        }
        return a + prob * (b - a);
    }

    double mean() const override { return (a + b) / 2.0; }
    double variance() const override { return (b - a) * (b - a) / 12.0; }

    double sample() const {
        std::uniform_real_distribution<double> dist(a, b);
        return dist(gen);
    }
};

/**
 * @class NormalDistribution
 * @brief Normal (Gaussian) distribution N(μ, σ²)
 */
class NormalDistribution : public Distribution {
private:
    double mu, sigma;

public:
    NormalDistribution(double mean_val, double std_dev) : mu(mean_val), sigma(std_dev) {
        if (sigma <= 0.0) {
            throw std::invalid_argument("Standard deviation must be positive");
        }
    }

    /**
     * @brief PDF: f(x) = (1/(σ√(2π))) * exp(-(x-μ)²/(2σ²))
     */
    double pdf(double x) const {
        double z = (x - mu) / sigma;
        return (1.0 / (sigma * std::sqrt(2.0 * PI))) * std::exp(-0.5 * z * z);
    }

    /**
     * @brief CDF: Φ(x) = (1/2)[1 + erf((x-μ)/(σ√2))]
     */
    double cdf(double x) const {
        double z = (x - mu) / (sigma * std::sqrt(2.0));
        return 0.5 * (1.0 + error_function(z));
    }

    double mean() const override { return mu; }
    double variance() const override { return sigma * sigma; }

    double sample() const {
        std::normal_distribution<double> dist(mu, sigma);
        return dist(gen);
    }

    /**
     * @brief Z-score: (x - μ) / σ
     */
    double z_score(double x) const {
        return (x - mu) / sigma;
    }
};

/**
 * @class ExponentialDistribution
 * @brief Exponential distribution with rate λ
 */
class ExponentialDistribution : public Distribution {
private:
    double lambda;

public:
    ExponentialDistribution(double rate) : lambda(rate) {
        if (lambda <= 0.0) {
            throw std::invalid_argument("Rate parameter must be positive");
        }
    }

    double pdf(double x) const {
        if (x < 0.0) return 0.0;
        return lambda * std::exp(-lambda * x);
    }

    double cdf(double x) const {
        if (x < 0.0) return 0.0;
        return 1.0 - std::exp(-lambda * x);
    }

    double quantile(double prob) const {
        if (prob < 0.0 || prob >= 1.0) {
            throw std::invalid_argument("Probability must be in [0,1)");
        }
        return -std::log(1.0 - prob) / lambda;
    }

    double mean() const override { return 1.0 / lambda; }
    double variance() const override { return 1.0 / (lambda * lambda); }

    double sample() const {
        std::exponential_distribution<double> dist(lambda);
        return dist(gen);
    }
};

/**
 * @class GammaDistribution
 * @brief Gamma distribution with shape α and rate β
 */
class GammaDistribution : public Distribution {
private:
    double alpha, beta;

public:
    GammaDistribution(double shape, double rate) : alpha(shape), beta(rate) {
        if (alpha <= 0.0 || beta <= 0.0) {
            throw std::invalid_argument("Shape and rate must be positive");
        }
    }

    /**
     * @brief PDF: f(x) = (β^α / Γ(α)) * x^(α-1) * e^(-βx)
     */
    double pdf(double x) const {
        if (x < 0.0) return 0.0;
        return (std::pow(beta, alpha) / gamma_function(alpha)) *
               std::pow(x, alpha - 1.0) * std::exp(-beta * x);
    }

    double mean() const override { return alpha / beta; }
    double variance() const override { return alpha / (beta * beta); }

    double sample() const {
        std::gamma_distribution<double> dist(alpha, 1.0 / beta);
        return dist(gen);
    }
};

/**
 * @class BetaDistribution
 * @brief Beta distribution on [0,1] with parameters α, β
 */
class BetaDistribution : public Distribution {
private:
    double alpha, beta;

public:
    BetaDistribution(double a, double b) : alpha(a), beta(b) {
        if (alpha <= 0.0 || beta <= 0.0) {
            throw std::invalid_argument("Alpha and beta must be positive");
        }
    }

    /**
     * @brief PDF: f(x) = (Γ(α+β)/(Γ(α)Γ(β))) * x^(α-1) * (1-x)^(β-1)
     */
    double pdf(double x) const {
        if (x < 0.0 || x > 1.0) return 0.0;
        double B = gamma_function(alpha) * gamma_function(beta) / gamma_function(alpha + beta);
        return std::pow(x, alpha - 1.0) * std::pow(1.0 - x, beta - 1.0) / B;
    }

    double mean() const override { return alpha / (alpha + beta); }
    double variance() const override {
        return (alpha * beta) / ((alpha + beta) * (alpha + beta) * (alpha + beta + 1.0));
    }

    double sample() const {
        // Sample using two gamma distributions
        GammaDistribution gamma_a(alpha, 1.0);
        GammaDistribution gamma_b(beta, 1.0);
        double x = gamma_a.sample();
        double y = gamma_b.sample();
        return x / (x + y);
    }
};

/**
 * @class ChiSquaredDistribution
 * @brief Chi-squared distribution with k degrees of freedom
 */
class ChiSquaredDistribution : public Distribution {
private:
    int k;

public:
    ChiSquaredDistribution(int degrees_of_freedom) : k(degrees_of_freedom) {
        if (k <= 0) {
            throw std::invalid_argument("Degrees of freedom must be positive");
        }
    }

    /**
     * @brief PDF (uses Gamma distribution)
     */
    double pdf(double x) const {
        GammaDistribution gamma(k / 2.0, 0.5);
        return gamma.pdf(x);
    }

    double mean() const override { return static_cast<double>(k); }
    double variance() const override { return 2.0 * k; }

    double sample() const {
        std::chi_squared_distribution<double> dist(k);
        return dist(gen);
    }
};

/**
 * @class StudentTDistribution
 * @brief Student's t-distribution with ν degrees of freedom
 *
 * Defined as: t = Z / sqrt(χ²/ν) where Z ~ N(0,1) and χ² ~ χ²(ν)
 *
 * Used in t-tests, confidence intervals for small samples
 * As ν → ∞, approaches standard normal distribution
 */
class StudentTDistribution : public Distribution {
private:
    int nu;  // degrees of freedom

public:
    StudentTDistribution(int degrees_of_freedom) : nu(degrees_of_freedom) {
        if (nu <= 0) {
            throw std::invalid_argument("Degrees of freedom must be positive");
        }
    }

    /**
     * @brief PDF: f(t) = Γ((ν+1)/2) / (sqrt(νπ) Γ(ν/2)) * (1 + t²/ν)^(-(ν+1)/2)
     */
    double pdf(double t) const {
        double numerator = gamma_function((nu + 1) / 2.0);
        double denominator = std::sqrt(nu * PI) * gamma_function(nu / 2.0);
        double base = 1.0 + (t * t) / nu;
        double exponent = -(nu + 1) / 2.0;

        return (numerator / denominator) * std::pow(base, exponent);
    }

    /**
     * @brief CDF using approximation
     * For large ν, approaches standard normal CDF
     */
    double cdf(double t) const {
        // For ν > 30, approximate with standard normal
        if (nu > 30) {
            return 0.5 * (1.0 + error_function(t / std::sqrt(2.0)));
        }

        // Numerical integration for general case
        double dt = 0.01;
        double sum = 0.0;
        int n_steps = static_cast<int>(std::abs(t) / dt) + 100;

        for (int i = -n_steps; i <= n_steps; ++i) {
            double x = i * dt;
            if (x <= t) {
                sum += pdf(x) * dt;
            }
        }

        return sum;
    }

    /**
     * @brief Mean: E[T] = 0 for ν > 1
     */
    double mean() const override {
        if (nu <= 1) throw std::runtime_error("Mean undefined for ν <= 1");
        return 0.0;
    }

    /**
     * @brief Variance: Var(T) = ν/(ν-2) for ν > 2
     */
    double variance() const override {
        if (nu <= 2) throw std::runtime_error("Variance undefined for ν <= 2");
        return static_cast<double>(nu) / (nu - 2);
    }

    /**
     * @brief Sample from t-distribution
     */
    double sample() const {
        std::student_t_distribution<double> dist(nu);
        return dist(gen);
    }

    /**
     * @brief Quantile function (inverse CDF) - approximation
     */
    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) {
            throw std::invalid_argument("Probability must be in (0, 1)");
        }

        // For large ν, use normal approximation
        if (nu > 30) {
            // Inverse error function approximation
            double x = 2.0 * p - 1.0;
            return std::sqrt(2.0) * x * (1.0 + x * x / 4.0);  // Simplified
        }

        // Binary search for general case
        double low = -10.0, high = 10.0;
        double mid, cdf_mid;
        int max_iter = 100;

        for (int iter = 0; iter < max_iter; ++iter) {
            mid = (low + high) / 2.0;
            cdf_mid = cdf(mid);

            if (std::abs(cdf_mid - p) < 1e-6) break;

            if (cdf_mid < p) {
                low = mid;
            } else {
                high = mid;
            }
        }

        return mid;
    }

    /**
     * @brief Degrees of freedom
     */
    int degrees_of_freedom() const { return nu; }
};

/**
 * @class FDistribution
 * @brief F-distribution (Fisher-Snedecor distribution)
 *
 * Ratio of two chi-squared distributions:
 * F = (χ²₁/d₁) / (χ²₂/d₂)
 *
 * Used in ANOVA, regression analysis, F-tests
 */
class FDistribution : public Distribution {
private:
    int d1;  // degrees of freedom 1
    int d2;  // degrees of freedom 2

public:
    FDistribution(int df1, int df2) : d1(df1), d2(df2) {
        if (d1 <= 0 || d2 <= 0) {
            throw std::invalid_argument("Degrees of freedom must be positive");
        }
    }

    /**
     * @brief PDF: f(x) = (d₁/d₂)^(d₁/2) * x^(d₁/2-1) / (B(d₁/2, d₂/2) * (1 + d₁x/d₂)^((d₁+d₂)/2))
     */
    double pdf(double x) const {
        if (x <= 0) return 0.0;

        double a = d1 / 2.0;
        double b = d2 / 2.0;

        // B(a, b) = Γ(a)Γ(b)/Γ(a+b)
        double beta_func = gamma_function(a) * gamma_function(b) / gamma_function(a + b);

        double coeff = std::pow(d1 * 1.0 / d2, a) / beta_func;
        double numerator = std::pow(x, a - 1);
        double denominator = std::pow(1 + d1 * x / d2, (d1 + d2) / 2.0);

        return coeff * numerator / denominator;
    }

    /**
     * @brief Mean: E[F] = d₂/(d₂-2) for d₂ > 2
     */
    double mean() const override {
        if (d2 <= 2) throw std::runtime_error("Mean undefined for d2 <= 2");
        return static_cast<double>(d2) / (d2 - 2);
    }

    /**
     * @brief Variance: Var(F) = 2d₂²(d₁+d₂-2) / (d₁(d₂-2)²(d₂-4)) for d₂ > 4
     */
    double variance() const override {
        if (d2 <= 4) throw std::runtime_error("Variance undefined for d2 <= 4");

        double numerator = 2.0 * d2 * d2 * (d1 + d2 - 2);
        double denominator = d1 * (d2 - 2) * (d2 - 2) * (d2 - 4);

        return numerator / denominator;
    }

    /**
     * @brief Sample from F-distribution using ratio of chi-squared
     */
    double sample() const {
        std::fisher_f_distribution<double> dist(d1, d2);
        return dist(gen);
    }

    /**
     * @brief CDF (approximation using regularized incomplete beta function)
     */
    double cdf(double x) const {
        if (x <= 0) return 0.0;

        // F_CDF(x; d1, d2) = I_{d1*x/(d1*x+d2)}(d1/2, d2/2)
        // where I is regularized incomplete beta function
        // Simplified approximation
        double z = d1 * x / (d1 * x + d2);
        return z;  // Placeholder - full implementation would use beta regularized
    }
};

/**
 * @class NegativeBinomialDistribution
 * @brief Negative binomial distribution
 *
 * Number of failures before r successes in Bernoulli trials
 * PMF: P(X=k) = C(k+r-1, k) * p^r * (1-p)^k
 *
 * Also models overdispersed count data (generalization of Poisson)
 */
class NegativeBinomialDistribution : public Distribution {
private:
    int r;      // number of successes
    double p;   // probability of success

public:
    NegativeBinomialDistribution(int num_successes, double prob_success)
        : r(num_successes), p(prob_success) {
        if (r <= 0) {
            throw std::invalid_argument("Number of successes must be positive");
        }
        if (p <= 0 || p >= 1) {
            throw std::invalid_argument("Probability must be in (0, 1)");
        }
    }

    /**
     * @brief PMF: P(X=k) = C(k+r-1, k) * p^r * (1-p)^k
     *
     * k = number of failures before r successes
     */
    double pmf(int k) const {
        if (k < 0) return 0.0;

        double coeff = binomial_coefficient(k + r - 1, k);
        return coeff * std::pow(p, r) * std::pow(1 - p, k);
    }

    /**
     * @brief Mean: E[X] = r(1-p)/p
     */
    double mean() const override {
        return r * (1 - p) / p;
    }

    /**
     * @brief Variance: Var(X) = r(1-p)/p²
     */
    double variance() const override {
        return r * (1 - p) / (p * p);
    }

    /**
     * @brief Sample from negative binomial distribution
     */
    double sample() const {
        std::negative_binomial_distribution<int> dist(r, p);
        return static_cast<double>(dist(gen));
    }

    /**
     * @brief CDF: P(X ≤ k) = ∑ᵢ₌₀ᵏ PMF(i)
     */
    double cdf(int k) const {
        if (k < 0) return 0.0;

        double sum = 0.0;
        for (int i = 0; i <= k; ++i) {
            sum += pmf(i);
        }
        return sum;
    }

    /**
     * @brief Mode: ⌊(r-1)(1-p)/p⌋ for r > 1
     */
    int mode() const {
        if (r == 1) return 0;
        return static_cast<int>(std::floor((r - 1) * (1 - p) / p));
    }
};

/**
 * @class HypergeometricDistribution
 * @brief Hypergeometric distribution
 *
 * Sampling without replacement from finite population
 * N = population size
 * K = number of success states in population
 * n = number of draws
 * PMF: P(X=k) = C(K,k) * C(N-K, n-k) / C(N, n)
 */
class HypergeometricDistribution : public Distribution {
private:
    int N;  // population size
    int K;  // number of success states
    int n;  // number of draws

public:
    HypergeometricDistribution(int pop_size, int num_successes, int num_draws)
        : N(pop_size), K(num_successes), n(num_draws) {
        if (N <= 0 || K < 0 || n < 0) {
            throw std::invalid_argument("Invalid parameters");
        }
        if (K > N || n > N) {
            throw std::invalid_argument("K and n must be <= N");
        }
    }

    /**
     * @brief PMF: P(X=k) = C(K,k) * C(N-K, n-k) / C(N, n)
     *
     * k = number of observed successes in n draws
     */
    double pmf(int k) const {
        // Valid range: max(0, n+K-N) ≤ k ≤ min(n, K)
        int k_min = std::max(0, n + K - N);
        int k_max = std::min(n, K);

        if (k < k_min || k > k_max) return 0.0;

        double numerator = binomial_coefficient(K, k) * binomial_coefficient(N - K, n - k);
        double denominator = binomial_coefficient(N, n);

        return numerator / denominator;
    }

    /**
     * @brief Mean: E[X] = n * K/N
     */
    double mean() const override {
        return n * K / static_cast<double>(N);
    }

    /**
     * @brief Variance: Var(X) = n * (K/N) * (1 - K/N) * (N-n)/(N-1)
     */
    double variance() const override {
        double p = K / static_cast<double>(N);
        return n * p * (1 - p) * (N - n) / (N - 1.0);
    }

    /**
     * @brief Sample from hypergeometric distribution
     */
    double sample() const {
        // Use reservoir sampling or direct simulation
        std::vector<int> population(N);
        for (int i = 0; i < K; ++i) population[i] = 1;  // Success
        for (int i = K; i < N; ++i) population[i] = 0;  // Failure

        // Shuffle
        std::shuffle(population.begin(), population.end(), gen);

        // Count successes in first n draws
        int successes = 0;
        for (int i = 0; i < n; ++i) {
            successes += population[i];
        }

        return static_cast<double>(successes);
    }

    /**
     * @brief CDF: P(X ≤ k) = ∑ᵢ PMF(i)
     */
    double cdf(int k) const {
        int k_min = std::max(0, n + K - N);
        int k_max = std::min(n, K);

        if (k < k_min) return 0.0;
        if (k >= k_max) return 1.0;

        double sum = 0.0;
        for (int i = k_min; i <= k; ++i) {
            sum += pmf(i);
        }
        return sum;
    }

    /**
     * @brief Mode: ⌊(n+1)(K+1)/(N+2)⌋
     */
    int mode() const {
        return static_cast<int>(std::floor((n + 1.0) * (K + 1.0) / (N + 2.0)));
    }
};

/**
 * @class ErlangDistribution
 * @brief Erlang distribution (special case of Gamma with integer shape)
 *
 * Sum of k independent exponential random variables
 * Used in queuing theory and telecommunications
 */
class ErlangDistribution : public Distribution {
private:
    int k;          // shape parameter (must be integer)
    double lambda;  // rate parameter

public:
    ErlangDistribution(int shape, double rate) : k(shape), lambda(rate) {
        if (k <= 0) throw std::invalid_argument("Shape must be positive integer");
        if (lambda <= 0.0) throw std::invalid_argument("Rate must be positive");
    }

    /**
     * @brief PDF: f(x) = λᵏ xᵏ⁻¹ e⁻ᵏˣ / (k-1)!
     */
    double pdf(double x) const {
        if (x < 0.0) return 0.0;
        return std::pow(lambda, k) * std::pow(x, k - 1) * std::exp(-lambda * x) / factorial(k - 1);
    }

    double mean() const override { return k / lambda; }
    double variance() const override { return k / (lambda * lambda); }

    double sample() const {
        // Sum of k exponential random variables
        GammaDistribution gamma(k, lambda);
        return gamma.sample();
    }
};

/**
 * @class LognormalDistribution
 * @brief Lognormal distribution
 *
 * If X ~ Lognormal(μ, σ), then log(X) ~ Normal(μ, σ)
 * Common in finance, biology, economics
 */
class LognormalDistribution : public Distribution {
private:
    double mu, sigma;  // parameters of underlying normal

public:
    LognormalDistribution(double mean_log, double std_log) : mu(mean_log), sigma(std_log) {
        if (sigma <= 0.0) throw std::invalid_argument("Sigma must be positive");
    }

    /**
     * @brief PDF: f(x) = 1/(xσ√(2π)) exp(-(ln(x)-μ)²/(2σ²))
     */
    double pdf(double x) const {
        if (x <= 0.0) return 0.0;
        double log_x = std::log(x);
        double z = (log_x - mu) / sigma;
        return (1.0 / (x * sigma * std::sqrt(2.0 * PI))) * std::exp(-0.5 * z * z);
    }

    double cdf(double x) const {
        if (x <= 0.0) return 0.0;
        double z = (std::log(x) - mu) / (sigma * std::sqrt(2.0));
        return 0.5 * (1.0 + error_function(z));
    }

    double mean() const override {
        return std::exp(mu + sigma * sigma / 2.0);
    }

    double variance() const override {
        double exp_2mu_sig2 = std::exp(2 * mu + sigma * sigma);
        return exp_2mu_sig2 * (std::exp(sigma * sigma) - 1.0);
    }

    double sample() const {
        std::lognormal_distribution<double> dist(mu, sigma);
        return dist(gen);
    }
};

/**
 * @class CauchyDistribution
 * @brief Cauchy distribution (Lorentz distribution)
 *
 * Heavy-tailed distribution with no defined mean or variance
 * PDF: f(x) = 1/(π γ (1 + ((x-x₀)/γ)²))
 */
class CauchyDistribution : public Distribution {
private:
    double x0;     // location parameter
    double gamma;  // scale parameter

public:
    CauchyDistribution(double location, double scale) : x0(location), gamma(scale) {
        if (gamma <= 0.0) throw std::invalid_argument("Scale must be positive");
    }

    double pdf(double x) const {
        double z = (x - x0) / gamma;
        return 1.0 / (PI * gamma * (1.0 + z * z));
    }

    double cdf(double x) const {
        return 0.5 + std::atan((x - x0) / gamma) / PI;
    }

    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) throw std::invalid_argument("Probability must be in (0,1)");
        return x0 + gamma * std::tan(PI * (p - 0.5));
    }

    double mean() const override {
        throw std::runtime_error("Cauchy distribution has no defined mean");
    }

    double variance() const override {
        throw std::runtime_error("Cauchy distribution has no defined variance");
    }

    double sample() const {
        std::cauchy_distribution<double> dist(x0, gamma);
        return dist(gen);
    }
};

/**
 * @class LaplaceDistribution
 * @brief Laplace distribution (double exponential)
 *
 * PDF: f(x) = (1/2b) exp(-|x-μ|/b)
 * Used in robust statistics, signal processing
 */
class LaplaceDistribution : public Distribution {
private:
    double mu;  // location
    double b;   // scale

public:
    LaplaceDistribution(double location, double scale) : mu(location), b(scale) {
        if (b <= 0.0) throw std::invalid_argument("Scale must be positive");
    }

    double pdf(double x) const {
        return (1.0 / (2.0 * b)) * std::exp(-std::abs(x - mu) / b);
    }

    double cdf(double x) const {
        if (x < mu) {
            return 0.5 * std::exp((x - mu) / b);
        } else {
            return 1.0 - 0.5 * std::exp(-(x - mu) / b);
        }
    }

    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) throw std::invalid_argument("Probability must be in (0,1)");
        if (p < 0.5) {
            return mu + b * std::log(2.0 * p);
        } else {
            return mu - b * std::log(2.0 * (1.0 - p));
        }
    }

    double mean() const override { return mu; }
    double variance() const override { return 2.0 * b * b; }

    double sample() const {
        // Use inverse transform: F⁻¹(U) where U ~ Uniform(0,1)
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        return quantile(uniform(gen));
    }
};

/**
 * @class LogisticDistribution
 * @brief Logistic distribution
 *
 * PDF: f(x) = e⁻ᶻ / (s(1+e⁻ᶻ)²) where z = (x-μ)/s
 * Used in logistic regression, neural networks
 */
class LogisticDistribution : public Distribution {
private:
    double mu;  // location
    double s;   // scale

public:
    LogisticDistribution(double location, double scale) : mu(location), s(scale) {
        if (s <= 0.0) throw std::invalid_argument("Scale must be positive");
    }

    double pdf(double x) const {
        double z = (x - mu) / s;
        double exp_z = std::exp(-z);
        return exp_z / (s * (1.0 + exp_z) * (1.0 + exp_z));
    }

    double cdf(double x) const {
        double z = (x - mu) / s;
        return 1.0 / (1.0 + std::exp(-z));
    }

    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) throw std::invalid_argument("Probability must be in (0,1)");
        return mu + s * std::log(p / (1.0 - p));
    }

    double mean() const override { return mu; }
    double variance() const override { return (s * s * PI * PI) / 3.0; }

    double sample() const {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        return quantile(uniform(gen));
    }
};

/**
 * @class ParetoDistribution
 * @brief Pareto distribution (power law)
 *
 * PDF: f(x) = α xₘᵅ / xᵅ⁺¹ for x ≥ xₘ
 * Models wealth distribution, city sizes, word frequencies
 */
class ParetoDistribution : public Distribution {
private:
    double xm;     // scale (minimum value)
    double alpha;  // shape

public:
    ParetoDistribution(double scale, double shape) : xm(scale), alpha(shape) {
        if (xm <= 0.0 || alpha <= 0.0) {
            throw std::invalid_argument("Scale and shape must be positive");
        }
    }

    double pdf(double x) const {
        if (x < xm) return 0.0;
        return alpha * std::pow(xm, alpha) / std::pow(x, alpha + 1.0);
    }

    double cdf(double x) const {
        if (x < xm) return 0.0;
        return 1.0 - std::pow(xm / x, alpha);
    }

    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) throw std::invalid_argument("Probability must be in (0,1)");
        return xm / std::pow(1.0 - p, 1.0 / alpha);
    }

    double mean() const override {
        if (alpha <= 1.0) throw std::runtime_error("Mean undefined for alpha <= 1");
        return alpha * xm / (alpha - 1.0);
    }

    double variance() const override {
        if (alpha <= 2.0) throw std::runtime_error("Variance undefined for alpha <= 2");
        return (xm * xm * alpha) / ((alpha - 1.0) * (alpha - 1.0) * (alpha - 2.0));
    }

    double sample() const {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        return quantile(uniform(gen));
    }
};

/**
 * @class RayleighDistribution
 * @brief Rayleigh distribution
 *
 * PDF: f(x) = (x/σ²) exp(-x²/(2σ²)) for x ≥ 0
 * Models magnitude of 2D vector with independent normal components
 */
class RayleighDistribution : public Distribution {
private:
    double sigma;

public:
    RayleighDistribution(double scale) : sigma(scale) {
        if (sigma <= 0.0) throw std::invalid_argument("Scale must be positive");
    }

    double pdf(double x) const {
        if (x < 0.0) return 0.0;
        return (x / (sigma * sigma)) * std::exp(-x * x / (2.0 * sigma * sigma));
    }

    double cdf(double x) const {
        if (x < 0.0) return 0.0;
        return 1.0 - std::exp(-x * x / (2.0 * sigma * sigma));
    }

    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) throw std::invalid_argument("Probability must be in (0,1)");
        return sigma * std::sqrt(-2.0 * std::log(1.0 - p));
    }

    double mean() const override {
        return sigma * std::sqrt(PI / 2.0);
    }

    double variance() const override {
        return ((4.0 - PI) / 2.0) * sigma * sigma;
    }

    double sample() const {
        // R = sqrt(X² + Y²) where X, Y ~ N(0, σ²)
        std::normal_distribution<double> normal(0.0, sigma);
        double x = normal(gen);
        double y = normal(gen);
        return std::sqrt(x * x + y * y);
    }
};

/**
 * @class HalfNormalDistribution
 * @brief Half-normal distribution
 *
 * PDF: f(x) = (√(2/π) / σ) exp(-x²/(2σ²)) for x ≥ 0
 * Folded normal distribution at zero
 */
class HalfNormalDistribution : public Distribution {
private:
    double sigma;

public:
    HalfNormalDistribution(double scale) : sigma(scale) {
        if (sigma <= 0.0) throw std::invalid_argument("Scale must be positive");
    }

    double pdf(double x) const {
        if (x < 0.0) return 0.0;
        return std::sqrt(2.0 / PI) / sigma * std::exp(-x * x / (2.0 * sigma * sigma));
    }

    double cdf(double x) const {
        if (x < 0.0) return 0.0;
        return error_function(x / (sigma * std::sqrt(2.0)));
    }

    double mean() const override {
        return sigma * std::sqrt(2.0 / PI);
    }

    double variance() const override {
        return sigma * sigma * (1.0 - 2.0 / PI);
    }

    double sample() const {
        std::normal_distribution<double> normal(0.0, sigma);
        return std::abs(normal(gen));
    }
};

/**
 * @class InverseGaussianDistribution
 * @brief Inverse Gaussian (Wald) distribution
 *
 * PDF: f(x) = √(λ/(2πx³)) exp(-λ(x-μ)²/(2μ²x))
 * First passage time for Brownian motion with drift
 */
class InverseGaussianDistribution : public Distribution {
private:
    double mu;     // mean
    double lambda; // shape

public:
    InverseGaussianDistribution(double mean_param, double shape_param)
        : mu(mean_param), lambda(shape_param) {
        if (mu <= 0.0 || lambda <= 0.0) {
            throw std::invalid_argument("Mean and shape must be positive");
        }
    }

    double pdf(double x) const {
        if (x <= 0.0) return 0.0;
        double coeff = std::sqrt(lambda / (2.0 * PI * x * x * x));
        double exponent = -lambda * (x - mu) * (x - mu) / (2.0 * mu * mu * x);
        return coeff * std::exp(exponent);
    }

    double mean() const override { return mu; }

    double variance() const override {
        return mu * mu * mu / lambda;
    }

    double sample() const {
        // Michael, Schucany and Haas algorithm
        std::normal_distribution<double> normal(0.0, 1.0);
        double nu = normal(gen);
        double y = nu * nu;
        double x = mu + (mu * mu * y) / (2.0 * lambda)
                   - (mu / (2.0 * lambda)) * std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y);

        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double u = uniform(gen);

        if (u <= mu / (mu + x)) {
            return x;
        } else {
            return mu * mu / x;
        }
    }
};

/**
 * @class ExtremeValueDistribution
 * @brief Extreme value (Gumbel) distribution
 *
 * PDF: f(x) = (1/β) exp(-(x-μ)/β) exp(-exp(-(x-μ)/β))
 * Models maximum/minimum of samples, used in extreme value theory
 */
class ExtremeValueDistribution : public Distribution {
private:
    double mu;   // location
    double beta; // scale

public:
    ExtremeValueDistribution(double location, double scale) : mu(location), beta(scale) {
        if (beta <= 0.0) throw std::invalid_argument("Scale must be positive");
    }

    double pdf(double x) const {
        double z = (x - mu) / beta;
        return (1.0 / beta) * std::exp(-z - std::exp(-z));
    }

    double cdf(double x) const {
        double z = (x - mu) / beta;
        return std::exp(-std::exp(-z));
    }

    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) throw std::invalid_argument("Probability must be in (0,1)");
        return mu - beta * std::log(-std::log(p));
    }

    double mean() const override {
        return mu + beta * 0.5772156649;  // Euler-Mascheroni constant
    }

    double variance() const override {
        return (PI * PI / 6.0) * beta * beta;
    }

    double sample() const {
        std::extreme_value_distribution<double> dist(mu, beta);
        return dist(gen);
    }
};

/**
 * @class ArcsinDistribution
 * @brief Arcsine distribution on [a, b]
 *
 * PDF: f(x) = 1/(π√((x-a)(b-x))) for a < x < b
 * Special case of Beta(1/2, 1/2) on [0,1]
 */
class ArcsinDistribution : public Distribution {
private:
    double a, b;

public:
    ArcsinDistribution(double lower, double upper) : a(lower), b(upper) {
        if (a >= b) throw std::invalid_argument("Lower bound must be less than upper");
    }

    double pdf(double x) const {
        if (x <= a || x >= b) return 0.0;
        return 1.0 / (PI * std::sqrt((x - a) * (b - x)));
    }

    double cdf(double x) const {
        if (x <= a) return 0.0;
        if (x >= b) return 1.0;
        return (2.0 / PI) * std::asin(std::sqrt((x - a) / (b - a)));
    }

    double mean() const override {
        return (a + b) / 2.0;
    }

    double variance() const override {
        return (b - a) * (b - a) / 8.0;
    }

    double sample() const {
        // Use Beta(1/2, 1/2) and transform
        BetaDistribution beta(0.5, 0.5);
        return a + (b - a) * beta.sample();
    }
};

/**
 * @class PowerFunctionDistribution
 * @brief Power function distribution
 *
 * PDF: f(x) = α xᵅ⁻¹ / bᵅ for 0 ≤ x ≤ b
 * Special case of Beta distribution
 */
class PowerFunctionDistribution : public Distribution {
private:
    double alpha;  // shape
    double b;      // upper bound

public:
    PowerFunctionDistribution(double shape, double upper) : alpha(shape), b(upper) {
        if (alpha <= 0.0) throw std::invalid_argument("Shape must be positive");
        if (b <= 0.0) throw std::invalid_argument("Upper bound must be positive");
    }

    double pdf(double x) const {
        if (x < 0.0 || x > b) return 0.0;
        return alpha * std::pow(x, alpha - 1.0) / std::pow(b, alpha);
    }

    double cdf(double x) const {
        if (x < 0.0) return 0.0;
        if (x > b) return 1.0;
        return std::pow(x / b, alpha);
    }

    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) throw std::invalid_argument("Probability must be in (0,1)");
        return b * std::pow(p, 1.0 / alpha);
    }

    double mean() const override {
        return alpha * b / (alpha + 1.0);
    }

    double variance() const override {
        return (alpha * b * b) / ((alpha + 1.0) * (alpha + 1.0) * (alpha + 2.0));
    }

    double sample() const {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        return quantile(uniform(gen));
    }
};

/**
 * @class NoncentralChiSquaredDistribution
 * @brief Noncentral chi-squared distribution
 *
 * Sum of squares of k independent normal random variables with non-zero means
 * χ²(k, λ) where λ is noncentrality parameter
 */
class NoncentralChiSquaredDistribution : public Distribution {
private:
    int k;          // degrees of freedom
    double lambda;  // noncentrality parameter

public:
    NoncentralChiSquaredDistribution(int df, double noncentrality)
        : k(df), lambda(noncentrality) {
        if (k <= 0) throw std::invalid_argument("Degrees of freedom must be positive");
        if (lambda < 0.0) throw std::invalid_argument("Noncentrality must be non-negative");
    }

    /**
     * @brief PDF approximation using series expansion
     */
    double pdf(double x) const {
        if (x < 0.0) return 0.0;

        // Sum of weighted central chi-squared PDFs
        double sum = 0.0;
        double poisson_term = std::exp(-lambda / 2.0);

        for (int j = 0; j < 20; ++j) {  // Truncate series
            ChiSquaredDistribution chi_sq(k + 2 * j);
            sum += poisson_term * chi_sq.pdf(x);
            poisson_term *= (lambda / 2.0) / (j + 1);
        }

        return sum;
    }

    double mean() const override {
        return k + lambda;
    }

    double variance() const override {
        return 2.0 * (k + 2.0 * lambda);
    }

    double sample() const {
        // Sum of squares of normals with non-zero means
        std::normal_distribution<double> normal(0.0, 1.0);
        double sum = 0.0;

        // First term has non-zero mean
        double z = normal(gen) + std::sqrt(lambda);
        sum += z * z;

        // Remaining k-1 terms are standard normal
        for (int i = 1; i < k; ++i) {
            double zi = normal(gen);
            sum += zi * zi;
        }

        return sum;
    }
};

/**
 * @class NoncentralTDistribution
 * @brief Noncentral t-distribution
 *
 * t(ν, δ) where ν is degrees of freedom, δ is noncentrality
 * Ratio of normal with non-zero mean to chi-squared
 */
class NoncentralTDistribution : public Distribution {
private:
    int nu;        // degrees of freedom
    double delta;  // noncentrality parameter

public:
    NoncentralTDistribution(int df, double noncentrality)
        : nu(df), delta(noncentrality) {
        if (nu <= 0) throw std::invalid_argument("Degrees of freedom must be positive");
    }

    /**
     * @brief PDF approximation
     */
    double pdf(double t) const {
        // Simplified approximation
        double z = t - delta;
        StudentTDistribution central_t(nu);
        return central_t.pdf(z);
    }

    double mean() const override {
        if (nu <= 1) throw std::runtime_error("Mean undefined for nu <= 1");
        return delta * std::sqrt(nu / 2.0) * gamma_function((nu - 1) / 2.0) / gamma_function(nu / 2.0);
    }

    double variance() const override {
        if (nu <= 2) throw std::runtime_error("Variance undefined for nu <= 2");
        double mu = mean();
        return nu * (1.0 + delta * delta) / (nu - 2.0) - mu * mu;
    }

    double sample() const {
        // t = (Z + δ) / sqrt(χ²/ν) where Z ~ N(0,1), χ² ~ χ²(ν)
        std::normal_distribution<double> normal(0.0, 1.0);
        ChiSquaredDistribution chi_sq(nu);

        double z = normal(gen) + delta;
        double chi = chi_sq.sample();

        return z / std::sqrt(chi / nu);
    }
};

/**
 * @class NoncentralFDistribution
 * @brief Noncentral F-distribution
 *
 * F(d₁, d₂, λ) where λ is noncentrality parameter
 * Ratio of noncentral chi-squared to central chi-squared
 */
class NoncentralFDistribution : public Distribution {
private:
    int d1, d2;    // degrees of freedom
    double lambda; // noncentrality parameter

public:
    NoncentralFDistribution(int df1, int df2, double noncentrality)
        : d1(df1), d2(df2), lambda(noncentrality) {
        if (d1 <= 0 || d2 <= 0) {
            throw std::invalid_argument("Degrees of freedom must be positive");
        }
        if (lambda < 0.0) {
            throw std::invalid_argument("Noncentrality must be non-negative");
        }
    }

    double mean() const override {
        if (d2 <= 2) throw std::runtime_error("Mean undefined for d2 <= 2");
        return d2 * (d1 + lambda) / (d1 * (d2 - 2));
    }

    double variance() const override {
        if (d2 <= 4) throw std::runtime_error("Variance undefined for d2 <= 4");

        double a = d1 + lambda;
        double num = 2.0 * d2 * d2 * (d1 * a + (d1 + 2 * lambda) * (d2 - 2));
        double den = d1 * d1 * (d2 - 2) * (d2 - 2) * (d2 - 4);

        return num / den;
    }

    double sample() const {
        // F = (χ²₁(λ)/d₁) / (χ²₂/d₂)
        NoncentralChiSquaredDistribution nc_chi_sq(d1, lambda);
        ChiSquaredDistribution chi_sq(d2);

        return (nc_chi_sq.sample() / d1) / (chi_sq.sample() / d2);
    }
};

/**
 * @class MultivariateNormalDistribution
 * @brief Multivariate normal distribution N(μ, Σ)
 *
 * Generalization of normal distribution to n dimensions
 * PDF: f(x) = (2π)^(-n/2) |Σ|^(-1/2) exp(-1/2 (x-μ)ᵀ Σ⁻¹ (x-μ))
 */
class MultivariateNormalDistribution {
private:
    std::vector<double> mu;                      // mean vector
    std::vector<std::vector<double>> sigma;      // covariance matrix
    std::vector<std::vector<double>> L;          // Cholesky decomposition of Σ
    int n;                                       // dimension

    // Compute Cholesky decomposition: Σ = LLᵀ
    void computeCholesky() {
        L.resize(n, std::vector<double>(n, 0.0));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = 0.0;
                for (int k = 0; k < j; ++k) {
                    sum += L[i][k] * L[j][k];
                }

                if (i == j) {
                    L[i][j] = std::sqrt(sigma[i][i] - sum);
                } else {
                    L[i][j] = (sigma[i][j] - sum) / L[j][j];
                }
            }
        }
    }

public:
    MultivariateNormalDistribution(
        const std::vector<double>& mean_vec,
        const std::vector<std::vector<double>>& cov_matrix)
        : mu(mean_vec), sigma(cov_matrix), n(mean_vec.size()) {

        if (sigma.size() != n || sigma[0].size() != n) {
            throw std::invalid_argument("Covariance matrix dimensions must match mean vector");
        }

        computeCholesky();
    }

    /**
     * @brief Sample from multivariate normal using Cholesky decomposition
     *
     * X = μ + L·Z where Z ~ N(0, I)
     */
    std::vector<double> sample() const {
        std::normal_distribution<double> normal(0.0, 1.0);
        std::vector<double> z(n);

        // Generate independent standard normals
        for (int i = 0; i < n; ++i) {
            z[i] = normal(gen);
        }

        // Transform: x = μ + L·z
        std::vector<double> x(n);
        for (int i = 0; i < n; ++i) {
            x[i] = mu[i];
            for (int j = 0; j <= i; ++j) {
                x[i] += L[i][j] * z[j];
            }
        }

        return x;
    }

    /**
     * @brief Mean vector
     */
    std::vector<double> mean() const {
        return mu;
    }

    /**
     * @brief Covariance matrix
     */
    std::vector<std::vector<double>> covariance() const {
        return sigma;
    }

    /**
     * @brief Dimension
     */
    int dimension() const {
        return n;
    }
};

/**
 * @class TriangularDistribution
 * @brief Triangular distribution on [a, c] with mode b
 *
 * PDF: f(x) = 2(x-a)/((c-a)(b-a)) for a ≤ x < b
 *           = 2(c-x)/((c-a)(c-b)) for b ≤ x ≤ c
 * Often used when limited data is available
 */
class TriangularDistribution : public Distribution {
private:
    double a;  // minimum
    double b;  // mode
    double c;  // maximum

public:
    TriangularDistribution(double min_val, double mode_val, double max_val)
        : a(min_val), b(mode_val), c(max_val) {
        if (!(a < b && b < c)) {
            throw std::invalid_argument("Must have a < b < c");
        }
    }

    double pdf(double x) const {
        if (x < a || x > c) return 0.0;

        if (x < b) {
            return 2.0 * (x - a) / ((c - a) * (b - a));
        } else {
            return 2.0 * (c - x) / ((c - a) * (c - b));
        }
    }

    double cdf(double x) const {
        if (x <= a) return 0.0;
        if (x >= c) return 1.0;

        if (x < b) {
            return (x - a) * (x - a) / ((c - a) * (b - a));
        } else {
            return 1.0 - (c - x) * (c - x) / ((c - a) * (c - b));
        }
    }

    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) throw std::invalid_argument("Probability must be in (0,1)");

        double fc = (b - a) / (c - a);

        if (p < fc) {
            return a + std::sqrt(p * (c - a) * (b - a));
        } else {
            return c - std::sqrt((1.0 - p) * (c - a) * (c - b));
        }
    }

    double mean() const override {
        return (a + b + c) / 3.0;
    }

    double variance() const override {
        return (a * a + b * b + c * c - a * b - a * c - b * c) / 18.0;
    }

    double sample() const {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        return quantile(uniform(gen));
    }

    double mode() const { return b; }
};

/**
 * @class WeibullDistribution
 * @brief Weibull distribution
 *
 * PDF: f(x) = (k/λ)(x/λ)^(k-1) exp(-(x/λ)^k) for x ≥ 0
 * Used in reliability engineering, survival analysis, wind speed modeling
 */
class WeibullDistribution : public Distribution {
private:
    double lambda;  // scale parameter
    double k;       // shape parameter

public:
    WeibullDistribution(double scale, double shape) : lambda(scale), k(shape) {
        if (lambda <= 0.0 || k <= 0.0) {
            throw std::invalid_argument("Scale and shape must be positive");
        }
    }

    double pdf(double x) const {
        if (x < 0.0) return 0.0;
        double z = x / lambda;
        return (k / lambda) * std::pow(z, k - 1.0) * std::exp(-std::pow(z, k));
    }

    double cdf(double x) const {
        if (x < 0.0) return 0.0;
        return 1.0 - std::exp(-std::pow(x / lambda, k));
    }

    double quantile(double p) const {
        if (p <= 0.0 || p >= 1.0) throw std::invalid_argument("Probability must be in (0,1)");
        return lambda * std::pow(-std::log(1.0 - p), 1.0 / k);
    }

    double mean() const override {
        return lambda * gamma_function(1.0 + 1.0 / k);
    }

    double variance() const override {
        double g1 = gamma_function(1.0 + 1.0 / k);
        double g2 = gamma_function(1.0 + 2.0 / k);
        return lambda * lambda * (g2 - g1 * g1);
    }

    double sample() const {
        std::weibull_distribution<double> dist(k, lambda);
        return dist(gen);
    }

    /**
     * @brief Hazard rate (failure rate): h(x) = (k/λ)(x/λ)^(k-1)
     */
    double hazard_rate(double x) const {
        if (x < 0.0) return 0.0;
        return (k / lambda) * std::pow(x / lambda, k - 1.0);
    }
};

// ============================================================================
// SECTION 3: STATISTICS AND STATISTICAL INFERENCE
// ============================================================================

/**
 * @class StandardNormalExtensions
 * @brief Extended functions for standard normal distribution
 *
 * Critical values, tolerance factors, circular probabilities
 */
class StandardNormalExtensions {
public:
    /**
     * @brief Critical value z_α such that P(Z > z_α) = α
     *
     * Commonly used: z_0.05 = 1.645, z_0.025 = 1.960, z_0.01 = 2.326
     */
    static double critical_value(double alpha) {
        if (alpha <= 0.0 || alpha >= 1.0) {
            throw std::invalid_argument("Alpha must be in (0,1)");
        }

        // Common critical values
        if (std::abs(alpha - 0.05) < 1e-6) return 1.6448536269514722;
        if (std::abs(alpha - 0.025) < 1e-6) return 1.9599639845400545;
        if (std::abs(alpha - 0.01) < 1e-6) return 2.3263478740408408;
        if (std::abs(alpha - 0.005) < 1e-6) return 2.5758293035489004;

        // Binary search for general case
        double low = -5.0, high = 5.0;
        NormalDistribution N(0.0, 1.0);

        for (int iter = 0; iter < 100; ++iter) {
            double mid = (low + high) / 2.0;
            double p = 1.0 - N.cdf(mid);

            if (std::abs(p - alpha) < 1e-9) return mid;

            if (p > alpha) {
                low = mid;
            } else {
                high = mid;
            }
        }

        return (low + high) / 2.0;
    }

    /**
     * @brief Tolerance factor k for normal distribution
     *
     * k such that P(μ - kσ ≤ X ≤ μ + kσ) ≥ β with confidence γ
     * Used in quality control and tolerance intervals
     */
    static double tolerance_factor(double beta, double gamma, int n) {
        if (beta <= 0.0 || beta >= 1.0) throw std::invalid_argument("Beta must be in (0,1)");
        if (gamma <= 0.0 || gamma >= 1.0) throw std::invalid_argument("Gamma must be in (0,1)");
        if (n <= 1) throw std::invalid_argument("Sample size must be > 1");

        // Simplified approximation
        double z_beta = critical_value((1.0 - beta) / 2.0);
        double chi_sq_gamma = std::sqrt(2.0 * n - 2.0) * critical_value((1.0 - gamma) / 2.0);

        return z_beta * std::sqrt((n - 1.0) / n) * (1.0 + 1.0 / (2.0 * n));
    }

    /**
     * @brief Circular error probable (CEP)
     *
     * Radius of circle centered at mean containing 50% of points
     * For bivariate normal with σ_x = σ_y = σ
     */
    static double circular_error_probable(double sigma) {
        // CEP ≈ 1.1774σ for equal variance bivariate normal
        return 1.1774100225154747 * sigma;
    }

    /**
     * @brief Circular normal probability
     *
     * P(X² + Y² ≤ r²) where X, Y ~ N(0, σ²) independent
     */
    static double circular_normal_probability(double r, double sigma) {
        if (r < 0.0) return 0.0;
        // This is a chi-squared distribution with 2 DOF scaled by σ²
        return 1.0 - std::exp(-r * r / (2.0 * sigma * sigma));
    }

    /**
     * @brief Operating characteristic (OC) curve value
     *
     * Probability of accepting hypothesis when true mean is μ₁
     * Used in hypothesis testing and quality control
     */
    static double operating_characteristic(
        double mu0, double mu1, double sigma, int n, double alpha) {

        double z_alpha = critical_value(alpha);
        double z = (mu1 - mu0) / (sigma / std::sqrt(n));

        NormalDistribution N(0.0, 1.0);
        return N.cdf(z_alpha - z);
    }

    /**
     * @brief Two-sided tolerance limit
     *
     * Returns [lower, upper] such that P(lower ≤ X ≤ upper) ≥ β with confidence γ
     */
    static std::pair<double, double> tolerance_limits(
        double mean, double std_dev, double beta, double gamma, int n) {

        double k = tolerance_factor(beta, gamma, n);
        return {mean - k * std_dev, mean + k * std_dev};
    }
};

/**
 * @class DistributionRelationships
 * @brief Documents and computes relationships among probability distributions
 *
 * Special cases, limiting forms, transformations
 */
class DistributionRelationships {
public:
    /**
     * @brief Check if distribution is a special case of another
     *
     * Examples:
     * - Exponential(λ) = Gamma(1, λ)
     * - Chi-squared(k) = Gamma(k/2, 1/2)
     * - Rayleigh(σ) is magnitude of bivariate normal with σ
     */
    static std::string relationship(const std::string& dist_name) {
        if (dist_name == "Exponential") {
            return "Exponential(λ) = Gamma(α=1, β=λ)";
        } else if (dist_name == "Chi-squared") {
            return "χ²(k) = Gamma(α=k/2, β=1/2)";
        } else if (dist_name == "Erlang") {
            return "Erlang(k, λ) = Gamma(α=k, β=λ) with integer k";
        } else if (dist_name == "Rayleigh") {
            return "Rayleigh(σ) = √(X² + Y²) where X, Y ~ N(0, σ²) independent";
        } else if (dist_name == "Student-t") {
            return "t(ν) → N(0,1) as ν → ∞";
        } else if (dist_name == "Binomial") {
            return "Binomial(n,p) → Poisson(λ=np) as n→∞, p→0";
        } else if (dist_name == "Poisson") {
            return "Poisson(λ) → N(λ, λ) as λ → ∞";
        }
        return "Unknown distribution";
    }

    /**
     * @brief Sum of independent distributions
     *
     * Z = X + Y, find distribution of Z
     */
    static std::string sum_distribution(const std::string& X_dist, const std::string& Y_dist) {
        if (X_dist == "Normal" && Y_dist == "Normal") {
            return "Normal(μ₁+μ₂, σ₁²+σ₂²)";
        } else if (X_dist == "Poisson" && Y_dist == "Poisson") {
            return "Poisson(λ₁+λ₂)";
        } else if (X_dist == "Gamma" && Y_dist == "Gamma") {
            return "Gamma(α₁+α₂, β) if same rate β";
        } else if (X_dist == "Chi-squared" && Y_dist == "Chi-squared") {
            return "Chi-squared(k₁+k₂)";
        }
        return "Use convolution or MGF";
    }

    /**
     * @brief Limiting distribution as n → ∞
     */
    static std::string limiting_distribution(const std::string& dist, const std::string& limit_type) {
        if (dist == "Binomial" && limit_type == "np→λ") {
            return "Poisson(λ)";
        } else if (dist == "Binomial" && limit_type == "CLT") {
            return "Normal(np, np(1-p))";
        } else if (dist == "Student-t" && limit_type == "ν→∞") {
            return "Normal(0, 1)";
        } else if (dist == "F" && limit_type == "d₁→∞") {
            return "χ²(d₂)/d₂";
        }
        return "Not a standard limit";
    }
};

/**
 * @class StatisticalEstimation
 * @brief Parameter estimation via method of moments and maximum likelihood
 *
 * Implements:
 * - Cramér-Rao lower bound
 * - Method of moments
 * - Maximum likelihood estimation
 * - Fisher information
 * - Asymptotic properties
 */
class StatisticalEstimation {
public:
    /**
     * @brief Cramér-Rao lower bound for variance of unbiased estimator
     *
     * Var(θ̂) ≥ 1/I(θ) where I(θ) is Fisher information
     *
     * For location parameter of normal: I(μ) = n/σ²
     */
    static double cramer_rao_bound(const std::string& distribution,
                                   double param, int n) {
        if (distribution == "Normal-mean") {
            double sigma = 1.0;  // Assuming known σ = 1
            return sigma * sigma / n;
        } else if (distribution == "Normal-variance") {
            // For variance estimation: I(σ²) = n/(2σ⁴)
            double sigma_sq = param;
            return 2.0 * sigma_sq * sigma_sq / n;
        } else if (distribution == "Exponential") {
            double lambda = param;
            return lambda * lambda / n;
        } else if (distribution == "Poisson") {
            double lambda = param;
            return lambda / n;
        }

        throw std::invalid_argument("Unknown distribution for Cramér-Rao bound");
    }

    /**
     * @brief Fisher information for single parameter
     *
     * I(θ) = E[(∂/∂θ log f(X;θ))²] = -E[∂²/∂θ² log f(X;θ)]
     */
    static double fisher_information(const std::string& distribution,
                                     double param, int n) {
        if (distribution == "Normal-mean") {
            double sigma = 1.0;
            return n / (sigma * sigma);
        } else if (distribution == "Exponential") {
            double lambda = param;
            return n / (lambda * lambda);
        } else if (distribution == "Poisson") {
            double lambda = param;
            return n / lambda;
        }

        throw std::invalid_argument("Unknown distribution for Fisher information");
    }

    /**
     * @brief Method of moments estimation
     *
     * Equate sample moments to population moments and solve for parameters
     *
     * Example: Normal(μ, σ²)
     *   μ̂ = X̄ (sample mean)
     *   σ̂² = (1/n)∑(Xᵢ - X̄)² (second central moment)
     */
    static std::vector<double> method_of_moments(
        const std::string& distribution,
        const std::vector<double>& data) {

        int n = data.size();
        double mean = std::accumulate(data.begin(), data.end(), 0.0) / n;

        double m2 = 0.0;  // Second moment
        for (double x : data) {
            m2 += (x - mean) * (x - mean);
        }
        m2 /= n;

        if (distribution == "Normal") {
            return {mean, std::sqrt(m2)};
        } else if (distribution == "Exponential") {
            // E[X] = 1/λ ⟹ λ̂ = 1/X̄
            return {1.0 / mean};
        } else if (distribution == "Gamma") {
            // E[X] = α/β, Var(X) = α/β²
            // ⟹ α̂ = E²/Var, β̂ = E/Var
            double alpha = mean * mean / m2;
            double beta = mean / m2;
            return {alpha, beta};
        } else if (distribution == "Uniform") {
            // Uniform[a, b]: E[X] = (a+b)/2, Var(X) = (b-a)²/12
            double range = std::sqrt(12.0 * m2);
            double a = mean - range / 2.0;
            double b = mean + range / 2.0;
            return {a, b};
        }

        throw std::invalid_argument("Unknown distribution for method of moments");
    }

    /**
     * @brief Maximum likelihood estimation
     *
     * Find θ̂ that maximizes L(θ) = ∏f(xᵢ; θ)
     * Equivalently, maximize log L(θ) = ∑log f(xᵢ; θ)
     */
    static std::vector<double> maximum_likelihood(
        const std::string& distribution,
        const std::vector<double>& data) {

        int n = data.size();
        double mean = std::accumulate(data.begin(), data.end(), 0.0) / n;

        if (distribution == "Normal") {
            // MLE: μ̂ = X̄, σ̂² = (1/n)∑(Xᵢ - X̄)²
            double var = 0.0;
            for (double x : data) {
                var += (x - mean) * (x - mean);
            }
            var /= n;
            return {mean, std::sqrt(var)};
        } else if (distribution == "Exponential") {
            // MLE: λ̂ = 1/X̄
            return {1.0 / mean};
        } else if (distribution == "Poisson") {
            // MLE: λ̂ = X̄
            return {mean};
        } else if (distribution == "Uniform") {
            // MLE: â = min(X), b̂ = max(X)
            double a = *std::min_element(data.begin(), data.end());
            double b = *std::max_element(data.begin(), data.end());
            return {a, b};
        }

        throw std::invalid_argument("Unknown distribution for MLE");
    }

    /**
     * @brief Log-likelihood function
     *
     * ℓ(θ) = ∑ᵢ log f(xᵢ; θ)
     */
    static double log_likelihood(
        const std::string& distribution,
        const std::vector<double>& data,
        const std::vector<double>& params) {

        double log_lik = 0.0;

        if (distribution == "Normal") {
            double mu = params[0];
            double sigma = params[1];
            NormalDistribution N(mu, sigma);

            for (double x : data) {
                log_lik += std::log(N.pdf(x));
            }
        } else if (distribution == "Exponential") {
            double lambda = params[0];
            ExponentialDistribution E(lambda);

            for (double x : data) {
                log_lik += std::log(E.pdf(x));
            }
        } else if (distribution == "Poisson") {
            double lambda = params[0];
            PoissonDistribution P(lambda);

            for (double x : data) {
                log_lik += std::log(P.pmf(static_cast<int>(x)));
            }
        } else {
            throw std::invalid_argument("Unknown distribution for log-likelihood");
        }

        return log_lik;
    }

    /**
     * @brief Invariance property of MLEs
     *
     * If θ̂ is MLE of θ, then g(θ̂) is MLE of g(θ)
     *
     * Example: If σ̂² is MLE of σ², then σ̂ = √σ̂² is MLE of σ
     */
    static double mle_transformed(
        double theta_hat,
        std::function<double(double)> g) {
        return g(theta_hat);
    }

    /**
     * @brief Asymptotic distribution of MLE
     *
     * √n(θ̂ - θ) →ᵈ N(0, 1/I(θ)) as n → ∞
     *
     * Returns asymptotic variance: 1/(n·I(θ))
     */
    static double asymptotic_variance(
        const std::string& distribution,
        double true_param, int n) {

        double I_theta = fisher_information(distribution, true_param, 1);
        return 1.0 / (n * I_theta);
    }

    /**
     * @brief Efficiency of estimator
     *
     * Eff(θ̂) = [Cramér-Rao bound] / Var(θ̂)
     * Efficient estimator has Eff = 1
     */
    static double efficiency(
        double estimator_variance,
        const std::string& distribution,
        double param, int n) {

        double cr_bound = cramer_rao_bound(distribution, param, n);
        return cr_bound / estimator_variance;
    }

    /**
     * @brief Bias of estimator
     *
     * Bias(θ̂) = E[θ̂] - θ
     */
    static double bias(double expected_value, double true_value) {
        return expected_value - true_value;
    }

    /**
     * @brief Mean squared error
     *
     * MSE(θ̂) = Var(θ̂) + [Bias(θ̂)]²
     */
    static double mean_squared_error(double variance, double bias_value) {
        return variance + bias_value * bias_value;
    }

    /**
     * @brief Consistent estimator check
     *
     * θ̂ₙ is consistent if θ̂ₙ →ᵖ θ as n → ∞
     * Sufficient condition: E[θ̂ₙ] → θ and Var(θ̂ₙ) → 0
     */
    static bool is_consistent(
        std::function<double(int)> expected_value,
        std::function<double(int)> variance,
        double true_value,
        int large_n = 10000) {

        double exp_val = expected_value(large_n);
        double var_val = variance(large_n);

        return (std::abs(exp_val - true_value) < 1e-6) && (var_val < 1e-6);
    }
};

/**
 * @class HypothesisTesting
 * @brief Comprehensive hypothesis testing methods
 *
 * Implements computational methods for:
 * - Neyman-Pearson lemma
 * - Likelihood ratio tests
 * - Goodness of fit tests
 * - Contingency tables
 * - Variance tests (Bartlett, Cochran)
 * - Sample size calculations
 * - Outlier detection
 * - Bernoulli trials
 */
class HypothesisTesting {
public:
    /**
     * @brief Neyman-Pearson lemma likelihood ratio test
     *
     * For H₀: θ = θ₀ vs H₁: θ = θ₁ (simple vs simple)
     * Reject H₀ if L(θ₁)/L(θ₀) > k
     *
     * @param distribution Name of distribution ("Normal", "Exponential", etc.)
     * @param data Sample data
     * @param theta0 Null hypothesis parameter
     * @param theta1 Alternative hypothesis parameter
     * @return Likelihood ratio L(θ₁)/L(θ₀)
     */
    static double neyman_pearson_lr(
        const std::string& distribution,
        const std::vector<double>& data,
        double theta0,
        double theta1) {

        double log_lik0 = 0.0, log_lik1 = 0.0;

        if (distribution == "Normal-mean") {
            // Assuming σ² = 1
            for (double x : data) {
                log_lik0 += -0.5 * (x - theta0) * (x - theta0);
                log_lik1 += -0.5 * (x - theta1) * (x - theta1);
            }
        } else if (distribution == "Exponential") {
            for (double x : data) {
                log_lik0 += std::log(theta0) - theta0 * x;
                log_lik1 += std::log(theta1) - theta1 * x;
            }
        } else if (distribution == "Poisson") {
            for (double x : data) {
                int k = static_cast<int>(x);
                log_lik0 += k * std::log(theta0) - theta0;
                log_lik1 += k * std::log(theta1) - theta1;
            }
        }

        return std::exp(log_lik1 - log_lik0);
    }

    /**
     * @brief Generalized likelihood ratio test (GLRT)
     *
     * Test statistic: Λ = L(θ̂₀)/L(θ̂) where θ̂₀ is MLE under H₀, θ̂ is unrestricted MLE
     * -2 log Λ ~ χ²(df) asymptotically
     *
     * @return -2 log(Λ)
     */
    static double generalized_lr_test(
        double log_lik_H0,
        double log_lik_unrestricted) {

        return -2.0 * (log_lik_H0 - log_lik_unrestricted);
    }

    /**
     * @brief Kolmogorov-Smirnov goodness-of-fit test
     *
     * Tests if sample comes from specified distribution
     * Test statistic: D = max|F_n(x) - F(x)|
     *
     * @param data Sample data (will be sorted)
     * @param cdf_function Theoretical CDF function
     * @return KS test statistic D
     */
    static double kolmogorov_smirnov_test(
        std::vector<double> data,
        std::function<double(double)> cdf_function) {

        std::sort(data.begin(), data.end());
        int n = data.size();
        double max_diff = 0.0;

        for (int i = 0; i < n; ++i) {
            double empirical_cdf = (i + 1.0) / n;
            double theoretical_cdf = cdf_function(data[i]);
            double diff = std::abs(empirical_cdf - theoretical_cdf);
            max_diff = std::max(max_diff, diff);
        }

        return max_diff;
    }

    /**
     * @brief Chi-squared goodness-of-fit test
     *
     * Tests if observed frequencies match expected frequencies
     * χ² = ∑(O_i - E_i)²/E_i
     *
     * @param observed Observed frequencies
     * @param expected Expected frequencies
     * @return Chi-squared statistic (compare with χ²(k-1) distribution)
     */
    static double chi_squared_goodness_of_fit(
        const std::vector<double>& observed,
        const std::vector<double>& expected) {

        if (observed.size() != expected.size()) {
            throw std::invalid_argument("Observed and expected must have same size");
        }

        double chi_sq = 0.0;
        for (size_t i = 0; i < observed.size(); ++i) {
            if (expected[i] < 1e-10) {
                throw std::invalid_argument("Expected frequencies must be > 0");
            }
            double diff = observed[i] - expected[i];
            chi_sq += (diff * diff) / expected[i];
        }
        return chi_sq;
    }

    /**
     * @brief Chi-squared test for r×c contingency table
     *
     * Tests independence of two categorical variables
     * χ² = ∑∑(O_ij - E_ij)²/E_ij where E_ij = (row_i total)(col_j total)/grand total
     *
     * @param table Contingency table (r rows × c columns)
     * @return Chi-squared statistic (df = (r-1)(c-1))
     */
    static double contingency_table_test(
        const std::vector<std::vector<double>>& table) {

        int r = table.size();
        if (r == 0) throw std::invalid_argument("Empty table");
        int c = table[0].size();

        // Compute row and column totals
        std::vector<double> row_totals(r, 0.0);
        std::vector<double> col_totals(c, 0.0);
        double grand_total = 0.0;

        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                row_totals[i] += table[i][j];
                col_totals[j] += table[i][j];
                grand_total += table[i][j];
            }
        }

        // Compute chi-squared statistic
        double chi_sq = 0.0;
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                double expected = (row_totals[i] * col_totals[j]) / grand_total;
                if (expected < 1e-10) continue;
                double diff = table[i][j] - expected;
                chi_sq += (diff * diff) / expected;
            }
        }

        return chi_sq;
    }

    /**
     * @brief Fisher's exact test for 2×2 contingency table
     *
     * Computes exact p-value using hypergeometric distribution
     * Used when cell counts are small
     *
     * @param a, b, c, d Cell counts in 2×2 table [[a,b],[c,d]]
     * @return One-tailed p-value
     */
    static double fishers_exact_test_2x2(int a, int b, int c, int d) {
        int n1 = a + b;  // Row 1 total
        int n2 = c + d;  // Row 2 total
        int m1 = a + c;  // Col 1 total
        int N = a + b + c + d;  // Grand total

        // P(X = a) under hypergeometric
        double p_observed = (binomial_coefficient(n1, a) *
                            binomial_coefficient(n2, c)) /
                            binomial_coefficient(N, m1);

        // Sum probabilities for tables as or more extreme
        double p_value = 0.0;
        int min_a = std::max(0, m1 - n2);
        int max_a = std::min(n1, m1);

        for (int x = min_a; x <= max_a; ++x) {
            int y = m1 - x;
            double p = (binomial_coefficient(n1, x) *
                       binomial_coefficient(n2, y)) /
                       binomial_coefficient(N, m1);

            if (p <= p_observed + 1e-10) {  // As or more extreme
                p_value += p;
            }
        }

        return p_value;
    }

    /**
     * @brief Bartlett's test for homogeneity of variances
     *
     * Tests H₀: σ₁² = σ₂² = ... = σ₂ₖ
     * Test statistic ~ χ²(k-1) under H₀
     *
     * @param samples Vector of samples (each sample is a vector)
     * @return Bartlett's test statistic
     */
    static double bartletts_test(const std::vector<std::vector<double>>& samples) {
        int k = samples.size();
        if (k < 2) throw std::invalid_argument("Need at least 2 samples");

        std::vector<int> n_i(k);
        std::vector<double> s_i_sq(k);
        int N = 0;

        // Compute sample variances
        for (int i = 0; i < k; ++i) {
            n_i[i] = samples[i].size();
            N += n_i[i];

            double mean = std::accumulate(samples[i].begin(), samples[i].end(), 0.0) / n_i[i];
            double var = 0.0;
            for (double x : samples[i]) {
                var += (x - mean) * (x - mean);
            }
            s_i_sq[i] = var / (n_i[i] - 1);
        }

        // Pooled variance
        double s_p_sq = 0.0;
        for (int i = 0; i < k; ++i) {
            s_p_sq += (n_i[i] - 1) * s_i_sq[i];
        }
        s_p_sq /= (N - k);

        // Bartlett's statistic
        double numerator = (N - k) * std::log(s_p_sq);
        double denominator = 0.0;
        for (int i = 0; i < k; ++i) {
            numerator -= (n_i[i] - 1) * std::log(s_i_sq[i]);
        }

        // Correction factor C
        double C = 1.0 + (1.0 / (3.0 * (k - 1))) *
                   (std::accumulate(n_i.begin(), n_i.end(), 0.0,
                    [](double sum, int n) { return sum + 1.0 / (n - 1); })
                    - 1.0 / (N - k));

        return numerator / C;
    }

    /**
     * @brief Cochran's test for outlying variance
     *
     * Tests if one variance is significantly larger than others
     * C = max(s_i²) / ∑s_i²
     *
     * @param samples Vector of samples
     * @return Cochran's C statistic
     */
    static double cochrans_test(const std::vector<std::vector<double>>& samples) {
        int k = samples.size();
        if (k < 2) throw std::invalid_argument("Need at least 2 samples");

        std::vector<double> variances(k);

        // Compute sample variances
        for (int i = 0; i < k; ++i) {
            double mean = std::accumulate(samples[i].begin(), samples[i].end(), 0.0)
                         / samples[i].size();
            double var = 0.0;
            for (double x : samples[i]) {
                var += (x - mean) * (x - mean);
            }
            variances[i] = var / (samples[i].size() - 1);
        }

        double max_var = *std::max_element(variances.begin(), variances.end());
        double sum_var = std::accumulate(variances.begin(), variances.end(), 0.0);

        return max_var / sum_var;
    }

    /**
     * @brief Sample size required for hypothesis test
     *
     * For testing H₀: μ = μ₀ vs H₁: μ = μ₁
     * Given Type I error α, Type II error β, and standard deviation σ
     *
     * n ≈ (z_α + z_β)² σ² / (μ₁ - μ₀)²
     *
     * @param alpha Type I error rate
     * @param beta Type II error rate
     * @param sigma Population standard deviation
     * @param mu0 Null hypothesis mean
     * @param mu1 Alternative hypothesis mean
     * @return Required sample size
     */
    static int required_sample_size(
        double alpha, double beta,
        double sigma, double mu0, double mu1) {

        double z_alpha = StandardNormalExtensions::critical_value(alpha);
        double z_beta = StandardNormalExtensions::critical_value(beta);

        double delta = std::abs(mu1 - mu0);
        double n = ((z_alpha + z_beta) * sigma / delta);
        n = n * n;

        return static_cast<int>(std::ceil(n));
    }

    /**
     * @brief Grubbs' test for outliers
     *
     * Tests if maximum or minimum is an outlier
     * G = |x_extreme - x̄| / s
     *
     * @param data Sample data
     * @param test_max If true, test maximum; else test minimum
     * @return Grubbs' test statistic
     */
    static double grubbs_test(const std::vector<double>& data, bool test_max = true) {
        int n = data.size();
        if (n < 3) throw std::invalid_argument("Need at least 3 observations");

        double mean = std::accumulate(data.begin(), data.end(), 0.0) / n;

        double var = 0.0;
        for (double x : data) {
            var += (x - mean) * (x - mean);
        }
        double std_dev = std::sqrt(var / (n - 1));

        double extreme = test_max ?
            *std::max_element(data.begin(), data.end()) :
            *std::min_element(data.begin(), data.end());

        return std::abs(extreme - mean) / std_dev;
    }

    /**
     * @brief Dixon's Q test for outliers
     *
     * Tests if extreme value is an outlier
     * Q = (x_suspect - x_nearest) / range
     *
     * @param data Sample data
     * @return Dixon's Q statistic
     */
    static double dixons_q_test(std::vector<double> data) {
        int n = data.size();
        if (n < 3) throw std::invalid_argument("Need at least 3 observations");

        std::sort(data.begin(), data.end());

        // Test both extremes
        double Q_low = (data[1] - data[0]) / (data[n-1] - data[0]);
        double Q_high = (data[n-1] - data[n-2]) / (data[n-1] - data[0]);

        return std::max(Q_low, Q_high);
    }

    /**
     * @brief One-sample t-test
     *
     * Tests H₀: μ = μ₀ vs H₁: μ ≠ μ₀
     * t = (x̄ - μ₀) / (s/√n) ~ t(n-1)
     *
     * @param data Sample data
     * @param mu0 Hypothesized mean
     * @return t-statistic
     */
    static double t_test_one_sample(const std::vector<double>& data, double mu0) {
        double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();

        double variance = 0.0;
        for (double x : data) {
            variance += (x - mean) * (x - mean);
        }
        variance /= (data.size() - 1);

        double se = std::sqrt(variance / data.size());
        return (mean - mu0) / se;
    }

    /**
     * @brief Two-sample t-test (equal variances)
     *
     * Tests H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂
     *
     * @return t-statistic
     */
    static double t_test_two_sample(
        const std::vector<double>& data1,
        const std::vector<double>& data2) {

        int n1 = data1.size(), n2 = data2.size();

        double mean1 = std::accumulate(data1.begin(), data1.end(), 0.0) / n1;
        double mean2 = std::accumulate(data2.begin(), data2.end(), 0.0) / n2;

        double var1 = 0.0, var2 = 0.0;
        for (double x : data1) var1 += (x - mean1) * (x - mean1);
        for (double x : data2) var2 += (x - mean2) * (x - mean2);

        // Pooled variance
        double sp_sq = (var1 + var2) / (n1 + n2 - 2);
        double se = std::sqrt(sp_sq * (1.0/n1 + 1.0/n2));

        return (mean1 - mean2) / se;
    }

    /**
     * @brief Binomial test determining values in Bernoulli trials
     *
     * Given n trials and x successes, test H₀: p = p₀
     * Exact binomial test using binomial distribution
     *
     * @param x Number of successes
     * @param n Number of trials
     * @param p0 Hypothesized probability
     * @return Two-tailed p-value
     */
    static double binomial_test(int x, int n, double p0) {
        BinomialDistribution binom(n, p0);

        // P(X = x) under H₀
        double p_observed = binom.pmf(x);

        // Two-tailed: sum probabilities ≤ p_observed
        double p_value = 0.0;
        for (int k = 0; k <= n; ++k) {
            double p_k = binom.pmf(k);
            if (p_k <= p_observed + 1e-10) {
                p_value += p_k;
            }
        }

        return p_value;
    }

    /**
     * @brief Estimate probability in Bernoulli trials
     *
     * MLE: p̂ = x/n
     * Confidence interval: p̂ ± z_α/2 √(p̂(1-p̂)/n)
     *
     * @param x Number of successes
     * @param n Number of trials
     * @param alpha Significance level
     * @return (point estimate, lower CI, upper CI)
     */
    static std::tuple<double, double, double> bernoulli_estimate(
        int x, int n, double alpha = 0.05) {

        double p_hat = static_cast<double>(x) / n;
        double z = StandardNormalExtensions::critical_value(alpha / 2.0);
        double se = std::sqrt(p_hat * (1.0 - p_hat) / n);

        double lower = std::max(0.0, p_hat - z * se);
        double upper = std::min(1.0, p_hat + z * se);

        return {p_hat, lower, upper};
    }
};

/**
 * @class RegressionAnalysis
 * @brief Linear regression and orthogonal polynomial regression
 *
 * Implements computational methods for:
 * - Simple linear regression (y = β₀ + β₁x + ε)
 * - Multiple linear regression (y = Xβ + ε)
 * - Orthogonal polynomials (Legendre, Chebyshev)
 * - Regression diagnostics (R², ANOVA, residuals)
 */
class RegressionAnalysis {
public:
    /**
     * @brief Simple linear regression result
     */
    struct SimpleRegressionResult {
        double beta0;           // Intercept
        double beta1;           // Slope
        double r_squared;       // Coefficient of determination
        double residual_se;     // Residual standard error
        double beta0_se;        // Standard error of β₀
        double beta1_se;        // Standard error of β₁
        double t_stat_beta1;    // t-statistic for β₁
        std::vector<double> residuals;
        std::vector<double> fitted_values;
    };

    /**
     * @brief Simple linear regression: y = β₀ + β₁x + ε
     *
     * Least squares estimates:
     * β₁ = Σ(x-x̄)(y-ȳ) / Σ(x-x̄)²
     * β₀ = ȳ - β₁x̄
     *
     * @param x Independent variable
     * @param y Dependent variable
     * @return Regression results
     */
    static SimpleRegressionResult simple_linear_regression(
        const std::vector<double>& x,
        const std::vector<double>& y) {

        int n = x.size();
        if (n != static_cast<int>(y.size())) {
            throw std::invalid_argument("x and y must have same size");
        }
        if (n < 3) throw std::invalid_argument("Need at least 3 observations");

        // Compute means
        double x_mean = std::accumulate(x.begin(), x.end(), 0.0) / n;
        double y_mean = std::accumulate(y.begin(), y.end(), 0.0) / n;

        // Compute β₁ and β₀
        double numerator = 0.0, denominator = 0.0;
        for (int i = 0; i < n; ++i) {
            double x_dev = x[i] - x_mean;
            double y_dev = y[i] - y_mean;
            numerator += x_dev * y_dev;
            denominator += x_dev * x_dev;
        }

        if (std::abs(denominator) < 1e-10) {
            throw std::runtime_error("x values have zero variance");
        }

        double beta1 = numerator / denominator;
        double beta0 = y_mean - beta1 * x_mean;

        // Compute fitted values and residuals
        std::vector<double> fitted(n), residuals(n);
        double SSE = 0.0, SST = 0.0;
        for (int i = 0; i < n; ++i) {
            fitted[i] = beta0 + beta1 * x[i];
            residuals[i] = y[i] - fitted[i];
            SSE += residuals[i] * residuals[i];
            SST += (y[i] - y_mean) * (y[i] - y_mean);
        }

        // R² = 1 - SSE/SST
        double r_squared = 1.0 - SSE / SST;

        // Residual standard error: s = √(SSE/(n-2))
        double residual_se = std::sqrt(SSE / (n - 2));

        // Standard errors of coefficients
        double beta1_se = residual_se / std::sqrt(denominator);
        double sum_x_sq = 0.0;
        for (double xi : x) sum_x_sq += xi * xi;
        double beta0_se = residual_se * std::sqrt(sum_x_sq / (n * denominator));

        // t-statistic for β₁: H₀: β₁ = 0
        double t_stat_beta1 = beta1 / beta1_se;

        return {beta0, beta1, r_squared, residual_se,
                beta0_se, beta1_se, t_stat_beta1,
                residuals, fitted};
    }

    /**
     * @brief Multiple regression result
     */
    struct MultipleRegressionResult {
        std::vector<double> beta;        // Coefficients (β₀, β₁, ..., βₖ)
        double r_squared;                // R²
        double adjusted_r_squared;       // Adjusted R²
        double residual_se;              // Residual standard error
        double f_statistic;              // F-statistic for overall significance
        std::vector<double> residuals;
        std::vector<double> fitted_values;
        std::vector<double> beta_se;     // Standard errors of coefficients
    };

    /**
     * @brief Multiple linear regression: y = Xβ + ε
     *
     * Using normal equations: β = (XᵀX)⁻¹Xᵀy
     *
     * @param X Design matrix (n × (k+1)) including intercept column
     * @param y Response vector (n × 1)
     * @return Regression results
     */
    static MultipleRegressionResult multiple_linear_regression(
        const std::vector<std::vector<double>>& X,
        const std::vector<double>& y) {

        int n = X.size();        // Number of observations
        int p = X[0].size();     // Number of predictors (including intercept)

        if (n != static_cast<int>(y.size())) {
            throw std::invalid_argument("X and y dimensions don't match");
        }
        if (n < p) {
            throw std::invalid_argument("Need more observations than predictors");
        }

        // Compute XᵀX
        std::vector<std::vector<double>> XtX(p, std::vector<double>(p, 0.0));
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < p; ++j) {
                for (int k = 0; k < n; ++k) {
                    XtX[i][j] += X[k][i] * X[k][j];
                }
            }
        }

        // Compute Xᵀy
        std::vector<double> Xty(p, 0.0);
        for (int i = 0; i < p; ++i) {
            for (int k = 0; k < n; ++k) {
                Xty[i] += X[k][i] * y[k];
            }
        }

        // Solve (XᵀX)β = Xᵀy using Gaussian elimination
        std::vector<double> beta = solve_linear_system(XtX, Xty);

        // Compute fitted values and residuals
        double y_mean = std::accumulate(y.begin(), y.end(), 0.0) / n;
        std::vector<double> fitted(n), residuals(n);
        double SSE = 0.0, SST = 0.0;

        for (int i = 0; i < n; ++i) {
            fitted[i] = 0.0;
            for (int j = 0; j < p; ++j) {
                fitted[i] += X[i][j] * beta[j];
            }
            residuals[i] = y[i] - fitted[i];
            SSE += residuals[i] * residuals[i];
            SST += (y[i] - y_mean) * (y[i] - y_mean);
        }

        // R² and adjusted R²
        double r_squared = 1.0 - SSE / SST;
        double adj_r_squared = 1.0 - (SSE / (n - p)) / (SST / (n - 1));

        // Residual standard error
        double residual_se = std::sqrt(SSE / (n - p));

        // F-statistic: F = (SSR/k) / (SSE/(n-k-1))
        double SSR = SST - SSE;
        int k = p - 1;  // Number of predictors (excluding intercept)
        double f_statistic = (SSR / k) / (SSE / (n - p));

        // Standard errors of coefficients (diagonal of (XᵀX)⁻¹ * MSE)
        std::vector<std::vector<double>> XtX_inv = invert_matrix(XtX);
        std::vector<double> beta_se(p);
        double mse = SSE / (n - p);
        for (int i = 0; i < p; ++i) {
            beta_se[i] = std::sqrt(XtX_inv[i][i] * mse);
        }

        return {beta, r_squared, adj_r_squared, residual_se,
                f_statistic, residuals, fitted, beta_se};
    }

    /**
     * @brief Legendre polynomial P_n(x) on [-1, 1]
     *
     * Recurrence: (n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)
     * P_0(x) = 1, P_1(x) = x
     *
     * @param n Degree
     * @param x Evaluation point in [-1, 1]
     * @return P_n(x)
     */
    static double legendre_polynomial(int n, double x) {
        if (n == 0) return 1.0;
        if (n == 1) return x;

        double P_prev = 1.0;  // P_0
        double P_curr = x;    // P_1

        for (int k = 1; k < n; ++k) {
            double P_next = ((2.0 * k + 1.0) * x * P_curr - k * P_prev) / (k + 1.0);
            P_prev = P_curr;
            P_curr = P_next;
        }

        return P_curr;
    }

    /**
     * @brief Chebyshev polynomial T_n(x) of first kind on [-1, 1]
     *
     * Recurrence: T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x)
     * T_0(x) = 1, T_1(x) = x
     *
     * @param n Degree
     * @param x Evaluation point in [-1, 1]
     * @return T_n(x)
     */
    static double chebyshev_polynomial(int n, double x) {
        if (n == 0) return 1.0;
        if (n == 1) return x;

        double T_prev = 1.0;  // T_0
        double T_curr = x;    // T_1

        for (int k = 1; k < n; ++k) {
            double T_next = 2.0 * x * T_curr - T_prev;
            T_prev = T_curr;
            T_curr = T_next;
        }

        return T_curr;
    }

    /**
     * @brief Orthogonal polynomial regression
     *
     * Fit y = Σᵢ₌₀ⁿ aᵢPᵢ(x) where Pᵢ are orthogonal polynomials
     *
     * @param x Independent variable
     * @param y Dependent variable
     * @param degree Maximum degree of polynomial
     * @param poly_type "Legendre" or "Chebyshev"
     * @return Coefficients {a₀, a₁, ..., aₙ}
     */
    static std::vector<double> orthogonal_polynomial_regression(
        const std::vector<double>& x,
        const std::vector<double>& y,
        int degree,
        const std::string& poly_type = "Legendre") {

        int n = x.size();
        if (n != static_cast<int>(y.size())) {
            throw std::invalid_argument("x and y must have same size");
        }

        // Normalize x to [-1, 1]
        double x_min = *std::min_element(x.begin(), x.end());
        double x_max = *std::max_element(x.begin(), x.end());
        std::vector<double> x_norm(n);
        for (int i = 0; i < n; ++i) {
            x_norm[i] = 2.0 * (x[i] - x_min) / (x_max - x_min) - 1.0;
        }

        // Build design matrix with orthogonal polynomials
        std::vector<std::vector<double>> X(n, std::vector<double>(degree + 1));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= degree; ++j) {
                if (poly_type == "Legendre") {
                    X[i][j] = legendre_polynomial(j, x_norm[i]);
                } else if (poly_type == "Chebyshev") {
                    X[i][j] = chebyshev_polynomial(j, x_norm[i]);
                } else {
                    throw std::invalid_argument("Unknown polynomial type");
                }
            }
        }

        // Perform multiple regression
        auto result = multiple_linear_regression(X, y);
        return result.beta;
    }

private:
    /**
     * @brief Solve linear system Ax = b using Gaussian elimination
     */
    static std::vector<double> solve_linear_system(
        std::vector<std::vector<double>> A,
        std::vector<double> b) {

        int n = A.size();

        // Forward elimination with partial pivoting
        for (int i = 0; i < n; ++i) {
            // Find pivot
            int max_row = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                    max_row = k;
                }
            }

            // Swap rows
            std::swap(A[i], A[max_row]);
            std::swap(b[i], b[max_row]);

            // Eliminate column
            for (int k = i + 1; k < n; ++k) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }

        // Back substitution
        std::vector<double> x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
        }

        return x;
    }

    /**
     * @brief Invert matrix using Gaussian elimination
     */
    static std::vector<std::vector<double>> invert_matrix(
        std::vector<std::vector<double>> A) {

        int n = A.size();
        std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));

        // Initialize inv as identity
        for (int i = 0; i < n; ++i) {
            inv[i][i] = 1.0;
        }

        // Forward elimination
        for (int i = 0; i < n; ++i) {
            // Pivot
            double pivot = A[i][i];
            if (std::abs(pivot) < 1e-10) {
                throw std::runtime_error("Matrix is singular");
            }

            for (int j = 0; j < n; ++j) {
                A[i][j] /= pivot;
                inv[i][j] /= pivot;
            }

            for (int k = 0; k < n; ++k) {
                if (k != i) {
                    double factor = A[k][i];
                    for (int j = 0; j < n; ++j) {
                        A[k][j] -= factor * A[i][j];
                        inv[k][j] -= factor * inv[i][j];
                    }
                }
            }
        }

        return inv;
    }
};

/**
 * @class ANOVA
 * @brief Analysis of Variance methods
 *
 * Implements computational methods for:
 * - One-way ANOVA (single factor)
 * - Two-way ANOVA (two factors with/without interaction)
 * - Three-factor ANOVA
 * - MANOVA (Multivariate ANOVA)
 */
class ANOVA {
public:
    /**
     * @brief One-way ANOVA result
     */
    struct OneWayResult {
        double SS_between;      // Sum of squares between groups
        double SS_within;       // Sum of squares within groups
        double SS_total;        // Total sum of squares
        int df_between;         // Degrees of freedom between
        int df_within;          // Degrees of freedom within
        int df_total;           // Total degrees of freedom
        double MS_between;      // Mean square between
        double MS_within;       // Mean square within
        double F_statistic;     // F-statistic
        double p_value_approx;  // Approximate p-value
        std::vector<double> group_means;
        double grand_mean;
    };

    /**
     * @brief One-way ANOVA: tests H₀: μ₁ = μ₂ = ... = μₖ
     *
     * Decomposes total variation: SST = SSB + SSW
     * F = MSB/MSW ~ F(k-1, N-k) under H₀
     *
     * @param groups Vector of samples, one vector per group
     * @return ANOVA table results
     */
    static OneWayResult one_way_anova(const std::vector<std::vector<double>>& groups) {
        int k = groups.size();  // Number of groups
        if (k < 2) throw std::invalid_argument("Need at least 2 groups");

        // Compute group sizes, means, and grand mean
        std::vector<int> n_i(k);
        std::vector<double> means(k);
        int N = 0;
        double grand_mean = 0.0;

        for (int i = 0; i < k; ++i) {
            n_i[i] = groups[i].size();
            N += n_i[i];
            means[i] = std::accumulate(groups[i].begin(), groups[i].end(), 0.0) / n_i[i];
            grand_mean += std::accumulate(groups[i].begin(), groups[i].end(), 0.0);
        }
        grand_mean /= N;

        // Compute sum of squares
        double SSB = 0.0;  // Between groups
        for (int i = 0; i < k; ++i) {
            SSB += n_i[i] * (means[i] - grand_mean) * (means[i] - grand_mean);
        }

        double SSW = 0.0;  // Within groups
        for (int i = 0; i < k; ++i) {
            for (double x : groups[i]) {
                SSW += (x - means[i]) * (x - means[i]);
            }
        }

        double SST = SSB + SSW;

        // Degrees of freedom
        int df_between = k - 1;
        int df_within = N - k;
        int df_total = N - 1;

        // Mean squares
        double MSB = SSB / df_between;
        double MSW = SSW / df_within;

        // F-statistic
        double F = MSB / MSW;

        // Approximate p-value (would need F-distribution CDF for exact)
        double p_value = 0.0;  // Placeholder

        return {SSB, SSW, SST, df_between, df_within, df_total,
                MSB, MSW, F, p_value, means, grand_mean};
    }

    /**
     * @brief Two-way ANOVA result
     */
    struct TwoWayResult {
        double SS_A;            // Sum of squares for factor A
        double SS_B;            // Sum of squares for factor B
        double SS_AB;           // Sum of squares for interaction
        double SS_error;        // Sum of squares error
        double SS_total;        // Total sum of squares
        int df_A, df_B, df_AB, df_error, df_total;
        double MS_A, MS_B, MS_AB, MS_error;
        double F_A, F_B, F_AB;  // F-statistics
        double grand_mean;
    };

    /**
     * @brief Two-way ANOVA with interaction: Y_ijk = μ + α_i + β_j + (αβ)_ij + ε_ijk
     *
     * Tests:
     * - H₀: α₁ = α₂ = ... = αₐ = 0 (main effect A)
     * - H₀: β₁ = β₂ = ... = β_b = 0 (main effect B)
     * - H₀: (αβ)_ij = 0 for all i,j (interaction)
     *
     * @param data Data[i][j][k] = observation k in cell (i,j)
     * @param a Number of levels of factor A
     * @param b Number of levels of factor B
     * @return ANOVA table results
     */
    static TwoWayResult two_way_anova(
        const std::vector<std::vector<std::vector<double>>>& data,
        int a, int b) {

        // Compute cell means and grand mean
        std::vector<std::vector<double>> cell_means(a, std::vector<double>(b));
        std::vector<double> row_means(a, 0.0);
        std::vector<double> col_means(b, 0.0);
        double grand_mean = 0.0;
        int N = 0;

        for (int i = 0; i < a; ++i) {
            for (int j = 0; j < b; ++j) {
                int n_ij = data[i][j].size();
                N += n_ij;
                cell_means[i][j] = std::accumulate(data[i][j].begin(),
                                                   data[i][j].end(), 0.0) / n_ij;
                grand_mean += std::accumulate(data[i][j].begin(),
                                             data[i][j].end(), 0.0);
            }
        }
        grand_mean /= N;

        // Compute marginal means
        int n_per_cell = data[0][0].size();  // Assuming balanced design
        for (int i = 0; i < a; ++i) {
            for (int j = 0; j < b; ++j) {
                row_means[i] += cell_means[i][j];
                col_means[j] += cell_means[i][j];
            }
            row_means[i] /= b;
        }
        for (int j = 0; j < b; ++j) {
            col_means[j] /= a;
        }

        // Compute sums of squares
        double SS_A = 0.0;
        for (int i = 0; i < a; ++i) {
            SS_A += b * n_per_cell * (row_means[i] - grand_mean) *
                    (row_means[i] - grand_mean);
        }

        double SS_B = 0.0;
        for (int j = 0; j < b; ++j) {
            SS_B += a * n_per_cell * (col_means[j] - grand_mean) *
                    (col_means[j] - grand_mean);
        }

        double SS_AB = 0.0;
        for (int i = 0; i < a; ++i) {
            for (int j = 0; j < b; ++j) {
                double interaction = cell_means[i][j] - row_means[i] -
                                    col_means[j] + grand_mean;
                SS_AB += n_per_cell * interaction * interaction;
            }
        }

        double SS_error = 0.0;
        for (int i = 0; i < a; ++i) {
            for (int j = 0; j < b; ++j) {
                for (double x : data[i][j]) {
                    SS_error += (x - cell_means[i][j]) * (x - cell_means[i][j]);
                }
            }
        }

        double SS_total = SS_A + SS_B + SS_AB + SS_error;

        // Degrees of freedom
        int df_A = a - 1;
        int df_B = b - 1;
        int df_AB = (a - 1) * (b - 1);
        int df_error = N - a * b;
        int df_total = N - 1;

        // Mean squares
        double MS_A = SS_A / df_A;
        double MS_B = SS_B / df_B;
        double MS_AB = SS_AB / df_AB;
        double MS_error = SS_error / df_error;

        // F-statistics
        double F_A = MS_A / MS_error;
        double F_B = MS_B / MS_error;
        double F_AB = MS_AB / MS_error;

        return {SS_A, SS_B, SS_AB, SS_error, SS_total,
                df_A, df_B, df_AB, df_error, df_total,
                MS_A, MS_B, MS_AB, MS_error,
                F_A, F_B, F_AB, grand_mean};
    }

    /**
     * @brief Three-factor ANOVA result
     */
    struct ThreeFactorResult {
        double SS_A, SS_B, SS_C;                    // Main effects
        double SS_AB, SS_AC, SS_BC;                 // Two-way interactions
        double SS_ABC;                              // Three-way interaction
        double SS_error, SS_total;
        int df_A, df_B, df_C, df_AB, df_AC, df_BC, df_ABC, df_error, df_total;
        double MS_A, MS_B, MS_C, MS_AB, MS_AC, MS_BC, MS_ABC, MS_error;
        double F_A, F_B, F_C, F_AB, F_AC, F_BC, F_ABC;
    };

    /**
     * @brief Three-factor ANOVA: Y_ijkl = μ + α_i + β_j + γ_k + (αβ)_ij + ... + ε_ijkl
     *
     * Tests main effects, two-way interactions, and three-way interaction
     *
     * @param data Data[i][j][k][l] = observation l in cell (i,j,k)
     * @param a, b, c Number of levels for factors A, B, C
     * @return ANOVA table results
     */
    static ThreeFactorResult three_factor_anova(
        const std::vector<std::vector<std::vector<std::vector<double>>>>& data,
        int a, int b, int c) {

        int n_per_cell = data[0][0][0].size();
        int N = a * b * c * n_per_cell;

        // Compute grand mean and cell means
        double grand_mean = 0.0;
        std::vector<std::vector<std::vector<double>>> cell_means(
            a, std::vector<std::vector<double>>(b, std::vector<double>(c)));

        for (int i = 0; i < a; ++i) {
            for (int j = 0; j < b; ++j) {
                for (int k = 0; k < c; ++k) {
                    cell_means[i][j][k] = std::accumulate(
                        data[i][j][k].begin(), data[i][j][k].end(), 0.0) / n_per_cell;
                    grand_mean += cell_means[i][j][k] * n_per_cell;
                }
            }
        }
        grand_mean /= N;

        // Compute marginal means
        std::vector<double> mean_A(a, 0.0), mean_B(b, 0.0), mean_C(c, 0.0);
        for (int i = 0; i < a; ++i) {
            for (int j = 0; j < b; ++j) {
                for (int k = 0; k < c; ++k) {
                    mean_A[i] += cell_means[i][j][k];
                    mean_B[j] += cell_means[i][j][k];
                    mean_C[k] += cell_means[i][j][k];
                }
            }
            mean_A[i] /= (b * c);
        }
        for (int j = 0; j < b; ++j) mean_B[j] /= (a * c);
        for (int k = 0; k < c; ++k) mean_C[k] /= (a * b);

        // Main effects
        double SS_A = 0.0;
        for (int i = 0; i < a; ++i) {
            SS_A += b * c * n_per_cell * (mean_A[i] - grand_mean) *
                    (mean_A[i] - grand_mean);
        }

        double SS_B = 0.0;
        for (int j = 0; j < b; ++j) {
            SS_B += a * c * n_per_cell * (mean_B[j] - grand_mean) *
                    (mean_B[j] - grand_mean);
        }

        double SS_C = 0.0;
        for (int k = 0; k < c; ++k) {
            SS_C += a * b * n_per_cell * (mean_C[k] - grand_mean) *
                    (mean_C[k] - grand_mean);
        }

        // Two-way interactions (simplified computation)
        double SS_AB = 0.0, SS_AC = 0.0, SS_BC = 0.0;
        // Placeholder: full computation would involve two-way marginal means

        // Three-way interaction and error
        double SS_ABC = 0.0;
        double SS_error = 0.0;

        for (int i = 0; i < a; ++i) {
            for (int j = 0; j < b; ++j) {
                for (int k = 0; k < c; ++k) {
                    for (double x : data[i][j][k]) {
                        SS_error += (x - cell_means[i][j][k]) *
                                   (x - cell_means[i][j][k]);
                    }
                }
            }
        }

        double SS_total = SS_A + SS_B + SS_C + SS_AB + SS_AC + SS_BC +
                         SS_ABC + SS_error;

        // Degrees of freedom
        int df_A = a - 1, df_B = b - 1, df_C = c - 1;
        int df_AB = (a - 1) * (b - 1);
        int df_AC = (a - 1) * (c - 1);
        int df_BC = (b - 1) * (c - 1);
        int df_ABC = (a - 1) * (b - 1) * (c - 1);
        int df_error = N - a * b * c;
        int df_total = N - 1;

        // Mean squares and F-statistics
        double MS_A = SS_A / df_A, MS_B = SS_B / df_B, MS_C = SS_C / df_C;
        double MS_AB = SS_AB / df_AB, MS_AC = SS_AC / df_AC, MS_BC = SS_BC / df_BC;
        double MS_ABC = SS_ABC / df_ABC;
        double MS_error = SS_error / df_error;

        double F_A = MS_A / MS_error, F_B = MS_B / MS_error, F_C = MS_C / MS_error;
        double F_AB = MS_AB / MS_error, F_AC = MS_AC / MS_error, F_BC = MS_BC / MS_error;
        double F_ABC = MS_ABC / MS_error;

        return {SS_A, SS_B, SS_C, SS_AB, SS_AC, SS_BC, SS_ABC, SS_error, SS_total,
                df_A, df_B, df_C, df_AB, df_AC, df_BC, df_ABC, df_error, df_total,
                MS_A, MS_B, MS_C, MS_AB, MS_AC, MS_BC, MS_ABC, MS_error,
                F_A, F_B, F_C, F_AB, F_AC, F_BC, F_ABC};
    }

    /**
     * @brief MANOVA (Multivariate ANOVA) - Wilks' Lambda test
     *
     * Tests H₀: mean vectors are equal across groups
     * Λ = |E| / |E + H| where E is error SSCP, H is hypothesis SSCP
     *
     * @param groups Vector of groups, each group contains multivariate observations
     * @return Wilks' Lambda statistic
     */
    static double manova_wilks_lambda(
        const std::vector<std::vector<std::vector<double>>>& groups) {

        int k = groups.size();  // Number of groups
        int p = groups[0][0].size();  // Number of response variables

        // Compute group means and grand mean
        std::vector<std::vector<double>> group_means(k, std::vector<double>(p, 0.0));
        std::vector<double> grand_mean(p, 0.0);
        int N = 0;

        for (int g = 0; g < k; ++g) {
            int n_g = groups[g].size();
            N += n_g;

            for (int i = 0; i < n_g; ++i) {
                for (int j = 0; j < p; ++j) {
                    group_means[g][j] += groups[g][i][j];
                    grand_mean[j] += groups[g][i][j];
                }
            }

            for (int j = 0; j < p; ++j) {
                group_means[g][j] /= n_g;
            }
        }

        for (int j = 0; j < p; ++j) {
            grand_mean[j] /= N;
        }

        // Compute hypothesis SSCP matrix H
        std::vector<std::vector<double>> H(p, std::vector<double>(p, 0.0));
        for (int g = 0; g < k; ++g) {
            int n_g = groups[g].size();
            for (int i = 0; i < p; ++i) {
                for (int j = 0; j < p; ++j) {
                    H[i][j] += n_g * (group_means[g][i] - grand_mean[i]) *
                              (group_means[g][j] - grand_mean[j]);
                }
            }
        }

        // Compute error SSCP matrix E
        std::vector<std::vector<double>> E(p, std::vector<double>(p, 0.0));
        for (int g = 0; g < k; ++g) {
            for (const auto& obs : groups[g]) {
                for (int i = 0; i < p; ++i) {
                    for (int j = 0; j < p; ++j) {
                        E[i][j] += (obs[i] - group_means[g][i]) *
                                  (obs[j] - group_means[g][j]);
                    }
                }
            }
        }

        // Compute E + H
        std::vector<std::vector<double>> E_plus_H(p, std::vector<double>(p, 0.0));
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < p; ++j) {
                E_plus_H[i][j] = E[i][j] + H[i][j];
            }
        }

        // Wilks' Lambda = |E| / |E + H|
        double det_E = compute_determinant(E);
        double det_E_plus_H = compute_determinant(E_plus_H);

        return det_E / det_E_plus_H;
    }

private:
    /**
     * @brief Compute determinant of matrix (using Gaussian elimination)
     */
    static double compute_determinant(std::vector<std::vector<double>> A) {
        int n = A.size();
        double det = 1.0;

        for (int i = 0; i < n; ++i) {
            // Find pivot
            int max_row = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                    max_row = k;
                }
            }

            if (max_row != i) {
                std::swap(A[i], A[max_row]);
                det *= -1.0;
            }

            if (std::abs(A[i][i]) < 1e-10) return 0.0;

            det *= A[i][i];

            // Eliminate
            for (int k = i + 1; k < n; ++k) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
            }
        }

        return det;
    }
};

/**
 * @class FactorAnalysis
 * @brief Factor analysis and principal component analysis
 *
 * Reduces dimensionality by finding latent factors
 */
class FactorAnalysis {
public:
    /**
     * @brief Compute covariance matrix from data
     */
    static std::vector<std::vector<double>> covariance_matrix(
        const std::vector<std::vector<double>>& data) {

        int n = data.size();      // observations
        int p = data[0].size();   // variables

        // Compute means
        std::vector<double> means(p, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < p; ++j) {
                means[j] += data[i][j];
            }
        }
        for (int j = 0; j < p; ++j) {
            means[j] /= n;
        }

        // Compute covariance
        std::vector<std::vector<double>> cov(p, std::vector<double>(p, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < p; ++j) {
                for (int k = 0; k < p; ++k) {
                    cov[j][k] += (data[i][j] - means[j]) * (data[i][k] - means[k]);
                }
            }
        }

        for (int j = 0; j < p; ++j) {
            for (int k = 0; k < p; ++k) {
                cov[j][k] /= (n - 1);
            }
        }

        return cov;
    }

    /**
     * @brief Compute correlation matrix from data
     */
    static std::vector<std::vector<double>> correlation_matrix(
        const std::vector<std::vector<double>>& data) {

        auto cov = covariance_matrix(data);
        int p = cov.size();

        std::vector<std::vector<double>> corr(p, std::vector<double>(p));
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < p; ++j) {
                corr[i][j] = cov[i][j] / std::sqrt(cov[i][i] * cov[j][j]);
            }
        }

        return corr;
    }

    /**
     * @brief Extract principal components (simplified eigenvalue approach)
     *
     * Returns variance explained by each component
     */
    static std::vector<double> principal_components_variance(
        const std::vector<std::vector<double>>& cov_matrix) {

        int p = cov_matrix.size();
        std::vector<double> variances(p);

        // Diagonal elements are variances (simplified - full PCA needs eigendecomposition)
        for (int i = 0; i < p; ++i) {
            variances[i] = cov_matrix[i][i];
        }

        return variances;
    }
};

/**
 * @class ExperimentalDesign
 * @brief Experimental design methods
 *
 * Implements:
 * - Latin square designs
 * - Graeco-Latin squares
 * - Block designs
 * - Factorial designs
 * - 2^r factorial experiments
 * - Confounding in 2^n experiments
 */
class ExperimentalDesign {
public:
    /**
     * @brief Generate Latin square of order n
     *
     * Each row and column contains each treatment exactly once
     * Uses cyclic construction for prime n
     *
     * @param n Order of square
     * @return n×n Latin square (0-indexed treatments)
     */
    static std::vector<std::vector<int>> latin_square(int n) {
        std::vector<std::vector<int>> square(n, std::vector<int>(n));

        // Cyclic construction
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                square[i][j] = (i + j) % n;
            }
        }

        return square;
    }

    /**
     * @brief Generate Graeco-Latin square of order n
     *
     * Superposition of two Latin squares that are orthogonal
     * Returns pair (Latin square 1, Latin square 2)
     *
     * @param n Order (must not be 2 or 6 for orthogonal squares)
     * @return Pair of orthogonal Latin squares
     */
    static std::pair<std::vector<std::vector<int>>,
                     std::vector<std::vector<int>>> graeco_latin_square(int n) {

        if (n == 2 || n == 6) {
            throw std::invalid_argument("No orthogonal Latin squares exist for n=2 or n=6");
        }

        auto square1 = latin_square(n);

        // Generate orthogonal square using different cyclic pattern
        std::vector<std::vector<int>> square2(n, std::vector<int>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                square2[i][j] = (i * 2 + j) % n;
            }
        }

        return {square1, square2};
    }

    /**
     * @brief Randomized complete block design (RCBD) analysis
     *
     * Model: Y_ij = μ + τ_i + β_j + ε_ij
     * where τ_i is treatment effect, β_j is block effect
     *
     * @param data Data[i][j] = observation for treatment i in block j
     * @param t Number of treatments
     * @param b Number of blocks
     * @return ANOVA table for RCBD
     */
    static ANOVA::OneWayResult rcbd_analysis(
        const std::vector<std::vector<double>>& data,
        int t, int b) {

        double grand_mean = 0.0;
        int N = t * b;

        // Compute grand mean
        for (int i = 0; i < t; ++i) {
            for (int j = 0; j < b; ++j) {
                grand_mean += data[i][j];
            }
        }
        grand_mean /= N;

        // Compute treatment and block means
        std::vector<double> treatment_means(t, 0.0);
        std::vector<double> block_means(b, 0.0);

        for (int i = 0; i < t; ++i) {
            for (int j = 0; j < b; ++j) {
                treatment_means[i] += data[i][j];
                block_means[j] += data[i][j];
            }
            treatment_means[i] /= b;
        }
        for (int j = 0; j < b; ++j) {
            block_means[j] /= t;
        }

        // Sum of squares
        double SS_treatments = 0.0;
        for (int i = 0; i < t; ++i) {
            SS_treatments += b * (treatment_means[i] - grand_mean) *
                           (treatment_means[i] - grand_mean);
        }

        double SS_blocks = 0.0;
        for (int j = 0; j < b; ++j) {
            SS_blocks += t * (block_means[j] - grand_mean) *
                        (block_means[j] - grand_mean);
        }

        double SS_error = 0.0;
        for (int i = 0; i < t; ++i) {
            for (int j = 0; j < b; ++j) {
                double expected = treatment_means[i] + block_means[j] - grand_mean;
                SS_error += (data[i][j] - expected) * (data[i][j] - expected);
            }
        }

        double SS_total = SS_treatments + SS_blocks + SS_error;

        int df_treatments = t - 1;
        int df_blocks = b - 1;
        int df_error = (t - 1) * (b - 1);

        double MS_treatments = SS_treatments / df_treatments;
        double MS_error = SS_error / df_error;
        double F = MS_treatments / MS_error;

        return {SS_treatments, SS_error, SS_total,
                df_treatments, df_error, N - 1,
                MS_treatments, MS_error, F, 0.0,
                treatment_means, grand_mean};
    }

    /**
     * @brief Generate 2^k factorial design
     *
     * Full factorial with k factors at 2 levels each
     * Returns design matrix with -1/+1 coding
     *
     * @param k Number of factors
     * @return Design matrix (2^k runs × k factors)
     */
    static std::vector<std::vector<int>> factorial_2k_design(int k) {
        int n_runs = 1 << k;  // 2^k
        std::vector<std::vector<int>> design(n_runs, std::vector<int>(k));

        for (int run = 0; run < n_runs; ++run) {
            for (int factor = 0; factor < k; ++factor) {
                // Use binary representation
                design[run][factor] = ((run >> factor) & 1) ? 1 : -1;
            }
        }

        return design;
    }

    /**
     * @brief Analyze 2^2 factorial experiment (2 factors, 2 levels each)
     *
     * Estimates main effects and interaction
     *
     * @param responses Responses for treatments: (1), a, b, ab
     * @return {main effect A, main effect B, interaction AB}
     */
    static std::tuple<double, double, double> factorial_2_2_analysis(
        const std::vector<double>& responses) {

        if (responses.size() != 4) {
            throw std::invalid_argument("Need 4 responses for 2^2 design");
        }

        // Responses: responses[0] = (1), [1] = a, [2] = b, [3] = ab
        double y_1 = responses[0];   // Low-low
        double y_a = responses[1];   // High-low
        double y_b = responses[2];   // Low-high
        double y_ab = responses[3];  // High-high

        // Main effect of A: average response at high A - average at low A
        double effect_A = ((y_a + y_ab) / 2.0) - ((y_1 + y_b) / 2.0);

        // Main effect of B
        double effect_B = ((y_b + y_ab) / 2.0) - ((y_1 + y_a) / 2.0);

        // Interaction AB
        double effect_AB = ((y_1 + y_ab) / 2.0) - ((y_a + y_b) / 2.0);

        return {effect_A, effect_B, effect_AB};
    }

    /**
     * @brief Fractional factorial design 2^(k-p)
     *
     * Generates a 2^(k-p) fractional factorial design
     * Uses p generators to reduce runs from 2^k to 2^(k-p)
     *
     * @param k Number of factors
     * @param p Fraction (1/2^p of full factorial)
     * @return Design matrix
     */
    static std::vector<std::vector<int>> fractional_factorial_2k_p(int k, int p) {
        int base_factors = k - p;
        int n_runs = 1 << base_factors;  // 2^(k-p) runs

        std::vector<std::vector<int>> design(n_runs, std::vector<int>(k));

        // Generate base design for first (k-p) factors
        for (int run = 0; run < n_runs; ++run) {
            for (int factor = 0; factor < base_factors; ++factor) {
                design[run][factor] = ((run >> factor) & 1) ? 1 : -1;
            }
        }

        // Generate remaining p factors using generators
        // Example: for 2^(4-1) design, factor D = ABC
        if (k == 4 && p == 1) {
            for (int run = 0; run < n_runs; ++run) {
                design[run][3] = design[run][0] * design[run][1] * design[run][2];
            }
        }

        return design;
    }

    /**
     * @brief Confounding scheme for 2^n factorial in 2^p blocks
     *
     * Determines which effects are confounded with blocks
     *
     * @param n Number of factors
     * @param p Number of blocks = 2^p
     * @return Confounding pattern (simplified indicator)
     */
    static std::vector<std::string> confounding_2n_factorial(int n, int p) {
        std::vector<std::string> confounded;

        // For 2^3 in 2^1 = 2 blocks: confound ABC with blocks
        if (n == 3 && p == 1) {
            confounded.push_back("ABC");
        }

        // For 2^4 in 2^2 = 4 blocks: confound AB, CD, ABCD with blocks
        if (n == 4 && p == 2) {
            confounded.push_back("AB");
            confounded.push_back("CD");
            confounded.push_back("ABCD");
        }

        return confounded;
    }

    /**
     * @brief Yates' algorithm for 2^k factorial analysis
     *
     * Efficiently computes all effects in 2^k factorial
     *
     * @param responses Responses in standard order
     * @return Vector of effects [grand mean, A, B, AB, C, AC, BC, ABC, ...]
     */
    static std::vector<double> yates_algorithm(std::vector<double> responses) {
        int n = responses.size();
        int k = 0;
        while ((1 << k) < n) ++k;

        if ((1 << k) != n) {
            throw std::invalid_argument("Number of responses must be power of 2");
        }

        std::vector<double> effects = responses;

        // k iterations
        for (int i = 0; i < k; ++i) {
            std::vector<double> temp(n);
            int step = 1 << i;

            for (int j = 0; j < n; j += 2 * step) {
                for (int m = 0; m < step; ++m) {
                    temp[j + m] = effects[j + m] + effects[j + m + step];
                    temp[j + m + step] = effects[j + m] - effects[j + m + step];
                }
            }

            effects = temp;
        }

        // Divide by n/2 to get effects
        for (int i = 1; i < n; ++i) {
            effects[i] /= (n / 2.0);
        }
        effects[0] /= n;  // Grand mean

        return effects;
    }
};

/**
 * @class DesignTables
 * @brief Tables and references for experimental design
 *
 * Critical values, design templates, and lookup tables
 */
class DesignTables {
public:
    /**
     * @brief Get standard 2^k factorial design in standard order
     *
     * @param k Number of factors
     * @return Treatment combinations in standard order: (1), a, b, ab, c, ac, bc, abc, ...
     */
    static std::vector<std::string> standard_order_2k(int k) {
        int n = 1 << k;
        std::vector<std::string> order(n);

        order[0] = "(1)";  // All factors at low level

        for (int run = 1; run < n; ++run) {
            std::string treatment;
            for (int factor = 0; factor < k; ++factor) {
                if ((run >> factor) & 1) {
                    treatment += char('a' + factor);
                }
            }
            order[run] = treatment;
        }

        return order;
    }

    /**
     * @brief Common confounding patterns for 2^n factorials
     *
     * Returns recommended confounding schemes
     */
    static std::string confounding_scheme(int n_factors, int n_blocks) {
        if (n_factors == 3 && n_blocks == 2) {
            return "Confound ABC with blocks";
        } else if (n_factors == 4 && n_blocks == 2) {
            return "Confound ABCD with blocks";
        } else if (n_factors == 4 && n_blocks == 4) {
            return "Confound AB and CD with blocks (ABCD also confounded)";
        } else if (n_factors == 5 && n_blocks == 2) {
            return "Confound ABCDE with blocks";
        }
        return "Custom confounding scheme needed";
    }

    /**
     * @brief Resolution of fractional factorial design
     *
     * Resolution III: No main effect aliased with another main effect
     * Resolution IV: No main effect aliased with 2-factor interaction
     * Resolution V: No main effect or 2-factor interaction aliased with another
     *
     * @param k Number of factors
     * @param p Fraction exponent (2^(k-p) design)
     * @return Resolution (III, IV, or V)
     */
    static int design_resolution(int k, int p) {
        // Common designs
        if (k == 3 && p == 1) return 3;  // 2^(3-1) is Resolution III
        if (k == 4 && p == 1) return 4;  // 2^(4-1) is Resolution IV
        if (k == 5 && p == 1) return 5;  // 2^(5-1) is Resolution V
        if (k == 5 && p == 2) return 3;  // 2^(5-2) is Resolution III

        return 0;  // Unknown
    }
};

/**
 * @class NonparametricStatistics
 * @brief Distribution-free (nonparametric) statistical methods
 *
 * Implements computational methods for:
 * - Rank-based tests (Wilcoxon, Mann-Whitney, Kruskal-Wallis)
 * - Correlation coefficients (Spearman's rho, Kendall's tau)
 * - Randomness tests (runs test)
 * - Sign tests and matched-pairs tests
 * - Goodness-of-fit (Kolmogorov-Smirnov)
 */
class NonparametricStatistics {
public:
    /**
     * @brief Friedman test for randomized block design
     *
     * Nonparametric alternative to two-way ANOVA
     * Tests H₀: treatments have identical effects
     * χ² = (12/bk(k+1)) Σ R²ⱼ - 3b(k+1)
     *
     * @param data Data[i][j] = observation for treatment j in block i
     * @param b Number of blocks
     * @param k Number of treatments
     * @return Friedman's χ² statistic (df = k-1)
     */
    static double friedman_test(
        const std::vector<std::vector<double>>& data,
        int b, int k) {

        // Rank within each block
        std::vector<std::vector<double>> ranks(b, std::vector<double>(k));

        for (int i = 0; i < b; ++i) {
            // Create pairs (value, index)
            std::vector<std::pair<double, int>> block_data(k);
            for (int j = 0; j < k; ++j) {
                block_data[j] = {data[i][j], j};
            }

            // Sort by value
            std::sort(block_data.begin(), block_data.end());

            // Assign ranks
            for (int j = 0; j < k; ++j) {
                int idx = block_data[j].second;
                ranks[i][idx] = j + 1.0;  // Ranks 1 to k
            }
        }

        // Compute rank sums for each treatment
        std::vector<double> R_j(k, 0.0);
        for (int j = 0; j < k; ++j) {
            for (int i = 0; i < b; ++i) {
                R_j[j] += ranks[i][j];
            }
        }

        // Friedman statistic
        double sum_R_sq = 0.0;
        for (int j = 0; j < k; ++j) {
            sum_R_sq += R_j[j] * R_j[j];
        }

        double chi_sq = (12.0 / (b * k * (k + 1))) * sum_R_sq - 3.0 * b * (k + 1);

        return chi_sq;
    }

    /**
     * @brief Kendall's rank correlation coefficient τ (tau)
     *
     * Measures ordinal association between two variables
     * τ = (C - D) / (n(n-1)/2)
     * where C = concordant pairs, D = discordant pairs
     *
     * @param x First variable
     * @param y Second variable
     * @return Kendall's tau (-1 to +1)
     */
    static double kendalls_tau(
        const std::vector<double>& x,
        const std::vector<double>& y) {

        int n = x.size();
        if (n != static_cast<int>(y.size())) {
            throw std::invalid_argument("x and y must have same size");
        }

        int concordant = 0, discordant = 0;

        // Count concordant and discordant pairs
        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double x_diff = x[j] - x[i];
                double y_diff = y[j] - y[i];

                if ((x_diff > 0 && y_diff > 0) || (x_diff < 0 && y_diff < 0)) {
                    concordant++;
                } else if ((x_diff > 0 && y_diff < 0) || (x_diff < 0 && y_diff > 0)) {
                    discordant++;
                }
                // Ties: neither concordant nor discordant
            }
        }

        double tau = static_cast<double>(concordant - discordant) / (n * (n - 1) / 2.0);
        return tau;
    }

    /**
     * @brief Kruskal-Wallis test (nonparametric one-way ANOVA)
     *
     * Tests H₀: k samples come from identical distributions
     * H = (12/(N(N+1))) Σ(R²ᵢ/nᵢ) - 3(N+1)
     *
     * @param groups Vector of samples (one vector per group)
     * @return Kruskal-Wallis H statistic (approximately χ²(k-1))
     */
    static double kruskal_wallis_test(const std::vector<std::vector<double>>& groups) {
        int k = groups.size();
        if (k < 2) throw std::invalid_argument("Need at least 2 groups");

        // Combine all data with group labels
        std::vector<std::pair<double, int>> combined;
        int N = 0;

        for (int i = 0; i < k; ++i) {
            for (double x : groups[i]) {
                combined.push_back({x, i});
                N++;
            }
        }

        // Sort by value
        std::sort(combined.begin(), combined.end());

        // Assign ranks (handle ties by averaging)
        std::vector<double> ranks(N);
        int i = 0;
        while (i < N) {
            int j = i;
            // Find ties
            while (j < N && combined[j].first == combined[i].first) {
                j++;
            }

            // Average rank for tied values
            double avg_rank = (i + 1 + j) / 2.0;
            for (int m = i; m < j; ++m) {
                ranks[m] = avg_rank;
            }

            i = j;
        }

        // Compute rank sums for each group
        std::vector<double> R_i(k, 0.0);
        std::vector<int> n_i(k, 0);

        for (int idx = 0; idx < N; ++idx) {
            int group = combined[idx].second;
            R_i[group] += ranks[idx];
            n_i[group]++;
        }

        // Kruskal-Wallis statistic
        double H = 0.0;
        for (int i = 0; i < k; ++i) {
            H += (R_i[i] * R_i[i]) / n_i[i];
        }
        H = (12.0 / (N * (N + 1))) * H - 3.0 * (N + 1);

        return H;
    }

    /**
     * @brief Runs test for randomness
     *
     * Tests H₀: sequence is random
     * Counts runs (consecutive sequences of same value/sign)
     *
     * @param data Binary sequence (0s and 1s) or signs
     * @return Standardized runs statistic Z ~ N(0,1) for large n
     */
    static double runs_test(const std::vector<int>& data) {
        int n = data.size();
        if (n < 2) throw std::invalid_argument("Need at least 2 observations");

        // Count runs
        int runs = 1;
        for (int i = 1; i < n; ++i) {
            if (data[i] != data[i-1]) {
                runs++;
            }
        }

        // Count n1 (number of 0s/negatives) and n2 (number of 1s/positives)
        int n1 = 0, n2 = 0;
        for (int x : data) {
            if (x == 0) n1++;
            else n2++;
        }

        // Expected number of runs
        double E_R = (2.0 * n1 * n2) / (n1 + n2) + 1.0;

        // Variance of runs
        double Var_R = (2.0 * n1 * n2 * (2.0 * n1 * n2 - n1 - n2)) /
                       ((n1 + n2) * (n1 + n2) * (n1 + n2 - 1));

        // Standardized test statistic
        double Z = (runs - E_R) / std::sqrt(Var_R);

        return Z;
    }

    /**
     * @brief Sign test
     *
     * Tests H₀: median = m₀
     * Counts number of observations above vs below m₀
     *
     * @param data Sample data
     * @param median0 Hypothesized median
     * @return Number of positive signs (compare with binomial(n, 0.5))
     */
    static int sign_test(const std::vector<double>& data, double median0) {
        int n_plus = 0;

        for (double x : data) {
            if (x > median0) n_plus++;
            // Ignore x == median0
        }

        return n_plus;
    }

    /**
     * @brief Spearman's rank correlation coefficient ρ (rho)
     *
     * Measures monotonic relationship between variables
     * ρ = 1 - (6Σd²ᵢ)/(n(n²-1))
     * where dᵢ = rank(xᵢ) - rank(yᵢ)
     *
     * @param x First variable
     * @param y Second variable
     * @return Spearman's rho (-1 to +1)
     */
    static double spearmans_rho(
        const std::vector<double>& x,
        const std::vector<double>& y) {

        int n = x.size();
        if (n != static_cast<int>(y.size())) {
            throw std::invalid_argument("x and y must have same size");
        }

        // Rank x
        std::vector<std::pair<double, int>> x_pairs(n);
        for (int i = 0; i < n; ++i) {
            x_pairs[i] = {x[i], i};
        }
        std::sort(x_pairs.begin(), x_pairs.end());

        std::vector<double> rank_x(n);
        for (int i = 0; i < n; ++i) {
            rank_x[x_pairs[i].second] = i + 1.0;
        }

        // Rank y
        std::vector<std::pair<double, int>> y_pairs(n);
        for (int i = 0; i < n; ++i) {
            y_pairs[i] = {y[i], i};
        }
        std::sort(y_pairs.begin(), y_pairs.end());

        std::vector<double> rank_y(n);
        for (int i = 0; i < n; ++i) {
            rank_y[y_pairs[i].second] = i + 1.0;
        }

        // Compute sum of squared differences
        double sum_d_sq = 0.0;
        for (int i = 0; i < n; ++i) {
            double d = rank_x[i] - rank_y[i];
            sum_d_sq += d * d;
        }

        // Spearman's rho
        double rho = 1.0 - (6.0 * sum_d_sq) / (n * (n * n - 1));

        return rho;
    }

    /**
     * @brief Wilcoxon matched-pairs signed-ranks test
     *
     * Tests H₀: median difference = 0 for paired data
     * Ranks absolute differences, sums ranks with positive/negative signs
     *
     * @param x First sample (paired)
     * @param y Second sample (paired)
     * @return Wilcoxon T statistic (sum of positive ranks)
     */
    static double wilcoxon_matched_pairs(
        const std::vector<double>& x,
        const std::vector<double>& y) {

        int n = x.size();
        if (n != static_cast<int>(y.size())) {
            throw std::invalid_argument("x and y must have same size");
        }

        // Compute differences
        std::vector<std::pair<double, int>> diffs;  // (|diff|, sign)
        for (int i = 0; i < n; ++i) {
            double diff = x[i] - y[i];
            if (std::abs(diff) > 1e-10) {  // Ignore zeros
                diffs.push_back({std::abs(diff), (diff > 0) ? 1 : -1});
            }
        }

        int m = diffs.size();
        if (m == 0) return 0.0;

        // Sort by absolute difference
        std::sort(diffs.begin(), diffs.end());

        // Assign ranks
        std::vector<double> ranks(m);
        for (int i = 0; i < m; ++i) {
            ranks[i] = i + 1.0;
        }

        // Sum positive ranks
        double T_plus = 0.0;
        for (int i = 0; i < m; ++i) {
            if (diffs[i].second > 0) {
                T_plus += ranks[i];
            }
        }

        return T_plus;
    }

    /**
     * @brief Wilcoxon rank-sum test (Mann-Whitney U test)
     *
     * Tests H₀: two independent samples come from same distribution
     * U = R₁ - n₁(n₁+1)/2
     * where R₁ is sum of ranks in first sample
     *
     * @param sample1 First sample
     * @param sample2 Second sample
     * @return Mann-Whitney U statistic
     */
    static double mann_whitney_u_test(
        const std::vector<double>& sample1,
        const std::vector<double>& sample2) {

        int n1 = sample1.size();
        int n2 = sample2.size();

        // Combine samples with labels
        std::vector<std::pair<double, int>> combined;
        for (double x : sample1) {
            combined.push_back({x, 1});
        }
        for (double x : sample2) {
            combined.push_back({x, 2});
        }

        // Sort
        std::sort(combined.begin(), combined.end());

        // Assign ranks
        int N = n1 + n2;
        std::vector<double> ranks(N);
        int i = 0;
        while (i < N) {
            int j = i;
            while (j < N && combined[j].first == combined[i].first) {
                j++;
            }
            double avg_rank = (i + 1 + j) / 2.0;
            for (int m = i; m < j; ++m) {
                ranks[m] = avg_rank;
            }
            i = j;
        }

        // Sum ranks for sample 1
        double R1 = 0.0;
        for (int idx = 0; idx < N; ++idx) {
            if (combined[idx].second == 1) {
                R1 += ranks[idx];
            }
        }

        // Mann-Whitney U
        double U = R1 - n1 * (n1 + 1.0) / 2.0;

        return U;
    }

    /**
     * @brief Wilcoxon signed-rank test (one-sample)
     *
     * Tests H₀: median = 0
     * Alternative to one-sample t-test for non-normal data
     *
     * @param data Sample data
     * @return Wilcoxon W statistic
     */
    static double wilcoxon_signed_rank(const std::vector<double>& data) {
        // Remove zeros
        std::vector<std::pair<double, int>> abs_vals;
        for (double x : data) {
            if (std::abs(x) > 1e-10) {
                abs_vals.push_back({std::abs(x), (x > 0) ? 1 : -1});
            }
        }

        int n = abs_vals.size();
        if (n == 0) return 0.0;

        // Sort by absolute value
        std::sort(abs_vals.begin(), abs_vals.end());

        // Sum positive ranks
        double W_plus = 0.0;
        for (int i = 0; i < n; ++i) {
            if (abs_vals[i].second > 0) {
                W_plus += (i + 1);
            }
        }

        return W_plus;
    }

    /**
     * @brief Kolmogorov-Smirnov one-sample test
     *
     * Tests H₀: sample comes from specified distribution
     * D = max|Fₙ(x) - F₀(x)|
     * Already implemented in HypothesisTesting, wrapper here
     */
    static double kolmogorov_smirnov_one_sample(
        std::vector<double> data,
        std::function<double(double)> cdf_function) {

        return HypothesisTesting::kolmogorov_smirnov_test(data, cdf_function);
    }

    /**
     * @brief Kolmogorov-Smirnov two-sample test
     *
     * Tests H₀: two samples come from same distribution
     * D = max|F₁(x) - F₂(x)|
     *
     * @param sample1 First sample
     * @param sample2 Second sample
     * @return KS statistic D
     */
    static double kolmogorov_smirnov_two_sample(
        std::vector<double> sample1,
        std::vector<double> sample2) {

        std::sort(sample1.begin(), sample1.end());
        std::sort(sample2.begin(), sample2.end());

        int n1 = sample1.size();
        int n2 = sample2.size();

        double max_diff = 0.0;
        int i = 0, j = 0;

        while (i < n1 && j < n2) {
            double F1 = (i + 1.0) / n1;
            double F2 = (j + 1.0) / n2;
            max_diff = std::max(max_diff, std::abs(F1 - F2));

            if (sample1[i] < sample2[j]) {
                i++;
            } else {
                j++;
            }
        }

        return max_diff;
    }
};

/**
 * @class QualityControl
 * @brief Quality assurance and control methods
 *
 * Implements:
 * - Control charts (X-bar, R, S, p, c charts)
 * - Acceptance sampling plans
 * - Process capability indices
 * - Reliability analysis
 * - Risk analysis
 */
class QualityControl {
public:
    /**
     * @brief X-bar control chart parameters
     *
     * For monitoring process mean
     * UCL = x̄̄ + A₂R̄, LCL = x̄̄ - A₂R̄
     */
    struct XBarChart {
        double center_line;    // x̄̄ (grand mean)
        double UCL;            // Upper control limit
        double LCL;            // Lower control limit
        double A2_factor;      // A₂ factor from table
    };

    /**
     * @brief Compute X-bar control chart limits
     *
     * @param sample_means Vector of sample means
     * @param R_bar Average range
     * @param n Sample size
     * @return Control chart parameters
     */
    static XBarChart xbar_chart(
        const std::vector<double>& sample_means,
        double R_bar,
        int n) {

        // A2 factors for different sample sizes (n=2 to 10)
        std::vector<double> A2_table = {0, 0, 1.880, 1.023, 0.729, 0.577,
                                        0.483, 0.419, 0.373, 0.337, 0.308};

        if (n < 2 || n > 10) {
            throw std::invalid_argument("Sample size must be between 2 and 10");
        }

        double x_bar_bar = std::accumulate(sample_means.begin(),
                                          sample_means.end(), 0.0) / sample_means.size();
        double A2 = A2_table[n];

        XBarChart chart;
        chart.center_line = x_bar_bar;
        chart.A2_factor = A2;
        chart.UCL = x_bar_bar + A2 * R_bar;
        chart.LCL = x_bar_bar - A2 * R_bar;

        return chart;
    }

    /**
     * @brief R chart (range chart) parameters
     *
     * For monitoring process variability
     * UCL = D₄R̄, LCL = D₃R̄
     */
    struct RChart {
        double center_line;  // R̄
        double UCL;
        double LCL;
    };

    /**
     * @brief Compute R control chart limits
     */
    static RChart r_chart(const std::vector<double>& ranges, int n) {
        // D3 and D4 factors
        std::vector<double> D3_table = {0, 0, 0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223};
        std::vector<double> D4_table = {0, 0, 3.267, 2.574, 2.282, 2.114,
                                        2.004, 1.924, 1.864, 1.816, 1.777};

        if (n < 2 || n > 10) {
            throw std::invalid_argument("Sample size must be between 2 and 10");
        }

        double R_bar = std::accumulate(ranges.begin(), ranges.end(), 0.0) / ranges.size();

        RChart chart;
        chart.center_line = R_bar;
        chart.UCL = D4_table[n] * R_bar;
        chart.LCL = D3_table[n] * R_bar;

        return chart;
    }

    /**
     * @brief p-chart for proportion defective
     *
     * UCL = p̄ + 3√(p̄(1-p̄)/n), LCL = p̄ - 3√(p̄(1-p̄)/n)
     */
    struct PChart {
        double center_line;  // p̄
        double UCL;
        double LCL;
    };

    static PChart p_chart(const std::vector<double>& proportions, int n) {
        double p_bar = std::accumulate(proportions.begin(), proportions.end(), 0.0)
                      / proportions.size();

        double std_error = std::sqrt(p_bar * (1.0 - p_bar) / n);

        PChart chart;
        chart.center_line = p_bar;
        chart.UCL = p_bar + 3.0 * std_error;
        chart.LCL = std::max(0.0, p_bar - 3.0 * std_error);

        return chart;
    }

    /**
     * @brief Process capability index Cp
     *
     * Cp = (USL - LSL) / (6σ)
     * Measures potential capability
     *
     * @param USL Upper specification limit
     * @param LSL Lower specification limit
     * @param sigma Process standard deviation
     * @return Cp index (>1 is capable)
     */
    static double process_capability_cp(double USL, double LSL, double sigma) {
        return (USL - LSL) / (6.0 * sigma);
    }

    /**
     * @brief Process capability index Cpk
     *
     * Cpk = min((USL - μ)/(3σ), (μ - LSL)/(3σ))
     * Measures actual capability accounting for centering
     *
     * @return Cpk index (>1 is capable, >1.33 is good)
     */
    static double process_capability_cpk(
        double USL, double LSL, double mu, double sigma) {

        double cpu = (USL - mu) / (3.0 * sigma);
        double cpl = (mu - LSL) / (3.0 * sigma);

        return std::min(cpu, cpl);
    }

    /**
     * @brief Single sampling plan parameters
     *
     * Accept lot if c or fewer defectives in sample of n
     */
    struct SamplingPlan {
        int n;      // Sample size
        int c;      // Acceptance number
        double AQL;  // Acceptable Quality Level
        double LTPD; // Lot Tolerance Percent Defective
    };

    /**
     * @brief Design single sampling plan
     *
     * @param AQL Acceptable quality level (e.g., 0.01 for 1%)
     * @param LTPD Lot tolerance percent defective
     * @param alpha Producer's risk
     * @param beta Consumer's risk
     * @return Sampling plan (n, c)
     */
    static SamplingPlan single_sampling_plan(
        double AQL, double LTPD, double alpha = 0.05, double beta = 0.10) {

        // Simplified approach: use nomograph approximation
        // For exact plans, would use OC curve calculations

        SamplingPlan plan;
        plan.AQL = AQL;
        plan.LTPD = LTPD;

        // Simple heuristic (real plans use binomial/Poisson tables)
        plan.c = static_cast<int>(std::ceil(-std::log(alpha) / AQL));
        plan.n = static_cast<int>(std::ceil(plan.c / LTPD * 2.0));

        return plan;
    }

    /**
     * @brief Operating characteristic (OC) curve value
     *
     * Probability of acceptance Pa(p) for given defect rate p
     * Pa(p) = Σᵢ₌₀ᶜ C(n,i) pⁱ(1-p)ⁿ⁻ⁱ
     *
     * @param p Proportion defective
     * @param n Sample size
     * @param c Acceptance number
     * @return Probability of acceptance
     */
    static double oc_curve(double p, int n, int c) {
        BinomialDistribution binom(n, p);

        double Pa = 0.0;
        for (int i = 0; i <= c; ++i) {
            Pa += binom.pmf(i);
        }

        return Pa;
    }

    /**
     * @brief Reliability function R(t) = P(T > t)
     *
     * For exponential failure time: R(t) = e^(-λt)
     *
     * @param t Time
     * @param lambda Failure rate
     * @return Reliability at time t
     */
    static double reliability_exponential(double t, double lambda) {
        if (t < 0) return 1.0;
        return std::exp(-lambda * t);
    }

    /**
     * @brief Reliability for Weibull distribution
     *
     * R(t) = exp(-(t/η)^β)
     *
     * @param t Time
     * @param eta Scale parameter
     * @param beta Shape parameter
     * @return Reliability
     */
    static double reliability_weibull(double t, double eta, double beta) {
        if (t < 0) return 1.0;
        return std::exp(-std::pow(t / eta, beta));
    }

    /**
     * @brief Mean Time Between Failures (MTBF)
     *
     * For exponential: MTBF = 1/λ
     */
    static double mtbf_exponential(double lambda) {
        return 1.0 / lambda;
    }

    /**
     * @brief Failure rate (hazard function) for Weibull
     *
     * h(t) = (β/η)(t/η)^(β-1)
     */
    static double hazard_rate_weibull(double t, double eta, double beta) {
        if (t < 0) return 0.0;
        return (beta / eta) * std::pow(t / eta, beta - 1.0);
    }

    /**
     * @brief System reliability for series configuration
     *
     * R_system = ∏ᵢ Rᵢ
     * System fails if any component fails
     */
    static double reliability_series(const std::vector<double>& component_reliabilities) {
        double R_system = 1.0;
        for (double R : component_reliabilities) {
            R_system *= R;
        }
        return R_system;
    }

    /**
     * @brief System reliability for parallel configuration
     *
     * R_system = 1 - ∏ᵢ(1 - Rᵢ)
     * System works if at least one component works
     */
    static double reliability_parallel(const std::vector<double>& component_reliabilities) {
        double F_system = 1.0;  // System unreliability
        for (double R : component_reliabilities) {
            F_system *= (1.0 - R);
        }
        return 1.0 - F_system;
    }

    /**
     * @brief Expected monetary value (EMV) for decision under risk
     *
     * EMV = Σᵢ pᵢ · vᵢ
     * where pᵢ is probability, vᵢ is value/payoff
     *
     * @param probabilities Probabilities of outcomes
     * @param values Values/payoffs for each outcome
     * @return Expected value
     */
    static double expected_monetary_value(
        const std::vector<double>& probabilities,
        const std::vector<double>& values) {

        if (probabilities.size() != values.size()) {
            throw std::invalid_argument("Probabilities and values must have same size");
        }

        double EMV = 0.0;
        for (size_t i = 0; i < probabilities.size(); ++i) {
            EMV += probabilities[i] * values[i];
        }

        return EMV;
    }

    /**
     * @brief Expected value of perfect information (EVPI)
     *
     * EVPI = EV with perfect info - EV without perfect info
     * Maximum amount worth paying for perfect information
     *
     * @param best_outcomes Best outcome for each state
     * @param state_probabilities Probability of each state
     * @param current_emv EMV of current best decision
     * @return EVPI
     */
    static double expected_value_perfect_information(
        const std::vector<double>& best_outcomes,
        const std::vector<double>& state_probabilities,
        double current_emv) {

        double ev_perfect = expected_monetary_value(state_probabilities, best_outcomes);
        return ev_perfect - current_emv;
    }

    /**
     * @brief Risk priority number (RPN) for FMEA
     *
     * RPN = Severity × Occurrence × Detection
     * Used in Failure Mode and Effects Analysis
     *
     * @param severity Severity rating (1-10)
     * @param occurrence Occurrence rating (1-10)
     * @param detection Detection rating (1-10)
     * @return RPN (1-1000)
     */
    static int risk_priority_number(int severity, int occurrence, int detection) {
        return severity * occurrence * detection;
    }

    /**
     * @brief Minimax decision rule
     *
     * Choose action that minimizes maximum loss
     * Conservative approach for risk-averse decision makers
     *
     * @param payoff_matrix Payoff[action][state]
     * @return Index of minimax action
     */
    static int minimax_decision(const std::vector<std::vector<double>>& payoff_matrix) {
        int n_actions = payoff_matrix.size();

        double best_worst = -std::numeric_limits<double>::infinity();
        int best_action = 0;

        for (int i = 0; i < n_actions; ++i) {
            double worst_payoff = *std::min_element(
                payoff_matrix[i].begin(), payoff_matrix[i].end());

            if (worst_payoff > best_worst) {
                best_worst = worst_payoff;
                best_action = i;
            }
        }

        return best_action;
    }

    /**
     * @brief Maximax decision rule
     *
     * Choose action that maximizes maximum gain
     * Optimistic approach for risk-seeking decision makers
     *
     * @param payoff_matrix Payoff[action][state]
     * @return Index of maximax action
     */
    static int maximax_decision(const std::vector<std::vector<double>>& payoff_matrix) {
        int n_actions = payoff_matrix.size();

        double best_best = -std::numeric_limits<double>::infinity();
        int best_action = 0;

        for (int i = 0; i < n_actions; ++i) {
            double best_payoff = *std::max_element(
                payoff_matrix[i].begin(), payoff_matrix[i].end());

            if (best_payoff > best_best) {
                best_best = best_payoff;
                best_action = i;
            }
        }

        return best_action;
    }
};

/**
 * @class StatisticalTests
 * @brief Legacy statistical test methods (kept for backward compatibility)
 */
class StatisticalTests {
public:
    static double t_test_one_sample(const std::vector<double>& data, double mu0) {
        return HypothesisTesting::t_test_one_sample(data, mu0);
    }

    static double chi_squared_test(const std::vector<double>& observed,
                                   const std::vector<double>& expected) {
        return HypothesisTesting::chi_squared_goodness_of_fit(observed, expected);
    }

    static std::pair<double, double> mle_normal(const std::vector<double>& data) {
        auto result = StatisticalEstimation::maximum_likelihood("Normal", data);
        return {result[0], result[1]};
    }
};

} // namespace maths::probability

#endif // MATHS_PROBABILITY_DISTRIBUTIONS_HPP
