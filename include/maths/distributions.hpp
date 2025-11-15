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
 * @class StatisticalTests
 * @brief Hypothesis testing and statistical inference
 */
class StatisticalTests {
public:
    /**
     * @brief One-sample t-test
     *
     * Tests H₀: μ = μ₀ vs H₁: μ ≠ μ₀
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
     * @brief Chi-squared goodness-of-fit test
     *
     * @param observed Observed frequencies
     * @param expected Expected frequencies
     * @return Chi-squared statistic
     */
    static double chi_squared_test(const std::vector<double>& observed,
                                   const std::vector<double>& expected) {
        if (observed.size() != expected.size()) {
            throw std::invalid_argument("Observed and expected must have same size");
        }

        double chi_sq = 0.0;
        for (size_t i = 0; i < observed.size(); ++i) {
            double diff = observed[i] - expected[i];
            chi_sq += (diff * diff) / expected[i];
        }
        return chi_sq;
    }

    /**
     * @brief Maximum likelihood estimation for normal distribution
     *
     * @param data Sample data
     * @return Pair (mean, standard_deviation)
     */
    static std::pair<double, double> mle_normal(const std::vector<double>& data) {
        double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();

        double variance = 0.0;
        for (double x : data) {
            variance += (x - mean) * (x - mean);
        }
        variance /= data.size();

        return {mean, std::sqrt(variance)};
    }
};

} // namespace maths::probability

#endif // MATHS_PROBABILITY_DISTRIBUTIONS_HPP
