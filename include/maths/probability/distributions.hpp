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
