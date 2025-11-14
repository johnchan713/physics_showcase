#ifndef MATHS_PROBABILITY_DISTRIBUTIONS_HPP
#define MATHS_PROBABILITY_DISTRIBUTIONS_HPP

#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>
#include <random>

/**
 * @file distributions.hpp
 * @brief Probability distributions and statistical functions
 *
 * Implements:
 * - Discrete distributions (Bernoulli, Binomial, Poisson)
 * - Continuous distributions (Uniform, Normal, Exponential)
 * - Distribution properties (mean, variance, PDF, CDF)
 * - Statistical functions
 */

namespace maths::probability {

/**
 * @class BernoulliDistribution
 * @brief Bernoulli distribution (single trial, success/failure)
 *
 * P(X = 1) = p, P(X = 0) = 1 - p
 */
class BernoulliDistribution {
public:
    /**
     * @brief Probability mass function
     */
    static double pmf(int x, double p) {
        if (p < 0.0 || p > 1.0) {
            throw std::invalid_argument("Probability must be in [0, 1]");
        }
        if (x == 0) return 1.0 - p;
        if (x == 1) return p;
        return 0.0;
    }

    /**
     * @brief Mean (expected value)
     */
    static double mean(double p) {
        return p;
    }

    /**
     * @brief Variance
     */
    static double variance(double p) {
        return p * (1.0 - p);
    }

    /**
     * @brief Standard deviation
     */
    static double stddev(double p) {
        return std::sqrt(variance(p));
    }
};

/**
 * @class BinomialDistribution
 * @brief Binomial distribution (n independent Bernoulli trials)
 *
 * P(X = k) = C(n, k) * p^k * (1-p)^(n-k)
 */
class BinomialDistribution {
public:
    /**
     * @brief Binomial coefficient C(n, k) = n! / (k!(n-k)!)
     */
    static double binomialCoefficient(int n, int k) {
        if (k < 0 || k > n) return 0.0;
        if (k == 0 || k == n) return 1.0;

        // Use symmetry: C(n,k) = C(n,n-k)
        if (k > n - k) k = n - k;

        double result = 1.0;
        for (int i = 0; i < k; ++i) {
            result *= (n - i);
            result /= (i + 1);
        }
        return result;
    }

    /**
     * @brief Probability mass function
     */
    static double pmf(int k, int n, double p) {
        if (p < 0.0 || p > 1.0) {
            throw std::invalid_argument("Probability must be in [0, 1]");
        }
        if (k < 0 || k > n) return 0.0;

        double coef = binomialCoefficient(n, k);
        return coef * std::pow(p, k) * std::pow(1.0 - p, n - k);
    }

    /**
     * @brief Mean
     */
    static double mean(int n, double p) {
        return n * p;
    }

    /**
     * @brief Variance
     */
    static double variance(int n, double p) {
        return n * p * (1.0 - p);
    }

    /**
     * @brief Standard deviation
     */
    static double stddev(int n, double p) {
        return std::sqrt(variance(n, p));
    }
};

/**
 * @class PoissonDistribution
 * @brief Poisson distribution (events in fixed interval)
 *
 * P(X = k) = (λ^k * e^(-λ)) / k!
 */
class PoissonDistribution {
public:
    /**
     * @brief Factorial
     */
    static double factorial(int n) {
        if (n < 0) throw std::invalid_argument("Factorial undefined for negative numbers");
        if (n == 0 || n == 1) return 1.0;

        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }

    /**
     * @brief Probability mass function
     */
    static double pmf(int k, double lambda) {
        if (lambda <= 0.0) {
            throw std::invalid_argument("Lambda must be positive");
        }
        if (k < 0) return 0.0;

        return std::pow(lambda, k) * std::exp(-lambda) / factorial(k);
    }

    /**
     * @brief Mean
     */
    static double mean(double lambda) {
        return lambda;
    }

    /**
     * @brief Variance
     */
    static double variance(double lambda) {
        return lambda;
    }

    /**
     * @brief Standard deviation
     */
    static double stddev(double lambda) {
        return std::sqrt(lambda);
    }
};

/**
 * @class UniformDistribution
 * @brief Continuous uniform distribution on [a, b]
 */
class UniformDistribution {
public:
    /**
     * @brief Probability density function
     */
    static double pdf(double x, double a, double b) {
        if (a >= b) {
            throw std::invalid_argument("Must have a < b");
        }
        if (x < a || x > b) return 0.0;
        return 1.0 / (b - a);
    }

    /**
     * @brief Cumulative distribution function
     */
    static double cdf(double x, double a, double b) {
        if (a >= b) {
            throw std::invalid_argument("Must have a < b");
        }
        if (x < a) return 0.0;
        if (x > b) return 1.0;
        return (x - a) / (b - a);
    }

    /**
     * @brief Mean
     */
    static double mean(double a, double b) {
        return (a + b) / 2.0;
    }

    /**
     * @brief Variance
     */
    static double variance(double a, double b) {
        return (b - a) * (b - a) / 12.0;
    }

    /**
     * @brief Standard deviation
     */
    static double stddev(double a, double b) {
        return std::sqrt(variance(a, b));
    }
};

/**
 * @class NormalDistribution
 * @brief Normal (Gaussian) distribution
 *
 * PDF: f(x) = (1/√(2πσ²)) exp(-(x-μ)²/(2σ²))
 */
class NormalDistribution {
public:
    /**
     * @brief Probability density function
     */
    static double pdf(double x, double mu, double sigma) {
        if (sigma <= 0.0) {
            throw std::invalid_argument("Sigma must be positive");
        }

        double coefficient = 1.0 / (sigma * std::sqrt(2.0 * M_PI));
        double exponent = -0.5 * std::pow((x - mu) / sigma, 2.0);
        return coefficient * std::exp(exponent);
    }

    /**
     * @brief Cumulative distribution function (approximate)
     *
     * Uses error function approximation
     */
    static double cdf(double x, double mu, double sigma) {
        if (sigma <= 0.0) {
            throw std::invalid_argument("Sigma must be positive");
        }

        // Use std::erf for error function
        double z = (x - mu) / (sigma * std::sqrt(2.0));
        return 0.5 * (1.0 + std::erf(z));
    }

    /**
     * @brief Standard normal PDF (μ=0, σ=1)
     */
    static double standardPDF(double z) {
        return pdf(z, 0.0, 1.0);
    }

    /**
     * @brief Standard normal CDF
     */
    static double standardCDF(double z) {
        return cdf(z, 0.0, 1.0);
    }

    /**
     * @brief Mean
     */
    static double mean(double mu, double sigma) {
        return mu;
    }

    /**
     * @brief Variance
     */
    static double variance(double mu, double sigma) {
        return sigma * sigma;
    }

    /**
     * @brief Standard deviation
     */
    static double stddev(double mu, double sigma) {
        return sigma;
    }

    /**
     * @brief 68-95-99.7 rule (empirical rule)
     */
    static std::string empiricalRule() {
        return "Empirical Rule (68-95-99.7):\n"
               "\n"
               "For normal distribution:\n"
               "~68% of data within μ ± σ\n"
               "~95% of data within μ ± 2σ\n"
               "~99.7% of data within μ ± 3σ";
    }
};

/**
 * @class ExponentialDistribution
 * @brief Exponential distribution (time between events)
 *
 * PDF: f(x) = λ * e^(-λx) for x ≥ 0
 */
class ExponentialDistribution {
public:
    /**
     * @brief Probability density function
     */
    static double pdf(double x, double lambda) {
        if (lambda <= 0.0) {
            throw std::invalid_argument("Lambda must be positive");
        }
        if (x < 0.0) return 0.0;
        return lambda * std::exp(-lambda * x);
    }

    /**
     * @brief Cumulative distribution function
     */
    static double cdf(double x, double lambda) {
        if (lambda <= 0.0) {
            throw std::invalid_argument("Lambda must be positive");
        }
        if (x < 0.0) return 0.0;
        return 1.0 - std::exp(-lambda * x);
    }

    /**
     * @brief Mean
     */
    static double mean(double lambda) {
        return 1.0 / lambda;
    }

    /**
     * @brief Variance
     */
    static double variance(double lambda) {
        return 1.0 / (lambda * lambda);
    }

    /**
     * @brief Standard deviation
     */
    static double stddev(double lambda) {
        return 1.0 / lambda;
    }

    /**
     * @brief Memoryless property
     */
    static std::string memorylessProperty() {
        return "Memoryless Property:\n"
               "\n"
               "P(X > s + t | X > s) = P(X > t)\n"
               "\n"
               "\"The distribution doesn't remember the past\"\n"
               "\n"
               "Only continuous distribution with this property!";
    }
};

/**
 * @class StatisticalFunctions
 * @brief Common statistical functions
 */
class StatisticalFunctions {
public:
    /**
     * @brief Sample mean
     */
    static double mean(const std::vector<double>& data) {
        if (data.empty()) {
            throw std::invalid_argument("Data cannot be empty");
        }

        double sum = 0.0;
        for (double x : data) {
            sum += x;
        }
        return sum / data.size();
    }

    /**
     * @brief Sample variance (unbiased estimator)
     */
    static double variance(const std::vector<double>& data) {
        if (data.size() < 2) {
            throw std::invalid_argument("Need at least 2 data points");
        }

        double m = mean(data);
        double sum_sq = 0.0;
        for (double x : data) {
            double diff = x - m;
            sum_sq += diff * diff;
        }
        return sum_sq / (data.size() - 1);  // Bessel's correction
    }

    /**
     * @brief Sample standard deviation
     */
    static double stddev(const std::vector<double>& data) {
        return std::sqrt(variance(data));
    }

    /**
     * @brief Median
     */
    static double median(std::vector<double> data) {
        if (data.empty()) {
            throw std::invalid_argument("Data cannot be empty");
        }

        std::sort(data.begin(), data.end());
        size_t n = data.size();

        if (n % 2 == 0) {
            return (data[n/2 - 1] + data[n/2]) / 2.0;
        } else {
            return data[n/2];
        }
    }

    /**
     * @brief Covariance of two datasets
     */
    static double covariance(const std::vector<double>& x,
                            const std::vector<double>& y) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Datasets must have same size");
        }
        if (x.size() < 2) {
            throw std::invalid_argument("Need at least 2 data points");
        }

        double mean_x = mean(x);
        double mean_y = mean(y);

        double sum = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            sum += (x[i] - mean_x) * (y[i] - mean_y);
        }

        return sum / (x.size() - 1);
    }

    /**
     * @brief Correlation coefficient (Pearson's r)
     */
    static double correlation(const std::vector<double>& x,
                             const std::vector<double>& y) {
        double cov = covariance(x, y);
        double std_x = stddev(x);
        double std_y = stddev(y);

        if (std_x == 0.0 || std_y == 0.0) {
            throw std::runtime_error("Cannot compute correlation with zero variance");
        }

        return cov / (std_x * std_y);
    }

    /**
     * @brief Z-score (standardize value)
     */
    static double zScore(double x, double mu, double sigma) {
        if (sigma == 0.0) {
            throw std::invalid_argument("Sigma cannot be zero");
        }
        return (x - mu) / sigma;
    }

    /**
     * @brief Central Limit Theorem
     */
    static std::string centralLimitTheorem() {
        return "Central Limit Theorem:\n"
               "\n"
               "For large n, the sample mean X̄ of n iid random variables\n"
               "is approximately normally distributed:\n"
               "\n"
               "X̄ ~ N(μ, σ²/n)\n"
               "\n"
               "where μ = E[X], σ² = Var(X)\n"
               "\n"
               "Holds regardless of original distribution!\n"
               "(Usually n ≥ 30 is sufficient)";
    }

    /**
     * @brief Law of Large Numbers
     */
    static std::string lawOfLargeNumbers() {
        return "Law of Large Numbers:\n"
               "\n"
               "As n → ∞, sample mean X̄ₙ → μ (expected value)\n"
               "\n"
               "Weak LLN: Convergence in probability\n"
               "Strong LLN: Almost sure convergence\n"
               "\n"
               "Justifies using sample mean to estimate μ";
    }
};

/**
 * @class BayesianInference
 * @brief Bayesian probability and inference
 */
class BayesianInference {
public:
    /**
     * @brief Bayes' Theorem
     *
     * P(A|B) = P(B|A) * P(A) / P(B)
     */
    static double bayesTheorem(double p_b_given_a, double p_a, double p_b) {
        if (p_b == 0.0) {
            throw std::invalid_argument("P(B) cannot be zero");
        }
        return (p_b_given_a * p_a) / p_b;
    }

    /**
     * @brief Statement of Bayes' Theorem
     */
    static std::string statement() {
        return "Bayes' Theorem:\n"
               "\n"
               "P(A|B) = P(B|A) * P(A) / P(B)\n"
               "\n"
               "where:\n"
               "P(A|B) = posterior probability\n"
               "P(B|A) = likelihood\n"
               "P(A) = prior probability\n"
               "P(B) = marginal likelihood (evidence)\n"
               "\n"
               "In words: \"Update beliefs based on new evidence\"";
    }

    /**
     * @brief Law of Total Probability
     */
    static std::string totalProbability() {
        return "Law of Total Probability:\n"
               "\n"
               "If {A₁, A₂, ..., Aₙ} partition the sample space:\n"
               "\n"
               "P(B) = Σᵢ P(B|Aᵢ) * P(Aᵢ)\n"
               "\n"
               "Used to compute P(B) in Bayes' Theorem";
    }
};

} // namespace maths::probability

#endif // MATHS_PROBABILITY_DISTRIBUTIONS_HPP
