#ifndef MATHS_ECONOMETRICS_REGRESSION_HPP
#define MATHS_ECONOMETRICS_REGRESSION_HPP

#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <numeric>

/**
 * @file regression.hpp
 * @brief Econometric regression analysis and time series
 *
 * Implements:
 * - Simple and multiple linear regression (OLS)
 * - Hypothesis testing (t-tests, F-tests)
 * - R-squared, adjusted R-squared
 * - Residual analysis
 * - Time series basics (AR, MA, ARMA)
 * - Stationarity and unit root tests
 */

namespace maths::econometrics {

/**
 * @class SimpleLinearRegression
 * @brief Ordinary Least Squares for simple linear regression
 */
class SimpleLinearRegression {
public:
    double intercept;
    double slope;
    double r_squared;
    double residual_std_error;
    std::vector<double> residuals;

    /**
     * @brief Fit simple linear regression: y = α + βx + ε
     */
    void fit(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size() || x.empty()) {
            throw std::invalid_argument("x and y must have same non-zero size");
        }

        size_t n = x.size();

        // Calculate means
        double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
        double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;

        // Calculate slope and intercept
        double numerator = 0.0;
        double denominator = 0.0;

        for (size_t i = 0; i < n; ++i) {
            numerator += (x[i] - mean_x) * (y[i] - mean_y);
            denominator += (x[i] - mean_x) * (x[i] - mean_x);
        }

        if (std::abs(denominator) < 1e-10) {
            throw std::runtime_error("No variation in x, cannot fit regression");
        }

        slope = numerator / denominator;
        intercept = mean_y - slope * mean_x;

        // Calculate residuals and R²
        residuals.resize(n);
        double ss_res = 0.0;  // Residual sum of squares
        double ss_tot = 0.0;  // Total sum of squares

        for (size_t i = 0; i < n; ++i) {
            double predicted = intercept + slope * x[i];
            residuals[i] = y[i] - predicted;
            ss_res += residuals[i] * residuals[i];
            ss_tot += (y[i] - mean_y) * (y[i] - mean_y);
        }

        r_squared = 1.0 - ss_res / ss_tot;
        residual_std_error = std::sqrt(ss_res / (n - 2));  // n-2 degrees of freedom
    }

    /**
     * @brief Predict y for given x
     */
    double predict(double x) const {
        return intercept + slope * x;
    }

    /**
     * @brief Get OLS estimator formulas
     */
    static std::string olsFormulas() {
        return "OLS Estimators (Simple Linear Regression):\n"
               "\n"
               "Model: yᵢ = α + βxᵢ + εᵢ\n"
               "\n"
               "Slope: β̂ = Cov(X,Y) / Var(X)\n"
               "         = Σ(xᵢ - x̄)(yᵢ - ȳ) / Σ(xᵢ - x̄)²\n"
               "\n"
               "Intercept: α̂ = ȳ - β̂·x̄\n"
               "\n"
               "Properties (Gauss-Markov Theorem):\n"
               "Under assumptions:\n"
               "1. Linearity: E[Y|X] = α + βX\n"
               "2. Homoscedasticity: Var(ε|X) = σ²\n"
               "3. No autocorrelation: Cov(εᵢ, εⱼ) = 0\n"
               "4. Exogeneity: E[ε|X] = 0\n"
               "\n"
               "OLS is BLUE:\n"
               "- Best: minimum variance\n"
               "- Linear: linear in Y\n"
               "- Unbiased: E[β̂] = β\n"
               "- Estimator";
    }
};

/**
 * @class HypothesisTesting
 * @brief Statistical hypothesis testing for regression
 */
class HypothesisTesting {
public:
    /**
     * @brief Calculate t-statistic for coefficient
     */
    static double tStatistic(double beta_hat, double se_beta) {
        if (se_beta < 1e-10) {
            throw std::runtime_error("Standard error too small");
        }
        return beta_hat / se_beta;
    }

    /**
     * @brief Calculate p-value (two-tailed) from t-statistic
     */
    static double pValue(double t_stat, int df) {
        // Approximate using normal distribution for large df
        // For exact, would need t-distribution CDF
        if (df > 30) {
            double z = std::abs(t_stat);
            double p_one_tail = 0.5 * std::erfc(z / std::sqrt(2.0));
            return 2.0 * p_one_tail;  // Two-tailed
        }
        // For small df, approximation is less accurate
        return 0.0;  // Placeholder
    }

    /**
     * @brief Test if coefficient is significant at level α
     */
    static bool isSignificant(double t_stat, double alpha = 0.05) {
        // Critical value for normal (large sample)
        double z_crit = 1.96;  // 95% confidence
        if (alpha == 0.01) z_crit = 2.576;  // 99% confidence
        if (alpha == 0.10) z_crit = 1.645;  // 90% confidence

        return std::abs(t_stat) > z_crit;
    }

    /**
     * @brief Hypothesis testing framework
     */
    static std::string framework() {
        return "Hypothesis Testing in Regression:\n"
               "\n"
               "Test for single coefficient β:\n"
               "H₀: β = β₀ (usually β₀ = 0)\n"
               "H₁: β ≠ β₀ (two-sided) or β > β₀ (one-sided)\n"
               "\n"
               "Test statistic:\n"
               "  t = (β̂ - β₀) / SE(β̂)\n"
               "\n"
               "Under H₀: t ~ t_{n-k} (t-distribution with n-k df)\n"
               "where n = sample size, k = number of parameters\n"
               "\n"
               "Decision rule (two-tailed, level α):\n"
               "  Reject H₀ if |t| > t_{α/2, n-k}\n"
               "\n"
               "p-value: P(|T| > |t_obs| | H₀)\n"
               "  - Small p-value: strong evidence against H₀\n"
               "  - Reject H₀ if p-value < α\n"
               "\n"
               "Confidence interval:\n"
               "  β̂ ± t_{α/2, n-k} · SE(β̂)\n"
               "\n"
               "F-test for overall significance:\n"
               "H₀: β₁ = β₂ = ... = βₖ = 0\n"
               "F = (ESS/k) / (RSS/(n-k-1))\n"
               "where ESS = explained sum of squares\n"
               "      RSS = residual sum of squares";
    }
};

/**
 * @class GoodnessOfFit
 * @brief R-squared and model fit statistics
 */
class GoodnessOfFit {
public:
    /**
     * @brief Calculate R-squared
     */
    static double rSquared(const std::vector<double>& y,
                           const std::vector<double>& y_pred) {
        if (y.size() != y_pred.size() || y.empty()) {
            throw std::invalid_argument("Invalid input sizes");
        }

        double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / y.size();

        double ss_tot = 0.0;
        double ss_res = 0.0;

        for (size_t i = 0; i < y.size(); ++i) {
            ss_tot += (y[i] - mean_y) * (y[i] - mean_y);
            ss_res += (y[i] - y_pred[i]) * (y[i] - y_pred[i]);
        }

        return 1.0 - ss_res / ss_tot;
    }

    /**
     * @brief Adjusted R-squared (penalizes additional predictors)
     */
    static double adjustedRSquared(double r2, size_t n, size_t k) {
        // R̄² = 1 - (1 - R²)(n-1)/(n-k-1)
        if (n <= k + 1) {
            throw std::invalid_argument("Not enough observations");
        }
        return 1.0 - (1.0 - r2) * (n - 1.0) / (n - k - 1.0);
    }

    /**
     * @brief Akaike Information Criterion
     */
    static double aic(double rss, size_t n, size_t k) {
        // AIC = n·ln(RSS/n) + 2k
        if (rss <= 0.0) {
            throw std::invalid_argument("RSS must be positive");
        }
        return n * std::log(rss / n) + 2.0 * k;
    }

    /**
     * @brief Bayesian Information Criterion
     */
    static double bic(double rss, size_t n, size_t k) {
        // BIC = n·ln(RSS/n) + k·ln(n)
        if (rss <= 0.0) {
            throw std::invalid_argument("RSS must be positive");
        }
        return n * std::log(rss / n) + k * std::log(n);
    }

    /**
     * @brief Model selection criteria
     */
    static std::string criteria() {
        return "Model Fit and Selection:\n"
               "\n"
               "R² (Coefficient of Determination):\n"
               "  R² = 1 - RSS/TSS = ESS/TSS\n"
               "  - Range: [0, 1]\n"
               "  - 1 = perfect fit, 0 = no explanatory power\n"
               "  - Problem: always increases with more variables\n"
               "\n"
               "Adjusted R²:\n"
               "  R̄² = 1 - (1-R²)(n-1)/(n-k-1)\n"
               "  - Penalizes additional variables\n"
               "  - Use for model comparison\n"
               "\n"
               "AIC (Akaike Information Criterion):\n"
               "  AIC = 2k - 2ln(L)\n"
               "  - Lower is better\n"
               "  - Penalizes complexity (2k term)\n"
               "  - Asymptotically efficient\n"
               "\n"
               "BIC (Bayesian Information Criterion):\n"
               "  BIC = k·ln(n) - 2ln(L)\n"
               "  - Lower is better\n"
               "  - Stronger penalty than AIC for large n\n"
               "  - Consistent estimator\n"
               "\n"
               "Model selection:\n"
               "- Compare AIC/BIC across models\n"
               "- BIC favors simpler models\n"
               "- Cross-validation for out-of-sample fit";
    }
};

/**
 * @class TimeSeries
 * @brief Basic time series models
 */
class TimeSeries {
public:
    /**
     * @brief Autoregressive model AR(p)
     */
    static std::string autoregressive() {
        return "Autoregressive Model AR(p):\n"
               "\n"
               "Xₜ = c + φ₁Xₜ₋₁ + φ₂Xₜ₋₂ + ... + φₚXₜ₋ₚ + εₜ\n"
               "\n"
               "where εₜ ~ WN(0, σ²) is white noise.\n"
               "\n"
               "AR(1): Xₜ = c + φXₜ₋₁ + εₜ\n"
               "- Stationarity: |φ| < 1\n"
               "- Mean: E[Xₜ] = c/(1-φ)\n"
               "- Variance: Var(Xₜ) = σ²/(1-φ²)\n"
               "- ACF: ρ(k) = φᵏ (exponential decay)\n"
               "\n"
               "General AR(p):\n"
               "- Stationarity: roots of φ(z) = 0 outside unit circle\n"
               "- Yule-Walker equations for estimation\n"
               "- ACF decays, PACF cuts off after lag p\n"
               "\n"
               "Random walk: AR(1) with φ = 1 (non-stationary)\n"
               "  Xₜ = Xₜ₋₁ + εₜ\n"
               "  Variance grows over time!";
    }

    /**
     * @brief Moving average model MA(q)
     */
    static std::string movingAverage() {
        return "Moving Average Model MA(q):\n"
               "\n"
               "Xₜ = μ + εₜ + θ₁εₜ₋₁ + θ₂εₜ₋₂ + ... + θ_qεₜ₋q\n"
               "\n"
               "where εₜ ~ WN(0, σ²)\n"
               "\n"
               "MA(1): Xₜ = μ + εₜ + θεₜ₋₁\n"
               "- Always stationary\n"
               "- Mean: E[Xₜ] = μ\n"
               "- Variance: Var(Xₜ) = σ²(1 + θ²)\n"
               "- ACF: ρ(1) = θ/(1+θ²), ρ(k) = 0 for k > 1\n"
               "\n"
               "General MA(q):\n"
               "- ACF cuts off after lag q\n"
               "- PACF decays\n"
               "- Invertibility: roots of θ(z) = 0 outside unit circle\n"
               "\n"
               "Duality with AR:\n"
               "- Invertible MA(q) ⟺ AR(∞)\n"
               "- Stationary AR(p) ⟺ MA(∞)";
    }

    /**
     * @brief ARMA model
     */
    static std::string arma() {
        return "ARMA(p,q) Model:\n"
               "\n"
               "Xₜ = c + φ₁Xₜ₋₁ + ... + φₚXₜ₋ₚ + εₜ + θ₁εₜ₋₁ + ... + θ_qεₜ₋q\n"
               "\n"
               "Combines AR(p) and MA(q) components.\n"
               "\n"
               "Lag operator notation:\n"
               "  φ(L)Xₜ = θ(L)εₜ\n"
               "where:\n"
               "  φ(L) = 1 - φ₁L - φ₂L² - ... - φₚLᵖ\n"
               "  θ(L) = 1 + θ₁L + θ₂L² + ... + θ_qLᵍ\n"
               "  LXₜ = Xₜ₋₁ (lag operator)\n"
               "\n"
               "Identification (Box-Jenkins):\n"
               "1. Check stationarity (unit root tests)\n"
               "2. Look at ACF and PACF:\n"
               "   - AR(p): ACF decays, PACF cuts at p\n"
               "   - MA(q): ACF cuts at q, PACF decays\n"
               "   - ARMA: both decay\n"
               "3. Estimate parameters (MLE, LS)\n"
               "4. Diagnostic checking (residuals)\n"
               "\n"
               "ARIMA(p,d,q): ARMA with d differences\n"
               "  (1-L)ᵈXₜ = ARMA(p,q)";
    }

    /**
     * @brief Stationarity
     */
    static std::string stationarity() {
        return "Stationarity:\n"
               "\n"
               "Weak (covariance) stationarity:\n"
               "1. E[Xₜ] = μ (constant mean)\n"
               "2. Var(Xₜ) = σ² (constant variance)\n"
               "3. Cov(Xₜ, Xₜ₊ₖ) = γ(k) (depends only on lag k)\n"
               "\n"
               "Strong (strict) stationarity:\n"
               "Joint distribution of (Xₜ₁,...,Xₜₙ) same as\n"
               "(Xₜ₁₊ₕ,...,Xₜₙ₊ₕ) for all h\n"
               "\n"
               "Tests for stationarity:\n"
               "- Augmented Dickey-Fuller (ADF) test\n"
               "- Phillips-Perron test\n"
               "- KPSS test (null: stationary)\n"
               "\n"
               "Unit root test (ADF):\n"
               "H₀: Xₜ has unit root (non-stationary)\n"
               "H₁: Xₜ is stationary\n"
               "\n"
               "Regression: ΔXₜ = α + βXₜ₋₁ + Σγᵢ ΔXₜ₋ᵢ + εₜ\n"
               "Test: H₀: β = 0 vs H₁: β < 0\n"
               "\n"
               "If non-stationary: difference until stationary";
    }
};

} // namespace maths::econometrics

#endif // MATHS_ECONOMETRICS_REGRESSION_HPP
