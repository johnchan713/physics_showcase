#ifndef MATHS_FINANCE_BLACK_SCHOLES_HPP
#define MATHS_FINANCE_BLACK_SCHOLES_HPP

#include <cmath>
#include <functional>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

/**
 * @file black_scholes.hpp
 * @brief Black-Scholes model and derivatives pricing
 *
 * Implements:
 * - Black-Scholes formula for European options
 * - The Greeks (delta, gamma, vega, theta, rho)
 * - Put-Call parity
 * - Implied volatility
 * - Risk-neutral valuation
 * - Binomial tree models
 */

namespace maths::finance {

/**
 * @class BlackScholes
 * @brief Black-Scholes option pricing model
 */
class BlackScholes {
public:
    /**
     * @brief Standard normal cumulative distribution function
     */
    static double normalCDF(double x) {
        return 0.5 * std::erfc(-x * M_SQRT1_2);
    }

    /**
     * @brief Standard normal probability density function
     */
    static double normalPDF(double x) {
        return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI);
    }

    /**
     * @brief Calculate d1 parameter in Black-Scholes
     *
     * d1 = [ln(S/K) + (r + σ²/2)T] / (σ√T)
     */
    static double d1(double S, double K, double r, double sigma, double T) {
        if (T <= 0.0 || sigma <= 0.0) {
            throw std::invalid_argument("T and sigma must be positive");
        }
        return (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    }

    /**
     * @brief Calculate d2 parameter in Black-Scholes
     *
     * d2 = d1 - σ√T
     */
    static double d2(double S, double K, double r, double sigma, double T) {
        return d1(S, K, r, sigma, T) - sigma * std::sqrt(T);
    }

    /**
     * @brief European call option price
     *
     * C = S·N(d1) - K·e^(-rT)·N(d2)
     *
     * @param S Current stock price
     * @param K Strike price
     * @param r Risk-free rate (continuous compounding)
     * @param sigma Volatility (annualized)
     * @param T Time to maturity (years)
     */
    static double callPrice(double S, double K, double r, double sigma, double T) {
        if (S <= 0.0 || K <= 0.0) {
            throw std::invalid_argument("S and K must be positive");
        }
        if (T <= 0.0) return std::max(S - K, 0.0);  // Intrinsic value at expiry

        double d1_val = d1(S, K, r, sigma, T);
        double d2_val = d2(S, K, r, sigma, T);

        return S * normalCDF(d1_val) - K * std::exp(-r * T) * normalCDF(d2_val);
    }

    /**
     * @brief European put option price
     *
     * P = K·e^(-rT)·N(-d2) - S·N(-d1)
     */
    static double putPrice(double S, double K, double r, double sigma, double T) {
        if (S <= 0.0 || K <= 0.0) {
            throw std::invalid_argument("S and K must be positive");
        }
        if (T <= 0.0) return std::max(K - S, 0.0);  // Intrinsic value at expiry

        double d1_val = d1(S, K, r, sigma, T);
        double d2_val = d2(S, K, r, sigma, T);

        return K * std::exp(-r * T) * normalCDF(-d2_val) - S * normalCDF(-d1_val);
    }

    /**
     * @brief Put-Call parity relationship
     *
     * C - P = S - K·e^(-rT)
     */
    static std::string putCallParity() {
        return "Put-Call Parity:\n"
               "\n"
               "C - P = S - K·e^(-rT)\n"
               "\n"
               "where:\n"
               "C = European call price\n"
               "P = European put price\n"
               "S = Current stock price\n"
               "K = Strike price\n"
               "r = Risk-free rate\n"
               "T = Time to maturity\n"
               "\n"
               "This relationship holds by no-arbitrage:\n"
               "Portfolio 1: Long call + cash K·e^(-rT)\n"
               "Portfolio 2: Long put + long stock\n"
               "Both worth max(S_T, K) at maturity";
    }

    /**
     * @brief Black-Scholes assumptions
     */
    static std::string assumptions() {
        return "Black-Scholes Model Assumptions:\n"
               "\n"
               "1. Stock price follows geometric Brownian motion:\n"
               "   dS = μS dt + σS dW\n"
               "\n"
               "2. No dividends paid during option's life\n"
               "\n"
               "3. Markets are frictionless:\n"
               "   - No transaction costs\n"
               "   - No taxes\n"
               "   - Assets infinitely divisible\n"
               "\n"
               "4. Risk-free rate r is constant\n"
               "\n"
               "5. Volatility σ is constant\n"
               "\n"
               "6. No arbitrage opportunities\n"
               "\n"
               "7. Trading is continuous\n"
               "\n"
               "8. Short selling allowed with full use of proceeds\n"
               "\n"
               "In practice, these assumptions are violated, leading to:\n"
               "- Volatility smile/skew\n"
               "- Jump diffusion models\n"
               "- Stochastic volatility models (Heston, SABR)";
    }
};

/**
 * @class Greeks
 * @brief Option sensitivities (The Greeks)
 */
class Greeks {
public:
    /**
     * @brief Delta: ∂V/∂S (sensitivity to stock price)
     *
     * Call: Δ_c = N(d1)
     * Put:  Δ_p = N(d1) - 1 = -N(-d1)
     */
    static double callDelta(double S, double K, double r, double sigma, double T) {
        double d1_val = BlackScholes::d1(S, K, r, sigma, T);
        return BlackScholes::normalCDF(d1_val);
    }

    static double putDelta(double S, double K, double r, double sigma, double T) {
        double d1_val = BlackScholes::d1(S, K, r, sigma, T);
        return BlackScholes::normalCDF(d1_val) - 1.0;
    }

    /**
     * @brief Gamma: ∂²V/∂S² (rate of change of delta)
     *
     * Γ = N'(d1) / (S·σ·√T)
     *
     * Same for calls and puts!
     */
    static double gamma(double S, double K, double r, double sigma, double T) {
        if (T <= 0.0) return 0.0;

        double d1_val = BlackScholes::d1(S, K, r, sigma, T);
        return BlackScholes::normalPDF(d1_val) / (S * sigma * std::sqrt(T));
    }

    /**
     * @brief Vega: ∂V/∂σ (sensitivity to volatility)
     *
     * ν = S·N'(d1)·√T
     *
     * Same for calls and puts!
     * Note: Vega is not a Greek letter (should be kappa)
     */
    static double vega(double S, double K, double r, double sigma, double T) {
        if (T <= 0.0) return 0.0;

        double d1_val = BlackScholes::d1(S, K, r, sigma, T);
        return S * BlackScholes::normalPDF(d1_val) * std::sqrt(T);
    }

    /**
     * @brief Theta: ∂V/∂T (time decay)
     *
     * Call: Θ_c = -[S·N'(d1)·σ/(2√T)] - r·K·e^(-rT)·N(d2)
     * Put:  Θ_p = -[S·N'(d1)·σ/(2√T)] + r·K·e^(-rT)·N(-d2)
     *
     * Usually negative (options lose value as time passes)
     */
    static double callTheta(double S, double K, double r, double sigma, double T) {
        if (T <= 0.0) return 0.0;

        double d1_val = BlackScholes::d1(S, K, r, sigma, T);
        double d2_val = BlackScholes::d2(S, K, r, sigma, T);

        double term1 = -S * BlackScholes::normalPDF(d1_val) * sigma / (2.0 * std::sqrt(T));
        double term2 = -r * K * std::exp(-r * T) * BlackScholes::normalCDF(d2_val);

        return term1 + term2;
    }

    static double putTheta(double S, double K, double r, double sigma, double T) {
        if (T <= 0.0) return 0.0;

        double d1_val = BlackScholes::d1(S, K, r, sigma, T);
        double d2_val = BlackScholes::d2(S, K, r, sigma, T);

        double term1 = -S * BlackScholes::normalPDF(d1_val) * sigma / (2.0 * std::sqrt(T));
        double term2 = r * K * std::exp(-r * T) * BlackScholes::normalCDF(-d2_val);

        return term1 + term2;
    }

    /**
     * @brief Rho: ∂V/∂r (sensitivity to interest rate)
     *
     * Call: ρ_c = K·T·e^(-rT)·N(d2)
     * Put:  ρ_p = -K·T·e^(-rT)·N(-d2)
     */
    static double callRho(double S, double K, double r, double sigma, double T) {
        if (T <= 0.0) return 0.0;

        double d2_val = BlackScholes::d2(S, K, r, sigma, T);
        return K * T * std::exp(-r * T) * BlackScholes::normalCDF(d2_val);
    }

    static double putRho(double S, double K, double r, double sigma, double T) {
        if (T <= 0.0) return 0.0;

        double d2_val = BlackScholes::d2(S, K, r, sigma, T);
        return -K * T * std::exp(-r * T) * BlackScholes::normalCDF(-d2_val);
    }

    /**
     * @brief Description of Greeks
     */
    static std::string description() {
        return "The Greeks - Option Sensitivities:\n"
               "\n"
               "Delta (Δ): ∂V/∂S\n"
               "- Hedge ratio: shares per option for delta-neutral portfolio\n"
               "- Call: 0 < Δ < 1 (positive exposure to stock)\n"
               "- Put: -1 < Δ < 0 (negative exposure to stock)\n"
               "- ATM options: Δ ≈ ±0.5\n"
               "\n"
               "Gamma (Γ): ∂²V/∂S² = ∂Δ/∂S\n"
               "- Curvature of option value\n"
               "- Measures delta hedging error\n"
               "- Highest for ATM options\n"
               "- Always positive for long options\n"
               "\n"
               "Vega (ν): ∂V/∂σ\n"
               "- Sensitivity to implied volatility\n"
               "- Always positive for long options\n"
               "- Highest for ATM options\n"
               "- Volatility trading: buy low vol, sell high vol\n"
               "\n"
               "Theta (Θ): ∂V/∂t\n"
               "- Time decay (typically negative)\n"
               "- Options lose value as expiration approaches\n"
               "- Highest (in magnitude) for ATM options\n"
               "\n"
               "Rho (ρ): ∂V/∂r\n"
               "- Sensitivity to interest rates\n"
               "- Usually less important than other Greeks\n"
               "- Call: positive, Put: negative";
    }
};

/**
 * @class ImpliedVolatility
 * @brief Calculate implied volatility from option prices
 */
class ImpliedVolatility {
public:
    /**
     * @brief Calculate implied volatility using Newton-Raphson
     *
     * Solves: marketPrice = BlackScholes(σ) for σ
     */
    static double calculate(double marketPrice, double S, double K, double r, double T,
                           bool isCall = true, double tolerance = 1e-6, int maxIter = 100) {
        if (marketPrice <= 0.0) {
            throw std::invalid_argument("Market price must be positive");
        }

        // Initial guess: use Brenner-Subrahmanyam approximation
        double sigma = std::sqrt(2.0 * M_PI / T) * marketPrice / S;
        sigma = std::max(0.01, std::min(sigma, 5.0));  // Bound initial guess

        for (int i = 0; i < maxIter; ++i) {
            double price = isCall ?
                BlackScholes::callPrice(S, K, r, sigma, T) :
                BlackScholes::putPrice(S, K, r, sigma, T);

            double diff = price - marketPrice;

            if (std::abs(diff) < tolerance) {
                return sigma;
            }

            // Vega for Newton-Raphson
            double vega_val = Greeks::vega(S, K, r, sigma, T);

            if (vega_val < 1e-10) {
                throw std::runtime_error("Vega too small, cannot compute IV");
            }

            // Newton-Raphson update
            sigma -= diff / vega_val;

            // Keep sigma in reasonable range
            sigma = std::max(0.001, std::min(sigma, 5.0));
        }

        throw std::runtime_error("Implied volatility did not converge");
    }

    /**
     * @brief Volatility smile/skew
     */
    static std::string volatilitySmile() {
        return "Implied Volatility Smile/Skew:\n"
               "\n"
               "In practice, implied volatility varies with strike:\n"
               "- Black-Scholes assumes constant σ\n"
               "- Market prices reveal σ_impl(K) depends on K\n"
               "\n"
               "Equity markets: Volatility skew\n"
               "- OTM puts have higher IV than ATM\n"
               "- Reflects crash risk (left tail is fatter)\n"
               "- Post-1987 crash phenomenon\n"
               "\n"
               "FX markets: Volatility smile\n"
               "- Both OTM puts and calls have higher IV\n"
               "- Symmetric around ATM\n"
               "- Reflects fat tails in both directions\n"
               "\n"
               "Causes:\n"
               "- Jump risk (non-continuous price movements)\n"
               "- Stochastic volatility\n"
               "- Supply/demand for protection\n"
               "\n"
               "Models addressing smile:\n"
               "- Local volatility (Dupire)\n"
               "- Stochastic volatility (Heston, SABR)\n"
               "- Jump diffusion (Merton)";
    }
};

/**
 * @class BinomialTree
 * @brief Binomial tree option pricing (Cox-Ross-Rubinstein)
 */
class BinomialTree {
public:
    /**
     * @brief Price American or European option using binomial tree
     *
     * @param S Current stock price
     * @param K Strike price
     * @param r Risk-free rate
     * @param sigma Volatility
     * @param T Time to maturity
     * @param steps Number of time steps
     * @param isCall true for call, false for put
     * @param isAmerican true for American, false for European
     */
    static double price(double S, double K, double r, double sigma, double T,
                       int steps, bool isCall, bool isAmerican) {
        if (steps <= 0) {
            throw std::invalid_argument("Number of steps must be positive");
        }

        double dt = T / steps;
        double u = std::exp(sigma * std::sqrt(dt));  // Up factor
        double d = 1.0 / u;                           // Down factor
        double p = (std::exp(r * dt) - d) / (u - d); // Risk-neutral probability

        // Stock prices at maturity
        std::vector<double> stockPrices(steps + 1);
        for (int i = 0; i <= steps; ++i) {
            stockPrices[i] = S * std::pow(u, steps - i) * std::pow(d, i);
        }

        // Option values at maturity
        std::vector<double> optionValues(steps + 1);
        for (int i = 0; i <= steps; ++i) {
            if (isCall) {
                optionValues[i] = std::max(stockPrices[i] - K, 0.0);
            } else {
                optionValues[i] = std::max(K - stockPrices[i], 0.0);
            }
        }

        // Backward induction
        for (int step = steps - 1; step >= 0; --step) {
            for (int i = 0; i <= step; ++i) {
                // Risk-neutral valuation
                double holdValue = std::exp(-r * dt) * (p * optionValues[i] + (1 - p) * optionValues[i + 1]);

                if (isAmerican) {
                    // Early exercise value
                    double S_current = S * std::pow(u, step - i) * std::pow(d, i);
                    double exerciseValue = isCall ?
                        std::max(S_current - K, 0.0) :
                        std::max(K - S_current, 0.0);

                    optionValues[i] = std::max(holdValue, exerciseValue);
                } else {
                    optionValues[i] = holdValue;
                }
            }
        }

        return optionValues[0];
    }

    /**
     * @brief Description of binomial tree method
     */
    static std::string description() {
        return "Binomial Tree Model (Cox-Ross-Rubinstein):\n"
               "\n"
               "Discrete-time approximation to continuous process:\n"
               "- Divide time [0, T] into N steps of length Δt = T/N\n"
               "- At each step, stock moves up or down:\n"
               "  u = e^(σ√Δt) (up factor)\n"
               "  d = 1/u = e^(-σ√Δt) (down factor)\n"
               "\n"
               "Risk-neutral probability:\n"
               "  p = (e^(rΔt) - d) / (u - d)\n"
               "\n"
               "Backward induction:\n"
               "1. Calculate terminal payoffs\n"
               "2. Work backwards: V_t = e^(-rΔt)[pV_u + (1-p)V_d]\n"
               "3. For American: V_t = max(hold, exercise)\n"
               "\n"
               "Advantages:\n"
               "- Can price American options\n"
               "- Handles dividends easily\n"
               "- Intuitive and flexible\n"
               "\n"
               "Convergence: As N → ∞, converges to Black-Scholes\n"
               "for European options";
    }
};

/**
 * @class RiskNeutralValuation
 * @brief Risk-neutral pricing framework
 */
class RiskNeutralValuation {
public:
    /**
     * @brief Risk-neutral valuation principle
     */
    static std::string principle() {
        return "Risk-Neutral Valuation:\n"
               "\n"
               "Fundamental theorem of asset pricing:\n"
               "Option price = e^(-rT) E^Q[Payoff]\n"
               "\n"
               "where E^Q is expectation under risk-neutral measure Q.\n"
               "\n"
               "Key insights:\n"
               "1. Actual probabilities don't matter!\n"
               "2. Under Q, all assets grow at rate r\n"
               "3. Discounting at r gives arbitrage-free price\n"
               "\n"
               "Risk-neutral world:\n"
               "- Investors are risk-neutral (μ doesn't appear)\n"
               "- Expected return on all assets is r\n"
               "- Stock price: S_T = S_0 e^((r - σ²/2)T + σW_T)\n"
               "\n"
               "Change of measure (Girsanov):\n"
               "- Physical measure P: drift μ\n"
               "- Risk-neutral measure Q: drift r\n"
               "- dW^Q = dW^P + ((μ-r)/σ)dt (market price of risk)\n"
               "\n"
               "Why it works:\n"
               "- Delta hedging eliminates risk\n"
               "- Replicating portfolio grows at r\n"
               "- No-arbitrage forces option price to match";
    }

    /**
     * @brief Martingale property
     */
    static std::string martingale() {
        return "Martingale Property:\n"
               "\n"
               "Under risk-neutral measure Q:\n"
               "  S_t/B_t is a martingale\n"
               "\n"
               "where B_t = e^(rt) is money market account.\n"
               "\n"
               "This means: E^Q[S_T/B_T | F_t] = S_t/B_t\n"
               "\n"
               "Equivalently: e^(-rt)S_t = E^Q[e^(-rT)S_T | F_t]\n"
               "\n"
               "Discounted stock price is martingale under Q!\n"
               "\n"
               "For derivatives:\n"
               "  V_t = e^(-r(T-t)) E^Q[Payoff(S_T) | F_t]\n"
               "\n"
               "This is the pricing formula for all derivatives.";
    }
};

} // namespace maths::finance

#endif // MATHS_FINANCE_BLACK_SCHOLES_HPP
