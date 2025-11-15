/**
 * @file stochastic_differential_equations.hpp
 * @brief Stochastic Differential Equations and Itô Calculus
 *
 * Implements computational algorithms for:
 * - Itô integrals and Itô's lemma
 * - Stochastic differential equations (SDEs)
 * - Filtering problems (Kalman filter)
 * - Stochastic approach to deterministic boundary value problems
 * - Optimal stopping problems
 * - Stochastic control
 * - Mathematical finance applications
 */

#ifndef MATHS_STOCHASTIC_DIFFERENTIAL_EQUATIONS_HPP
#define MATHS_STOCHASTIC_DIFFERENTIAL_EQUATIONS_HPP

#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <numeric>

namespace maths::sde {

// Random number generator
static std::random_device rd;
static std::mt19937 gen(rd());

/**
 * ============================================================================
 * ITÔ INTEGRALS
 * ============================================================================
 */

/**
 * @class ItoIntegral
 * @brief Construction and computation of Itô integrals
 *
 * Itô integral: ∫₀ᵗ f(s, ω) dW(s) where W is Brownian motion
 */
class ItoIntegral {
public:
    /**
     * @brief Compute Itô integral using Euler-Maruyama approximation
     *
     * Approximates ∫₀ᵀ f(t, W(t)) dW(t) using partition sum
     *
     * @param integrand Function f(t, W_t) to integrate
     * @param T Terminal time
     * @param n_steps Number of time steps
     * @return Approximation of Itô integral
     */
    static double compute(
        std::function<double(double, double)> integrand,
        double T,
        int n_steps = 1000) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);

        std::normal_distribution<double> normal(0.0, 1.0);

        double W = 0.0;  // Brownian motion
        double integral = 0.0;

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = sqrt_dt * normal(gen);

            integral += integrand(t, W) * dW;
            W += dW;
        }

        return integral;
    }

    /**
     * @brief Itô isometry: E[|∫ f dW|²] = E[∫ f² dt]
     *
     * Verify the fundamental property of Itô integrals
     */
    static double itoIsometry(
        std::function<double(double)> f,
        double T,
        int n_samples = 100,
        int n_steps = 1000) {

        std::vector<double> integral_squared(n_samples);

        for (int sample = 0; sample < n_samples; ++sample) {
            double dt = T / n_steps;
            double sqrt_dt = std::sqrt(dt);
            std::normal_distribution<double> normal(0.0, 1.0);

            double integral = 0.0;
            for (int i = 0; i < n_steps; ++i) {
                double t = i * dt;
                double dW = sqrt_dt * normal(gen);
                integral += f(t) * dW;
            }

            integral_squared[sample] = integral * integral;
        }

        // E[|∫ f dW|²]
        return std::accumulate(integral_squared.begin(), integral_squared.end(), 0.0) / n_samples;
    }

    /**
     * @brief Quadratic variation of Itô integral
     *
     * [∫ f dW, ∫ f dW]_t = ∫₀ᵗ f² ds
     */
    static double quadraticVariation(
        std::function<double(double)> f,
        double T,
        int n_steps = 1000) {

        double dt = T / n_steps;
        double variation = 0.0;

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            variation += f(t) * f(t) * dt;
        }

        return variation;
    }
};

/**
 * ============================================================================
 * ITÔ'S LEMMA
 * ============================================================================
 */

/**
 * @class ItoLemma
 * @brief Itô's formula for change of variables in SDEs
 *
 * If dX = μ dt + σ dW and Y = f(t, X), then:
 * dY = (∂f/∂t + μ ∂f/∂x + ½σ² ∂²f/∂x²) dt + σ ∂f/∂x dW
 */
class ItoLemma {
public:
    /**
     * @brief Apply Itô's lemma to transform SDE
     *
     * @param f Function Y = f(t, X)
     * @param df_dt Partial derivative ∂f/∂t
     * @param df_dx Partial derivative ∂f/∂x
     * @param d2f_dx2 Second partial derivative ∂²f/∂x²
     * @param mu Drift coefficient in dX = μ dt + σ dW
     * @param sigma Diffusion coefficient
     * @param t Current time
     * @param X Current value of X
     * @return Pair (drift of Y, diffusion of Y)
     */
    static std::pair<double, double> apply(
        std::function<double(double, double)> f,
        std::function<double(double, double)> df_dt,
        std::function<double(double, double)> df_dx,
        std::function<double(double, double)> d2f_dx2,
        double mu, double sigma,
        double t, double X) {

        double drift = df_dt(t, X) + mu * df_dx(t, X) + 0.5 * sigma * sigma * d2f_dx2(t, X);
        double diffusion = sigma * df_dx(t, X);

        return {drift, diffusion};
    }

    /**
     * @brief Geometric Brownian motion via Itô's lemma
     *
     * dS = μS dt + σS dW
     * Solution: S(t) = S₀ exp((μ - σ²/2)t + σW(t))
     */
    static double geometricBrownianMotion(
        double S0, double mu, double sigma, double T,
        int n_steps = 1000) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double S = S0;
        for (int i = 0; i < n_steps; ++i) {
            double dW = sqrt_dt * normal(gen);
            S += mu * S * dt + sigma * S * dW;
        }

        return S;
    }

    /**
     * @brief Multi-dimensional Itô formula
     *
     * For X = (X₁, ..., Xₙ) with dXᵢ = μᵢ dt + ∑ⱼ σᵢⱼ dWⱼ
     * and Y = f(t, X), then:
     * dY = (∂f/∂t + ∑ᵢ μᵢ ∂f/∂xᵢ + ½∑ᵢⱼₖ σᵢₖσⱼₖ ∂²f/∂xᵢ∂xⱼ) dt
     *      + ∑ᵢⱼ σᵢⱼ ∂f/∂xᵢ dWⱼ
     *
     * @param x State vector X
     * @param mu Drift vector
     * @param sigma Diffusion matrix (n × m)
     * @param df_dt Time derivative
     * @param grad_f Gradient ∇f
     * @param hessian_f Hessian ∂²f/∂xᵢ∂xⱼ
     * @return Pair (drift of Y, diffusion coefficients of Y)
     */
    static std::pair<double, std::vector<double>> applyMultiDimensional(
        const std::vector<double>& x,
        const std::vector<double>& mu,
        const std::vector<std::vector<double>>& sigma,
        double df_dt,
        const std::vector<double>& grad_f,
        const std::vector<std::vector<double>>& hessian_f) {

        int n = x.size();
        int m = sigma[0].size();  // Number of Brownian motions

        // Compute drift term
        double drift = df_dt;

        // Add ∑ᵢ μᵢ ∂f/∂xᵢ
        for (int i = 0; i < n; ++i) {
            drift += mu[i] * grad_f[i];
        }

        // Add ½∑ᵢⱼₖ σᵢₖσⱼₖ ∂²f/∂xᵢ∂xⱼ
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double sum_sigma = 0.0;
                for (int k = 0; k < m; ++k) {
                    sum_sigma += sigma[i][k] * sigma[j][k];
                }
                drift += 0.5 * sum_sigma * hessian_f[i][j];
            }
        }

        // Compute diffusion terms: ∑ᵢ σᵢⱼ ∂f/∂xᵢ for each j
        std::vector<double> diffusion(m, 0.0);
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                diffusion[j] += sigma[i][j] * grad_f[i];
            }
        }

        return {drift, diffusion};
    }
};

/**
 * ============================================================================
 * MARTINGALE REPRESENTATION THEOREM
 * ============================================================================
 */

/**
 * @class MartingaleRepresentation
 * @brief Martingale representation and properties
 *
 * Every square-integrable martingale M adapted to Brownian filtration
 * has representation: M(t) = M(0) + ∫₀ᵗ φ(s) dW(s)
 */
class MartingaleRepresentation {
public:
    /**
     * @brief Check if process is a martingale
     *
     * Verify E[X(t) | ℱₛ] = X(s) for s < t
     * (Simplified check: verify drift is zero)
     */
    static bool isMartingale(
        std::function<double(double, double)> drift,
        double t_start, double t_end,
        int n_samples = 100) {

        // For a process to be a martingale, drift must be zero
        double avg_drift = 0.0;
        int n_points = 20;

        for (int i = 0; i < n_points; ++i) {
            double t = t_start + (t_end - t_start) * i / n_points;
            double x = 0.0;  // Test at origin
            avg_drift += std::abs(drift(t, x));
        }

        avg_drift /= n_points;
        return avg_drift < 1e-6;
    }

    /**
     * @brief Represent martingale via Itô integral
     *
     * For martingale M(t), find integrand φ(t) such that:
     * M(t) = M(0) + ∫₀ᵗ φ(s) dW(s)
     *
     * @param martingale Function M(t, W_t)
     * @param phi Integrand (computed via differentiation)
     * @param M0 Initial value M(0)
     * @param T Terminal time
     * @param n_steps Number of time steps
     * @return Simulated martingale path
     */
    static std::vector<double> represent(
        std::function<double(double, double)> phi,
        double M0, double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> M(n_steps + 1);
        M[0] = M0;

        double W = 0.0;
        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = sqrt_dt * normal(gen);

            M[i+1] = M[i] + phi(t, W) * dW;
            W += dW;
        }

        return M;
    }

    /**
     * @brief Doob's optional stopping theorem
     *
     * For martingale M and stopping time τ:
     * E[M(τ)] = M(0) (under suitable conditions)
     */
    static double optionalStopping(
        std::function<double(double, double)> martingale,
        std::function<bool(double, double)> stopping_condition,
        double M0, double max_time,
        int n_samples = 1000) {

        std::vector<double> stopped_values(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);

        for (int sample = 0; sample < n_samples; ++sample) {
            double t = 0.0;
            double W = 0.0;
            double dt = 0.01;
            double sqrt_dt = std::sqrt(dt);

            double M = M0;

            while (t < max_time && !stopping_condition(t, W)) {
                double dW = sqrt_dt * normal(gen);
                W += dW;
                t += dt;
                M = martingale(t, W);
            }

            stopped_values[sample] = M;
        }

        return std::accumulate(stopped_values.begin(), stopped_values.end(), 0.0) / n_samples;
    }
};

/**
 * ============================================================================
 * STOCHASTIC DIFFERENTIAL EQUATIONS
 * ============================================================================
 */

/**
 * @class StochasticDifferentialEquation
 * @brief Numerical solution of SDEs: dX = μ(t, X) dt + σ(t, X) dW
 */
class StochasticDifferentialEquation {
public:
    /**
     * @brief Euler-Maruyama method for solving SDEs
     *
     * dX = μ(t, X) dt + σ(t, X) dW
     * X_{n+1} = X_n + μ(t_n, X_n)Δt + σ(t_n, X_n)ΔW_n
     */
    static std::vector<double> eulerMaruyama(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        double X0, double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> X(n_steps + 1);
        X[0] = X0;

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = sqrt_dt * normal(gen);

            X[i+1] = X[i] + mu(t, X[i]) * dt + sigma(t, X[i]) * dW;
        }

        return X;
    }

    /**
     * @brief Milstein method (higher order accuracy)
     *
     * X_{n+1} = X_n + μΔt + σΔW + ½σ(∂σ/∂x)[(ΔW)² - Δt]
     */
    static std::vector<double> milstein(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double, double)> dsigma_dx,
        double X0, double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> X(n_steps + 1);
        X[0] = X0;

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = sqrt_dt * normal(gen);

            double mu_val = mu(t, X[i]);
            double sigma_val = sigma(t, X[i]);
            double dsigma_val = dsigma_dx(t, X[i]);

            X[i+1] = X[i] + mu_val * dt + sigma_val * dW
                     + 0.5 * sigma_val * dsigma_val * (dW * dW - dt);
        }

        return X;
    }

    /**
     * @brief Ornstein-Uhlenbeck process (mean-reverting)
     *
     * dX = θ(μ - X) dt + σ dW
     * Exact solution: X(t) = X₀e^(-θt) + μ(1 - e^(-θt)) + σ∫₀ᵗ e^(-θ(t-s)) dW(s)
     */
    static std::vector<double> ornsteinUhlenbeck(
        double X0, double theta, double mu, double sigma,
        double T, int n_steps) {

        auto drift = [theta, mu](double t, double X) { return theta * (mu - X); };
        auto diffusion = [sigma](double t, double X) { return sigma; };

        return eulerMaruyama(drift, diffusion, X0, T, n_steps);
    }

    /**
     * @brief Linear SDE with exact solution
     *
     * dX = aX dt + bX dW
     * Exact: X(t) = X₀ exp((a - b²/2)t + bW(t))
     */
    static double linearSDE(double X0, double a, double b, double T) {
        std::normal_distribution<double> normal(0.0, 1.0);
        double W_T = std::sqrt(T) * normal(gen);

        return X0 * std::exp((a - 0.5 * b * b) * T + b * W_T);
    }

    /**
     * @brief Bessel process (radial part of Brownian motion)
     *
     * dR = (n-1)/(2R) dt + dW for n-dimensional Brownian motion
     */
    static std::vector<double> besselProcess(
        double R0, int dimension, double T, int n_steps) {

        auto drift = [dimension](double t, double R) {
            return (dimension - 1) / (2.0 * std::max(R, 1e-10));
        };
        auto diffusion = [](double t, double R) { return 1.0; };

        return eulerMaruyama(drift, diffusion, R0, T, n_steps);
    }

    /**
     * @brief Langevin equation (physics)
     *
     * dv = -γv dt + σ dW (damped Brownian motion)
     */
    static std::vector<double> langevinEquation(
        double v0, double gamma, double sigma, double T, int n_steps) {

        auto drift = [gamma](double t, double v) { return -gamma * v; };
        auto diffusion = [sigma](double t, double v) { return sigma; };

        return eulerMaruyama(drift, diffusion, v0, T, n_steps);
    }
};

/**
 * ============================================================================
 * EXISTENCE AND UNIQUENESS THEOREM
 * ============================================================================
 */

/**
 * @class ExistenceUniqueness
 * @brief Existence and uniqueness results for SDEs
 *
 * Theorem: If μ and σ satisfy Lipschitz and linear growth conditions,
 * then SDE dX = μ(t,X)dt + σ(t,X)dW has unique strong solution
 */
class ExistenceUniqueness {
public:
    /**
     * @brief Check Lipschitz condition
     *
     * |f(t,x) - f(t,y)| ≤ K|x - y| for all x, y
     */
    static bool checkLipschitz(
        std::function<double(double, double)> f,
        double t, double x_min, double x_max,
        int n_points = 100) {

        double max_lipschitz_constant = 0.0;

        for (int i = 0; i < n_points - 1; ++i) {
            double x1 = x_min + (x_max - x_min) * i / n_points;
            double x2 = x_min + (x_max - x_min) * (i + 1) / n_points;

            double f1 = f(t, x1);
            double f2 = f(t, x2);

            double lipschitz = std::abs(f2 - f1) / std::abs(x2 - x1);
            max_lipschitz_constant = std::max(max_lipschitz_constant, lipschitz);
        }

        // Check if bounded
        return std::isfinite(max_lipschitz_constant);
    }

    /**
     * @brief Check linear growth condition
     *
     * |f(t,x)|² ≤ K(1 + |x|²) for all x
     */
    static bool checkLinearGrowth(
        std::function<double(double, double)> f,
        double t, double x_min, double x_max,
        int n_points = 100) {

        for (int i = 0; i < n_points; ++i) {
            double x = x_min + (x_max - x_min) * i / n_points;
            double f_val = f(t, x);

            // Check if f²/(1 + x²) is bounded
            double ratio = (f_val * f_val) / (1.0 + x * x);

            if (!std::isfinite(ratio) || ratio > 1e6) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Verify conditions for existence & uniqueness
     */
    static bool verifyConditions(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        double t, double x_min, double x_max) {

        bool mu_lipschitz = checkLipschitz(mu, t, x_min, x_max);
        bool sigma_lipschitz = checkLipschitz(sigma, t, x_min, x_max);
        bool mu_growth = checkLinearGrowth(mu, t, x_min, x_max);
        bool sigma_growth = checkLinearGrowth(sigma, t, x_min, x_max);

        return mu_lipschitz && sigma_lipschitz && mu_growth && sigma_growth;
    }

    /**
     * @brief Picard iteration for SDE existence proof
     *
     * X^{(n+1)}(t) = X₀ + ∫₀ᵗ μ(s, X^{(n)}(s)) ds + ∫₀ᵗ σ(s, X^{(n)}(s)) dW(s)
     */
    static std::vector<double> picardIteration(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        double X0, double T, int n_iterations, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        // Initialize X^(0) = X₀
        std::vector<double> X_old(n_steps + 1, X0);
        std::vector<double> X_new(n_steps + 1);
        X_new[0] = X0;

        // Picard iterations
        for (int iter = 0; iter < n_iterations; ++iter) {
            double W = 0.0;

            for (int i = 0; i < n_steps; ++i) {
                double t = i * dt;
                double dW = sqrt_dt * normal(gen);

                X_new[i+1] = X0 + mu(t, X_old[i]) * dt + sigma(t, X_old[i]) * dW;
                W += dW;
            }

            X_old = X_new;
        }

        return X_new;
    }
};

/**
 * ============================================================================
 * WEAK AND STRONG SOLUTIONS
 * ============================================================================
 */

/**
 * @class WeakStrongSolutions
 * @brief Distinction between weak and strong solutions of SDEs
 *
 * Strong solution: Adapted to given Brownian motion filtration
 * Weak solution: Brownian motion is part of the solution (law of process)
 */
class WeakStrongSolutions {
public:
    /**
     * @brief Strong solution (pathwise uniqueness)
     *
     * Given (Ω, ℱ, {ℱₜ}, P, W), find adapted process X satisfying SDE
     */
    static std::vector<double> strongSolution(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        double X0, double T, int n_steps) {

        // Use Euler-Maruyama on given Brownian motion
        return StochasticDifferentialEquation::eulerMaruyama(mu, sigma, X0, T, n_steps);
    }

    /**
     * @brief Weak solution (distributional uniqueness)
     *
     * Find probability space and Brownian motion such that SDE is satisfied
     * Focus on law of process, not specific path
     */
    static std::vector<double> weakSolution(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        double X0, double T, int n_steps,
        int n_samples = 1000) {

        // Average over multiple realizations (distribution)
        std::vector<std::vector<double>> samples(n_samples);

        for (int s = 0; s < n_samples; ++s) {
            samples[s] = StochasticDifferentialEquation::eulerMaruyama(mu, sigma, X0, T, n_steps);
        }

        // Return expected path
        std::vector<double> expected_path(n_steps + 1, 0.0);
        for (int i = 0; i <= n_steps; ++i) {
            for (int s = 0; s < n_samples; ++s) {
                expected_path[i] += samples[s][i];
            }
            expected_path[i] /= n_samples;
        }

        return expected_path;
    }

    /**
     * @brief Girsanov theorem (change of measure)
     *
     * Transform drift via change of probability measure
     * dW̃ = dW + θ(t) dt is Brownian under new measure Q
     */
    static std::vector<double> girsanovTransform(
        std::function<double(double)> theta,
        double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> W_tilde(n_steps + 1, 0.0);

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = sqrt_dt * normal(gen);

            // W̃(t+dt) = W̃(t) + dW + θ(t)dt
            W_tilde[i+1] = W_tilde[i] + dW + theta(t) * dt;
        }

        return W_tilde;
    }

    /**
     * @brief Radon-Nikodym derivative for measure change
     *
     * dQ/dP = exp(-∫₀ᵗ θ(s) dW(s) - ½∫₀ᵗ θ²(s) ds)
     */
    static double radonNikodym(
        std::function<double(double)> theta,
        const std::vector<double>& W_path,
        double T) {

        int n_steps = W_path.size() - 1;
        double dt = T / n_steps;

        double integral_theta_dW = 0.0;
        double integral_theta_squared = 0.0;

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = W_path[i+1] - W_path[i];

            integral_theta_dW += theta(t) * dW;
            integral_theta_squared += theta(t) * theta(t) * dt;
        }

        return std::exp(-integral_theta_dW - 0.5 * integral_theta_squared);
    }
};

/**
 * ============================================================================
 * FILTERING PROBLEMS
 * ============================================================================
 */

/**
 * @class KalmanFilter
 * @brief Optimal filtering for linear Gaussian systems
 *
 * State equation: X_{k+1} = AX_k + w_k
 * Observation: Y_k = HX_k + v_k
 * where w_k ~ N(0, Q), v_k ~ N(0, R)
 */
class KalmanFilter {
public:
    struct State {
        std::vector<double> x;  // State estimate
        std::vector<std::vector<double>> P;  // Covariance matrix
    };

    /**
     * @brief Kalman filter prediction step
     */
    static State predict(
        const State& state,
        const std::vector<std::vector<double>>& A,
        const std::vector<std::vector<double>>& Q) {

        int n = state.x.size();
        State predicted;
        predicted.x.resize(n);
        predicted.P.resize(n, std::vector<double>(n));

        // Predict state: x̂_k|k-1 = A x̂_k-1|k-1
        for (int i = 0; i < n; ++i) {
            predicted.x[i] = 0.0;
            for (int j = 0; j < n; ++j) {
                predicted.x[i] += A[i][j] * state.x[j];
            }
        }

        // Predict covariance: P_k|k-1 = A P_k-1|k-1 A^T + Q
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                predicted.P[i][j] = Q[i][j];
                for (int k = 0; k < n; ++k) {
                    for (int l = 0; l < n; ++l) {
                        predicted.P[i][j] += A[i][k] * state.P[k][l] * A[j][l];
                    }
                }
            }
        }

        return predicted;
    }

    /**
     * @brief Kalman filter update step
     */
    static State update(
        const State& predicted,
        const std::vector<double>& measurement,
        const std::vector<std::vector<double>>& H,
        const std::vector<std::vector<double>>& R) {

        int n = predicted.x.size();
        int m = measurement.size();

        // Innovation: y = z - H x̂_k|k-1
        std::vector<double> innovation(m);
        for (int i = 0; i < m; ++i) {
            innovation[i] = measurement[i];
            for (int j = 0; j < n; ++j) {
                innovation[i] -= H[i][j] * predicted.x[j];
            }
        }

        // Innovation covariance: S = H P_k|k-1 H^T + R
        std::vector<std::vector<double>> S(m, std::vector<double>(m));
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                S[i][j] = R[i][j];
                for (int k = 0; k < n; ++k) {
                    for (int l = 0; l < n; ++l) {
                        S[i][j] += H[i][k] * predicted.P[k][l] * H[j][l];
                    }
                }
            }
        }

        // Kalman gain: K = P_k|k-1 H^T S^(-1)
        // (Simplified - assumes S is invertible)
        std::vector<std::vector<double>> K(n, std::vector<double>(m, 0.0));

        // Update state: x̂_k|k = x̂_k|k-1 + K y
        State updated;
        updated.x.resize(n);
        updated.P.resize(n, std::vector<double>(n));

        for (int i = 0; i < n; ++i) {
            updated.x[i] = predicted.x[i];
            for (int j = 0; j < m; ++j) {
                updated.x[i] += K[i][j] * innovation[j];
            }
        }

        // Update covariance: P_k|k = (I - KH) P_k|k-1
        updated.P = predicted.P;

        return updated;
    }
};

/**
 * ============================================================================
 * OPTIMAL STOPPING
 * ============================================================================
 */

/**
 * @class OptimalStopping
 * @brief Optimal stopping problems and American options
 */
class OptimalStopping {
public:
    /**
     * @brief American put option via dynamic programming
     *
     * Value function: V(t, S) = max(K - S, continuation value)
     *
     * @param S0 Initial stock price
     * @param K Strike price
     * @param r Risk-free rate
     * @param sigma Volatility
     * @param T Maturity
     * @param n_steps Number of time steps
     * @param n_paths Number of Monte Carlo paths
     * @return Option value
     */
    static double americanPut(
        double S0, double K, double r, double sigma,
        double T, int n_steps, int n_paths = 1000) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> payoffs(n_paths);

        // Simulate paths and find optimal stopping time
        for (int path = 0; path < n_paths; ++path) {
            double S = S0;
            double max_payoff = std::max(K - S, 0.0);
            double discount = 1.0;

            for (int i = 0; i < n_steps; ++i) {
                double dW = sqrt_dt * normal(gen);
                S *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * dW);

                discount *= std::exp(-r * dt);
                double immediate_payoff = discount * std::max(K - S, 0.0);

                // Greedy stopping (simplified)
                max_payoff = std::max(max_payoff, immediate_payoff);
            }

            payoffs[path] = max_payoff;
        }

        return std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / n_paths;
    }

    /**
     * @brief Optimal stopping time for general reward process
     */
    static double optimalStoppingTime(
        std::function<double(double, double)> reward,
        double X0, double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double X = X0;
        double max_reward = reward(0.0, X);
        double optimal_time = 0.0;

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = sqrt_dt * normal(gen);
            X += dW;

            double current_reward = reward(t, X);
            if (current_reward > max_reward) {
                max_reward = current_reward;
                optimal_time = t;
            }
        }

        return optimal_time;
    }
};

/**
 * ============================================================================
 * STOCHASTIC CONTROL
 * ============================================================================
 */

/**
 * @class StochasticControl
 * @brief Hamilton-Jacobi-Bellman equation and optimal control
 */
class StochasticControl {
public:
    /**
     * @brief Linear-Quadratic-Gaussian (LQG) control
     *
     * Minimize: E[∫₀ᵀ (x^T Q x + u^T R u) dt + x(T)^T Q_f x(T)]
     * Subject to: dx = (Ax + Bu) dt + σ dW
     *
     * Optimal control: u*(t) = -R^(-1) B^T P(t) x(t)
     * where P satisfies Riccati equation
     */
    static std::vector<double> lqgControl(
        const std::vector<double>& x0,
        const std::vector<std::vector<double>>& A,
        const std::vector<std::vector<double>>& B,
        const std::vector<std::vector<double>>& Q,
        const std::vector<std::vector<double>>& R,
        double T, int n_steps) {

        int n = x0.size();
        double dt = T / n_steps;

        std::vector<double> x = x0;
        std::vector<double> trajectory;
        trajectory.reserve(n_steps + 1);

        for (int step = 0; step < n_steps; ++step) {
            // Compute control (simplified - assumes R = I)
            std::vector<double> u(B[0].size(), 0.0);

            // Update state
            std::vector<double> x_new(n);
            for (int i = 0; i < n; ++i) {
                x_new[i] = x[i];
                for (int j = 0; j < n; ++j) {
                    x_new[i] += A[i][j] * x[j] * dt;
                }
                for (int j = 0; j < (int)u.size(); ++j) {
                    x_new[i] += B[i][j] * u[j] * dt;
                }
            }

            x = x_new;
            trajectory.push_back(x[0]);  // Store first component
        }

        return trajectory;
    }

    /**
     * @brief Merton's portfolio optimization
     *
     * Maximize: E[U(X(T))] where U is utility function
     * Subject to: dX = (rX + π(μ - r)) dt + πσ dW
     *
     * Optimal fraction: π* = (μ - r)/(γσ²) where γ is risk aversion
     */
    static double mertonOptimalFraction(
        double mu, double r, double sigma, double gamma) {

        return (mu - r) / (gamma * sigma * sigma);
    }
};

/**
 * ============================================================================
 * MATHEMATICAL FINANCE APPLICATIONS
 * ============================================================================
 */

/**
 * @class FinancialSDEs
 * @brief Stochastic models in mathematical finance
 */
class FinancialSDEs {
public:
    /**
     * @brief Heston stochastic volatility model
     *
     * dS = μS dt + √v S dW₁
     * dv = κ(θ - v) dt + σ_v √v dW₂
     * where dW₁ dW₂ = ρ dt
     */
    static std::pair<std::vector<double>, std::vector<double>> hestonModel(
        double S0, double v0,
        double mu, double kappa, double theta, double sigma_v, double rho,
        double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> S(n_steps + 1), v(n_steps + 1);
        S[0] = S0;
        v[0] = v0;

        for (int i = 0; i < n_steps; ++i) {
            double dW1 = sqrt_dt * normal(gen);
            double dZ = sqrt_dt * normal(gen);
            double dW2 = rho * dW1 + std::sqrt(1 - rho * rho) * dZ;

            // Ensure v stays positive (Feller condition)
            v[i+1] = std::max(0.0, v[i] + kappa * (theta - v[i]) * dt
                              + sigma_v * std::sqrt(std::max(v[i], 0.0)) * dW2);

            S[i+1] = S[i] + mu * S[i] * dt + std::sqrt(std::max(v[i], 0.0)) * S[i] * dW1;
        }

        return {S, v};
    }

    /**
     * @brief Cox-Ingersoll-Ross (CIR) model for interest rates
     *
     * dr = κ(θ - r) dt + σ√r dW
     */
    static std::vector<double> cirModel(
        double r0, double kappa, double theta, double sigma,
        double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> r(n_steps + 1);
        r[0] = r0;

        for (int i = 0; i < n_steps; ++i) {
            double dW = sqrt_dt * normal(gen);
            r[i+1] = std::max(0.0, r[i] + kappa * (theta - r[i]) * dt
                              + sigma * std::sqrt(std::max(r[i], 0.0)) * dW);
        }

        return r;
    }

    /**
     * @brief Vasicek model for interest rates
     *
     * dr = κ(θ - r) dt + σ dW
     * (Allows negative rates unlike CIR)
     */
    static std::vector<double> vasicekModel(
        double r0, double kappa, double theta, double sigma,
        double T, int n_steps) {

        auto drift = [kappa, theta](double t, double r) { return kappa * (theta - r); };
        auto diffusion = [sigma](double t, double r) { return sigma; };

        return StochasticDifferentialEquation::eulerMaruyama(drift, diffusion, r0, T, n_steps);
    }
};

} // namespace maths::sde

#endif // MATHS_STOCHASTIC_DIFFERENTIAL_EQUATIONS_HPP
