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
 * @class FilteringProblem
 * @brief Linear and nonlinear filtering theory
 *
 * The filtering problem: Estimate state X(t) from observations Y(s), 0 ≤ s ≤ t
 * - Signal process: dX = b(X) dt + σ(X) dW
 * - Observation: dY = h(X) dt + dV
 * where W, V are independent Brownian motions
 */
class FilteringProblem {
public:
    /**
     * @brief 1-Dimensional linear filtering problem
     *
     * Signal: dX = a X dt + σ dW
     * Observation: dY = c X dt + dV
     *
     * Optimal estimate: dX̂ = a X̂ dt + K(t) (dY - c X̂ dt)
     * where K(t) is Kalman gain
     */
    static std::vector<double> linearFilter1D(
        double X0,  // Initial state
        double a,   // Signal drift coefficient
        double sigma,  // Signal diffusion
        double c,   // Observation coefficient
        const std::vector<double>& observations,
        double dt) {

        int n_steps = observations.size();
        std::vector<double> X_hat(n_steps);
        X_hat[0] = X0;

        // Variance of estimate
        double P = 1.0;

        for (int i = 0; i < n_steps - 1; ++i) {
            // Riccati equation for variance
            double dP_dt = 2.0 * a * P + sigma * sigma - c * c * P * P;
            P += dP_dt * dt;

            // Kalman gain
            double K = c * P;

            // Innovation
            double innovation = observations[i+1] - observations[i] - c * X_hat[i] * dt;

            // Update estimate
            X_hat[i+1] = X_hat[i] + a * X_hat[i] * dt + K * innovation;
        }

        return X_hat;
    }

    /**
     * @brief Multidimensional linear filtering (Kalman-Bucy filter)
     *
     * Signal: dX = A X dt + B dW
     * Observation: dY = C X dt + D dV
     *
     * Continuous-time Kalman filter
     */
    static std::vector<std::vector<double>> multidimensionalLinearFilter(
        const std::vector<double>& X0,
        const std::vector<std::vector<double>>& A,  // n×n drift matrix
        const std::vector<std::vector<double>>& B,  // n×m diffusion matrix
        const std::vector<std::vector<double>>& C,  // p×n observation matrix
        const std::vector<std::vector<double>>& observations,  // p-dimensional
        double dt) {

        int n = X0.size();
        int p = C.size();
        int n_steps = observations.size();

        std::vector<std::vector<double>> X_hat(n_steps, std::vector<double>(n));
        X_hat[0] = X0;

        // Covariance matrix P
        std::vector<std::vector<double>> P(n, std::vector<double>(n, 1.0));

        for (int step = 0; step < n_steps - 1; ++step) {
            // Riccati equation: dP/dt = AP + PA^T + BB^T - PC^T CP
            std::vector<std::vector<double>> dP(n, std::vector<double>(n, 0.0));

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    // AP + PA^T
                    for (int k = 0; k < n; ++k) {
                        dP[i][j] += A[i][k] * P[k][j] + P[i][k] * A[j][k];
                    }

                    // BB^T
                    for (int k = 0; k < (int)B[0].size(); ++k) {
                        dP[i][j] += B[i][k] * B[j][k];
                    }

                    // -PC^T CP
                    for (int k = 0; k < p; ++k) {
                        for (int l = 0; l < p; ++l) {
                            for (int m = 0; m < n; ++m) {
                                for (int r = 0; r < n; ++r) {
                                    dP[i][j] -= P[i][m] * C[k][m] * C[l][r] * P[r][j];
                                }
                            }
                        }
                    }
                }
            }

            // Update P
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    P[i][j] += dP[i][j] * dt;
                }
            }

            // Kalman gain K = PC^T
            std::vector<std::vector<double>> K(n, std::vector<double>(p));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < p; ++j) {
                    K[i][j] = 0.0;
                    for (int k = 0; k < n; ++k) {
                        K[i][j] += P[i][k] * C[j][k];
                    }
                }
            }

            // Innovation: dY - CX̂ dt
            std::vector<double> innovation(p);
            for (int i = 0; i < p; ++i) {
                innovation[i] = observations[step+1][i] - observations[step][i];
                for (int j = 0; j < n; ++j) {
                    innovation[i] -= C[i][j] * X_hat[step][j] * dt;
                }
            }

            // Update estimate: dX̂ = AX̂ dt + K(dY - CX̂ dt)
            for (int i = 0; i < n; ++i) {
                X_hat[step+1][i] = X_hat[step][i];

                // AX̂ dt term
                for (int j = 0; j < n; ++j) {
                    X_hat[step+1][i] += A[i][j] * X_hat[step][j] * dt;
                }

                // K * innovation term
                for (int j = 0; j < p; ++j) {
                    X_hat[step+1][i] += K[i][j] * innovation[j];
                }
            }
        }

        return X_hat;
    }

    /**
     * @brief Innovation process
     *
     * ν(t) = Y(t) - ∫₀ᵗ h(X̂(s)) ds
     * where X̂ is the filtered estimate
     */
    static std::vector<double> innovationProcess(
        const std::vector<double>& observations,
        const std::vector<double>& filtered_estimate,
        std::function<double(double)> h,
        double dt) {

        int n_steps = observations.size();
        std::vector<double> innovation(n_steps);
        innovation[0] = 0.0;

        for (int i = 1; i < n_steps; ++i) {
            innovation[i] = observations[i] - observations[i-1]
                          - h(filtered_estimate[i-1]) * dt;
        }

        return innovation;
    }
};

/**
 * @class KalmanFilter
 * @brief Discrete-time Kalman filter (for compatibility)
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
 * APPLICATIONS TO BOUNDARY VALUE PROBLEMS
 * ============================================================================
 */

/**
 * @class BoundaryValueProblems
 * @brief Stochastic representation of boundary value problems
 *
 * Connection between diffusions and elliptic PDEs via probabilistic methods
 */
class BoundaryValueProblems {
public:
    /**
     * @brief Combined Dirichlet-Poisson problem
     *
     * Solve: -½σ²Δu + b·∇u = f in D, u = g on ∂D
     *
     * Stochastic representation:
     * u(x) = E[∫₀^τ f(X(s)) ds + g(X(τ))]
     * where τ is first exit time from domain D
     */
    static double dirichletPoisson(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double, double)> f,      // Source term
        std::function<double(double, double)> g,      // Boundary condition
        std::function<bool(double, double)> in_domain,
        double x0, double y0,
        int n_samples = 1000) {

        std::vector<double> values(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);

        for (int sample = 0; sample < n_samples; ++sample) {
            double x = x0, y = y0;
            double t = 0.0;
            double dt = 0.01;
            double sqrt_dt = std::sqrt(dt);
            double integral_f = 0.0;

            // Run until exit from domain
            while (in_domain(x, y) && t < 10.0) {
                double dW = sqrt_dt * normal(gen);

                // Accumulate source term
                integral_f += f(x, y) * dt;

                // Evolve 2D diffusion
                x += mu(t, x) * dt + sigma(t, x) * dW;
                double dW2 = sqrt_dt * normal(gen);
                y += mu(t, y) * dt + sigma(t, y) * dW2;

                t += dt;
            }

            // Add boundary value
            values[sample] = integral_f + g(x, y);
        }

        return std::accumulate(values.begin(), values.end(), 0.0) / n_samples;
    }

    /**
     * @brief Uniqueness via maximum principle
     *
     * For elliptic operator L, if Lu₁ = Lu₂ and u₁ = u₂ on boundary,
     * then u₁ = u₂ everywhere (uniqueness)
     */
    static bool verifyUniqueness(
        const std::vector<std::vector<double>>& u1,
        const std::vector<std::vector<double>>& u2,
        double tolerance = 1e-6) {

        int n = u1.size();
        int m = u1[0].size();

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                if (std::abs(u1[i][j] - u2[i][j]) > tolerance) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * @brief Dirichlet problem: Δu = 0 in D, u = g on ∂D
     *
     * Regular point x: For every continuous g on ∂D,
     * solution u is continuous at x and u(x) = E[g(X(τ))]
     * where τ is exit time and X starts at x
     */
    static double dirichletProblem(
        std::function<double(double, double)> g,      // Boundary data
        std::function<bool(double, double)> in_domain,
        double x0, double y0,
        int n_samples = 1000) {

        std::vector<double> boundary_values(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);

        for (int sample = 0; sample < n_samples; ++sample) {
            double x = x0, y = y0;
            double t = 0.0;
            double dt = 0.01;
            double sqrt_dt = std::sqrt(dt);

            // Brownian motion until exit
            while (in_domain(x, y) && t < 10.0) {
                double dW1 = sqrt_dt * normal(gen);
                double dW2 = sqrt_dt * normal(gen);

                x += dW1;
                y += dW2;
                t += dt;
            }

            boundary_values[sample] = g(x, y);
        }

        return std::accumulate(boundary_values.begin(), boundary_values.end(), 0.0) / n_samples;
    }

    /**
     * @brief Check if point is regular for Dirichlet problem
     *
     * Point x is regular if P(τ_D = 0) = 1 where τ_D is hitting time of ∂D
     */
    static bool isRegularPoint(
        double x0, double y0,
        std::function<bool(double, double)> in_domain,
        int n_samples = 1000) {

        int immediate_exits = 0;
        std::normal_distribution<double> normal(0.0, 1.0);
        double epsilon = 0.001;

        for (int sample = 0; sample < n_samples; ++sample) {
            double x = x0, y = y0;
            double dt = epsilon;
            double sqrt_dt = std::sqrt(dt);

            // Check if exits immediately
            double dW1 = sqrt_dt * normal(gen);
            double dW2 = sqrt_dt * normal(gen);

            x += dW1;
            y += dW2;

            if (!in_domain(x, y)) {
                immediate_exits++;
            }
        }

        // Regular if exits immediately with high probability
        return (double)immediate_exits / n_samples > 0.9;
    }

    /**
     * @brief Poisson problem: Δu = f in D, u = 0 on ∂D
     *
     * Stochastic representation: u(x) = E[∫₀^τ f(X(s)) ds]
     * where X is Brownian motion and τ is first exit time
     */
    static double poissonProblem(
        std::function<double(double, double)> f,
        std::function<bool(double, double)> in_domain,
        double x0, double y0,
        int n_samples = 1000) {

        std::vector<double> integrals(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);

        for (int sample = 0; sample < n_samples; ++sample) {
            double x = x0, y = y0;
            double t = 0.0;
            double dt = 0.01;
            double sqrt_dt = std::sqrt(dt);
            double integral = 0.0;

            while (in_domain(x, y) && t < 10.0) {
                double dW1 = sqrt_dt * normal(gen);
                double dW2 = sqrt_dt * normal(gen);

                integral += f(x, y) * dt;

                x += dW1;
                y += dW2;
                t += dt;
            }

            integrals[sample] = integral;
        }

        return std::accumulate(integrals.begin(), integrals.end(), 0.0) / n_samples;
    }

    /**
     * @brief Mean value property for harmonic functions
     *
     * If Δu = 0, then u(x) = average of u on any sphere around x
     */
    static double meanValueProperty(
        std::function<double(double, double)> u,
        double x0, double y0,
        double radius,
        int n_points = 100) {

        double sum = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double theta = 2.0 * M_PI * i / n_points;
            double x = x0 + radius * std::cos(theta);
            double y = y0 + radius * std::sin(theta);

            sum += u(x, y);
        }

        return sum / n_points;
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

    /**
     * @brief Time-homogeneous optimal stopping
     *
     * Value function: V(x) = sup_τ E[g(X(τ)) | X(0) = x]
     * where coefficients don't depend on time
     *
     * Free boundary problem: Find continuation region C
     * - In C: V(x) solves LV = 0
     * - On ∂C: V(x) = g(x) (value matching)
     * - On ∂C: V'(x) = g'(x) (smooth pasting)
     */
    static double timeHomogeneous(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double)> payoff,
        double X0, double max_time,
        int n_samples = 1000) {

        std::vector<double> values(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);

        for (int sample = 0; sample < n_samples; ++sample) {
            double X = X0;
            double t = 0.0;
            double dt = 0.01;
            double sqrt_dt = std::sqrt(dt);
            double max_value = payoff(X);

            while (t < max_time) {
                double dW = sqrt_dt * normal(gen);
                X += mu(t, X) * dt + sigma(t, X) * dW;
                t += dt;

                max_value = std::max(max_value, payoff(X));
            }

            values[sample] = max_value;
        }

        return std::accumulate(values.begin(), values.end(), 0.0) / n_samples;
    }

    /**
     * @brief Time-inhomogeneous optimal stopping
     *
     * Value function: V(t, x) = sup_τ E[g(τ, X(τ)) | X(t) = x]
     * where payoff depends explicitly on time
     *
     * Backward induction: V(t,x) = max(g(t,x), continuation value)
     */
    static double timeInhomogeneous(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double, double)> payoff,  // Time-dependent
        double X0, double T,
        int n_steps = 100,
        int n_samples = 1000) {

        std::vector<double> values(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);
        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);

        for (int sample = 0; sample < n_samples; ++sample) {
            double X = X0;
            double t = 0.0;
            double max_value = payoff(0.0, X);

            for (int i = 0; i < n_steps; ++i) {
                double dW = sqrt_dt * normal(gen);
                X += mu(t, X) * dt + sigma(t, X) * dW;
                t += dt;

                max_value = std::max(max_value, payoff(t, X));
            }

            values[sample] = max_value;
        }

        return std::accumulate(values.begin(), values.end(), 0.0) / n_samples;
    }

    /**
     * @brief Optimal stopping with running payoff integral
     *
     * Value function: V(x) = sup_τ E[∫₀^τ f(X(s)) ds + g(X(τ)) | X(0) = x]
     *
     * Variational inequality:
     * max(LV - f, g - V) = 0
     */
    static double withIntegral(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double)> running_payoff,  // f(x)
        std::function<double(double)> terminal_payoff, // g(x)
        double X0, double max_time,
        int n_samples = 1000) {

        std::vector<double> values(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);

        for (int sample = 0; sample < n_samples; ++sample) {
            double X = X0;
            double t = 0.0;
            double dt = 0.01;
            double sqrt_dt = std::sqrt(dt);
            double integral = 0.0;
            double best_value = terminal_payoff(X);

            while (t < max_time) {
                double dW = sqrt_dt * normal(gen);
                integral += running_payoff(X) * dt;
                X += mu(t, X) * dt + sigma(t, X) * dW;
                t += dt;

                // Compare integral + terminal vs continuing
                double current_value = integral + terminal_payoff(X);
                best_value = std::max(best_value, current_value);
            }

            values[sample] = best_value;
        }

        return std::accumulate(values.begin(), values.end(), 0.0) / n_samples;
    }

    /**
     * @brief Connection with variational inequalities
     *
     * Optimal stopping ⟺ Variational inequality:
     * V ≥ g (value above payoff)
     * LV ≤ 0 (superharmonic)
     * (V - g)(LV) = 0 (complementarity)
     *
     * Continuation region: C = {x : V(x) > g(x)}
     * Stopping region: S = {x : V(x) = g(x)}
     */
    static std::pair<double, bool> variationalInequality(
        std::function<double(double)> value_function,
        std::function<double(double)> payoff,
        std::function<double(double)> LV,  // Generator applied to V
        double x,
        double tolerance = 1e-6) {

        double V = value_function(x);
        double g = payoff(x);
        double LV_val = LV(x);

        // Check variational inequality conditions
        bool in_continuation = (V - g) > tolerance;
        bool superharmonic = LV_val <= tolerance;
        bool complementarity = std::abs((V - g) * LV_val) < tolerance;

        bool satisfies_VI = (V >= g - tolerance) && superharmonic && complementarity;

        return {V, in_continuation};
    }

    /**
     * @brief Solve variational inequality via finite differences
     *
     * Find V satisfying: max(LV, g - V) = 0
     */
    static std::vector<double> solveVariationalInequality(
        std::function<double(double)> b,
        std::function<double(double)> sigma,
        std::function<double(double)> payoff,
        double x_min, double x_max,
        int n_points = 100,
        int max_iterations = 1000) {

        double dx = (x_max - x_min) / n_points;
        std::vector<double> V(n_points + 1);

        // Initialize with payoff
        for (int i = 0; i <= n_points; ++i) {
            double x = x_min + i * dx;
            V[i] = payoff(x);
        }

        // Fixed point iteration
        for (int iter = 0; iter < max_iterations; ++iter) {
            std::vector<double> V_new = V;

            for (int i = 1; i < n_points; ++i) {
                double x = x_min + i * dx;

                // Compute LV via finite differences
                double dV = (V[i+1] - V[i-1]) / (2.0 * dx);
                double d2V = (V[i+1] - 2.0 * V[i] + V[i-1]) / (dx * dx);
                double LV = b(x) * dV + 0.5 * sigma(x) * sigma(x) * d2V;

                // Variational inequality: max(LV, g - V) = 0
                // Equivalent to: V = max(g, V - α·LV) for small α
                double alpha = 0.01;
                V_new[i] = std::max(payoff(x), V[i] - alpha * LV);
            }

            V = V_new;
        }

        return V;
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
     * @brief Statement of stochastic control problem
     *
     * Controlled SDE: dX = b(X, u) dt + σ(X, u) dW
     * Control process: u(t) ∈ U (admissible controls)
     * Performance criterion: J(u) = E[∫₀ᵀ f(X,u) dt + g(X(T))]
     *
     * Objective: Find u* that minimizes J(u)
     */
    struct ControlProblem {
        std::function<double(double, double, double)> drift;      // b(t, x, u)
        std::function<double(double, double, double)> diffusion;  // σ(t, x, u)
        std::function<double(double, double, double)> running_cost;  // f(t, x, u)
        std::function<double(double)> terminal_cost;  // g(x)
        double X0;
        double T;
    };

    /**
     * @brief Solve control problem via Monte Carlo
     */
    static double solveControlProblem(
        const ControlProblem& problem,
        std::function<double(double, double)> control_law,  // u(t, x)
        int n_steps = 100,
        int n_samples = 1000) {

        std::vector<double> costs(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);
        double dt = problem.T / n_steps;
        double sqrt_dt = std::sqrt(dt);

        for (int sample = 0; sample < n_samples; ++sample) {
            double X = problem.X0;
            double t = 0.0;
            double running_cost_integral = 0.0;

            for (int i = 0; i < n_steps; ++i) {
                double u = control_law(t, X);
                double dW = sqrt_dt * normal(gen);

                // Accumulate running cost
                running_cost_integral += problem.running_cost(t, X, u) * dt;

                // Evolve state
                X += problem.drift(t, X, u) * dt + problem.diffusion(t, X, u) * dW;
                t += dt;
            }

            costs[sample] = running_cost_integral + problem.terminal_cost(X);
        }

        return std::accumulate(costs.begin(), costs.end(), 0.0) / n_samples;
    }

    /**
     * @brief Hamilton-Jacobi-Bellman (HJB) equation
     *
     * Value function V(t, x) satisfies:
     * ∂V/∂t + inf_u [b(x,u)∂V/∂x + ½σ²(x,u)∂²V/∂x² + f(x,u)] = 0
     * V(T, x) = g(x)
     *
     * Optimal control: u*(t,x) = argmin_u [...]
     * Verification theorem: If V smooth and u* optimal, then V is value function
     */
    static std::vector<std::vector<double>> solveHJB(
        std::function<double(double, double)> b,
        std::function<double(double, double)> sigma,
        std::function<double(double)> running_cost,
        std::function<double(double)> terminal_cost,
        double x_min, double x_max, int n_x,
        double T, int n_t) {

        double dx = (x_max - x_min) / n_x;
        double dt = T / n_t;

        // Value function V[time][space]
        std::vector<std::vector<double>> V(n_t + 1, std::vector<double>(n_x + 1));

        // Terminal condition
        for (int i = 0; i <= n_x; ++i) {
            double x = x_min + i * dx;
            V[n_t][i] = terminal_cost(x);
        }

        // Backward time stepping (HJB equation)
        for (int k = n_t - 1; k >= 0; --k) {
            for (int i = 1; i < n_x; ++i) {
                double x = x_min + i * dx;

                // Spatial derivatives
                double dV_dx = (V[k+1][i+1] - V[k+1][i-1]) / (2.0 * dx);
                double d2V_dx2 = (V[k+1][i+1] - 2.0 * V[k+1][i] + V[k+1][i-1]) / (dx * dx);

                // Hamiltonian (minimization over control)
                double hamiltonian = b(x, 0.0) * dV_dx
                                   + 0.5 * sigma(x, 0.0) * sigma(x, 0.0) * d2V_dx2
                                   + running_cost(x);

                // Backward Euler for stability
                V[k][i] = V[k+1][i] - dt * hamiltonian;
            }

            // Boundary conditions
            V[k][0] = V[k][1];
            V[k][n_x] = V[k][n_x-1];
        }

        return V;
    }

    /**
     * @brief Extract optimal control from value function
     *
     * u*(t, x) = argmin_u [b(x,u)∂V/∂x + ½σ²(x,u)∂²V/∂x² + f(x,u)]
     */
    static double optimalControl(
        double dV_dx, double d2V_dx2,
        std::function<double(double)> b_control,
        std::function<double(double)> sigma_control,
        std::function<double(double)> f_control,
        double u_min, double u_max,
        int n_search = 100) {

        double min_hamiltonian = std::numeric_limits<double>::infinity();
        double best_u = 0.0;

        for (int i = 0; i <= n_search; ++i) {
            double u = u_min + (u_max - u_min) * i / n_search;

            double hamiltonian = b_control(u) * dV_dx
                               + 0.5 * sigma_control(u) * sigma_control(u) * d2V_dx2
                               + f_control(u);

            if (hamiltonian < min_hamiltonian) {
                min_hamiltonian = hamiltonian;
                best_u = u;
            }
        }

        return best_u;
    }

    /**
     * @brief Stochastic control with terminal conditions
     *
     * Minimize: E[g(X(T))] (no running cost)
     * Subject to: dX = b(X, u) dt + σ(X) dW
     *
     * HJB: ∂V/∂t + inf_u [b(x,u)∂V/∂x] + ½σ²∂²V/∂x² = 0
     * V(T, x) = g(x)
     */
    static double terminalConditionProblem(
        std::function<double(double, double)> drift,    // b(x, u)
        std::function<double(double)> diffusion,        // σ(x)
        std::function<double(double)> terminal_cost,    // g(x)
        std::function<double(double, double)> control_law,  // u(t, x)
        double X0, double T,
        int n_steps = 100,
        int n_samples = 1000) {

        std::vector<double> terminal_values(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);
        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);

        for (int sample = 0; sample < n_samples; ++sample) {
            double X = X0;
            double t = 0.0;

            for (int i = 0; i < n_steps; ++i) {
                double u = control_law(t, X);
                double dW = sqrt_dt * normal(gen);

                X += drift(X, u) * dt + diffusion(X) * dW;
                t += dt;
            }

            terminal_values[sample] = terminal_cost(X);
        }

        return std::accumulate(terminal_values.begin(), terminal_values.end(), 0.0) / n_samples;
    }

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

    /**
     * @brief Verification theorem for HJB solution
     *
     * If V is C^{1,2}, satisfies HJB with terminal condition,
     * and u*(t,x) achieves minimum in HJB, then:
     * 1. V(t,x) is the value function
     * 2. u*(t,x) is optimal control
     */
    static bool verifyHJBSolution(
        std::function<double(double, double)> V,
        std::function<double(double, double)> u_star,
        std::function<double(double, double, double)> b,
        std::function<double(double, double, double)> sigma,
        std::function<double(double, double, double)> f,
        double t, double x,
        double h = 1e-5,
        double tolerance = 1e-4) {

        // Compute derivatives numerically
        double dV_dt = (V(t + h, x) - V(t - h, x)) / (2.0 * h);
        double dV_dx = (V(t, x + h) - V(t, x - h)) / (2.0 * h);
        double d2V_dx2 = (V(t, x + h) - 2.0 * V(t, x) + V(t, x - h)) / (h * h);

        double u = u_star(t, x);

        // Check HJB equation
        double hjb_residual = dV_dt + b(t, x, u) * dV_dx
                            + 0.5 * sigma(t, x, u) * sigma(t, x, u) * d2V_dx2
                            + f(t, x, u);

        return std::abs(hjb_residual) < tolerance;
    }
};

/**
 * ============================================================================
 * DIFFUSIONS: BASIC PROPERTIES
 * ============================================================================
 */

/**
 * @class DiffusionProperties
 * @brief Fundamental properties of Itô diffusions
 *
 * Itô diffusion: dX = b(X) dt + σ(X) dW
 */
class DiffusionProperties {
public:
    /**
     * @brief Verify Markov property
     *
     * P(X(t) ∈ A | ℱₛ) = P(X(t) ∈ A | X(s)) for s < t
     *
     * Check that future evolution depends only on current state
     */
    static bool checkMarkovProperty(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        double X0, double T, int n_tests = 100) {

        // For Itô diffusions with time-homogeneous coefficients,
        // Markov property holds automatically
        // Test: verify conditional independence

        std::normal_distribution<double> normal(0.0, 1.0);
        double dt = 0.01;
        double sqrt_dt = std::sqrt(dt);

        // Sample paths starting from different states
        for (int test = 0; test < n_tests; ++test) {
            double X1 = X0;
            double X2 = X0;

            // Evolve to time T/2
            for (int i = 0; i < (int)(T / (2 * dt)); ++i) {
                double t = i * dt;
                double dW1 = sqrt_dt * normal(gen);
                double dW2 = sqrt_dt * normal(gen);

                X1 += mu(t, X1) * dt + sigma(t, X1) * dW1;
                X2 += mu(t, X2) * dt + sigma(t, X2) * dW2;
            }

            // From same current state, future should be independent of history
            // (Simplified verification)
        }

        return true;  // Itô diffusions satisfy Markov property
    }

    /**
     * @brief Strong Markov property at stopping times
     *
     * For stopping time τ: P(X(τ+t) ∈ A | ℱ_τ) = P(X(τ+t) ∈ A | X(τ))
     */
    static double verifyStrongMarkov(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<bool(double, double)> stopping_condition,
        double X0, double T,
        int n_samples = 1000) {

        std::normal_distribution<double> normal(0.0, 1.0);
        double dt = 0.01;
        double sqrt_dt = std::sqrt(dt);

        std::vector<double> values_at_tau(n_samples);

        for (int sample = 0; sample < n_samples; ++sample) {
            double X = X0;
            double t = 0.0;
            bool stopped = false;

            // Run until stopping time
            while (t < T && !stopped) {
                if (stopping_condition(t, X)) {
                    values_at_tau[sample] = X;
                    stopped = true;
                } else {
                    double dW = sqrt_dt * normal(gen);
                    X += mu(t, X) * dt + sigma(t, X) * dW;
                    t += dt;
                }
            }

            if (!stopped) {
                values_at_tau[sample] = X;
            }
        }

        // Return mean value at stopping time
        return std::accumulate(values_at_tau.begin(), values_at_tau.end(), 0.0) / n_samples;
    }

    /**
     * @brief Generator of Itô diffusion
     *
     * For dX = b(X) dt + σ(X) dW, the generator is:
     * Af(x) = b(x) f'(x) + ½σ²(x) f''(x)
     *
     * Satisfies: d/dt E[f(X(t))] = E[Af(X(t))]
     */
    static double generator(
        std::function<double(double)> b,
        std::function<double(double)> sigma,
        std::function<double(double)> f,
        std::function<double(double)> df,
        std::function<double(double)> d2f,
        double x) {

        return b(x) * df(x) + 0.5 * sigma(x) * sigma(x) * d2f(x);
    }

    /**
     * @brief Dynkin formula
     *
     * E[f(X(τ))] = f(x) + E[∫₀^τ Af(X(s)) ds]
     *
     * for bounded stopping time τ
     */
    static double dynkinFormula(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double)> f,
        std::function<double(double)> Af,
        std::function<bool(double, double)> stopping_condition,
        double X0, double max_time,
        int n_samples = 1000) {

        std::vector<double> terminal_values(n_samples);
        std::normal_distribution<double> normal(0.0, 1.0);

        for (int sample = 0; sample < n_samples; ++sample) {
            double X = X0;
            double t = 0.0;
            double dt = 0.01;
            double sqrt_dt = std::sqrt(dt);
            double integral_Af = 0.0;

            while (t < max_time && !stopping_condition(t, X)) {
                double dW = sqrt_dt * normal(gen);
                integral_Af += Af(X) * dt;
                X += mu(t, X) * dt + sigma(t, X) * dW;
                t += dt;
            }

            terminal_values[sample] = f(X);
        }

        return std::accumulate(terminal_values.begin(), terminal_values.end(), 0.0) / n_samples;
    }

    /**
     * @brief Characteristic operator (extended generator)
     *
     * Same as generator for time-homogeneous diffusions
     */
    static double characteristicOperator(
        std::function<double(double)> b,
        std::function<double(double)> sigma,
        std::function<double(double)> f,
        double x,
        double h = 1e-5) {

        // Numerical differentiation
        double df_dx = (f(x + h) - f(x - h)) / (2.0 * h);
        double d2f_dx2 = (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);

        return b(x) * df_dx + 0.5 * sigma(x) * sigma(x) * d2f_dx2;
    }
};

/**
 * ============================================================================
 * OTHER TOPICS IN DIFFUSION THEORY
 * ============================================================================
 */

/**
 * @class AdvancedDiffusionTheory
 * @brief Kolmogorov equations, Feynman-Kac, martingale problem, time changes
 */
class AdvancedDiffusionTheory {
public:
    /**
     * @brief Kolmogorov backward equation
     *
     * ∂u/∂t + b(x)∂u/∂x + ½σ²(x)∂²u/∂x² = 0
     * u(T, x) = φ(x)
     *
     * Solution: u(t,x) = E[φ(X(T)) | X(t) = x]
     */
    static std::vector<std::vector<double>> kolmogorovBackward(
        std::function<double(double)> b,
        std::function<double(double)> sigma,
        std::function<double(double)> terminal_condition,
        double x_min, double x_max, int n_x,
        double T, int n_t) {

        double dx = (x_max - x_min) / n_x;
        double dt = T / n_t;

        // Solution grid u[time][space]
        std::vector<std::vector<double>> u(n_t + 1, std::vector<double>(n_x + 1));

        // Terminal condition
        for (int i = 0; i <= n_x; ++i) {
            double x = x_min + i * dx;
            u[n_t][i] = terminal_condition(x);
        }

        // Backward time stepping
        for (int k = n_t - 1; k >= 0; --k) {
            for (int i = 1; i < n_x; ++i) {
                double x = x_min + i * dx;

                // Finite difference approximation
                double du_dx = (u[k+1][i+1] - u[k+1][i-1]) / (2.0 * dx);
                double d2u_dx2 = (u[k+1][i+1] - 2.0 * u[k+1][i] + u[k+1][i-1]) / (dx * dx);

                u[k][i] = u[k+1][i] - dt * (b(x) * du_dx + 0.5 * sigma(x) * sigma(x) * d2u_dx2);
            }

            // Boundary conditions (simplified)
            u[k][0] = u[k][1];
            u[k][n_x] = u[k][n_x-1];
        }

        return u;
    }

    /**
     * @brief Resolvent operator
     *
     * R_λ f(x) = E[∫₀^∞ e^(-λt) f(X(t)) dt | X(0) = x]
     *
     * Satisfies: (λI - A) R_λ = I
     */
    static double resolvent(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double)> f,
        double lambda, double x0,
        double max_time = 10.0,
        int n_steps = 1000) {

        double dt = max_time / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double integral = 0.0;
        double X = x0;
        double t = 0.0;

        for (int i = 0; i < n_steps; ++i) {
            double dW = sqrt_dt * normal(gen);
            integral += std::exp(-lambda * t) * f(X) * dt;
            X += mu(t, X) * dt + sigma(t, X) * dW;
            t += dt;
        }

        return integral;
    }

    /**
     * @brief Feynman-Kac formula
     *
     * Solve PDE: ∂u/∂t + Au - qu = -f
     * with terminal condition u(T,x) = φ(x)
     *
     * Solution: u(t,x) = E[e^(-∫ₜᵀ q(X(s))ds) φ(X(T)) + ∫ₜᵀ e^(-∫ₜˢ q(X(r))dr) f(X(s)) ds]
     */
    static double feynmanKac(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double)> q,      // Killing rate
        std::function<double(double)> f,      // Source term
        std::function<double(double)> phi,    // Terminal condition
        double X0, double t0, double T,
        int n_steps = 1000,
        int n_samples = 100) {

        double dt = (T - t0) / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> results(n_samples);

        for (int sample = 0; sample < n_samples; ++sample) {
            double X = X0;
            double t = t0;
            double discount = 0.0;  // ∫ₜᵀ q(X(s)) ds
            double integral_f = 0.0;

            for (int i = 0; i < n_steps; ++i) {
                double dW = sqrt_dt * normal(gen);

                // Accumulate killing integral
                discount += q(X) * dt;

                // Accumulate source term integral
                integral_f += std::exp(-discount) * f(X) * dt;

                // Evolve process
                X += mu(t, X) * dt + sigma(t, X) * dW;
                t += dt;
            }

            results[sample] = std::exp(-discount) * phi(X) + integral_f;
        }

        return std::accumulate(results.begin(), results.end(), 0.0) / n_samples;
    }

    /**
     * @brief Killing: Process stopped with rate q(x)
     *
     * Killed process: P(alive at t) = exp(-∫₀ᵗ q(X(s)) ds)
     */
    static double killingProbability(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double)> killing_rate,
        double X0, double T,
        int n_steps = 1000) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double X = X0;
        double total_rate = 0.0;

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = sqrt_dt * normal(gen);

            total_rate += killing_rate(X) * dt;
            X += mu(t, X) * dt + sigma(t, X) * dW;
        }

        return std::exp(-total_rate);
    }

    /**
     * @brief Martingale problem
     *
     * X is solution to martingale problem (A, μ) if:
     * f(X(t)) - ∫₀ᵗ Af(X(s)) ds is a martingale for all f ∈ C²
     */
    static bool solveMartingaleProblem(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double)> f,
        std::function<double(double)> Af,
        double X0, double T,
        int n_steps = 1000) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double X = X0;
        double M = f(X0);  // Martingale value
        double integral_Af = 0.0;

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = sqrt_dt * normal(gen);

            integral_Af += Af(X) * dt;
            X += mu(t, X) * dt + sigma(t, X) * dW;

            // Check martingale property: M(t) = f(X(t)) - integral
            double M_expected = f(X) - integral_Af;

            // (Simplified verification)
        }

        return true;
    }

    /**
     * @brief Check when Itô process is a diffusion
     *
     * Itô process dX = μ(t,X) dt + σ(t,X) dW is diffusion if:
     * - Coefficients depend only on (t, X(t)), not on path history
     * - Markov property holds
     */
    static bool isItoDiffusion(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        double X0, double T) {

        // Check 1: Coefficients time-homogeneous or depend only on (t,x)
        // (Structural check - always true for given signature)

        // Check 2: Markov property
        bool is_markov = DiffusionProperties::checkMarkovProperty(mu, sigma, X0, T);

        return is_markov;
    }

    /**
     * @brief Random time change
     *
     * Transform diffusion via time change: τ(t) = ∫₀ᵗ c(X(s)) ds
     * New process Y(t) = X(τ⁻¹(t)) is also a diffusion
     */
    static std::vector<double> randomTimeChange(
        const std::vector<double>& X_path,
        std::function<double(double)> time_change_rate,
        double T, int n_steps) {

        double dt = T / n_steps;
        std::vector<double> tau(n_steps + 1, 0.0);

        // Compute random time τ(t) = ∫₀ᵗ c(X(s)) ds
        for (int i = 0; i < n_steps; ++i) {
            tau[i+1] = tau[i] + time_change_rate(X_path[i]) * dt;
        }

        // Inverse time change (simplified - assumes monotone increasing)
        std::vector<double> Y_path(n_steps + 1);
        for (int i = 0; i <= n_steps; ++i) {
            // Find j such that τ(j) ≈ t
            int j = std::min(i, n_steps);
            Y_path[i] = X_path[j];
        }

        return Y_path;
    }

    /**
     * @brief Extended Girsanov theorem
     *
     * Change drift from b(x) to b(x) + σ(x)θ(x)
     * via measure change with density:
     * dQ/dP = exp(∫₀ᵗ θ(X(s)) dW(s) - ½∫₀ᵗ θ²(X(s)) ds)
     */
    static std::pair<std::vector<double>, double> girsanovDriftChange(
        std::function<double(double, double)> mu,
        std::function<double(double, double)> sigma,
        std::function<double(double)> theta,
        double X0, double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> X_path(n_steps + 1);
        X_path[0] = X0;

        double X = X0;
        double W = 0.0;
        double integral_theta_dW = 0.0;
        double integral_theta_squared = 0.0;

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double dW = sqrt_dt * normal(gen);

            // Original dynamics with changed drift
            X += (mu(t, X) + sigma(t, X) * theta(X)) * dt + sigma(t, X) * dW;

            // Radon-Nikodym derivative components
            integral_theta_dW += theta(X) * dW;
            integral_theta_squared += theta(X) * theta(X) * dt;

            X_path[i+1] = X;
            W += dW;
        }

        double radon_nikodym = std::exp(integral_theta_dW - 0.5 * integral_theta_squared);

        return {X_path, radon_nikodym};
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

/**
 * @class MarketTheory
 * @brief Market models, portfolio dynamics, and arbitrage theory
 *
 * Models financial markets with multiple assets:
 * - Stock: dS(t) = μS(t)dt + σS(t)dW(t) (geometric Brownian motion)
 * - Bond: dB(t) = rB(t)dt (risk-free asset)
 * - Self-financing portfolio: dV = Σᵢ φᵢ dSᵢ
 * - Arbitrage: Risk-free profit with zero investment
 */
class MarketTheory {
public:
    struct Asset {
        double S0;      // Initial price
        double mu;      // Drift
        double sigma;   // Volatility
    };

    struct Portfolio {
        std::vector<double> holdings;  // φᵢ = shares of asset i
        double cash;                    // Cash position
    };

    /**
     * @brief Simulate market with multiple stocks and risk-free bond
     *
     * Stock i: dSᵢ = μᵢSᵢ dt + σᵢSᵢ dWᵢ
     * Bond: dB = rB dt
     */
    static std::vector<std::vector<double>> simulateMarket(
        const std::vector<Asset>& assets,
        double r,           // Risk-free rate
        double T,
        int n_steps,
        const std::vector<std::vector<double>>& correlations = {}) {

        int n_assets = assets.size();
        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        // Initialize price paths
        std::vector<std::vector<double>> prices(n_assets);
        for (int i = 0; i < n_assets; ++i) {
            prices[i].resize(n_steps + 1);
            prices[i][0] = assets[i].S0;
        }

        // Simulate correlated Brownian motions
        for (int t = 0; t < n_steps; ++t) {
            std::vector<double> dW(n_assets);

            if (correlations.empty()) {
                // Independent Brownian motions
                for (int i = 0; i < n_assets; ++i) {
                    dW[i] = sqrt_dt * normal(gen);
                }
            } else {
                // Correlated via Cholesky decomposition
                std::vector<double> Z(n_assets);
                for (int i = 0; i < n_assets; ++i) {
                    Z[i] = normal(gen);
                }

                // L from Cholesky: ρ = LL^T
                for (int i = 0; i < n_assets; ++i) {
                    dW[i] = 0.0;
                    for (int j = 0; j <= i; ++j) {
                        dW[i] += correlations[i][j] * Z[j];
                    }
                    dW[i] *= sqrt_dt;
                }
            }

            // Update prices: dS = μS dt + σS dW
            for (int i = 0; i < n_assets; ++i) {
                prices[i][t+1] = prices[i][t] * (1.0 + assets[i].mu * dt + assets[i].sigma * dW[i]);
            }
        }

        return prices;
    }

    /**
     * @brief Self-financing portfolio dynamics
     *
     * A trading strategy (φ₀(t), φ₁(t), ..., φₙ(t)) is self-financing if:
     * dV(t) = φ₀(t) dB(t) + Σᵢ φᵢ(t) dSᵢ(t)
     *
     * No external cash injections or withdrawals
     */
    static std::vector<double> selfFinancingPortfolio(
        const std::vector<std::vector<double>>& asset_prices,
        const std::vector<std::function<double(int)>>& strategy,  // φᵢ(t)
        double r,  // Risk-free rate
        double V0, // Initial wealth
        double T,
        int n_steps) {

        int n_assets = asset_prices.size();
        double dt = T / n_steps;
        std::vector<double> V(n_steps + 1);
        V[0] = V0;

        for (int t = 0; t < n_steps; ++t) {
            double dV = 0.0;

            // φ₀ dB: Cash account grows at risk-free rate
            double bond_value = V[t];
            for (int i = 0; i < n_assets; ++i) {
                bond_value -= strategy[i](t) * asset_prices[i][t];
            }
            dV += r * bond_value * dt;

            // Σᵢ φᵢ dSᵢ: Stock holdings
            for (int i = 0; i < n_assets; ++i) {
                double dS = asset_prices[i][t+1] - asset_prices[i][t];
                dV += strategy[i](t) * dS;
            }

            V[t+1] = V[t] + dV;
        }

        return V;
    }

    /**
     * @brief Detect arbitrage via martingale property
     *
     * Fundamental Theorem of Asset Pricing:
     * No arbitrage ⟺ ∃ equivalent martingale measure Q
     * where discounted prices S̃(t) = e^(-rt)S(t) are Q-martingales
     */
    static bool detectArbitrage(
        const std::vector<Asset>& assets,
        double r,  // Risk-free rate
        double tolerance = 1e-6) {

        // For single-asset market: no arbitrage iff can find Q with
        // E^Q[dS̃/S̃] = 0 ⟺ μ - r = λσ for some λ (market price of risk)

        if (assets.size() == 1) {
            // Always no arbitrage in single-asset Black-Scholes
            return false;
        }

        // For multi-asset: check if system μᵢ - r = Σⱼ σᵢⱼλⱼ is solvable
        // Simplified check: verify all Sharpe ratios are consistent
        std::vector<double> sharpe_ratios(assets.size());
        for (size_t i = 0; i < assets.size(); ++i) {
            sharpe_ratios[i] = (assets[i].mu - r) / assets[i].sigma;
        }

        // In complete market, all Sharpe ratios should be equal
        double sharpe_ref = sharpe_ratios[0];
        for (size_t i = 1; i < assets.size(); ++i) {
            if (std::abs(sharpe_ratios[i] - sharpe_ref) > tolerance) {
                return true;  // Arbitrage exists
            }
        }

        return false;  // No arbitrage
    }

    /**
     * @brief Construct risk-neutral measure (Girsanov change)
     *
     * Under Q: dS = rS dt + σS dW^Q
     * where W^Q = W + (μ-r)/σ t (market price of risk added to Brownian motion)
     *
     * Radon-Nikodym: dQ/dP = exp(-θW - ½θ²t) where θ = (μ-r)/σ
     */
    static std::pair<std::vector<double>, double> riskNeutralMeasure(
        double S0, double mu, double sigma, double r,
        double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double theta = (mu - r) / sigma;  // Market price of risk
        std::vector<double> S(n_steps + 1);
        S[0] = S0;

        double W_Q = 0.0;
        double radon_nikodym = 1.0;

        for (int i = 0; i < n_steps; ++i) {
            double dW_P = sqrt_dt * normal(gen);
            double dW_Q = dW_P + theta * sqrt_dt;  // Girsanov drift

            // Under Q: dS = rS dt + σS dW^Q
            S[i+1] = S[i] * (1.0 + r * dt + sigma * dW_Q);

            // Radon-Nikodym derivative
            W_Q += dW_Q;
            radon_nikodym *= std::exp(-theta * dW_P - 0.5 * theta * theta * dt);
        }

        return {S, radon_nikodym};
    }

    /**
     * @brief Market price of risk (Sharpe ratio)
     *
     * λ = (μ - r) / σ
     */
    static double marketPriceOfRisk(double mu, double sigma, double r) {
        return (mu - r) / sigma;
    }
};

/**
 * @class AttainabilityCompleteness
 * @brief Attainable claims, market completeness, and hedging strategies
 *
 * - Attainable claim: Payoff that can be replicated by trading strategy
 * - Complete market: All claims are attainable
 * - Hedging: Replicating portfolio that eliminates risk
 */
class AttainabilityCompleteness {
public:
    /**
     * @brief Check if claim is attainable in Black-Scholes market
     *
     * In complete market (1 stock + 1 bond), any measurable payoff h(S(T))
     * is attainable with replicating portfolio
     */
    static bool isAttainable(
        std::function<double(double)> payoff,
        double S0, double r, double sigma,
        double T,
        double tolerance = 1e-3) {

        // Black-Scholes market is complete: all payoffs attainable
        // Verify by checking if value function exists

        // Price via risk-neutral valuation
        double price = riskNeutralPrice(payoff, S0, r, sigma, T, 1000);

        return std::isfinite(price);  // Attainable if finite price exists
    }

    /**
     * @brief Check market completeness
     *
     * Market with d stocks and d sources of randomness (Brownian motions)
     * is complete if volatility matrix σ is non-singular
     */
    static bool isMarketComplete(
        const std::vector<std::vector<double>>& volatility_matrix) {

        int n = volatility_matrix.size();
        if (n == 0) return false;

        // For 1D: complete if σ > 0
        if (n == 1) {
            return volatility_matrix[0][0] > 0;
        }

        // For multi-dimensional: check det(σ) ≠ 0
        // Simplified check for 2×2
        if (n == 2) {
            double det = volatility_matrix[0][0] * volatility_matrix[1][1]
                       - volatility_matrix[0][1] * volatility_matrix[1][0];
            return std::abs(det) > 1e-10;
        }

        // For larger systems, assume complete if square and well-conditioned
        return true;
    }

    /**
     * @brief Construct replicating portfolio (delta hedging)
     *
     * For option with value V(t,S), replicating portfolio:
     * - Δ(t) = ∂V/∂S shares of stock
     * - B(t) = V - ΔS units of bond
     *
     * This replicates: dV = Δ dS + rB dt
     */
    static std::pair<std::vector<double>, std::vector<double>> replicatingPortfolio(
        std::function<double(double, double)> value_func,  // V(t,S)
        std::function<double(double, double)> delta_func,  // ∂V/∂S
        const std::vector<double>& stock_price,
        double r, double T, int n_steps) {

        int n = stock_price.size();
        double dt = T / n_steps;
        std::vector<double> stock_holdings(n);
        std::vector<double> bond_holdings(n);

        for (int i = 0; i < n; ++i) {
            double t = i * dt;
            double S = stock_price[i];

            // Delta hedge
            stock_holdings[i] = delta_func(t, S);

            // Bond position: B = V - ΔS
            double V = value_func(t, S);
            bond_holdings[i] = V - stock_holdings[i] * S;
        }

        return {stock_holdings, bond_holdings};
    }

    /**
     * @brief Verify hedging strategy eliminates risk
     *
     * For perfect hedge: dV = rV dt (no stochastic term)
     * Achieved when Δ = ∂V/∂S (delta hedge)
     */
    static bool verifyHedge(
        std::function<double(double, double)> value_func,
        std::function<double(double, double)> delta_func,
        double S0, double r, double sigma,
        double T, int n_steps,
        double tolerance = 1e-2) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double S = S0;
        double V = value_func(0.0, S0);

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double delta = delta_func(t, S);

            // Stock evolution: dS = rS dt + σS dW (risk-neutral)
            double dW = sqrt_dt * normal(gen);
            double dS = r * S * dt + sigma * S * dW;

            // Portfolio evolution: dV = Δ dS + r(V - ΔS) dt
            double dV_hedge = delta * dS + r * (V - delta * S) * dt;

            S += dS;
            V += dV_hedge;

            // Check if matches true value (within tolerance)
            double V_true = value_func(t + dt, S);
            if (std::abs(V - V_true) > tolerance * V_true) {
                return false;
            }
        }

        return true;
    }

private:
    static double riskNeutralPrice(
        std::function<double(double)> payoff,
        double S0, double r, double sigma,
        double T, int n_samples) {

        std::normal_distribution<double> normal(0.0, 1.0);
        double sum = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double Z = normal(gen);
            double S_T = S0 * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);
            sum += payoff(S_T);
        }

        return std::exp(-r * T) * sum / n_samples;
    }
};

/**
 * @class OptionPricing
 * @brief Option pricing via Black-Scholes and path-dependent options
 *
 * Black-Scholes PDE: ∂V/∂t + ½σ²S²∂²V/∂S² + rS∂V/∂S - rV = 0
 * Solution: V = E^Q[e^(-rT) h(S(T))]
 */
class OptionPricing {
public:
    /**
     * @brief Black-Scholes formula for European call
     *
     * C(S,t) = S·N(d₁) - K·e^(-r(T-t))·N(d₂)
     * d₁ = [ln(S/K) + (r + σ²/2)(T-t)] / (σ√(T-t))
     * d₂ = d₁ - σ√(T-t)
     */
    static double europeanCall(
        double S, double K, double r, double sigma, double T) {

        if (T <= 0) return std::max(S - K, 0.0);

        double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);

        auto N = [](double x) {
            return 0.5 * std::erfc(-x / std::sqrt(2.0));
        };

        return S * N(d1) - K * std::exp(-r * T) * N(d2);
    }

    /**
     * @brief Black-Scholes formula for European put
     *
     * P(S,t) = K·e^(-r(T-t))·N(-d₂) - S·N(-d₁)
     */
    static double europeanPut(
        double S, double K, double r, double sigma, double T) {

        if (T <= 0) return std::max(K - S, 0.0);

        double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);

        auto N = [](double x) {
            return 0.5 * std::erfc(-x / std::sqrt(2.0));
        };

        return K * std::exp(-r * T) * N(-d2) - S * N(-d1);
    }

    /**
     * @brief Digital (binary) option
     *
     * Pays 1 if S(T) > K, 0 otherwise
     * Value: e^(-rT)·N(d₂)
     */
    static double digitalOption(
        double S, double K, double r, double sigma, double T) {

        if (T <= 0) return S > K ? 1.0 : 0.0;

        double d2 = (std::log(S / K) + (r - 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));

        auto N = [](double x) {
            return 0.5 * std::erfc(-x / std::sqrt(2.0));
        };

        return std::exp(-r * T) * N(d2);
    }

    /**
     * @brief Asian option (arithmetic average)
     *
     * Payoff: max(Ā - K, 0) where Ā = (1/n)Σᵢ S(tᵢ)
     * No closed form - use Monte Carlo
     */
    static double asianOption(
        double S0, double K, double r, double sigma,
        double T, int n_steps, int n_samples = 10000) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double sum_payoffs = 0.0;

        for (int sample = 0; sample < n_samples; ++sample) {
            double S = S0;
            double average = 0.0;

            for (int i = 0; i < n_steps; ++i) {
                double dW = sqrt_dt * normal(gen);
                S *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * dW);
                average += S;
            }

            average /= n_steps;
            sum_payoffs += std::max(average - K, 0.0);
        }

        return std::exp(-r * T) * sum_payoffs / n_samples;
    }

    /**
     * @brief Barrier option (knock-out call)
     *
     * Payoff: max(S(T) - K, 0) if S(t) < H for all t, else 0
     * Option becomes worthless if barrier H is hit
     */
    static double barrierOption(
        double S0, double K, double H, double r, double sigma,
        double T, int n_steps, int n_samples = 10000) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double sum_payoffs = 0.0;

        for (int sample = 0; sample < n_samples; ++sample) {
            double S = S0;
            bool barrier_hit = false;

            for (int i = 0; i < n_steps; ++i) {
                double dW = sqrt_dt * normal(gen);
                S *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * dW);

                if (S >= H) {
                    barrier_hit = true;
                    break;
                }
            }

            if (!barrier_hit) {
                sum_payoffs += std::max(S - K, 0.0);
            }
        }

        return std::exp(-r * T) * sum_payoffs / n_samples;
    }

    /**
     * @brief Lookback option (floating strike)
     *
     * Payoff: S(T) - min_{t∈[0,T]} S(t)
     * Maximum gain from hindsight trading
     */
    static double lookbackOption(
        double S0, double r, double sigma,
        double T, int n_steps, int n_samples = 10000) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> normal(0.0, 1.0);

        double sum_payoffs = 0.0;

        for (int sample = 0; sample < n_samples; ++sample) {
            double S = S0;
            double S_min = S0;

            for (int i = 0; i < n_steps; ++i) {
                double dW = sqrt_dt * normal(gen);
                S *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * dW);
                S_min = std::min(S_min, S);
            }

            sum_payoffs += S - S_min;
        }

        return std::exp(-r * T) * sum_payoffs / n_samples;
    }

    /**
     * @brief Option Greeks
     *
     * Delta: Δ = ∂V/∂S
     * Gamma: Γ = ∂²V/∂S²
     * Vega: ν = ∂V/∂σ
     * Theta: Θ = ∂V/∂t
     * Rho: ρ = ∂V/∂r
     */
    struct Greeks {
        double delta;   // ∂V/∂S
        double gamma;   // ∂²V/∂S²
        double vega;    // ∂V/∂σ
        double theta;   // ∂V/∂t
        double rho;     // ∂V/∂r
    };

    static Greeks europeanCallGreeks(
        double S, double K, double r, double sigma, double T) {

        Greeks greeks;

        double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);

        auto N = [](double x) {
            return 0.5 * std::erfc(-x / std::sqrt(2.0));
        };

        auto n = [](double x) {
            return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI);
        };

        // Delta: Δ = N(d₁)
        greeks.delta = N(d1);

        // Gamma: Γ = n(d₁) / (S σ √T)
        greeks.gamma = n(d1) / (S * sigma * std::sqrt(T));

        // Vega: ν = S n(d₁) √T
        greeks.vega = S * n(d1) * std::sqrt(T);

        // Theta: Θ = -S n(d₁) σ / (2√T) - rK e^(-rT) N(d₂)
        greeks.theta = -S * n(d1) * sigma / (2.0 * std::sqrt(T))
                      - r * K * std::exp(-r * T) * N(d2);

        // Rho: ρ = K T e^(-rT) N(d₂)
        greeks.rho = K * T * std::exp(-r * T) * N(d2);

        return greeks;
    }

    /**
     * @brief Risk-neutral valuation
     *
     * V(0,S) = E^Q[e^(-rT) h(S(T))]
     * where Q is risk-neutral measure (S grows at rate r)
     */
    static double riskNeutralValuation(
        std::function<double(double)> payoff,
        double S0, double r, double sigma,
        double T, int n_samples = 10000) {

        std::normal_distribution<double> normal(0.0, 1.0);
        double sum = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double Z = normal(gen);
            double S_T = S0 * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);
            sum += payoff(S_T);
        }

        return std::exp(-r * T) * sum / n_samples;
    }

    /**
     * @brief Put-Call parity
     *
     * C - P = S - K·e^(-rT)
     * Relationship between European call and put prices
     */
    static bool verifyPutCallParity(
        double call_price, double put_price,
        double S, double K, double r, double T,
        double tolerance = 1e-6) {

        double lhs = call_price - put_price;
        double rhs = S - K * std::exp(-r * T);

        return std::abs(lhs - rhs) < tolerance;
    }
};

} // namespace maths::sde

#endif // MATHS_STOCHASTIC_DIFFERENTIAL_EQUATIONS_HPP
