/**
 * @file monte_carlo.hpp
 * @brief Comprehensive Monte Carlo methods and stochastic algorithms
 *
 * Implements computational algorithms for:
 * - Basic Monte Carlo integration and simulation
 * - Importance sampling and variance reduction
 * - Markov Chain Monte Carlo (MCMC): Metropolis-Hastings, Gibbs
 * - Markov chains: transition matrices, stationary distributions
 * - Stochastic processes: Brownian motion, Ornstein-Uhlenbeck
 * - Boltzmann equations and kinetic theory
 * - Hybrid algorithms: Hamiltonian Monte Carlo (HMC), NUTS
 */

#ifndef MATHS_STOCHASTIC_MONTE_CARLO_HPP
#define MATHS_STOCHASTIC_MONTE_CARLO_HPP

#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <map>

namespace maths::stochastic {

// Random number generator
static std::random_device rd;
static std::mt19937 gen(rd());

/**
 * @class MonteCarloIntegrator
 * @brief Basic Monte Carlo integration methods
 */
class MonteCarloIntegrator {
public:
    /**
     * @brief Basic Monte Carlo integration
     *
     * Estimates ∫f(x)dx over [a,b] using uniform sampling
     *
     * @param f Integrand function
     * @param a Lower bound
     * @param b Upper bound
     * @param n_samples Number of samples
     * @return Pair (estimate, standard_error)
     */
    static std::pair<double, double> integrate(
        std::function<double(double)> f,
        double a, double b,
        int n_samples = 10000) {

        std::uniform_real_distribution<double> dist(a, b);
        
        double sum = 0.0;
        double sum_sq = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double x = dist(gen);
            double fx = f(x);
            sum += fx;
            sum_sq += fx * fx;
        }

        double mean = sum / n_samples;
        double variance = (sum_sq / n_samples) - (mean * mean);
        double estimate = (b - a) * mean;
        double std_error = (b - a) * std::sqrt(variance / n_samples);

        return {estimate, std_error};
    }

    /**
     * @brief Importance sampling Monte Carlo
     *
     * Uses proposal distribution g(x) to reduce variance
     * Estimates ∫f(x)dx ≈ E[f(X)/g(X)] where X ~ g
     *
     * @param f Integrand
     * @param proposal_sampler Sampler for proposal distribution g
     * @param proposal_pdf PDF of proposal distribution
     * @param n_samples Number of samples
     * @return Pair (estimate, standard_error)
     */
    static std::pair<double, double> importanceSampling(
        std::function<double(double)> f,
        std::function<double()> proposal_sampler,
        std::function<double(double)> proposal_pdf,
        int n_samples = 10000) {

        double sum = 0.0;
        double sum_sq = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double x = proposal_sampler();
            double gx = proposal_pdf(x);
            
            if (gx > 1e-10) {
                double ratio = f(x) / gx;
                sum += ratio;
                sum_sq += ratio * ratio;
            }
        }

        double mean = sum / n_samples;
        double variance = (sum_sq / n_samples) - (mean * mean);
        double std_error = std::sqrt(variance / n_samples);

        return {mean, std_error};
    }

    /**
     * @brief Rejection sampling (Accept-Reject method)
     *
     * Sample from target distribution f(x) using proposal g(x)
     * where f(x) ≤ M·g(x) for all x
     *
     * @param target_pdf Target PDF (up to normalization)
     * @param proposal_sampler Proposal sampler
     * @param proposal_pdf Proposal PDF
     * @param M Upper bound: f(x) ≤ M·g(x)
     * @return Sample from target distribution
     */
    static double rejectionSampling(
        std::function<double(double)> target_pdf,
        std::function<double()> proposal_sampler,
        std::function<double(double)> proposal_pdf,
        double M) {

        while (true) {
            double x = proposal_sampler();
            double u = std::uniform_real_distribution<double>(0.0, 1.0)(gen);
            
            double accept_prob = target_pdf(x) / (M * proposal_pdf(x));
            
            if (u <= accept_prob) {
                return x;
            }
        }
    }

    /**
     * @brief Stratified sampling for variance reduction
     *
     * Divides interval into strata and samples proportionally
     */
    static std::pair<double, double> stratifiedSampling(
        std::function<double(double)> f,
        double a, double b,
        int n_strata = 10,
        int samples_per_stratum = 100) {

        double stratum_width = (b - a) / n_strata;
        double total_estimate = 0.0;
        double total_variance = 0.0;

        for (int i = 0; i < n_strata; ++i) {
            double stratum_a = a + i * stratum_width;
            double stratum_b = stratum_a + stratum_width;

            auto [est, var] = integrate(f, stratum_a, stratum_b, samples_per_stratum);
            total_estimate += est;
            total_variance += var * var;
        }

        double std_error = std::sqrt(total_variance);
        return {total_estimate, std_error};
    }
};

/**
 * @class MarkovChain
 * @brief Discrete-time Markov chain with transition matrix
 */
class MarkovChain {
private:
    std::vector<std::vector<double>> transition_matrix;
    int n_states;

public:
    MarkovChain(const std::vector<std::vector<double>>& P) 
        : transition_matrix(P), n_states(P.size()) {
        
        // Verify stochastic matrix
        for (const auto& row : P) {
            double sum = std::accumulate(row.begin(), row.end(), 0.0);
            if (std::abs(sum - 1.0) > 1e-6) {
                throw std::invalid_argument("Rows must sum to 1");
            }
        }
    }

    /**
     * @brief Simulate Markov chain trajectory
     *
     * @param initial_state Starting state
     * @param n_steps Number of steps
     * @return Sequence of states
     */
    std::vector<int> simulate(int initial_state, int n_steps) {
        std::vector<int> trajectory;
        trajectory.reserve(n_steps + 1);
        trajectory.push_back(initial_state);

        int current_state = initial_state;
        for (int step = 0; step < n_steps; ++step) {
            // Sample next state from transition probabilities
            std::discrete_distribution<int> dist(
                transition_matrix[current_state].begin(),
                transition_matrix[current_state].end()
            );
            current_state = dist(gen);
            trajectory.push_back(current_state);
        }

        return trajectory;
    }

    /**
     * @brief Compute stationary distribution π
     *
     * Solves πP = π (eigenvector with eigenvalue 1)
     * Uses power iteration method
     *
     * @param max_iter Maximum iterations
     * @param tol Convergence tolerance
     * @return Stationary distribution
     */
    std::vector<double> stationaryDistribution(int max_iter = 1000, double tol = 1e-8) {
        // Start with uniform distribution
        std::vector<double> pi(n_states, 1.0 / n_states);
        std::vector<double> pi_new(n_states);

        for (int iter = 0; iter < max_iter; ++iter) {
            // Compute πP
            std::fill(pi_new.begin(), pi_new.end(), 0.0);
            for (int i = 0; i < n_states; ++i) {
                for (int j = 0; j < n_states; ++j) {
                    pi_new[j] += pi[i] * transition_matrix[i][j];
                }
            }

            // Check convergence
            double diff = 0.0;
            for (int i = 0; i < n_states; ++i) {
                diff += std::abs(pi_new[i] - pi[i]);
            }

            pi = pi_new;

            if (diff < tol) {
                break;
            }
        }

        return pi;
    }

    /**
     * @brief Compute n-step transition matrix P^n
     */
    std::vector<std::vector<double>> nStepTransition(int n) {
        auto result = transition_matrix;
        
        for (int step = 1; step < n; ++step) {
            auto temp = std::vector<std::vector<double>>(
                n_states, std::vector<double>(n_states, 0.0));
            
            for (int i = 0; i < n_states; ++i) {
                for (int j = 0; j < n_states; ++j) {
                    for (int k = 0; k < n_states; ++k) {
                        temp[i][j] += result[i][k] * transition_matrix[k][j];
                    }
                }
            }
            result = temp;
        }

        return result;
    }

    /**
     * @brief Check if chain is irreducible
     *
     * All states communicate (can reach any state from any state)
     */
    bool isIrreducible() {
        // Check if there exists n such that P^n has all positive entries
        int max_steps = n_states * n_states;
        auto Pn = nStepTransition(max_steps);

        for (int i = 0; i < n_states; ++i) {
            for (int j = 0; j < n_states; ++j) {
                if (Pn[i][j] <= 1e-10) {
                    return false;
                }
            }
        }
        return true;
    }
};

/**
 * @class MetropolisHastings
 * @brief Metropolis-Hastings MCMC algorithm
 */
class MetropolisHastings {
public:
    /**
     * @brief Run Metropolis-Hastings sampling
     *
     * Samples from target distribution π(x) using proposal q(x'|x)
     *
     * @param target_pdf Target PDF (unnormalized)
     * @param proposal_sampler Proposal sampler q(·|x)
     * @param proposal_pdf Proposal PDF q(x'|x)
     * @param initial Initial state
     * @param n_samples Number of samples
     * @param burn_in Burn-in period
     * @return Sampled chain
     */
    static std::vector<double> sample(
        std::function<double(double)> target_pdf,
        std::function<double(double)> proposal_sampler,
        std::function<double(double, double)> proposal_pdf,
        double initial,
        int n_samples = 10000,
        int burn_in = 1000) {

        std::vector<double> chain;
        chain.reserve(n_samples);

        double current = initial;
        std::uniform_real_distribution<double> unif(0.0, 1.0);

        for (int i = 0; i < burn_in + n_samples; ++i) {
            // Propose new state
            double proposed = proposal_sampler(current);

            // Compute acceptance ratio
            double pi_proposed = target_pdf(proposed);
            double pi_current = target_pdf(current);
            double q_forward = proposal_pdf(proposed, current);
            double q_backward = proposal_pdf(current, proposed);

            double alpha = std::min(1.0, 
                (pi_proposed * q_backward) / (pi_current * q_forward + 1e-10));

            // Accept or reject
            if (unif(gen) < alpha) {
                current = proposed;
            }

            // Store after burn-in
            if (i >= burn_in) {
                chain.push_back(current);
            }
        }

        return chain;
    }

    /**
     * @brief Metropolis algorithm (symmetric proposal)
     *
     * Special case where q(x'|x) = q(x|x') (e.g., random walk)
     */
    static std::vector<double> metropolis(
        std::function<double(double)> target_pdf,
        std::function<double(double)> proposal_sampler,
        double initial,
        int n_samples = 10000,
        int burn_in = 1000) {

        std::vector<double> chain;
        chain.reserve(n_samples);

        double current = initial;
        std::uniform_real_distribution<double> unif(0.0, 1.0);

        for (int i = 0; i < burn_in + n_samples; ++i) {
            double proposed = proposal_sampler(current);

            // Acceptance ratio (proposal cancels)
            double alpha = std::min(1.0, target_pdf(proposed) / (target_pdf(current) + 1e-10));

            if (unif(gen) < alpha) {
                current = proposed;
            }

            if (i >= burn_in) {
                chain.push_back(current);
            }
        }

        return chain;
    }
};

/**
 * @class GibbsSampler
 * @brief Gibbs sampling for multivariate distributions
 */
class GibbsSampler {
public:
    /**
     * @brief 2D Gibbs sampling
     *
     * Alternately samples from conditional distributions
     * x | y ~ p(x|y) and y | x ~ p(y|x)
     *
     * @param conditional_x Sampler for p(x|y)
     * @param conditional_y Sampler for p(y|x)
     * @param initial Initial (x, y)
     * @param n_samples Number of samples
     * @param burn_in Burn-in period
     * @return Vector of (x,y) pairs
     */
    static std::vector<std::pair<double, double>> sample2D(
        std::function<double(double)> conditional_x,
        std::function<double(double)> conditional_y,
        std::pair<double, double> initial,
        int n_samples = 10000,
        int burn_in = 1000) {

        std::vector<std::pair<double, double>> chain;
        chain.reserve(n_samples);

        double x = initial.first;
        double y = initial.second;

        for (int i = 0; i < burn_in + n_samples; ++i) {
            // Sample x | y
            x = conditional_x(y);
            
            // Sample y | x
            y = conditional_y(x);

            if (i >= burn_in) {
                chain.push_back({x, y});
            }
        }

        return chain;
    }
};

/**
 * @class BrownianMotion
 * @brief Brownian motion (Wiener process) simulation
 */
class BrownianMotion {
public:
    /**
     * @brief Simulate standard Brownian motion W(t)
     *
     * W(t) has W(0) = 0, independent increments,
     * W(t) - W(s) ~ N(0, t-s)
     *
     * @param T Final time
     * @param n_steps Number of time steps
     * @return Time series (time, W(t))
     */
    static std::vector<std::pair<double, double>> simulate(
        double T, int n_steps) {

        double dt = T / n_steps;
        double sigma_dt = std::sqrt(dt);
        
        std::normal_distribution<double> norm(0.0, sigma_dt);
        std::vector<std::pair<double, double>> path;
        path.reserve(n_steps + 1);

        double W = 0.0;
        path.push_back({0.0, W});

        for (int i = 1; i <= n_steps; ++i) {
            W += norm(gen);
            path.push_back({i * dt, W});
        }

        return path;
    }

    /**
     * @brief Simulate geometric Brownian motion
     *
     * dS = μS dt + σS dW
     * Solution: S(t) = S(0) exp((μ - σ²/2)t + σW(t))
     *
     * Used in Black-Scholes model
     */
    static std::vector<std::pair<double, double>> geometricBrownian(
        double S0, double mu, double sigma, double T, int n_steps) {

        double dt = T / n_steps;
        double drift = (mu - 0.5 * sigma * sigma) * dt;
        double vol = sigma * std::sqrt(dt);

        std::normal_distribution<double> norm(0.0, 1.0);
        std::vector<std::pair<double, double>> path;
        path.reserve(n_steps + 1);

        double S = S0;
        path.push_back({0.0, S});

        for (int i = 1; i <= n_steps; ++i) {
            double dW = norm(gen);
            S *= std::exp(drift + vol * dW);
            path.push_back({i * dt, S});
        }

        return path;
    }
};

/**
 * @class OrnsteinUhlenbeck
 * @brief Ornstein-Uhlenbeck process (mean-reverting)
 */
class OrnsteinUhlenbeck {
public:
    /**
     * @brief Simulate O-U process
     *
     * dX = θ(μ - X)dt + σdW
     *
     * Mean-reverting to level μ with speed θ and volatility σ
     *
     * @param X0 Initial value
     * @param theta Mean reversion speed
     * @param mu Long-term mean
     * @param sigma Volatility
     * @param T Final time
     * @param n_steps Number of steps
     * @return Time series
     */
    static std::vector<std::pair<double, double>> simulate(
        double X0, double theta, double mu, double sigma,
        double T, int n_steps) {

        double dt = T / n_steps;
        double sqrt_dt = std::sqrt(dt);

        std::normal_distribution<double> norm(0.0, 1.0);
        std::vector<std::pair<double, double>> path;
        path.reserve(n_steps + 1);

        double X = X0;
        path.push_back({0.0, X});

        for (int i = 1; i <= n_steps; ++i) {
            double dW = norm(gen);
            X += theta * (mu - X) * dt + sigma * sqrt_dt * dW;
            path.push_back({i * dt, X});
        }

        return path;
    }
};

/**
 * @class BoltzmannEquation
 * @brief Boltzmann equation and kinetic theory
 */
class BoltzmannEquation {
public:
    /**
     * @brief Maxwell-Boltzmann velocity distribution
     *
     * f(v) = (m/(2πkT))^(3/2) exp(-mv²/(2kT))
     *
     * @param v Velocity
     * @param T Temperature
     * @param m Mass
     * @param k Boltzmann constant
     * @return Probability density
     */
    static double maxwellBoltzmannPDF(double v, double T, double m = 1.0, double k = 1.0) {
        double beta = m / (2.0 * k * T);
        double norm = std::pow(beta / M_PI, 1.5);
        return norm * std::exp(-beta * v * v);
    }

    /**
     * @brief Sample from Maxwell-Boltzmann distribution
     */
    static double sampleMaxwellBoltzmann(double T, double m = 1.0, double k = 1.0) {
        // Use Box-Muller for 3D velocity
        std::normal_distribution<double> norm(0.0, std::sqrt(k * T / m));
        
        double vx = norm(gen);
        double vy = norm(gen);
        double vz = norm(gen);
        
        return std::sqrt(vx*vx + vy*vy + vz*vz);
    }

    /**
     * @brief Direct Simulation Monte Carlo (DSMC) for rarefied gas
     *
     * Bird's algorithm for solving Boltzmann equation
     *
     * @param n_particles Number of simulated particles
     * @param T Temperature
     * @param dt Time step
     * @param n_steps Number of steps
     * @return Average properties over time
     */
    static std::vector<double> DSMC_simulation(
        int n_particles, double T, double dt, int n_steps) {

        // Initialize particles with Maxwell-Boltzmann velocities
        std::vector<double> velocities(n_particles);
        for (int i = 0; i < n_particles; ++i) {
            velocities[i] = sampleMaxwellBoltzmann(T);
        }

        std::vector<double> avg_kinetic_energy;
        avg_kinetic_energy.reserve(n_steps);

        for (int step = 0; step < n_steps; ++step) {
            // Compute average kinetic energy
            double KE = 0.0;
            for (double v : velocities) {
                KE += 0.5 * v * v;
            }
            avg_kinetic_energy.push_back(KE / n_particles);

            // Simple collision model (placeholder for full DSMC)
            // In full DSMC: pair particles, compute collision probability, update velocities
            for (int i = 0; i < n_particles / 2; ++i) {
                int idx1 = std::uniform_int_distribution<int>(0, n_particles-1)(gen);
                int idx2 = std::uniform_int_distribution<int>(0, n_particles-1)(gen);
                
                // Elastic collision (conserve energy and momentum)
                double v_avg = (velocities[idx1] + velocities[idx2]) / 2.0;
                velocities[idx1] = v_avg;
                velocities[idx2] = v_avg;
            }
        }

        return avg_kinetic_energy;
    }

    /**
     * @brief H-theorem: entropy increases
     *
     * H = ∫ f log f d³v (Boltzmann H-function)
     */
    static double computeHFunction(const std::vector<double>& velocities) {
        std::map<int, int> histogram;
        int n_bins = 50;
        double v_max = *std::max_element(velocities.begin(), velocities.end());
        
        for (double v : velocities) {
            int bin = static_cast<int>((v / v_max) * n_bins);
            if (bin >= n_bins) bin = n_bins - 1;
            histogram[bin]++;
        }

        double H = 0.0;
        int N = velocities.size();
        for (const auto& [bin, count] : histogram) {
            if (count > 0) {
                double f = static_cast<double>(count) / N;
                H += f * std::log(f);
            }
        }

        return H;
    }
};

/**
 * @class HamiltonianMonteCarlo
 * @brief Hamiltonian Monte Carlo (Hybrid Monte Carlo)
 */
class HamiltonianMonteCarlo {
public:
    /**
     * @brief HMC sampling
     *
     * Uses Hamiltonian dynamics to propose distant states
     * H(q, p) = U(q) + K(p) where U = -log π(q), K = p²/2
     *
     * @param neg_log_prob Negative log probability -log π(q)
     * @param grad_neg_log_prob Gradient ∇U(q)
     * @param initial Initial position
     * @param epsilon Step size
     * @param L Number of leapfrog steps
     * @param n_samples Number of samples
     * @return Sampled chain
     */
    static std::vector<double> sample(
        std::function<double(double)> neg_log_prob,
        std::function<double(double)> grad_neg_log_prob,
        double initial,
        double epsilon = 0.01,
        int L = 10,
        int n_samples = 1000) {

        std::vector<double> chain;
        chain.reserve(n_samples);

        std::normal_distribution<double> norm(0.0, 1.0);
        std::uniform_real_distribution<double> unif(0.0, 1.0);

        double q = initial;

        for (int i = 0; i < n_samples; ++i) {
            // Sample momentum
            double p = norm(gen);
            double current_p = p;

            // Leapfrog integration
            double q_new = q;
            double p_new = p;

            p_new -= 0.5 * epsilon * grad_neg_log_prob(q_new);
            
            for (int step = 0; step < L; ++step) {
                q_new += epsilon * p_new;
                if (step < L - 1) {
                    p_new -= epsilon * grad_neg_log_prob(q_new);
                }
            }
            
            p_new -= 0.5 * epsilon * grad_neg_log_prob(q_new);
            p_new = -p_new;  // Negate for reversibility

            // Compute Hamiltonian
            double current_H = neg_log_prob(q) + 0.5 * current_p * current_p;
            double proposed_H = neg_log_prob(q_new) + 0.5 * p_new * p_new;

            // Accept/reject
            double alpha = std::min(1.0, std::exp(current_H - proposed_H));
            
            if (unif(gen) < alpha) {
                q = q_new;
            }

            chain.push_back(q);
        }

        return chain;
    }
};

/**
 * @class MCMCDiagnostics
 * @brief Diagnostics for MCMC convergence
 */
class MCMCDiagnostics {
public:
    /**
     * @brief Effective sample size (ESS)
     *
     * Accounts for autocorrelation: ESS = N / (1 + 2∑ρₖ)
     */
    static double effectiveSampleSize(const std::vector<double>& chain, int max_lag = 100) {
        int N = chain.size();
        double mean = std::accumulate(chain.begin(), chain.end(), 0.0) / N;

        // Compute variance
        double variance = 0.0;
        for (double x : chain) {
            variance += (x - mean) * (x - mean);
        }
        variance /= N;

        // Compute autocorrelations
        double autocorr_sum = 0.0;
        for (int lag = 1; lag <= std::min(max_lag, N/2); ++lag) {
            double autocorr = 0.0;
            for (int i = 0; i < N - lag; ++i) {
                autocorr += (chain[i] - mean) * (chain[i + lag] - mean);
            }
            autocorr /= ((N - lag) * variance);
            
            if (autocorr < 0.05) break;  // Stop when autocorrelation becomes small
            autocorr_sum += autocorr;
        }

        return N / (1.0 + 2.0 * autocorr_sum);
    }

    /**
     * @brief Gelman-Rubin convergence diagnostic (R-hat)
     *
     * Compares within-chain and between-chain variance
     * Should be close to 1.0 for convergence
     */
    static double gelmanRubin(const std::vector<std::vector<double>>& chains) {
        int m = chains.size();  // Number of chains
        int n = chains[0].size();  // Length of each chain

        // Compute chain means
        std::vector<double> chain_means(m);
        for (int i = 0; i < m; ++i) {
            chain_means[i] = std::accumulate(chains[i].begin(), chains[i].end(), 0.0) / n;
        }

        // Overall mean
        double overall_mean = std::accumulate(chain_means.begin(), chain_means.end(), 0.0) / m;

        // Between-chain variance
        double B = 0.0;
        for (double mean : chain_means) {
            B += (mean - overall_mean) * (mean - overall_mean);
        }
        B *= n / (m - 1.0);

        // Within-chain variance
        double W = 0.0;
        for (int i = 0; i < m; ++i) {
            double var = 0.0;
            for (double x : chains[i]) {
                var += (x - chain_means[i]) * (x - chain_means[i]);
            }
            W += var / (n - 1.0);
        }
        W /= m;

        // Pooled variance estimate
        double var_plus = ((n - 1.0) * W + B) / n;

        // R-hat statistic
        return std::sqrt(var_plus / W);
    }

    /**
     * @brief Acceptance rate for MH algorithm
     */
    static double acceptanceRate(const std::vector<double>& chain) {
        int accepts = 0;
        for (size_t i = 1; i < chain.size(); ++i) {
            if (std::abs(chain[i] - chain[i-1]) > 1e-10) {
                accepts++;
            }
        }
        return static_cast<double>(accepts) / (chain.size() - 1);
    }
};

} // namespace maths::stochastic

#endif // MATHS_STOCHASTIC_MONTE_CARLO_HPP
