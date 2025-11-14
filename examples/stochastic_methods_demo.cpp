/**
 * @file stochastic_methods_demo.cpp
 * @brief Comprehensive demonstration of stochastic methods and Monte Carlo algorithms
 *
 * Demonstrates:
 * - Monte Carlo integration (basic, importance, rejection, stratified)
 * - Markov chains (discrete-time, transition matrices, stationary distributions)
 * - MCMC sampling (Metropolis-Hastings, Gibbs sampling)
 * - Stochastic processes (Brownian motion, geometric Brownian, Ornstein-Uhlenbeck)
 * - Boltzmann equation and kinetic theory (DSMC, H-theorem)
 * - Hamiltonian Monte Carlo
 * - MCMC convergence diagnostics
 */

#include "maths/stochastic/monte_carlo.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace maths::stochastic;

void printSeparator(const std::string& title) {
    std::cout << "\n=== " << title << " ===" << std::endl;
}

// Demo 1: Monte Carlo Integration
void demoMonteCarloIntegration() {
    printSeparator("Demo 1: Monte Carlo Integration");

    // Integrate sin(x) from 0 to π (exact answer = 2)
    auto f = [](double x) { return std::sin(x); };
    double a = 0.0, b = M_PI;
    int n_samples = 10000;

    std::cout << "Integrating sin(x) from 0 to π (exact = 2.0)" << std::endl;
    std::cout << "Number of samples: " << n_samples << "\n" << std::endl;

    // Basic Monte Carlo
    auto [result_basic, error_basic] = MonteCarloIntegrator::integrate(f, a, b, n_samples);
    std::cout << "Basic Monte Carlo:" << std::endl;
    std::cout << "  Result: " << std::setprecision(6) << result_basic << std::endl;
    std::cout << "  Error estimate: ± " << error_basic << std::endl;
    std::cout << "  Actual error: " << std::abs(result_basic - 2.0) << std::endl;

    // Stratified sampling (better convergence)
    auto [result_strat, error_strat] = MonteCarloIntegrator::stratifiedSampling(f, a, b, n_samples, 10);
    std::cout << "\nStratified Sampling (10 strata):" << std::endl;
    std::cout << "  Result: " << result_strat << std::endl;
    std::cout << "  Error estimate: ± " << error_strat << std::endl;
    std::cout << "  Actual error: " << std::abs(result_strat - 2.0) << std::endl;
    std::cout << "  (Reduced variance compared to basic MC)" << std::endl;
}

// Demo 2: Discrete Markov Chains
void demoMarkovChains() {
    printSeparator("Demo 2: Discrete-Time Markov Chains");

    // Simple 3-state weather model: Sunny, Cloudy, Rainy
    std::vector<std::vector<double>> P = {
        {0.7, 0.2, 0.1},  // From Sunny
        {0.3, 0.4, 0.3},  // From Cloudy
        {0.2, 0.3, 0.5}   // From Rainy
    };

    MarkovChain weather(P);

    std::cout << "Weather Model (3 states: Sunny, Cloudy, Rainy)" << std::endl;
    std::cout << "\nTransition Matrix:" << std::endl;
    std::cout << "       S    C    R" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << (i == 0 ? "  S " : (i == 1 ? "  C " : "  R ")) << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << std::setw(4) << P[i][j] << " ";
        }
        std::cout << "]" << std::endl;
    }

    // Simulate chain
    auto states = weather.simulate(0, 10);  // Start sunny, 10 steps
    std::cout << "\nSimulation (starting from Sunny):" << std::endl;
    std::cout << "  States: ";
    const char* labels[] = {"S", "C", "R"};
    for (int s : states) {
        std::cout << labels[s] << " ";
    }
    std::cout << std::endl;

    // Stationary distribution
    auto stationary = weather.stationaryDistribution();
    std::cout << "\nStationary Distribution (long-run probabilities):" << std::endl;
    std::cout << "  Sunny:  " << std::setprecision(4) << stationary[0] << std::endl;
    std::cout << "  Cloudy: " << stationary[1] << std::endl;
    std::cout << "  Rainy:  " << stationary[2] << std::endl;

    // Check irreducibility
    bool irreducible = weather.isIrreducible();
    std::cout << "\nIs irreducible? " << (irreducible ? "Yes" : "No") << std::endl;
}

// Demo 3: Metropolis-Hastings MCMC
void demoMetropolisHastings() {
    printSeparator("Demo 3: Metropolis-Hastings MCMC");

    // Target: Standard normal distribution
    auto target_pdf = [](double x) {
        return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI);
    };

    std::cout << "Sampling from standard normal N(0,1)" << std::endl;
    std::cout << "Using symmetric proposal (random walk)" << std::endl;

    // Proposal: random walk with fixed step size
    double sigma = 1.0;
    auto proposal_sampler = [sigma](double x) {
        static std::mt19937 rng(std::random_device{}());
        std::normal_distribution<double> dist(x, sigma);
        return dist(rng);
    };

    // Sample using symmetric Metropolis
    auto samples = MetropolisHastings::metropolis(target_pdf, proposal_sampler, 0.0, 5000, 1000);

    // Compute sample statistics
    double mean = 0.0;
    for (double x : samples) mean += x;
    mean /= samples.size();

    double variance = 0.0;
    for (double x : samples) variance += (x - mean) * (x - mean);
    variance /= samples.size();

    std::cout << "\nResults (after burn-in of 1000):" << std::endl;
    std::cout << "  Samples collected: " << samples.size() << std::endl;
    std::cout << "  Sample mean: " << std::setprecision(4) << mean << " (expected 0)" << std::endl;
    std::cout << "  Sample variance: " << variance << " (expected 1)" << std::endl;

    // Diagnostics
    double acceptance = MCMCDiagnostics::acceptanceRate(samples);
    double ess = MCMCDiagnostics::effectiveSampleSize(samples, 50);

    std::cout << "\nMCMC Diagnostics:" << std::endl;
    std::cout << "  Acceptance rate: " << acceptance << " (target ~0.234 for 1D)" << std::endl;
    std::cout << "  Effective sample size: " << ess << std::endl;
}

// Demo 4: Gibbs Sampling
void demoGibbsSampling() {
    printSeparator("Demo 4: Gibbs Sampling for 2D Distribution");

    // Sample from bivariate normal (independent components for simplicity)
    auto conditional_x = [](double y) {
        static std::mt19937 rng(std::random_device{}());
        std::normal_distribution<double> dist(0.0, 1.0);
        return dist(rng);
    };

    auto conditional_y = [](double x) {
        static std::mt19937 rng(std::random_device{}());
        std::normal_distribution<double> dist(0.0, 1.0);
        return dist(rng);
    };

    std::cout << "Sampling from 2D standard normal" << std::endl;
    std::cout << "Using Gibbs sampling (alternating conditionals)" << std::endl;

    auto samples = GibbsSampler::sample2D(conditional_x, conditional_y, {0.0, 0.0}, 5000, 1000);

    // Compute sample statistics
    double mean_x = 0.0, mean_y = 0.0;
    for (auto [x, y] : samples) {
        mean_x += x;
        mean_y += y;
    }
    mean_x /= samples.size();
    mean_y /= samples.size();

    std::cout << "\nResults:" << std::endl;
    std::cout << "  Samples collected: " << samples.size() << std::endl;
    std::cout << "  Mean X: " << std::setprecision(4) << mean_x << " (expected 0)" << std::endl;
    std::cout << "  Mean Y: " << mean_y << " (expected 0)" << std::endl;
}

// Demo 5: Brownian Motion
void demoBrownianMotion() {
    printSeparator("Demo 5: Brownian Motion and Stochastic Processes");

    double T = 1.0;  // Time horizon
    int n_steps = 1000;

    // Standard Brownian motion (Wiener process)
    auto wiener = BrownianMotion::simulate(T, n_steps);

    std::cout << "Standard Brownian Motion W(t):" << std::endl;
    std::cout << "  Time horizon: " << T << std::endl;
    std::cout << "  Number of steps: " << n_steps << std::endl;
    std::cout << "  Final position W(" << T << "): " << wiener.back().second << std::endl;
    std::cout << "  (Expected variance: " << T << ")" << std::endl;

    // Geometric Brownian motion (stock prices)
    double S0 = 100.0;  // Initial price
    double mu = 0.05;   // Drift (5% annual return)
    double sigma = 0.2; // Volatility (20%)

    auto gbm = BrownianMotion::geometricBrownian(S0, mu, sigma, T, n_steps);

    std::cout << "\nGeometric Brownian Motion (Stock Price Model):" << std::endl;
    std::cout << "  Initial price S(0): $" << S0 << std::endl;
    std::cout << "  Drift μ: " << mu * 100 << "%" << std::endl;
    std::cout << "  Volatility σ: " << sigma * 100 << "%" << std::endl;
    std::cout << "  Final price S(" << T << "): $" << std::setprecision(2) << std::fixed
              << gbm.back().second << std::endl;
    std::cout << std::defaultfloat << std::setprecision(6);

    // Ornstein-Uhlenbeck (mean-reverting)
    double X0 = 5.0;      // Initial value
    double theta = 2.0;   // Mean reversion speed
    double mu_ou = 0.0;   // Long-term mean
    double sigma_ou = 1.0; // Volatility

    auto ou = OrnsteinUhlenbeck::simulate(X0, theta, mu_ou, sigma_ou, T, n_steps);

    std::cout << "\nOrnstein-Uhlenbeck Process (Mean-Reverting):" << std::endl;
    std::cout << "  Initial value: " << X0 << std::endl;
    std::cout << "  Long-term mean: " << mu_ou << std::endl;
    std::cout << "  Mean reversion speed θ: " << theta << std::endl;
    std::cout << "  Final value: " << std::setprecision(4) << ou.back().second << std::endl;
    std::cout << "  (Process reverts to mean " << mu_ou << ")" << std::endl;
}

// Demo 6: Boltzmann Equation and Kinetic Theory
void demoBoltzmannEquation() {
    printSeparator("Demo 6: Boltzmann Equation and Kinetic Theory");

    double T = 300.0;  // Temperature (K)
    double m = 1.0;    // Particle mass (arbitrary units)
    double k = 1.0;    // Boltzmann constant (arbitrary units)

    std::cout << "Maxwell-Boltzmann Velocity Distribution" << std::endl;
    std::cout << "  Temperature T: " << T << " K" << std::endl;
    std::cout << "  Particle mass m: " << m << std::endl;

    // Sample velocities
    std::cout << "\nSampled velocities:" << std::endl;
    for (int i = 0; i < 5; i++) {
        double v = BoltzmannEquation::sampleMaxwellBoltzmann(T, m, k);
        double pdf = BoltzmannEquation::maxwellBoltzmannPDF(v, T, m, k);
        std::cout << "  v = " << std::setw(8) << std::setprecision(4) << v
                  << "  (PDF = " << pdf << ")" << std::endl;
    }

    // DSMC simulation
    std::cout << "\nDirect Simulation Monte Carlo (DSMC):" << std::endl;
    int n_particles = 1000;
    double dt = 0.01;
    int n_steps = 100;

    auto velocities = BoltzmannEquation::DSMC_simulation(n_particles, T, dt, n_steps);

    // Compute statistics
    double mean_v = 0.0;
    for (double v : velocities) mean_v += v;
    mean_v /= velocities.size();

    double mean_v2 = 0.0;
    for (double v : velocities) mean_v2 += v * v;
    mean_v2 /= velocities.size();

    std::cout << "  Particles: " << n_particles << std::endl;
    std::cout << "  Time steps: " << n_steps << std::endl;
    std::cout << "  Mean velocity: " << mean_v << std::endl;
    std::cout << "  Mean v²: " << mean_v2 << std::endl;

    // H-theorem
    double H = BoltzmannEquation::computeHFunction(velocities);
    std::cout << "\nH-function (entropy): " << H << std::endl;
    std::cout << "  (H decreases toward equilibrium - H-theorem)" << std::endl;
}

// Demo 7: Hamiltonian Monte Carlo
void demoHamiltonianMonteCarlo() {
    printSeparator("Demo 7: Hamiltonian Monte Carlo (HMC)");

    // Sample from standard normal using HMC
    auto neg_log_prob = [](double x) {
        return 0.5 * x * x;  // -log(N(0,1)) ∝ x²/2
    };

    auto grad_neg_log_prob = [](double x) {
        return x;  // gradient
    };

    std::cout << "Sampling from standard normal N(0,1)" << std::endl;
    std::cout << "Using Hamiltonian Monte Carlo (no random walk!)" << std::endl;

    double epsilon = 0.1;  // Step size
    int L = 10;            // Number of leapfrog steps
    int n_samples = 2000;

    auto samples = HamiltonianMonteCarlo::sample(
        neg_log_prob, grad_neg_log_prob, 0.0, epsilon, L, n_samples);

    // Statistics
    double mean = 0.0;
    for (double x : samples) mean += x;
    mean /= samples.size();

    double variance = 0.0;
    for (double x : samples) variance += (x - mean) * (x - mean);
    variance /= samples.size();

    std::cout << "\nHMC Parameters:" << std::endl;
    std::cout << "  Step size ε: " << epsilon << std::endl;
    std::cout << "  Leapfrog steps L: " << L << std::endl;

    std::cout << "\nResults:" << std::endl;
    std::cout << "  Samples: " << samples.size() << std::endl;
    std::cout << "  Sample mean: " << std::setprecision(4) << mean << " (expected 0)" << std::endl;
    std::cout << "  Sample variance: " << variance << " (expected 1)" << std::endl;

    // Diagnostics
    double acceptance = MCMCDiagnostics::acceptanceRate(samples);
    double ess = MCMCDiagnostics::effectiveSampleSize(samples, 50);

    std::cout << "\nDiagnostics:" << std::endl;
    std::cout << "  Acceptance rate: " << acceptance << " (HMC typically ~0.65-0.9)" << std::endl;
    std::cout << "  Effective sample size: " << ess << std::endl;
    std::cout << "  ESS/N ratio: " << ess / samples.size() << " (HMC has low autocorr)" << std::endl;
}

// Demo 8: MCMC Convergence Diagnostics
void demoMCMCDiagnostics() {
    printSeparator("Demo 8: MCMC Convergence Diagnostics");

    // Run multiple chains for Gelman-Rubin diagnostic
    auto target_pdf = [](double x) {
        return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI);
    };

    std::cout << "Running 4 independent MCMC chains" << std::endl;
    std::cout << "Target: Standard normal N(0,1)" << std::endl;

    std::vector<std::vector<double>> chains;
    std::vector<double> starting_points = {-2.0, -1.0, 1.0, 2.0};

    // Proposal: random walk
    double sigma = 1.0;
    auto proposal_sampler = [sigma](double x) {
        static std::mt19937 rng(std::random_device{}());
        std::normal_distribution<double> dist(x, sigma);
        return dist(rng);
    };

    for (double start : starting_points) {
        auto chain = MetropolisHastings::metropolis(target_pdf, proposal_sampler, start, 3000, 500);
        chains.push_back(chain);
    }

    std::cout << "\nChain Statistics:" << std::endl;
    for (size_t i = 0; i < chains.size(); i++) {
        double mean = 0.0;
        for (double x : chains[i]) mean += x;
        mean /= chains[i].size();

        double ess = MCMCDiagnostics::effectiveSampleSize(chains[i], 50);

        std::cout << "  Chain " << i + 1 << ": mean = " << std::setprecision(4) << mean
                  << ", ESS = " << ess << std::endl;
    }

    // Gelman-Rubin diagnostic
    double R_hat = MCMCDiagnostics::gelmanRubin(chains);
    std::cout << "\nGelman-Rubin R̂ statistic: " << R_hat << std::endl;
    std::cout << "  (R̂ < 1.1 indicates convergence)" << std::endl;
    std::cout << "  Status: " << (R_hat < 1.1 ? "CONVERGED ✓" : "NOT CONVERGED ✗") << std::endl;
}

// Demo 9: Applications
void demoApplications() {
    printSeparator("Demo 9: Practical Applications");

    std::cout << "Monte Carlo Methods in Science and Engineering:\n" << std::endl;

    std::cout << "1. Financial Mathematics:" << std::endl;
    std::cout << "   - Option pricing (Black-Scholes via GBM)" << std::endl;
    std::cout << "   - Risk assessment (VaR, CVaR)" << std::endl;
    std::cout << "   - Portfolio optimization" << std::endl;

    std::cout << "\n2. Statistical Physics:" << std::endl;
    std::cout << "   - Boltzmann equation for gas dynamics" << std::endl;
    std::cout << "   - Ising model simulations (Metropolis)" << std::endl;
    std::cout << "   - Molecular dynamics (DSMC)" << std::endl;

    std::cout << "\n3. Bayesian Inference:" << std::endl;
    std::cout << "   - Parameter estimation (MCMC)" << std::endl;
    std::cout << "   - Model comparison (Bayes factors)" << std::endl;
    std::cout << "   - Hierarchical models (Gibbs sampling)" << std::endl;

    std::cout << "\n4. Machine Learning:" << std::endl;
    std::cout << "   - Bayesian neural networks (HMC)" << std::endl;
    std::cout << "   - Variational inference" << std::endl;
    std::cout << "   - Reinforcement learning (policy gradients)" << std::endl;

    std::cout << "\n5. Numerical Integration:" << std::endl;
    std::cout << "   - High-dimensional integrals" << std::endl;
    std::cout << "   - Path integrals (quantum mechanics)" << std::endl;
    std::cout << "   - Expectation computation" << std::endl;
}

int main() {
    std::cout << "============================================" << std::endl;
    std::cout << "Stochastic Methods and Monte Carlo Demo" << std::endl;
    std::cout << "Comprehensive Computational Statistics" << std::endl;
    std::cout << "============================================" << std::endl;

    try {
        demoMonteCarloIntegration();
        demoMarkovChains();
        demoMetropolisHastings();
        demoGibbsSampling();
        demoBrownianMotion();
        demoBoltzmannEquation();
        demoHamiltonianMonteCarlo();
        demoMCMCDiagnostics();
        demoApplications();

        std::cout << "\n============================================" << std::endl;
        std::cout << "All stochastic methods demos completed!" << std::endl;
        std::cout << "============================================" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
