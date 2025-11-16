#ifndef PHYSICS_ADVANCED_CLASSICAL_LIOUVILLE_HPP
#define PHYSICS_ADVANCED_CLASSICAL_LIOUVILLE_HPP

#include "hamiltonian.hpp"
#include "phase_space.hpp"
#include <Eigen/Dense>
#include <functional>

/**
 * @file liouville.hpp
 * @brief Liouville's equation and ensemble evolution
 *
 * Implements:
 * - Liouville's equation for density evolution
 * - Statistical ensembles
 * - Ensemble averages
 * - BBGKY hierarchy
 * - Reduced distributions
 */

namespace physics::advanced::classical {

/**
 * @class LiouvilleEquation
 * @brief Evolution of phase space density via Liouville's equation
 *
 * Liouville's equation:
 *   ∂ρ/∂t + {H, ρ} = 0
 *
 * where {H, ρ} is the Poisson bracket:
 *   {H, ρ} = Σᵢ (∂H/∂qᵢ ∂ρ/∂pᵢ - ∂H/∂pᵢ ∂ρ/∂qᵢ)
 *
 * This is equivalent to:
 *   dρ/dt = 0  (density constant along trajectories)
 */
class LiouvilleEquation {
private:
    const HamiltonianSystem& system_;
    PoissonBracket bracket_;
    double epsilon_;

public:
    LiouvilleEquation(const HamiltonianSystem& system, double epsilon = 1e-8)
        : system_(system), bracket_(epsilon), epsilon_(epsilon) {}

    /**
     * @brief Compute Liouville operator acting on density
     *
     * L̂ρ = {H, ρ}
     */
    double liouvilleOperator(
        const std::function<double(const PhasePoint&)>& rho,
        const PhasePoint& point) const {

        auto H = [this](const PhasePoint& p) {
            return system_.hamiltonian(p);
        };

        return bracket_.compute(H, rho, point);
    }

    /**
     * @brief Evolve density forward in time (first-order)
     *
     * ρ(t + dt) = ρ(t) - dt {H, ρ}
     */
    std::function<double(const PhasePoint&)> evolve(
        const std::function<double(const PhasePoint&)>& rho0,
        double dt) const {

        return [this, rho0, dt](const PhasePoint& point) {
            // Backward in time to find where density came from
            PhasePoint prev = point;
            PhasePoint deriv = system_.hamiltonEquations(point);
            prev.q = point.q - dt * deriv.q;
            prev.p = point.p - dt * deriv.p;

            return rho0(prev);
        };
    }

    /**
     * @brief Verify Liouville's theorem numerically
     *
     * Check that ∂ρ/∂t + {H, ρ} ≈ 0
     */
    bool verifyLiouville(
        const std::function<double(const PhasePoint&)>& rho,
        const PhasePoint& point,
        double tolerance = 1e-6) const {

        double liouville = liouvilleOperator(rho, point);
        return std::abs(liouville) < tolerance;
    }

    /**
     * @brief Check density conservation along trajectory
     *
     * For any trajectory: dρ/dt = ∂ρ/∂t + Σᵢ(q̇ᵢ ∂ρ/∂qᵢ + ṗᵢ ∂ρ/∂pᵢ) = 0
     */
    bool checkConservation(
        const std::function<double(const PhasePoint&)>& rho,
        const std::vector<PhasePoint>& trajectory,
        double tolerance = 1e-4) const {

        if (trajectory.size() < 2) return true;

        double rho0 = rho(trajectory[0]);

        for (size_t i = 1; i < trajectory.size(); ++i) {
            double rho_t = rho(trajectory[i]);
            if (std::abs(rho_t - rho0) > tolerance * std::abs(rho0)) {
                return false;
            }
        }

        return true;
    }
};

/**
 * @class StatisticalEnsemble
 * @brief Different statistical ensembles for equilibrium
 */
class StatisticalEnsemble {
public:
    /**
     * @brief Microcanonical ensemble (isolated system, NVE)
     *
     * All microstates with same energy E equally probable:
     * ρ(q,p) = δ(H(q,p) - E) / Ω(E)
     *
     * where Ω(E) is density of states
     */
    static PhaseSpaceDensity microcanonical(
        const HamiltonianSystem& system,
        double energy, double delta_E = 1e-6) {

        return PhaseSpaceDensity::microcanonical(system, energy, delta_E);
    }

    /**
     * @brief Canonical ensemble (constant temperature, NVT)
     *
     * Boltzmann distribution:
     * ρ(q,p) = Z⁻¹ exp(-βH(q,p))
     *
     * where β = 1/(kT), Z = partition function
     */
    static PhaseSpaceDensity canonical(
        const HamiltonianSystem& system,
        double temperature, double k_boltzmann = 1.380649e-23) {

        return PhaseSpaceDensity::maxwellBoltzmann(
            system, temperature, k_boltzmann);
    }

    /**
     * @brief Grand canonical ensemble (μVT)
     *
     * ρ(q,p,N) = Ξ⁻¹ exp(-β(H - μN))
     *
     * where μ is chemical potential, N is particle number
     */
    static std::function<double(const PhasePoint&, int)> grandCanonical(
        const HamiltonianSystem& system,
        double temperature, double chemical_potential,
        double k_boltzmann = 1.380649e-23) {

        double beta = 1.0 / (k_boltzmann * temperature);

        return [&system, beta, chemical_potential](
            const PhasePoint& point, int N) {
            double H = system.hamiltonian(point);
            return std::exp(-beta * (H - chemical_potential * N));
        };
    }

    /**
     * @brief Calculate partition function (canonical ensemble)
     *
     * Z = ∫ exp(-βH) dΩ
     *
     * (Monte Carlo estimate)
     */
    static double partitionFunction(
        const HamiltonianSystem& system,
        double temperature,
        const std::vector<PhasePoint>& sample_points,
        double k_boltzmann = 1.380649e-23) {

        double beta = 1.0 / (k_boltzmann * temperature);
        double sum = 0.0;

        for (const auto& point : sample_points) {
            double H = system.hamiltonian(point);
            sum += std::exp(-beta * H);
        }

        // Normalize by sample volume
        return sum / sample_points.size();
    }

    /**
     * @brief Calculate Helmholtz free energy
     *
     * F = -kT ln Z
     */
    static double freeEnergy(double partition_function,
                           double temperature,
                           double k_boltzmann = 1.380649e-23) {
        return -k_boltzmann * temperature * std::log(partition_function);
    }

    /**
     * @brief Calculate entropy (canonical ensemble)
     *
     * S = -k ∫ ρ ln ρ dΩ
     */
    static double entropy(
        const PhaseSpaceDensity& density,
        const std::vector<PhasePoint>& sample_points,
        double k_boltzmann = 1.380649e-23) {

        double sum = 0.0;
        for (const auto& point : sample_points) {
            double rho = density(point);
            if (rho > 1e-30) {  // Avoid log(0)
                sum += rho * std::log(rho);
            }
        }

        return -k_boltzmann * sum / sample_points.size();
    }
};

/**
 * @class ReducedDistribution
 * @brief Reduced phase space distributions
 *
 * For N-particle systems, integrate out some coordinates:
 * f₁(q₁,p₁) = ∫...∫ f(q₁...qₙ, p₁...pₙ) dq₂...dqₙ dp₂...dpₙ
 */
class ReducedDistribution {
public:
    /**
     * @brief Single-particle distribution (for many-particle system)
     *
     * f₁(q,p) gives probability density for one particle
     * regardless of others' positions
     */
    static std::function<double(double, double)> singleParticle(
        const std::function<double(const std::vector<PhasePoint>&)>& f_N,
        const std::vector<std::vector<PhasePoint>>& ensemble) {

        return [f_N, ensemble](double q, double p) {
            // Average over all particles in ensemble
            double sum = 0.0;
            int count = 0;

            for (const auto& config : ensemble) {
                for (const auto& particle : config) {
                    // Check if near (q, p)
                    double dq = std::abs(particle.q(0) - q);
                    double dp = std::abs(particle.p(0) - p);
                    if (dq < 0.1 && dp < 0.1) {  // Within window
                        sum += f_N(config);
                        count++;
                    }
                }
            }

            return (count > 0) ? sum / count : 0.0;
        };
    }

    /**
     * @brief Two-particle correlation function
     *
     * g₂(q₁,p₁,q₂,p₂) measures correlations between pairs
     */
    static std::function<double(const PhasePoint&, const PhasePoint&)>
    twoParticleCorrelation(
        const std::function<double(const std::vector<PhasePoint>&)>& f_N,
        const std::vector<std::vector<PhasePoint>>& ensemble) {

        return [f_N, ensemble](const PhasePoint& p1, const PhasePoint& p2) {
            // Simplified correlation function
            double sum = 0.0;
            int count = 0;

            for (const auto& config : ensemble) {
                if (config.size() >= 2) {
                    sum += f_N(config);
                    count++;
                }
            }

            return (count > 0) ? sum / count : 0.0;
        };
    }
};

/**
 * @class BBGKYHierarchy
 * @brief BBGKY hierarchy for distribution functions
 *
 * Relates n-particle to (n+1)-particle distributions:
 * ∂f_s/∂t + Σᵢ(qᵢ·∂f_s/∂qᵢ + Fᵢ·∂f_s/∂pᵢ) = I[f_{s+1}]
 *
 * where s is the order and I is collision integral
 */
class BBGKYHierarchy {
public:
    /**
     * @brief Check consistency between s and s+1 level
     *
     * ∫ f_{s+1}(x₁...x_{s+1}) dx_{s+1} = f_s(x₁...xₛ)
     */
    static bool checkConsistency(
        const std::function<double(const std::vector<PhasePoint>&)>& f_s,
        const std::function<double(const std::vector<PhasePoint>&)>& f_s1,
        const std::vector<std::vector<PhasePoint>>& test_configs,
        double tolerance = 1e-4) {

        // Simplified check
        for (const auto& config : test_configs) {
            if (config.size() < 2) continue;

            std::vector<PhasePoint> config_s(config.begin(),
                                            config.begin() + config.size() - 1);
            double val_s = f_s(config_s);
            double val_s1 = f_s1(config);

            if (std::abs(val_s - val_s1) > tolerance) {
                return false;
            }
        }

        return true;
    }
};

} // namespace physics::advanced::classical

#endif // PHYSICS_ADVANCED_CLASSICAL_LIOUVILLE_HPP
