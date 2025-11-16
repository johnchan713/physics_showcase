#ifndef PHYSICS_ADVANCED_CLASSICAL_PHASE_SPACE_HPP
#define PHYSICS_ADVANCED_CLASSICAL_PHASE_SPACE_HPP

#include "hamiltonian.hpp"
#include <Eigen/Dense>
#include <functional>
#include <cmath>

/**
 * @file phase_space.hpp
 * @brief Phase space structure and analysis
 *
 * Implements:
 * - Phase space volume
 * - Phase space density
 * - Poincaré sections
 * - Phase portraits
 * - Lyapunov exponents
 * - Symplectic structure
 */

namespace physics::advanced::classical {

/**
 * @class PhaseSpaceVolume
 * @brief Calculate and track phase space volumes
 */
class PhaseSpaceVolume {
public:
    /**
     * @brief Calculate volume element in phase space
     *
 * dΩ = dq₁...dqₙ dp₁...dpₙ
     */
    static double volumeElement(const std::vector<PhasePoint>& points) {
        if (points.empty()) return 0.0;

        int dim = points[0].dimension();
        if (dim == 1) {
            // 2D phase space: area = Δq × Δp
            double q_min = points[0].q(0), q_max = points[0].q(0);
            double p_min = points[0].p(0), p_max = points[0].p(0);

            for (const auto& pt : points) {
                q_min = std::min(q_min, pt.q(0));
                q_max = std::max(q_max, pt.q(0));
                p_min = std::min(p_min, pt.p(0));
                p_max = std::max(p_max, pt.p(0));
            }

            return (q_max - q_min) * (p_max - p_min);
        }

        // For higher dimensions, use convex hull volume (simplified)
        return 0.0;  // Would need more sophisticated algorithm
    }

    /**
     * @brief Verify Liouville's theorem: phase space volume preservation
     *
     * In Hamiltonian systems, phase space volume is conserved:
     * dΩ(t) = dΩ(0)
     */
    static bool verifyLiouvilleTheorem(
        const HamiltonianSystem& system,
        const std::vector<PhasePoint>& initial_region,
        double time, double dt,
        double tolerance = 1e-6) {

        // Integrate all points in region
        int num_steps = static_cast<int>(time / dt);
        std::vector<PhasePoint> final_region;
        final_region.reserve(initial_region.size());

        for (const auto& point : initial_region) {
            auto trajectory = system.integrate(point, dt, num_steps);
            final_region.push_back(trajectory.back());
        }

        // Compare volumes
        double V0 = volumeElement(initial_region);
        double Vt = volumeElement(final_region);

        return std::abs(Vt - V0) < tolerance * V0;
    }
};

/**
 * @class PhaseSpaceDensity
 * @brief Distribution function in phase space ρ(q, p, t)
 */
class PhaseSpaceDensity {
private:
    int dim_;
    std::function<double(const PhasePoint&)> rho_;

public:
    PhaseSpaceDensity(int dimension,
                     std::function<double(const PhasePoint&)> density)
        : dim_(dimension), rho_(density) {}

    /**
     * @brief Evaluate density at phase point
     */
    double operator()(const PhasePoint& point) const {
        return rho_(point);
    }

    /**
     * @brief Maxwell-Boltzmann distribution
     *
     * ρ(q,p) = Z⁻¹ exp(-βH(q,p))
     * where β = 1/(kT) and Z is partition function
     */
    static PhaseSpaceDensity maxwellBoltzmann(
        const HamiltonianSystem& system,
        double temperature, double k_boltzmann = 1.380649e-23) {

        double beta = 1.0 / (k_boltzmann * temperature);

        return PhaseSpaceDensity(system.dimension(),
            [&system, beta](const PhasePoint& point) {
                double H = system.hamiltonian(point);
                return std::exp(-beta * H);
            });
    }

    /**
     * @brief Microcanonical ensemble (constant energy surface)
     *
     * ρ(q,p) = δ(H(q,p) - E) / Ω(E)
     */
    static PhaseSpaceDensity microcanonical(
        const HamiltonianSystem& system,
        double energy, double delta_E = 1e-6) {

        return PhaseSpaceDensity(system.dimension(),
            [&system, energy, delta_E](const PhasePoint& point) {
                double H = system.hamiltonian(point);
                double dE = std::abs(H - energy);
                return (dE < delta_E) ? 1.0 / delta_E : 0.0;
            });
    }

    /**
     * @brief Calculate average of observable f
     *
     * ⟨f⟩ = ∫ f(q,p) ρ(q,p) dΩ / ∫ ρ(q,p) dΩ
     *
     * (Simplified Monte Carlo estimate)
     */
    double average(const std::function<double(const PhasePoint&)>& observable,
                  const std::vector<PhasePoint>& sample_points) const {
        double sum_f = 0.0;
        double sum_rho = 0.0;

        for (const auto& point : sample_points) {
            double rho_val = rho_(point);
            sum_f += observable(point) * rho_val;
            sum_rho += rho_val;
        }

        return sum_f / sum_rho;
    }
};

/**
 * @class PoincareSection
 * @brief Poincaré surface of section analysis
 *
 * Reduces continuous dynamics to discrete map by recording
 * crossings of a surface in phase space
 */
class PoincareSection {
private:
    std::function<double(const PhasePoint&)> surface_;
    std::vector<PhasePoint> crossings_;

public:
    /**
     * @brief Define Poincaré section
     * @param surface Function defining surface: S(q,p) = 0
     */
    PoincareSection(std::function<double(const PhasePoint&)> surface)
        : surface_(surface) {}

    /**
     * @brief Check if trajectory crosses section
     */
    bool crosses(const PhasePoint& p1, const PhasePoint& p2) const {
        double s1 = surface_(p1);
        double s2 = surface_(p2);
        // Detect sign change (crossing)
        return (s1 * s2 < 0.0);
    }

    /**
     * @brief Record crossings from trajectory
     */
    void recordCrossings(const std::vector<PhasePoint>& trajectory) {
        crossings_.clear();

        for (size_t i = 1; i < trajectory.size(); ++i) {
            if (crosses(trajectory[i-1], trajectory[i])) {
                // Linear interpolation to find crossing point
                PhasePoint crossing = interpolateCrossing(
                    trajectory[i-1], trajectory[i]);
                crossings_.push_back(crossing);
            }
        }
    }

    /**
     * @brief Get recorded crossings
     */
    const std::vector<PhasePoint>& getCrossings() const {
        return crossings_;
    }

    /**
     * @brief Analyze periodicity from crossings
     */
    int detectPeriod(double tolerance = 1e-4) const {
        if (crossings_.size() < 2) return 0;

        const PhasePoint& first = crossings_[0];

        for (size_t i = 1; i < crossings_.size(); ++i) {
            double dist = (crossings_[i].q - first.q).norm() +
                         (crossings_[i].p - first.p).norm();
            if (dist < tolerance) {
                return static_cast<int>(i);
            }
        }

        return 0;  // No period detected
    }

private:
    PhasePoint interpolateCrossing(const PhasePoint& p1,
                                   const PhasePoint& p2) const {
        double s1 = surface_(p1);
        double s2 = surface_(p2);
        double alpha = -s1 / (s2 - s1);

        PhasePoint crossing = p1;
        crossing.q = p1.q + alpha * (p2.q - p1.q);
        crossing.p = p1.p + alpha * (p2.p - p1.p);
        return crossing;
    }
};

/**
 * @class LyapunovExponent
 * @brief Calculate Lyapunov exponents to detect chaos
 *
 * Measures exponential divergence of nearby trajectories:
 * λ = lim_{t→∞} (1/t) ln(δ(t)/δ(0))
 */
class LyapunovExponent {
private:
    const HamiltonianSystem& system_;

public:
    LyapunovExponent(const HamiltonianSystem& system) : system_(system) {}

    /**
     * @brief Compute largest Lyapunov exponent
     *
     * @param initial Initial phase point
     * @param dt Time step
     * @param num_steps Number of integration steps
     * @param delta0 Initial separation
     * @return Largest Lyapunov exponent λ
     */
    double compute(const PhasePoint& initial,
                  double dt, int num_steps,
                  double delta0 = 1e-8) const {

        // Create nearby trajectory
        PhasePoint perturbed = initial;
        perturbed.q(0) += delta0;

        auto traj1 = system_.integrate(initial, dt, num_steps);
        auto traj2 = system_.integrate(perturbed, dt, num_steps);

        // Calculate final separation
        const PhasePoint& final1 = traj1.back();
        const PhasePoint& final2 = traj2.back();

        double delta_q = (final2.q - final1.q).norm();
        double delta_p = (final2.p - final1.p).norm();
        double delta_final = std::sqrt(delta_q * delta_q + delta_p * delta_p);

        double time = dt * num_steps;
        return std::log(delta_final / delta0) / time;
    }

    /**
     * @brief Classify dynamics based on Lyapunov exponent
     */
    enum class DynamicsType {
        STABLE,       // λ < 0 (fixed point)
        NEUTRAL,      // λ ≈ 0 (quasiperiodic)
        CHAOTIC       // λ > 0 (chaos)
    };

    DynamicsType classify(double lambda, double tolerance = 1e-6) const {
        if (lambda < -tolerance) return DynamicsType::STABLE;
        if (lambda > tolerance) return DynamicsType::CHAOTIC;
        return DynamicsType::NEUTRAL;
    }
};

/**
 * @brief Symplectic matrix Ω for canonical structure
 *
 * For n degrees of freedom:
 *     ⎛  0   I ⎞
 * Ω = ⎜        ⎟
 *     ⎝ -I   0 ⎠
 *
 * where I is n×n identity matrix
 */
class SymplecticMatrix {
private:
    int dim_;
    Eigen::MatrixXd omega_;

public:
    SymplecticMatrix(int dimension) : dim_(dimension) {
        int n = 2 * dimension;
        omega_ = Eigen::MatrixXd::Zero(n, n);

        // Upper right block: +I
        omega_.block(0, dimension, dimension, dimension) =
            Eigen::MatrixXd::Identity(dimension, dimension);

        // Lower left block: -I
        omega_.block(dimension, 0, dimension, dimension) =
            -Eigen::MatrixXd::Identity(dimension, dimension);
    }

    const Eigen::MatrixXd& matrix() const { return omega_; }

    /**
     * @brief Verify symplectic property of transformation M
     *
     * M is symplectic if: M^T Ω M = Ω
     */
    bool isSymplectic(const Eigen::MatrixXd& M, double tolerance = 1e-10) const {
        Eigen::MatrixXd result = M.transpose() * omega_ * M;
        return (result - omega_).norm() < tolerance;
    }
};

} // namespace physics::advanced::classical

#endif // PHYSICS_ADVANCED_CLASSICAL_PHASE_SPACE_HPP
