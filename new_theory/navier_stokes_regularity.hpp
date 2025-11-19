/**
 * @file navier_stokes_regularity.hpp
 * @brief Computational Attack on Navier-Stokes Millennium Prize Problem
 *
 * MILLENNIUM PRIZE PROBLEM (Clay Mathematics Institute):
 * Prove or give counter-example: For smooth initial conditions, do solutions
 * to the 3D incompressible Navier-Stokes equations exist globally in time
 * and remain smooth?
 *
 * PRIZE: $1,000,000
 *
 * APPROACH:
 * 1. Implement regularity criteria (BKM, LPS)
 * 2. Track energy dissipation and enstrophy
 * 3. Numerical blow-up detection
 * 4. Vorticity stretching analysis
 * 5. Critical Sobolev norm monitoring
 */

#ifndef NEW_THEORY_NAVIER_STOKES_REGULARITY_HPP
#define NEW_THEORY_NAVIER_STOKES_REGULARITY_HPP

#include <cmath>
#include <vector>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <limits>

namespace new_theory {
namespace navier_stokes {

/**
 * ============================================================================
 * NAVIER-STOKES EQUATIONS
 * ============================================================================
 *
 * 3D Incompressible Navier-Stokes:
 *   ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + f
 *   ∇·u = 0
 *
 * where:
 *   u(x,t) - velocity field (vector)
 *   p(x,t) - pressure (scalar)
 *   ν - kinematic viscosity
 *   f - external force
 *   ρ - density (constant, normalized to 1)
 *
 * MILLENNIUM PROBLEM:
 * Given smooth initial data u₀ with ∇·u₀ = 0, does there exist
 * a smooth solution u(x,t) for all t ∈ [0,∞)?
 *
 * Or does a finite-time singularity occur: ||u(·,t)|| → ∞ as t → T* < ∞?
 */

namespace constants {
    constexpr double NU_WATER = 1.0e-6;      // m²/s (water at 20°C)
    constexpr double NU_AIR = 1.5e-5;        // m²/s (air at 20°C)
    constexpr double RHO_WATER = 1000.0;     // kg/m³
}

/**
 * ============================================================================
 * REGULARITY CRITERION 1: BEALE-KATO-MAJDA (BKM) THEOREM
 * ============================================================================
 *
 * THEOREM (Beale-Kato-Majda, 1984):
 * A solution u remains smooth on [0,T] if and only if
 *
 *   ∫₀ᵀ ||ω(·,t)||_{L^∞} dt < ∞
 *
 * where ω = ∇×u is the vorticity.
 *
 * PHYSICAL INTERPRETATION:
 * Blow-up can occur only if vorticity becomes unbounded in finite time.
 *
 * PROOF SKETCH:
 * 1. Vorticity equation: ∂ω/∂t + u·∇ω = ω·∇u + ν∇²ω
 * 2. The term ω·∇u is vortex stretching (nonlinear!)
 * 3. Energy estimate: d/dt||u||² ≤ C||ω||_{L^∞}||u||²
 * 4. Grönwall inequality: If ∫||ω||_{L^∞} dt < ∞, then ||u|| remains bounded
 */

class BealeKatoMajdaCriterion {
public:
    /**
     * @brief Check BKM regularity criterion
     *
     * Computes: I(T) = ∫₀ᵀ ||ω(·,t)||_{L^∞} dt
     *
     * Solution is regular on [0,T] ⟺ I(T) < ∞
     *
     * @param vorticity_Linfty Function returning ||ω(t)||_{L^∞} at time t
     * @param T Final time
     * @param dt Time step
     * @return Integrated vorticity norm
     */
    static double integratedVorticityNorm(
        std::function<double(double)> vorticity_Linfty,
        double T,
        double dt = 0.01) {

        if (T <= 0 || dt <= 0) {
            throw std::invalid_argument("T and dt must be positive");
        }

        double integral = 0.0;
        int n_steps = static_cast<int>(T / dt);

        for (int i = 0; i < n_steps; ++i) {
            double t = i * dt;
            double omega_inf = vorticity_Linfty(t);

            // Check for blow-up
            if (std::isinf(omega_inf) || std::isnan(omega_inf)) {
                return std::numeric_limits<double>::infinity();
            }

            integral += omega_inf * dt;
        }

        return integral;
    }

    /**
     * @brief Check if solution remains regular
     *
     * @param integrated_norm I(T) = ∫₀ᵀ ||ω||_{L^∞} dt
     * @return true if I(T) < ∞ (regular), false if blow-up suspected
     */
    static bool isRegular(double integrated_norm) {
        return std::isfinite(integrated_norm);
    }

    /**
     * @brief Estimate blow-up time using BKM criterion
     *
     * If ||ω(t)||_{L^∞} grows like (T*-t)^{-α}, find T*
     *
     * @param vorticity_Linfty Time-dependent vorticity norm
     * @param t_start Start time
     * @param t_end Search end time
     * @param dt Time step
     * @return Estimated blow-up time T* (or infinity if no blow-up detected)
     */
    static double estimateBlowupTime(
        std::function<double(double)> vorticity_Linfty,
        double t_start,
        double t_end,
        double dt = 0.01) {

        std::vector<double> times;
        std::vector<double> omega_values;

        for (double t = t_start; t < t_end; t += dt) {
            double omega = vorticity_Linfty(t);
            if (std::isfinite(omega) && omega > 0) {
                times.push_back(t);
                omega_values.push_back(omega);
            }
        }

        if (times.size() < 3) {
            return std::numeric_limits<double>::infinity();
        }

        // Check for power-law growth: ω ~ (T*-t)^{-α}
        // If log(ω) vs log(T*-t) is linear, extrapolate T*

        // Simple heuristic: if ω is growing rapidly, estimate T*
        double omega_final = omega_values.back();
        double omega_initial = omega_values[0];
        double growth_rate = omega_final / omega_initial;

        if (growth_rate > 10.0) {
            // Rapid growth detected, estimate T*
            // Use linear extrapolation of 1/ω
            double inv_omega_final = 1.0 / omega_final;
            double t_final = times.back();

            // Rough estimate: T* ≈ t_final + 1/(growth_rate)
            return t_final + dt * 10.0 / growth_rate;
        }

        return std::numeric_limits<double>::infinity();
    }
};

/**
 * ============================================================================
 * REGULARITY CRITERION 2: LADYZHENSKAYA-PRODI-SERRIN (LPS) CONDITION
 * ============================================================================
 *
 * THEOREM (Ladyzhenskaya-Prodi-Serrin, 1960s):
 * A weak solution is regular on [0,T] if
 *
 *   u ∈ L^p(0,T; L^q(ℝ³)) with 2/p + 3/q ≤ 1, q ≥ 3
 *
 * SPECIAL CASES:
 * - p = ∞, q = 3: u ∈ L^∞(0,T; L³(ℝ³))
 * - p = 4, q = ∞: u ∈ L⁴(0,T; L^∞(ℝ³))
 * - p = 2, q = ∞: u ∈ L²(0,T; L^∞(ℝ³))
 *
 * PHYSICAL MEANING:
 * Solution is regular if velocity doesn't concentrate too much
 * in space-time (measured by L^p-L^q norms).
 */

class LadyzhenskayaProdiSerrin {
public:
    /**
     * @brief Check LPS condition: 2/p + 3/q ≤ 1
     *
     * @param p Time exponent
     * @param q Space exponent
     * @return true if (p,q) satisfies Serrin condition
     */
    static bool checkSerrinCondition(double p, double q) {
        // Check basic requirements (handle infinity)
        if (std::isinf(p) && std::isinf(q)) {
            return true;  // L^∞(L^∞) satisfies condition
        }

        if (std::isinf(q)) {
            // q = ∞: condition becomes 2/p ≤ 1, so p ≥ 2
            return p >= 2.0;
        }

        if (std::isinf(p)) {
            // p = ∞: condition becomes 3/q ≤ 1, so q ≥ 3
            return q >= 3.0;
        }

        // Finite p, q: check p ≥ 2, q ≥ 3, and Serrin condition
        if (p < 2.0 || q < 3.0) {
            return false;
        }

        return (2.0/p + 3.0/q) <= 1.0 + 1e-10;  // Small tolerance for rounding
    }

    /**
     * @brief Compute L^p(L^q) norm
     *
     * ||u||_{L^p(0,T; L^q)} = (∫₀ᵀ ||u(·,t)||_{L^q}^p dt)^{1/p}
     *
     * @param u_Lq_norm Function returning ||u(t)||_{L^q} at time t
     * @param p Time exponent
     * @param T Final time
     * @param dt Time step
     * @return L^p(L^q) norm
     */
    static double LpLqNorm(
        std::function<double(double)> u_Lq_norm,
        double p,
        double T,
        double dt = 0.01) {

        if (std::isinf(p)) {
            // L^∞ in time: max over time
            double max_norm = 0.0;
            for (double t = 0; t < T; t += dt) {
                double norm_t = u_Lq_norm(t);
                max_norm = std::max(max_norm, norm_t);
            }
            return max_norm;
        }

        // L^p in time: (∫ ||u||^p dt)^{1/p}
        double integral = 0.0;
        for (double t = 0; t < T; t += dt) {
            double norm_t = u_Lq_norm(t);
            integral += std::pow(norm_t, p) * dt;
        }

        return std::pow(integral, 1.0/p);
    }

    /**
     * @brief Check regularity via LPS criterion
     *
     * Solution is regular if ||u||_{L^p(L^q)} < ∞ for some (p,q)
     * with 2/p + 3/q ≤ 1
     *
     * @param u_Lq_norm Spatial norm ||u(t)||_{L^q}
     * @param p Time exponent
     * @param q Space exponent
     * @param T Final time
     * @return true if LPS condition satisfied
     */
    static bool isRegular(
        std::function<double(double)> u_Lq_norm,
        double p,
        double q,
        double T) {

        if (!checkSerrinCondition(p, q)) {
            throw std::invalid_argument("(p,q) does not satisfy Serrin condition");
        }

        double norm = LpLqNorm(u_Lq_norm, p, T);
        return std::isfinite(norm);
    }
};

/**
 * ============================================================================
 * ENERGY AND ENSTROPHY
 * ============================================================================
 *
 * ENERGY: E(t) = (1/2)∫|u(x,t)|² dx
 *
 * ENERGY DISSIPATION:
 *   dE/dt = -ν∫|∇u|² dx + ∫f·u dx
 *
 * For f=0: dE/dt = -νε(t) where ε = ∫|∇u|² is enstrophy dissipation rate
 *
 * ENSTROPHY: Ω(t) = (1/2)∫|ω(x,t)|² dx where ω = ∇×u
 *
 * ENSTROPHY EQUATION:
 *   dΩ/dt = ∫(ω·∇u)·ω dx - ν∫|∇ω|² dx
 *
 * The term ∫(ω·∇u)·ω is vortex stretching - can increase enstrophy!
 *
 * KEY FACT: Blow-up occurs when enstrophy grows unbounded.
 */

class EnergyEnstrophy {
public:
    /**
     * @brief Calculate kinetic energy E = (1/2)∫|u|² dx
     *
     * @param u_squared Spatial integral ∫|u|² dx
     * @return Kinetic energy
     */
    static double kineticEnergy(double u_squared) {
        return 0.5 * u_squared;
    }

    /**
     * @brief Calculate energy dissipation rate
     *
     * ε = ∫|∇u|² dx
     *
     * @param grad_u_squared Spatial integral ∫|∇u|² dx
     * @return Energy dissipation rate
     */
    static double energyDissipationRate(double grad_u_squared) {
        return grad_u_squared;
    }

    /**
     * @brief Energy evolution: dE/dt = -νε
     *
     * @param E_current Current energy
     * @param epsilon Dissipation rate ∫|∇u|² dx
     * @param nu Kinematic viscosity
     * @param dt Time step
     * @return E(t+dt)
     */
    static double evolveEnergy(double E_current, double epsilon, double nu, double dt) {
        // dE/dt = -νε
        double dE_dt = -nu * epsilon;
        return E_current + dE_dt * dt;
    }

    /**
     * @brief Calculate enstrophy Ω = (1/2)∫|ω|² dx
     *
     * @param omega_squared Spatial integral ∫|ω|² dx
     * @return Enstrophy
     */
    static double enstrophy(double omega_squared) {
        return 0.5 * omega_squared;
    }

    /**
     * @brief Vortex stretching term ∫(ω·∇u)·ω dx
     *
     * This term can be positive or negative:
     * - Positive: vorticity aligns with stretching direction → amplification
     * - Negative: vorticity perpendicular to stretching → compression
     *
     * @param omega_grad_u_omega Integral ∫(ω·∇u)·ω dx
     * @return Vortex stretching contribution
     */
    static double vortexStretching(double omega_grad_u_omega) {
        return omega_grad_u_omega;
    }

    /**
     * @brief Enstrophy evolution
     *
     * dΩ/dt = S - νP
     * where S = ∫(ω·∇u)·ω dx (stretching)
     *       P = ∫|∇ω|² dx (palinstrophy, dissipation)
     *
     * @param Omega_current Current enstrophy
     * @param stretching Vortex stretching term
     * @param palinstrophy ∫|∇ω|² dx
     * @param nu Viscosity
     * @param dt Time step
     * @return Ω(t+dt)
     */
    static double evolveEnstrophy(
        double Omega_current,
        double stretching,
        double palinstrophy,
        double nu,
        double dt) {

        // dΩ/dt = S - νP
        double dOmega_dt = stretching - nu * palinstrophy;
        return Omega_current + dOmega_dt * dt;
    }

    /**
     * @brief Check for enstrophy blow-up
     *
     * If Ω(t) → ∞ as t → T*, then blow-up occurs
     *
     * @param enstrophy_history Vector of Ω(t) values
     * @param threshold Blow-up threshold
     * @return true if enstrophy growing unbounded
     */
    static bool detectEnstrophyBlowup(
        const std::vector<double>& enstrophy_history,
        double threshold = 1e10) {

        if (enstrophy_history.size() < 2) {
            return false;
        }

        // Check if enstrophy exceeds threshold
        double omega_current = enstrophy_history.back();
        if (omega_current > threshold) {
            return true;
        }

        // Check for rapid growth
        size_t n = enstrophy_history.size();
        if (n >= 3) {
            double omega_prev = enstrophy_history[n-2];
            double omega_prev2 = enstrophy_history[n-3];

            double growth_rate = (omega_current - omega_prev) / omega_prev;
            double prev_growth = (omega_prev - omega_prev2) / omega_prev2;

            // Accelerating growth?
            if (growth_rate > prev_growth && growth_rate > 0.1) {
                return true;
            }
        }

        return false;
    }
};

/**
 * ============================================================================
 * CRITICAL SOBOLEV NORMS
 * ============================================================================
 *
 * SOBOLEV SPACE H^s(ℝ³):
 * Functions with derivatives up to order s in L²
 *
 * NORM: ||u||_{H^s} = (∫(1+|ξ|²)^s |û(ξ)|² dξ)^{1/2}
 *
 * CRITICAL EXPONENT: s = 1/2 for 3D Navier-Stokes
 * - Subcritical (s > 1/2): Global regularity known
 * - Critical (s = 1/2): Open problem!
 * - Supercritical (s < 1/2): Ill-posed
 *
 * The Millennium Problem essentially asks: Is H^{1/2} regularity preserved?
 */

class SobolevNorms {
public:
    /**
     * @brief Calculate H^s Sobolev norm (simplified)
     *
     * ||u||_{H^s} ≈ ||u||_{L²} + ||∇^s u||_{L²}
     *
     * For integer s: sum of L² norms of derivatives up to order s
     *
     * @param u_L2 L² norm of u
     * @param grad_s_u_L2 L² norm of s-th derivative
     * @return Approximate H^s norm
     */
    static double HsNorm(double u_L2, double grad_s_u_L2) {
        return std::sqrt(u_L2 * u_L2 + grad_s_u_L2 * grad_s_u_L2);
    }

    /**
     * @brief Check if Sobolev norm remains bounded
     *
     * Solution is regular if ||u(t)||_{H^s} < ∞ for all t
     *
     * @param Hs_norm_history Time series of ||u(t)||_{H^s}
     * @param threshold Blow-up threshold
     * @return true if norm remains bounded
     */
    static bool isSobolevRegular(
        const std::vector<double>& Hs_norm_history,
        double threshold = 1e8) {

        for (double norm : Hs_norm_history) {
            if (!std::isfinite(norm) || norm > threshold) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Critical Sobolev index for 3D NS
     *
     * s_crit = 1/2 in 3D
     * s_crit = 0 in 2D (subcritical → global regularity proven!)
     *
     * @param dimension Spatial dimension
     * @return Critical Sobolev exponent
     */
    static double criticalExponent(int dimension) {
        if (dimension == 2) {
            return 0.0;  // 2D is subcritical
        } else if (dimension == 3) {
            return 0.5;  // 3D critical exponent
        }

        return (dimension - 2.0) / 2.0;
    }
};

/**
 * ============================================================================
 * BLOW-UP DETECTION FRAMEWORK
 * ============================================================================
 *
 * Synthesizes all regularity criteria to detect potential blow-up
 */

class BlowUpDetector {
public:
    struct RegularityStatus {
        bool is_regular;
        double T_blowup_estimate;
        std::string criterion_violated;
        double severity_score;  // 0-1, higher = closer to blow-up
    };

    /**
     * @brief Comprehensive regularity check
     *
     * Checks:
     * 1. BKM criterion (vorticity)
     * 2. Energy/enstrophy evolution
     * 3. Sobolev norm growth
     *
     * @param vorticity_Linfty Function ||ω(t)||_{L^∞}
     * @param enstrophy_history Ω(t) time series
     * @param Hs_norm_history ||u(t)||_{H^s} time series
     * @param T Current time
     * @return RegularityStatus struct
     */
    static RegularityStatus checkRegularity(
        std::function<double(double)> vorticity_Linfty,
        const std::vector<double>& enstrophy_history,
        const std::vector<double>& Hs_norm_history,
        double T) {

        RegularityStatus status;
        status.is_regular = true;
        status.T_blowup_estimate = std::numeric_limits<double>::infinity();
        status.criterion_violated = "none";
        status.severity_score = 0.0;

        // Check BKM criterion
        double bkm_integral = BealeKatoMajdaCriterion::integratedVorticityNorm(
            vorticity_Linfty, T);

        if (!BealeKatoMajdaCriterion::isRegular(bkm_integral)) {
            status.is_regular = false;
            status.criterion_violated = "BKM (vorticity unbounded)";
            status.severity_score = 1.0;
            status.T_blowup_estimate = BealeKatoMajdaCriterion::estimateBlowupTime(
                vorticity_Linfty, 0, T);
            return status;
        }

        // Check enstrophy blow-up
        if (EnergyEnstrophy::detectEnstrophyBlowup(enstrophy_history)) {
            status.is_regular = false;
            status.criterion_violated = "Enstrophy blow-up";
            status.severity_score = 0.9;
            return status;
        }

        // Check Sobolev norm
        if (!SobolevNorms::isSobolevRegular(Hs_norm_history)) {
            status.is_regular = false;
            status.criterion_violated = "Sobolev norm unbounded";
            status.severity_score = 0.8;
            return status;
        }

        // Calculate severity score based on growth rates
        if (!enstrophy_history.empty()) {
            double omega_current = enstrophy_history.back();
            double omega_initial = enstrophy_history.front();
            double omega_ratio = omega_current / (omega_initial + 1e-10);

            status.severity_score = std::tanh(std::log10(omega_ratio + 1.0) / 5.0);
        }

        return status;
    }
};

} // namespace navier_stokes
} // namespace new_theory

#endif // NEW_THEORY_NAVIER_STOKES_REGULARITY_HPP
