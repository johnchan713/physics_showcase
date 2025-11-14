#ifndef PHYSICS_SPECIAL_RELATIVITY_HPP
#define PHYSICS_SPECIAL_RELATIVITY_HPP

#include <cmath>
#include <stdexcept>
#include <array>

/**
 * @file special_relativity.hpp
 * @brief Special relativity: Lorentz transformations and relativistic mechanics
 *
 * This module implements:
 * - Lorentz transformations
 * - Time dilation and length contraction
 * - Relativistic velocity addition
 * - Relativistic energy and momentum
 * - Mass-energy equivalence
 * - Relativistic Doppler effect
 * - Four-vectors and invariants
 *
 * All calculations use SI units with c = 299,792,458 m/s.
 */

namespace physics {
namespace special_relativity {

/**
 * @namespace constants
 * @brief Physical constants for special relativity
 */
namespace constants {
    constexpr double SPEED_OF_LIGHT = 299792458.0;  // Speed of light in vacuum (m/s)
}

// ============================================================================
// LORENTZ FACTOR
// ============================================================================

/**
 * @brief Calculate Lorentz factor γ
 *
 * γ = 1/√(1 - v²/c²)
 *
 * @param velocity Velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Lorentz factor γ (dimensionless, γ ≥ 1)
 * @throws std::invalid_argument if velocity >= c
 */
inline double lorentzFactor(double velocity, double speedOfLight = constants::SPEED_OF_LIGHT) {
    double beta = velocity / speedOfLight;
    double beta2 = beta * beta;

    if (beta2 >= 1.0) {
        throw std::invalid_argument("Velocity must be less than speed of light");
    }

    return 1.0 / std::sqrt(1.0 - beta2);
}

/**
 * @brief Calculate velocity from Lorentz factor
 *
 * v = c√(1 - 1/γ²)
 *
 * @param gamma Lorentz factor γ
 * @param speedOfLight Speed of light c (m/s)
 * @return Velocity (m/s)
 * @throws std::invalid_argument if γ < 1
 */
inline double velocityFromLorentzFactor(double gamma, double speedOfLight = constants::SPEED_OF_LIGHT) {
    if (gamma < 1.0) {
        throw std::invalid_argument("Lorentz factor must be >= 1");
    }
    return speedOfLight * std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

/**
 * @brief Calculate β = v/c
 *
 * @param velocity Velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return β (dimensionless, 0 ≤ β < 1)
 * @throws std::invalid_argument if velocity >= c
 */
inline double beta(double velocity, double speedOfLight = constants::SPEED_OF_LIGHT) {
    double b = velocity / speedOfLight;
    if (std::abs(b) >= 1.0) {
        throw std::invalid_argument("Velocity magnitude must be less than speed of light");
    }
    return b;
}

// ============================================================================
// TIME DILATION
// ============================================================================

/**
 * @brief Calculate time dilation
 *
 * Δt = γΔt₀
 *
 * Time interval in moving frame appears longer in stationary frame
 *
 * @param properTime Proper time Δt₀ (time in rest frame) (s)
 * @param velocity Relative velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Dilated time Δt (s)
 */
inline double timeDilation(double properTime, double velocity,
                          double speedOfLight = constants::SPEED_OF_LIGHT) {
    double gamma = lorentzFactor(velocity, speedOfLight);
    return gamma * properTime;
}

/**
 * @brief Calculate proper time from dilated time
 *
 * Δt₀ = Δt/γ = Δt√(1 - v²/c²)
 *
 * @param dilatedTime Time in stationary frame Δt (s)
 * @param velocity Relative velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Proper time Δt₀ (s)
 */
inline double properTime(double dilatedTime, double velocity,
                        double speedOfLight = constants::SPEED_OF_LIGHT) {
    double gamma = lorentzFactor(velocity, speedOfLight);
    return dilatedTime / gamma;
}

// ============================================================================
// LENGTH CONTRACTION
// ============================================================================

/**
 * @brief Calculate length contraction
 *
 * L = L₀/γ = L₀√(1 - v²/c²)
 *
 * Length in direction of motion appears contracted
 *
 * @param properLength Rest length L₀ (m)
 * @param velocity Relative velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Contracted length L (m)
 */
inline double lengthContraction(double properLength, double velocity,
                                double speedOfLight = constants::SPEED_OF_LIGHT) {
    double gamma = lorentzFactor(velocity, speedOfLight);
    return properLength / gamma;
}

/**
 * @brief Calculate proper length from contracted length
 *
 * L₀ = γL
 *
 * @param contractedLength Length in moving frame L (m)
 * @param velocity Relative velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Proper length L₀ (m)
 */
inline double properLength(double contractedLength, double velocity,
                          double speedOfLight = constants::SPEED_OF_LIGHT) {
    double gamma = lorentzFactor(velocity, speedOfLight);
    return gamma * contractedLength;
}

// ============================================================================
// RELATIVISTIC VELOCITY ADDITION
// ============================================================================

/**
 * @brief Calculate relativistic velocity addition (1D)
 *
 * u = (v + u')/[1 + (vu')/(c²)]
 *
 * @param velocityFrame Velocity of frame v (m/s)
 * @param velocityInFrame Velocity in moving frame u' (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Velocity in stationary frame u (m/s)
 */
inline double velocityAddition(double velocityFrame, double velocityInFrame,
                               double speedOfLight = constants::SPEED_OF_LIGHT) {
    double c2 = speedOfLight * speedOfLight;
    double numerator = velocityFrame + velocityInFrame;
    double denominator = 1.0 + (velocityFrame * velocityInFrame) / c2;

    if (std::abs(denominator) < 1e-15) {
        throw std::invalid_argument("Invalid velocity configuration");
    }

    return numerator / denominator;
}

/**
 * @brief Verify that velocity addition never exceeds c
 *
 * For any v, u' < c: (v + u')/(1 + vu'/c²) < c
 *
 * @param result Result of velocity addition (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return true if result < c
 */
inline bool isValidRelativistic Velocity(double result, double speedOfLight = constants::SPEED_OF_LIGHT) {
    return std::abs(result) < speedOfLight;
}

// ============================================================================
// RELATIVISTIC ENERGY
// ============================================================================

/**
 * @brief Calculate relativistic total energy
 *
 * E = γmc²
 *
 * @param mass Rest mass m (kg)
 * @param velocity Velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Total energy E (J)
 */
inline double relativisticEnergy(double mass, double velocity,
                                 double speedOfLight = constants::SPEED_OF_LIGHT) {
    if (mass < 0.0) {
        throw std::invalid_argument("Mass cannot be negative");
    }
    double gamma = lorentzFactor(velocity, speedOfLight);
    return gamma * mass * speedOfLight * speedOfLight;
}

/**
 * @brief Calculate rest mass energy
 *
 * E₀ = mc²
 *
 * Einstein's mass-energy equivalence
 *
 * @param mass Rest mass m (kg)
 * @param speedOfLight Speed of light c (m/s)
 * @return Rest energy E₀ (J)
 */
inline double restEnergy(double mass, double speedOfLight = constants::SPEED_OF_LIGHT) {
    if (mass < 0.0) {
        throw std::invalid_argument("Mass cannot be negative");
    }
    return mass * speedOfLight * speedOfLight;
}

/**
 * @brief Calculate relativistic kinetic energy
 *
 * K = (γ - 1)mc² = E - mc²
 *
 * @param mass Rest mass m (kg)
 * @param velocity Velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Kinetic energy K (J)
 */
inline double relativisticKineticEnergy(double mass, double velocity,
                                        double speedOfLight = constants::SPEED_OF_LIGHT) {
    double E_total = relativisticEnergy(mass, velocity, speedOfLight);
    double E_rest = restEnergy(mass, speedOfLight);
    return E_total - E_rest;
}

/**
 * @brief Calculate velocity from kinetic energy
 *
 * v = c√[1 - (mc²/(K + mc²))²]
 *
 * @param kineticEnergy Kinetic energy K (J)
 * @param mass Rest mass m (kg)
 * @param speedOfLight Speed of light c (m/s)
 * @return Velocity (m/s)
 */
inline double velocityFromKineticEnergy(double kineticEnergy, double mass,
                                        double speedOfLight = constants::SPEED_OF_LIGHT) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    double E_rest = restEnergy(mass, speedOfLight);
    double E_total = kineticEnergy + E_rest;
    double ratio = E_rest / E_total;
    return speedOfLight * std::sqrt(1.0 - ratio * ratio);
}

// ============================================================================
// RELATIVISTIC MOMENTUM
// ============================================================================

/**
 * @brief Calculate relativistic momentum
 *
 * p = γmv
 *
 * @param mass Rest mass m (kg)
 * @param velocity Velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Momentum p (kg⋅m/s)
 */
inline double relativisticMomentum(double mass, double velocity,
                                   double speedOfLight = constants::SPEED_OF_LIGHT) {
    if (mass < 0.0) {
        throw std::invalid_argument("Mass cannot be negative");
    }
    double gamma = lorentzFactor(velocity, speedOfLight);
    return gamma * mass * velocity;
}

/**
 * @brief Calculate energy-momentum relation
 *
 * E² = (pc)² + (mc²)²
 *
 * Relativistic invariant
 *
 * @param momentum Momentum p (kg⋅m/s)
 * @param mass Rest mass m (kg)
 * @param speedOfLight Speed of light c (m/s)
 * @return Total energy E (J)
 */
inline double energyFromMomentum(double momentum, double mass,
                                 double speedOfLight = constants::SPEED_OF_LIGHT) {
    double pc = momentum * speedOfLight;
    double mc2 = mass * speedOfLight * speedOfLight;
    return std::sqrt(pc * pc + mc2 * mc2);
}

/**
 * @brief Calculate momentum from energy
 *
 * p = √(E² - (mc²)²)/c
 *
 * @param energy Total energy E (J)
 * @param mass Rest mass m (kg)
 * @param speedOfLight Speed of light c (m/s)
 * @return Momentum p (kg⋅m/s)
 */
inline double momentumFromEnergy(double energy, double mass,
                                 double speedOfLight = constants::SPEED_OF_LIGHT) {
    double mc2 = mass * speedOfLight * speedOfLight;
    double E2_minus_m2c4 = energy * energy - mc2 * mc2;

    if (E2_minus_m2c4 < 0.0) {
        throw std::invalid_argument("Energy less than rest energy");
    }

    return std::sqrt(E2_minus_m2c4) / speedOfLight;
}

/**
 * @brief Calculate velocity from momentum
 *
 * v = pc²/E = pc²/√[(pc)² + (mc²)²]
 *
 * @param momentum Momentum p (kg⋅m/s)
 * @param mass Rest mass m (kg)
 * @param speedOfLight Speed of light c (m/s)
 * @return Velocity (m/s)
 */
inline double velocityFromMomentum(double momentum, double mass,
                                   double speedOfLight = constants::SPEED_OF_LIGHT) {
    double energy = energyFromMomentum(momentum, mass, speedOfLight);
    return momentum * speedOfLight * speedOfLight / energy;
}

// ============================================================================
// MASSLESS PARTICLES (PHOTONS)
// ============================================================================

/**
 * @brief Calculate photon energy
 *
 * E = pc (for m = 0)
 *
 * @param momentum Photon momentum p (kg⋅m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Photon energy (J)
 */
inline double photonEnergy(double momentum, double speedOfLight = constants::SPEED_OF_LIGHT) {
    return momentum * speedOfLight;
}

/**
 * @brief Calculate photon momentum from energy
 *
 * p = E/c
 *
 * @param energy Photon energy E (J)
 * @param speedOfLight Speed of light c (m/s)
 * @return Photon momentum (kg⋅m/s)
 */
inline double photonMomentum(double energy, double speedOfLight = constants::SPEED_OF_LIGHT) {
    return energy / speedOfLight;
}

// ============================================================================
// RELATIVISTIC DOPPLER EFFECT
// ============================================================================

/**
 * @brief Calculate relativistic Doppler shift (longitudinal)
 *
 * f' = f√[(1 - β)/(1 + β)]  (receding)
 * f' = f√[(1 + β)/(1 - β)]  (approaching)
 *
 * @param frequency Source frequency f (Hz)
 * @param velocity Relative velocity v (positive = approaching) (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Observed frequency f' (Hz)
 */
inline double relativisticDoppler(double frequency, double velocity,
                                  double speedOfLight = constants::SPEED_OF_LIGHT) {
    double b = beta(velocity, speedOfLight);

    if (velocity > 0.0) {
        // Approaching: f' = f√[(1 + β)/(1 - β)]
        return frequency * std::sqrt((1.0 + b) / (1.0 - b));
    } else {
        // Receding: f' = f√[(1 - β)/(1 + β)]
        double bAbs = std::abs(b);
        return frequency * std::sqrt((1.0 - bAbs) / (1.0 + bAbs));
    }
}

/**
 * @brief Calculate transverse Doppler effect
 *
 * f' = f/γ
 *
 * Pure time dilation effect when motion is perpendicular
 *
 * @param frequency Source frequency f (Hz)
 * @param velocity Transverse velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Observed frequency f' (Hz)
 */
inline double transverseDoppler(double frequency, double velocity,
                                double speedOfLight = constants::SPEED_OF_LIGHT) {
    double gamma = lorentzFactor(velocity, speedOfLight);
    return frequency / gamma;
}

/**
 * @brief Calculate redshift parameter z
 *
 * z = Δλ/λ₀ = (λ_obs - λ_em)/λ_em
 *
 * For recession: z = √[(1 + β)/(1 - β)] - 1
 *
 * @param velocity Recession velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Redshift z
 */
inline double cosmologicalRedshift(double velocity, double speedOfLight = constants::SPEED_OF_LIGHT) {
    double b = beta(velocity, speedOfLight);
    return std::sqrt((1.0 + b) / (1.0 - b)) - 1.0;
}

// ============================================================================
// FOUR-VECTORS AND INVARIANTS
// ============================================================================

/**
 * @brief Calculate spacetime interval (Minkowski metric)
 *
 * s² = c²t² - x² - y² - z² = c²t² - r²
 *
 * Lorentz invariant
 *
 * @param time Time coordinate t (s)
 * @param spatialDistance Spatial distance r (m)
 * @param speedOfLight Speed of light c (m/s)
 * @return Spacetime interval squared s² (m²)
 */
inline double spacetimeInterval(double time, double spatialDistance,
                                double speedOfLight = constants::SPEED_OF_LIGHT) {
    double ct = speedOfLight * time;
    return ct * ct - spatialDistance * spatialDistance;
}

/**
 * @brief Calculate proper time along worldline
 *
 * τ = √(s²)/c
 *
 * @param spacetimeInterval2 Spacetime interval squared s² (m²)
 * @param speedOfLight Speed of light c (m/s)
 * @return Proper time τ (s)
 * @throws std::invalid_argument if s² < 0 (spacelike interval)
 */
inline double properTimeFromInterval(double spacetimeInterval2,
                                     double speedOfLight = constants::SPEED_OF_LIGHT) {
    if (spacetimeInterval2 < 0.0) {
        throw std::invalid_argument("Spacelike interval has no proper time");
    }
    return std::sqrt(spacetimeInterval2) / speedOfLight;
}

/**
 * @brief Calculate four-momentum magnitude (invariant mass)
 *
 * p_μp^μ = E²/c² - p² = m²c²
 *
 * @param energy Total energy E (J)
 * @param momentum Three-momentum magnitude p (kg⋅m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Invariant mass m (kg)
 */
inline double invariantMass(double energy, double momentum,
                            double speedOfLight = constants::SPEED_OF_LIGHT) {
    double E_over_c = energy / speedOfLight;
    double p_mu_squared = E_over_c * E_over_c - momentum * momentum;

    if (p_mu_squared < 0.0) {
        throw std::invalid_argument("Invalid four-momentum");
    }

    return std::sqrt(p_mu_squared) / speedOfLight;
}

/**
 * @brief Calculate rapidity
 *
 * η = arctanh(v/c) = (1/2)ln[(1 + β)/(1 - β)]
 *
 * Rapidity adds linearly under boosts
 *
 * @param velocity Velocity v (m/s)
 * @param speedOfLight Speed of light c (m/s)
 * @return Rapidity η (dimensionless)
 */
inline double rapidity(double velocity, double speedOfLight = constants::SPEED_OF_LIGHT) {
    double b = beta(velocity, speedOfLight);
    return 0.5 * std::log((1.0 + b) / (1.0 - b));
}

/**
 * @brief Calculate velocity from rapidity
 *
 * v = c tanh(η)
 *
 * @param eta Rapidity η
 * @param speedOfLight Speed of light c (m/s)
 * @return Velocity (m/s)
 */
inline double velocityFromRapidity(double eta, double speedOfLight = constants::SPEED_OF_LIGHT) {
    return speedOfLight * std::tanh(eta);
}

} // namespace special_relativity
} // namespace physics

#endif // PHYSICS_SPECIAL_RELATIVITY_HPP
