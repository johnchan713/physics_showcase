#ifndef PHYSICS_ENERGY_MOMENTUM_HPP
#define PHYSICS_ENERGY_MOMENTUM_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace energy_momentum {

/**
 * @brief Energy and Momentum - Comparison and Relationships
 *
 * This namespace explores the relationship between kinetic energy (KE)
 * and momentum (p), two fundamental quantities in mechanics.
 *
 * Key formulas:
 * - Kinetic Energy: KE = (1/2)mv²
 * - Momentum: p = mv
 * - Relationship: KE = p²/(2m)
 * - Also: p = √(2m·KE)
 */

// ============================================================================
// Kinetic Energy Calculations
// ============================================================================

/**
 * @brief Calculate kinetic energy from mass and velocity
 *
 * KE = (1/2)mv²
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param velocity Velocity of object (in m/s)
 * @return Kinetic energy (in Joules)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateKineticEnergy(double mass, double velocity) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return 0.5 * mass * velocity * velocity;
}

/**
 * @brief Calculate velocity from kinetic energy and mass
 *
 * From KE = (1/2)mv², we get v = √(2·KE/m)
 *
 * @param kineticEnergy Kinetic energy (in Joules, must be >= 0)
 * @param mass Mass of object (in kilograms, must be > 0)
 * @return Velocity magnitude (in m/s)
 * @throws std::invalid_argument if KE < 0 or mass <= 0
 */
inline double calculateVelocityFromKE(double kineticEnergy, double mass) {
    if (kineticEnergy < 0) {
        throw std::invalid_argument("Kinetic energy cannot be negative");
    }
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return std::sqrt(2.0 * kineticEnergy / mass);
}

// ============================================================================
// Momentum Calculations
// ============================================================================

/**
 * @brief Calculate momentum from mass and velocity
 *
 * p = mv
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param velocity Velocity of object (in m/s)
 * @return Momentum (in kg·m/s)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateMomentum(double mass, double velocity) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return mass * velocity;
}

/**
 * @brief Calculate velocity from momentum and mass
 *
 * From p = mv, we get v = p/m
 *
 * @param momentum Momentum (in kg·m/s)
 * @param mass Mass of object (in kilograms, must be > 0)
 * @return Velocity (in m/s)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateVelocityFromMomentum(double momentum, double mass) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return momentum / mass;
}

// ============================================================================
// Energy-Momentum Relationships
// ============================================================================

/**
 * @brief Calculate kinetic energy from momentum
 *
 * Relationship: KE = p²/(2m)
 *
 * This shows that for equal momentum, lighter objects have more kinetic energy.
 *
 * @param momentum Momentum (in kg·m/s)
 * @param mass Mass of object (in kilograms, must be > 0)
 * @return Kinetic energy (in Joules)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateKEFromMomentum(double momentum, double mass) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return (momentum * momentum) / (2.0 * mass);
}

/**
 * @brief Calculate momentum from kinetic energy
 *
 * Relationship: p = √(2m·KE)
 *
 * This shows that for equal kinetic energy, heavier objects have more momentum.
 *
 * @param kineticEnergy Kinetic energy (in Joules, must be >= 0)
 * @param mass Mass of object (in kilograms, must be > 0)
 * @return Momentum magnitude (in kg·m/s)
 * @throws std::invalid_argument if KE < 0 or mass <= 0
 */
inline double calculateMomentumFromKE(double kineticEnergy, double mass) {
    if (kineticEnergy < 0) {
        throw std::invalid_argument("Kinetic energy cannot be negative");
    }
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return std::sqrt(2.0 * mass * kineticEnergy);
}

// ============================================================================
// Comparisons and Analysis
// ============================================================================

/**
 * @brief Compare KE/momentum ratio for two objects
 *
 * For object with momentum p and mass m:
 * KE/p = p/(2m) = v/2
 *
 * This ratio tells us the "energy efficiency" per unit momentum.
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param velocity Velocity of object (in m/s)
 * @return KE/momentum ratio (in m/s)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateKEToMomentumRatio(double mass, double velocity) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return velocity / 2.0;
}

/**
 * @brief Calculate how mass affects KE for constant momentum
 *
 * If momentum is constant: KE₂/KE₁ = m₁/m₂
 *
 * Doubling mass halves the kinetic energy (for same momentum).
 *
 * @param mass1 First mass (in kilograms, must be > 0)
 * @param mass2 Second mass (in kilograms, must be > 0)
 * @return Ratio of kinetic energies (KE₂/KE₁)
 * @throws std::invalid_argument if either mass <= 0
 */
inline double kineticEnergyRatioConstantMomentum(double mass1, double mass2) {
    if (mass1 <= 0 || mass2 <= 0) {
        throw std::invalid_argument("Both masses must be greater than zero");
    }
    return mass1 / mass2;
}

/**
 * @brief Calculate how mass affects momentum for constant KE
 *
 * If kinetic energy is constant: p₂/p₁ = √(m₂/m₁)
 *
 * Doubling mass increases momentum by factor of √2 (for same KE).
 *
 * @param mass1 First mass (in kilograms, must be > 0)
 * @param mass2 Second mass (in kilograms, must be > 0)
 * @return Ratio of momenta (p₂/p₁)
 * @throws std::invalid_argument if either mass <= 0
 */
inline double momentumRatioConstantKE(double mass1, double mass2) {
    if (mass1 <= 0 || mass2 <= 0) {
        throw std::invalid_argument("Both masses must be greater than zero");
    }
    return std::sqrt(mass2 / mass1);
}

// ============================================================================
// Potential Energy
// ============================================================================

/**
 * @brief Calculate gravitational potential energy
 *
 * PE = mgh (near Earth's surface)
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param height Height above reference point (in meters)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Potential energy (in Joules)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculatePotentialEnergy(double mass, double height, double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return mass * gravity * height;
}

/**
 * @brief Calculate elastic potential energy in a spring
 *
 * PE_spring = (1/2)kx²
 * where k is spring constant and x is displacement from equilibrium
 *
 * @param springConstant Spring constant (in N/m, must be > 0)
 * @param displacement Displacement from equilibrium (in meters)
 * @return Elastic potential energy (in Joules)
 * @throws std::invalid_argument if springConstant <= 0
 */
inline double calculateSpringPotentialEnergy(double springConstant, double displacement) {
    if (springConstant <= 0) {
        throw std::invalid_argument("Spring constant must be positive");
    }
    return 0.5 * springConstant * displacement * displacement;
}

// ============================================================================
// Energy Conservation
// ============================================================================

/**
 * @brief Calculate final velocity from energy conservation (falling object)
 *
 * Initial: PE = mgh, KE = 0
 * Final: PE = 0, KE = (1/2)mv²
 * Conservation: mgh = (1/2)mv²
 * Therefore: v = √(2gh)
 *
 * @param height Height fallen (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Final velocity (in m/s)
 * @throws std::invalid_argument if height <= 0
 */
inline double velocityFromFall(double height, double gravity = 9.81) {
    if (height <= 0) {
        throw std::invalid_argument("Height must be positive");
    }
    return std::sqrt(2.0 * gravity * height);
}

/**
 * @brief Calculate maximum height from initial upward velocity
 *
 * Using energy conservation:
 * Initial: KE = (1/2)mv², PE = 0
 * At max height: KE = 0, PE = mgh
 * Therefore: h = v²/(2g)
 *
 * @param initialVelocity Initial upward velocity (in m/s, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Maximum height (in meters)
 * @throws std::invalid_argument if velocity <= 0
 */
inline double maxHeightFromVelocity(double initialVelocity, double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    return (initialVelocity * initialVelocity) / (2.0 * gravity);
}

// ============================================================================
// Impulse and Change in Momentum
// ============================================================================

/**
 * @brief Calculate impulse (change in momentum)
 *
 * Impulse = Δp = m(v_f - v_i)
 * Also equals: Impulse = F·Δt (impulse-momentum theorem)
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param initialVelocity Initial velocity (in m/s)
 * @param finalVelocity Final velocity (in m/s)
 * @return Impulse (in N·s or kg·m/s)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateImpulse(double mass, double initialVelocity, double finalVelocity) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return mass * (finalVelocity - initialVelocity);
}

/**
 * @brief Calculate change in kinetic energy
 *
 * ΔKE = (1/2)m(v_f² - v_i²)
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param initialVelocity Initial velocity (in m/s)
 * @param finalVelocity Final velocity (in m/s)
 * @return Change in kinetic energy (in Joules)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateChangeInKE(double mass, double initialVelocity, double finalVelocity) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return 0.5 * mass * (finalVelocity * finalVelocity - initialVelocity * initialVelocity);
}

/**
 * @brief Calculate average force from impulse and time
 *
 * From Impulse = F·Δt, we get F_avg = Impulse/Δt
 *
 * @param impulse Impulse (in N·s, magnitude)
 * @param timeInterval Time interval (in seconds, must be > 0)
 * @return Average force (in Newtons)
 * @throws std::invalid_argument if timeInterval <= 0
 */
inline double calculateAverageForceFromImpulse(double impulse, double timeInterval) {
    if (timeInterval <= 0) {
        throw std::invalid_argument("Time interval must be positive");
    }
    return impulse / timeInterval;
}

// ============================================================================
// Power
// ============================================================================

/**
 * @brief Calculate power from work and time
 *
 * Power = Work / Time
 *
 * @param work Work done (in Joules)
 * @param time Time interval (in seconds, must be > 0)
 * @return Power (in Watts)
 * @throws std::invalid_argument if time <= 0
 */
inline double calculatePowerFromWork(double work, double time) {
    if (time <= 0) {
        throw std::invalid_argument("Time must be positive");
    }
    return work / time;
}

/**
 * @brief Calculate power from force and velocity
 *
 * Power = F·v (for force parallel to velocity)
 *
 * @param force Force applied (in Newtons)
 * @param velocity Velocity (in m/s)
 * @return Power (in Watts)
 */
inline double calculatePowerFromForceVelocity(double force, double velocity) {
    return force * velocity;
}

/**
 * @brief Calculate power to lift object at constant velocity
 *
 * Power = mgh/t = mgv (where v is lifting speed)
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param liftingSpeed Constant upward velocity (in m/s, must be >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Power required (in Watts)
 * @throws std::invalid_argument if mass <= 0 or liftingSpeed < 0
 */
inline double calculatePowerToLift(double mass, double liftingSpeed, double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (liftingSpeed < 0) {
        throw std::invalid_argument("Lifting speed must be non-negative");
    }
    return mass * gravity * liftingSpeed;
}

} // namespace energy_momentum
} // namespace physics

#endif // PHYSICS_ENERGY_MOMENTUM_HPP
