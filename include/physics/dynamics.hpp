#ifndef PHYSICS_DYNAMICS_HPP
#define PHYSICS_DYNAMICS_HPP

#include <cmath>
#include <vector>
#include <stdexcept>

namespace physics {
namespace dynamics {

/**
 * @brief Dynamics: Force Causing Rectilinear Motion with Constant Acceleration
 *
 * This namespace contains functions that combine Newton's Second Law (F = ma)
 * with kinematic equations to analyze motion caused by forces.
 *
 * These functions bridge the gap between forces (dynamics) and motion (kinematics),
 * showing how forces cause acceleration and resulting motion.
 */

// ============================================================================
// Force-Acceleration Relationships
// ============================================================================

/**
 * @brief Calculate net force from multiple forces acting on an object
 *
 * For one-dimensional motion, forces in the positive direction are positive,
 * and forces in the negative direction are negative.
 *
 * @param forces Vector of forces acting on the object (in Newtons)
 * @return Net force (in Newtons)
 */
inline double calculateNetForce(const std::vector<double>& forces) {
    double netForce = 0.0;
    for (const auto& force : forces) {
        netForce += force;
    }
    return netForce;
}

/**
 * @brief Calculate acceleration caused by a net force
 *
 * Using Newton's Second Law: a = F_net / m
 *
 * @param netForce Net force acting on the object (in Newtons)
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @return Acceleration (in m/s²)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateAccelerationFromForce(double netForce, double mass) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return netForce / mass;
}

/**
 * @brief Calculate force required to produce a desired acceleration
 *
 * Using Newton's Second Law: F = ma
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param desiredAcceleration Desired acceleration (in m/s²)
 * @return Required force (in Newtons)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateRequiredForce(double mass, double desiredAcceleration) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return mass * desiredAcceleration;
}

// ============================================================================
// Force-Motion Integration: Velocity Changes
// ============================================================================

/**
 * @brief Calculate final velocity after force acts for a given time
 *
 * Combines F = ma with v = v₀ + at
 * First calculates a = F/m, then uses v = v₀ + at
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param netForce Net force acting on the object (in Newtons)
 * @param initialVelocity Initial velocity (in m/s)
 * @param time Time force acts (in seconds, must be >= 0)
 * @return Final velocity (in m/s)
 * @throws std::invalid_argument if mass <= 0 or time < 0
 */
inline double calculateFinalVelocityFromForce(double mass, double netForce,
                                               double initialVelocity, double time) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }

    double acceleration = netForce / mass;
    return initialVelocity + (acceleration * time);
}

/**
 * @brief Calculate velocity change (Δv) due to force acting for a time
 *
 * Δv = at = (F/m)t
 * This is the impulse-momentum relationship: Δv = (F·Δt) / m
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param netForce Net force acting on the object (in Newtons)
 * @param time Time force acts (in seconds, must be >= 0)
 * @return Change in velocity (in m/s)
 * @throws std::invalid_argument if mass <= 0 or time < 0
 */
inline double calculateVelocityChange(double mass, double netForce, double time) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }

    double acceleration = netForce / mass;
    return acceleration * time;
}

/**
 * @brief Calculate time required for force to produce a velocity change
 *
 * From Δv = (F/m)t, we get t = m·Δv / F
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param netForce Net force acting on the object (in Newtons, must be != 0)
 * @param velocityChange Desired change in velocity (in m/s)
 * @return Time required (in seconds)
 * @throws std::invalid_argument if mass <= 0 or netForce == 0
 */
inline double calculateTimeForVelocityChange(double mass, double netForce, double velocityChange) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (std::abs(netForce) < 1e-10) {
        throw std::invalid_argument("Net force must be non-zero");
    }

    return (mass * velocityChange) / netForce;
}

// ============================================================================
// Force-Motion Integration: Displacement
// ============================================================================

/**
 * @brief Calculate displacement when force acts on object for given time
 *
 * Combines F = ma with s = v₀t + (1/2)at²
 * First calculates a = F/m, then uses displacement equation
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param netForce Net force acting on the object (in Newtons)
 * @param initialVelocity Initial velocity (in m/s)
 * @param time Time force acts (in seconds, must be >= 0)
 * @return Displacement (in meters)
 * @throws std::invalid_argument if mass <= 0 or time < 0
 */
inline double calculateDisplacementFromForce(double mass, double netForce,
                                              double initialVelocity, double time) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }

    double acceleration = netForce / mass;
    return (initialVelocity * time) + (0.5 * acceleration * time * time);
}

/**
 * @brief Calculate final velocity after force causes displacement
 *
 * Combines F = ma with v² = v₀² + 2as
 * First calculates a = F/m, then uses v = √(v₀² + 2as)
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param netForce Net force acting on the object (in Newtons)
 * @param initialVelocity Initial velocity (in m/s)
 * @param displacement Displacement (in meters)
 * @return Final velocity magnitude (in m/s)
 * @throws std::invalid_argument if mass <= 0 or if v² < 0 (no real solution)
 */
inline double calculateFinalVelocityFromForceAndDisplacement(double mass, double netForce,
                                                              double initialVelocity,
                                                              double displacement) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }

    double acceleration = netForce / mass;
    double vSquared = (initialVelocity * initialVelocity) + (2.0 * acceleration * displacement);

    if (vSquared < 0) {
        throw std::invalid_argument("No real solution: v² cannot be negative");
    }

    return std::sqrt(vSquared);
}

// ============================================================================
// Force Analysis: Stopping Problems
// ============================================================================

/**
 * @brief Calculate braking force required to stop object in given distance
 *
 * Uses v² = v₀² + 2as with v = 0 (stopping)
 * a = -v₀² / (2s), then F = ma
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param stoppingDistance Distance available for stopping (in meters, must be > 0)
 * @return Braking force required (in Newtons, negative indicating opposite to motion)
 * @throws std::invalid_argument if mass <= 0, initialVelocity <= 0, or stoppingDistance <= 0
 */
inline double calculateBrakingForce(double mass, double initialVelocity, double stoppingDistance) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (stoppingDistance <= 0) {
        throw std::invalid_argument("Stopping distance must be positive");
    }

    // Deceleration needed: a = -v₀² / (2s)
    double acceleration = -(initialVelocity * initialVelocity) / (2.0 * stoppingDistance);
    return mass * acceleration;
}

/**
 * @brief Calculate stopping distance for given braking force
 *
 * Uses F = ma to find a, then s = v₀² / (2|a|) to find stopping distance
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param brakingForce Braking force magnitude (in Newtons, must be > 0, always positive)
 * @return Stopping distance (in meters)
 * @throws std::invalid_argument if mass <= 0, initialVelocity <= 0, or brakingForce <= 0
 */
inline double calculateStoppingDistanceFromForce(double mass, double initialVelocity,
                                                  double brakingForce) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (brakingForce <= 0) {
        throw std::invalid_argument("Braking force must be positive");
    }

    double deceleration = brakingForce / mass;
    return (initialVelocity * initialVelocity) / (2.0 * deceleration);
}

/**
 * @brief Calculate stopping time for given braking force
 *
 * Uses F = ma to find a, then t = v₀ / |a| to find stopping time
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param brakingForce Braking force magnitude (in Newtons, must be > 0, always positive)
 * @return Stopping time (in seconds)
 * @throws std::invalid_argument if mass <= 0, initialVelocity <= 0, or brakingForce <= 0
 */
inline double calculateStoppingTimeFromForce(double mass, double initialVelocity,
                                              double brakingForce) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (brakingForce <= 0) {
        throw std::invalid_argument("Braking force must be positive");
    }

    double deceleration = brakingForce / mass;
    return initialVelocity / deceleration;
}

// ============================================================================
// Friction and Motion
// ============================================================================

/**
 * @brief Calculate frictional force
 *
 * Friction force: f = μN, where μ is coefficient of friction and N is normal force
 * For horizontal surfaces, N = mg (weight)
 *
 * @param coefficientOfFriction Coefficient of friction (dimensionless, must be >= 0)
 * @param normalForce Normal force (in Newtons, must be >= 0)
 * @return Frictional force (in Newtons)
 * @throws std::invalid_argument if coefficientOfFriction < 0 or normalForce < 0
 */
inline double calculateFrictionForce(double coefficientOfFriction, double normalForce) {
    if (coefficientOfFriction < 0) {
        throw std::invalid_argument("Coefficient of friction must be non-negative");
    }
    if (normalForce < 0) {
        throw std::invalid_argument("Normal force must be non-negative");
    }

    return coefficientOfFriction * normalForce;
}

/**
 * @brief Calculate acceleration of object on horizontal surface with friction
 *
 * For object on horizontal surface: N = mg
 * Net force: F_net = F_applied - f = F_applied - μmg
 * Acceleration: a = F_net / m = (F_applied / m) - μg
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param appliedForce Applied force (in Newtons)
 * @param coefficientOfFriction Coefficient of friction (dimensionless, must be >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Acceleration (in m/s²)
 * @throws std::invalid_argument if mass <= 0 or coefficientOfFriction < 0
 */
inline double calculateAccelerationWithFriction(double mass, double appliedForce,
                                                double coefficientOfFriction,
                                                double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (coefficientOfFriction < 0) {
        throw std::invalid_argument("Coefficient of friction must be non-negative");
    }

    double frictionForce = coefficientOfFriction * mass * gravity;
    double netForce = appliedForce - frictionForce;
    return netForce / mass;
}

/**
 * @brief Calculate minimum force to overcome static friction and start motion
 *
 * For object to start moving: F_applied >= f_static_max = μ_s × N
 * On horizontal surface: N = mg
 * Therefore: F_min = μ_s × mg
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param staticFrictionCoefficient Coefficient of static friction (dimensionless, must be >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Minimum force to overcome friction (in Newtons)
 * @throws std::invalid_argument if mass <= 0 or staticFrictionCoefficient < 0
 */
inline double calculateMinimumForceToOvercomeFriction(double mass,
                                                      double staticFrictionCoefficient,
                                                      double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (staticFrictionCoefficient < 0) {
        throw std::invalid_argument("Static friction coefficient must be non-negative");
    }

    return staticFrictionCoefficient * mass * gravity;
}

// ============================================================================
// Work-Energy Considerations (preview for future extensions)
// ============================================================================

/**
 * @brief Calculate work done by constant force over displacement
 *
 * Work: W = F · s (for force parallel to displacement)
 * Work represents energy transferred to/from the object
 *
 * @param force Force applied (in Newtons)
 * @param displacement Displacement (in meters)
 * @return Work done (in Joules)
 */
inline double calculateWork(double force, double displacement) {
    return force * displacement;
}

/**
 * @brief Calculate power (rate of work) for constant force and velocity
 *
 * Power: P = F · v (for force parallel to velocity)
 * Power represents rate of energy transfer
 *
 * @param force Force applied (in Newtons)
 * @param velocity Velocity (in m/s)
 * @return Power (in Watts)
 */
inline double calculatePower(double force, double velocity) {
    return force * velocity;
}

} // namespace dynamics
} // namespace physics

#endif // PHYSICS_DYNAMICS_HPP
