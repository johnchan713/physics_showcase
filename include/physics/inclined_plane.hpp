#ifndef PHYSICS_INCLINED_PLANE_HPP
#define PHYSICS_INCLINED_PLANE_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace inclined_plane {

/**
 * @brief Inclined Plane Mechanics
 *
 * Functions for analyzing motion on inclined planes, including:
 * - Force components parallel and perpendicular to the plane
 * - Acceleration down the plane
 * - Velocity at the foot of the plane
 * - Effects of friction
 *
 * Standard notation:
 * - θ (theta): angle of inclination (in radians)
 * - m: mass of object
 * - g: gravitational acceleration
 * - μ (mu): coefficient of friction
 */

// ============================================================================
// Force Components on Inclined Plane
// ============================================================================

/**
 * @brief Calculate component of weight parallel to inclined plane
 *
 * The parallel component causes the object to slide down: F_parallel = mg sin(θ)
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Force parallel to plane (in Newtons)
 * @throws std::invalid_argument if mass <= 0 or angle out of range
 */
inline double calculateParallelForce(double mass, double angleRadians, double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }
    return mass * gravity * std::sin(angleRadians);
}

/**
 * @brief Calculate component of weight perpendicular to inclined plane
 *
 * The perpendicular component equals the normal force: F_perpendicular = mg cos(θ)
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Force perpendicular to plane (Normal force in Newtons)
 * @throws std::invalid_argument if mass <= 0 or angle out of range
 */
inline double calculateNormalForce(double mass, double angleRadians, double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }
    return mass * gravity * std::cos(angleRadians);
}

// ============================================================================
// Acceleration on Inclined Plane
// ============================================================================

/**
 * @brief Calculate acceleration down frictionless inclined plane
 *
 * Without friction: a = g sin(θ)
 * Acceleration is independent of mass!
 *
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Acceleration down the plane (in m/s²)
 * @throws std::invalid_argument if angle out of range
 */
inline double calculateAccelerationFrictionless(double angleRadians, double gravity = 9.81) {
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }
    return gravity * std::sin(angleRadians);
}

/**
 * @brief Calculate acceleration down inclined plane with friction
 *
 * With friction: a = g(sin(θ) - μ cos(θ))
 * where μ is the coefficient of kinetic friction
 *
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param frictionCoefficient Coefficient of kinetic friction (dimensionless, >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Acceleration down the plane (in m/s²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateAccelerationWithFriction(double angleRadians,
                                                double frictionCoefficient,
                                                double gravity = 9.81) {
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }
    if (frictionCoefficient < 0) {
        throw std::invalid_argument("Friction coefficient must be non-negative");
    }
    return gravity * (std::sin(angleRadians) - frictionCoefficient * std::cos(angleRadians));
}

/**
 * @brief Calculate minimum angle for object to slide (with static friction)
 *
 * Object starts sliding when: tan(θ) > μ_s
 * Therefore: θ_min = arctan(μ_s)
 *
 * @param staticFrictionCoefficient Coefficient of static friction (dimensionless, must be > 0)
 * @return Minimum angle for sliding (in radians)
 * @throws std::invalid_argument if coefficient <= 0
 */
inline double calculateMinimumAngleToSlide(double staticFrictionCoefficient) {
    if (staticFrictionCoefficient <= 0) {
        throw std::invalid_argument("Static friction coefficient must be positive");
    }
    return std::atan(staticFrictionCoefficient);
}

// ============================================================================
// Velocity at Foot of Inclined Plane
// ============================================================================

/**
 * @brief Calculate velocity at foot of frictionless inclined plane
 *
 * Object starts from rest at height h and slides down.
 * Using energy conservation: mgh = (1/2)mv²
 * Therefore: v = √(2gh)
 *
 * Alternatively, using kinematics with a = g sin(θ) and s = h/sin(θ):
 * v² = 2as = 2g sin(θ) × (h/sin(θ)) = 2gh
 *
 * @param height Vertical height of incline (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Velocity at bottom (in m/s)
 * @throws std::invalid_argument if height <= 0
 */
inline double calculateVelocityAtFootFromHeight(double height, double gravity = 9.81) {
    if (height <= 0) {
        throw std::invalid_argument("Height must be positive");
    }
    return std::sqrt(2.0 * gravity * height);
}

/**
 * @brief Calculate velocity at foot of inclined plane from length and angle
 *
 * Object starts from rest and slides down length L at angle θ.
 * Using kinematics: v² = 2as, where a = g sin(θ) and s = L
 * Therefore: v = √(2gL sin(θ))
 *
 * @param length Length along the plane (in meters, must be > 0)
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Velocity at bottom (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateVelocityAtFootFrictionless(double length, double angleRadians,
                                                   double gravity = 9.81) {
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }

    double acceleration = calculateAccelerationFrictionless(angleRadians, gravity);
    return std::sqrt(2.0 * acceleration * length);
}

/**
 * @brief Calculate velocity at foot of inclined plane with friction
 *
 * Object starts from rest and slides down length L with friction.
 * Using kinematics: v² = 2as, where a = g(sin(θ) - μ cos(θ)) and s = L
 * Therefore: v = √(2gL(sin(θ) - μ cos(θ)))
 *
 * @param length Length along the plane (in meters, must be > 0)
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param frictionCoefficient Coefficient of kinetic friction (dimensionless, >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Velocity at bottom (in m/s)
 * @throws std::invalid_argument if parameters out of range or if friction too high
 */
inline double calculateVelocityAtFootWithFriction(double length, double angleRadians,
                                                  double frictionCoefficient,
                                                  double gravity = 9.81) {
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }
    if (frictionCoefficient < 0) {
        throw std::invalid_argument("Friction coefficient must be non-negative");
    }

    double acceleration = calculateAccelerationWithFriction(angleRadians, frictionCoefficient, gravity);

    if (acceleration <= 0) {
        throw std::invalid_argument("Friction too high: object won't slide");
    }

    return std::sqrt(2.0 * acceleration * length);
}

// ============================================================================
// Time to Slide Down Inclined Plane
// ============================================================================

/**
 * @brief Calculate time to slide down frictionless inclined plane
 *
 * Starting from rest: s = (1/2)at²
 * Therefore: t = √(2s/a) = √(2L/(g sin(θ)))
 *
 * @param length Length along the plane (in meters, must be > 0)
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Time to reach bottom (in seconds)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateTimeToSlide(double length, double angleRadians, double gravity = 9.81) {
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }

    double acceleration = calculateAccelerationFrictionless(angleRadians, gravity);
    return std::sqrt(2.0 * length / acceleration);
}

/**
 * @brief Calculate time to slide down inclined plane with friction
 *
 * Starting from rest: t = √(2L/a), where a = g(sin(θ) - μ cos(θ))
 *
 * @param length Length along the plane (in meters, must be > 0)
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param frictionCoefficient Coefficient of kinetic friction (dimensionless, >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Time to reach bottom (in seconds)
 * @throws std::invalid_argument if parameters out of range or if friction too high
 */
inline double calculateTimeToSlideWithFriction(double length, double angleRadians,
                                               double frictionCoefficient,
                                               double gravity = 9.81) {
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }
    if (frictionCoefficient < 0) {
        throw std::invalid_argument("Friction coefficient must be non-negative");
    }

    double acceleration = calculateAccelerationWithFriction(angleRadians, frictionCoefficient, gravity);

    if (acceleration <= 0) {
        throw std::invalid_argument("Friction too high: object won't slide");
    }

    return std::sqrt(2.0 * length / acceleration);
}

// ============================================================================
// Energy Considerations
// ============================================================================

/**
 * @brief Calculate work done against friction on inclined plane
 *
 * Work against friction: W_friction = μ × N × L = μ mg cos(θ) × L
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param length Length along the plane (in meters, must be > 0)
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param frictionCoefficient Coefficient of kinetic friction (dimensionless, >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Work done against friction (in Joules, positive value)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateWorkAgainstFriction(double mass, double length, double angleRadians,
                                           double frictionCoefficient, double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }
    if (frictionCoefficient < 0) {
        throw std::invalid_argument("Friction coefficient must be non-negative");
    }

    double normalForce = calculateNormalForce(mass, angleRadians, gravity);
    return frictionCoefficient * normalForce * length;
}

/**
 * @brief Calculate energy loss to friction as fraction of initial potential energy
 *
 * Energy loss ratio = (Work against friction) / (Initial PE)
 *                   = (μ mg cos(θ) L) / (mgh)
 *                   = (μ mg cos(θ) L) / (mg L sin(θ))
 *                   = μ cot(θ)
 *
 * @param angleRadians Angle of inclination (in radians, 0 < θ < π/2)
 * @param frictionCoefficient Coefficient of kinetic friction (dimensionless, >= 0)
 * @return Fraction of energy lost to friction (0 to 1)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateEnergyLossFraction(double angleRadians, double frictionCoefficient) {
    if (angleRadians <= 0 || angleRadians >= M_PI / 2.0) {
        throw std::invalid_argument("Angle must be between 0 and π/2 radians");
    }
    if (frictionCoefficient < 0) {
        throw std::invalid_argument("Friction coefficient must be non-negative");
    }

    return frictionCoefficient / std::tan(angleRadians);
}

} // namespace inclined_plane
} // namespace physics

#endif // PHYSICS_INCLINED_PLANE_HPP
