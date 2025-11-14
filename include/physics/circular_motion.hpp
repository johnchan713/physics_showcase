#ifndef PHYSICS_CIRCULAR_MOTION_HPP
#define PHYSICS_CIRCULAR_MOTION_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace circular_motion {

/**
 * @brief Motion in a Circle with Constant Speed (Uniform Circular Motion)
 *
 * Functions for analyzing uniform circular motion where speed is constant
 * but velocity direction changes continuously.
 *
 * Key concepts:
 * - Speed is constant, velocity is not (velocity changes direction)
 * - Centripetal acceleration directed toward center: a_c = v²/r
 * - Centripetal force required: F_c = mv²/r
 * - Angular velocity: ω = v/r
 * - Period: T = 2πr/v
 *
 * Standard notation:
 * - r: radius of circular path
 * - v: tangential (linear) speed
 * - ω (omega): angular velocity
 * - T: period (time for one complete revolution)
 * - f: frequency (revolutions per unit time)
 */

// ============================================================================
// Basic Circular Motion Parameters
// ============================================================================

/**
 * @brief Calculate centripetal acceleration
 *
 * Centripetal acceleration: a_c = v²/r = ω²r
 *
 * This acceleration is always directed toward the center of the circle.
 *
 * @param velocity Tangential speed (in m/s, must be >= 0)
 * @param radius Radius of circular path (in meters, must be > 0)
 * @return Centripetal acceleration (in m/s²)
 * @throws std::invalid_argument if velocity < 0 or radius <= 0
 */
inline double calculateCentripetalAcceleration(double velocity, double radius) {
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return (velocity * velocity) / radius;
}

/**
 * @brief Calculate centripetal force required
 *
 * Centripetal force: F_c = mv²/r = mω²r
 *
 * This is the net force toward center needed to maintain circular motion.
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param velocity Tangential speed (in m/s, must be >= 0)
 * @param radius Radius of circular path (in meters, must be > 0)
 * @return Centripetal force (in Newtons)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateCentripetalForce(double mass, double velocity, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return mass * (velocity * velocity) / radius;
}

/**
 * @brief Calculate velocity for given centripetal acceleration
 *
 * From a_c = v²/r, we get v = √(a_c × r)
 *
 * @param centripetalAccel Centripetal acceleration (in m/s², must be >= 0)
 * @param radius Radius of circular path (in meters, must be > 0)
 * @return Tangential velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateVelocityFromAcceleration(double centripetalAccel, double radius) {
    if (centripetalAccel < 0) {
        throw std::invalid_argument("Centripetal acceleration must be non-negative");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return std::sqrt(centripetalAccel * radius);
}

// ============================================================================
// Angular Motion
// ============================================================================

/**
 * @brief Calculate angular velocity
 *
 * Angular velocity: ω = v/r (in radians per second)
 *
 * Also: ω = 2π/T = 2πf
 *
 * @param velocity Tangential speed (in m/s, must be >= 0)
 * @param radius Radius of circular path (in meters, must be > 0)
 * @return Angular velocity (in radians/second)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateAngularVelocity(double velocity, double radius) {
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return velocity / radius;
}

/**
 * @brief Calculate tangential velocity from angular velocity
 *
 * Tangential velocity: v = ωr
 *
 * @param angularVelocity Angular velocity (in rad/s, must be >= 0)
 * @param radius Radius of circular path (in meters, must be > 0)
 * @return Tangential velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateTangentialVelocity(double angularVelocity, double radius) {
    if (angularVelocity < 0) {
        throw std::invalid_argument("Angular velocity must be non-negative");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return angularVelocity * radius;
}

/**
 * @brief Calculate centripetal acceleration from angular velocity
 *
 * a_c = ω²r
 *
 * @param angularVelocity Angular velocity (in rad/s, must be >= 0)
 * @param radius Radius of circular path (in meters, must be > 0)
 * @return Centripetal acceleration (in m/s²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateCentripetalAccelFromAngular(double angularVelocity, double radius) {
    if (angularVelocity < 0) {
        throw std::invalid_argument("Angular velocity must be non-negative");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return angularVelocity * angularVelocity * radius;
}

// ============================================================================
// Period and Frequency
// ============================================================================

/**
 * @brief Calculate period (time for one complete revolution)
 *
 * Period: T = 2πr/v = 2π/ω
 *
 * @param velocity Tangential speed (in m/s, must be > 0)
 * @param radius Radius of circular path (in meters, must be > 0)
 * @return Period (in seconds)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculatePeriod(double velocity, double radius) {
    if (velocity <= 0) {
        throw std::invalid_argument("Velocity must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return (2.0 * M_PI * radius) / velocity;
}

/**
 * @brief Calculate period from angular velocity
 *
 * Period: T = 2π/ω
 *
 * @param angularVelocity Angular velocity (in rad/s, must be > 0)
 * @return Period (in seconds)
 * @throws std::invalid_argument if angularVelocity <= 0
 */
inline double calculatePeriodFromAngular(double angularVelocity) {
    if (angularVelocity <= 0) {
        throw std::invalid_argument("Angular velocity must be positive");
    }
    return (2.0 * M_PI) / angularVelocity;
}

/**
 * @brief Calculate frequency (revolutions per second)
 *
 * Frequency: f = 1/T = v/(2πr) = ω/(2π)
 *
 * @param velocity Tangential speed (in m/s, must be > 0)
 * @param radius Radius of circular path (in meters, must be > 0)
 * @return Frequency (in Hz or rev/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateFrequency(double velocity, double radius) {
    if (velocity <= 0) {
        throw std::invalid_argument("Velocity must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return velocity / (2.0 * M_PI * radius);
}

/**
 * @brief Calculate angular velocity from period
 *
 * ω = 2π/T
 *
 * @param period Period (in seconds, must be > 0)
 * @return Angular velocity (in rad/s)
 * @throws std::invalid_argument if period <= 0
 */
inline double calculateAngularVelocityFromPeriod(double period) {
    if (period <= 0) {
        throw std::invalid_argument("Period must be positive");
    }
    return (2.0 * M_PI) / period;
}

/**
 * @brief Calculate angular velocity from frequency
 *
 * ω = 2πf
 *
 * @param frequency Frequency (in Hz, must be > 0)
 * @return Angular velocity (in rad/s)
 * @throws std::invalid_argument if frequency <= 0
 */
inline double calculateAngularVelocityFromFrequency(double frequency) {
    if (frequency <= 0) {
        throw std::invalid_argument("Frequency must be positive");
    }
    return 2.0 * M_PI * frequency;
}

// ============================================================================
// Energy in Circular Motion
// ============================================================================

/**
 * @brief Calculate kinetic energy in circular motion
 *
 * KE = (1/2)mv² (same as linear motion)
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param velocity Tangential speed (in m/s, must be >= 0)
 * @return Kinetic energy (in Joules)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateKineticEnergy(double mass, double velocity) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    return 0.5 * mass * velocity * velocity;
}

/**
 * @brief Calculate rotational kinetic energy using angular velocity
 *
 * For a point mass: KE = (1/2)mv² = (1/2)m(ωr)² = (1/2)(mr²)ω²
 *
 * Note: mr² is the moment of inertia for a point mass
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param angularVelocity Angular velocity (in rad/s, must be >= 0)
 * @param radius Radius of circular path (in meters, must be > 0)
 * @return Kinetic energy (in Joules)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateKineticEnergyAngular(double mass, double angularVelocity, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (angularVelocity < 0) {
        throw std::invalid_argument("Angular velocity must be non-negative");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return 0.5 * mass * radius * radius * angularVelocity * angularVelocity;
}

// ============================================================================
// Banking and Vertical Circles
// ============================================================================

/**
 * @brief Calculate banking angle for frictionless circular turn
 *
 * For a car on a banked curve with no friction:
 * tan(θ) = v²/(rg)
 *
 * @param velocity Speed of vehicle (in m/s, must be > 0)
 * @param radius Radius of curve (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Banking angle (in radians)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateBankingAngle(double velocity, double radius, double gravity = 9.81) {
    if (velocity <= 0) {
        throw std::invalid_argument("Velocity must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return std::atan((velocity * velocity) / (radius * gravity));
}

/**
 * @brief Calculate minimum velocity at top of vertical circle
 *
 * For object on inside of vertical circle (like loop-the-loop):
 * At top, minimum velocity to maintain contact: v_min = √(gr)
 *
 * Below this speed, the object loses contact with the track.
 *
 * @param radius Radius of vertical circle (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Minimum velocity at top (in m/s)
 * @throws std::invalid_argument if radius <= 0
 */
inline double calculateMinVelocityTopOfLoop(double radius, double gravity = 9.81) {
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return std::sqrt(gravity * radius);
}

/**
 * @brief Calculate tension in string for object in vertical circle
 *
 * At top of circle: T = mv²/r - mg (tension can be zero if v = √(gr))
 * At bottom: T = mv²/r + mg
 *
 * This function calculates tension at the bottom.
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param velocity Speed at bottom (in m/s, must be >= 0)
 * @param radius Radius of circle (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Tension at bottom (in Newtons)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateTensionAtBottom(double mass, double velocity, double radius,
                                       double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }

    double centripetalForce = calculateCentripetalForce(mass, velocity, radius);
    double weight = mass * gravity;
    return centripetalForce + weight;
}

/**
 * @brief Calculate tension in string at top of vertical circle
 *
 * At top: T = mv²/r - mg
 *
 * If T < 0, the object has lost contact (only valid if v >= √(gr))
 *
 * @param mass Mass of object (in kilograms, must be > 0)
 * @param velocity Speed at top (in m/s, must be >= √(gr))
 * @param radius Radius of circle (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Tension at top (in Newtons)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateTensionAtTop(double mass, double velocity, double radius,
                                    double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }

    double centripetalForce = calculateCentripetalForce(mass, velocity, radius);
    double weight = mass * gravity;
    return centripetalForce - weight;
}

// ============================================================================
// Conversions
// ============================================================================

/**
 * @brief Convert RPM (revolutions per minute) to radians per second
 *
 * ω(rad/s) = RPM × (2π/60)
 *
 * @param rpm Revolutions per minute (must be >= 0)
 * @return Angular velocity (in rad/s)
 * @throws std::invalid_argument if rpm < 0
 */
inline double rpmToRadPerSec(double rpm) {
    if (rpm < 0) {
        throw std::invalid_argument("RPM must be non-negative");
    }
    return rpm * (2.0 * M_PI / 60.0);
}

/**
 * @brief Convert radians per second to RPM
 *
 * RPM = ω(rad/s) × (60/2π)
 *
 * @param radPerSec Angular velocity (in rad/s, must be >= 0)
 * @return Revolutions per minute
 * @throws std::invalid_argument if radPerSec < 0
 */
inline double radPerSecToRpm(double radPerSec) {
    if (radPerSec < 0) {
        throw std::invalid_argument("Angular velocity must be non-negative");
    }
    return radPerSec * (60.0 / (2.0 * M_PI));
}

} // namespace circular_motion
} // namespace physics

#endif // PHYSICS_CIRCULAR_MOTION_HPP
