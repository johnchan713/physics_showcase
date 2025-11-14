#ifndef PHYSICS_KINEMATICS_HPP
#define PHYSICS_KINEMATICS_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace kinematics {

/**
 * @brief Motion in a Straight Line with Constant Acceleration
 *
 * This namespace contains functions for analyzing one-dimensional motion
 * with constant acceleration using the kinematic equations.
 *
 * Standard notation:
 * - s or Δx: displacement (in meters)
 * - v₀ or u: initial velocity (in m/s)
 * - v: final velocity (in m/s)
 * - a: acceleration (in m/s²)
 * - t: time (in seconds)
 */

// ============================================================================
// First Kinematic Equation: v = v₀ + at
// ============================================================================

/**
 * @brief Calculate final velocity using: v = v₀ + at
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param acceleration Constant acceleration (in m/s²)
 * @param time Time interval (in seconds, must be >= 0)
 * @return Final velocity (in m/s)
 * @throws std::invalid_argument if time < 0
 */
inline double calculateFinalVelocity(double initialVelocity, double acceleration, double time) {
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return initialVelocity + (acceleration * time);
}

/**
 * @brief Calculate acceleration using: a = (v - v₀) / t
 *
 * Derived from v = v₀ + at
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param finalVelocity Final velocity (in m/s)
 * @param time Time interval (in seconds, must be > 0)
 * @return Acceleration (in m/s²)
 * @throws std::invalid_argument if time <= 0
 */
inline double calculateAccelerationFromVelocities(double initialVelocity, double finalVelocity, double time) {
    if (time <= 0) {
        throw std::invalid_argument("Time must be positive");
    }
    return (finalVelocity - initialVelocity) / time;
}

/**
 * @brief Calculate time required to reach final velocity: t = (v - v₀) / a
 *
 * Derived from v = v₀ + at
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param finalVelocity Final velocity (in m/s)
 * @param acceleration Constant acceleration (in m/s², must be != 0)
 * @return Time interval (in seconds)
 * @throws std::invalid_argument if acceleration == 0
 */
inline double calculateTimeFromVelocities(double initialVelocity, double finalVelocity, double acceleration) {
    if (std::abs(acceleration) < 1e-10) {
        throw std::invalid_argument("Acceleration must be non-zero");
    }
    return (finalVelocity - initialVelocity) / acceleration;
}

// ============================================================================
// Second Kinematic Equation: s = v₀t + (1/2)at²
// ============================================================================

/**
 * @brief Calculate displacement using: s = v₀t + (1/2)at²
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param acceleration Constant acceleration (in m/s²)
 * @param time Time interval (in seconds, must be >= 0)
 * @return Displacement (in meters)
 * @throws std::invalid_argument if time < 0
 */
inline double calculateDisplacement(double initialVelocity, double acceleration, double time) {
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return (initialVelocity * time) + (0.5 * acceleration * time * time);
}

/**
 * @brief Calculate acceleration from displacement: a = 2(s - v₀t) / t²
 *
 * Derived from s = v₀t + (1/2)at²
 *
 * @param displacement Displacement (in meters)
 * @param initialVelocity Initial velocity (in m/s)
 * @param time Time interval (in seconds, must be > 0)
 * @return Acceleration (in m/s²)
 * @throws std::invalid_argument if time <= 0
 */
inline double calculateAccelerationFromDisplacement(double displacement, double initialVelocity, double time) {
    if (time <= 0) {
        throw std::invalid_argument("Time must be positive");
    }
    return (2.0 * (displacement - (initialVelocity * time))) / (time * time);
}

// ============================================================================
// Third Kinematic Equation: v² = v₀² + 2as
// ============================================================================

/**
 * @brief Calculate final velocity using: v = √(v₀² + 2as)
 *
 * Note: Returns the positive root. For objects moving in negative direction,
 * the result should be negated.
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param acceleration Constant acceleration (in m/s²)
 * @param displacement Displacement (in meters)
 * @return Final velocity magnitude (in m/s)
 * @throws std::invalid_argument if v₀² + 2as < 0 (no real solution)
 */
inline double calculateFinalVelocityFromDisplacement(double initialVelocity, double acceleration, double displacement) {
    double vSquared = (initialVelocity * initialVelocity) + (2.0 * acceleration * displacement);
    if (vSquared < 0) {
        throw std::invalid_argument("No real solution: v² cannot be negative");
    }
    return std::sqrt(vSquared);
}

/**
 * @brief Calculate acceleration using: a = (v² - v₀²) / (2s)
 *
 * Derived from v² = v₀² + 2as
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param finalVelocity Final velocity (in m/s)
 * @param displacement Displacement (in meters, must be != 0)
 * @return Acceleration (in m/s²)
 * @throws std::invalid_argument if displacement == 0
 */
inline double calculateAccelerationFromVelocitySquared(double initialVelocity, double finalVelocity, double displacement) {
    if (std::abs(displacement) < 1e-10) {
        throw std::invalid_argument("Displacement must be non-zero");
    }
    return ((finalVelocity * finalVelocity) - (initialVelocity * initialVelocity)) / (2.0 * displacement);
}

/**
 * @brief Calculate displacement using: s = (v² - v₀²) / (2a)
 *
 * Derived from v² = v₀² + 2as
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param finalVelocity Final velocity (in m/s)
 * @param acceleration Constant acceleration (in m/s², must be != 0)
 * @return Displacement (in meters)
 * @throws std::invalid_argument if acceleration == 0
 */
inline double calculateDisplacementFromVelocities(double initialVelocity, double finalVelocity, double acceleration) {
    if (std::abs(acceleration) < 1e-10) {
        throw std::invalid_argument("Acceleration must be non-zero");
    }
    return ((finalVelocity * finalVelocity) - (initialVelocity * initialVelocity)) / (2.0 * acceleration);
}

// ============================================================================
// Fourth Kinematic Equation: s = ((v + v₀) / 2) * t
// ============================================================================

/**
 * @brief Calculate displacement using average velocity: s = ((v + v₀) / 2) * t
 *
 * This equation is valid only for constant acceleration, where average
 * velocity equals (initial velocity + final velocity) / 2
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param finalVelocity Final velocity (in m/s)
 * @param time Time interval (in seconds, must be >= 0)
 * @return Displacement (in meters)
 * @throws std::invalid_argument if time < 0
 */
inline double calculateDisplacementFromAverageVelocity(double initialVelocity, double finalVelocity, double time) {
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return ((initialVelocity + finalVelocity) / 2.0) * time;
}

/**
 * @brief Calculate average velocity for constant acceleration
 *
 * For motion with constant acceleration: v_avg = (v₀ + v) / 2
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param finalVelocity Final velocity (in m/s)
 * @return Average velocity (in m/s)
 */
inline double calculateAverageVelocity(double initialVelocity, double finalVelocity) {
    return (initialVelocity + finalVelocity) / 2.0;
}

/**
 * @brief Calculate time from displacement and average velocity: t = s / v_avg
 *
 * @param displacement Displacement (in meters)
 * @param initialVelocity Initial velocity (in m/s)
 * @param finalVelocity Final velocity (in m/s)
 * @return Time interval (in seconds)
 * @throws std::invalid_argument if average velocity == 0
 */
inline double calculateTimeFromAverageVelocity(double displacement, double initialVelocity, double finalVelocity) {
    double avgVelocity = calculateAverageVelocity(initialVelocity, finalVelocity);
    if (std::abs(avgVelocity) < 1e-10) {
        throw std::invalid_argument("Average velocity must be non-zero");
    }
    return displacement / avgVelocity;
}

// ============================================================================
// Additional Utility Functions
// ============================================================================

/**
 * @brief Calculate distance traveled (always positive)
 *
 * Unlike displacement which can be negative, distance is the total
 * path length traveled, always positive.
 *
 * For constant acceleration in one direction, distance = |displacement|
 *
 * @param initialVelocity Initial velocity (in m/s)
 * @param acceleration Constant acceleration (in m/s²)
 * @param time Time interval (in seconds, must be >= 0)
 * @return Distance traveled (in meters, always >= 0)
 * @throws std::invalid_argument if time < 0
 */
inline double calculateDistance(double initialVelocity, double acceleration, double time) {
    return std::abs(calculateDisplacement(initialVelocity, acceleration, time));
}

/**
 * @brief Calculate stopping distance for a decelerating object
 *
 * Calculates the distance required to bring an object to rest
 * Uses: s = -v₀² / (2a), where a is negative (deceleration)
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param deceleration Deceleration magnitude (in m/s², must be > 0, always positive)
 * @return Stopping distance (in meters)
 * @throws std::invalid_argument if initialVelocity <= 0 or deceleration <= 0
 */
inline double calculateStoppingDistance(double initialVelocity, double deceleration) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (deceleration <= 0) {
        throw std::invalid_argument("Deceleration must be positive");
    }
    return (initialVelocity * initialVelocity) / (2.0 * deceleration);
}

/**
 * @brief Calculate stopping time for a decelerating object
 *
 * Calculates the time required to bring an object to rest
 * Uses: t = v₀ / a, where a is the deceleration magnitude
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param deceleration Deceleration magnitude (in m/s², must be > 0, always positive)
 * @return Stopping time (in seconds)
 * @throws std::invalid_argument if initialVelocity <= 0 or deceleration <= 0
 */
inline double calculateStoppingTime(double initialVelocity, double deceleration) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (deceleration <= 0) {
        throw std::invalid_argument("Deceleration must be positive");
    }
    return initialVelocity / deceleration;
}

} // namespace kinematics
} // namespace physics

#endif // PHYSICS_KINEMATICS_HPP
