#ifndef PHYSICS_PROJECTILE_HPP
#define PHYSICS_PROJECTILE_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace projectile {

/**
 * @brief Projectile Motion
 *
 * Functions for analyzing 2D projectile motion under constant gravitational acceleration.
 * Assumes:
 * - No air resistance
 * - Constant gravitational acceleration
 * - Motion in a vertical plane
 *
 * Standard notation:
 * - v₀: initial velocity (magnitude)
 * - θ: launch angle (measured from horizontal)
 * - g: gravitational acceleration
 * - Horizontal: x-direction (no acceleration)
 * - Vertical: y-direction (acceleration = -g)
 */

// ============================================================================
// Initial Velocity Components
// ============================================================================

/**
 * @brief Calculate horizontal component of initial velocity
 *
 * v₀ₓ = v₀ cos(θ)
 *
 * @param initialVelocity Initial velocity magnitude (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @return Horizontal velocity component (in m/s)
 * @throws std::invalid_argument if initialVelocity <= 0
 */
inline double calculateHorizontalVelocity(double initialVelocity, double angleRadians) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    return initialVelocity * std::cos(angleRadians);
}

/**
 * @brief Calculate vertical component of initial velocity
 *
 * v₀ᵧ = v₀ sin(θ)
 *
 * @param initialVelocity Initial velocity magnitude (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @return Vertical velocity component (in m/s)
 * @throws std::invalid_argument if initialVelocity <= 0
 */
inline double calculateVerticalVelocity(double initialVelocity, double angleRadians) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    return initialVelocity * std::sin(angleRadians);
}

// ============================================================================
// Time of Flight
// ============================================================================

/**
 * @brief Calculate total time of flight (landing at same height as launch)
 *
 * Time of flight: T = 2v₀ sin(θ) / g
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Total flight time (in seconds)
 * @throws std::invalid_argument if initialVelocity <= 0
 */
inline double calculateTimeOfFlight(double initialVelocity, double angleRadians,
                                    double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    double verticalVel = calculateVerticalVelocity(initialVelocity, angleRadians);
    return (2.0 * verticalVel) / gravity;
}

/**
 * @brief Calculate time to reach maximum height
 *
 * Time to max height: t_max = v₀ sin(θ) / g
 * (This is half the total flight time)
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Time to maximum height (in seconds)
 * @throws std::invalid_argument if initialVelocity <= 0
 */
inline double calculateTimeToMaxHeight(double initialVelocity, double angleRadians,
                                       double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    double verticalVel = calculateVerticalVelocity(initialVelocity, angleRadians);
    return verticalVel / gravity;
}

// ============================================================================
// Range and Height
// ============================================================================

/**
 * @brief Calculate horizontal range (landing at same height as launch)
 *
 * Range: R = v₀² sin(2θ) / g
 *
 * Maximum range occurs at 45° angle.
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Horizontal range (in meters)
 * @throws std::invalid_argument if initialVelocity <= 0
 */
inline double calculateRange(double initialVelocity, double angleRadians, double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    return (initialVelocity * initialVelocity * std::sin(2.0 * angleRadians)) / gravity;
}

/**
 * @brief Calculate maximum height reached
 *
 * Max height: H = (v₀ sin(θ))² / (2g)
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Maximum height (in meters)
 * @throws std::invalid_argument if initialVelocity <= 0
 */
inline double calculateMaxHeight(double initialVelocity, double angleRadians, double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    double verticalVel = calculateVerticalVelocity(initialVelocity, angleRadians);
    return (verticalVel * verticalVel) / (2.0 * gravity);
}

// ============================================================================
// Position at Time t
// ============================================================================

/**
 * @brief Calculate horizontal position at time t
 *
 * x(t) = v₀ cos(θ) × t
 * (No acceleration in horizontal direction)
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @param time Time since launch (in seconds, must be >= 0)
 * @return Horizontal position (in meters)
 * @throws std::invalid_argument if initialVelocity <= 0 or time < 0
 */
inline double calculateHorizontalPosition(double initialVelocity, double angleRadians, double time) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    double horizontalVel = calculateHorizontalVelocity(initialVelocity, angleRadians);
    return horizontalVel * time;
}

/**
 * @brief Calculate vertical position at time t
 *
 * y(t) = v₀ sin(θ) × t - (1/2)gt²
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @param time Time since launch (in seconds, must be >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Vertical position (in meters)
 * @throws std::invalid_argument if initialVelocity <= 0 or time < 0
 */
inline double calculateVerticalPosition(double initialVelocity, double angleRadians,
                                        double time, double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    double verticalVel = calculateVerticalVelocity(initialVelocity, angleRadians);
    return (verticalVel * time) - (0.5 * gravity * time * time);
}

// ============================================================================
// Velocity at Time t
// ============================================================================

/**
 * @brief Calculate horizontal velocity at time t
 *
 * vₓ(t) = v₀ cos(θ) (constant throughout flight)
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @return Horizontal velocity (in m/s)
 * @throws std::invalid_argument if initialVelocity <= 0
 */
inline double getHorizontalVelocityAtTime(double initialVelocity, double angleRadians) {
    return calculateHorizontalVelocity(initialVelocity, angleRadians);
}

/**
 * @brief Calculate vertical velocity at time t
 *
 * vᵧ(t) = v₀ sin(θ) - gt
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @param time Time since launch (in seconds, must be >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Vertical velocity (in m/s, positive = upward)
 * @throws std::invalid_argument if initialVelocity <= 0 or time < 0
 */
inline double getVerticalVelocityAtTime(double initialVelocity, double angleRadians,
                                        double time, double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    double verticalVel = calculateVerticalVelocity(initialVelocity, angleRadians);
    return verticalVel - (gravity * time);
}

/**
 * @brief Calculate speed (velocity magnitude) at time t
 *
 * speed = √(vₓ² + vᵧ²)
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @param time Time since launch (in seconds, must be >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Speed (in m/s)
 * @throws std::invalid_argument if initialVelocity <= 0 or time < 0
 */
inline double calculateSpeedAtTime(double initialVelocity, double angleRadians,
                                   double time, double gravity = 9.81) {
    double vx = getHorizontalVelocityAtTime(initialVelocity, angleRadians);
    double vy = getVerticalVelocityAtTime(initialVelocity, angleRadians, time, gravity);
    return std::sqrt(vx * vx + vy * vy);
}

// ============================================================================
// Trajectory Equation
// ============================================================================

/**
 * @brief Calculate height y at horizontal position x (trajectory equation)
 *
 * y = x tan(θ) - (gx²)/(2v₀²cos²(θ))
 *
 * This gives the parabolic trajectory of the projectile.
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param angleRadians Launch angle from horizontal (in radians)
 * @param horizontalDistance Horizontal position (in meters, must be >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Height at horizontal position (in meters)
 * @throws std::invalid_argument if initialVelocity <= 0 or horizontalDistance < 0
 */
inline double calculateHeightAtDistance(double initialVelocity, double angleRadians,
                                        double horizontalDistance, double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (horizontalDistance < 0) {
        throw std::invalid_argument("Horizontal distance must be non-negative");
    }

    double tanTheta = std::tan(angleRadians);
    double cosTheta = std::cos(angleRadians);
    double v0Squared = initialVelocity * initialVelocity;

    return (horizontalDistance * tanTheta) -
           ((gravity * horizontalDistance * horizontalDistance) / (2.0 * v0Squared * cosTheta * cosTheta));
}

// ============================================================================
// Special Cases and Analysis
// ============================================================================

/**
 * @brief Calculate launch angle for maximum range
 *
 * Maximum range occurs at 45° (π/4 radians) for level ground.
 *
 * @return Optimal angle for maximum range (in radians) = π/4
 */
inline double getAngleForMaxRange() {
    return M_PI / 4.0; // 45 degrees
}

/**
 * @brief Calculate maximum possible range for given initial velocity
 *
 * Max range = v₀² / g (at 45° angle)
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Maximum possible range (in meters)
 * @throws std::invalid_argument if initialVelocity <= 0
 */
inline double calculateMaximumRange(double initialVelocity, double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    return (initialVelocity * initialVelocity) / gravity;
}

/**
 * @brief Calculate two possible launch angles for a given range
 *
 * For a given range R < R_max, there are two angles that achieve it:
 * θ₁ = (1/2) arcsin(gR/v₀²)  (lower angle, flatter trajectory)
 * θ₂ = 90° - θ₁               (higher angle, higher trajectory)
 *
 * This function returns the lower angle.
 *
 * @param initialVelocity Initial velocity (in m/s, must be > 0)
 * @param range Desired range (in meters, must be > 0 and <= max range)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Lower launch angle (in radians)
 * @throws std::invalid_argument if parameters invalid or range too large
 */
inline double calculateLaunchAngleForRange(double initialVelocity, double range,
                                           double gravity = 9.81) {
    if (initialVelocity <= 0) {
        throw std::invalid_argument("Initial velocity must be positive");
    }
    if (range <= 0) {
        throw std::invalid_argument("Range must be positive");
    }

    double maxRange = calculateMaximumRange(initialVelocity, gravity);
    if (range > maxRange) {
        throw std::invalid_argument("Range exceeds maximum possible range");
    }

    double sinValue = (gravity * range) / (initialVelocity * initialVelocity);
    return 0.5 * std::asin(sinValue);
}

/**
 * @brief Calculate velocity required to hit target at given range and height
 *
 * For target at (R, H), minimum velocity is:
 * v₀ = √(g(R² + 4H²)/(2H + R sin(2θ)))
 *
 * Simplified for horizontal target (H = 0):
 * v₀ = √(gR/sin(2θ))
 *
 * @param range Horizontal distance to target (in meters, must be > 0)
 * @param angleRadians Launch angle (in radians)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Required initial velocity (in m/s)
 * @throws std::invalid_argument if range <= 0 or angle invalid
 */
inline double calculateRequiredVelocity(double range, double angleRadians, double gravity = 9.81) {
    if (range <= 0) {
        throw std::invalid_argument("Range must be positive");
    }

    double sin2Theta = std::sin(2.0 * angleRadians);
    if (std::abs(sin2Theta) < 1e-10) {
        throw std::invalid_argument("Invalid angle: sin(2θ) too close to zero");
    }

    return std::sqrt((gravity * range) / sin2Theta);
}

} // namespace projectile
} // namespace physics

#endif // PHYSICS_PROJECTILE_HPP
