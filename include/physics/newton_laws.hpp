#ifndef PHYSICS_NEWTON_LAWS_HPP
#define PHYSICS_NEWTON_LAWS_HPP

#include <cmath>
#include <vector>

namespace physics {
namespace newton {

/**
 * @brief Newton's First Law of Motion (Law of Inertia)
 *
 * An object at rest stays at rest and an object in motion stays in motion
 * with the same speed and in the same direction unless acted upon by an
 * unbalanced force.
 *
 * This function checks if an object is in equilibrium (net force = 0)
 *
 * @param forces Vector of forces acting on the object (in Newtons)
 * @param tolerance Acceptable error margin for equilibrium check (default: 1e-6)
 * @return true if the object is in equilibrium, false otherwise
 */
inline bool isInEquilibrium(const std::vector<double>& forces, double tolerance = 1e-6) {
    double netForce = 0.0;
    for (const auto& force : forces) {
        netForce += force;
    }
    return std::abs(netForce) < tolerance;
}

/**
 * @brief Calculate net force acting on an object (First Law application)
 *
 * Sums all forces acting on an object to determine the net force.
 * If net force is zero, the object is in equilibrium.
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
 * @brief Newton's Second Law of Motion: F = ma
 *
 * The acceleration of an object is directly proportional to the net force
 * acting on it and inversely proportional to its mass.
 *
 * Calculate force given mass and acceleration
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param acceleration Acceleration of the object (in m/s²)
 * @return Force acting on the object (in Newtons)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateForce(double mass, double acceleration) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return mass * acceleration;
}

/**
 * @brief Newton's Second Law: Calculate acceleration from force and mass
 *
 * Derives acceleration using F = ma, therefore a = F/m
 *
 * @param force Net force acting on the object (in Newtons)
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @return Acceleration of the object (in m/s²)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateAcceleration(double force, double mass) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return force / mass;
}

/**
 * @brief Newton's Second Law: Calculate mass from force and acceleration
 *
 * Derives mass using F = ma, therefore m = F/a
 *
 * @param force Net force acting on the object (in Newtons)
 * @param acceleration Acceleration of the object (in m/s², must be != 0)
 * @return Mass of the object (in kilograms)
 * @throws std::invalid_argument if acceleration == 0
 */
inline double calculateMass(double force, double acceleration) {
    if (std::abs(acceleration) < 1e-10) {
        throw std::invalid_argument("Acceleration must be non-zero");
    }
    return force / acceleration;
}

/**
 * @brief Newton's Third Law of Motion: Action-Reaction Pairs
 *
 * For every action, there is an equal and opposite reaction.
 * When one object exerts a force on another, the second object exerts
 * an equal force in the opposite direction on the first object.
 *
 * This function calculates the reaction force given an action force
 *
 * @param actionForce The force exerted by object A on object B (in Newtons)
 * @return Reaction force exerted by object B on object A (in Newtons, opposite direction)
 */
inline double calculateReactionForce(double actionForce) {
    return -actionForce;
}

/**
 * @brief Verify Newton's Third Law for a pair of forces
 *
 * Checks if two forces form a valid action-reaction pair
 * (equal in magnitude, opposite in direction)
 *
 * @param force1 First force (in Newtons)
 * @param force2 Second force (in Newtons)
 * @param tolerance Acceptable error margin (default: 1e-6)
 * @return true if forces form an action-reaction pair, false otherwise
 */
inline bool verifyActionReactionPair(double force1, double force2, double tolerance = 1e-6) {
    return std::abs(force1 + force2) < tolerance;
}

/**
 * @brief Calculate weight force (application of Second Law)
 *
 * Weight is the gravitational force acting on an object: W = mg
 * where g is the acceleration due to gravity
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param gravity Acceleration due to gravity (in m/s², default: 9.81 m/s² on Earth)
 * @return Weight force (in Newtons)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateWeight(double mass, double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be greater than zero");
    }
    return mass * gravity;
}

} // namespace newton
} // namespace physics

#endif // PHYSICS_NEWTON_LAWS_HPP
