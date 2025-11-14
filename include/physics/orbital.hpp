#ifndef PHYSICS_ORBITAL_HPP
#define PHYSICS_ORBITAL_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace orbital {

/**
 * @brief Orbital Mechanics - Motion Around the Earth
 *
 * Functions for analyzing circular orbits around Earth (or other bodies).
 *
 * Key concepts:
 * - Gravitational force provides centripetal force: GMm/r² = mv²/r
 * - Orbital velocity: v = √(GM/r)
 * - Orbital period: T = 2π√(r³/GM)
 * - Escape velocity: v_esc = √(2GM/r)
 *
 * Constants for Earth:
 * - Mass: M_earth ≈ 5.972 × 10^24 kg
 * - Radius: R_earth ≈ 6.371 × 10^6 m
 * - GM_earth ≈ 3.986 × 10^14 m³/s²
 *
 * Standard notation:
 * - G: gravitational constant (6.674 × 10^-11 N⋅m²/kg²)
 * - M: mass of central body (e.g., Earth)
 * - m: mass of orbiting object
 * - r: orbital radius (from center of central body)
 * - h: altitude above surface
 */

// ============================================================================
// Physical Constants
// ============================================================================

namespace constants {
    constexpr double G = 6.674e-11;              // Gravitational constant (N⋅m²/kg²)
    constexpr double EARTH_MASS = 5.972e24;      // kg
    constexpr double EARTH_RADIUS = 6.371e6;     // meters
    constexpr double EARTH_GM = 3.986e14;        // m³/s²
    constexpr double MOON_ORBITAL_RADIUS = 3.844e8; // meters (Earth-Moon distance)
}

// ============================================================================
// Gravitational Force
// ============================================================================

/**
 * @brief Calculate gravitational force between two masses
 *
 * F = GMm/r²
 *
 * @param mass1 First mass (in kg, must be > 0)
 * @param mass2 Second mass (in kg, must be > 0)
 * @param distance Distance between centers (in meters, must be > 0)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Gravitational force (in Newtons)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateGravitationalForce(double mass1, double mass2, double distance,
                                          double G = constants::G) {
    if (mass1 <= 0 || mass2 <= 0) {
        throw std::invalid_argument("Masses must be positive");
    }
    if (distance <= 0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return (G * mass1 * mass2) / (distance * distance);
}

/**
 * @brief Calculate gravitational acceleration at distance from Earth's center
 *
 * g = GM/r²
 *
 * At Earth's surface (r = R_earth): g ≈ 9.81 m/s²
 *
 * @param distance Distance from Earth's center (in meters, must be > 0)
 * @param GM Product of G and Earth's mass (default: 3.986e14 m³/s²)
 * @return Gravitational acceleration (in m/s²)
 * @throws std::invalid_argument if distance <= 0
 */
inline double calculateGravitationalAcceleration(double distance, double GM = constants::EARTH_GM) {
    if (distance <= 0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return GM / (distance * distance);
}

// ============================================================================
// Orbital Velocity
// ============================================================================

/**
 * @brief Calculate orbital velocity for circular orbit
 *
 * For circular orbit: Gravitational force = Centripetal force
 * GMm/r² = mv²/r
 * Therefore: v = √(GM/r)
 *
 * Note: Orbital velocity is independent of the orbiting object's mass!
 *
 * @param orbitalRadius Distance from center of Earth (in meters, must be > 0)
 * @param GM Product of G and central body mass (default: Earth's GM)
 * @return Orbital velocity (in m/s)
 * @throws std::invalid_argument if orbitalRadius <= 0
 */
inline double calculateOrbitalVelocity(double orbitalRadius, double GM = constants::EARTH_GM) {
    if (orbitalRadius <= 0) {
        throw std::invalid_argument("Orbital radius must be positive");
    }
    return std::sqrt(GM / orbitalRadius);
}

/**
 * @brief Calculate orbital velocity from altitude above surface
 *
 * Orbital radius r = R_earth + altitude
 * Then: v = √(GM/r)
 *
 * @param altitude Height above Earth's surface (in meters, must be >= 0)
 * @param earthRadius Radius of Earth (default: 6.371e6 m)
 * @param GM Product of G and Earth's mass (default: 3.986e14 m³/s²)
 * @return Orbital velocity (in m/s)
 * @throws std::invalid_argument if altitude < 0
 */
inline double calculateOrbitalVelocityFromAltitude(double altitude,
                                                   double earthRadius = constants::EARTH_RADIUS,
                                                   double GM = constants::EARTH_GM) {
    if (altitude < 0) {
        throw std::invalid_argument("Altitude must be non-negative");
    }
    double orbitalRadius = earthRadius + altitude;
    return calculateOrbitalVelocity(orbitalRadius, GM);
}

// ============================================================================
// Orbital Period
// ============================================================================

/**
 * @brief Calculate orbital period (Kepler's Third Law)
 *
 * For circular orbit: T = 2πr/v = 2π√(r³/GM)
 *
 * This is Kepler's Third Law: T² ∝ r³
 *
 * @param orbitalRadius Distance from center of Earth (in meters, must be > 0)
 * @param GM Product of G and central body mass (default: Earth's GM)
 * @return Orbital period (in seconds)
 * @throws std::invalid_argument if orbitalRadius <= 0
 */
inline double calculateOrbitalPeriod(double orbitalRadius, double GM = constants::EARTH_GM) {
    if (orbitalRadius <= 0) {
        throw std::invalid_argument("Orbital radius must be positive");
    }
    return 2.0 * M_PI * std::sqrt((orbitalRadius * orbitalRadius * orbitalRadius) / GM);
}

/**
 * @brief Calculate orbital period from altitude above surface
 *
 * @param altitude Height above Earth's surface (in meters, must be >= 0)
 * @param earthRadius Radius of Earth (default: 6.371e6 m)
 * @param GM Product of G and Earth's mass (default: 3.986e14 m³/s²)
 * @return Orbital period (in seconds)
 * @throws std::invalid_argument if altitude < 0
 */
inline double calculateOrbitalPeriodFromAltitude(double altitude,
                                                 double earthRadius = constants::EARTH_RADIUS,
                                                 double GM = constants::EARTH_GM) {
    if (altitude < 0) {
        throw std::invalid_argument("Altitude must be non-negative");
    }
    double orbitalRadius = earthRadius + altitude;
    return calculateOrbitalPeriod(orbitalRadius, GM);
}

/**
 * @brief Calculate orbital radius for geostationary orbit
 *
 * Geostationary orbit: Period = 24 hours = 86400 seconds
 * From T² = (4π²/GM)r³, we get r = ∛(GMT²/(4π²))
 *
 * For Earth: r ≈ 42,164 km from center (≈ 35,786 km altitude)
 *
 * @param GM Product of G and Earth's mass (default: 3.986e14 m³/s²)
 * @return Orbital radius for geostationary orbit (in meters)
 */
inline double calculateGeostationaryOrbitRadius(double GM = constants::EARTH_GM) {
    double period = 86400.0; // 24 hours in seconds
    double periodSquared = period * period;
    return std::cbrt((GM * periodSquared) / (4.0 * M_PI * M_PI));
}

/**
 * @brief Calculate orbital radius for given period
 *
 * From Kepler's Third Law: r = ∛(GMT²/(4π²))
 *
 * @param period Desired orbital period (in seconds, must be > 0)
 * @param GM Product of G and central body mass (default: Earth's GM)
 * @return Orbital radius (in meters)
 * @throws std::invalid_argument if period <= 0
 */
inline double calculateOrbitalRadiusFromPeriod(double period, double GM = constants::EARTH_GM) {
    if (period <= 0) {
        throw std::invalid_argument("Period must be positive");
    }
    return std::cbrt((GM * period * period) / (4.0 * M_PI * M_PI));
}

// ============================================================================
// Escape Velocity
// ============================================================================

/**
 * @brief Calculate escape velocity from a distance
 *
 * Escape velocity is the minimum velocity needed to escape gravitational pull.
 * v_esc = √(2GM/r)
 *
 * Note: Escape velocity = √2 × orbital velocity
 *
 * @param distance Distance from center of body (in meters, must be > 0)
 * @param GM Product of G and central body mass (default: Earth's GM)
 * @return Escape velocity (in m/s)
 * @throws std::invalid_argument if distance <= 0
 */
inline double calculateEscapeVelocity(double distance, double GM = constants::EARTH_GM) {
    if (distance <= 0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return std::sqrt(2.0 * GM / distance);
}

/**
 * @brief Calculate escape velocity from Earth's surface
 *
 * @param earthRadius Radius of Earth (default: 6.371e6 m)
 * @param GM Product of G and Earth's mass (default: 3.986e14 m³/s²)
 * @return Escape velocity from surface (in m/s, ≈ 11,186 m/s)
 */
inline double calculateEscapeVelocityFromSurface(double earthRadius = constants::EARTH_RADIUS,
                                                 double GM = constants::EARTH_GM) {
    return calculateEscapeVelocity(earthRadius, GM);
}

/**
 * @brief Verify relationship: escape velocity = √2 × orbital velocity
 *
 * @param distance Distance from center (in meters, must be > 0)
 * @param GM Product of G and central body mass (default: Earth's GM)
 * @return Ratio of escape to orbital velocity (should be √2 ≈ 1.414)
 * @throws std::invalid_argument if distance <= 0
 */
inline double escapeToOrbitalVelocityRatio(double distance, double GM = constants::EARTH_GM) {
    if (distance <= 0) {
        throw std::invalid_argument("Distance must be positive");
    }
    double vOrbital = calculateOrbitalVelocity(distance, GM);
    double vEscape = calculateEscapeVelocity(distance, GM);
    return vEscape / vOrbital; // Should return √2
}

// ============================================================================
// Orbital Energy
// ============================================================================

/**
 * @brief Calculate gravitational potential energy in orbit
 *
 * PE = -GMm/r (negative because it's bound to Earth)
 *
 * @param mass Mass of orbiting object (in kg, must be > 0)
 * @param distance Distance from center of Earth (in meters, must be > 0)
 * @param GM Product of G and Earth's mass (default: 3.986e14 m³/s²)
 * @return Gravitational potential energy (in Joules, negative)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateOrbitalPotentialEnergy(double mass, double distance,
                                              double GM = constants::EARTH_GM) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (distance <= 0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return -(GM * mass) / distance;
}

/**
 * @brief Calculate kinetic energy in orbit
 *
 * KE = (1/2)mv² = (1/2)m(GM/r) = GMm/(2r)
 *
 * @param mass Mass of orbiting object (in kg, must be > 0)
 * @param orbitalRadius Distance from center (in meters, must be > 0)
 * @param GM Product of G and central body mass (default: Earth's GM)
 * @return Kinetic energy (in Joules)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateOrbitalKineticEnergy(double mass, double orbitalRadius,
                                            double GM = constants::EARTH_GM) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (orbitalRadius <= 0) {
        throw std::invalid_argument("Orbital radius must be positive");
    }
    return (GM * mass) / (2.0 * orbitalRadius);
}

/**
 * @brief Calculate total mechanical energy in orbit
 *
 * E_total = KE + PE = GMm/(2r) - GMm/r = -GMm/(2r)
 *
 * Total energy is negative (bound orbit) and equals -KE
 *
 * @param mass Mass of orbiting object (in kg, must be > 0)
 * @param orbitalRadius Distance from center (in meters, must be > 0)
 * @param GM Product of G and central body mass (default: Earth's GM)
 * @return Total mechanical energy (in Joules, negative for bound orbit)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateTotalOrbitalEnergy(double mass, double orbitalRadius,
                                         double GM = constants::EARTH_GM) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (orbitalRadius <= 0) {
        throw std::invalid_argument("Orbital radius must be positive");
    }
    return -(GM * mass) / (2.0 * orbitalRadius);
}

// ============================================================================
// Special Orbits
// ============================================================================

/**
 * @brief Calculate parameters for Low Earth Orbit (LEO)
 *
 * LEO typically at altitude 200-2000 km above Earth's surface
 * Returns orbital velocity for given LEO altitude
 *
 * @param altitude Altitude above Earth's surface (in meters, typically 200-2000 km)
 * @param earthRadius Radius of Earth (default: 6.371e6 m)
 * @param GM Product of G and Earth's mass (default: 3.986e14 m³/s²)
 * @return Orbital velocity (in m/s)
 */
inline double calculateLEOVelocity(double altitude,
                                   double earthRadius = constants::EARTH_RADIUS,
                                   double GM = constants::EARTH_GM) {
    return calculateOrbitalVelocityFromAltitude(altitude, earthRadius, GM);
}

/**
 * @brief Calculate centripetal acceleration in orbit
 *
 * a_c = GM/r² = v²/r
 *
 * This is also the gravitational acceleration at that altitude.
 *
 * @param orbitalRadius Distance from center (in meters, must be > 0)
 * @param GM Product of G and central body mass (default: Earth's GM)
 * @return Centripetal acceleration (in m/s²)
 * @throws std::invalid_argument if orbitalRadius <= 0
 */
inline double calculateOrbitalAcceleration(double orbitalRadius, double GM = constants::EARTH_GM) {
    return calculateGravitationalAcceleration(orbitalRadius, GM);
}

/**
 * @brief Calculate weight in orbit (apparent weight)
 *
 * In orbit, apparent weight is zero (weightlessness) even though
 * gravitational force still acts. This is because the satellite
 * and everything in it are in free fall.
 *
 * This function returns the actual gravitational force (not apparent weight).
 *
 * @param mass Mass of object (in kg, must be > 0)
 * @param orbitalRadius Distance from Earth's center (in meters, must be > 0)
 * @param GM Product of G and Earth's mass (default: 3.986e14 m³/s²)
 * @return Gravitational force (in Newtons)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateWeightInOrbit(double mass, double orbitalRadius,
                                     double GM = constants::EARTH_GM) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (orbitalRadius <= 0) {
        throw std::invalid_argument("Orbital radius must be positive");
    }
    double g = calculateGravitationalAcceleration(orbitalRadius, GM);
    return mass * g;
}

} // namespace orbital
} // namespace physics

#endif // PHYSICS_ORBITAL_HPP
