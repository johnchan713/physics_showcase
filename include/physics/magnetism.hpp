#ifndef PHYSICS_MAGNETISM_HPP
#define PHYSICS_MAGNETISM_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file magnetism.hpp
 * @brief Comprehensive implementation of magnetism and magnetic fields
 *
 * This module implements:
 * - Magnetic fields and forces
 * - Magnetic moments
 * - Unit magnetic pole
 * - Magnetic field intensity
 * - Force between magnets
 * - Gauss' method for measuring field intensity
 * - Magnet-armature attraction
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace magnetism {

/**
 * @namespace constants
 * @brief Physical constants for magnetism
 */
namespace constants {
    // Permeability of free space
    constexpr double MU_0 = 1.25663706212e-6;  // H/m or T·m/A or N/A²

    // Earth's magnetic field (approximate)
    constexpr double EARTH_MAGNETIC_FIELD = 5e-5;  // T (50 μT)

    // Electron charge
    constexpr double ELECTRON_CHARGE = 1.602176634e-19;  // C
}

// ============================================================================
// MAGNETIC FORCE ON MOVING CHARGES
// ============================================================================

/**
 * @brief Calculate magnetic force on moving charge (Lorentz force)
 *
 * F = qvB sin(θ)
 *
 * where θ is angle between velocity and magnetic field
 *
 * @param charge Charge (C)
 * @param velocity Speed of charge (m/s)
 * @param magneticField Magnetic field strength (T)
 * @param angle Angle between velocity and field (radians)
 * @return Magnetic force (N)
 */
inline double magneticForceOnCharge(double charge, double velocity,
                                    double magneticField, double angle) {
    return std::abs(charge) * velocity * magneticField * std::sin(angle);
}

/**
 * @brief Calculate magnetic force when perpendicular (maximum force)
 *
 * F = qvB (when v ⊥ B)
 *
 * @param charge Charge (C)
 * @param velocity Speed (m/s)
 * @param magneticField Magnetic field strength (T)
 * @return Maximum magnetic force (N)
 */
inline double maxMagneticForce(double charge, double velocity, double magneticField) {
    return std::abs(charge) * velocity * magneticField;
}

/**
 * @brief Calculate radius of circular motion in magnetic field
 *
 * r = mv/(qB)
 *
 * When charged particle moves perpendicular to uniform B field
 *
 * @param mass Mass of particle (kg)
 * @param velocity Speed (m/s)
 * @param charge Charge (C)
 * @param magneticField Magnetic field strength (T)
 * @return Radius of circular path (m)
 * @throws std::invalid_argument if charge or field is zero
 */
inline double cyclotronRadius(double mass, double velocity, double charge,
                              double magneticField) {
    if (std::abs(charge) < 1e-30 || std::abs(magneticField) < 1e-15) {
        throw std::invalid_argument("Charge and magnetic field must be non-zero");
    }
    return (mass * velocity) / (std::abs(charge) * magneticField);
}

/**
 * @brief Calculate cyclotron frequency
 *
 * f = qB/(2πm)
 *
 * Frequency of circular motion in magnetic field
 *
 * @param charge Charge (C)
 * @param magneticField Magnetic field strength (T)
 * @param mass Mass of particle (kg)
 * @return Frequency (Hz)
 * @throws std::invalid_argument if mass is non-positive
 */
inline double cyclotronFrequency(double charge, double magneticField, double mass) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    return (std::abs(charge) * magneticField) / (2.0 * M_PI * mass);
}

// ============================================================================
// MAGNETIC FORCE ON CURRENT-CARRYING WIRE
// ============================================================================

/**
 * @brief Calculate force on current-carrying wire in magnetic field
 *
 * F = BIL sin(θ)
 *
 * where θ is angle between current direction and field
 *
 * @param magneticField Magnetic field strength (T)
 * @param current Current in wire (A)
 * @param length Length of wire in field (m)
 * @param angle Angle between current and field (radians)
 * @return Magnetic force on wire (N)
 */
inline double forceOnWire(double magneticField, double current, double length,
                          double angle) {
    if (length < 0.0) {
        throw std::invalid_argument("Length must be non-negative");
    }
    return magneticField * current * length * std::sin(angle);
}

/**
 * @brief Calculate maximum force on wire (perpendicular case)
 *
 * F = BIL (when I ⊥ B)
 *
 * @param magneticField Magnetic field strength (T)
 * @param current Current (A)
 * @param length Length of wire (m)
 * @return Maximum force (N)
 */
inline double maxForceOnWire(double magneticField, double current, double length) {
    if (length < 0.0) {
        throw std::invalid_argument("Length must be non-negative");
    }
    return magneticField * current * length;
}

// ============================================================================
// MAGNETIC FIELD FROM CURRENT
// ============================================================================

/**
 * @brief Calculate magnetic field from long straight wire (Ampere's law)
 *
 * B = μ₀I/(2πr)
 *
 * @param current Current in wire (A)
 * @param distance Distance from wire (m)
 * @param mu0 Permeability of free space (H/m)
 * @return Magnetic field strength (T)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double fieldFromStraightWire(double current, double distance,
                                    double mu0 = constants::MU_0) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return (mu0 * current) / (2.0 * M_PI * distance);
}

/**
 * @brief Calculate magnetic field at center of circular loop
 *
 * B = μ₀I/(2R)
 *
 * @param current Current in loop (A)
 * @param radius Radius of loop (m)
 * @param mu0 Permeability of free space (H/m)
 * @return Magnetic field at center (T)
 * @throws std::invalid_argument if radius is non-positive
 */
inline double fieldAtCenterOfLoop(double current, double radius,
                                  double mu0 = constants::MU_0) {
    if (radius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return (mu0 * current) / (2.0 * radius);
}

/**
 * @brief Calculate magnetic field inside solenoid
 *
 * B = μ₀nI
 *
 * where n is number of turns per unit length
 *
 * @param turnsPerLength Number of turns per meter (turns/m)
 * @param current Current (A)
 * @param mu0 Permeability of free space (H/m)
 * @return Magnetic field inside solenoid (T)
 */
inline double fieldInsideSolenoid(double turnsPerLength, double current,
                                  double mu0 = constants::MU_0) {
    if (turnsPerLength < 0.0) {
        throw std::invalid_argument("Turns per length must be non-negative");
    }
    return mu0 * turnsPerLength * current;
}

/**
 * @brief Calculate turns per length for solenoid
 *
 * n = N/L
 *
 * @param totalTurns Total number of turns
 * @param length Length of solenoid (m)
 * @return Turns per unit length (turns/m)
 * @throws std::invalid_argument if length is non-positive
 */
inline double calculateTurnsPerLength(int totalTurns, double length) {
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    return totalTurns / length;
}

// ============================================================================
// MAGNETIC MOMENT
// ============================================================================

/**
 * @brief Calculate magnetic moment of current loop
 *
 * μ = IA
 *
 * where A is area of loop
 *
 * @param current Current in loop (A)
 * @param area Area enclosed by loop (m²)
 * @return Magnetic moment (A·m²)
 */
inline double magneticMoment(double current, double area) {
    if (area < 0.0) {
        throw std::invalid_argument("Area must be non-negative");
    }
    return current * area;
}

/**
 * @brief Calculate magnetic moment of circular loop
 *
 * μ = I × πR²
 *
 * @param current Current (A)
 * @param radius Radius of loop (m)
 * @return Magnetic moment (A·m²)
 * @throws std::invalid_argument if radius is negative
 */
inline double magneticMomentCircularLoop(double current, double radius) {
    if (radius < 0.0) {
        throw std::invalid_argument("Radius must be non-negative");
    }
    return current * M_PI * radius * radius;
}

/**
 * @brief Calculate torque on magnetic dipole in field
 *
 * τ = μB sin(θ)
 *
 * where θ is angle between moment and field
 *
 * @param magneticMoment Magnetic moment (A·m²)
 * @param magneticField Magnetic field strength (T)
 * @param angle Angle between moment and field (radians)
 * @return Torque (N·m)
 */
inline double torqueOnDipole(double magneticMoment, double magneticField, double angle) {
    return magneticMoment * magneticField * std::sin(angle);
}

/**
 * @brief Calculate potential energy of magnetic dipole
 *
 * U = -μB cos(θ)
 *
 * @param magneticMoment Magnetic moment (A·m²)
 * @param magneticField Magnetic field strength (T)
 * @param angle Angle between moment and field (radians)
 * @return Potential energy (J)
 */
inline double dipoleEnergy(double magneticMoment, double magneticField, double angle) {
    return -magneticMoment * magneticField * std::cos(angle);
}

// ============================================================================
// MAGNETIC POLES AND BAR MAGNETS
// ============================================================================

/**
 * @brief Calculate force between two magnetic poles (analogous to Coulomb)
 *
 * F = (μ₀/(4π)) × (m₁m₂/r²)
 *
 * where m₁, m₂ are pole strengths in A·m
 *
 * @param poleStrength1 First pole strength (A·m)
 * @param poleStrength2 Second pole strength (A·m)
 * @param distance Distance between poles (m)
 * @param mu0 Permeability of free space (H/m)
 * @return Force between poles (N)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double forceBetweenPoles(double poleStrength1, double poleStrength2,
                                double distance, double mu0 = constants::MU_0) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return (mu0 / (4.0 * M_PI)) * std::abs(poleStrength1 * poleStrength2) /
           (distance * distance);
}

/**
 * @brief Calculate magnetic field strength from magnetic pole
 *
 * H = m/(4πr²)
 *
 * H is magnetic field intensity (A/m)
 *
 * @param poleStrength Pole strength (A·m)
 * @param distance Distance from pole (m)
 * @return Magnetic field intensity (A/m)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double fieldIntensityFromPole(double poleStrength, double distance) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return std::abs(poleStrength) / (4.0 * M_PI * distance * distance);
}

/**
 * @brief Calculate magnetic field B from field intensity H
 *
 * B = μ₀H (in vacuum/air)
 *
 * @param fieldIntensity Magnetic field intensity H (A/m)
 * @param mu0 Permeability of free space (H/m)
 * @return Magnetic field B (T)
 */
inline double magneticFieldFromIntensity(double fieldIntensity,
                                         double mu0 = constants::MU_0) {
    return mu0 * fieldIntensity;
}

/**
 * @brief Calculate field at point near bar magnet (on axis)
 *
 * For short magnet: B ≈ (μ₀/(4π)) × (2M/r³)
 *
 * where M is magnetic moment
 *
 * @param magneticMoment Magnetic moment of magnet (A·m²)
 * @param distance Distance along axis from center (m)
 * @param mu0 Permeability of free space (H/m)
 * @return Magnetic field on axis (T)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double fieldOnAxisOfMagnet(double magneticMoment, double distance,
                                  double mu0 = constants::MU_0) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return (mu0 / (4.0 * M_PI)) * (2.0 * magneticMoment) /
           (distance * distance * distance);
}

/**
 * @brief Calculate field at equatorial point of bar magnet
 *
 * For short magnet: B ≈ (μ₀/(4π)) × (M/r³)
 *
 * @param magneticMoment Magnetic moment of magnet (A·m²)
 * @param distance Distance perpendicular to axis (m)
 * @param mu0 Permeability of free space (H/m)
 * @return Magnetic field at equator (T)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double fieldAtEquatorOfMagnet(double magneticMoment, double distance,
                                     double mu0 = constants::MU_0) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return (mu0 / (4.0 * M_PI)) * magneticMoment /
           (distance * distance * distance);
}

// ============================================================================
// FORCE WITH WHICH MAGNET ATTRACTS ARMATURE
// ============================================================================

/**
 * @brief Calculate attractive force between magnet and ferromagnetic armature
 *
 * F = B²A/(2μ₀)
 *
 * where B is field at pole face, A is contact area
 *
 * @param magneticField Magnetic field at pole face (T)
 * @param area Contact area (m²)
 * @param mu0 Permeability of free space (H/m)
 * @return Attractive force (N)
 * @throws std::invalid_argument if area is negative
 */
inline double magnetArmatureForce(double magneticField, double area,
                                  double mu0 = constants::MU_0) {
    if (area < 0.0) {
        throw std::invalid_argument("Area must be non-negative");
    }
    return (magneticField * magneticField * area) / (2.0 * mu0);
}

/**
 * @brief Calculate magnetic pressure at interface
 *
 * P = B²/(2μ₀)
 *
 * Pressure (force per unit area) due to magnetic field
 *
 * @param magneticField Magnetic field strength (T)
 * @param mu0 Permeability of free space (H/m)
 * @return Magnetic pressure (Pa)
 */
inline double magneticPressure(double magneticField, double mu0 = constants::MU_0) {
    return (magneticField * magneticField) / (2.0 * mu0);
}

// ============================================================================
// GAUSS' METHOD FOR MEASURING HORIZONTAL INTENSITY
// ============================================================================

/**
 * @brief Calculate horizontal component of Earth's magnetic field (Gauss method)
 *
 * Using oscillation method: B_H = 4π²I/(MT²)
 *
 * where I is moment of inertia, M is magnetic moment, T is period
 *
 * @param momentOfInertia Moment of inertia of magnet (kg·m²)
 * @param magneticMoment Magnetic moment of magnet (A·m²)
 * @param period Period of oscillation (s)
 * @return Horizontal component of Earth's field (T)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double horizontalFieldGaussMethod(double momentOfInertia,
                                         double magneticMoment, double period) {
    if (momentOfInertia <= 0.0 || magneticMoment <= 0.0 || period <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return (4.0 * M_PI * M_PI * momentOfInertia) /
           (magneticMoment * period * period);
}

/**
 * @brief Calculate period of oscillation of magnet in Earth's field
 *
 * T = 2π√(I/(MB_H))
 *
 * @param momentOfInertia Moment of inertia (kg·m²)
 * @param magneticMoment Magnetic moment (A·m²)
 * @param horizontalField Horizontal component of field (T)
 * @return Period of oscillation (s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double magnetOscillationPeriod(double momentOfInertia, double magneticMoment,
                                      double horizontalField) {
    if (momentOfInertia <= 0.0 || magneticMoment <= 0.0 || horizontalField <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return 2.0 * M_PI * std::sqrt(momentOfInertia / (magneticMoment * horizontalField));
}

// ============================================================================
// MAGNETIC FLUX
// ============================================================================

/**
 * @brief Calculate magnetic flux through surface
 *
 * Φ = BA cos(θ)
 *
 * where θ is angle between field and surface normal
 *
 * @param magneticField Magnetic field strength (T)
 * @param area Area of surface (m²)
 * @param angle Angle between field and normal (radians)
 * @return Magnetic flux (Wb = T·m²)
 */
inline double magneticFlux(double magneticField, double area, double angle) {
    if (area < 0.0) {
        throw std::invalid_argument("Area must be non-negative");
    }
    return magneticField * area * std::cos(angle);
}

/**
 * @brief Calculate maximum flux (perpendicular case)
 *
 * Φ_max = BA (when B ⊥ surface)
 *
 * @param magneticField Magnetic field strength (T)
 * @param area Area (m²)
 * @return Maximum flux (Wb)
 */
inline double maxMagneticFlux(double magneticField, double area) {
    if (area < 0.0) {
        throw std::invalid_argument("Area must be non-negative");
    }
    return magneticField * area;
}

} // namespace magnetism
} // namespace physics

#endif // PHYSICS_MAGNETISM_HPP
