#ifndef PHYSICS_ELECTROSTATICS_HPP
#define PHYSICS_ELECTROSTATICS_HPP

#include <cmath>
#include <stdexcept>
#include <vector>

/**
 * @file electrostatics.hpp
 * @brief Comprehensive implementation of electrostatics
 *
 * This module implements:
 * - Coulomb's law (electric force between charges)
 * - Electric field calculations
 * - Electric potential and potential energy
 * - Capacitance (isolated conductors, parallel plate, spherical, cylindrical)
 * - Energy stored in capacitors
 * - Equipotential surfaces
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace electrostatics {

/**
 * @namespace constants
 * @brief Physical constants for electrostatics
 */
namespace constants {
    // Coulomb's constant
    constexpr double K_E = 8.9875517923e9;                  // N·m²/C²

    // Permittivity of free space
    constexpr double EPSILON_0 = 8.8541878128e-12;          // F/m (C²/(N·m²))

    // Elementary charge
    constexpr double ELEMENTARY_CHARGE = 1.602176634e-19;    // C

    // Electron mass
    constexpr double ELECTRON_MASS = 9.1093837015e-31;       // kg
}

// ============================================================================
// COULOMB'S LAW
// ============================================================================

/**
 * @brief Calculate electrostatic force between two point charges (Coulomb's law)
 *
 * F = k|q₁q₂|/r²
 *
 * where k = 1/(4πε₀) ≈ 8.99×10⁹ N·m²/C²
 *
 * @param charge1 First charge (C)
 * @param charge2 Second charge (C)
 * @param distance Distance between charges (m)
 * @param k Coulomb's constant (N·m²/C²)
 * @return Magnitude of electrostatic force (N)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double coulombForce(double charge1, double charge2, double distance,
                           double k = constants::K_E) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return k * std::abs(charge1 * charge2) / (distance * distance);
}

/**
 * @brief Determine if force is attractive or repulsive
 *
 * @param charge1 First charge (C)
 * @param charge2 Second charge (C)
 * @return true if attractive (opposite signs), false if repulsive (same signs)
 */
inline bool isAttractive(double charge1, double charge2) {
    return (charge1 * charge2) < 0.0;
}

/**
 * @brief Calculate electric field magnitude from point charge
 *
 * E = k|q|/r²
 *
 * @param charge Source charge (C)
 * @param distance Distance from charge (m)
 * @param k Coulomb's constant (N·m²/C²)
 * @return Electric field magnitude (N/C or V/m)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double electricField(double charge, double distance, double k = constants::K_E) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return k * std::abs(charge) / (distance * distance);
}

/**
 * @brief Calculate force on test charge in electric field
 *
 * F = qE
 *
 * @param charge Test charge (C)
 * @param electricField Electric field strength (N/C)
 * @return Force on charge (N)
 */
inline double forceInField(double charge, double electricField) {
    return std::abs(charge * electricField);
}

// ============================================================================
// ELECTRIC POTENTIAL
// ============================================================================

/**
 * @brief Calculate electric potential from point charge
 *
 * V = kq/r
 *
 * Potential is a scalar (sign matters)
 *
 * @param charge Source charge (C)
 * @param distance Distance from charge (m)
 * @param k Coulomb's constant (N·m²/C²)
 * @return Electric potential (V)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double electricPotential(double charge, double distance, double k = constants::K_E) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return k * charge / distance;
}

/**
 * @brief Calculate potential energy between two point charges
 *
 * U = kq₁q₂/r
 *
 * @param charge1 First charge (C)
 * @param charge2 Second charge (C)
 * @param distance Distance between charges (m)
 * @param k Coulomb's constant (N·m²/C²)
 * @return Potential energy (J)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double electricPotentialEnergy(double charge1, double charge2, double distance,
                                      double k = constants::K_E) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return k * charge1 * charge2 / distance;
}

/**
 * @brief Calculate work done in moving charge in electric field
 *
 * W = q(V_b - V_a) = qΔV
 *
 * @param charge Charge being moved (C)
 * @param potentialDifference Potential difference (V)
 * @return Work done (J), positive = work done by field
 */
inline double workInElectricField(double charge, double potentialDifference) {
    return charge * potentialDifference;
}

/**
 * @brief Calculate potential at point due to multiple charges
 *
 * V_total = Σ(kq_i/r_i)
 *
 * Superposition principle for potential
 *
 * @param charges Vector of charges (C)
 * @param distances Vector of distances from point (m)
 * @param k Coulomb's constant (N·m²/C²)
 * @return Total electric potential (V)
 * @throws std::invalid_argument if vectors have different sizes or invalid distances
 */
inline double totalPotential(const std::vector<double>& charges,
                             const std::vector<double>& distances,
                             double k = constants::K_E) {
    if (charges.size() != distances.size()) {
        throw std::invalid_argument("Charges and distances vectors must have same size");
    }

    double total = 0.0;
    for (size_t i = 0; i < charges.size(); ++i) {
        if (distances[i] <= 0.0) {
            throw std::invalid_argument("All distances must be positive");
        }
        total += k * charges[i] / distances[i];
    }
    return total;
}

/**
 * @brief Calculate electric field from potential gradient
 *
 * E = -dV/dr
 *
 * Electric field points in direction of decreasing potential
 *
 * @param potentialDifference Potential difference (V)
 * @param distance Distance over which potential changes (m)
 * @return Electric field magnitude (V/m)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double fieldFromPotentialGradient(double potentialDifference, double distance) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return -potentialDifference / distance;
}

// ============================================================================
// CAPACITANCE
// ============================================================================

/**
 * @brief Calculate capacitance from charge and potential
 *
 * C = Q/V
 *
 * @param charge Charge stored (C)
 * @param potential Potential difference (V)
 * @return Capacitance (F)
 * @throws std::invalid_argument if potential is zero
 */
inline double capacitanceFromCharge(double charge, double potential) {
    if (std::abs(potential) < 1e-15) {
        throw std::invalid_argument("Potential must be non-zero");
    }
    return std::abs(charge / potential);
}

/**
 * @brief Calculate charge stored on capacitor
 *
 * Q = CV
 *
 * @param capacitance Capacitance (F)
 * @param voltage Voltage across capacitor (V)
 * @return Charge stored (C)
 * @throws std::invalid_argument if capacitance is non-positive
 */
inline double chargeOnCapacitor(double capacitance, double voltage) {
    if (capacitance <= 0.0) {
        throw std::invalid_argument("Capacitance must be positive");
    }
    return capacitance * std::abs(voltage);
}

/**
 * @brief Calculate capacitance of isolated conducting sphere
 *
 * C = 4πε₀R
 *
 * @param radius Radius of sphere (m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Capacitance (F)
 * @throws std::invalid_argument if radius is non-positive
 */
inline double sphereCapacitance(double radius, double epsilon0 = constants::EPSILON_0) {
    if (radius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return 4.0 * M_PI * epsilon0 * radius;
}

/**
 * @brief Calculate capacitance of parallel plate capacitor
 *
 * C = ε₀A/d
 *
 * @param area Area of plates (m²)
 * @param separation Distance between plates (m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Capacitance (F)
 * @throws std::invalid_argument if area or separation is non-positive
 */
inline double parallelPlateCapacitance(double area, double separation,
                                       double epsilon0 = constants::EPSILON_0) {
    if (area <= 0.0 || separation <= 0.0) {
        throw std::invalid_argument("Area and separation must be positive");
    }
    return epsilon0 * area / separation;
}

/**
 * @brief Calculate capacitance with dielectric material
 *
 * C = κε₀A/d
 *
 * where κ is dielectric constant (relative permittivity)
 *
 * @param area Area of plates (m²)
 * @param separation Distance between plates (m)
 * @param dielectricConstant Relative permittivity (dimensionless)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Capacitance (F)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double capacitanceWithDielectric(double area, double separation,
                                        double dielectricConstant,
                                        double epsilon0 = constants::EPSILON_0) {
    if (area <= 0.0 || separation <= 0.0 || dielectricConstant <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return dielectricConstant * epsilon0 * area / separation;
}

/**
 * @brief Calculate capacitance of two concentric spheres
 *
 * C = 4πε₀ × (R₁R₂)/(R₂ - R₁)
 *
 * where R₁ < R₂
 *
 * @param innerRadius Inner sphere radius (m)
 * @param outerRadius Outer sphere radius (m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Capacitance (F)
 * @throws std::invalid_argument if radii are invalid
 */
inline double concentricSpheresCapacitance(double innerRadius, double outerRadius,
                                           double epsilon0 = constants::EPSILON_0) {
    if (innerRadius <= 0.0 || outerRadius <= 0.0) {
        throw std::invalid_argument("Radii must be positive");
    }
    if (outerRadius <= innerRadius) {
        throw std::invalid_argument("Outer radius must be greater than inner radius");
    }
    return 4.0 * M_PI * epsilon0 * innerRadius * outerRadius / (outerRadius - innerRadius);
}

/**
 * @brief Calculate capacitance of cylindrical capacitor
 *
 * C = 2πε₀L/ln(R₂/R₁)
 *
 * @param length Length of cylinders (m)
 * @param innerRadius Inner cylinder radius (m)
 * @param outerRadius Outer cylinder radius (m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Capacitance (F)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double cylindricalCapacitance(double length, double innerRadius, double outerRadius,
                                     double epsilon0 = constants::EPSILON_0) {
    if (length <= 0.0 || innerRadius <= 0.0 || outerRadius <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    if (outerRadius <= innerRadius) {
        throw std::invalid_argument("Outer radius must be greater than inner radius");
    }
    return 2.0 * M_PI * epsilon0 * length / std::log(outerRadius / innerRadius);
}

// ============================================================================
// CAPACITOR COMBINATIONS
// ============================================================================

/**
 * @brief Calculate equivalent capacitance for series combination
 *
 * 1/C_eq = 1/C₁ + 1/C₂ + 1/C₃ + ...
 *
 * @param capacitances Vector of capacitances (F)
 * @return Equivalent capacitance (F)
 * @throws std::invalid_argument if any capacitance is non-positive
 */
inline double seriesCapacitance(const std::vector<double>& capacitances) {
    if (capacitances.empty()) {
        throw std::invalid_argument("Capacitances vector cannot be empty");
    }

    double reciprocalSum = 0.0;
    for (double c : capacitances) {
        if (c <= 0.0) {
            throw std::invalid_argument("All capacitances must be positive");
        }
        reciprocalSum += 1.0 / c;
    }
    return 1.0 / reciprocalSum;
}

/**
 * @brief Calculate equivalent capacitance for parallel combination
 *
 * C_eq = C₁ + C₂ + C₃ + ...
 *
 * @param capacitances Vector of capacitances (F)
 * @return Equivalent capacitance (F)
 * @throws std::invalid_argument if any capacitance is negative
 */
inline double parallelCapacitance(const std::vector<double>& capacitances) {
    if (capacitances.empty()) {
        throw std::invalid_argument("Capacitances vector cannot be empty");
    }

    double total = 0.0;
    for (double c : capacitances) {
        if (c < 0.0) {
            throw std::invalid_argument("Capacitances must be non-negative");
        }
        total += c;
    }
    return total;
}

// ============================================================================
// ENERGY IN CAPACITORS
// ============================================================================

/**
 * @brief Calculate energy stored in capacitor
 *
 * U = (1/2)CV² = (1/2)QV = Q²/(2C)
 *
 * @param capacitance Capacitance (F)
 * @param voltage Voltage across capacitor (V)
 * @return Energy stored (J)
 * @throws std::invalid_argument if capacitance is non-positive
 */
inline double capacitorEnergy(double capacitance, double voltage) {
    if (capacitance <= 0.0) {
        throw std::invalid_argument("Capacitance must be positive");
    }
    return 0.5 * capacitance * voltage * voltage;
}

/**
 * @brief Calculate energy stored from charge
 *
 * U = Q²/(2C)
 *
 * @param charge Charge stored (C)
 * @param capacitance Capacitance (F)
 * @return Energy stored (J)
 * @throws std::invalid_argument if capacitance is non-positive
 */
inline double capacitorEnergyFromCharge(double charge, double capacitance) {
    if (capacitance <= 0.0) {
        throw std::invalid_argument("Capacitance must be positive");
    }
    return (charge * charge) / (2.0 * capacitance);
}

/**
 * @brief Calculate energy density in electric field
 *
 * u = (1/2)ε₀E²
 *
 * Energy per unit volume
 *
 * @param electricField Electric field strength (V/m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Energy density (J/m³)
 */
inline double electricEnergyDensity(double electricField, double epsilon0 = constants::EPSILON_0) {
    return 0.5 * epsilon0 * electricField * electricField;
}

/**
 * @brief Calculate electric field between parallel plates
 *
 * E = V/d = σ/ε₀
 *
 * where σ is surface charge density
 *
 * @param voltage Voltage across plates (V)
 * @param separation Distance between plates (m)
 * @return Electric field strength (V/m)
 * @throws std::invalid_argument if separation is non-positive
 */
inline double fieldBetweenPlates(double voltage, double separation) {
    if (separation <= 0.0) {
        throw std::invalid_argument("Separation must be positive");
    }
    return std::abs(voltage) / separation;
}

/**
 * @brief Calculate surface charge density on capacitor plate
 *
 * σ = Q/A
 *
 * @param charge Total charge on plate (C)
 * @param area Area of plate (m²)
 * @return Surface charge density (C/m²)
 * @throws std::invalid_argument if area is non-positive
 */
inline double surfaceChargeDensity(double charge, double area) {
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    return std::abs(charge) / area;
}

// ============================================================================
// UNIT CHARGE AND CONVERSIONS
// ============================================================================

/**
 * @brief Calculate number of elementary charges
 *
 * n = Q/e
 *
 * @param charge Total charge (C)
 * @param elementaryCharge Elementary charge (C)
 * @return Number of elementary charges
 */
inline double numberOfCharges(double charge, double elementaryCharge = constants::ELEMENTARY_CHARGE) {
    return charge / elementaryCharge;
}

/**
 * @brief Calculate charge from number of electrons/protons
 *
 * Q = ne
 *
 * @param number Number of charges
 * @param elementaryCharge Elementary charge (C)
 * @return Total charge (C)
 */
inline double chargeFromNumber(double number, double elementaryCharge = constants::ELEMENTARY_CHARGE) {
    return number * elementaryCharge;
}

} // namespace electrostatics
} // namespace physics

#endif // PHYSICS_ELECTROSTATICS_HPP
