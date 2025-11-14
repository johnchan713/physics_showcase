#ifndef PHYSICS_MAXWELL_EQUATIONS_HPP
#define PHYSICS_MAXWELL_EQUATIONS_HPP

#include <cmath>
#include <stdexcept>
#include <array>

/**
 * @file maxwell_equations.hpp
 * @brief Maxwell's equations and electromagnetic field theory
 *
 * This module implements:
 * - Maxwell's equations in various forms
 * - Electromagnetic field energy and momentum
 * - Poynting vector and energy flow
 * - Electromagnetic potentials (scalar and vector)
 * - Gauge transformations
 * - Electromagnetic stress tensor
 * - Radiation from accelerating charges
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace maxwell {

/**
 * @namespace constants
 * @brief Physical constants for electromagnetism
 */
namespace constants {
    constexpr double EPSILON_0 = 8.854187817e-12;   // Permittivity of free space (F/m)
    constexpr double MU_0 = 4.0 * M_PI * 1e-7;      // Permeability of free space (H/m)
    constexpr double SPEED_OF_LIGHT = 299792458.0;  // Speed of light (m/s)
}

// ============================================================================
// ELECTRIC FIELD AND GAUSS'S LAW
// ============================================================================

/**
 * @brief Calculate electric field from charge density (Gauss's law differential form)
 *
 * ∇·E = ρ/ε₀
 *
 * For spherical symmetry: E = (Q/(4πε₀r²)) r̂
 *
 * @param charge Total enclosed charge Q (C)
 * @param distance Distance from center r (m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Electric field magnitude (V/m)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double electricFieldFromCharge(double charge, double distance,
                                      double epsilon0 = constants::EPSILON_0) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return charge / (4.0 * M_PI * epsilon0 * distance * distance);
}

/**
 * @brief Calculate electric flux through closed surface (Gauss's law integral form)
 *
 * Φ_E = ∮ E·dA = Q_enc/ε₀
 *
 * @param enclosedCharge Charge enclosed by surface (C)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Electric flux (V⋅m)
 */
inline double electricFlux(double enclosedCharge, double epsilon0 = constants::EPSILON_0) {
    return enclosedCharge / epsilon0;
}

/**
 * @brief Calculate electric field between parallel plates
 *
 * E = σ/ε₀ = Q/(Aε₀)
 *
 * @param surfaceChargeDensity Surface charge density σ (C/m²)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Electric field (V/m)
 */
inline double fieldBetweenPlates(double surfaceChargeDensity,
                                 double epsilon0 = constants::EPSILON_0) {
    return surfaceChargeDensity / epsilon0;
}

// ============================================================================
// MAGNETIC FIELD AND AMPERE'S LAW
// ============================================================================

/**
 * @brief Calculate magnetic field from current (Ampere's law for long straight wire)
 *
 * ∇×B = μ₀J + μ₀ε₀(∂E/∂t)
 *
 * For steady currents in wire: B = μ₀I/(2πr)
 *
 * @param current Current I (A)
 * @param distance Distance from wire r (m)
 * @param mu0 Permeability of free space (H/m)
 * @return Magnetic field magnitude (T)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double magneticFieldFromCurrent(double current, double distance,
                                       double mu0 = constants::MU_0) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return mu0 * current / (2.0 * M_PI * distance);
}

/**
 * @brief Calculate displacement current
 *
 * I_d = ε₀(dΦ_E/dt)
 *
 * Maxwell's correction to Ampere's law
 *
 * @param electricFluxRate Rate of change of electric flux (V⋅m/s)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Displacement current (A)
 */
inline double displacementCurrent(double electricFluxRate,
                                  double epsilon0 = constants::EPSILON_0) {
    return epsilon0 * electricFluxRate;
}

/**
 * @brief Calculate magnetic field in solenoid
 *
 * B = μ₀nI
 *
 * @param turnsPerLength Number of turns per unit length n (turns/m)
 * @param current Current I (A)
 * @param mu0 Permeability of free space (H/m)
 * @return Magnetic field inside solenoid (T)
 */
inline double solenoidField(double turnsPerLength, double current,
                           double mu0 = constants::MU_0) {
    return mu0 * turnsPerLength * current;
}

// ============================================================================
// FARADAY'S LAW
// ============================================================================

/**
 * @brief Calculate induced EMF from changing magnetic flux (Faraday's law)
 *
 * ε = -dΦ_B/dt = -d(BA)/dt
 *
 * @param magneticFluxRate Rate of change of magnetic flux (Wb/s = V)
 * @return Induced EMF magnitude (V)
 */
inline double inducedEMF(double magneticFluxRate) {
    return std::abs(magneticFluxRate);
}

/**
 * @brief Calculate induced electric field from changing magnetic flux
 *
 * ∇×E = -∂B/∂t
 *
 * For cylindrical symmetry: E = r(dB/dt)/2
 *
 * @param radius Distance from axis (m)
 * @param magneticFieldRate Rate of change of B (T/s)
 * @return Induced electric field (V/m)
 */
inline double inducedElectricField(double radius, double magneticFieldRate) {
    return radius * std::abs(magneticFieldRate) / 2.0;
}

// ============================================================================
// ELECTROMAGNETIC ENERGY AND POYNTING VECTOR
// ============================================================================

/**
 * @brief Calculate energy density of electric field
 *
 * u_E = (ε₀/2)E²
 *
 * @param electricField Electric field magnitude (V/m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Energy density (J/m³)
 */
inline double electricEnergyDensity(double electricField,
                                    double epsilon0 = constants::EPSILON_0) {
    return 0.5 * epsilon0 * electricField * electricField;
}

/**
 * @brief Calculate energy density of magnetic field
 *
 * u_B = B²/(2μ₀)
 *
 * @param magneticField Magnetic field magnitude (T)
 * @param mu0 Permeability of free space (H/m)
 * @return Energy density (J/m³)
 */
inline double magneticEnergyDensity(double magneticField,
                                    double mu0 = constants::MU_0) {
    return magneticField * magneticField / (2.0 * mu0);
}

/**
 * @brief Calculate total electromagnetic energy density
 *
 * u = (ε₀/2)E² + B²/(2μ₀)
 *
 * @param electricField Electric field magnitude (V/m)
 * @param magneticField Magnetic field magnitude (T)
 * @param epsilon0 Permittivity of free space (F/m)
 * @param mu0 Permeability of free space (H/m)
 * @return Total energy density (J/m³)
 */
inline double totalEnergyDensity(double electricField, double magneticField,
                                 double epsilon0 = constants::EPSILON_0,
                                 double mu0 = constants::MU_0) {
    return electricEnergyDensity(electricField, epsilon0) +
           magneticEnergyDensity(magneticField, mu0);
}

/**
 * @brief Calculate Poynting vector magnitude (energy flux)
 *
 * S = (E × B)/μ₀
 *
 * For perpendicular E and B: |S| = EB/μ₀
 *
 * @param electricField Electric field magnitude (V/m)
 * @param magneticField Magnetic field magnitude (T)
 * @param mu0 Permeability of free space (H/m)
 * @return Poynting vector magnitude (W/m²)
 */
inline double poyntingVector(double electricField, double magneticField,
                            double mu0 = constants::MU_0) {
    return (electricField * magneticField) / mu0;
}

/**
 * @brief Calculate electromagnetic momentum density
 *
 * g = S/c² = ε₀(E × B)
 *
 * @param electricField Electric field magnitude (V/m)
 * @param magneticField Magnetic field magnitude (T)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Momentum density magnitude (kg/(m²⋅s))
 */
inline double momentumDensity(double electricField, double magneticField,
                              double epsilon0 = constants::EPSILON_0) {
    return epsilon0 * electricField * magneticField;
}

/**
 * @brief Calculate radiation pressure from EM wave
 *
 * P = u = S/c (for perfect absorption)
 * P = 2u = 2S/c (for perfect reflection)
 *
 * @param intensity Wave intensity S (W/m²)
 * @param speedOfLight Speed of light c (m/s)
 * @param isReflection true for perfect reflection, false for absorption
 * @return Radiation pressure (Pa)
 */
inline double radiationPressure(double intensity, double speedOfLight = constants::SPEED_OF_LIGHT,
                                bool isReflection = false) {
    double pressure = intensity / speedOfLight;
    return isReflection ? 2.0 * pressure : pressure;
}

// ============================================================================
// ELECTROMAGNETIC POTENTIALS
// ============================================================================

/**
 * @brief Calculate electric potential from point charge
 *
 * φ = kQ/r = Q/(4πε₀r)
 *
 * @param charge Charge Q (C)
 * @param distance Distance r (m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Electric potential (V)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double scalarPotential(double charge, double distance,
                              double epsilon0 = constants::EPSILON_0) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return charge / (4.0 * M_PI * epsilon0 * distance);
}

/**
 * @brief Calculate magnetic vector potential for magnetic dipole
 *
 * A = (μ₀/4π)(m × r̂)/r²
 *
 * Magnitude for perpendicular geometry
 *
 * @param magneticMoment Magnetic dipole moment m (A⋅m²)
 * @param distance Distance r (m)
 * @param mu0 Permeability of free space (H/m)
 * @return Vector potential magnitude (T⋅m)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double vectorPotentialDipole(double magneticMoment, double distance,
                                    double mu0 = constants::MU_0) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return (mu0 * magneticMoment) / (4.0 * M_PI * distance * distance);
}

/**
 * @brief Verify Coulomb gauge condition
 *
 * ∇·A = 0
 *
 * Check if divergence of vector potential is zero
 *
 * @param divergenceA Divergence of vector potential
 * @param tolerance Tolerance for zero check
 * @return true if Coulomb gauge is satisfied
 */
inline bool isCoulombGauge(double divergenceA, double tolerance = 1e-10) {
    return std::abs(divergenceA) < tolerance;
}

/**
 * @brief Verify Lorenz gauge condition
 *
 * ∇·A + (1/c²)(∂φ/∂t) = 0
 *
 * @param divergenceA Divergence of A
 * @param timeDerivativePhii Time derivative of φ
 * @param speedOfLight Speed of light c
 * @param tolerance Tolerance for zero check
 * @return true if Lorenz gauge is satisfied
 */
inline bool isLorenzGauge(double divergenceA, double timeDerivativePhi,
                          double speedOfLight = constants::SPEED_OF_LIGHT,
                          double tolerance = 1e-10) {
    double gauge = divergenceA + timeDerivativePhi / (speedOfLight * speedOfLight);
    return std::abs(gauge) < tolerance;
}

// ============================================================================
// ELECTROMAGNETIC WAVES IN VACUUM
// ============================================================================

/**
 * @brief Calculate E/B ratio in electromagnetic wave
 *
 * E/B = c
 *
 * @param speedOfLight Speed of light (m/s)
 * @return E/B ratio (m/s)
 */
inline double eBRatio(double speedOfLight = constants::SPEED_OF_LIGHT) {
    return speedOfLight;
}

/**
 * @brief Calculate characteristic impedance of free space
 *
 * Z₀ = √(μ₀/ε₀) = μ₀c ≈ 377 Ω
 *
 * @param mu0 Permeability of free space (H/m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Impedance of free space (Ω)
 */
inline double impedanceOfFreeSpace(double mu0 = constants::MU_0,
                                   double epsilon0 = constants::EPSILON_0) {
    return std::sqrt(mu0 / epsilon0);
}

/**
 * @brief Calculate wave intensity from field amplitude
 *
 * I = (ε₀c/2)E₀² = cB₀²/(2μ₀)
 *
 * @param electricFieldAmplitude E₀ (V/m)
 * @param speedOfLight c (m/s)
 * @param epsilon0 Permittivity (F/m)
 * @return Intensity (W/m²)
 */
inline double waveIntensity(double electricFieldAmplitude,
                           double speedOfLight = constants::SPEED_OF_LIGHT,
                           double epsilon0 = constants::EPSILON_0) {
    return 0.5 * speedOfLight * epsilon0 * electricFieldAmplitude * electricFieldAmplitude;
}

// ============================================================================
// RADIATION FROM ACCELERATING CHARGES
// ============================================================================

/**
 * @brief Calculate Larmor formula for radiated power
 *
 * P = (μ₀q²a²)/(6πc)
 *
 * Power radiated by non-relativistic accelerating charge
 *
 * @param charge Charge q (C)
 * @param acceleration Acceleration a (m/s²)
 * @param mu0 Permeability of free space (H/m)
 * @param speedOfLight Speed of light (m/s)
 * @return Radiated power (W)
 */
inline double larmorPower(double charge, double acceleration,
                         double mu0 = constants::MU_0,
                         double speedOfLight = constants::SPEED_OF_LIGHT) {
    double q2 = charge * charge;
    double a2 = acceleration * acceleration;
    return (mu0 * q2 * a2) / (6.0 * M_PI * speedOfLight);
}

/**
 * @brief Calculate classical electron radius
 *
 * r_e = e²/(4πε₀m_e c²)
 *
 * @param electronCharge Elementary charge e (C)
 * @param electronMass Electron mass m_e (kg)
 * @param speedOfLight Speed of light (m/s)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Classical electron radius (m)
 */
inline double classicalElectronRadius(double electronCharge, double electronMass,
                                      double speedOfLight = constants::SPEED_OF_LIGHT,
                                      double epsilon0 = constants::EPSILON_0) {
    double e2 = electronCharge * electronCharge;
    double c2 = speedOfLight * speedOfLight;
    return e2 / (4.0 * M_PI * epsilon0 * electronMass * c2);
}

/**
 * @brief Calculate radiation reaction force (Abraham-Lorentz)
 *
 * F_rad = (μ₀q²)/(6πc) × (da/dt)
 *
 * Self-force on accelerating charge due to radiation
 *
 * @param charge Charge q (C)
 * @param jerk Rate of change of acceleration (m/s³)
 * @param mu0 Permeability of free space (H/m)
 * @param speedOfLight Speed of light (m/s)
 * @return Radiation reaction force (N)
 */
inline double radiationReactionForce(double charge, double jerk,
                                     double mu0 = constants::MU_0,
                                     double speedOfLight = constants::SPEED_OF_LIGHT) {
    double q2 = charge * charge;
    return (mu0 * q2 * jerk) / (6.0 * M_PI * speedOfLight);
}

// ============================================================================
// MULTIPOLE MOMENTS
// ============================================================================

/**
 * @brief Calculate electric dipole moment
 *
 * p = qd
 *
 * @param charge Charge magnitude (C)
 * @param separation Separation distance (m)
 * @return Dipole moment (C⋅m)
 */
inline double electricDipoleMoment(double charge, double separation) {
    return charge * separation;
}

/**
 * @brief Calculate electric field from dipole (on axis)
 *
 * E = (1/(4πε₀)) × (2p/r³)
 *
 * Far-field approximation along dipole axis
 *
 * @param dipoleMoment Dipole moment p (C⋅m)
 * @param distance Distance r (m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Electric field (V/m)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double dipoleFieldOnAxis(double dipoleMoment, double distance,
                                double epsilon0 = constants::EPSILON_0) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return (2.0 * dipoleMoment) / (4.0 * M_PI * epsilon0 * distance * distance * distance);
}

/**
 * @brief Calculate electric quadrupole moment (simplified)
 *
 * Q = q(3z² - r²)
 *
 * For point charge at position
 *
 * @param charge Charge (C)
 * @param z Position along axis (m)
 * @param r Radial distance (m)
 * @return Quadrupole moment component (C⋅m²)
 */
inline double electricQuadrupoleMoment(double charge, double z, double r) {
    return charge * (3.0 * z * z - r * r);
}

// ============================================================================
// BOUNDARY CONDITIONS
// ============================================================================

/**
 * @brief Calculate discontinuity in normal E across surface charge
 *
 * E_⊥(above) - E_⊥(below) = σ/ε₀
 *
 * @param surfaceChargeDensity σ (C/m²)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Discontinuity in normal E component (V/m)
 */
inline double normalEDiscontinuity(double surfaceChargeDensity,
                                   double epsilon0 = constants::EPSILON_0) {
    return surfaceChargeDensity / epsilon0;
}

/**
 * @brief Calculate discontinuity in tangential B across surface current
 *
 * B_∥(above) - B_∥(below) = μ₀K
 *
 * where K is surface current density
 *
 * @param surfaceCurrentDensity K (A/m)
 * @param mu0 Permeability of free space (H/m)
 * @return Discontinuity in tangential B component (T)
 */
inline double tangentialBDiscontinuity(double surfaceCurrentDensity,
                                       double mu0 = constants::MU_0) {
    return mu0 * surfaceCurrentDensity;
}

} // namespace maxwell
} // namespace physics

#endif // PHYSICS_MAXWELL_EQUATIONS_HPP
