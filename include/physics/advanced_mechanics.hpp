#ifndef PHYSICS_ADVANCED_MECHANICS_HPP
#define PHYSICS_ADVANCED_MECHANICS_HPP

#include <cmath>
#include <stdexcept>
#include <array>
#include <vector>

/**
 * @file advanced_mechanics.hpp
 * @brief Advanced classical mechanics including orbital dynamics and variational methods
 *
 * This module implements:
 * - Polar coordinates (conversions, velocities, accelerations)
 * - Relative motion between reference frames
 * - Point-dynamics in fixed coordinate systems
 * - Conservative force fields and potential energy
 * - Orbital equations and dynamics
 * - Kepler's laws and orbital elements
 * - Virial theorem
 * - Basic tensor operations
 * - Variational calculus (Euler-Lagrange equations)
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace advanced_mechanics {

/**
 * @namespace constants
 * @brief Astronomical and physical constants
 */
namespace constants {
    constexpr double GRAVITATIONAL_CONSTANT = 6.67430e-11;  // G (m³ kg⁻¹ s⁻²)
    constexpr double SOLAR_MASS = 1.989e30;                 // kg
    constexpr double EARTH_MASS = 5.972e24;                 // kg
    constexpr double AU = 1.496e11;                         // Astronomical Unit (m)
}

// ============================================================================
// POLAR COORDINATES
// ============================================================================

/**
 * @brief Convert Cartesian to polar coordinates
 *
 * r = √(x² + y²)
 * θ = arctan(y/x)
 *
 * @param x X coordinate (m)
 * @param y Y coordinate (m)
 * @return [r, θ] where r is radial distance (m) and θ is angle (radians)
 */
inline std::array<double, 2> cartesianToPolar(double x, double y) {
    double r = std::sqrt(x * x + y * y);
    double theta = std::atan2(y, x);
    return {r, theta};
}

/**
 * @brief Convert polar to Cartesian coordinates
 *
 * x = r cos(θ)
 * y = r sin(θ)
 *
 * @param r Radial distance (m)
 * @param theta Angle (radians)
 * @return [x, y] Cartesian coordinates (m)
 */
inline std::array<double, 2> polarToCartesian(double r, double theta) {
    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    return {x, y};
}

/**
 * @brief Calculate radial velocity in polar coordinates
 *
 * v_r = dr/dt
 *
 * @param radialRate Rate of change of r (m/s)
 * @return Radial velocity (m/s)
 */
inline double radialVelocity(double radialRate) {
    return radialRate;
}

/**
 * @brief Calculate tangential velocity in polar coordinates
 *
 * v_θ = r × dθ/dt
 *
 * @param r Radial distance (m)
 * @param angularVelocity Angular velocity (rad/s)
 * @return Tangential velocity (m/s)
 */
inline double tangentialVelocity(double r, double angularVelocity) {
    return r * angularVelocity;
}

/**
 * @brief Calculate radial acceleration in polar coordinates
 *
 * a_r = d²r/dt² - r(dθ/dt)²
 *
 * Includes centripetal acceleration term
 *
 * @param radialAccel Second derivative of r (m/s²)
 * @param r Radial distance (m)
 * @param angularVelocity Angular velocity (rad/s)
 * @return Radial acceleration (m/s²)
 */
inline double radialAcceleration(double radialAccel, double r, double angularVelocity) {
    return radialAccel - r * angularVelocity * angularVelocity;
}

/**
 * @brief Calculate tangential acceleration in polar coordinates
 *
 * a_θ = r × d²θ/dt² + 2(dr/dt)(dθ/dt)
 *
 * Includes Coriolis-like term
 *
 * @param r Radial distance (m)
 * @param angularAccel Angular acceleration (rad/s²)
 * @param radialVel Radial velocity (m/s)
 * @param angularVel Angular velocity (rad/s)
 * @return Tangential acceleration (m/s²)
 */
inline double tangentialAcceleration(double r, double angularAccel,
                                     double radialVel, double angularVel) {
    return r * angularAccel + 2.0 * radialVel * angularVel;
}

// ============================================================================
// RELATIVE MOTION
// ============================================================================

/**
 * @brief Calculate relative velocity between two frames
 *
 * v_rel = v_A - v_B
 *
 * @param velocityA Velocity in frame A (m/s)
 * @param velocityB Velocity in frame B (m/s)
 * @return Relative velocity (m/s)
 */
inline double relativeVelocity(double velocityA, double velocityB) {
    return velocityA - velocityB;
}

/**
 * @brief Calculate velocity in moving frame (Galilean transformation)
 *
 * v' = v - V
 *
 * where V is velocity of moving frame
 *
 * @param velocityLab Velocity in lab frame (m/s)
 * @param frameVelocity Velocity of moving frame (m/s)
 * @return Velocity in moving frame (m/s)
 */
inline double velocityInMovingFrame(double velocityLab, double frameVelocity) {
    return velocityLab - frameVelocity;
}

/**
 * @brief Calculate relative acceleration (in inertial frames)
 *
 * a_rel = a_A - a_B
 *
 * @param accelA Acceleration in frame A (m/s²)
 * @param accelB Acceleration in frame B (m/s²)
 * @return Relative acceleration (m/s²)
 */
inline double relativeAcceleration(double accelA, double accelB) {
    return accelA - accelB;
}

// ============================================================================
// CONSERVATIVE FORCE FIELDS
// ============================================================================

/**
 * @brief Calculate force from potential energy (1D)
 *
 * F = -dU/dx
 *
 * @param potentialGradient Gradient of potential (J/m)
 * @return Force (N)
 */
inline double forceFromPotential(double potentialGradient) {
    return -potentialGradient;
}

/**
 * @brief Calculate gravitational potential energy
 *
 * U = -GMm/r
 *
 * @param mass1 First mass (kg)
 * @param mass2 Second mass (kg)
 * @param distance Separation (m)
 * @param G Gravitational constant (m³ kg⁻¹ s⁻²)
 * @return Potential energy (J)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double gravitationalPotential(double mass1, double mass2, double distance,
                                     double G = constants::GRAVITATIONAL_CONSTANT) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return -G * mass1 * mass2 / distance;
}

/**
 * @brief Calculate escape velocity from gravitational field
 *
 * v_esc = √(2GM/r)
 *
 * @param mass Mass of central body (kg)
 * @param distance Distance from center (m)
 * @param G Gravitational constant (m³ kg⁻¹ s⁻²)
 * @return Escape velocity (m/s)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double escapeVelocity(double mass, double distance,
                            double G = constants::GRAVITATIONAL_CONSTANT) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return std::sqrt(2.0 * G * mass / distance);
}

/**
 * @brief Check if force field is conservative (from curl = 0)
 *
 * For 2D field: ∂F_y/∂x = ∂F_x/∂y
 *
 * @param dFy_dx Partial derivative of F_y with respect to x
 * @param dFx_dy Partial derivative of F_x with respect to y
 * @param tolerance Tolerance for comparison
 * @return true if field is conservative
 */
inline bool isConservative(double dFy_dx, double dFx_dy, double tolerance = 1e-10) {
    return std::abs(dFy_dx - dFx_dy) < tolerance;
}

/**
 * @brief Calculate mechanical energy (total energy)
 *
 * E = K + U
 *
 * @param kineticEnergy Kinetic energy (J)
 * @param potentialEnergy Potential energy (J)
 * @return Total mechanical energy (J)
 */
inline double mechanicalEnergy(double kineticEnergy, double potentialEnergy) {
    return kineticEnergy + potentialEnergy;
}

// ============================================================================
// ORBITAL EQUATIONS
// ============================================================================

/**
 * @brief Calculate orbital angular momentum
 *
 * L = m × r × v_⊥
 *
 * @param mass Mass of orbiting body (kg)
 * @param distance Distance from center (m)
 * @param tangentialVelocity Tangential velocity (m/s)
 * @return Angular momentum (kg⋅m²/s)
 */
inline double orbitalAngularMomentum(double mass, double distance, double tangentialVelocity) {
    return mass * distance * tangentialVelocity;
}

/**
 * @brief Calculate specific angular momentum
 *
 * h = r × v = L/m
 *
 * @param distance Distance (m)
 * @param tangentialVelocity Tangential velocity (m/s)
 * @return Specific angular momentum (m²/s)
 */
inline double specificAngularMomentum(double distance, double tangentialVelocity) {
    return distance * tangentialVelocity;
}

/**
 * @brief Calculate orbital velocity for circular orbit
 *
 * v = √(GM/r)
 *
 * @param centralMass Mass of central body (kg)
 * @param radius Orbital radius (m)
 * @param G Gravitational constant (m³ kg⁻¹ s⁻²)
 * @return Orbital velocity (m/s)
 * @throws std::invalid_argument if radius is non-positive
 */
inline double circularOrbitVelocity(double centralMass, double radius,
                                    double G = constants::GRAVITATIONAL_CONSTANT) {
    if (radius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return std::sqrt(G * centralMass / radius);
}

/**
 * @brief Calculate orbital period for circular orbit
 *
 * T = 2π√(r³/GM)
 *
 * @param radius Orbital radius (m)
 * @param centralMass Mass of central body (kg)
 * @param G Gravitational constant (m³ kg⁻¹ s⁻²)
 * @return Orbital period (s)
 * @throws std::invalid_argument if radius is non-positive
 */
inline double circularOrbitPeriod(double radius, double centralMass,
                                  double G = constants::GRAVITATIONAL_CONSTANT) {
    if (radius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return 2.0 * M_PI * std::sqrt((radius * radius * radius) / (G * centralMass));
}

/**
 * @brief Calculate specific orbital energy
 *
 * ε = -GM/(2a)
 *
 * where a is semi-major axis (negative for bound orbits)
 *
 * @param centralMass Mass of central body (kg)
 * @param semiMajorAxis Semi-major axis (m)
 * @param G Gravitational constant (m³ kg⁻¹ s⁻²)
 * @return Specific orbital energy (J/kg)
 * @throws std::invalid_argument if semi-major axis is non-positive
 */
inline double specificOrbitalEnergy(double centralMass, double semiMajorAxis,
                                    double G = constants::GRAVITATIONAL_CONSTANT) {
    if (semiMajorAxis <= 0.0) {
        throw std::invalid_argument("Semi-major axis must be positive");
    }
    return -G * centralMass / (2.0 * semiMajorAxis);
}

// ============================================================================
// KEPLER'S ORBITAL EQUATIONS
// ============================================================================

/**
 * @brief Kepler's Third Law: T² ∝ a³
 *
 * T² = (4π²/GM) × a³
 *
 * @param semiMajorAxis Semi-major axis (m)
 * @param centralMass Mass of central body (kg)
 * @param G Gravitational constant (m³ kg⁻¹ s⁻²)
 * @return Orbital period (s)
 * @throws std::invalid_argument if semi-major axis is non-positive
 */
inline double keplersThirdLaw(double semiMajorAxis, double centralMass,
                              double G = constants::GRAVITATIONAL_CONSTANT) {
    if (semiMajorAxis <= 0.0) {
        throw std::invalid_argument("Semi-major axis must be positive");
    }
    double a3 = semiMajorAxis * semiMajorAxis * semiMajorAxis;
    return 2.0 * M_PI * std::sqrt(a3 / (G * centralMass));
}

/**
 * @brief Calculate semi-major axis from period (inverse of Kepler's Third Law)
 *
 * a = ∛[(GMT²)/(4π²)]
 *
 * @param period Orbital period (s)
 * @param centralMass Mass of central body (kg)
 * @param G Gravitational constant (m³ kg⁻¹ s⁻²)
 * @return Semi-major axis (m)
 * @throws std::invalid_argument if period is non-positive
 */
inline double semiMajorAxisFromPeriod(double period, double centralMass,
                                      double G = constants::GRAVITATIONAL_CONSTANT) {
    if (period <= 0.0) {
        throw std::invalid_argument("Period must be positive");
    }
    double T2 = period * period;
    return std::cbrt((G * centralMass * T2) / (4.0 * M_PI * M_PI));
}

/**
 * @brief Calculate orbital eccentricity from perihelion and aphelion
 *
 * e = (r_a - r_p)/(r_a + r_p)
 *
 * @param aphelion Maximum distance (m)
 * @param perihelion Minimum distance (m)
 * @return Eccentricity (0 ≤ e < 1 for ellipse)
 * @throws std::invalid_argument if distances are non-positive
 */
inline double eccentricity(double aphelion, double perihelion) {
    if (aphelion <= 0.0 || perihelion <= 0.0) {
        throw std::invalid_argument("Distances must be positive");
    }
    if (perihelion > aphelion) {
        throw std::invalid_argument("Perihelion must be less than aphelion");
    }
    return (aphelion - perihelion) / (aphelion + perihelion);
}

/**
 * @brief Calculate semi-major axis from perihelion and aphelion
 *
 * a = (r_a + r_p)/2
 *
 * @param aphelion Maximum distance (m)
 * @param perihelion Minimum distance (m)
 * @return Semi-major axis (m)
 * @throws std::invalid_argument if distances are non-positive
 */
inline double semiMajorAxis(double aphelion, double perihelion) {
    if (aphelion <= 0.0 || perihelion <= 0.0) {
        throw std::invalid_argument("Distances must be positive");
    }
    return (aphelion + perihelion) / 2.0;
}

/**
 * @brief Calculate semi-minor axis
 *
 * b = a√(1 - e²)
 *
 * @param semiMajorAxis Semi-major axis (m)
 * @param eccentricity Eccentricity
 * @return Semi-minor axis (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double semiMinorAxis(double semiMajorAxis, double eccentricity) {
    if (semiMajorAxis <= 0.0) {
        throw std::invalid_argument("Semi-major axis must be positive");
    }
    if (eccentricity < 0.0 || eccentricity >= 1.0) {
        throw std::invalid_argument("Eccentricity must be in range [0, 1) for ellipse");
    }
    return semiMajorAxis * std::sqrt(1.0 - eccentricity * eccentricity);
}

/**
 * @brief Calculate perihelion distance
 *
 * r_p = a(1 - e)
 *
 * @param semiMajorAxis Semi-major axis (m)
 * @param eccentricity Eccentricity
 * @return Perihelion distance (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double perihelionDistance(double semiMajorAxis, double eccentricity) {
    if (semiMajorAxis <= 0.0) {
        throw std::invalid_argument("Semi-major axis must be positive");
    }
    if (eccentricity < 0.0 || eccentricity >= 1.0) {
        throw std::invalid_argument("Eccentricity must be in range [0, 1)");
    }
    return semiMajorAxis * (1.0 - eccentricity);
}

/**
 * @brief Calculate aphelion distance
 *
 * r_a = a(1 + e)
 *
 * @param semiMajorAxis Semi-major axis (m)
 * @param eccentricity Eccentricity
 * @return Aphelion distance (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double aphelionDistance(double semiMajorAxis, double eccentricity) {
    if (semiMajorAxis <= 0.0) {
        throw std::invalid_argument("Semi-major axis must be positive");
    }
    if (eccentricity < 0.0 || eccentricity >= 1.0) {
        throw std::invalid_argument("Eccentricity must be in range [0, 1)");
    }
    return semiMajorAxis * (1.0 + eccentricity);
}

/**
 * @brief Calculate orbital radius at given true anomaly
 *
 * r = a(1 - e²)/(1 + e cos(ν))
 *
 * @param semiMajorAxis Semi-major axis (m)
 * @param eccentricity Eccentricity
 * @param trueAnomaly True anomaly (radians)
 * @return Orbital radius (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double orbitalRadius(double semiMajorAxis, double eccentricity, double trueAnomaly) {
    if (semiMajorAxis <= 0.0) {
        throw std::invalid_argument("Semi-major axis must be positive");
    }
    if (eccentricity < 0.0 || eccentricity >= 1.0) {
        throw std::invalid_argument("Eccentricity must be in range [0, 1)");
    }

    double numerator = semiMajorAxis * (1.0 - eccentricity * eccentricity);
    double denominator = 1.0 + eccentricity * std::cos(trueAnomaly);

    if (denominator <= 0.0) {
        throw std::invalid_argument("Invalid true anomaly for this eccentricity");
    }

    return numerator / denominator;
}

/**
 * @brief Calculate velocity at given point in elliptical orbit (vis-viva equation)
 *
 * v² = GM(2/r - 1/a)
 *
 * @param centralMass Mass of central body (kg)
 * @param distance Current distance from center (m)
 * @param semiMajorAxis Semi-major axis (m)
 * @param G Gravitational constant (m³ kg⁻¹ s⁻²)
 * @return Orbital velocity (m/s)
 * @throws std::invalid_argument if distances are non-positive
 */
inline double visVivaEquation(double centralMass, double distance, double semiMajorAxis,
                              double G = constants::GRAVITATIONAL_CONSTANT) {
    if (distance <= 0.0 || semiMajorAxis <= 0.0) {
        throw std::invalid_argument("Distances must be positive");
    }
    double v2 = G * centralMass * (2.0/distance - 1.0/semiMajorAxis);
    if (v2 < 0.0) {
        throw std::invalid_argument("Invalid orbital parameters");
    }
    return std::sqrt(v2);
}

// ============================================================================
// VIRIAL THEOREM
// ============================================================================

/**
 * @brief Apply virial theorem for gravitational systems
 *
 * 2⟨K⟩ + ⟨U⟩ = 0
 * ⟨K⟩ = -⟨U⟩/2 = E (for bound systems)
 *
 * For inverse-square force: time-averaged kinetic energy equals
 * half the magnitude of time-averaged potential energy
 *
 * @param potentialEnergy Time-averaged potential energy (J)
 * @return Time-averaged kinetic energy (J)
 */
inline double virialTheoremKinetic(double potentialEnergy) {
    return -potentialEnergy / 2.0;
}

/**
 * @brief Calculate total energy from virial theorem
 *
 * E = ⟨K⟩ + ⟨U⟩ = ⟨U⟩/2
 *
 * @param averagePotentialEnergy Time-averaged potential energy (J)
 * @return Total energy (J)
 */
inline double virialTheoremEnergy(double averagePotentialEnergy) {
    return averagePotentialEnergy / 2.0;
}

/**
 * @brief Verify virial theorem condition
 *
 * Check if 2K + U ≈ 0
 *
 * @param kineticEnergy Kinetic energy (J)
 * @param potentialEnergy Potential energy (J)
 * @param tolerance Tolerance for verification
 * @return true if virial theorem is satisfied
 */
inline bool verifyVirialTheorem(double kineticEnergy, double potentialEnergy,
                               double tolerance = 1e-6) {
    double virial = 2.0 * kineticEnergy + potentialEnergy;
    double scale = std::abs(potentialEnergy);
    if (scale < 1e-15) scale = 1.0;
    return std::abs(virial / scale) < tolerance;
}

// ============================================================================
// TENSOR OPERATIONS (BASIC)
// ============================================================================

/**
 * @brief Calculate moment of inertia tensor element (2D)
 *
 * I_ij = m(r²δ_ij - r_i r_j)
 *
 * For point mass, simplified calculation
 *
 * @param mass Mass (kg)
 * @param x X coordinate (m)
 * @param y Y coordinate (m)
 * @param i First index (0 or 1)
 * @param j Second index (0 or 1)
 * @return Tensor element (kg⋅m²)
 */
inline double inertiatensorElement2D(double mass, double x, double y, int i, int j) {
    if (i < 0 || i > 1 || j < 0 || j > 1) {
        throw std::invalid_argument("Indices must be 0 or 1 for 2D");
    }

    double coords[2] = {x, y};
    double r2 = x*x + y*y;

    if (i == j) {
        return mass * (r2 - coords[i] * coords[i]);
    } else {
        return -mass * coords[i] * coords[j];
    }
}

/**
 * @brief Calculate trace of 2D inertia tensor
 *
 * Tr(I) = I_00 + I_11
 *
 * @param I00 Element (0,0)
 * @param I11 Element (1,1)
 * @return Trace
 */
inline double tensorTrace2D(double I00, double I11) {
    return I00 + I11;
}

// ============================================================================
// VARIATIONAL CALCULUS
// ============================================================================

/**
 * @brief Calculate action integral for free particle
 *
 * S = ∫ L dt = ∫ (1/2)mv² dt
 *
 * For constant velocity: S = (1/2)mv²t
 *
 * @param mass Mass (kg)
 * @param velocity Velocity (m/s)
 * @param time Time interval (s)
 * @return Action (J⋅s)
 * @throws std::invalid_argument if time is negative
 */
inline double actionFreeParticle(double mass, double velocity, double time) {
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    return 0.5 * mass * velocity * velocity * time;
}

/**
 * @brief Calculate Lagrangian for free particle
 *
 * L = T - V = (1/2)mv²
 *
 * (V = 0 for free particle)
 *
 * @param mass Mass (kg)
 * @param velocity Velocity (m/s)
 * @return Lagrangian (J)
 */
inline double lagrangianFreeParticle(double mass, double velocity) {
    return 0.5 * mass * velocity * velocity;
}

/**
 * @brief Calculate Lagrangian for particle in potential
 *
 * L = T - V = (1/2)mv² - U
 *
 * @param kineticEnergy Kinetic energy (J)
 * @param potentialEnergy Potential energy (J)
 * @return Lagrangian (J)
 */
inline double lagrangian(double kineticEnergy, double potentialEnergy) {
    return kineticEnergy - potentialEnergy;
}

/**
 * @brief Calculate Hamiltonian
 *
 * H = T + V = p²/(2m) + U
 *
 * @param kineticEnergy Kinetic energy (J)
 * @param potentialEnergy Potential energy (J)
 * @return Hamiltonian (J)
 */
inline double hamiltonian(double kineticEnergy, double potentialEnergy) {
    return kineticEnergy + potentialEnergy;
}

/**
 * @brief Calculate canonical momentum
 *
 * p = ∂L/∂v = mv (for standard kinetic energy)
 *
 * @param mass Mass (kg)
 * @param velocity Velocity (m/s)
 * @return Momentum (kg⋅m/s)
 */
inline double canonicalMomentum(double mass, double velocity) {
    return mass * velocity;
}

/**
 * @brief Verify conservation of Hamiltonian (energy)
 *
 * dH/dt = 0 for time-independent systems
 *
 * @param hamiltonian1 Hamiltonian at time t1 (J)
 * @param hamiltonian2 Hamiltonian at time t2 (J)
 * @param tolerance Tolerance for comparison
 * @return true if energy is conserved
 */
inline bool isHamiltonianConserved(double hamiltonian1, double hamiltonian2,
                                  double tolerance = 1e-6) {
    double scale = std::abs(hamiltonian1);
    if (scale < 1e-15) scale = 1.0;
    return std::abs((hamiltonian2 - hamiltonian1) / scale) < tolerance;
}

/**
 * @brief Calculate reduced mass for two-body problem
 *
 * μ = (m₁ × m₂)/(m₁ + m₂)
 *
 * @param mass1 First mass (kg)
 * @param mass2 Second mass (kg)
 * @return Reduced mass (kg)
 * @throws std::invalid_argument if masses are non-positive
 */
inline double reducedMass(double mass1, double mass2) {
    if (mass1 <= 0.0 || mass2 <= 0.0) {
        throw std::invalid_argument("Masses must be positive");
    }
    return (mass1 * mass2) / (mass1 + mass2);
}

/**
 * @brief Calculate effective one-body problem mass
 *
 * M_eff = m₁ + m₂
 *
 * @param mass1 First mass (kg)
 * @param mass2 Second mass (kg)
 * @return Effective mass (kg)
 */
inline double effectiveMass(double mass1, double mass2) {
    return mass1 + mass2;
}

} // namespace advanced_mechanics
} // namespace physics

#endif // PHYSICS_ADVANCED_MECHANICS_HPP
