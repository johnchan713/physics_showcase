#ifndef PHYSICS_QUANTUM_BASICS_HPP
#define PHYSICS_QUANTUM_BASICS_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file quantum_basics.hpp
 * @brief Basic quantum mechanics formulas and relations
 *
 * This module implements:
 * - De Broglie wavelength and matter waves
 * - Compton scattering
 * - Heisenberg uncertainty principle
 * - Bohr model of hydrogen
 * - Energy levels in simple quantum systems
 * - Wave-particle duality relations
 * - Photoelectric effect
 * - Quantum tunneling probability
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace quantum_basics {

/**
 * @namespace constants
 * @brief Physical constants for quantum mechanics
 */
namespace constants {
    constexpr double PLANCK_H = 6.62607015e-34;        // Planck's constant (J⋅s)
    constexpr double HBAR = 1.054571817e-34;           // Reduced Planck's constant ℏ = h/(2π) (J⋅s)
    constexpr double ELECTRON_MASS = 9.1093837015e-31; // Electron mass (kg)
    constexpr double PROTON_MASS = 1.67262192369e-27;  // Proton mass (kg)
    constexpr double ELEMENTARY_CHARGE = 1.602176634e-19; // Elementary charge (C)
    constexpr double SPEED_OF_LIGHT = 299792458.0;     // Speed of light (m/s)
    constexpr double BOHR_RADIUS = 5.29177210903e-11;  // Bohr radius (m)
    constexpr double RYDBERG_ENERGY = 13.605693122994; // Rydberg energy (eV)
    constexpr double FINE_STRUCTURE = 1.0/137.035999084; // Fine structure constant α
}

// ============================================================================
// DE BROGLIE WAVELENGTH
// ============================================================================

/**
 * @brief Calculate de Broglie wavelength from momentum
 *
 * λ = h/p
 *
 * @param momentum Momentum p (kg⋅m/s)
 * @param planckH Planck's constant (J⋅s)
 * @return De Broglie wavelength (m)
 * @throws std::invalid_argument if momentum is zero
 */
inline double deBroglieWavelength(double momentum, double planckH = constants::PLANCK_H) {
    if (std::abs(momentum) < 1e-50) {
        throw std::invalid_argument("Momentum must be non-zero");
    }
    return planckH / std::abs(momentum);
}

/**
 * @brief Calculate de Broglie wavelength from velocity (non-relativistic)
 *
 * λ = h/(mv)
 *
 * @param mass Mass (kg)
 * @param velocity Velocity (m/s)
 * @param planckH Planck's constant (J⋅s)
 * @return De Broglie wavelength (m)
 * @throws std::invalid_argument if mass or velocity is zero
 */
inline double deBroglieWavelengthFromVelocity(double mass, double velocity,
                                               double planckH = constants::PLANCK_H) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (std::abs(velocity) < 1e-15) {
        throw std::invalid_argument("Velocity must be non-zero");
    }
    return planckH / (mass * velocity);
}

/**
 * @brief Calculate de Broglie wavelength from kinetic energy (non-relativistic)
 *
 * λ = h/√(2mE)
 *
 * @param mass Mass (kg)
 * @param kineticEnergy Kinetic energy (J)
 * @param planckH Planck's constant (J⋅s)
 * @return De Broglie wavelength (m)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double deBroglieWavelengthFromEnergy(double mass, double kineticEnergy,
                                             double planckH = constants::PLANCK_H) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (kineticEnergy <= 0.0) {
        throw std::invalid_argument("Kinetic energy must be positive");
    }
    return planckH / std::sqrt(2.0 * mass * kineticEnergy);
}

/**
 * @brief Calculate momentum from de Broglie wavelength
 *
 * p = h/λ
 *
 * @param wavelength De Broglie wavelength (m)
 * @param planckH Planck's constant (J⋅s)
 * @return Momentum (kg⋅m/s)
 * @throws std::invalid_argument if wavelength is non-positive
 */
inline double momentumFromWavelength(double wavelength, double planckH = constants::PLANCK_H) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return planckH / wavelength;
}

// ============================================================================
// COMPTON SCATTERING
// ============================================================================

/**
 * @brief Calculate Compton wavelength shift
 *
 * Δλ = λ' - λ = (h/(m_e c))(1 - cos θ)
 *
 * @param scatteringAngle Scattering angle θ (radians)
 * @param electronMass Electron mass (kg)
 * @param planckH Planck's constant (J⋅s)
 * @param speedOfLight Speed of light (m/s)
 * @return Wavelength shift Δλ (m)
 */
inline double comptonShift(double scatteringAngle,
                           double electronMass = constants::ELECTRON_MASS,
                           double planckH = constants::PLANCK_H,
                           double speedOfLight = constants::SPEED_OF_LIGHT) {
    double comptonWavelength = planckH / (electronMass * speedOfLight);
    return comptonWavelength * (1.0 - std::cos(scatteringAngle));
}

/**
 * @brief Calculate Compton wavelength
 *
 * λ_C = h/(m_e c)
 *
 * Characteristic wavelength for Compton scattering
 *
 * @param electronMass Electron mass (kg)
 * @param planckH Planck's constant (J⋅s)
 * @param speedOfLight Speed of light (m/s)
 * @return Compton wavelength (m)
 */
inline double comptonWavelength(double electronMass = constants::ELECTRON_MASS,
                                double planckH = constants::PLANCK_H,
                                double speedOfLight = constants::SPEED_OF_LIGHT) {
    return planckH / (electronMass * speedOfLight);
}

/**
 * @brief Calculate scattered photon wavelength in Compton effect
 *
 * λ' = λ + Δλ
 *
 * @param incidentWavelength Initial photon wavelength λ (m)
 * @param scatteringAngle Scattering angle θ (radians)
 * @param electronMass Electron mass (kg)
 * @param planckH Planck's constant (J⋅s)
 * @param speedOfLight Speed of light (m/s)
 * @return Scattered wavelength λ' (m)
 */
inline double scatteredWavelength(double incidentWavelength, double scatteringAngle,
                                   double electronMass = constants::ELECTRON_MASS,
                                   double planckH = constants::PLANCK_H,
                                   double speedOfLight = constants::SPEED_OF_LIGHT) {
    double shift = comptonShift(scatteringAngle, electronMass, planckH, speedOfLight);
    return incidentWavelength + shift;
}

/**
 * @brief Calculate energy of scattered photon in Compton effect
 *
 * E' = E/(1 + (E/(m_e c²))(1 - cos θ))
 *
 * @param incidentEnergy Initial photon energy (J)
 * @param scatteringAngle Scattering angle θ (radians)
 * @param electronMass Electron mass (kg)
 * @param speedOfLight Speed of light (m/s)
 * @return Scattered photon energy (J)
 */
inline double comptonScatteredEnergy(double incidentEnergy, double scatteringAngle,
                                      double electronMass = constants::ELECTRON_MASS,
                                      double speedOfLight = constants::SPEED_OF_LIGHT) {
    double restEnergy = electronMass * speedOfLight * speedOfLight;
    double denominator = 1.0 + (incidentEnergy / restEnergy) * (1.0 - std::cos(scatteringAngle));
    return incidentEnergy / denominator;
}

// ============================================================================
// HEISENBERG UNCERTAINTY PRINCIPLE
// ============================================================================

/**
 * @brief Calculate minimum position uncertainty from momentum uncertainty
 *
 * Δx ≥ ℏ/(2Δp)
 *
 * @param momentumUncertainty Δp (kg⋅m/s)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Minimum position uncertainty Δx (m)
 * @throws std::invalid_argument if momentum uncertainty is non-positive
 */
inline double positionUncertainty(double momentumUncertainty, double hbar = constants::HBAR) {
    if (momentumUncertainty <= 0.0) {
        throw std::invalid_argument("Momentum uncertainty must be positive");
    }
    return hbar / (2.0 * momentumUncertainty);
}

/**
 * @brief Calculate minimum momentum uncertainty from position uncertainty
 *
 * Δp ≥ ℏ/(2Δx)
 *
 * @param positionUncertainty Δx (m)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Minimum momentum uncertainty Δp (kg⋅m/s)
 * @throws std::invalid_argument if position uncertainty is non-positive
 */
inline double momentumUncertainty(double positionUncertainty, double hbar = constants::HBAR) {
    if (positionUncertainty <= 0.0) {
        throw std::invalid_argument("Position uncertainty must be positive");
    }
    return hbar / (2.0 * positionUncertainty);
}

/**
 * @brief Calculate minimum energy-time uncertainty product
 *
 * ΔE⋅Δt ≥ ℏ/2
 *
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Minimum uncertainty product (J⋅s)
 */
inline double energyTimeUncertainty(double hbar = constants::HBAR) {
    return hbar / 2.0;
}

/**
 * @brief Calculate minimum time uncertainty from energy uncertainty
 *
 * Δt ≥ ℏ/(2ΔE)
 *
 * @param energyUncertainty ΔE (J)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Minimum time uncertainty Δt (s)
 * @throws std::invalid_argument if energy uncertainty is non-positive
 */
inline double timeUncertainty(double energyUncertainty, double hbar = constants::HBAR) {
    if (energyUncertainty <= 0.0) {
        throw std::invalid_argument("Energy uncertainty must be positive");
    }
    return hbar / (2.0 * energyUncertainty);
}

// ============================================================================
// BOHR MODEL
// ============================================================================

/**
 * @brief Calculate energy level in Bohr hydrogen atom
 *
 * E_n = -13.6 eV / n² = -R_∞ / n²
 *
 * @param principalQuantumNumber n (positive integer)
 * @param rydbergEnergy Rydberg energy (eV)
 * @return Energy of level n (eV, negative for bound states)
 * @throws std::invalid_argument if n <= 0
 */
inline double bohrEnergyLevel(int principalQuantumNumber,
                               double rydbergEnergy = constants::RYDBERG_ENERGY) {
    if (principalQuantumNumber <= 0) {
        throw std::invalid_argument("Principal quantum number must be positive");
    }
    return -rydbergEnergy / (principalQuantumNumber * principalQuantumNumber);
}

/**
 * @brief Calculate Bohr radius for nth orbit
 *
 * r_n = n² a₀
 *
 * where a₀ is the Bohr radius
 *
 * @param principalQuantumNumber n
 * @param bohrRadius Bohr radius a₀ (m)
 * @return Orbital radius (m)
 * @throws std::invalid_argument if n <= 0
 */
inline double bohrOrbitalRadius(int principalQuantumNumber,
                                 double bohrRadius = constants::BOHR_RADIUS) {
    if (principalQuantumNumber <= 0) {
        throw std::invalid_argument("Principal quantum number must be positive");
    }
    return principalQuantumNumber * principalQuantumNumber * bohrRadius;
}

/**
 * @brief Calculate photon energy for transition between levels
 *
 * ΔE = R_∞(1/n_f² - 1/n_i²)
 *
 * @param initialLevel Initial quantum number n_i
 * @param finalLevel Final quantum number n_f
 * @param rydbergEnergy Rydberg energy (eV)
 * @return Photon energy (eV, positive for emission)
 * @throws std::invalid_argument if levels are non-positive
 */
inline double hydrogenTransitionEnergy(int initialLevel, int finalLevel,
                                        double rydbergEnergy = constants::RYDBERG_ENERGY) {
    if (initialLevel <= 0 || finalLevel <= 0) {
        throw std::invalid_argument("Quantum numbers must be positive");
    }
    double invFinal = 1.0 / (finalLevel * finalLevel);
    double invInitial = 1.0 / (initialLevel * initialLevel);
    return rydbergEnergy * (invFinal - invInitial);
}

/**
 * @brief Calculate wavelength of emitted photon in hydrogen transition
 *
 * 1/λ = R(1/n_f² - 1/n_i²)
 *
 * Rydberg formula
 *
 * @param initialLevel Initial quantum number
 * @param finalLevel Final quantum number
 * @param rydbergConstant Rydberg constant (1/m), default for hydrogen
 * @return Wavelength (m)
 * @throws std::invalid_argument if invalid transition
 */
inline double rydbergWavelength(int initialLevel, int finalLevel,
                                 double rydbergConstant = 1.0973731568160e7) {
    if (initialLevel <= finalLevel) {
        throw std::invalid_argument("Initial level must be greater than final for emission");
    }
    if (finalLevel <= 0) {
        throw std::invalid_argument("Final level must be positive");
    }

    double invFinal = 1.0 / (finalLevel * finalLevel);
    double invInitial = 1.0 / (initialLevel * initialLevel);
    double reciprocalWavelength = rydbergConstant * (invFinal - invInitial);

    return 1.0 / reciprocalWavelength;
}

/**
 * @brief Calculate ionization energy for hydrogen
 *
 * E_ion = 13.6 eV (from ground state n=1)
 *
 * @param rydbergEnergy Rydberg energy (eV)
 * @return Ionization energy (eV)
 */
inline double hydrogenIonizationEnergy(double rydbergEnergy = constants::RYDBERG_ENERGY) {
    return rydbergEnergy;
}

// ============================================================================
// PHOTOELECTRIC EFFECT
// ============================================================================

/**
 * @brief Calculate maximum kinetic energy of photoelectrons
 *
 * KE_max = hf - φ = hc/λ - φ
 *
 * Einstein's photoelectric equation
 *
 * @param photonEnergy Energy of incident photon (J)
 * @param workFunction Work function φ of material (J)
 * @return Maximum kinetic energy of ejected electron (J)
 * @throws std::invalid_argument if photon energy < work function
 */
inline double photoelectronKineticEnergy(double photonEnergy, double workFunction) {
    if (photonEnergy < workFunction) {
        throw std::invalid_argument("Photon energy must exceed work function");
    }
    return photonEnergy - workFunction;
}

/**
 * @brief Calculate threshold frequency for photoelectric effect
 *
 * f_0 = φ/h
 *
 * @param workFunction Work function φ (J)
 * @param planckH Planck's constant (J⋅s)
 * @return Threshold frequency (Hz)
 * @throws std::invalid_argument if work function is non-positive
 */
inline double thresholdFrequency(double workFunction, double planckH = constants::PLANCK_H) {
    if (workFunction <= 0.0) {
        throw std::invalid_argument("Work function must be positive");
    }
    return workFunction / planckH;
}

/**
 * @brief Calculate stopping potential
 *
 * eV_s = KE_max = hf - φ
 *
 * @param photonEnergy Photon energy (J)
 * @param workFunction Work function (J)
 * @param elementaryCharge Elementary charge e (C)
 * @return Stopping potential (V)
 */
inline double stoppingPotential(double photonEnergy, double workFunction,
                                double elementaryCharge = constants::ELEMENTARY_CHARGE) {
    if (photonEnergy < workFunction) {
        return 0.0; // No photoelectron emission
    }
    return (photonEnergy - workFunction) / elementaryCharge;
}

// ============================================================================
// PARTICLE IN A BOX
// ============================================================================

/**
 * @brief Calculate energy levels for particle in 1D infinite square well
 *
 * E_n = n²π²ℏ²/(2mL²) = n²h²/(8mL²)
 *
 * @param quantumNumber n (positive integer)
 * @param mass Particle mass (kg)
 * @param boxLength Length of box L (m)
 * @param planckH Planck's constant (J⋅s)
 * @return Energy level E_n (J)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double particleInBoxEnergy(int quantumNumber, double mass, double boxLength,
                                   double planckH = constants::PLANCK_H) {
    if (quantumNumber <= 0) {
        throw std::invalid_argument("Quantum number must be positive");
    }
    if (mass <= 0.0 || boxLength <= 0.0) {
        throw std::invalid_argument("Mass and box length must be positive");
    }

    double n2 = quantumNumber * quantumNumber;
    double L2 = boxLength * boxLength;
    return (n2 * planckH * planckH) / (8.0 * mass * L2);
}

/**
 * @brief Calculate zero-point energy for particle in box
 *
 * E_0 = h²/(8mL²)
 *
 * Minimum energy (n=1)
 *
 * @param mass Particle mass (kg)
 * @param boxLength Length of box L (m)
 * @param planckH Planck's constant (J⋅s)
 * @return Zero-point energy (J)
 */
inline double zeroPointEnergy(double mass, double boxLength,
                               double planckH = constants::PLANCK_H) {
    return particleInBoxEnergy(1, mass, boxLength, planckH);
}

// ============================================================================
// QUANTUM HARMONIC OSCILLATOR
// ============================================================================

/**
 * @brief Calculate energy levels of quantum harmonic oscillator
 *
 * E_n = ℏω(n + 1/2)
 *
 * @param quantumNumber n (non-negative integer)
 * @param angularFrequency Angular frequency ω (rad/s)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Energy level E_n (J)
 * @throws std::invalid_argument if n < 0
 */
inline double harmonicOscillatorEnergy(int quantumNumber, double angularFrequency,
                                        double hbar = constants::HBAR) {
    if (quantumNumber < 0) {
        throw std::invalid_argument("Quantum number must be non-negative");
    }
    return hbar * angularFrequency * (quantumNumber + 0.5);
}

/**
 * @brief Calculate zero-point energy of harmonic oscillator
 *
 * E_0 = ℏω/2
 *
 * @param angularFrequency Angular frequency ω (rad/s)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Zero-point energy (J)
 */
inline double harmonicOscillatorZeroPoint(double angularFrequency,
                                           double hbar = constants::HBAR) {
    return 0.5 * hbar * angularFrequency;
}

// ============================================================================
// QUANTUM TUNNELING
// ============================================================================

/**
 * @brief Calculate tunneling probability (WKB approximation, rectangular barrier)
 *
 * T ≈ exp(-2κL)
 * where κ = √(2m(V-E))/ℏ
 *
 * @param energy Particle energy E (J)
 * @param barrierHeight Barrier potential V (J)
 * @param barrierWidth Barrier width L (m)
 * @param mass Particle mass (kg)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Tunneling probability (0 to 1)
 * @throws std::invalid_argument if E >= V (classical over-barrier motion)
 */
inline double tunnelingProbability(double energy, double barrierHeight,
                                    double barrierWidth, double mass,
                                    double hbar = constants::HBAR) {
    if (energy >= barrierHeight) {
        throw std::invalid_argument("Energy must be less than barrier for tunneling");
    }
    if (mass <= 0.0 || barrierWidth <= 0.0) {
        throw std::invalid_argument("Mass and barrier width must be positive");
    }

    double kappa = std::sqrt(2.0 * mass * (barrierHeight - energy)) / hbar;
    return std::exp(-2.0 * kappa * barrierWidth);
}

/**
 * @brief Calculate decay constant κ for exponential decay in barrier
 *
 * κ = √[2m(V-E)]/ℏ
 *
 * @param energy Particle energy E (J)
 * @param barrierHeight Barrier potential V (J)
 * @param mass Particle mass (kg)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Decay constant κ (1/m)
 */
inline double tunnelingDecayConstant(double energy, double barrierHeight, double mass,
                                      double hbar = constants::HBAR) {
    if (energy >= barrierHeight) {
        throw std::invalid_argument("Energy must be less than barrier");
    }
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }

    return std::sqrt(2.0 * mass * (barrierHeight - energy)) / hbar;
}

// ============================================================================
// QUANTUM ANGULAR MOMENTUM
// ============================================================================

/**
 * @brief Calculate magnitude of orbital angular momentum
 *
 * L = √[l(l+1)] ℏ
 *
 * @param orbitalQuantumNumber l (non-negative integer)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Angular momentum magnitude (J⋅s)
 * @throws std::invalid_argument if l < 0
 */
inline double orbitalAngularMomentum(int orbitalQuantumNumber, double hbar = constants::HBAR) {
    if (orbitalQuantumNumber < 0) {
        throw std::invalid_argument("Orbital quantum number must be non-negative");
    }
    return hbar * std::sqrt(orbitalQuantumNumber * (orbitalQuantumNumber + 1.0));
}

/**
 * @brief Calculate z-component of angular momentum
 *
 * L_z = m_l ℏ
 *
 * @param magneticQuantumNumber m_l (integer, |m_l| ≤ l)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Z-component of angular momentum (J⋅s)
 */
inline double angularMomentumZComponent(int magneticQuantumNumber, double hbar = constants::HBAR) {
    return magneticQuantumNumber * hbar;
}

/**
 * @brief Calculate spin angular momentum magnitude
 *
 * S = √[s(s+1)] ℏ
 *
 * For electron: s = 1/2, so S = (√3/2)ℏ
 *
 * @param spinQuantumNumber s (typically 1/2 for fermions)
 * @param hbar Reduced Planck's constant (J⋅s)
 * @return Spin angular momentum magnitude (J⋅s)
 */
inline double spinAngularMomentum(double spinQuantumNumber, double hbar = constants::HBAR) {
    return hbar * std::sqrt(spinQuantumNumber * (spinQuantumNumber + 1.0));
}

} // namespace quantum_basics
} // namespace physics

#endif // PHYSICS_QUANTUM_BASICS_HPP
