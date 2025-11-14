#ifndef PHYSICS_ELECTROMAGNETIC_INDUCTION_HPP
#define PHYSICS_ELECTROMAGNETIC_INDUCTION_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file electromagnetic_induction.hpp
 * @brief Comprehensive implementation of electromagnetic induction
 *
 * This module implements:
 * - Faraday's law of electromagnetic induction
 * - Lenz's law
 * - Induced EMF and current
 * - Motional EMF
 * - Self-inductance and mutual inductance
 * - Energy in inductors
 * - Transformers
 * - Electric motors and generators
 * - Eddy currents
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace electromagnetic_induction {

/**
 * @namespace constants
 * @brief Physical constants for electromagnetic induction
 */
namespace constants {
    constexpr double MU_0 = 4.0 * M_PI * 1e-7;  // Permeability of free space (H/m)
}

// ============================================================================
// FARADAY'S LAW AND INDUCED EMF
// ============================================================================

/**
 * @brief Calculate induced EMF from rate of flux change (Faraday's law)
 *
 * ε = -N × dΦ/dt
 *
 * The negative sign represents Lenz's law (direction of induced EMF opposes change)
 *
 * @param numberOfTurns Number of turns in coil
 * @param fluxChange Change in magnetic flux (Wb)
 * @param timeInterval Time interval (s)
 * @return Magnitude of induced EMF (V)
 * @throws std::invalid_argument if time interval is zero or number of turns is negative
 */
inline double inducedEMF(int numberOfTurns, double fluxChange, double timeInterval) {
    if (numberOfTurns < 0) {
        throw std::invalid_argument("Number of turns cannot be negative");
    }
    if (std::abs(timeInterval) < 1e-15) {
        throw std::invalid_argument("Time interval must be non-zero");
    }
    return std::abs(numberOfTurns * fluxChange / timeInterval);
}

/**
 * @brief Calculate induced EMF from changing magnetic field
 *
 * ε = N × A × dB/dt
 *
 * When flux change is due to changing magnetic field through fixed area
 *
 * @param numberOfTurns Number of turns in coil
 * @param area Area of coil (m²)
 * @param fieldChange Change in magnetic field (T)
 * @param timeInterval Time interval (s)
 * @return Magnitude of induced EMF (V)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double inducedEMFFromFieldChange(int numberOfTurns, double area,
                                        double fieldChange, double timeInterval) {
    if (numberOfTurns < 0) {
        throw std::invalid_argument("Number of turns cannot be negative");
    }
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    if (std::abs(timeInterval) < 1e-15) {
        throw std::invalid_argument("Time interval must be non-zero");
    }
    return std::abs(numberOfTurns * area * fieldChange / timeInterval);
}

/**
 * @brief Calculate magnetic flux through surface
 *
 * Φ = B × A × cos(θ)
 *
 * @param magneticField Magnetic field strength (T)
 * @param area Area (m²)
 * @param angle Angle between field and area normal (radians)
 * @return Magnetic flux (Wb)
 * @throws std::invalid_argument if area is negative
 */
inline double magneticFlux(double magneticField, double area, double angle = 0.0) {
    if (area < 0.0) {
        throw std::invalid_argument("Area cannot be negative");
    }
    return magneticField * area * std::cos(angle);
}

/**
 * @brief Calculate flux linkage
 *
 * Λ = N × Φ
 *
 * Total flux linkage through coil with N turns
 *
 * @param numberOfTurns Number of turns
 * @param flux Magnetic flux through one turn (Wb)
 * @return Flux linkage (Wb-turns)
 * @throws std::invalid_argument if number of turns is negative
 */
inline double fluxLinkage(int numberOfTurns, double flux) {
    if (numberOfTurns < 0) {
        throw std::invalid_argument("Number of turns cannot be negative");
    }
    return numberOfTurns * flux;
}

// ============================================================================
// MOTIONAL EMF
// ============================================================================

/**
 * @brief Calculate motional EMF in rod moving through magnetic field
 *
 * ε = B × L × v
 *
 * When rod of length L moves with velocity v perpendicular to field B
 *
 * @param magneticField Magnetic field strength (T)
 * @param length Length of rod (m)
 * @param velocity Velocity of rod (m/s)
 * @return Induced EMF (V)
 * @throws std::invalid_argument if length is non-positive
 */
inline double motionalEMF(double magneticField, double length, double velocity) {
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    return magneticField * length * velocity;
}

/**
 * @brief Calculate motional EMF with angle
 *
 * ε = B × L × v × sin(θ)
 *
 * General case where velocity makes angle θ with magnetic field
 *
 * @param magneticField Magnetic field strength (T)
 * @param length Length of conductor (m)
 * @param velocity Velocity (m/s)
 * @param angle Angle between velocity and field (radians)
 * @return Induced EMF (V)
 * @throws std::invalid_argument if length is non-positive
 */
inline double motionalEMFWithAngle(double magneticField, double length,
                                   double velocity, double angle) {
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    return magneticField * length * velocity * std::sin(angle);
}

/**
 * @brief Calculate EMF in rotating coil
 *
 * ε = N × B × A × ω × sin(ωt)
 *
 * For coil rotating in magnetic field (AC generator principle)
 *
 * @param numberOfTurns Number of turns
 * @param magneticField Magnetic field (T)
 * @param area Area of coil (m²)
 * @param angularVelocity Angular velocity (rad/s)
 * @param time Time (s)
 * @return Instantaneous EMF (V)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double rotatingCoilEMF(int numberOfTurns, double magneticField,
                              double area, double angularVelocity, double time) {
    if (numberOfTurns < 0) {
        throw std::invalid_argument("Number of turns cannot be negative");
    }
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    return numberOfTurns * magneticField * area * angularVelocity *
           std::sin(angularVelocity * time);
}

/**
 * @brief Calculate maximum EMF in rotating coil
 *
 * ε_max = N × B × A × ω
 *
 * @param numberOfTurns Number of turns
 * @param magneticField Magnetic field (T)
 * @param area Area of coil (m²)
 * @param angularVelocity Angular velocity (rad/s)
 * @return Maximum EMF (V)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double maxRotatingCoilEMF(int numberOfTurns, double magneticField,
                                 double area, double angularVelocity) {
    if (numberOfTurns < 0) {
        throw std::invalid_argument("Number of turns cannot be negative");
    }
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    return numberOfTurns * magneticField * area * angularVelocity;
}

// ============================================================================
// INDUCED CURRENT
// ============================================================================

/**
 * @brief Calculate induced current from induced EMF
 *
 * I = ε / R
 *
 * @param inducedEMF Induced EMF (V)
 * @param resistance Total circuit resistance (Ω)
 * @return Induced current (A)
 * @throws std::invalid_argument if resistance is non-positive
 */
inline double inducedCurrent(double inducedEMF, double resistance) {
    if (resistance <= 0.0) {
        throw std::invalid_argument("Resistance must be positive");
    }
    return inducedEMF / resistance;
}

/**
 * @brief Calculate total charge flow during flux change
 *
 * Q = N × ΔΦ / R
 *
 * Total charge that flows when flux changes
 *
 * @param numberOfTurns Number of turns
 * @param fluxChange Change in flux (Wb)
 * @param resistance Circuit resistance (Ω)
 * @return Total charge (C)
 * @throws std::invalid_argument if resistance is non-positive
 */
inline double totalChargeFlow(int numberOfTurns, double fluxChange, double resistance) {
    if (numberOfTurns < 0) {
        throw std::invalid_argument("Number of turns cannot be negative");
    }
    if (resistance <= 0.0) {
        throw std::invalid_argument("Resistance must be positive");
    }
    return std::abs(numberOfTurns * fluxChange / resistance);
}

// ============================================================================
// SELF-INDUCTANCE
// ============================================================================

/**
 * @brief Calculate self-inductance from flux linkage
 *
 * L = N × Φ / I = Λ / I
 *
 * @param fluxLinkage Flux linkage (Wb-turns)
 * @param current Current (A)
 * @return Self-inductance (H)
 * @throws std::invalid_argument if current is zero
 */
inline double selfInductance(double fluxLinkage, double current) {
    if (std::abs(current) < 1e-15) {
        throw std::invalid_argument("Current must be non-zero");
    }
    return fluxLinkage / current;
}

/**
 * @brief Calculate self-inductance of solenoid
 *
 * L = μ₀ × n² × A × l
 *
 * where n is turns per unit length
 *
 * @param turnsPerLength Turns per unit length (turns/m)
 * @param area Cross-sectional area (m²)
 * @param length Length of solenoid (m)
 * @param mu0 Permeability of free space (H/m)
 * @return Self-inductance (H)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double solenoidInductance(double turnsPerLength, double area, double length,
                                 double mu0 = constants::MU_0) {
    if (turnsPerLength <= 0.0) {
        throw std::invalid_argument("Turns per length must be positive");
    }
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    return mu0 * turnsPerLength * turnsPerLength * area * length;
}

/**
 * @brief Calculate self-inductance of solenoid from total turns
 *
 * L = μ₀ × N² × A / l
 *
 * @param totalTurns Total number of turns
 * @param area Cross-sectional area (m²)
 * @param length Length of solenoid (m)
 * @param mu0 Permeability of free space (H/m)
 * @return Self-inductance (H)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double solenoidInductanceFromTurns(int totalTurns, double area, double length,
                                          double mu0 = constants::MU_0) {
    if (totalTurns < 0) {
        throw std::invalid_argument("Total turns cannot be negative");
    }
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    return mu0 * totalTurns * totalTurns * area / length;
}

/**
 * @brief Calculate back EMF due to self-inductance
 *
 * ε = -L × dI/dt
 *
 * Self-induced EMF opposes change in current
 *
 * @param inductance Self-inductance (H)
 * @param currentChange Change in current (A)
 * @param timeInterval Time interval (s)
 * @return Magnitude of back EMF (V)
 * @throws std::invalid_argument if inductance is negative or time is zero
 */
inline double backEMF(double inductance, double currentChange, double timeInterval) {
    if (inductance < 0.0) {
        throw std::invalid_argument("Inductance cannot be negative");
    }
    if (std::abs(timeInterval) < 1e-15) {
        throw std::invalid_argument("Time interval must be non-zero");
    }
    return std::abs(inductance * currentChange / timeInterval);
}

/**
 * @brief Calculate energy stored in inductor
 *
 * U = (1/2) × L × I²
 *
 * @param inductance Inductance (H)
 * @param current Current (A)
 * @return Energy stored (J)
 * @throws std::invalid_argument if inductance is negative
 */
inline double inductorEnergy(double inductance, double current) {
    if (inductance < 0.0) {
        throw std::invalid_argument("Inductance cannot be negative");
    }
    return 0.5 * inductance * current * current;
}

/**
 * @brief Calculate energy density in magnetic field
 *
 * u = B² / (2μ₀)
 *
 * @param magneticField Magnetic field (T)
 * @param mu0 Permeability of free space (H/m)
 * @return Energy density (J/m³)
 */
inline double magneticEnergyDensity(double magneticField, double mu0 = constants::MU_0) {
    return (magneticField * magneticField) / (2.0 * mu0);
}

// ============================================================================
// MUTUAL INDUCTANCE
// ============================================================================

/**
 * @brief Calculate mutual inductance
 *
 * M = N₂ × Φ₂₁ / I₁
 *
 * Flux through coil 2 due to current in coil 1
 *
 * @param turns2 Number of turns in second coil
 * @param flux21 Flux through coil 2 due to coil 1 (Wb)
 * @param current1 Current in first coil (A)
 * @return Mutual inductance (H)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double mutualInductance(int turns2, double flux21, double current1) {
    if (turns2 < 0) {
        throw std::invalid_argument("Number of turns cannot be negative");
    }
    if (std::abs(current1) < 1e-15) {
        throw std::invalid_argument("Current must be non-zero");
    }
    return std::abs(turns2 * flux21 / current1);
}

/**
 * @brief Calculate mutually induced EMF
 *
 * ε₂ = -M × dI₁/dt
 *
 * EMF induced in coil 2 by changing current in coil 1
 *
 * @param mutualInductance Mutual inductance (H)
 * @param currentChange Change in current in first coil (A)
 * @param timeInterval Time interval (s)
 * @return Magnitude of induced EMF (V)
 * @throws std::invalid_argument if mutual inductance is negative or time is zero
 */
inline double mutuallyInducedEMF(double mutualInductance, double currentChange,
                                 double timeInterval) {
    if (mutualInductance < 0.0) {
        throw std::invalid_argument("Mutual inductance cannot be negative");
    }
    if (std::abs(timeInterval) < 1e-15) {
        throw std::invalid_argument("Time interval must be non-zero");
    }
    return std::abs(mutualInductance * currentChange / timeInterval);
}

// ============================================================================
// MOTORS AND GENERATORS
// ============================================================================

/**
 * @brief Calculate EMF of armature (generator)
 *
 * ε = N × B × A × ω
 *
 * Average EMF for DC generator with commutator
 *
 * @param numberOfTurns Number of turns in armature
 * @param magneticField Magnetic field (T)
 * @param area Area of each turn (m²)
 * @param angularVelocity Angular velocity (rad/s)
 * @return EMF (V)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double generatorEMF(int numberOfTurns, double magneticField,
                           double area, double angularVelocity) {
    if (numberOfTurns < 0) {
        throw std::invalid_argument("Number of turns cannot be negative");
    }
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    return numberOfTurns * magneticField * area * angularVelocity;
}

/**
 * @brief Calculate mechanical power output of generator
 *
 * P_mech = τ × ω
 *
 * @param torque Torque (N⋅m)
 * @param angularVelocity Angular velocity (rad/s)
 * @return Mechanical power (W)
 */
inline double generatorMechanicalPower(double torque, double angularVelocity) {
    return torque * angularVelocity;
}

/**
 * @brief Calculate electrical power output of generator
 *
 * P_elec = ε × I
 *
 * @param emf Generated EMF (V)
 * @param current Output current (A)
 * @return Electrical power (W)
 */
inline double generatorElectricalPower(double emf, double current) {
    return emf * current;
}

/**
 * @brief Calculate energy spent in motor
 *
 * E = (V - ε_back) × I × t
 *
 * Energy converted to mechanical work
 *
 * @param appliedVoltage Applied voltage (V)
 * @param backEMF Back EMF (V)
 * @param current Current (A)
 * @param time Time (s)
 * @return Energy (J)
 * @throws std::invalid_argument if time is negative
 */
inline double motorEnergy(double appliedVoltage, double backEMF,
                         double current, double time) {
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    return (appliedVoltage - backEMF) * current * time;
}

/**
 * @brief Calculate mechanical power developed by motor
 *
 * P_mech = ε_back × I
 *
 * @param backEMF Back EMF (V)
 * @param current Current (A)
 * @return Mechanical power (W)
 */
inline double motorMechanicalPower(double backEMF, double current) {
    return backEMF * current;
}

/**
 * @brief Calculate electrical power input to motor
 *
 * P_input = V × I
 *
 * @param voltage Applied voltage (V)
 * @param current Current (A)
 * @return Input power (W)
 */
inline double motorInputPower(double voltage, double current) {
    return voltage * current;
}

/**
 * @brief Calculate motor efficiency
 *
 * η = P_mech / P_input = ε_back / V
 *
 * @param backEMF Back EMF (V)
 * @param appliedVoltage Applied voltage (V)
 * @return Efficiency (0 to 1)
 * @throws std::invalid_argument if applied voltage is zero
 */
inline double motorEfficiency(double backEMF, double appliedVoltage) {
    if (std::abs(appliedVoltage) < 1e-15) {
        throw std::invalid_argument("Applied voltage must be non-zero");
    }
    return backEMF / appliedVoltage;
}

/**
 * @brief Calculate torque developed by motor
 *
 * τ = P_mech / ω
 *
 * @param mechanicalPower Mechanical power (W)
 * @param angularVelocity Angular velocity (rad/s)
 * @return Torque (N⋅m)
 * @throws std::invalid_argument if angular velocity is zero
 */
inline double motorTorque(double mechanicalPower, double angularVelocity) {
    if (std::abs(angularVelocity) < 1e-15) {
        throw std::invalid_argument("Angular velocity must be non-zero");
    }
    return mechanicalPower / angularVelocity;
}

// ============================================================================
// TRANSFORMERS
// ============================================================================

/**
 * @brief Calculate transformer voltage ratio
 *
 * V_s / V_p = N_s / N_p
 *
 * @param primaryVoltage Primary voltage (V)
 * @param primaryTurns Primary turns
 * @param secondaryTurns Secondary turns
 * @return Secondary voltage (V)
 * @throws std::invalid_argument if primary turns is zero
 */
inline double transformerSecondaryVoltage(double primaryVoltage,
                                          int primaryTurns, int secondaryTurns) {
    if (primaryTurns == 0) {
        throw std::invalid_argument("Primary turns must be non-zero");
    }
    return primaryVoltage * static_cast<double>(secondaryTurns) / primaryTurns;
}

/**
 * @brief Calculate transformer current ratio
 *
 * I_s / I_p = N_p / N_s (for ideal transformer)
 *
 * @param primaryCurrent Primary current (A)
 * @param primaryTurns Primary turns
 * @param secondaryTurns Secondary turns
 * @return Secondary current (A)
 * @throws std::invalid_argument if secondary turns is zero
 */
inline double transformerSecondaryCurrent(double primaryCurrent,
                                          int primaryTurns, int secondaryTurns) {
    if (secondaryTurns == 0) {
        throw std::invalid_argument("Secondary turns must be non-zero");
    }
    return primaryCurrent * static_cast<double>(primaryTurns) / secondaryTurns;
}

/**
 * @brief Calculate transformer turns ratio
 *
 * n = N_s / N_p
 *
 * @param primaryTurns Primary turns
 * @param secondaryTurns Secondary turns
 * @return Turns ratio
 * @throws std::invalid_argument if primary turns is zero
 */
inline double transformerTurnsRatio(int primaryTurns, int secondaryTurns) {
    if (primaryTurns == 0) {
        throw std::invalid_argument("Primary turns must be non-zero");
    }
    return static_cast<double>(secondaryTurns) / primaryTurns;
}

/**
 * @brief Calculate transformer efficiency
 *
 * η = P_out / P_in
 *
 * @param outputPower Output power (W)
 * @param inputPower Input power (W)
 * @return Efficiency (0 to 1)
 * @throws std::invalid_argument if input power is zero
 */
inline double transformerEfficiency(double outputPower, double inputPower) {
    if (std::abs(inputPower) < 1e-15) {
        throw std::invalid_argument("Input power must be non-zero");
    }
    return outputPower / inputPower;
}

// ============================================================================
// EDDY CURRENTS
// ============================================================================

/**
 * @brief Calculate power loss due to eddy currents
 *
 * P ∝ B² × f² × t² × V
 *
 * Qualitative relationship for eddy current losses
 *
 * @param magneticField Magnetic field (T)
 * @param frequency Frequency (Hz)
 * @param thickness Conductor thickness (m)
 * @param volume Volume of conductor (m³)
 * @param proportionalityConstant Material-dependent constant
 * @return Power loss (W)
 * @throws std::invalid_argument if parameters are negative
 */
inline double eddyCurrentLoss(double magneticField, double frequency,
                              double thickness, double volume,
                              double proportionalityConstant = 1.0) {
    if (frequency < 0.0 || thickness < 0.0 || volume < 0.0) {
        throw std::invalid_argument("Parameters cannot be negative");
    }
    return proportionalityConstant * magneticField * magneticField *
           frequency * frequency * thickness * thickness * volume;
}

// ============================================================================
// TIME CONSTANT AND RL CIRCUITS
// ============================================================================

/**
 * @brief Calculate time constant of RL circuit
 *
 * τ = L / R
 *
 * @param inductance Inductance (H)
 * @param resistance Resistance (Ω)
 * @return Time constant (s)
 * @throws std::invalid_argument if resistance is non-positive
 */
inline double rlTimeConstant(double inductance, double resistance) {
    if (resistance <= 0.0) {
        throw std::invalid_argument("Resistance must be positive");
    }
    if (inductance < 0.0) {
        throw std::invalid_argument("Inductance cannot be negative");
    }
    return inductance / resistance;
}

/**
 * @brief Calculate current growth in RL circuit
 *
 * I(t) = I_max × (1 - e^(-t/τ))
 * where I_max = V/R and τ = L/R
 *
 * @param voltage Applied voltage (V)
 * @param resistance Resistance (Ω)
 * @param inductance Inductance (H)
 * @param time Time (s)
 * @return Current at time t (A)
 * @throws std::invalid_argument if resistance is non-positive or time is negative
 */
inline double rlCurrentGrowth(double voltage, double resistance,
                              double inductance, double time) {
    if (resistance <= 0.0) {
        throw std::invalid_argument("Resistance must be positive");
    }
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    double maxCurrent = voltage / resistance;
    double tau = inductance / resistance;
    return maxCurrent * (1.0 - std::exp(-time / tau));
}

/**
 * @brief Calculate current decay in RL circuit
 *
 * I(t) = I₀ × e^(-t/τ)
 *
 * @param initialCurrent Initial current (A)
 * @param resistance Resistance (Ω)
 * @param inductance Inductance (H)
 * @param time Time (s)
 * @return Current at time t (A)
 * @throws std::invalid_argument if resistance is non-positive or time is negative
 */
inline double rlCurrentDecay(double initialCurrent, double resistance,
                             double inductance, double time) {
    if (resistance <= 0.0) {
        throw std::invalid_argument("Resistance must be positive");
    }
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    double tau = inductance / resistance;
    return initialCurrent * std::exp(-time / tau);
}

} // namespace electromagnetic_induction
} // namespace physics

#endif // PHYSICS_ELECTROMAGNETIC_INDUCTION_HPP
