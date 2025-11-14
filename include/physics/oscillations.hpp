#ifndef PHYSICS_OSCILLATIONS_HPP
#define PHYSICS_OSCILLATIONS_HPP

#include <cmath>
#include <stdexcept>
#include <complex>

/**
 * @file oscillations.hpp
 * @brief Advanced oscillation theory including damped, forced, and coupled oscillators
 *
 * This module implements:
 * - Damped harmonic oscillations (underdamped, critically damped, overdamped)
 * - Forced oscillations and resonance
 * - Quality factor and energy dissipation
 * - Electrical oscillations (LC, RLC circuits)
 * - Coupled oscillators and normal modes
 * - Transient response
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace oscillations {

// ============================================================================
// DAMPED HARMONIC OSCILLATIONS
// ============================================================================

/**
 * @brief Calculate damping coefficient from damping ratio
 *
 * γ = 2ζω₀
 *
 * where ζ is damping ratio and ω₀ is natural frequency
 *
 * @param dampingRatio Damping ratio ζ (dimensionless)
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @return Damping coefficient γ (1/s)
 */
inline double dampingCoefficient(double dampingRatio, double naturalFrequency) {
    return 2.0 * dampingRatio * naturalFrequency;
}

/**
 * @brief Calculate damping ratio from system parameters
 *
 * ζ = c/(2√(km)) = c/(2mω₀)
 *
 * @param dampingConstant Damping constant c (kg/s)
 * @param mass Mass (kg)
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @return Damping ratio ζ (dimensionless)
 * @throws std::invalid_argument if mass or frequency is non-positive
 */
inline double dampingRatio(double dampingConstant, double mass, double naturalFrequency) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (naturalFrequency <= 0.0) {
        throw std::invalid_argument("Natural frequency must be positive");
    }
    return dampingConstant / (2.0 * mass * naturalFrequency);
}

/**
 * @brief Calculate damped angular frequency
 *
 * ω_d = ω₀√(1 - ζ²)
 *
 * Valid for underdamped case (ζ < 1)
 *
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @param dampingRatio Damping ratio ζ
 * @return Damped frequency ω_d (rad/s)
 * @throws std::invalid_argument if damping ratio >= 1 (not underdamped)
 */
inline double dampedFrequency(double naturalFrequency, double dampingRatio) {
    if (dampingRatio >= 1.0) {
        throw std::invalid_argument("Damping ratio must be < 1 for underdamped oscillation");
    }
    return naturalFrequency * std::sqrt(1.0 - dampingRatio * dampingRatio);
}

/**
 * @brief Calculate amplitude decay envelope for underdamped oscillator
 *
 * A(t) = A₀ e^(-ζω₀t)
 *
 * @param initialAmplitude Initial amplitude A₀ (m)
 * @param dampingRatio Damping ratio ζ
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @param time Time (s)
 * @return Amplitude at time t (m)
 * @throws std::invalid_argument if time is negative
 */
inline double dampedAmplitude(double initialAmplitude, double dampingRatio,
                              double naturalFrequency, double time) {
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    return initialAmplitude * std::exp(-dampingRatio * naturalFrequency * time);
}

/**
 * @brief Calculate displacement for underdamped oscillator
 *
 * x(t) = A₀ e^(-ζω₀t) cos(ω_d t + φ)
 *
 * @param amplitude Initial amplitude (m)
 * @param dampingRatio Damping ratio ζ
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @param time Time (s)
 * @param phase Phase angle φ (radians)
 * @return Displacement (m)
 */
inline double underdampedDisplacement(double amplitude, double dampingRatio,
                                      double naturalFrequency, double time, double phase = 0.0) {
    double dampedFreq = dampedFrequency(naturalFrequency, dampingRatio);
    double envelope = dampedAmplitude(amplitude, dampingRatio, naturalFrequency, time);
    return envelope * std::cos(dampedFreq * time + phase);
}

/**
 * @brief Calculate displacement for critically damped oscillator
 *
 * x(t) = (A + Bt) e^(-ω₀t)
 *
 * For initial conditions x(0) = x₀, v(0) = v₀:
 * A = x₀, B = v₀ + ω₀x₀
 *
 * @param initialDisplacement x₀ (m)
 * @param initialVelocity v₀ (m/s)
 * @param naturalFrequency ω₀ (rad/s)
 * @param time Time (s)
 * @return Displacement (m)
 * @throws std::invalid_argument if time is negative
 */
inline double criticallyDampedDisplacement(double initialDisplacement, double initialVelocity,
                                           double naturalFrequency, double time) {
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    double A = initialDisplacement;
    double B = initialVelocity + naturalFrequency * initialDisplacement;
    return (A + B * time) * std::exp(-naturalFrequency * time);
}

/**
 * @brief Calculate relaxation time (time constant) for damped oscillator
 *
 * τ = 1/(ζω₀)
 *
 * Time for amplitude to decay to 1/e of initial value
 *
 * @param dampingRatio Damping ratio ζ
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @return Relaxation time τ (s)
 * @throws std::invalid_argument if product is zero
 */
inline double relaxationTime(double dampingRatio, double naturalFrequency) {
    double product = dampingRatio * naturalFrequency;
    if (std::abs(product) < 1e-15) {
        throw std::invalid_argument("Damping ratio times frequency must be non-zero");
    }
    return 1.0 / product;
}

/**
 * @brief Calculate logarithmic decrement
 *
 * δ = ln(A_n/A_{n+1}) = 2πζ/√(1-ζ²)
 *
 * Ratio of successive peak amplitudes
 *
 * @param dampingRatio Damping ratio ζ
 * @return Logarithmic decrement δ
 * @throws std::invalid_argument if damping ratio >= 1
 */
inline double logarithmicDecrement(double dampingRatio) {
    if (dampingRatio >= 1.0) {
        throw std::invalid_argument("Damping ratio must be < 1");
    }
    return (2.0 * M_PI * dampingRatio) / std::sqrt(1.0 - dampingRatio * dampingRatio);
}

// ============================================================================
// QUALITY FACTOR AND ENERGY
// ============================================================================

/**
 * @brief Calculate quality factor
 *
 * Q = 1/(2ζ) = ω₀/(2γ)
 *
 * Measure of oscillator's underdamping
 *
 * @param dampingRatio Damping ratio ζ
 * @return Quality factor Q
 * @throws std::invalid_argument if damping ratio is zero
 */
inline double qualityFactor(double dampingRatio) {
    if (std::abs(dampingRatio) < 1e-15) {
        throw std::invalid_argument("Damping ratio must be non-zero");
    }
    return 1.0 / (2.0 * dampingRatio);
}

/**
 * @brief Calculate damping ratio from quality factor
 *
 * ζ = 1/(2Q)
 *
 * @param qualityFactor Quality factor Q
 * @return Damping ratio ζ
 * @throws std::invalid_argument if Q is zero
 */
inline double dampingRatioFromQ(double qualityFactor) {
    if (std::abs(qualityFactor) < 1e-15) {
        throw std::invalid_argument("Quality factor must be non-zero");
    }
    return 1.0 / (2.0 * qualityFactor);
}

/**
 * @brief Calculate energy decay in damped oscillator
 *
 * E(t) = E₀ e^(-2ζω₀t) = E₀ e^(-t/τ_E)
 *
 * where τ_E = 1/(2ζω₀) is energy relaxation time
 *
 * @param initialEnergy Initial energy E₀ (J)
 * @param dampingRatio Damping ratio ζ
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @param time Time (s)
 * @return Energy at time t (J)
 * @throws std::invalid_argument if time is negative
 */
inline double energyDecay(double initialEnergy, double dampingRatio,
                         double naturalFrequency, double time) {
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    return initialEnergy * std::exp(-2.0 * dampingRatio * naturalFrequency * time);
}

/**
 * @brief Calculate power dissipated in damped oscillator
 *
 * P = -dE/dt = 2ζω₀E
 *
 * @param energy Current energy (J)
 * @param dampingRatio Damping ratio ζ
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @return Power dissipated (W)
 */
inline double powerDissipated(double energy, double dampingRatio, double naturalFrequency) {
    return 2.0 * dampingRatio * naturalFrequency * energy;
}

// ============================================================================
// FORCED OSCILLATIONS AND RESONANCE
// ============================================================================

/**
 * @brief Calculate steady-state amplitude for forced oscillator
 *
 * A(ω) = F₀/(m√[(ω₀² - ω²)² + (2ζω₀ω)²])
 *
 * @param drivingForceAmplitude F₀ (N)
 * @param mass Mass (kg)
 * @param naturalFrequency ω₀ (rad/s)
 * @param drivingFrequency ω (rad/s)
 * @param dampingRatio ζ
 * @return Amplitude (m)
 * @throws std::invalid_argument if mass is non-positive
 */
inline double forcedAmplitude(double drivingForceAmplitude, double mass,
                             double naturalFrequency, double drivingFrequency,
                             double dampingRatio) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }

    double omega0_sq = naturalFrequency * naturalFrequency;
    double omega_sq = drivingFrequency * drivingFrequency;
    double dampingTerm = 2.0 * dampingRatio * naturalFrequency * drivingFrequency;

    double denomSq = (omega0_sq - omega_sq) * (omega0_sq - omega_sq) + dampingTerm * dampingTerm;
    return drivingForceAmplitude / (mass * std::sqrt(denomSq));
}

/**
 * @brief Calculate resonance frequency for amplitude
 *
 * ω_r = ω₀√(1 - 2ζ²)
 *
 * Frequency at which amplitude is maximum (valid for ζ < 1/√2)
 *
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @param dampingRatio Damping ratio ζ
 * @return Resonance frequency (rad/s)
 * @throws std::invalid_argument if ζ >= 1/√2
 */
inline double resonanceFrequency(double naturalFrequency, double dampingRatio) {
    double criticalValue = 1.0 / std::sqrt(2.0);
    if (dampingRatio >= criticalValue) {
        throw std::invalid_argument("Damping ratio must be < 1/√2 for amplitude resonance");
    }
    return naturalFrequency * std::sqrt(1.0 - 2.0 * dampingRatio * dampingRatio);
}

/**
 * @brief Calculate maximum amplitude at resonance
 *
 * A_max = F₀/(2mζω₀²) = F₀Q/(mω₀²)
 *
 * @param drivingForceAmplitude F₀ (N)
 * @param mass Mass (kg)
 * @param naturalFrequency ω₀ (rad/s)
 * @param dampingRatio ζ
 * @return Maximum amplitude (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double maxAmplitudeAtResonance(double drivingForceAmplitude, double mass,
                                      double naturalFrequency, double dampingRatio) {
    if (mass <= 0.0 || naturalFrequency <= 0.0 || dampingRatio <= 0.0) {
        throw std::invalid_argument("Mass, frequency, and damping ratio must be positive");
    }
    return drivingForceAmplitude / (2.0 * mass * dampingRatio * naturalFrequency * naturalFrequency);
}

/**
 * @brief Calculate phase lag between driving force and displacement
 *
 * tan(φ) = 2ζω₀ω/(ω₀² - ω²)
 *
 * @param naturalFrequency ω₀ (rad/s)
 * @param drivingFrequency ω (rad/s)
 * @param dampingRatio ζ
 * @return Phase lag φ (radians, 0 to π)
 */
inline double phaseLag(double naturalFrequency, double drivingFrequency, double dampingRatio) {
    double numerator = 2.0 * dampingRatio * naturalFrequency * drivingFrequency;
    double denominator = naturalFrequency * naturalFrequency - drivingFrequency * drivingFrequency;
    return std::atan2(numerator, denominator);
}

/**
 * @brief Calculate bandwidth (full width at half maximum)
 *
 * Δω = 2ζω₀ = ω₀/Q
 *
 * Frequency range where power is above half maximum
 *
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @param dampingRatio Damping ratio ζ
 * @return Bandwidth Δω (rad/s)
 */
inline double bandwidth(double naturalFrequency, double dampingRatio) {
    return 2.0 * dampingRatio * naturalFrequency;
}

/**
 * @brief Calculate average power absorbed at steady state
 *
 * ⟨P⟩ = F₀²ω²ζω₀/(m[(ω₀² - ω²)² + (2ζω₀ω)²])
 *
 * @param drivingForceAmplitude F₀ (N)
 * @param mass Mass (kg)
 * @param naturalFrequency ω₀ (rad/s)
 * @param drivingFrequency ω (rad/s)
 * @param dampingRatio ζ
 * @return Average power (W)
 */
inline double averagePowerAbsorbed(double drivingForceAmplitude, double mass,
                                   double naturalFrequency, double drivingFrequency,
                                   double dampingRatio) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }

    double omega0_sq = naturalFrequency * naturalFrequency;
    double omega_sq = drivingFrequency * drivingFrequency;
    double dampingTerm = 2.0 * dampingRatio * naturalFrequency * drivingFrequency;

    double denominator = (omega0_sq - omega_sq) * (omega0_sq - omega_sq) + dampingTerm * dampingTerm;
    double numerator = drivingForceAmplitude * drivingForceAmplitude * omega_sq * dampingRatio * naturalFrequency;

    return numerator / (mass * denominator);
}

// ============================================================================
// ELECTRICAL OSCILLATIONS (RLC CIRCUITS)
// ============================================================================

/**
 * @brief Calculate natural frequency of LC circuit
 *
 * ω₀ = 1/√(LC)
 *
 * @param inductance Inductance L (H)
 * @param capacitance Capacitance C (F)
 * @return Natural frequency (rad/s)
 * @throws std::invalid_argument if L or C is non-positive
 */
inline double lcNaturalFrequency(double inductance, double capacitance) {
    if (inductance <= 0.0 || capacitance <= 0.0) {
        throw std::invalid_argument("Inductance and capacitance must be positive");
    }
    return 1.0 / std::sqrt(inductance * capacitance);
}

/**
 * @brief Calculate damping ratio for RLC circuit
 *
 * ζ = R/(2√(L/C)) = R/(2Z₀)
 *
 * where Z₀ = √(L/C) is characteristic impedance
 *
 * @param resistance Resistance R (Ω)
 * @param inductance Inductance L (H)
 * @param capacitance Capacitance C (F)
 * @return Damping ratio ζ
 * @throws std::invalid_argument if L or C is non-positive
 */
inline double rlcDampingRatio(double resistance, double inductance, double capacitance) {
    if (inductance <= 0.0 || capacitance <= 0.0) {
        throw std::invalid_argument("Inductance and capacitance must be positive");
    }
    double characteristicImpedance = std::sqrt(inductance / capacitance);
    return resistance / (2.0 * characteristicImpedance);
}

/**
 * @brief Calculate characteristic impedance of LC circuit
 *
 * Z₀ = √(L/C)
 *
 * @param inductance Inductance L (H)
 * @param capacitance Capacitance C (F)
 * @return Characteristic impedance (Ω)
 * @throws std::invalid_argument if L or C is non-positive
 */
inline double characteristicImpedance(double inductance, double capacitance) {
    if (inductance <= 0.0 || capacitance <= 0.0) {
        throw std::invalid_argument("Inductance and capacitance must be positive");
    }
    return std::sqrt(inductance / capacitance);
}

/**
 * @brief Calculate quality factor of RLC circuit
 *
 * Q = ω₀L/R = 1/(ω₀RC) = (1/R)√(L/C)
 *
 * @param resistance Resistance R (Ω)
 * @param inductance Inductance L (H)
 * @param capacitance Capacitance C (F)
 * @return Quality factor Q
 * @throws std::invalid_argument if parameters are invalid
 */
inline double rlcQualityFactor(double resistance, double inductance, double capacitance) {
    if (resistance <= 0.0 || inductance <= 0.0 || capacitance <= 0.0) {
        throw std::invalid_argument("All circuit parameters must be positive");
    }
    return (1.0 / resistance) * std::sqrt(inductance / capacitance);
}

/**
 * @brief Calculate energy stored in LC oscillator
 *
 * E = (1/2)LI₀² = (1/2)CV₀²
 *
 * Total energy (constant in ideal LC circuit)
 *
 * @param inductance Inductance L (H)
 * @param currentAmplitude Peak current I₀ (A)
 * @return Energy (J)
 * @throws std::invalid_argument if inductance is non-positive
 */
inline double lcEnergy(double inductance, double currentAmplitude) {
    if (inductance <= 0.0) {
        throw std::invalid_argument("Inductance must be positive");
    }
    return 0.5 * inductance * currentAmplitude * currentAmplitude;
}

/**
 * @brief Calculate impedance of RLC circuit at frequency
 *
 * |Z(ω)| = √[R² + (ωL - 1/(ωC))²]
 *
 * @param resistance R (Ω)
 * @param inductance L (H)
 * @param capacitance C (F)
 * @param frequency ω (rad/s)
 * @return Impedance magnitude (Ω)
 * @throws std::invalid_argument if frequency is zero
 */
inline double rlcImpedance(double resistance, double inductance, double capacitance, double frequency) {
    if (std::abs(frequency) < 1e-15) {
        throw std::invalid_argument("Frequency must be non-zero");
    }

    double reactance = frequency * inductance - 1.0 / (frequency * capacitance);
    return std::sqrt(resistance * resistance + reactance * reactance);
}

// ============================================================================
// COUPLED OSCILLATORS
// ============================================================================

/**
 * @brief Calculate normal mode frequencies for two coupled oscillators
 *
 * ω± = √(ω₀² ± κω₀²)
 *
 * where κ is coupling strength
 *
 * @param naturalFrequency Natural frequency of uncoupled oscillator ω₀ (rad/s)
 * @param couplingStrength Coupling strength κ (dimensionless)
 * @return [ω₋, ω₊] Lower and upper mode frequencies (rad/s)
 * @throws std::invalid_argument if frequencies would be imaginary
 */
inline std::array<double, 2> normalModeFrequencies(double naturalFrequency, double couplingStrength) {
    if (naturalFrequency <= 0.0) {
        throw std::invalid_argument("Natural frequency must be positive");
    }
    if (couplingStrength < -1.0) {
        throw std::invalid_argument("Coupling strength too negative (frequencies would be imaginary)");
    }

    double omega0_sq = naturalFrequency * naturalFrequency;
    double omegaMinus = std::sqrt(omega0_sq * (1.0 - couplingStrength));
    double omegaPlus = std::sqrt(omega0_sq * (1.0 + couplingStrength));

    return {omegaMinus, omegaPlus};
}

/**
 * @brief Calculate beat frequency for weakly coupled oscillators
 *
 * ω_beat = |ω₊ - ω₋| ≈ κω₀ (for weak coupling)
 *
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @param couplingStrength Coupling strength κ
 * @return Beat frequency (rad/s)
 */
inline double beatFrequency(double naturalFrequency, double couplingStrength) {
    return std::abs(couplingStrength * naturalFrequency);
}

/**
 * @brief Calculate energy transfer time for coupled oscillators
 *
 * τ_transfer = π/ω_beat = π/(κω₀)
 *
 * Time for energy to transfer from one oscillator to the other
 *
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @param couplingStrength Coupling strength κ
 * @return Transfer time (s)
 * @throws std::invalid_argument if coupling is zero
 */
inline double energyTransferTime(double naturalFrequency, double couplingStrength) {
    if (std::abs(couplingStrength) < 1e-15) {
        throw std::invalid_argument("Coupling strength must be non-zero");
    }
    return M_PI / (std::abs(couplingStrength) * naturalFrequency);
}

// ============================================================================
// TRANSIENT RESPONSE
// ============================================================================

/**
 * @brief Calculate time to reach steady state (approximately)
 *
 * t_ss ≈ 5τ = 5/(ζω₀)
 *
 * Time for transients to decay to ~1% of initial value
 *
 * @param dampingRatio Damping ratio ζ
 * @param naturalFrequency Natural frequency ω₀ (rad/s)
 * @return Settling time (s)
 */
inline double settlingTime(double dampingRatio, double naturalFrequency) {
    return 5.0 / (dampingRatio * naturalFrequency);
}

/**
 * @brief Calculate number of oscillations before amplitude decays to 1/e
 *
 * N = Q/π
 *
 * @param qualityFactor Quality factor Q
 * @return Number of cycles
 */
inline double cyclesBeforeDecay(double qualityFactor) {
    return qualityFactor / M_PI;
}

} // namespace oscillations
} // namespace physics

#endif // PHYSICS_OSCILLATIONS_HPP
