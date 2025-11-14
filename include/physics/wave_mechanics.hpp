#ifndef PHYSICS_WAVE_MECHANICS_HPP
#define PHYSICS_WAVE_MECHANICS_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file wave_mechanics.hpp
 * @brief Comprehensive implementation of wave phenomena, oscillations, and acoustics
 *
 * This module implements:
 * - Wave fundamentals (wavelength, frequency, velocity)
 * - Sound waves and compressional waves
 * - Newton's formula for sound velocity
 * - Laplace's correction
 * - Sound intensity and decibels
 * - Doppler effect
 * - String vibrations and standing waves
 * - Wave energy and power
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace wave_mechanics {

/**
 * @namespace constants
 * @brief Physical constants for wave mechanics calculations
 */
namespace constants {
    // Speed of sound (m/s at 20°C, sea level)
    constexpr double SOUND_SPEED_AIR_20C = 343.0;
    constexpr double SOUND_SPEED_WATER = 1482.0;
    constexpr double SOUND_SPEED_STEEL = 5960.0;
    constexpr double SOUND_SPEED_ALUMINUM = 6420.0;

    // Reference values for sound intensity
    constexpr double REFERENCE_INTENSITY = 1e-12;  // W/m² (threshold of hearing)
    constexpr double REFERENCE_PRESSURE = 2e-5;    // Pa (threshold pressure)

    // Thermodynamic constants
    constexpr double GAMMA_AIR = 1.4;              // Adiabatic index for air
    constexpr double R = 8.314;                    // Universal gas constant (J/(mol·K))

    // Standard conditions
    constexpr double STANDARD_TEMPERATURE = 273.15; // K (0°C)
    constexpr double STANDARD_PRESSURE = 101325.0;  // Pa

    // Musical frequencies (Hz)
    constexpr double A4_FREQUENCY = 440.0;         // Concert pitch A
    constexpr double MIDDLE_C_FREQUENCY = 261.63;  // Middle C (C4)
}

// ============================================================================
// WAVE FUNDAMENTALS
// ============================================================================

/**
 * @brief Calculate wavelength from velocity and frequency
 *
 * λ = v/f
 *
 * @param velocity Wave velocity (m/s)
 * @param frequency Wave frequency (Hz)
 * @return Wavelength (m)
 * @throws std::invalid_argument if frequency is non-positive
 */
inline double calculateWavelength(double velocity, double frequency) {
    if (frequency <= 0.0) {
        throw std::invalid_argument("Frequency must be positive");
    }
    return velocity / frequency;
}

/**
 * @brief Calculate wave frequency from velocity and wavelength
 *
 * f = v/λ
 *
 * @param velocity Wave velocity (m/s)
 * @param wavelength Wavelength (m)
 * @return Frequency (Hz)
 * @throws std::invalid_argument if wavelength is non-positive
 */
inline double calculateFrequency(double velocity, double wavelength) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return velocity / wavelength;
}

/**
 * @brief Calculate wave velocity from frequency and wavelength
 *
 * v = f × λ
 *
 * @param frequency Wave frequency (Hz)
 * @param wavelength Wavelength (m)
 * @return Wave velocity (m/s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double calculateWaveVelocity(double frequency, double wavelength) {
    if (frequency <= 0.0 || wavelength <= 0.0) {
        throw std::invalid_argument("Frequency and wavelength must be positive");
    }
    return frequency * wavelength;
}

/**
 * @brief Calculate angular frequency from frequency
 *
 * ω = 2πf
 *
 * @param frequency Frequency (Hz)
 * @return Angular frequency (rad/s)
 */
inline double calculateAngularFrequency(double frequency) {
    return 2.0 * M_PI * frequency;
}

/**
 * @brief Calculate wave number (spatial frequency)
 *
 * k = 2π/λ
 *
 * @param wavelength Wavelength (m)
 * @return Wave number (rad/m)
 * @throws std::invalid_argument if wavelength is non-positive
 */
inline double calculateWaveNumber(double wavelength) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return (2.0 * M_PI) / wavelength;
}

/**
 * @brief Calculate period from frequency
 *
 * T = 1/f
 *
 * @param frequency Frequency (Hz)
 * @return Period (s)
 * @throws std::invalid_argument if frequency is non-positive
 */
inline double calculatePeriod(double frequency) {
    if (frequency <= 0.0) {
        throw std::invalid_argument("Frequency must be positive");
    }
    return 1.0 / frequency;
}

// ============================================================================
// VELOCITY OF COMPRESSIONAL WAVES (SOUND)
// ============================================================================

/**
 * @brief Newton's formula for velocity of sound in a medium
 *
 * v = √(E/ρ)
 *
 * where E is elastic modulus and ρ is density
 * (Note: This gives incorrect results for gases; use Laplace correction)
 *
 * @param elasticModulus Elastic modulus (Pa)
 * @param density Density (kg/m³)
 * @return Velocity of sound (m/s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double newtonsFormulaSound(double elasticModulus, double density) {
    if (elasticModulus <= 0.0 || density <= 0.0) {
        throw std::invalid_argument("Elastic modulus and density must be positive");
    }
    return std::sqrt(elasticModulus / density);
}

/**
 * @brief Laplace's correction for velocity of sound in gases
 *
 * v = √(γP/ρ) = √(γRT/M)
 *
 * where γ is adiabatic index (ratio of specific heats)
 * This corrects Newton's formula for gases by using adiabatic (not isothermal) process
 *
 * @param adiabaticIndex Ratio of specific heats γ = Cp/Cv (dimensionless)
 * @param pressure Gas pressure (Pa)
 * @param density Gas density (kg/m³)
 * @return Velocity of sound (m/s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double laplaceFormulaSound(double adiabaticIndex, double pressure,
                                  double density) {
    if (adiabaticIndex <= 0.0 || pressure <= 0.0 || density <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return std::sqrt((adiabaticIndex * pressure) / density);
}

/**
 * @brief Calculate velocity of sound in ideal gas from temperature
 *
 * v = √(γRT/M)
 *
 * where R is gas constant, T is temperature, M is molar mass
 *
 * @param adiabaticIndex Ratio of specific heats (dimensionless)
 * @param temperature Temperature (K)
 * @param molarMass Molar mass (kg/mol)
 * @param R Universal gas constant (J/(mol·K))
 * @return Velocity of sound (m/s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double soundVelocityFromTemperature(double adiabaticIndex, double temperature,
                                           double molarMass,
                                           double R = constants::R) {
    if (adiabaticIndex <= 0.0 || temperature <= 0.0 || molarMass <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return std::sqrt((adiabaticIndex * R * temperature) / molarMass);
}

/**
 * @brief Calculate velocity of sound in air at given temperature
 *
 * v ≈ 331.3 + 0.606T (empirical formula)
 *
 * where T is temperature in Celsius
 *
 * @param temperatureCelsius Temperature in Celsius
 * @return Velocity of sound in air (m/s)
 */
inline double soundVelocityInAir(double temperatureCelsius) {
    return 331.3 + 0.606 * temperatureCelsius;
}

/**
 * @brief Calculate velocity of sound in solids (longitudinal waves)
 *
 * v = √(E/ρ)
 *
 * where E is Young's modulus
 *
 * @param youngsModulus Young's modulus (Pa)
 * @param density Density (kg/m³)
 * @return Velocity of longitudinal waves (m/s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double soundVelocityInSolid(double youngsModulus, double density) {
    return newtonsFormulaSound(youngsModulus, density);
}

/**
 * @brief Calculate velocity of sound in liquids
 *
 * v = √(K/ρ)
 *
 * where K is bulk modulus
 *
 * @param bulkModulus Bulk modulus (Pa)
 * @param density Density (kg/m³)
 * @return Velocity of sound (m/s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double soundVelocityInLiquid(double bulkModulus, double density) {
    return newtonsFormulaSound(bulkModulus, density);
}

// ============================================================================
// SOUND INTENSITY AND DECIBELS
// ============================================================================

/**
 * @brief Calculate sound intensity from power and area
 *
 * I = P/A
 *
 * @param power Sound power (W)
 * @param area Area through which sound passes (m²)
 * @return Intensity (W/m²)
 * @throws std::invalid_argument if area is non-positive
 */
inline double calculateSoundIntensity(double power, double area) {
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    return power / area;
}

/**
 * @brief Calculate sound intensity level in decibels
 *
 * β = 10 log₁₀(I/I₀)
 *
 * where I₀ = 10⁻¹² W/m² is the reference intensity (threshold of hearing)
 *
 * @param intensity Sound intensity (W/m²)
 * @param referenceIntensity Reference intensity (W/m²)
 * @return Sound intensity level (dB)
 * @throws std::invalid_argument if intensities are non-positive
 */
inline double calculateSoundLevelDecibels(double intensity,
                                          double referenceIntensity = constants::REFERENCE_INTENSITY) {
    if (intensity <= 0.0 || referenceIntensity <= 0.0) {
        throw std::invalid_argument("Intensities must be positive");
    }
    return 10.0 * std::log10(intensity / referenceIntensity);
}

/**
 * @brief Calculate intensity from decibel level
 *
 * I = I₀ × 10^(β/10)
 *
 * @param decibelLevel Sound level (dB)
 * @param referenceIntensity Reference intensity (W/m²)
 * @return Sound intensity (W/m²)
 */
inline double intensityFromDecibels(double decibelLevel,
                                    double referenceIntensity = constants::REFERENCE_INTENSITY) {
    return referenceIntensity * std::pow(10.0, decibelLevel / 10.0);
}

/**
 * @brief Calculate sound intensity from pressure amplitude
 *
 * I = (ΔP²)/(2ρv)
 *
 * where ΔP is pressure amplitude, ρ is density, v is sound velocity
 *
 * @param pressureAmplitude Pressure amplitude (Pa)
 * @param density Medium density (kg/m³)
 * @param soundVelocity Velocity of sound in medium (m/s)
 * @return Sound intensity (W/m²)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double intensityFromPressure(double pressureAmplitude, double density,
                                    double soundVelocity) {
    if (density <= 0.0 || soundVelocity <= 0.0) {
        throw std::invalid_argument("Density and sound velocity must be positive");
    }
    return (pressureAmplitude * pressureAmplitude) / (2.0 * density * soundVelocity);
}

/**
 * @brief Calculate spherical wave intensity at distance from point source
 *
 * I = P/(4πr²)
 *
 * Intensity decreases as 1/r² for spherical waves
 *
 * @param power Total power of source (W)
 * @param distance Distance from source (m)
 * @return Intensity at distance (W/m²)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double sphericalWaveIntensity(double power, double distance) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return power / (4.0 * M_PI * distance * distance);
}

// ============================================================================
// DOPPLER EFFECT
// ============================================================================

/**
 * @brief Calculate observed frequency using Doppler effect
 *
 * f' = f × (v + v_observer)/(v - v_source)
 *
 * Positive velocities: observer/source moving toward each other
 * Negative velocities: moving apart
 *
 * @param sourceFrequency Frequency emitted by source (Hz)
 * @param soundVelocity Velocity of sound in medium (m/s)
 * @param observerVelocity Velocity of observer (m/s), positive toward source
 * @param sourceVelocity Velocity of source (m/s), positive toward observer
 * @return Observed frequency (Hz)
 * @throws std::invalid_argument if denominator is non-positive
 */
inline double dopplerFrequency(double sourceFrequency, double soundVelocity,
                               double observerVelocity, double sourceVelocity) {
    double denominator = soundVelocity - sourceVelocity;
    if (denominator <= 0.0) {
        throw std::invalid_argument("Source velocity must be less than sound velocity");
    }

    return sourceFrequency * (soundVelocity + observerVelocity) / denominator;
}

/**
 * @brief Calculate Doppler shift for source moving toward stationary observer
 *
 * f' = f × v/(v - v_s)
 *
 * @param sourceFrequency Emitted frequency (Hz)
 * @param soundVelocity Velocity of sound (m/s)
 * @param sourceVelocity Velocity of source toward observer (m/s)
 * @return Observed frequency (Hz)
 */
inline double dopplerSourceApproaching(double sourceFrequency, double soundVelocity,
                                       double sourceVelocity) {
    return dopplerFrequency(sourceFrequency, soundVelocity, 0.0, sourceVelocity);
}

/**
 * @brief Calculate Doppler shift for source moving away from stationary observer
 *
 * f' = f × v/(v + v_s)
 *
 * @param sourceFrequency Emitted frequency (Hz)
 * @param soundVelocity Velocity of sound (m/s)
 * @param sourceVelocity Velocity of source away from observer (m/s)
 * @return Observed frequency (Hz)
 */
inline double dopplerSourceReceding(double sourceFrequency, double soundVelocity,
                                    double sourceVelocity) {
    return dopplerFrequency(sourceFrequency, soundVelocity, 0.0, -sourceVelocity);
}

/**
 * @brief Calculate Doppler shift for observer moving toward stationary source
 *
 * f' = f × (v + v_o)/v
 *
 * @param sourceFrequency Emitted frequency (Hz)
 * @param soundVelocity Velocity of sound (m/s)
 * @param observerVelocity Velocity of observer toward source (m/s)
 * @return Observed frequency (Hz)
 */
inline double dopplerObserverApproaching(double sourceFrequency, double soundVelocity,
                                         double observerVelocity) {
    return dopplerFrequency(sourceFrequency, soundVelocity, observerVelocity, 0.0);
}

/**
 * @brief Calculate beat frequency from two interfering waves
 *
 * f_beat = |f₁ - f₂|
 *
 * @param frequency1 First frequency (Hz)
 * @param frequency2 Second frequency (Hz)
 * @return Beat frequency (Hz)
 */
inline double calculateBeatFrequency(double frequency1, double frequency2) {
    return std::abs(frequency1 - frequency2);
}

// ============================================================================
// TRANSVERSE WAVES IN STRINGS (CORDS)
// ============================================================================

/**
 * @brief Calculate velocity of transverse wave in a string
 *
 * v = √(T/μ)
 *
 * where T is tension and μ is linear mass density (mass per unit length)
 *
 * @param tension Tension in string (N)
 * @param linearDensity Mass per unit length (kg/m)
 * @return Wave velocity (m/s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double stringWaveVelocity(double tension, double linearDensity) {
    if (tension <= 0.0 || linearDensity <= 0.0) {
        throw std::invalid_argument("Tension and linear density must be positive");
    }
    return std::sqrt(tension / linearDensity);
}

/**
 * @brief Calculate linear mass density from total mass and length
 *
 * μ = m/L
 *
 * @param mass Total mass of string (kg)
 * @param length Length of string (m)
 * @return Linear mass density (kg/m)
 * @throws std::invalid_argument if length is non-positive
 */
inline double calculateLinearDensity(double mass, double length) {
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    return mass / length;
}

/**
 * @brief Calculate fundamental frequency of vibrating string
 *
 * f₁ = v/(2L) = (1/2L)√(T/μ)
 *
 * First harmonic (n = 1)
 *
 * @param length Length of string (m)
 * @param tension Tension in string (N)
 * @param linearDensity Mass per unit length (kg/m)
 * @return Fundamental frequency (Hz)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double stringFundamentalFrequency(double length, double tension,
                                         double linearDensity) {
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }

    double velocity = stringWaveVelocity(tension, linearDensity);
    return velocity / (2.0 * length);
}

/**
 * @brief Calculate nth harmonic frequency of vibrating string
 *
 * f_n = n × f₁ = (n/2L)√(T/μ)
 *
 * @param harmonicNumber Harmonic number (n = 1, 2, 3, ...)
 * @param length Length of string (m)
 * @param tension Tension in string (N)
 * @param linearDensity Mass per unit length (kg/m)
 * @return Frequency of nth harmonic (Hz)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double stringHarmonicFrequency(int harmonicNumber, double length,
                                      double tension, double linearDensity) {
    if (harmonicNumber < 1) {
        throw std::invalid_argument("Harmonic number must be at least 1");
    }

    double fundamentalFreq = stringFundamentalFrequency(length, tension, linearDensity);
    return harmonicNumber * fundamentalFreq;
}

/**
 * @brief Calculate wavelength of nth harmonic in a string
 *
 * λ_n = 2L/n
 *
 * @param length Length of string (m)
 * @param harmonicNumber Harmonic number (n = 1, 2, 3, ...)
 * @return Wavelength (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double stringHarmonicWavelength(double length, int harmonicNumber) {
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (harmonicNumber < 1) {
        throw std::invalid_argument("Harmonic number must be at least 1");
    }

    return (2.0 * length) / harmonicNumber;
}

/**
 * @brief Calculate required tension for desired fundamental frequency
 *
 * T = 4μL²f²
 *
 * @param targetFrequency Desired fundamental frequency (Hz)
 * @param length Length of string (m)
 * @param linearDensity Mass per unit length (kg/m)
 * @return Required tension (N)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double calculateRequiredTension(double targetFrequency, double length,
                                       double linearDensity) {
    if (targetFrequency <= 0.0 || length <= 0.0 || linearDensity <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }

    return 4.0 * linearDensity * length * length * targetFrequency * targetFrequency;
}

// ============================================================================
// STANDING WAVES AND RESONANCE
// ============================================================================

/**
 * @brief Calculate resonant frequencies in a tube open at both ends
 *
 * f_n = n × v/(2L)  (n = 1, 2, 3, ...)
 *
 * All harmonics are present
 *
 * @param harmonicNumber Harmonic number (n = 1, 2, 3, ...)
 * @param length Length of tube (m)
 * @param soundVelocity Velocity of sound (m/s)
 * @return Resonant frequency (Hz)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double openTubeResonance(int harmonicNumber, double length,
                                double soundVelocity) {
    if (harmonicNumber < 1) {
        throw std::invalid_argument("Harmonic number must be at least 1");
    }
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }

    return harmonicNumber * soundVelocity / (2.0 * length);
}

/**
 * @brief Calculate resonant frequencies in a tube closed at one end
 *
 * f_n = n × v/(4L)  (n = 1, 3, 5, ... odd only)
 *
 * Only odd harmonics are present
 *
 * @param oddHarmonicNumber Odd harmonic number (n = 1, 3, 5, ...)
 * @param length Length of tube (m)
 * @param soundVelocity Velocity of sound (m/s)
 * @return Resonant frequency (Hz)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double closedTubeResonance(int oddHarmonicNumber, double length,
                                  double soundVelocity) {
    if (oddHarmonicNumber < 1 || oddHarmonicNumber % 2 == 0) {
        throw std::invalid_argument("Harmonic number must be odd (1, 3, 5, ...)");
    }
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }

    return oddHarmonicNumber * soundVelocity / (4.0 * length);
}

// ============================================================================
// WAVE ENERGY AND POWER
// ============================================================================

/**
 * @brief Calculate energy density of a wave
 *
 * u = ½ρω²A²
 *
 * Energy per unit volume for a sinusoidal wave
 *
 * @param density Medium density (kg/m³)
 * @param angularFrequency Angular frequency (rad/s)
 * @param amplitude Wave amplitude (m)
 * @return Energy density (J/m³)
 * @throws std::invalid_argument if density is non-positive
 */
inline double waveEnergyDensity(double density, double angularFrequency,
                                double amplitude) {
    if (density <= 0.0) {
        throw std::invalid_argument("Density must be positive");
    }
    return 0.5 * density * angularFrequency * angularFrequency *
           amplitude * amplitude;
}

/**
 * @brief Calculate power transmitted by a wave
 *
 * P = ½ρvω²A² × Area
 *
 * @param density Medium density (kg/m³)
 * @param waveVelocity Wave velocity (m/s)
 * @param angularFrequency Angular frequency (rad/s)
 * @param amplitude Wave amplitude (m)
 * @param area Cross-sectional area (m²)
 * @return Power (W)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double wavePower(double density, double waveVelocity, double angularFrequency,
                        double amplitude, double area) {
    if (density <= 0.0 || area <= 0.0) {
        throw std::invalid_argument("Density and area must be positive");
    }

    double energyDensity = waveEnergyDensity(density, angularFrequency, amplitude);
    return energyDensity * waveVelocity * area;
}

} // namespace wave_mechanics
} // namespace physics

#endif // PHYSICS_WAVE_MECHANICS_HPP
