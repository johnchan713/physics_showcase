#ifndef PHYSICS_ELECTROMAGNETIC_WAVES_HPP
#define PHYSICS_ELECTROMAGNETIC_WAVES_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file electromagnetic_waves.hpp
 * @brief Comprehensive implementation of electromagnetic wave theory
 *
 * This module implements:
 * - Maxwell's equations and wave propagation
 * - Speed of electromagnetic waves
 * - Relationship between E and B fields
 * - Wave energy and intensity (Poynting vector)
 * - Energy density in EM waves
 * - Radiation pressure
 * - Wavelength and frequency relationships
 * - Wave propagation in different media
 * - Electromagnetic spectrum
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace electromagnetic_waves {

/**
 * @namespace constants
 * @brief Physical constants for electromagnetic waves
 */
namespace constants {
    constexpr double EPSILON_0 = 8.854187817e-12;   // Permittivity of free space (F/m)
    constexpr double MU_0 = 4.0 * M_PI * 1e-7;      // Permeability of free space (H/m)
    constexpr double SPEED_OF_LIGHT = 299792458.0;  // Speed of light in vacuum (m/s)
    constexpr double PLANCK_H = 6.62607015e-34;     // Planck's constant (J⋅s)
}

// ============================================================================
// SPEED OF ELECTROMAGNETIC WAVES
// ============================================================================

/**
 * @brief Calculate speed of EM waves in vacuum
 *
 * c = 1 / √(μ₀ε₀)
 *
 * @param epsilon0 Permittivity of free space (F/m)
 * @param mu0 Permeability of free space (H/m)
 * @return Speed of light (m/s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double speedOfLight(double epsilon0 = constants::EPSILON_0,
                          double mu0 = constants::MU_0) {
    if (epsilon0 <= 0.0 || mu0 <= 0.0) {
        throw std::invalid_argument("Permittivity and permeability must be positive");
    }
    return 1.0 / std::sqrt(epsilon0 * mu0);
}

/**
 * @brief Calculate speed of EM waves in medium
 *
 * v = c / n = 1 / √(με)
 *
 * @param permittivity Permittivity of medium (F/m)
 * @param permeability Permeability of medium (H/m)
 * @return Speed in medium (m/s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double speedInMedium(double permittivity, double permeability) {
    if (permittivity <= 0.0 || permeability <= 0.0) {
        throw std::invalid_argument("Permittivity and permeability must be positive");
    }
    return 1.0 / std::sqrt(permittivity * permeability);
}

/**
 * @brief Calculate refractive index from permittivity and permeability
 *
 * n = √(εᵣ × μᵣ) = √[(ε/ε₀) × (μ/μ₀)]
 *
 * @param relativePermittivity Relative permittivity (εᵣ)
 * @param relativePermeability Relative permeability (μᵣ)
 * @return Refractive index
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double refractiveIndex(double relativePermittivity,
                             double relativePermeability = 1.0) {
    if (relativePermittivity <= 0.0 || relativePermeability <= 0.0) {
        throw std::invalid_argument("Relative permittivity and permeability must be positive");
    }
    return std::sqrt(relativePermittivity * relativePermeability);
}

/**
 * @brief Calculate speed from refractive index
 *
 * v = c / n
 *
 * @param refractiveIndex Refractive index of medium
 * @param speedInVacuum Speed of light in vacuum (m/s)
 * @return Speed in medium (m/s)
 * @throws std::invalid_argument if refractive index is non-positive
 */
inline double speedFromRefractiveIndex(double refractiveIndex,
                                       double speedInVacuum = constants::SPEED_OF_LIGHT) {
    if (refractiveIndex <= 0.0) {
        throw std::invalid_argument("Refractive index must be positive");
    }
    return speedInVacuum / refractiveIndex;
}

// ============================================================================
// WAVELENGTH AND FREQUENCY RELATIONSHIPS
// ============================================================================

/**
 * @brief Calculate wavelength from frequency
 *
 * λ = c / f
 *
 * @param frequency Frequency (Hz)
 * @param speed Speed of wave (m/s)
 * @return Wavelength (m)
 * @throws std::invalid_argument if frequency is non-positive
 */
inline double wavelengthFromFrequency(double frequency,
                                      double speed = constants::SPEED_OF_LIGHT) {
    if (frequency <= 0.0) {
        throw std::invalid_argument("Frequency must be positive");
    }
    return speed / frequency;
}

/**
 * @brief Calculate frequency from wavelength
 *
 * f = c / λ
 *
 * @param wavelength Wavelength (m)
 * @param speed Speed of wave (m/s)
 * @return Frequency (Hz)
 * @throws std::invalid_argument if wavelength is non-positive
 */
inline double frequencyFromWavelength(double wavelength,
                                      double speed = constants::SPEED_OF_LIGHT) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return speed / wavelength;
}

/**
 * @brief Calculate angular frequency
 *
 * ω = 2πf
 *
 * @param frequency Frequency (Hz)
 * @return Angular frequency (rad/s)
 * @throws std::invalid_argument if frequency is negative
 */
inline double angularFrequency(double frequency) {
    if (frequency < 0.0) {
        throw std::invalid_argument("Frequency cannot be negative");
    }
    return 2.0 * M_PI * frequency;
}

/**
 * @brief Calculate wave number
 *
 * k = 2π / λ = ω / c
 *
 * @param wavelength Wavelength (m)
 * @return Wave number (rad/m)
 * @throws std::invalid_argument if wavelength is non-positive
 */
inline double waveNumber(double wavelength) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return 2.0 * M_PI / wavelength;
}

// ============================================================================
// RELATIONSHIP BETWEEN E AND B FIELDS
// ============================================================================

/**
 * @brief Calculate magnetic field amplitude from electric field
 *
 * B = E / c
 *
 * In EM wave, E and B are perpendicular and in phase
 *
 * @param electricField Electric field amplitude (V/m)
 * @param speed Speed of wave (m/s)
 * @return Magnetic field amplitude (T)
 * @throws std::invalid_argument if speed is non-positive
 */
inline double magneticFieldFromElectric(double electricField,
                                        double speed = constants::SPEED_OF_LIGHT) {
    if (speed <= 0.0) {
        throw std::invalid_argument("Speed must be positive");
    }
    return electricField / speed;
}

/**
 * @brief Calculate electric field amplitude from magnetic field
 *
 * E = c × B
 *
 * @param magneticField Magnetic field amplitude (T)
 * @param speed Speed of wave (m/s)
 * @return Electric field amplitude (V/m)
 */
inline double electricFieldFromMagnetic(double magneticField,
                                        double speed = constants::SPEED_OF_LIGHT) {
    return magneticField * speed;
}

/**
 * @brief Verify E/B ratio equals wave speed
 *
 * c = E / B
 *
 * @param electricField Electric field amplitude (V/m)
 * @param magneticField Magnetic field amplitude (T)
 * @return Wave speed (m/s)
 * @throws std::invalid_argument if magnetic field is zero
 */
inline double waveSpeedFromFields(double electricField, double magneticField) {
    if (std::abs(magneticField) < 1e-15) {
        throw std::invalid_argument("Magnetic field must be non-zero");
    }
    return electricField / magneticField;
}

// ============================================================================
// ENERGY AND INTENSITY
// ============================================================================

/**
 * @brief Calculate energy density in electric field
 *
 * u_E = (1/2) × ε₀ × E²
 *
 * @param electricField Electric field amplitude (V/m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Energy density (J/m³)
 */
inline double electricEnergyDensity(double electricField,
                                    double epsilon0 = constants::EPSILON_0) {
    return 0.5 * epsilon0 * electricField * electricField;
}

/**
 * @brief Calculate energy density in magnetic field
 *
 * u_B = (1/2) × B² / μ₀
 *
 * @param magneticField Magnetic field amplitude (T)
 * @param mu0 Permeability of free space (H/m)
 * @return Energy density (J/m³)
 */
inline double magneticEnergyDensity(double magneticField,
                                    double mu0 = constants::MU_0) {
    return 0.5 * magneticField * magneticField / mu0;
}

/**
 * @brief Calculate total energy density in EM wave
 *
 * u = u_E + u_B = ε₀E²
 *
 * (Since u_E = u_B in EM wave)
 *
 * @param electricField Electric field amplitude (V/m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Total energy density (J/m³)
 */
inline double totalEnergyDensity(double electricField,
                                 double epsilon0 = constants::EPSILON_0) {
    return epsilon0 * electricField * electricField;
}

/**
 * @brief Calculate instantaneous Poynting vector magnitude
 *
 * S = (1/μ₀) × E × B
 *
 * Represents energy flux (power per unit area)
 *
 * @param electricField Electric field (V/m)
 * @param magneticField Magnetic field (T)
 * @param mu0 Permeability of free space (H/m)
 * @return Poynting vector magnitude (W/m²)
 */
inline double poyntingVector(double electricField, double magneticField,
                            double mu0 = constants::MU_0) {
    return (electricField * magneticField) / mu0;
}

/**
 * @brief Calculate intensity (average power per unit area)
 *
 * I = (1/2) × c × ε₀ × E₀²
 *
 * For sinusoidal wave with amplitude E₀
 *
 * @param electricFieldAmplitude Peak electric field (V/m)
 * @param speed Speed of light (m/s)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Intensity (W/m²)
 */
inline double waveIntensity(double electricFieldAmplitude,
                           double speed = constants::SPEED_OF_LIGHT,
                           double epsilon0 = constants::EPSILON_0) {
    return 0.5 * speed * epsilon0 * electricFieldAmplitude * electricFieldAmplitude;
}

/**
 * @brief Calculate intensity from energy density
 *
 * I = u × c
 *
 * @param energyDensity Energy density (J/m³)
 * @param speed Speed of wave (m/s)
 * @return Intensity (W/m²)
 */
inline double intensityFromEnergyDensity(double energyDensity,
                                         double speed = constants::SPEED_OF_LIGHT) {
    return energyDensity * speed;
}

/**
 * @brief Calculate power from intensity and area
 *
 * P = I × A
 *
 * @param intensity Intensity (W/m²)
 * @param area Area (m²)
 * @return Power (W)
 * @throws std::invalid_argument if area is negative
 */
inline double powerFromIntensity(double intensity, double area) {
    if (area < 0.0) {
        throw std::invalid_argument("Area cannot be negative");
    }
    return intensity * area;
}

/**
 * @brief Calculate intensity at distance from point source
 *
 * I = P / (4πr²)
 *
 * Inverse square law for isotropic source
 *
 * @param power Total radiated power (W)
 * @param distance Distance from source (m)
 * @return Intensity (W/m²)
 * @throws std::invalid_argument if distance is non-positive
 */
inline double intensityAtDistance(double power, double distance) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return power / (4.0 * M_PI * distance * distance);
}

// ============================================================================
// RADIATION PRESSURE
// ============================================================================

/**
 * @brief Calculate radiation pressure (complete absorption)
 *
 * P = I / c
 *
 * Pressure exerted by EM wave on perfectly absorbing surface
 *
 * @param intensity Intensity of wave (W/m²)
 * @param speed Speed of light (m/s)
 * @return Radiation pressure (Pa)
 * @throws std::invalid_argument if speed is non-positive
 */
inline double radiationPressureAbsorption(double intensity,
                                          double speed = constants::SPEED_OF_LIGHT) {
    if (speed <= 0.0) {
        throw std::invalid_argument("Speed must be positive");
    }
    return intensity / speed;
}

/**
 * @brief Calculate radiation pressure (complete reflection)
 *
 * P = 2I / c
 *
 * Pressure exerted by EM wave on perfectly reflecting surface
 *
 * @param intensity Intensity of wave (W/m²)
 * @param speed Speed of light (m/s)
 * @return Radiation pressure (Pa)
 * @throws std::invalid_argument if speed is non-positive
 */
inline double radiationPressureReflection(double intensity,
                                          double speed = constants::SPEED_OF_LIGHT) {
    if (speed <= 0.0) {
        throw std::invalid_argument("Speed must be positive");
    }
    return 2.0 * intensity / speed;
}

/**
 * @brief Calculate radiation force on surface
 *
 * F = P × A
 *
 * @param radiationPressure Radiation pressure (Pa)
 * @param area Surface area (m²)
 * @return Force (N)
 * @throws std::invalid_argument if area is negative
 */
inline double radiationForce(double radiationPressure, double area) {
    if (area < 0.0) {
        throw std::invalid_argument("Area cannot be negative");
    }
    return radiationPressure * area;
}

// ============================================================================
// PHOTON ENERGY AND MOMENTUM
// ============================================================================

/**
 * @brief Calculate photon energy from frequency
 *
 * E = h × f
 *
 * @param frequency Frequency (Hz)
 * @param planckConstant Planck's constant (J⋅s)
 * @return Photon energy (J)
 * @throws std::invalid_argument if frequency is negative
 */
inline double photonEnergy(double frequency,
                          double planckConstant = constants::PLANCK_H) {
    if (frequency < 0.0) {
        throw std::invalid_argument("Frequency cannot be negative");
    }
    return planckConstant * frequency;
}

/**
 * @brief Calculate photon energy from wavelength
 *
 * E = hc / λ
 *
 * @param wavelength Wavelength (m)
 * @param planckConstant Planck's constant (J⋅s)
 * @param speedOfLight Speed of light (m/s)
 * @return Photon energy (J)
 * @throws std::invalid_argument if wavelength is non-positive
 */
inline double photonEnergyFromWavelength(double wavelength,
                                         double planckConstant = constants::PLANCK_H,
                                         double speedOfLight = constants::SPEED_OF_LIGHT) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return (planckConstant * speedOfLight) / wavelength;
}

/**
 * @brief Calculate photon momentum
 *
 * p = E / c = h / λ
 *
 * @param wavelength Wavelength (m)
 * @param planckConstant Planck's constant (J⋅s)
 * @return Photon momentum (kg⋅m/s)
 * @throws std::invalid_argument if wavelength is non-positive
 */
inline double photonMomentum(double wavelength,
                            double planckConstant = constants::PLANCK_H) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return planckConstant / wavelength;
}

/**
 * @brief Calculate number of photons per second
 *
 * N = P / E_photon = P × λ / (hc)
 *
 * @param power Total power (W)
 * @param photonEnergy Energy per photon (J)
 * @return Photons per second
 * @throws std::invalid_argument if photon energy is non-positive
 */
inline double photonsPerSecond(double power, double photonEnergy) {
    if (photonEnergy <= 0.0) {
        throw std::invalid_argument("Photon energy must be positive");
    }
    return power / photonEnergy;
}

// ============================================================================
// WAVE PROPAGATION IN DIFFERENT MEDIA
// ============================================================================

/**
 * @brief Calculate wavelength in medium
 *
 * λ_medium = λ_vacuum / n
 *
 * Frequency remains constant, wavelength changes
 *
 * @param vacuumWavelength Wavelength in vacuum (m)
 * @param refractiveIndex Refractive index of medium
 * @return Wavelength in medium (m)
 * @throws std::invalid_argument if refractive index is non-positive
 */
inline double wavelengthInMedium(double vacuumWavelength, double refractiveIndex) {
    if (refractiveIndex <= 0.0) {
        throw std::invalid_argument("Refractive index must be positive");
    }
    return vacuumWavelength / refractiveIndex;
}

/**
 * @brief Calculate phase velocity in medium
 *
 * v_p = c / n
 *
 * @param refractiveIndex Refractive index
 * @param speedInVacuum Speed in vacuum (m/s)
 * @return Phase velocity (m/s)
 * @throws std::invalid_argument if refractive index is non-positive
 */
inline double phaseVelocity(double refractiveIndex,
                           double speedInVacuum = constants::SPEED_OF_LIGHT) {
    if (refractiveIndex <= 0.0) {
        throw std::invalid_argument("Refractive index must be positive");
    }
    return speedInVacuum / refractiveIndex;
}

/**
 * @brief Calculate impedance of medium
 *
 * Z = √(μ/ε)
 *
 * For vacuum: Z₀ = √(μ₀/ε₀) ≈ 377 Ω
 *
 * @param permeability Permeability (H/m)
 * @param permittivity Permittivity (F/m)
 * @return Impedance (Ω)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double mediumImpedance(double permeability, double permittivity) {
    if (permeability <= 0.0 || permittivity <= 0.0) {
        throw std::invalid_argument("Permeability and permittivity must be positive");
    }
    return std::sqrt(permeability / permittivity);
}

/**
 * @brief Calculate impedance of free space
 *
 * Z₀ = √(μ₀/ε₀) ≈ 376.73 Ω
 *
 * @param mu0 Permeability of free space (H/m)
 * @param epsilon0 Permittivity of free space (F/m)
 * @return Impedance of free space (Ω)
 */
inline double freeSpaceImpedance(double mu0 = constants::MU_0,
                                double epsilon0 = constants::EPSILON_0) {
    return std::sqrt(mu0 / epsilon0);
}

// ============================================================================
// ELECTROMAGNETIC SPECTRUM
// ============================================================================

/**
 * @brief Check if wavelength is in visible spectrum
 *
 * Visible light: approximately 380-750 nm
 *
 * @param wavelength Wavelength (m)
 * @return true if in visible spectrum
 */
inline bool isVisible(double wavelength) {
    constexpr double MIN_VISIBLE = 380e-9;  // 380 nm
    constexpr double MAX_VISIBLE = 750e-9;  // 750 nm
    return wavelength >= MIN_VISIBLE && wavelength <= MAX_VISIBLE;
}

/**
 * @brief Determine region of EM spectrum
 *
 * Returns approximate region based on wavelength
 *
 * @param wavelength Wavelength (m)
 * @return String describing the region
 */
inline const char* spectrumRegion(double wavelength) {
    if (wavelength < 1e-11) return "Gamma rays";
    if (wavelength < 1e-8) return "X-rays";
    if (wavelength < 380e-9) return "Ultraviolet";
    if (wavelength < 750e-9) return "Visible light";
    if (wavelength < 1e-3) return "Infrared";
    if (wavelength < 1.0) return "Microwaves";
    return "Radio waves";
}

/**
 * @brief Convert wavelength to photon energy in eV
 *
 * E(eV) = 1240 / λ(nm)
 *
 * Common approximation for photon energy
 *
 * @param wavelengthNm Wavelength in nanometers
 * @return Energy in electron volts (eV)
 * @throws std::invalid_argument if wavelength is non-positive
 */
inline double photonEnergyEV(double wavelengthNm) {
    if (wavelengthNm <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return 1240.0 / wavelengthNm;
}

} // namespace electromagnetic_waves
} // namespace physics

#endif // PHYSICS_ELECTROMAGNETIC_WAVES_HPP
