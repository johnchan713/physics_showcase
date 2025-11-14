#ifndef PHYSICS_THERMODYNAMICS_HPP
#define PHYSICS_THERMODYNAMICS_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace thermodynamics {

/**
 * @brief Thermodynamics and Gas Laws
 *
 * Functions for analyzing gas behavior, pressure, temperature, and volume relationships.
 *
 * Key concepts:
 * - Ideal Gas Law: PV = nRT
 * - Boyle's Law: P₁V₁ = P₂V₂ (constant T)
 * - Charles's Law: V₁/T₁ = V₂/T₂ (constant P)
 * - Gay-Lussac's Law: P₁/T₁ = P₂/T₂ (constant V)
 * - Combined Gas Law: (P₁V₁)/T₁ = (P₂V₂)/T₂
 */

// ============================================================================
// Physical Constants
// ============================================================================

namespace constants {
    constexpr double R = 8.314;              // Universal gas constant (J/(mol·K))
    constexpr double k_B = 1.381e-23;        // Boltzmann constant (J/K)
    constexpr double N_A = 6.022e23;         // Avogadro's number (molecules/mol)
    constexpr double STP_PRESSURE = 101325;  // Standard pressure (Pa)
    constexpr double STP_TEMPERATURE = 273.15; // Standard temperature (K)
}

// ============================================================================
// Boyle's Law and Variations
// ============================================================================

/**
 * @brief Boyle's Law: Calculate final pressure at constant temperature
 *
 * P₁V₁ = P₂V₂ (at constant temperature)
 * Therefore: P₂ = P₁V₁/V₂
 *
 * Robert Boyle (1662): Pressure and volume are inversely proportional
 *
 * @param initialPressure Initial pressure (in Pa, must be > 0)
 * @param initialVolume Initial volume (in m³, must be > 0)
 * @param finalVolume Final volume (in m³, must be > 0)
 * @return Final pressure (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double boylesLaw(double initialPressure, double initialVolume, double finalVolume) {
    if (initialPressure <= 0) {
        throw std::invalid_argument("Initial pressure must be positive");
    }
    if (initialVolume <= 0 || finalVolume <= 0) {
        throw std::invalid_argument("Volumes must be positive");
    }
    return (initialPressure * initialVolume) / finalVolume;
}

/**
 * @brief Charles's Law: Calculate final volume at constant pressure
 *
 * V₁/T₁ = V₂/T₂ (at constant pressure)
 * Therefore: V₂ = V₁(T₂/T₁)
 *
 * Jacques Charles (1787): Volume and temperature are directly proportional
 *
 * @param initialVolume Initial volume (in m³, must be > 0)
 * @param initialTemp Initial temperature (in Kelvin, must be > 0)
 * @param finalTemp Final temperature (in Kelvin, must be > 0)
 * @return Final volume (in m³)
 * @throws std::invalid_argument if parameters out of range
 */
inline double charlesLaw(double initialVolume, double initialTemp, double finalTemp) {
    if (initialVolume <= 0) {
        throw std::invalid_argument("Initial volume must be positive");
    }
    if (initialTemp <= 0 || finalTemp <= 0) {
        throw std::invalid_argument("Temperatures must be positive (Kelvin)");
    }
    return initialVolume * (finalTemp / initialTemp);
}

/**
 * @brief Gay-Lussac's Law: Calculate final pressure at constant volume
 *
 * P₁/T₁ = P₂/T₂ (at constant volume)
 * Therefore: P₂ = P₁(T₂/T₁)
 *
 * Joseph Gay-Lussac (1802): Pressure and temperature are directly proportional
 *
 * @param initialPressure Initial pressure (in Pa, must be > 0)
 * @param initialTemp Initial temperature (in Kelvin, must be > 0)
 * @param finalTemp Final temperature (in Kelvin, must be > 0)
 * @return Final pressure (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double gayLussacsLaw(double initialPressure, double initialTemp, double finalTemp) {
    if (initialPressure <= 0) {
        throw std::invalid_argument("Initial pressure must be positive");
    }
    if (initialTemp <= 0 || finalTemp <= 0) {
        throw std::invalid_argument("Temperatures must be positive (Kelvin)");
    }
    return initialPressure * (finalTemp / initialTemp);
}

/**
 * @brief Combined Gas Law: General case with changing P, V, T
 *
 * (P₁V₁)/T₁ = (P₂V₂)/T₂
 * Therefore: P₂ = P₁(V₁/V₂)(T₂/T₁)
 *
 * Combines Boyle's, Charles's, and Gay-Lussac's laws
 *
 * @param initialPressure Initial pressure (in Pa, must be > 0)
 * @param initialVolume Initial volume (in m³, must be > 0)
 * @param initialTemp Initial temperature (in Kelvin, must be > 0)
 * @param finalVolume Final volume (in m³, must be > 0)
 * @param finalTemp Final temperature (in Kelvin, must be > 0)
 * @return Final pressure (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double combinedGasLaw(double initialPressure, double initialVolume, double initialTemp,
                             double finalVolume, double finalTemp) {
    if (initialPressure <= 0) {
        throw std::invalid_argument("Initial pressure must be positive");
    }
    if (initialVolume <= 0 || finalVolume <= 0) {
        throw std::invalid_argument("Volumes must be positive");
    }
    if (initialTemp <= 0 || finalTemp <= 0) {
        throw std::invalid_argument("Temperatures must be positive (Kelvin)");
    }
    return initialPressure * (initialVolume / finalVolume) * (finalTemp / initialTemp);
}

// ============================================================================
// Ideal Gas Law
// ============================================================================

/**
 * @brief Ideal Gas Law: PV = nRT
 *
 * Calculate pressure from amount, volume, and temperature
 *
 * @param moles Amount of gas (in moles, must be > 0)
 * @param volume Volume (in m³, must be > 0)
 * @param temperature Temperature (in Kelvin, must be > 0)
 * @param R Gas constant (default: 8.314 J/(mol·K))
 * @return Pressure (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double idealGasLawPressure(double moles, double volume, double temperature,
                                  double R = constants::R) {
    if (moles <= 0) {
        throw std::invalid_argument("Moles must be positive");
    }
    if (volume <= 0) {
        throw std::invalid_argument("Volume must be positive");
    }
    if (temperature <= 0) {
        throw std::invalid_argument("Temperature must be positive (Kelvin)");
    }
    return (moles * R * temperature) / volume;
}

/**
 * @brief Calculate volume from ideal gas law
 *
 * V = nRT/P
 *
 * @param moles Amount of gas (in moles, must be > 0)
 * @param temperature Temperature (in Kelvin, must be > 0)
 * @param pressure Pressure (in Pa, must be > 0)
 * @param R Gas constant (default: 8.314 J/(mol·K))
 * @return Volume (in m³)
 * @throws std::invalid_argument if parameters out of range
 */
inline double idealGasLawVolume(double moles, double temperature, double pressure,
                                double R = constants::R) {
    if (moles <= 0) {
        throw std::invalid_argument("Moles must be positive");
    }
    if (temperature <= 0) {
        throw std::invalid_argument("Temperature must be positive (Kelvin)");
    }
    if (pressure <= 0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    return (moles * R * temperature) / pressure;
}

/**
 * @brief Calculate temperature from ideal gas law
 *
 * T = PV/(nR)
 *
 * @param pressure Pressure (in Pa, must be > 0)
 * @param volume Volume (in m³, must be > 0)
 * @param moles Amount of gas (in moles, must be > 0)
 * @param R Gas constant (default: 8.314 J/(mol·K))
 * @return Temperature (in Kelvin)
 * @throws std::invalid_argument if parameters out of range
 */
inline double idealGasLawTemperature(double pressure, double volume, double moles,
                                     double R = constants::R) {
    if (pressure <= 0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    if (volume <= 0) {
        throw std::invalid_argument("Volume must be positive");
    }
    if (moles <= 0) {
        throw std::invalid_argument("Moles must be positive");
    }
    return (pressure * volume) / (moles * R);
}

/**
 * @brief Calculate number of moles from ideal gas law
 *
 * n = PV/(RT)
 *
 * @param pressure Pressure (in Pa, must be > 0)
 * @param volume Volume (in m³, must be > 0)
 * @param temperature Temperature (in Kelvin, must be > 0)
 * @param R Gas constant (default: 8.314 J/(mol·K))
 * @return Number of moles
 * @throws std::invalid_argument if parameters out of range
 */
inline double idealGasLawMoles(double pressure, double volume, double temperature,
                               double R = constants::R) {
    if (pressure <= 0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    if (volume <= 0) {
        throw std::invalid_argument("Volume must be positive");
    }
    if (temperature <= 0) {
        throw std::invalid_argument("Temperature must be positive (Kelvin)");
    }
    return (pressure * volume) / (R * temperature);
}

// ============================================================================
// Pressure of a Gas
// ============================================================================

/**
 * @brief Calculate pressure from molecular kinetic theory
 *
 * P = (1/3)(N/V)m<v²>
 * where N is number of molecules, m is molecular mass, <v²> is mean square velocity
 *
 * @param numberDensity Number density N/V (molecules per m³, must be > 0)
 * @param molecularMass Mass of one molecule (in kg, must be > 0)
 * @param meanSquareVelocity Mean square velocity (in m²/s², must be >= 0)
 * @return Pressure (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double pressureFromKineticTheory(double numberDensity, double molecularMass,
                                       double meanSquareVelocity) {
    if (numberDensity <= 0) {
        throw std::invalid_argument("Number density must be positive");
    }
    if (molecularMass <= 0) {
        throw std::invalid_argument("Molecular mass must be positive");
    }
    if (meanSquareVelocity < 0) {
        throw std::invalid_argument("Mean square velocity must be non-negative");
    }
    return (1.0 / 3.0) * numberDensity * molecularMass * meanSquareVelocity;
}

/**
 * @brief Calculate RMS (root mean square) velocity of gas molecules
 *
 * v_rms = √<v²> = √(3kT/m) = √(3RT/M)
 * where k is Boltzmann constant, T is temperature, m is molecular mass
 *
 * @param temperature Temperature (in Kelvin, must be > 0)
 * @param molecularMass Mass of one molecule (in kg, must be > 0)
 * @param k_B Boltzmann constant (default: 1.381e-23 J/K)
 * @return RMS velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double rmsVelocity(double temperature, double molecularMass, double k_B = constants::k_B) {
    if (temperature <= 0) {
        throw std::invalid_argument("Temperature must be positive (Kelvin)");
    }
    if (molecularMass <= 0) {
        throw std::invalid_argument("Molecular mass must be positive");
    }
    return std::sqrt(3.0 * k_B * temperature / molecularMass);
}

/**
 * @brief Calculate RMS velocity using molar mass
 *
 * v_rms = √(3RT/M)
 * where R is gas constant, M is molar mass
 *
 * @param temperature Temperature (in Kelvin, must be > 0)
 * @param molarMass Molar mass (in kg/mol, must be > 0)
 * @param R Gas constant (default: 8.314 J/(mol·K))
 * @return RMS velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double rmsVelocityMolar(double temperature, double molarMass, double R = constants::R) {
    if (temperature <= 0) {
        throw std::invalid_argument("Temperature must be positive (Kelvin)");
    }
    if (molarMass <= 0) {
        throw std::invalid_argument("Molar mass must be positive");
    }
    return std::sqrt(3.0 * R * temperature / molarMass);
}

// ============================================================================
// Elasticity of Gases
// ============================================================================

/**
 * @brief Calculate isothermal bulk modulus of gas
 *
 * For isothermal process (constant T):
 * K_T = -V(dP/dV) = P
 *
 * Bulk modulus equals the pressure for ideal gas at constant temperature
 *
 * @param pressure Pressure (in Pa, must be > 0)
 * @return Isothermal bulk modulus (in Pa)
 * @throws std::invalid_argument if pressure <= 0
 */
inline double isothermalBulkModulus(double pressure) {
    if (pressure <= 0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    return pressure;
}

/**
 * @brief Calculate adiabatic bulk modulus of gas
 *
 * For adiabatic process (no heat transfer):
 * K_S = -V(dP/dV) = γP
 * where γ is the heat capacity ratio (C_p/C_v)
 *
 * @param pressure Pressure (in Pa, must be > 0)
 * @param gamma Heat capacity ratio (dimensionless, must be > 1)
 * @return Adiabatic bulk modulus (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double adiabaticBulkModulus(double pressure, double gamma) {
    if (pressure <= 0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    if (gamma <= 1.0) {
        throw std::invalid_argument("Gamma must be greater than 1");
    }
    return gamma * pressure;
}

/**
 * @brief Calculate compressibility of gas
 *
 * Compressibility κ = 1/K = -1/V (dV/dP)
 * For isothermal: κ = 1/P
 *
 * @param pressure Pressure (in Pa, must be > 0)
 * @return Isothermal compressibility (in Pa⁻¹)
 * @throws std::invalid_argument if pressure <= 0
 */
inline double isothermalCompressibility(double pressure) {
    if (pressure <= 0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    return 1.0 / pressure;
}

// ============================================================================
// Temperature Conversions
// ============================================================================

/**
 * @brief Convert Celsius to Kelvin
 *
 * K = °C + 273.15
 *
 * @param celsius Temperature in Celsius
 * @return Temperature in Kelvin
 */
inline double celsiusToKelvin(double celsius) {
    return celsius + 273.15;
}

/**
 * @brief Convert Kelvin to Celsius
 *
 * °C = K - 273.15
 *
 * @param kelvin Temperature in Kelvin (must be >= 0)
 * @return Temperature in Celsius
 * @throws std::invalid_argument if kelvin < 0
 */
inline double kelvinToCelsius(double kelvin) {
    if (kelvin < 0) {
        throw std::invalid_argument("Temperature in Kelvin cannot be negative");
    }
    return kelvin - 273.15;
}

/**
 * @brief Convert Fahrenheit to Kelvin
 *
 * K = (°F + 459.67) × 5/9
 *
 * @param fahrenheit Temperature in Fahrenheit
 * @return Temperature in Kelvin
 */
inline double fahrenheitToKelvin(double fahrenheit) {
    return (fahrenheit + 459.67) * (5.0 / 9.0);
}

/**
 * @brief Convert Kelvin to Fahrenheit
 *
 * °F = K × 9/5 - 459.67
 *
 * @param kelvin Temperature in Kelvin (must be >= 0)
 * @return Temperature in Fahrenheit
 * @throws std::invalid_argument if kelvin < 0
 */
inline double kelvinToFahrenheit(double kelvin) {
    if (kelvin < 0) {
        throw std::invalid_argument("Temperature in Kelvin cannot be negative");
    }
    return kelvin * (9.0 / 5.0) - 459.67;
}

// ============================================================================
// Standard Conditions
// ============================================================================

/**
 * @brief Calculate molar volume at STP (Standard Temperature and Pressure)
 *
 * At STP (0°C, 1 atm): V_m = RT/P = 22.4 L/mol
 *
 * @param R Gas constant (default: 8.314 J/(mol·K))
 * @return Molar volume at STP (in m³/mol)
 */
inline double molarVolumeAtSTP(double R = constants::R) {
    return (R * constants::STP_TEMPERATURE) / constants::STP_PRESSURE;
}

/**
 * @brief Calculate density of gas at given conditions
 *
 * ρ = PM/(RT)
 * where M is molar mass
 *
 * @param pressure Pressure (in Pa, must be > 0)
 * @param molarMass Molar mass (in kg/mol, must be > 0)
 * @param temperature Temperature (in Kelvin, must be > 0)
 * @param R Gas constant (default: 8.314 J/(mol·K))
 * @return Density (in kg/m³)
 * @throws std::invalid_argument if parameters out of range
 */
inline double gasDensity(double pressure, double molarMass, double temperature,
                        double R = constants::R) {
    if (pressure <= 0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    if (molarMass <= 0) {
        throw std::invalid_argument("Molar mass must be positive");
    }
    if (temperature <= 0) {
        throw std::invalid_argument("Temperature must be positive (Kelvin)");
    }
    return (pressure * molarMass) / (R * temperature);
}

} // namespace thermodynamics
} // namespace physics

#endif // PHYSICS_THERMODYNAMICS_HPP
