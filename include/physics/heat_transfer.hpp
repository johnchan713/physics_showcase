#ifndef PHYSICS_HEAT_TRANSFER_HPP
#define PHYSICS_HEAT_TRANSFER_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file heat_transfer.hpp
 * @brief Comprehensive implementation of heat transfer mechanisms and thermodynamic processes
 *
 * This module implements:
 * - Thermal conductivity (Fourier's law)
 * - Convection heat transfer
 * - Thermal radiation (Stefan-Boltzmann law)
 * - Heat engines and efficiency (Carnot cycle)
 * - Blackbody radiation and quantum theory
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace heat_transfer {

/**
 * @namespace constants
 * @brief Physical constants for heat transfer calculations
 */
namespace constants {
    // Thermal conductivity (W/(m·K))
    constexpr double SILVER_CONDUCTIVITY = 429.0;
    constexpr double COPPER_CONDUCTIVITY = 401.0;
    constexpr double ALUMINUM_CONDUCTIVITY = 237.0;
    constexpr double STEEL_CONDUCTIVITY = 50.0;
    constexpr double GLASS_CONDUCTIVITY = 1.0;
    constexpr double CONCRETE_CONDUCTIVITY = 1.4;
    constexpr double WOOD_CONDUCTIVITY = 0.15;
    constexpr double AIR_CONDUCTIVITY = 0.024;
    constexpr double WATER_CONDUCTIVITY = 0.6;

    // Stefan-Boltzmann constant
    constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;  // W/(m²·K⁴)

    // Planck's constant
    constexpr double PLANCK_H = 6.62607015e-34;          // J·s

    // Speed of light
    constexpr double SPEED_OF_LIGHT = 299792458.0;       // m/s

    // Boltzmann constant
    constexpr double BOLTZMANN_K = 1.380649e-23;         // J/K
}

// ============================================================================
// THERMAL CONDUCTIVITY (CONDUCTION)
// ============================================================================

/**
 * @brief Calculate heat transfer rate through conduction (Fourier's law)
 *
 * Q/t = kA(T₁ - T₂)/d
 *
 * @param thermalConductivity Thermal conductivity of material (W/(m·K))
 * @param area Cross-sectional area (m²)
 * @param temperatureDifference Temperature difference across material (K or °C)
 * @param thickness Thickness of material (m)
 * @return Heat transfer rate (W)
 * @throws std::invalid_argument if area, conductivity, or thickness is non-positive
 */
inline double conductionHeatRate(double thermalConductivity, double area,
                                 double temperatureDifference, double thickness) {
    if (thermalConductivity <= 0.0 || area <= 0.0 || thickness <= 0.0) {
        throw std::invalid_argument("Conductivity, area, and thickness must be positive");
    }
    return (thermalConductivity * area * temperatureDifference) / thickness;
}

/**
 * @brief Calculate thermal resistance
 *
 * R = d/(kA)
 *
 * Thermal resistance is analogous to electrical resistance
 *
 * @param thickness Thickness of material (m)
 * @param thermalConductivity Thermal conductivity (W/(m·K))
 * @param area Cross-sectional area (m²)
 * @return Thermal resistance (K/W)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double thermalResistance(double thickness, double thermalConductivity, double area) {
    if (thickness <= 0.0 || thermalConductivity <= 0.0 || area <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return thickness / (thermalConductivity * area);
}

/**
 * @brief Calculate heat transfer rate using thermal resistance
 *
 * Q/t = ΔT/R
 *
 * @param temperatureDifference Temperature difference (K or °C)
 * @param thermalResistance Thermal resistance (K/W)
 * @return Heat transfer rate (W)
 * @throws std::invalid_argument if thermal resistance is non-positive
 */
inline double heatRateFromResistance(double temperatureDifference, double thermalResistance) {
    if (thermalResistance <= 0.0) {
        throw std::invalid_argument("Thermal resistance must be positive");
    }
    return temperatureDifference / thermalResistance;
}

/**
 * @brief Calculate thermal resistance for composite wall (series)
 *
 * R_total = R₁ + R₂ + R₃ + ...
 *
 * For layers in series (e.g., wall with multiple materials)
 *
 * @param resistances Vector of thermal resistances (K/W)
 * @return Total thermal resistance (K/W)
 */
inline double seriesThermalResistance(const std::vector<double>& resistances) {
    double total = 0.0;
    for (double r : resistances) {
        if (r < 0.0) {
            throw std::invalid_argument("Resistances must be non-negative");
        }
        total += r;
    }
    return total;
}

/**
 * @brief Calculate temperature gradient in steady-state conduction
 *
 * dT/dx = -Q/(kA)
 *
 * @param heatRate Heat transfer rate (W)
 * @param thermalConductivity Thermal conductivity (W/(m·K))
 * @param area Cross-sectional area (m²)
 * @return Temperature gradient (K/m)
 * @throws std::invalid_argument if conductivity or area is non-positive
 */
inline double temperatureGradient(double heatRate, double thermalConductivity, double area) {
    if (thermalConductivity <= 0.0 || area <= 0.0) {
        throw std::invalid_argument("Conductivity and area must be positive");
    }
    return -heatRate / (thermalConductivity * area);
}

// ============================================================================
// THERMAL RADIATION
// ============================================================================

/**
 * @brief Calculate radiated power using Stefan-Boltzmann law
 *
 * P = εσAT⁴
 *
 * where ε is emissivity (0 to 1), σ is Stefan-Boltzmann constant
 *
 * @param emissivity Emissivity of surface (0 = perfect reflector, 1 = blackbody)
 * @param area Surface area (m²)
 * @param temperature Absolute temperature (K)
 * @param sigma Stefan-Boltzmann constant (W/(m²·K⁴))
 * @return Radiated power (W)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double stefanBoltzmannRadiation(double emissivity, double area, double temperature,
                                       double sigma = constants::STEFAN_BOLTZMANN) {
    if (emissivity < 0.0 || emissivity > 1.0) {
        throw std::invalid_argument("Emissivity must be between 0 and 1");
    }
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    if (temperature <= 0.0) {
        throw std::invalid_argument("Temperature must be positive (in Kelvin)");
    }

    return emissivity * sigma * area * std::pow(temperature, 4);
}

/**
 * @brief Calculate net radiation between object and surroundings
 *
 * P_net = εσA(T⁴ - T_surr⁴)
 *
 * @param emissivity Emissivity of object
 * @param area Surface area (m²)
 * @param objectTemp Temperature of object (K)
 * @param surroundingTemp Temperature of surroundings (K)
 * @param sigma Stefan-Boltzmann constant (W/(m²·K⁴))
 * @return Net radiated power (W), positive = net emission
 * @throws std::invalid_argument if parameters are invalid
 */
inline double netRadiation(double emissivity, double area, double objectTemp,
                           double surroundingTemp,
                           double sigma = constants::STEFAN_BOLTZMANN) {
    if (emissivity < 0.0 || emissivity > 1.0) {
        throw std::invalid_argument("Emissivity must be between 0 and 1");
    }
    if (area <= 0.0 || objectTemp <= 0.0 || surroundingTemp <= 0.0) {
        throw std::invalid_argument("Area and temperatures must be positive");
    }

    return emissivity * sigma * area *
           (std::pow(objectTemp, 4) - std::pow(surroundingTemp, 4));
}

/**
 * @brief Calculate Wien's displacement law (peak wavelength)
 *
 * λ_max = b/T
 *
 * where b = 2.898×10⁻³ m·K (Wien's displacement constant)
 *
 * @param temperature Absolute temperature (K)
 * @return Peak wavelength (m)
 * @throws std::invalid_argument if temperature is non-positive
 */
inline double wienDisplacementLaw(double temperature) {
    if (temperature <= 0.0) {
        throw std::invalid_argument("Temperature must be positive");
    }
    constexpr double WIEN_CONSTANT = 2.897771955e-3;  // m·K
    return WIEN_CONSTANT / temperature;
}

/**
 * @brief Calculate Planck's blackbody radiation formula
 *
 * B(λ,T) = (2hc²/λ⁵) × 1/(e^(hc/λkT) - 1)
 *
 * Spectral radiance per unit wavelength
 *
 * @param wavelength Wavelength (m)
 * @param temperature Absolute temperature (K)
 * @return Spectral radiance (W/(m²·sr·m))
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double planckRadiation(double wavelength, double temperature) {
    if (wavelength <= 0.0 || temperature <= 0.0) {
        throw std::invalid_argument("Wavelength and temperature must be positive");
    }

    const double h = constants::PLANCK_H;
    const double c = constants::SPEED_OF_LIGHT;
    const double k = constants::BOLTZMANN_K;

    double lambda5 = std::pow(wavelength, 5);
    double numerator = 2.0 * h * c * c;
    double exponent = (h * c) / (wavelength * k * temperature);
    double denominator = lambda5 * (std::exp(exponent) - 1.0);

    return numerator / denominator;
}

/**
 * @brief Calculate photon energy
 *
 * E = hf = hc/λ
 *
 * Quantum theory: energy is quantized in photons
 *
 * @param wavelength Wavelength of photon (m)
 * @return Energy of photon (J)
 * @throws std::invalid_argument if wavelength is non-positive
 */
inline double photonEnergy(double wavelength) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return (constants::PLANCK_H * constants::SPEED_OF_LIGHT) / wavelength;
}

/**
 * @brief Calculate photon energy from frequency
 *
 * E = hf
 *
 * @param frequency Frequency of photon (Hz)
 * @return Energy of photon (J)
 * @throws std::invalid_argument if frequency is non-positive
 */
inline double photonEnergyFromFrequency(double frequency) {
    if (frequency <= 0.0) {
        throw std::invalid_argument("Frequency must be positive");
    }
    return constants::PLANCK_H * frequency;
}

// ============================================================================
// HEAT ENGINES AND EFFICIENCY
// ============================================================================

/**
 * @brief Calculate thermal efficiency of heat engine
 *
 * η = W/Q_h = (Q_h - Q_c)/Q_h = 1 - Q_c/Q_h
 *
 * where W is work output, Q_h is heat input, Q_c is heat rejected
 *
 * @param heatInput Heat absorbed from hot reservoir (J)
 * @param heatOutput Heat rejected to cold reservoir (J)
 * @return Efficiency (dimensionless, 0 to 1)
 * @throws std::invalid_argument if heat input is non-positive
 */
inline double engineEfficiency(double heatInput, double heatOutput) {
    if (heatInput <= 0.0) {
        throw std::invalid_argument("Heat input must be positive");
    }
    if (heatOutput < 0.0 || heatOutput > heatInput) {
        throw std::invalid_argument("Heat output must be between 0 and heat input");
    }
    return 1.0 - (heatOutput / heatInput);
}

/**
 * @brief Calculate Carnot efficiency (maximum possible efficiency)
 *
 * η_Carnot = 1 - T_c/T_h
 *
 * @param hotTemperature Temperature of hot reservoir (K)
 * @param coldTemperature Temperature of cold reservoir (K)
 * @return Carnot efficiency (dimensionless, 0 to 1)
 * @throws std::invalid_argument if temperatures are non-positive or T_c > T_h
 */
inline double carnotEfficiency(double hotTemperature, double coldTemperature) {
    if (hotTemperature <= 0.0 || coldTemperature <= 0.0) {
        throw std::invalid_argument("Temperatures must be positive (in Kelvin)");
    }
    if (coldTemperature >= hotTemperature) {
        throw std::invalid_argument("Cold temperature must be less than hot temperature");
    }
    return 1.0 - (coldTemperature / hotTemperature);
}

/**
 * @brief Calculate work output from heat engine
 *
 * W = Q_h - Q_c = ηQ_h
 *
 * @param heatInput Heat absorbed (J)
 * @param efficiency Engine efficiency (0 to 1)
 * @return Work output (J)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double engineWorkOutput(double heatInput, double efficiency) {
    if (heatInput <= 0.0) {
        throw std::invalid_argument("Heat input must be positive");
    }
    if (efficiency < 0.0 || efficiency > 1.0) {
        throw std::invalid_argument("Efficiency must be between 0 and 1");
    }
    return efficiency * heatInput;
}

/**
 * @brief Calculate heat rejected by engine
 *
 * Q_c = Q_h(1 - η)
 *
 * @param heatInput Heat absorbed (J)
 * @param efficiency Engine efficiency (0 to 1)
 * @return Heat rejected (J)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double engineHeatRejected(double heatInput, double efficiency) {
    if (heatInput <= 0.0) {
        throw std::invalid_argument("Heat input must be positive");
    }
    if (efficiency < 0.0 || efficiency > 1.0) {
        throw std::invalid_argument("Efficiency must be between 0 and 1");
    }
    return heatInput * (1.0 - efficiency);
}

/**
 * @brief Calculate coefficient of performance for refrigerator
 *
 * COP_refrigerator = Q_c/W = Q_c/(Q_h - Q_c)
 *
 * @param heatRemoved Heat removed from cold reservoir (J)
 * @param workInput Work input required (J)
 * @return Coefficient of performance (dimensionless, can be > 1)
 * @throws std::invalid_argument if work is non-positive
 */
inline double refrigeratorCOP(double heatRemoved, double workInput) {
    if (workInput <= 0.0) {
        throw std::invalid_argument("Work input must be positive");
    }
    if (heatRemoved < 0.0) {
        throw std::invalid_argument("Heat removed must be non-negative");
    }
    return heatRemoved / workInput;
}

/**
 * @brief Calculate Carnot COP for refrigerator (maximum possible)
 *
 * COP_Carnot = T_c/(T_h - T_c)
 *
 * @param coldTemperature Temperature of cold reservoir (K)
 * @param hotTemperature Temperature of hot reservoir (K)
 * @return Carnot COP (dimensionless)
 * @throws std::invalid_argument if temperatures are invalid
 */
inline double carnotCOPRefrigerator(double coldTemperature, double hotTemperature) {
    if (coldTemperature <= 0.0 || hotTemperature <= 0.0) {
        throw std::invalid_argument("Temperatures must be positive (in Kelvin)");
    }
    if (coldTemperature >= hotTemperature) {
        throw std::invalid_argument("Cold temperature must be less than hot temperature");
    }
    return coldTemperature / (hotTemperature - coldTemperature);
}

/**
 * @brief Calculate coefficient of performance for heat pump
 *
 * COP_heatpump = Q_h/W = Q_h/(Q_h - Q_c)
 *
 * @param heatDelivered Heat delivered to hot reservoir (J)
 * @param workInput Work input required (J)
 * @return Coefficient of performance (dimensionless, always > 1)
 * @throws std::invalid_argument if work is non-positive
 */
inline double heatPumpCOP(double heatDelivered, double workInput) {
    if (workInput <= 0.0) {
        throw std::invalid_argument("Work input must be positive");
    }
    if (heatDelivered <= 0.0) {
        throw std::invalid_argument("Heat delivered must be positive");
    }
    return heatDelivered / workInput;
}

/**
 * @brief Calculate Carnot COP for heat pump
 *
 * COP_Carnot = T_h/(T_h - T_c)
 *
 * @param hotTemperature Temperature of hot reservoir (K)
 * @param coldTemperature Temperature of cold reservoir (K)
 * @return Carnot COP for heat pump (dimensionless)
 * @throws std::invalid_argument if temperatures are invalid
 */
inline double carnotCOPHeatPump(double hotTemperature, double coldTemperature) {
    if (hotTemperature <= 0.0 || coldTemperature <= 0.0) {
        throw std::invalid_argument("Temperatures must be positive (in Kelvin)");
    }
    if (coldTemperature >= hotTemperature) {
        throw std::invalid_argument("Cold temperature must be less than hot temperature");
    }
    return hotTemperature / (hotTemperature - coldTemperature);
}

// ============================================================================
// CONVECTION HEAT TRANSFER
// ============================================================================

/**
 * @brief Calculate convection heat transfer rate (Newton's law of cooling)
 *
 * Q/t = hA(T_surface - T_fluid)
 *
 * where h is convection heat transfer coefficient
 *
 * @param convectionCoefficient Heat transfer coefficient (W/(m²·K))
 * @param area Surface area (m²)
 * @param surfaceTemp Surface temperature (K or °C)
 * @param fluidTemp Fluid temperature (K or °C)
 * @return Heat transfer rate (W)
 * @throws std::invalid_argument if coefficient or area is non-positive
 */
inline double convectionHeatRate(double convectionCoefficient, double area,
                                 double surfaceTemp, double fluidTemp) {
    if (convectionCoefficient <= 0.0 || area <= 0.0) {
        throw std::invalid_argument("Coefficient and area must be positive");
    }
    return convectionCoefficient * area * (surfaceTemp - fluidTemp);
}

/**
 * @brief Calculate cooling/heating time constant
 *
 * τ = (mc)/(hA)
 *
 * Time constant for exponential cooling: T(t) = T_∞ + (T₀ - T_∞)e^(-t/τ)
 *
 * @param mass Mass of object (kg)
 * @param specificHeat Specific heat capacity (J/(kg·K))
 * @param convectionCoefficient Heat transfer coefficient (W/(m²·K))
 * @param area Surface area (m²)
 * @return Time constant (s)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double coolingTimeConstant(double mass, double specificHeat,
                                  double convectionCoefficient, double area) {
    if (mass <= 0.0 || specificHeat <= 0.0 || convectionCoefficient <= 0.0 || area <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return (mass * specificHeat) / (convectionCoefficient * area);
}

/**
 * @brief Calculate temperature after cooling/heating time
 *
 * T(t) = T_∞ + (T₀ - T_∞)e^(-t/τ)
 *
 * Newton's law of cooling (exponential decay)
 *
 * @param initialTemp Initial temperature (K or °C)
 * @param ambientTemp Ambient/fluid temperature (K or °C)
 * @param time Time elapsed (s)
 * @param timeConstant Time constant τ (s)
 * @return Temperature at time t (same unit as input temps)
 * @throws std::invalid_argument if time constant is non-positive
 */
inline double temperatureAfterCooling(double initialTemp, double ambientTemp,
                                      double time, double timeConstant) {
    if (timeConstant <= 0.0) {
        throw std::invalid_argument("Time constant must be positive");
    }
    if (time < 0.0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return ambientTemp + (initialTemp - ambientTemp) * std::exp(-time / timeConstant);
}

} // namespace heat_transfer
} // namespace physics

#endif // PHYSICS_HEAT_TRANSFER_HPP
