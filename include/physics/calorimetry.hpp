#ifndef PHYSICS_CALORIMETRY_HPP
#define PHYSICS_CALORIMETRY_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file calorimetry.hpp
 * @brief Comprehensive implementation of calorimetry and heat transfer
 *
 * This module implements:
 * - Heat capacity and specific heat
 * - Method of mixtures
 * - Phase changes (latent heat)
 * - Electrical calorimeter
 * - Heat transfer applications
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace calorimetry {

/**
 * @namespace constants
 * @brief Physical constants for calorimetry calculations
 */

namespace constants {
    // Specific heat capacities (J/(kg·K))
    constexpr double WATER_SPECIFIC_HEAT = 4186.0;         // Water
    constexpr double ICE_SPECIFIC_HEAT = 2090.0;           // Ice
    constexpr double STEAM_SPECIFIC_HEAT = 2010.0;         // Steam
    constexpr double ALUMINUM_SPECIFIC_HEAT = 900.0;       // Aluminum
    constexpr double COPPER_SPECIFIC_HEAT = 385.0;         // Copper
    constexpr double IRON_SPECIFIC_HEAT = 450.0;           // Iron
    constexpr double STEEL_SPECIFIC_HEAT = 420.0;          // Steel
    constexpr double GLASS_SPECIFIC_HEAT = 840.0;          // Glass
    constexpr double CONCRETE_SPECIFIC_HEAT = 880.0;       // Concrete
    constexpr double ETHANOL_SPECIFIC_HEAT = 2440.0;       // Ethanol

    // Latent heats (J/kg)
    constexpr double WATER_LATENT_FUSION = 334000.0;       // Ice → Water (melting)
    constexpr double WATER_LATENT_VAPORIZATION = 2260000.0; // Water → Steam (boiling)
    constexpr double ETHANOL_LATENT_VAPORIZATION = 855000.0; // Ethanol boiling

    // Standard temperatures (K)
    constexpr double WATER_FREEZING_POINT = 273.15;        // 0°C
    constexpr double WATER_BOILING_POINT = 373.15;         // 100°C
}

// ============================================================================
// HEAT CAPACITY AND SPECIFIC HEAT
// ============================================================================

/**
 * @brief Calculate heat transferred
 *
 * Q = mcΔT
 *
 * @param mass Mass of substance (kg)
 * @param specificHeat Specific heat capacity (J/(kg·K))
 * @param temperatureChange Change in temperature (K or °C)
 * @return Heat transferred (J), positive = heat absorbed
 * @throws std::invalid_argument if mass or specific heat is non-positive
 */
inline double calculateHeat(double mass, double specificHeat, double temperatureChange) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (specificHeat <= 0.0) {
        throw std::invalid_argument("Specific heat must be positive");
    }
    return mass * specificHeat * temperatureChange;
}

/**
 * @brief Calculate specific heat capacity from measurements
 *
 * c = Q/(mΔT)
 *
 * @param heat Heat transferred (J)
 * @param mass Mass of substance (kg)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Specific heat capacity (J/(kg·K))
 * @throws std::invalid_argument if mass or temperature change is zero
 */
inline double calculateSpecificHeat(double heat, double mass, double temperatureChange) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (std::abs(temperatureChange) < 1e-10) {
        throw std::invalid_argument("Temperature change must be non-zero");
    }
    return heat / (mass * temperatureChange);
}

/**
 * @brief Calculate heat capacity of an object
 *
 * C = mc
 *
 * Heat capacity is the product of mass and specific heat
 *
 * @param mass Mass of object (kg)
 * @param specificHeat Specific heat capacity (J/(kg·K))
 * @return Heat capacity (J/K)
 * @throws std::invalid_argument if mass or specific heat is non-positive
 */
inline double calculateHeatCapacity(double mass, double specificHeat) {
    if (mass <= 0.0 || specificHeat <= 0.0) {
        throw std::invalid_argument("Mass and specific heat must be positive");
    }
    return mass * specificHeat;
}

/**
 * @brief Calculate final temperature after heat transfer
 *
 * T_f = T_i + Q/(mc)
 *
 * @param initialTemp Initial temperature (K or °C)
 * @param heat Heat transferred (J), positive = heat absorbed
 * @param mass Mass of substance (kg)
 * @param specificHeat Specific heat capacity (J/(kg·K))
 * @return Final temperature (same unit as initialTemp)
 * @throws std::invalid_argument if mass or specific heat is non-positive
 */
inline double calculateFinalTemperature(double initialTemp, double heat,
                                        double mass, double specificHeat) {
    if (mass <= 0.0 || specificHeat <= 0.0) {
        throw std::invalid_argument("Mass and specific heat must be positive");
    }
    return initialTemp + heat / (mass * specificHeat);
}

// ============================================================================
// METHOD OF MIXTURES
// ============================================================================

/**
 * @brief Calculate equilibrium temperature of two substances mixing
 *
 * Method of Mixtures: Heat lost = Heat gained
 * m₁c₁(T_eq - T₁) = -m₂c₂(T_eq - T₂)
 *
 * T_eq = (m₁c₁T₁ + m₂c₂T₂)/(m₁c₁ + m₂c₂)
 *
 * @param mass1 Mass of first substance (kg)
 * @param specificHeat1 Specific heat of first substance (J/(kg·K))
 * @param temp1 Initial temperature of first substance (K or °C)
 * @param mass2 Mass of second substance (kg)
 * @param specificHeat2 Specific heat of second substance (J/(kg·K))
 * @param temp2 Initial temperature of second substance (K or °C)
 * @return Equilibrium temperature (same unit as input temps)
 * @throws std::invalid_argument if masses or specific heats are non-positive
 */
inline double equilibriumTemperature(double mass1, double specificHeat1, double temp1,
                                     double mass2, double specificHeat2, double temp2) {
    if (mass1 <= 0.0 || mass2 <= 0.0) {
        throw std::invalid_argument("Masses must be positive");
    }
    if (specificHeat1 <= 0.0 || specificHeat2 <= 0.0) {
        throw std::invalid_argument("Specific heats must be positive");
    }

    double heatCap1 = mass1 * specificHeat1;
    double heatCap2 = mass2 * specificHeat2;

    return (heatCap1 * temp1 + heatCap2 * temp2) / (heatCap1 + heatCap2);
}

/**
 * @brief Calculate heat exchanged between two substances
 *
 * Q = m₁c₁(T_eq - T₁) = -m₂c₂(T_eq - T₂)
 *
 * @param mass Mass of substance (kg)
 * @param specificHeat Specific heat (J/(kg·K))
 * @param initialTemp Initial temperature (K or °C)
 * @param equilibriumTemp Equilibrium temperature (K or °C)
 * @return Heat transferred (J), positive = heat absorbed
 * @throws std::invalid_argument if mass or specific heat is non-positive
 */
inline double heatExchanged(double mass, double specificHeat,
                            double initialTemp, double equilibriumTemp) {
    if (mass <= 0.0 || specificHeat <= 0.0) {
        throw std::invalid_argument("Mass and specific heat must be positive");
    }
    return mass * specificHeat * (equilibriumTemp - initialTemp);
}

/**
 * @brief Calculate unknown specific heat using method of mixtures
 *
 * From heat balance: m₁c₁(T_eq - T₁) + m₂c₂(T_eq - T₂) = 0
 *
 * @param mass Mass of unknown substance (kg)
 * @param initialTemp Initial temperature of unknown (K or °C)
 * @param equilibriumTemp Final equilibrium temperature (K or °C)
 * @param knownMass Mass of known substance (kg)
 * @param knownSpecificHeat Specific heat of known substance (J/(kg·K))
 * @param knownInitialTemp Initial temperature of known substance (K or °C)
 * @return Specific heat of unknown substance (J/(kg·K))
 * @throws std::invalid_argument if parameters are invalid
 */
inline double findUnknownSpecificHeat(double mass, double initialTemp, double equilibriumTemp,
                                      double knownMass, double knownSpecificHeat,
                                      double knownInitialTemp) {
    if (mass <= 0.0 || knownMass <= 0.0) {
        throw std::invalid_argument("Masses must be positive");
    }
    if (knownSpecificHeat <= 0.0) {
        throw std::invalid_argument("Known specific heat must be positive");
    }

    double tempChange = equilibriumTemp - initialTemp;
    if (std::abs(tempChange) < 1e-10) {
        throw std::invalid_argument("Temperature change for unknown must be non-zero");
    }

    // Heat gained by unknown = -Heat lost by known
    double heatFromKnown = knownMass * knownSpecificHeat * (equilibriumTemp - knownInitialTemp);
    return -heatFromKnown / (mass * tempChange);
}

// ============================================================================
// PHASE CHANGES (LATENT HEAT)
// ============================================================================

/**
 * @brief Calculate heat required for phase change
 *
 * Q = mL
 *
 * where L is latent heat (fusion or vaporization)
 *
 * @param mass Mass of substance (kg)
 * @param latentHeat Latent heat of phase change (J/kg)
 * @return Heat required (J)
 * @throws std::invalid_argument if mass or latent heat is non-positive
 */
inline double calculateLatentHeat(double mass, double latentHeat) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (latentHeat <= 0.0) {
        throw std::invalid_argument("Latent heat must be positive");
    }
    return mass * latentHeat;
}

/**
 * @brief Calculate total heat to raise temperature and change phase
 *
 * Example: Ice at -10°C to water at 20°C
 * Q = Q_heating_ice + Q_melting + Q_heating_water
 *
 * @param mass Mass of substance (kg)
 * @param specificHeat1 Specific heat in initial phase (J/(kg·K))
 * @param tempChange1 Temperature change in initial phase (K or °C)
 * @param latentHeat Latent heat of phase transition (J/kg)
 * @param specificHeat2 Specific heat in final phase (J/(kg·K))
 * @param tempChange2 Temperature change in final phase (K or °C)
 * @return Total heat required (J)
 * @throws std::invalid_argument if mass or specific heats are non-positive
 */
inline double totalHeatWithPhaseChange(double mass,
                                       double specificHeat1, double tempChange1,
                                       double latentHeat,
                                       double specificHeat2, double tempChange2) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (specificHeat1 <= 0.0 || specificHeat2 <= 0.0) {
        throw std::invalid_argument("Specific heats must be positive");
    }

    double q1 = mass * specificHeat1 * tempChange1;    // Heat to reach phase transition temp
    double q2 = mass * latentHeat;                      // Heat for phase change
    double q3 = mass * specificHeat2 * tempChange2;    // Heat in new phase

    return q1 + q2 + q3;
}

// ============================================================================
// ELECTRICAL CALORIMETER
// ============================================================================

/**
 * @brief Calculate heat generated by electrical current
 *
 * Q = I²Rt = VIt = Pt
 *
 * Joule heating / Electrical work converted to heat
 *
 * @param current Electric current (A)
 * @param resistance Resistance of heater (Ω)
 * @param time Time duration (s)
 * @return Heat generated (J)
 * @throws std::invalid_argument if current, resistance, or time is negative
 */
inline double electricalHeatGenerated(double current, double resistance, double time) {
    if (current < 0.0 || resistance < 0.0 || time < 0.0) {
        throw std::invalid_argument("Current, resistance, and time must be non-negative");
    }
    return current * current * resistance * time;
}

/**
 * @brief Calculate heat generated from voltage and current
 *
 * Q = VIt
 *
 * @param voltage Voltage across heater (V)
 * @param current Electric current (A)
 * @param time Time duration (s)
 * @return Heat generated (J)
 * @throws std::invalid_argument if parameters are negative
 */
inline double electricalHeatFromVoltage(double voltage, double current, double time) {
    if (voltage < 0.0 || current < 0.0 || time < 0.0) {
        throw std::invalid_argument("Voltage, current, and time must be non-negative");
    }
    return voltage * current * time;
}

/**
 * @brief Calculate specific heat using electrical calorimeter
 *
 * Energy supplied electrically = Heat absorbed by substance
 * VIt = mcΔT
 *
 * @param voltage Voltage (V)
 * @param current Current (A)
 * @param time Time of heating (s)
 * @param mass Mass of substance (kg)
 * @param temperatureChange Temperature rise (K or °C)
 * @return Specific heat capacity (J/(kg·K))
 * @throws std::invalid_argument if parameters are invalid
 */
inline double specificHeatElectricalMethod(double voltage, double current, double time,
                                           double mass, double temperatureChange) {
    if (mass <= 0.0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (std::abs(temperatureChange) < 1e-10) {
        throw std::invalid_argument("Temperature change must be non-zero");
    }

    double energySupplied = voltage * current * time;
    return energySupplied / (mass * temperatureChange);
}

// ============================================================================
// CALORIMETER CALCULATIONS
// ============================================================================

/**
 * @brief Calculate water equivalent of calorimeter
 *
 * Water equivalent W = (mass of calorimeter) × (specific heat of calorimeter) / (specific heat of water)
 *
 * The mass of water that would absorb the same heat as the calorimeter
 *
 * @param calorimeterMass Mass of calorimeter (kg)
 * @param calorimeterSpecificHeat Specific heat of calorimeter material (J/(kg·K))
 * @param waterSpecificHeat Specific heat of water (J/(kg·K))
 * @return Water equivalent (kg)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double waterEquivalent(double calorimeterMass, double calorimeterSpecificHeat,
                              double waterSpecificHeat = constants::WATER_SPECIFIC_HEAT) {
    if (calorimeterMass <= 0.0 || calorimeterSpecificHeat <= 0.0 || waterSpecificHeat <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return (calorimeterMass * calorimeterSpecificHeat) / waterSpecificHeat;
}

/**
 * @brief Calculate equilibrium temperature accounting for calorimeter
 *
 * Includes heat absorbed by calorimeter:
 * T_eq = (m₁c₁T₁ + m₂c₂T₂ + m_cal·c_cal·T_cal)/(m₁c₁ + m₂c₂ + m_cal·c_cal)
 *
 * @param mass Mass of added substance (kg)
 * @param specificHeat Specific heat of substance (J/(kg·K))
 * @param temp Initial temperature of substance (K or °C)
 * @param waterMass Mass of water in calorimeter (kg)
 * @param waterTemp Initial water temperature (K or °C)
 * @param calorimeterMass Mass of calorimeter (kg)
 * @param calorimeterSpecificHeat Specific heat of calorimeter (J/(kg·K))
 * @param waterSpecificHeat Specific heat of water (J/(kg·K))
 * @return Equilibrium temperature (same unit as input temps)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double equilibriumWithCalorimeter(double mass, double specificHeat, double temp,
                                         double waterMass, double waterTemp,
                                         double calorimeterMass, double calorimeterSpecificHeat,
                                         double waterSpecificHeat = constants::WATER_SPECIFIC_HEAT) {
    if (mass <= 0.0 || waterMass <= 0.0 || calorimeterMass <= 0.0) {
        throw std::invalid_argument("Masses must be positive");
    }
    if (specificHeat <= 0.0 || waterSpecificHeat <= 0.0 || calorimeterSpecificHeat <= 0.0) {
        throw std::invalid_argument("Specific heats must be positive");
    }

    double heatCap1 = mass * specificHeat;
    double heatCap2 = waterMass * waterSpecificHeat;
    double heatCapCal = calorimeterMass * calorimeterSpecificHeat;

    return (heatCap1 * temp + heatCap2 * waterTemp + heatCapCal * waterTemp) /
           (heatCap1 + heatCap2 + heatCapCal);
}

// ============================================================================
// HEAT TRANSFER RATE
// ============================================================================

/**
 * @brief Calculate heat transfer rate (power)
 *
 * P = Q/t
 *
 * @param heat Total heat transferred (J)
 * @param time Time duration (s)
 * @return Heat transfer rate (W)
 * @throws std::invalid_argument if time is non-positive
 */
inline double heatTransferRate(double heat, double time) {
    if (time <= 0.0) {
        throw std::invalid_argument("Time must be positive");
    }
    return heat / time;
}

/**
 * @brief Calculate time required for heat transfer
 *
 * t = Q/P
 *
 * @param heat Total heat to transfer (J)
 * @param power Heat transfer rate (W)
 * @return Time required (s)
 * @throws std::invalid_argument if power is non-positive
 */
inline double timeForHeatTransfer(double heat, double power) {
    if (power <= 0.0) {
        throw std::invalid_argument("Power must be positive");
    }
    return heat / power;
}

} // namespace calorimetry
} // namespace physics

#endif // PHYSICS_CALORIMETRY_HPP
