#ifndef PHYSICS_ELECTRIC_CIRCUITS_HPP
#define PHYSICS_ELECTRIC_CIRCUITS_HPP

#include <cmath>
#include <stdexcept>
#include <vector>

/**
 * @file electric_circuits.hpp
 * @brief Comprehensive implementation of electric circuit theory and analysis
 *
 * This module implements:
 * - Ohm's law and circuit analysis
 * - Resistance calculations (wire resistance, specific resistance)
 * - Series and parallel combinations
 * - Power dissipation and heating effects
 * - Battery circuits (EMF, internal resistance)
 * - Energy in electric circuits
 * - Kirchhoff's laws
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace electric_circuits {

/**
 * @namespace constants
 * @brief Physical constants and resistivities for common materials
 */
namespace constants {
    // Resistivities (ρ) in Ω⋅m at 20°C
    constexpr double SILVER_RESISTIVITY = 1.59e-8;          // Silver (best conductor)
    constexpr double COPPER_RESISTIVITY = 1.68e-8;          // Copper
    constexpr double GOLD_RESISTIVITY = 2.44e-8;            // Gold
    constexpr double ALUMINUM_RESISTIVITY = 2.82e-8;        // Aluminum
    constexpr double TUNGSTEN_RESISTIVITY = 5.60e-8;        // Tungsten (used in bulbs)
    constexpr double IRON_RESISTIVITY = 1.0e-7;             // Iron
    constexpr double NICHROME_RESISTIVITY = 1.10e-6;        // Nichrome (heating elements)
    constexpr double CARBON_RESISTIVITY = 3.5e-5;           // Carbon

    // Temperature coefficient of resistance (α) in K⁻¹
    constexpr double COPPER_TEMP_COEFF = 0.00393;
    constexpr double ALUMINUM_TEMP_COEFF = 0.00429;
    constexpr double TUNGSTEN_TEMP_COEFF = 0.0045;
    constexpr double NICHROME_TEMP_COEFF = 0.0004;
}

// ============================================================================
// OHM'S LAW AND BASIC CIRCUIT RELATIONS
// ============================================================================

/**
 * @brief Calculate voltage using Ohm's law
 *
 * V = I × R
 *
 * @param current Current flowing through resistor (A)
 * @param resistance Resistance (Ω)
 * @return Voltage across resistor (V)
 * @throws std::invalid_argument if resistance is negative
 */
inline double ohmsLawVoltage(double current, double resistance) {
    if (resistance < 0.0) {
        throw std::invalid_argument("Resistance cannot be negative");
    }
    return current * resistance;
}

/**
 * @brief Calculate current using Ohm's law
 *
 * I = V / R
 *
 * @param voltage Voltage across resistor (V)
 * @param resistance Resistance (Ω)
 * @return Current through resistor (A)
 * @throws std::invalid_argument if resistance is non-positive
 */
inline double ohmsLawCurrent(double voltage, double resistance) {
    if (resistance <= 0.0) {
        throw std::invalid_argument("Resistance must be positive");
    }
    return voltage / resistance;
}

/**
 * @brief Calculate resistance using Ohm's law
 *
 * R = V / I
 *
 * @param voltage Voltage across resistor (V)
 * @param current Current through resistor (A)
 * @return Resistance (Ω)
 * @throws std::invalid_argument if current is zero
 */
inline double ohmsLawResistance(double voltage, double current) {
    if (std::abs(current) < 1e-15) {
        throw std::invalid_argument("Current must be non-zero");
    }
    return voltage / current;
}

// ============================================================================
// RESISTANCE OF WIRES AND SPECIFIC RESISTANCE
// ============================================================================

/**
 * @brief Calculate resistance of a wire
 *
 * R = ρL/A
 *
 * @param resistivity Resistivity of material (Ω⋅m)
 * @param length Length of wire (m)
 * @param crossSectionalArea Cross-sectional area (m²)
 * @return Resistance (Ω)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double wireResistance(double resistivity, double length, double crossSectionalArea) {
    if (resistivity < 0.0) {
        throw std::invalid_argument("Resistivity cannot be negative");
    }
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (crossSectionalArea <= 0.0) {
        throw std::invalid_argument("Cross-sectional area must be positive");
    }
    return resistivity * length / crossSectionalArea;
}

/**
 * @brief Calculate resistance of circular wire from diameter
 *
 * R = ρL/(πd²/4) = 4ρL/(πd²)
 *
 * @param resistivity Resistivity of material (Ω⋅m)
 * @param length Length of wire (m)
 * @param diameter Diameter of wire (m)
 * @return Resistance (Ω)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double circularWireResistance(double resistivity, double length, double diameter) {
    if (resistivity < 0.0) {
        throw std::invalid_argument("Resistivity cannot be negative");
    }
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (diameter <= 0.0) {
        throw std::invalid_argument("Diameter must be positive");
    }
    double area = M_PI * diameter * diameter / 4.0;
    return resistivity * length / area;
}

/**
 * @brief Calculate resistivity from resistance measurements
 *
 * ρ = RA/L
 *
 * @param resistance Measured resistance (Ω)
 * @param crossSectionalArea Cross-sectional area (m²)
 * @param length Length of wire (m)
 * @return Resistivity (Ω⋅m)
 * @throws std::invalid_argument if length is non-positive
 */
inline double calculateResistivity(double resistance, double crossSectionalArea, double length) {
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (crossSectionalArea <= 0.0) {
        throw std::invalid_argument("Cross-sectional area must be positive");
    }
    return resistance * crossSectionalArea / length;
}

/**
 * @brief Calculate resistance at different temperature
 *
 * R_T = R₀(1 + α(T - T₀))
 *
 * @param resistanceAtRef Resistance at reference temperature (Ω)
 * @param tempCoefficient Temperature coefficient of resistance (K⁻¹)
 * @param temperature Final temperature (K or °C)
 * @param referenceTemp Reference temperature (K or °C)
 * @return Resistance at final temperature (Ω)
 * @throws std::invalid_argument if reference resistance is negative
 */
inline double resistanceAtTemperature(double resistanceAtRef, double tempCoefficient,
                                      double temperature, double referenceTemp = 293.15) {
    if (resistanceAtRef < 0.0) {
        throw std::invalid_argument("Reference resistance cannot be negative");
    }
    return resistanceAtRef * (1.0 + tempCoefficient * (temperature - referenceTemp));
}

// ============================================================================
// SERIES AND PARALLEL COMBINATIONS
// ============================================================================

/**
 * @brief Calculate total resistance of resistors in series
 *
 * R_total = R₁ + R₂ + R₃ + ...
 *
 * @param resistances Vector of resistances (Ω)
 * @return Total resistance (Ω)
 * @throws std::invalid_argument if any resistance is negative or vector is empty
 */
inline double seriesResistance(const std::vector<double>& resistances) {
    if (resistances.empty()) {
        throw std::invalid_argument("Resistances vector cannot be empty");
    }
    double total = 0.0;
    for (double r : resistances) {
        if (r < 0.0) {
            throw std::invalid_argument("Resistance cannot be negative");
        }
        total += r;
    }
    return total;
}

/**
 * @brief Calculate total resistance of resistors in parallel
 *
 * 1/R_total = 1/R₁ + 1/R₂ + 1/R₃ + ...
 *
 * @param resistances Vector of resistances (Ω)
 * @return Total resistance (Ω)
 * @throws std::invalid_argument if any resistance is non-positive or vector is empty
 */
inline double parallelResistance(const std::vector<double>& resistances) {
    if (resistances.empty()) {
        throw std::invalid_argument("Resistances vector cannot be empty");
    }
    double reciprocalSum = 0.0;
    for (double r : resistances) {
        if (r <= 0.0) {
            throw std::invalid_argument("Resistance must be positive for parallel combination");
        }
        reciprocalSum += 1.0 / r;
    }
    return 1.0 / reciprocalSum;
}

/**
 * @brief Calculate total resistance of two resistors in parallel (simplified)
 *
 * R = (R₁ × R₂)/(R₁ + R₂)
 *
 * @param r1 First resistance (Ω)
 * @param r2 Second resistance (Ω)
 * @return Total resistance (Ω)
 * @throws std::invalid_argument if resistances are non-positive
 */
inline double twoResistorsParallel(double r1, double r2) {
    if (r1 <= 0.0 || r2 <= 0.0) {
        throw std::invalid_argument("Resistances must be positive");
    }
    return (r1 * r2) / (r1 + r2);
}

// ============================================================================
// POWER AND ENERGY
// ============================================================================

/**
 * @brief Calculate electrical power (P = VI)
 *
 * P = V × I
 *
 * @param voltage Voltage (V)
 * @param current Current (A)
 * @return Power (W)
 */
inline double electricalPower(double voltage, double current) {
    return voltage * current;
}

/**
 * @brief Calculate power dissipated in resistor (P = I²R)
 *
 * P = I² × R
 *
 * @param current Current through resistor (A)
 * @param resistance Resistance (Ω)
 * @return Power dissipated (W)
 * @throws std::invalid_argument if resistance is negative
 */
inline double powerFromCurrent(double current, double resistance) {
    if (resistance < 0.0) {
        throw std::invalid_argument("Resistance cannot be negative");
    }
    return current * current * resistance;
}

/**
 * @brief Calculate power dissipated in resistor (P = V²/R)
 *
 * P = V² / R
 *
 * @param voltage Voltage across resistor (V)
 * @param resistance Resistance (Ω)
 * @return Power dissipated (W)
 * @throws std::invalid_argument if resistance is non-positive
 */
inline double powerFromVoltage(double voltage, double resistance) {
    if (resistance <= 0.0) {
        throw std::invalid_argument("Resistance must be positive");
    }
    return (voltage * voltage) / resistance;
}

/**
 * @brief Calculate energy dissipated over time
 *
 * E = P × t = V × I × t
 *
 * @param power Power (W)
 * @param time Time duration (s)
 * @return Energy (J)
 * @throws std::invalid_argument if time is negative
 */
inline double energyDissipated(double power, double time) {
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    return power * time;
}

/**
 * @brief Calculate energy from current and resistance (E = I²Rt)
 *
 * E = I² × R × t
 *
 * Heating effect of current (Joule heating)
 *
 * @param current Current (A)
 * @param resistance Resistance (Ω)
 * @param time Time duration (s)
 * @return Energy dissipated as heat (J)
 * @throws std::invalid_argument if resistance or time is negative
 */
inline double jouleHeating(double current, double resistance, double time) {
    if (resistance < 0.0) {
        throw std::invalid_argument("Resistance cannot be negative");
    }
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    return current * current * resistance * time;
}

/**
 * @brief Calculate heat generated in calories
 *
 * H = 0.24 × I² × R × t (calories)
 *
 * Where 0.24 is the mechanical equivalent of heat (cal/J)
 *
 * @param current Current (A)
 * @param resistance Resistance (Ω)
 * @param time Time (s)
 * @return Heat in calories (cal)
 * @throws std::invalid_argument if resistance or time is negative
 */
inline double heatingEffectCalories(double current, double resistance, double time) {
    if (resistance < 0.0) {
        throw std::invalid_argument("Resistance cannot be negative");
    }
    if (time < 0.0) {
        throw std::invalid_argument("Time cannot be negative");
    }
    constexpr double JOULES_TO_CALORIES = 0.239006; // More precise value
    return JOULES_TO_CALORIES * current * current * resistance * time;
}

// ============================================================================
// CELLS AND BATTERIES
// ============================================================================

/**
 * @brief Calculate terminal voltage of cell with internal resistance
 *
 * V = ε - I × r
 *
 * where ε is EMF and r is internal resistance
 *
 * @param emf Electromotive force (V)
 * @param current Current drawn (A)
 * @param internalResistance Internal resistance of cell (Ω)
 * @return Terminal voltage (V)
 * @throws std::invalid_argument if internal resistance is negative
 */
inline double terminalVoltage(double emf, double current, double internalResistance) {
    if (internalResistance < 0.0) {
        throw std::invalid_argument("Internal resistance cannot be negative");
    }
    return emf - current * internalResistance;
}

/**
 * @brief Calculate current in simple circuit with cell
 *
 * I = ε / (R + r)
 *
 * @param emf Electromotive force (V)
 * @param externalResistance External resistance (Ω)
 * @param internalResistance Internal resistance (Ω)
 * @return Current (A)
 * @throws std::invalid_argument if total resistance is non-positive
 */
inline double cellCurrent(double emf, double externalResistance, double internalResistance) {
    if (externalResistance < 0.0 || internalResistance < 0.0) {
        throw std::invalid_argument("Resistances cannot be negative");
    }
    double totalResistance = externalResistance + internalResistance;
    if (totalResistance <= 0.0) {
        throw std::invalid_argument("Total resistance must be positive");
    }
    return emf / totalResistance;
}

/**
 * @brief Calculate EMF of cells in series
 *
 * ε_total = ε₁ + ε₂ + ε₃ + ...
 * r_total = r₁ + r₂ + r₃ + ...
 *
 * @param emfs Vector of EMFs (V)
 * @return Total EMF (V)
 * @throws std::invalid_argument if vector is empty
 */
inline double cellsInSeriesEMF(const std::vector<double>& emfs) {
    if (emfs.empty()) {
        throw std::invalid_argument("EMF vector cannot be empty");
    }
    double total = 0.0;
    for (double e : emfs) {
        total += e;
    }
    return total;
}

/**
 * @brief Calculate internal resistance of cells in series
 *
 * r_total = r₁ + r₂ + r₃ + ...
 *
 * @param internalResistances Vector of internal resistances (Ω)
 * @return Total internal resistance (Ω)
 * @throws std::invalid_argument if vector is empty
 */
inline double cellsInSeriesResistance(const std::vector<double>& internalResistances) {
    return seriesResistance(internalResistances);
}

/**
 * @brief Calculate current from n identical cells in series
 *
 * I = nε / (R + nr)
 *
 * @param numberOfCells Number of identical cells
 * @param cellEMF EMF of each cell (V)
 * @param cellInternalResistance Internal resistance of each cell (Ω)
 * @param externalResistance External resistance (Ω)
 * @return Current (A)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double identicalCellsSeriesCurrent(int numberOfCells, double cellEMF,
                                          double cellInternalResistance,
                                          double externalResistance) {
    if (numberOfCells <= 0) {
        throw std::invalid_argument("Number of cells must be positive");
    }
    if (cellInternalResistance < 0.0 || externalResistance < 0.0) {
        throw std::invalid_argument("Resistances cannot be negative");
    }
    double totalResistance = externalResistance + numberOfCells * cellInternalResistance;
    if (totalResistance <= 0.0) {
        throw std::invalid_argument("Total resistance must be positive");
    }
    return (numberOfCells * cellEMF) / totalResistance;
}

/**
 * @brief Calculate current from n identical cells in parallel
 *
 * I = ε / (R + r/n)
 *
 * @param numberOfCells Number of identical cells in parallel
 * @param cellEMF EMF of each cell (V)
 * @param cellInternalResistance Internal resistance of each cell (Ω)
 * @param externalResistance External resistance (Ω)
 * @return Current (A)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double identicalCellsParallelCurrent(int numberOfCells, double cellEMF,
                                            double cellInternalResistance,
                                            double externalResistance) {
    if (numberOfCells <= 0) {
        throw std::invalid_argument("Number of cells must be positive");
    }
    if (cellInternalResistance < 0.0 || externalResistance < 0.0) {
        throw std::invalid_argument("Resistances cannot be negative");
    }
    double effectiveInternalR = cellInternalResistance / numberOfCells;
    double totalResistance = externalResistance + effectiveInternalR;
    if (totalResistance <= 0.0) {
        throw std::invalid_argument("Total resistance must be positive");
    }
    return cellEMF / totalResistance;
}

/**
 * @brief Calculate current from m rows of n cells (mixed arrangement)
 *
 * m rows in parallel, each row has n cells in series
 * I = nε / (R + nr/m)
 *
 * @param rows Number of parallel rows (m)
 * @param cellsPerRow Number of cells in each row (n)
 * @param cellEMF EMF of each cell (V)
 * @param cellInternalResistance Internal resistance of each cell (Ω)
 * @param externalResistance External resistance (Ω)
 * @return Current (A)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double mixedCellsCurrent(int rows, int cellsPerRow, double cellEMF,
                                double cellInternalResistance, double externalResistance) {
    if (rows <= 0 || cellsPerRow <= 0) {
        throw std::invalid_argument("Number of rows and cells per row must be positive");
    }
    if (cellInternalResistance < 0.0 || externalResistance < 0.0) {
        throw std::invalid_argument("Resistances cannot be negative");
    }
    double totalEMF = cellsPerRow * cellEMF;
    double effectiveInternalR = (cellsPerRow * cellInternalResistance) / rows;
    double totalResistance = externalResistance + effectiveInternalR;
    if (totalResistance <= 0.0) {
        throw std::invalid_argument("Total resistance must be positive");
    }
    return totalEMF / totalResistance;
}

/**
 * @brief Calculate optimal number of cells in series for maximum power
 *
 * For maximum power transfer: R = nr
 * Optimal n = R/r
 *
 * @param externalResistance External resistance (Ω)
 * @param cellInternalResistance Internal resistance of each cell (Ω)
 * @return Optimal number of cells in series
 * @throws std::invalid_argument if internal resistance is non-positive
 */
inline int optimalCellsInSeries(double externalResistance, double cellInternalResistance) {
    if (cellInternalResistance <= 0.0) {
        throw std::invalid_argument("Cell internal resistance must be positive");
    }
    if (externalResistance < 0.0) {
        throw std::invalid_argument("External resistance cannot be negative");
    }
    return static_cast<int>(std::round(externalResistance / cellInternalResistance));
}

// ============================================================================
// KIRCHHOFF'S LAWS HELPERS
// ============================================================================

/**
 * @brief Calculate voltage drop across part of circuit
 *
 * V_drop = I × R
 *
 * Ohm's law applied to part of circuit
 *
 * @param current Current through the section (A)
 * @param resistance Resistance of the section (Ω)
 * @return Voltage drop (V)
 * @throws std::invalid_argument if resistance is negative
 */
inline double voltageDrop(double current, double resistance) {
    if (resistance < 0.0) {
        throw std::invalid_argument("Resistance cannot be negative");
    }
    return current * resistance;
}

/**
 * @brief Calculate power delivered by source
 *
 * P_source = ε × I
 *
 * @param emf EMF of source (V)
 * @param current Current supplied (A)
 * @return Power delivered (W)
 */
inline double sourcePower(double emf, double current) {
    return emf * current;
}

/**
 * @brief Calculate power dissipated in internal resistance
 *
 * P_internal = I² × r
 *
 * @param current Current (A)
 * @param internalResistance Internal resistance (Ω)
 * @return Power dissipated internally (W)
 * @throws std::invalid_argument if internal resistance is negative
 */
inline double internalPowerLoss(double current, double internalResistance) {
    if (internalResistance < 0.0) {
        throw std::invalid_argument("Internal resistance cannot be negative");
    }
    return current * current * internalResistance;
}

/**
 * @brief Calculate efficiency of power transfer
 *
 * η = P_external / P_total = R / (R + r)
 *
 * @param externalResistance External load resistance (Ω)
 * @param internalResistance Internal resistance (Ω)
 * @return Efficiency (0 to 1)
 * @throws std::invalid_argument if resistances are negative
 */
inline double powerTransferEfficiency(double externalResistance, double internalResistance) {
    if (externalResistance < 0.0 || internalResistance < 0.0) {
        throw std::invalid_argument("Resistances cannot be negative");
    }
    double totalResistance = externalResistance + internalResistance;
    if (totalResistance <= 0.0) {
        return 0.0;
    }
    return externalResistance / totalResistance;
}

/**
 * @brief Calculate maximum power delivered to external load
 *
 * P_max = ε² / (4r)
 *
 * Maximum power occurs when R = r (matched load)
 *
 * @param emf EMF of source (V)
 * @param internalResistance Internal resistance (Ω)
 * @return Maximum power that can be delivered (W)
 * @throws std::invalid_argument if internal resistance is non-positive
 */
inline double maximumPowerTransfer(double emf, double internalResistance) {
    if (internalResistance <= 0.0) {
        throw std::invalid_argument("Internal resistance must be positive");
    }
    return (emf * emf) / (4.0 * internalResistance);
}

} // namespace electric_circuits
} // namespace physics

#endif // PHYSICS_ELECTRIC_CIRCUITS_HPP
