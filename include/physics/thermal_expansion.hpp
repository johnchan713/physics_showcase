#ifndef PHYSICS_THERMAL_EXPANSION_HPP
#define PHYSICS_THERMAL_EXPANSION_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file thermal_expansion.hpp
 * @brief Comprehensive implementation of thermal expansion of solids and liquids
 *
 * This module implements:
 * - Linear expansion of solids
 * - Volume expansion (solids and liquids)
 * - Coefficient of expansion measurements
 * - Practical applications
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace thermal_expansion {

/**
 * @namespace constants
 * @brief Coefficients of expansion for common materials
 */
namespace constants {
    // Coefficients of linear expansion (α) in K⁻¹ or °C⁻¹
    constexpr double STEEL_LINEAR_EXPANSION = 11e-6;        // Steel
    constexpr double ALUMINUM_LINEAR_EXPANSION = 23e-6;     // Aluminum
    constexpr double COPPER_LINEAR_EXPANSION = 17e-6;       // Copper
    constexpr double BRASS_LINEAR_EXPANSION = 19e-6;        // Brass
    constexpr double GLASS_LINEAR_EXPANSION = 9e-6;         // Glass (typical)
    constexpr double CONCRETE_LINEAR_EXPANSION = 12e-6;     // Concrete
    constexpr double IRON_LINEAR_EXPANSION = 12e-6;         // Iron
    constexpr double INVAR_LINEAR_EXPANSION = 1.2e-6;       // Invar (low expansion alloy)

    // Coefficients of volume expansion (β) in K⁻¹ or °C⁻¹
    // For isotropic solids: β ≈ 3α
    constexpr double STEEL_VOLUME_EXPANSION = 33e-6;
    constexpr double ALUMINUM_VOLUME_EXPANSION = 69e-6;
    constexpr double COPPER_VOLUME_EXPANSION = 51e-6;

    // Liquids (volume expansion only)
    constexpr double WATER_VOLUME_EXPANSION = 207e-6;       // Water (at 20°C)
    constexpr double MERCURY_VOLUME_EXPANSION = 182e-6;     // Mercury
    constexpr double ETHANOL_VOLUME_EXPANSION = 1100e-6;    // Ethanol
    constexpr double GLYCERIN_VOLUME_EXPANSION = 485e-6;    // Glycerin
}

// ============================================================================
// LINEAR EXPANSION OF SOLIDS
// ============================================================================

/**
 * @brief Calculate change in length due to thermal expansion
 *
 * ΔL = α × L₀ × ΔT
 *
 * @param originalLength Original length at initial temperature (m)
 * @param coefficientLinear Coefficient of linear expansion (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Change in length (m)
 * @throws std::invalid_argument if original length is non-positive
 */
inline double calculateLengthChange(double originalLength, double coefficientLinear,
                                    double temperatureChange) {
    if (originalLength <= 0.0) {
        throw std::invalid_argument("Original length must be positive");
    }
    return coefficientLinear * originalLength * temperatureChange;
}

/**
 * @brief Calculate final length after thermal expansion
 *
 * L = L₀(1 + αΔT)
 *
 * @param originalLength Original length (m)
 * @param coefficientLinear Coefficient of linear expansion (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Final length (m)
 * @throws std::invalid_argument if original length is non-positive
 */
inline double calculateFinalLength(double originalLength, double coefficientLinear,
                                   double temperatureChange) {
    if (originalLength <= 0.0) {
        throw std::invalid_argument("Original length must be positive");
    }
    return originalLength * (1.0 + coefficientLinear * temperatureChange);
}

/**
 * @brief Calculate coefficient of linear expansion from measurements
 *
 * α = ΔL/(L₀ × ΔT)
 *
 * @param lengthChange Measured change in length (m)
 * @param originalLength Original length (m)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Coefficient of linear expansion (K⁻¹)
 * @throws std::invalid_argument if original length or temperature change is zero
 */
inline double calculateLinearExpansionCoefficient(double lengthChange, double originalLength,
                                                  double temperatureChange) {
    if (originalLength <= 0.0) {
        throw std::invalid_argument("Original length must be positive");
    }
    if (std::abs(temperatureChange) < 1e-10) {
        throw std::invalid_argument("Temperature change must be non-zero");
    }
    return lengthChange / (originalLength * temperatureChange);
}

/**
 * @brief Calculate change in area due to thermal expansion
 *
 * ΔA ≈ 2α × A₀ × ΔT (for small expansions)
 *
 * For isotropic materials, area expansion coefficient β_A ≈ 2α
 *
 * @param originalArea Original area (m²)
 * @param coefficientLinear Coefficient of linear expansion (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Change in area (m²)
 * @throws std::invalid_argument if original area is non-positive
 */
inline double calculateAreaChange(double originalArea, double coefficientLinear,
                                  double temperatureChange) {
    if (originalArea <= 0.0) {
        throw std::invalid_argument("Original area must be positive");
    }
    return 2.0 * coefficientLinear * originalArea * temperatureChange;
}

/**
 * @brief Calculate final area after thermal expansion
 *
 * A = A₀(1 + 2αΔT)
 *
 * @param originalArea Original area (m²)
 * @param coefficientLinear Coefficient of linear expansion (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Final area (m²)
 * @throws std::invalid_argument if original area is non-positive
 */
inline double calculateFinalArea(double originalArea, double coefficientLinear,
                                 double temperatureChange) {
    if (originalArea <= 0.0) {
        throw std::invalid_argument("Original area must be positive");
    }
    return originalArea * (1.0 + 2.0 * coefficientLinear * temperatureChange);
}

// ============================================================================
// VOLUME EXPANSION
// ============================================================================

/**
 * @brief Calculate change in volume due to thermal expansion
 *
 * ΔV = β × V₀ × ΔT
 *
 * For isotropic solids: β ≈ 3α
 * For liquids: β is measured directly
 *
 * @param originalVolume Original volume (m³)
 * @param coefficientVolume Coefficient of volume expansion (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Change in volume (m³)
 * @throws std::invalid_argument if original volume is non-positive
 */
inline double calculateVolumeChange(double originalVolume, double coefficientVolume,
                                    double temperatureChange) {
    if (originalVolume <= 0.0) {
        throw std::invalid_argument("Original volume must be positive");
    }
    return coefficientVolume * originalVolume * temperatureChange;
}

/**
 * @brief Calculate final volume after thermal expansion
 *
 * V = V₀(1 + βΔT)
 *
 * @param originalVolume Original volume (m³)
 * @param coefficientVolume Coefficient of volume expansion (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Final volume (m³)
 * @throws std::invalid_argument if original volume is non-positive
 */
inline double calculateFinalVolume(double originalVolume, double coefficientVolume,
                                   double temperatureChange) {
    if (originalVolume <= 0.0) {
        throw std::invalid_argument("Original volume must be positive");
    }
    return originalVolume * (1.0 + coefficientVolume * temperatureChange);
}

/**
 * @brief Calculate coefficient of volume expansion from measurements
 *
 * β = ΔV/(V₀ × ΔT)
 *
 * @param volumeChange Measured change in volume (m³)
 * @param originalVolume Original volume (m³)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Coefficient of volume expansion (K⁻¹)
 * @throws std::invalid_argument if original volume or temperature change is zero
 */
inline double calculateVolumeExpansionCoefficient(double volumeChange, double originalVolume,
                                                  double temperatureChange) {
    if (originalVolume <= 0.0) {
        throw std::invalid_argument("Original volume must be positive");
    }
    if (std::abs(temperatureChange) < 1e-10) {
        throw std::invalid_argument("Temperature change must be non-zero");
    }
    return volumeChange / (originalVolume * temperatureChange);
}

/**
 * @brief Calculate volume expansion coefficient from linear expansion coefficient
 *
 * β ≈ 3α (for isotropic materials)
 *
 * @param coefficientLinear Coefficient of linear expansion (K⁻¹)
 * @return Coefficient of volume expansion (K⁻¹)
 */
inline double volumeFromLinearCoefficient(double coefficientLinear) {
    return 3.0 * coefficientLinear;
}

/**
 * @brief Calculate linear expansion coefficient from volume expansion coefficient
 *
 * α ≈ β/3 (for isotropic materials)
 *
 * @param coefficientVolume Coefficient of volume expansion (K⁻¹)
 * @return Coefficient of linear expansion (K⁻¹)
 */
inline double linearFromVolumeCoefficient(double coefficientVolume) {
    return coefficientVolume / 3.0;
}

// ============================================================================
// DENSITY CHANGE WITH TEMPERATURE
// ============================================================================

/**
 * @brief Calculate density after thermal expansion
 *
 * ρ = ρ₀/(1 + βΔT)
 *
 * Density decreases as volume expands
 *
 * @param originalDensity Original density (kg/m³)
 * @param coefficientVolume Coefficient of volume expansion (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Final density (kg/m³)
 * @throws std::invalid_argument if original density is non-positive
 */
inline double calculateDensityAfterExpansion(double originalDensity, double coefficientVolume,
                                             double temperatureChange) {
    if (originalDensity <= 0.0) {
        throw std::invalid_argument("Original density must be positive");
    }
    return originalDensity / (1.0 + coefficientVolume * temperatureChange);
}

/**
 * @brief Calculate relative density change
 *
 * Δρ/ρ₀ ≈ -βΔT (for small expansions)
 *
 * @param coefficientVolume Coefficient of volume expansion (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Relative density change (dimensionless)
 */
inline double relativeDensityChange(double coefficientVolume, double temperatureChange) {
    return -coefficientVolume * temperatureChange;
}

// ============================================================================
// PRACTICAL APPLICATIONS
// ============================================================================

/**
 * @brief Calculate thermal stress in constrained rod
 *
 * σ = E × α × ΔT
 *
 * When a rod is prevented from expanding, thermal stress develops
 *
 * @param youngsModulus Young's modulus of material (Pa)
 * @param coefficientLinear Coefficient of linear expansion (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Thermal stress (Pa)
 * @throws std::invalid_argument if Young's modulus is non-positive
 */
inline double calculateThermalStress(double youngsModulus, double coefficientLinear,
                                     double temperatureChange) {
    if (youngsModulus <= 0.0) {
        throw std::invalid_argument("Young's modulus must be positive");
    }
    return youngsModulus * coefficientLinear * temperatureChange;
}

/**
 * @brief Calculate required gap for thermal expansion joint
 *
 * Gap = α × L × ΔT_max
 *
 * @param length Length of structure (m)
 * @param coefficientLinear Coefficient of linear expansion (K⁻¹)
 * @param maxTemperatureChange Maximum expected temperature change (K or °C)
 * @return Required gap width (m)
 * @throws std::invalid_argument if length is non-positive
 */
inline double calculateExpansionGap(double length, double coefficientLinear,
                                    double maxTemperatureChange) {
    if (length <= 0.0) {
        throw std::invalid_argument("Length must be positive");
    }
    return coefficientLinear * length * maxTemperatureChange;
}

/**
 * @brief Calculate temperature change needed for specific length change
 *
 * ΔT = ΔL/(α × L₀)
 *
 * @param desiredLengthChange Desired change in length (m)
 * @param originalLength Original length (m)
 * @param coefficientLinear Coefficient of linear expansion (K⁻¹)
 * @return Required temperature change (K or °C)
 * @throws std::invalid_argument if original length or coefficient is non-positive
 */
inline double temperatureChangeForLengthChange(double desiredLengthChange, double originalLength,
                                               double coefficientLinear) {
    if (originalLength <= 0.0) {
        throw std::invalid_argument("Original length must be positive");
    }
    if (std::abs(coefficientLinear) < 1e-15) {
        throw std::invalid_argument("Coefficient of linear expansion must be non-zero");
    }
    return desiredLengthChange / (coefficientLinear * originalLength);
}

/**
 * @brief Calculate overflow from heated liquid in container
 *
 * V_overflow = V₀ × (β_liquid - β_container) × ΔT
 *
 * Apparent expansion = true expansion of liquid - expansion of container
 *
 * @param liquidVolume Initial volume of liquid (m³)
 * @param liquidExpansion Coefficient of volume expansion of liquid (K⁻¹)
 * @param containerExpansion Coefficient of volume expansion of container (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Volume of liquid overflow (m³)
 * @throws std::invalid_argument if liquid volume is non-positive
 */
inline double calculateLiquidOverflow(double liquidVolume, double liquidExpansion,
                                      double containerExpansion, double temperatureChange) {
    if (liquidVolume <= 0.0) {
        throw std::invalid_argument("Liquid volume must be positive");
    }
    return liquidVolume * (liquidExpansion - containerExpansion) * temperatureChange;
}

/**
 * @brief Calculate apparent coefficient of expansion
 *
 * β_apparent = β_liquid - β_container
 *
 * When liquid and container both expand
 *
 * @param liquidExpansion True coefficient of volume expansion of liquid (K⁻¹)
 * @param containerExpansion Coefficient of volume expansion of container (K⁻¹)
 * @return Apparent coefficient of expansion (K⁻¹)
 */
inline double apparentExpansionCoefficient(double liquidExpansion, double containerExpansion) {
    return liquidExpansion - containerExpansion;
}

/**
 * @brief Calculate bimetallic strip curvature
 *
 * For bimetallic strip (two metals with different expansion coefficients):
 * Radius of curvature r ≈ t/(3(α₁ - α₂)ΔT)
 *
 * where t is total thickness
 *
 * @param thickness Total thickness of bimetallic strip (m)
 * @param alpha1 Coefficient of first metal (K⁻¹)
 * @param alpha2 Coefficient of second metal (K⁻¹)
 * @param temperatureChange Change in temperature (K or °C)
 * @return Approximate radius of curvature (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double bimetallicStripCurvature(double thickness, double alpha1, double alpha2,
                                       double temperatureChange) {
    if (thickness <= 0.0) {
        throw std::invalid_argument("Thickness must be positive");
    }
    double alphaDiff = alpha1 - alpha2;
    if (std::abs(alphaDiff) < 1e-15) {
        throw std::invalid_argument("Expansion coefficients must be different");
    }
    if (std::abs(temperatureChange) < 1e-10) {
        throw std::invalid_argument("Temperature change must be non-zero");
    }

    return thickness / (3.0 * alphaDiff * temperatureChange);
}

} // namespace thermal_expansion
} // namespace physics

#endif // PHYSICS_THERMAL_EXPANSION_HPP
