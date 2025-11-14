#ifndef PHYSICS_ELASTICITY_HPP
#define PHYSICS_ELASTICITY_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file elasticity.hpp
 * @brief Comprehensive implementation of elasticity, stress-strain relationships,
 *        and deformation of solid materials
 *
 * This module implements:
 * - Young's Modulus (Elastic/Tensile Modulus)
 * - Bulk Modulus (Volume Elasticity)
 * - Shear Modulus (Rigidity Modulus)
 * - Beam bending and deflection
 * - Poisson's Ratio
 * - Elastic energy storage
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace elasticity {

/**
 * @namespace constants
 * @brief Physical constants for elasticity calculations
 */
namespace constants {
    // Material properties (Young's Modulus in Pa)
    constexpr double STEEL_YOUNGS_MODULUS = 200e9;      // 200 GPa
    constexpr double ALUMINUM_YOUNGS_MODULUS = 70e9;    // 70 GPa
    constexpr double COPPER_YOUNGS_MODULUS = 130e9;     // 130 GPa
    constexpr double GLASS_YOUNGS_MODULUS = 70e9;       // 70 GPa
    constexpr double CONCRETE_YOUNGS_MODULUS = 30e9;    // 30 GPa
    constexpr double WOOD_YOUNGS_MODULUS = 11e9;        // 11 GPa (varies)
    constexpr double RUBBER_YOUNGS_MODULUS = 0.01e9;    // 0.01 GPa

    // Bulk Modulus (in Pa)
    constexpr double WATER_BULK_MODULUS = 2.2e9;        // 2.2 GPa
    constexpr double STEEL_BULK_MODULUS = 160e9;        // 160 GPa

    // Shear Modulus (in Pa)
    constexpr double STEEL_SHEAR_MODULUS = 80e9;        // 80 GPa
    constexpr double ALUMINUM_SHEAR_MODULUS = 26e9;     // 26 GPa

    // Poisson's ratios (dimensionless)
    constexpr double STEEL_POISSON_RATIO = 0.30;
    constexpr double ALUMINUM_POISSON_RATIO = 0.33;
    constexpr double RUBBER_POISSON_RATIO = 0.49;       // Nearly incompressible
}

// ============================================================================
// STRESS AND STRAIN
// ============================================================================

/**
 * @brief Calculate stress (force per unit area)
 *
 * Stress σ = F/A
 *
 * @param force Applied force (N)
 * @param area Cross-sectional area (m²)
 * @return Stress (Pa = N/m²)
 * @throws std::invalid_argument if area is non-positive
 */
inline double calculateStress(double force, double area) {
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    return force / area;
}

/**
 * @brief Calculate strain (relative deformation)
 *
 * Strain ε = ΔL/L₀
 *
 * @param deformation Change in length (m)
 * @param originalLength Original length (m)
 * @return Strain (dimensionless)
 * @throws std::invalid_argument if original length is non-positive
 */
inline double calculateStrain(double deformation, double originalLength) {
    if (originalLength <= 0.0) {
        throw std::invalid_argument("Original length must be positive");
    }
    return deformation / originalLength;
}

// ============================================================================
// YOUNG'S MODULUS (ELASTIC/TENSILE MODULUS)
// ============================================================================

/**
 * @brief Calculate Young's Modulus from stress-strain data
 *
 * Young's Modulus E = σ/ε = (F/A)/(ΔL/L₀)
 *
 * Measures the stiffness of a material under tensile or compressive stress.
 *
 * @param stress Tensile/compressive stress (Pa)
 * @param strain Corresponding strain (dimensionless)
 * @return Young's Modulus (Pa)
 * @throws std::invalid_argument if strain is zero
 */
inline double calculateYoungsModulus(double stress, double strain) {
    if (std::abs(strain) < 1e-15) {
        throw std::invalid_argument("Strain must be non-zero");
    }
    return stress / strain;
}

/**
 * @brief Calculate elongation from Young's Modulus
 *
 * ΔL = (F × L₀)/(A × E)
 *
 * @param force Applied force (N)
 * @param originalLength Original length (m)
 * @param area Cross-sectional area (m²)
 * @param youngsModulus Young's Modulus of material (Pa)
 * @return Elongation/compression (m)
 * @throws std::invalid_argument if area or Young's Modulus is non-positive
 */
inline double calculateElongation(double force, double originalLength,
                                  double area, double youngsModulus) {
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    if (youngsModulus <= 0.0) {
        throw std::invalid_argument("Young's Modulus must be positive");
    }
    return (force * originalLength) / (area * youngsModulus);
}

/**
 * @brief Calculate stress from Young's Modulus and strain
 *
 * σ = E × ε (Hooke's Law for linear elastic materials)
 *
 * @param youngsModulus Young's Modulus (Pa)
 * @param strain Strain (dimensionless)
 * @return Stress (Pa)
 */
inline double calculateStressFromStrain(double youngsModulus, double strain) {
    return youngsModulus * strain;
}

/**
 * @brief Calculate strain from stress using Young's Modulus
 *
 * ε = σ/E
 *
 * @param stress Stress (Pa)
 * @param youngsModulus Young's Modulus (Pa)
 * @return Strain (dimensionless)
 * @throws std::invalid_argument if Young's Modulus is non-positive
 */
inline double calculateStrainFromStress(double stress, double youngsModulus) {
    if (youngsModulus <= 0.0) {
        throw std::invalid_argument("Young's Modulus must be positive");
    }
    return stress / youngsModulus;
}

// ============================================================================
// BULK MODULUS (VOLUME ELASTICITY)
// ============================================================================

/**
 * @brief Calculate Bulk Modulus (Volume Elasticity)
 *
 * K = -V₀(ΔP/ΔV) = -ΔP/(ΔV/V₀)
 *
 * Measures resistance to uniform compression. Negative sign ensures K is positive
 * (volume decreases when pressure increases).
 *
 * @param pressureChange Change in pressure (Pa)
 * @param volumeChange Change in volume (m³)
 * @param originalVolume Original volume (m³)
 * @return Bulk Modulus (Pa)
 * @throws std::invalid_argument if original volume is non-positive
 */
inline double calculateBulkModulus(double pressureChange, double volumeChange,
                                   double originalVolume) {
    if (originalVolume <= 0.0) {
        throw std::invalid_argument("Original volume must be positive");
    }
    if (std::abs(volumeChange) < 1e-15) {
        throw std::invalid_argument("Volume change must be non-zero");
    }
    return -pressureChange * originalVolume / volumeChange;
}

/**
 * @brief Calculate volume change from bulk modulus
 *
 * ΔV = -V₀ × ΔP/K
 *
 * @param originalVolume Original volume (m³)
 * @param pressureChange Change in pressure (Pa)
 * @param bulkModulus Bulk Modulus (Pa)
 * @return Volume change (m³), negative means compression
 * @throws std::invalid_argument if bulk modulus is non-positive
 */
inline double calculateVolumeChange(double originalVolume, double pressureChange,
                                    double bulkModulus) {
    if (bulkModulus <= 0.0) {
        throw std::invalid_argument("Bulk Modulus must be positive");
    }
    return -(originalVolume * pressureChange) / bulkModulus;
}

/**
 * @brief Calculate compressibility (reciprocal of bulk modulus)
 *
 * β = 1/K
 *
 * @param bulkModulus Bulk Modulus (Pa)
 * @return Compressibility (Pa⁻¹)
 * @throws std::invalid_argument if bulk modulus is non-positive
 */
inline double calculateCompressibility(double bulkModulus) {
    if (bulkModulus <= 0.0) {
        throw std::invalid_argument("Bulk Modulus must be positive");
    }
    return 1.0 / bulkModulus;
}

/**
 * @brief Calculate volumetric strain
 *
 * ε_v = ΔV/V₀
 *
 * @param volumeChange Change in volume (m³)
 * @param originalVolume Original volume (m³)
 * @return Volumetric strain (dimensionless)
 * @throws std::invalid_argument if original volume is non-positive
 */
inline double calculateVolumetricStrain(double volumeChange, double originalVolume) {
    if (originalVolume <= 0.0) {
        throw std::invalid_argument("Original volume must be positive");
    }
    return volumeChange / originalVolume;
}

// ============================================================================
// SHEAR MODULUS (RIGIDITY MODULUS)
// ============================================================================

/**
 * @brief Calculate Shear Modulus (Rigidity Modulus)
 *
 * G = τ/γ = (F/A)/(Δx/h)
 *
 * where τ is shear stress and γ is shear strain
 *
 * @param shearStress Shear stress (Pa)
 * @param shearStrain Shear strain (radians)
 * @return Shear Modulus (Pa)
 * @throws std::invalid_argument if shear strain is zero
 */
inline double calculateShearModulus(double shearStress, double shearStrain) {
    if (std::abs(shearStrain) < 1e-15) {
        throw std::invalid_argument("Shear strain must be non-zero");
    }
    return shearStress / shearStrain;
}

/**
 * @brief Calculate shear stress
 *
 * τ = F/A (force parallel to surface divided by area)
 *
 * @param force Tangential force (N)
 * @param area Area parallel to force (m²)
 * @return Shear stress (Pa)
 * @throws std::invalid_argument if area is non-positive
 */
inline double calculateShearStress(double force, double area) {
    if (area <= 0.0) {
        throw std::invalid_argument("Area must be positive");
    }
    return force / area;
}

/**
 * @brief Calculate shear strain from displacement
 *
 * γ = Δx/h (tangent of shear angle for small angles)
 *
 * @param displacement Lateral displacement (m)
 * @param height Original height (m)
 * @return Shear strain (radians)
 * @throws std::invalid_argument if height is non-positive
 */
inline double calculateShearStrain(double displacement, double height) {
    if (height <= 0.0) {
        throw std::invalid_argument("Height must be positive");
    }
    return displacement / height;
}

// ============================================================================
// POISSON'S RATIO
// ============================================================================

/**
 * @brief Calculate Poisson's Ratio
 *
 * ν = -ε_transverse/ε_axial
 *
 * Ratio of transverse strain to axial strain. For most materials, 0 < ν < 0.5
 *
 * @param transverseStrain Strain perpendicular to applied force (dimensionless)
 * @param axialStrain Strain parallel to applied force (dimensionless)
 * @return Poisson's Ratio (dimensionless)
 * @throws std::invalid_argument if axial strain is zero
 */
inline double calculatePoissonsRatio(double transverseStrain, double axialStrain) {
    if (std::abs(axialStrain) < 1e-15) {
        throw std::invalid_argument("Axial strain must be non-zero");
    }
    return -transverseStrain / axialStrain;
}

/**
 * @brief Calculate transverse strain from Poisson's Ratio
 *
 * ε_transverse = -ν × ε_axial
 *
 * @param poissonsRatio Poisson's Ratio (dimensionless)
 * @param axialStrain Axial strain (dimensionless)
 * @return Transverse strain (dimensionless)
 */
inline double calculateTransverseStrain(double poissonsRatio, double axialStrain) {
    return -poissonsRatio * axialStrain;
}

/**
 * @brief Relate elastic moduli using Poisson's Ratio
 *
 * G = E/(2(1 + ν))
 *
 * @param youngsModulus Young's Modulus (Pa)
 * @param poissonsRatio Poisson's Ratio (dimensionless)
 * @return Shear Modulus (Pa)
 * @throws std::invalid_argument if (1 + ν) is non-positive
 */
inline double calculateShearModulusFromYoungs(double youngsModulus,
                                              double poissonsRatio) {
    if ((1.0 + poissonsRatio) <= 0.0) {
        throw std::invalid_argument("Invalid Poisson's ratio");
    }
    return youngsModulus / (2.0 * (1.0 + poissonsRatio));
}

/**
 * @brief Calculate bulk modulus from Young's modulus and Poisson's ratio
 *
 * K = E/(3(1 - 2ν))
 *
 * @param youngsModulus Young's Modulus (Pa)
 * @param poissonsRatio Poisson's Ratio (dimensionless)
 * @return Bulk Modulus (Pa)
 * @throws std::invalid_argument if denominator is non-positive
 */
inline double calculateBulkModulusFromYoungs(double youngsModulus,
                                             double poissonsRatio) {
    double denominator = 3.0 * (1.0 - 2.0 * poissonsRatio);
    if (std::abs(denominator) < 1e-15) {
        throw std::invalid_argument("Invalid Poisson's ratio (≈ 0.5 for incompressible)");
    }
    return youngsModulus / denominator;
}

// ============================================================================
// BEAM BENDING AND DEFLECTION
// ============================================================================

/**
 * @brief Calculate maximum deflection of a simply supported beam with center load
 *
 * δ_max = (F × L³)/(48 × E × I)
 *
 * For a beam supported at both ends with load at center.
 *
 * @param force Load at center (N)
 * @param length Beam length (m)
 * @param youngsModulus Young's Modulus (Pa)
 * @param momentOfInertia Second moment of area (m⁴)
 * @return Maximum deflection at center (m)
 * @throws std::invalid_argument if E or I is non-positive
 */
inline double calculateBeamDeflectionCenterLoad(double force, double length,
                                                double youngsModulus,
                                                double momentOfInertia) {
    if (youngsModulus <= 0.0 || momentOfInertia <= 0.0) {
        throw std::invalid_argument("Young's Modulus and Moment of Inertia must be positive");
    }
    return (force * std::pow(length, 3)) / (48.0 * youngsModulus * momentOfInertia);
}

/**
 * @brief Calculate maximum deflection of a cantilever beam with end load
 *
 * δ_max = (F × L³)/(3 × E × I)
 *
 * For a beam fixed at one end with load at free end.
 *
 * @param force Load at free end (N)
 * @param length Beam length (m)
 * @param youngsModulus Young's Modulus (Pa)
 * @param momentOfInertia Second moment of area (m⁴)
 * @return Maximum deflection at free end (m)
 * @throws std::invalid_argument if E or I is non-positive
 */
inline double calculateCantileverDeflection(double force, double length,
                                            double youngsModulus,
                                            double momentOfInertia) {
    if (youngsModulus <= 0.0 || momentOfInertia <= 0.0) {
        throw std::invalid_argument("Young's Modulus and Moment of Inertia must be positive");
    }
    return (force * std::pow(length, 3)) / (3.0 * youngsModulus * momentOfInertia);
}

/**
 * @brief Calculate deflection of simply supported beam with uniform load
 *
 * δ_max = (5 × w × L⁴)/(384 × E × I)
 *
 * @param loadPerLength Uniform load per unit length (N/m)
 * @param length Beam length (m)
 * @param youngsModulus Young's Modulus (Pa)
 * @param momentOfInertia Second moment of area (m⁴)
 * @return Maximum deflection at center (m)
 * @throws std::invalid_argument if E or I is non-positive
 */
inline double calculateBeamDeflectionUniformLoad(double loadPerLength, double length,
                                                 double youngsModulus,
                                                 double momentOfInertia) {
    if (youngsModulus <= 0.0 || momentOfInertia <= 0.0) {
        throw std::invalid_argument("Young's Modulus and Moment of Inertia must be positive");
    }
    return (5.0 * loadPerLength * std::pow(length, 4)) /
           (384.0 * youngsModulus * momentOfInertia);
}

/**
 * @brief Calculate second moment of area for rectangular cross-section
 *
 * I = (b × h³)/12
 *
 * For a rectangle with width b and height h (about neutral axis).
 *
 * @param width Width of rectangle (m)
 * @param height Height of rectangle (m)
 * @return Second moment of area (m⁴)
 * @throws std::invalid_argument if width or height is non-positive
 */
inline double calculateRectangularMomentOfInertia(double width, double height) {
    if (width <= 0.0 || height <= 0.0) {
        throw std::invalid_argument("Width and height must be positive");
    }
    return (width * std::pow(height, 3)) / 12.0;
}

/**
 * @brief Calculate second moment of area for circular cross-section
 *
 * I = (π × d⁴)/64 = (π × r⁴)/4
 *
 * @param diameter Diameter of circle (m)
 * @return Second moment of area (m⁴)
 * @throws std::invalid_argument if diameter is non-positive
 */
inline double calculateCircularMomentOfInertia(double diameter) {
    if (diameter <= 0.0) {
        throw std::invalid_argument("Diameter must be positive");
    }
    return (M_PI * std::pow(diameter, 4)) / 64.0;
}

/**
 * @brief Calculate bending stress in a beam
 *
 * σ = (M × y)/I
 *
 * where M is bending moment, y is distance from neutral axis, I is moment of inertia
 *
 * @param bendingMoment Bending moment (N·m)
 * @param distanceFromNeutralAxis Distance from neutral axis (m)
 * @param momentOfInertia Second moment of area (m⁴)
 * @return Bending stress (Pa)
 * @throws std::invalid_argument if moment of inertia is non-positive
 */
inline double calculateBendingStress(double bendingMoment,
                                     double distanceFromNeutralAxis,
                                     double momentOfInertia) {
    if (momentOfInertia <= 0.0) {
        throw std::invalid_argument("Moment of Inertia must be positive");
    }
    return (bendingMoment * distanceFromNeutralAxis) / momentOfInertia;
}

// ============================================================================
// ELASTIC ENERGY
// ============================================================================

/**
 * @brief Calculate elastic potential energy stored in stretched/compressed material
 *
 * U = (1/2) × k × (ΔL)² = (1/2) × (EA/L) × (ΔL)²
 *
 * where k = EA/L is the effective spring constant
 *
 * @param youngsModulus Young's Modulus (Pa)
 * @param area Cross-sectional area (m²)
 * @param originalLength Original length (m)
 * @param deformation Elongation/compression (m)
 * @return Elastic potential energy (J)
 * @throws std::invalid_argument if parameters are non-positive
 */
inline double calculateElasticEnergy(double youngsModulus, double area,
                                     double originalLength, double deformation) {
    if (youngsModulus <= 0.0 || area <= 0.0 || originalLength <= 0.0) {
        throw std::invalid_argument("Material parameters must be positive");
    }
    double springConstant = (youngsModulus * area) / originalLength;
    return 0.5 * springConstant * deformation * deformation;
}

/**
 * @brief Calculate energy density (energy per unit volume) in elastic material
 *
 * u = (1/2) × E × ε² = (1/2) × σ × ε = σ²/(2E)
 *
 * @param stress Stress (Pa)
 * @param youngsModulus Young's Modulus (Pa)
 * @return Energy density (J/m³)
 * @throws std::invalid_argument if Young's Modulus is non-positive
 */
inline double calculateEnergyDensity(double stress, double youngsModulus) {
    if (youngsModulus <= 0.0) {
        throw std::invalid_argument("Young's Modulus must be positive");
    }
    return (stress * stress) / (2.0 * youngsModulus);
}

/**
 * @brief Calculate elastic energy from stress and strain
 *
 * U = (1/2) × σ × ε × V
 *
 * @param stress Stress (Pa)
 * @param strain Strain (dimensionless)
 * @param volume Volume of material (m³)
 * @return Elastic energy (J)
 * @throws std::invalid_argument if volume is non-positive
 */
inline double calculateElasticEnergyFromStressStrain(double stress, double strain,
                                                     double volume) {
    if (volume <= 0.0) {
        throw std::invalid_argument("Volume must be positive");
    }
    return 0.5 * stress * strain * volume;
}

} // namespace elasticity
} // namespace physics

#endif // PHYSICS_ELASTICITY_HPP
