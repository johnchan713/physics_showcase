#ifndef PHYSICS_SURFACE_TENSION_HPP
#define PHYSICS_SURFACE_TENSION_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file surface_tension.hpp
 * @brief Comprehensive implementation of surface tension phenomena
 *
 * This module implements:
 * - Surface tension and surface energy
 * - Pressure due to curved surfaces (Young-Laplace equation)
 * - Capillary rise and depression
 * - Contact angle effects
 * - Meniscus formation
 * - Soap bubbles and droplets
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace surface_tension {

/**
 * @namespace constants
 * @brief Physical constants for surface tension calculations
 */
namespace constants {
    // Surface tension coefficients (N/m at 20°C)
    constexpr double WATER_SURFACE_TENSION = 0.0728;          // Water-air
    constexpr double MERCURY_SURFACE_TENSION = 0.4865;        // Mercury-air
    constexpr double ETHANOL_SURFACE_TENSION = 0.0223;        // Ethanol-air
    constexpr double GLYCERIN_SURFACE_TENSION = 0.0631;       // Glycerin-air
    constexpr double SOAP_SOLUTION_SURFACE_TENSION = 0.025;   // Soap solution-air
    constexpr double OLIVE_OIL_SURFACE_TENSION = 0.032;       // Olive oil-air

    // Densities (kg/m³)
    constexpr double WATER_DENSITY = 1000.0;
    constexpr double MERCURY_DENSITY = 13600.0;
    constexpr double ETHANOL_DENSITY = 789.0;

    // Standard gravity
    constexpr double GRAVITY = 9.81;  // m/s²

    // Contact angles (radians)
    constexpr double WATER_GLASS_CONTACT_ANGLE = 0.0;         // Perfect wetting
    constexpr double MERCURY_GLASS_CONTACT_ANGLE = 2.4435;    // ~140° (non-wetting)
}

// ============================================================================
// SURFACE TENSION FUNDAMENTALS
// ============================================================================

/**
 * @brief Calculate force due to surface tension
 *
 * F = γ × L
 *
 * where γ is surface tension and L is length of contact line
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param length Length of the contact line (m)
 * @return Force due to surface tension (N)
 * @throws std::invalid_argument if length is negative
 */
inline double calculateSurfaceTensionForce(double surfaceTension, double length) {
    if (length < 0.0) {
        throw std::invalid_argument("Length cannot be negative");
    }
    return surfaceTension * length;
}

/**
 * @brief Calculate surface energy
 *
 * E = γ × A
 *
 * Energy required to create a surface of area A
 *
 * @param surfaceTension Surface tension coefficient (N/m = J/m²)
 * @param area Surface area (m²)
 * @return Surface energy (J)
 * @throws std::invalid_argument if area is negative
 */
inline double calculateSurfaceEnergy(double surfaceTension, double area) {
    if (area < 0.0) {
        throw std::invalid_argument("Area cannot be negative");
    }
    return surfaceTension * area;
}

/**
 * @brief Calculate work done in stretching a liquid film
 *
 * W = γ × ΔA
 *
 * For a film with two surfaces (like soap film), W = 2γ × ΔA
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param areaChange Change in area (m²)
 * @param numberOfSurfaces Number of surfaces (1 for droplet, 2 for film)
 * @return Work done (J)
 * @throws std::invalid_argument if number of surfaces is non-positive
 */
inline double calculateWorkInStretching(double surfaceTension, double areaChange,
                                        int numberOfSurfaces = 1) {
    if (numberOfSurfaces <= 0) {
        throw std::invalid_argument("Number of surfaces must be positive");
    }
    return numberOfSurfaces * surfaceTension * areaChange;
}

// ============================================================================
// PRESSURE DUE TO CURVED SURFACES (YOUNG-LAPLACE EQUATION)
// ============================================================================

/**
 * @brief Calculate excess pressure inside a spherical droplet
 *
 * ΔP = 2γ/r
 *
 * Pressure inside exceeds outside pressure due to surface tension
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param radius Radius of droplet (m)
 * @return Excess pressure (Pa)
 * @throws std::invalid_argument if radius is non-positive
 */
inline double calculateDropletPressure(double surfaceTension, double radius) {
    if (radius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return (2.0 * surfaceTension) / radius;
}

/**
 * @brief Calculate excess pressure inside a soap bubble
 *
 * ΔP = 4γ/r
 *
 * Soap bubble has two surfaces (inside and outside), so pressure is doubled
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param radius Radius of bubble (m)
 * @return Excess pressure (Pa)
 * @throws std::invalid_argument if radius is non-positive
 */
inline double calculateBubblePressure(double surfaceTension, double radius) {
    if (radius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return (4.0 * surfaceTension) / radius;
}

/**
 * @brief Calculate excess pressure in cylindrical surface (like liquid jet)
 *
 * ΔP = γ/r
 *
 * For a cylinder, only one principal radius of curvature contributes
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param radius Radius of cylinder (m)
 * @return Excess pressure (Pa)
 * @throws std::invalid_argument if radius is non-positive
 */
inline double calculateCylindricalPressure(double surfaceTension, double radius) {
    if (radius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return surfaceTension / radius;
}

/**
 * @brief General Young-Laplace equation for arbitrary curved surface
 *
 * ΔP = γ(1/R₁ + 1/R₂)
 *
 * where R₁ and R₂ are principal radii of curvature
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param radius1 First principal radius of curvature (m)
 * @param radius2 Second principal radius of curvature (m)
 * @return Excess pressure (Pa)
 * @throws std::invalid_argument if any radius is non-positive
 */
inline double calculateYoungLaplacePressure(double surfaceTension,
                                            double radius1, double radius2) {
    if (radius1 <= 0.0 || radius2 <= 0.0) {
        throw std::invalid_argument("Radii must be positive");
    }
    return surfaceTension * (1.0 / radius1 + 1.0 / radius2);
}

// ============================================================================
// CAPILLARY RISE AND DEPRESSION
// ============================================================================

/**
 * @brief Calculate capillary rise (or depression) in a tube
 *
 * h = (2γ cos θ)/(ρgr)
 *
 * Positive h means rise (θ < 90°, wetting fluid like water in glass)
 * Negative h means depression (θ > 90°, non-wetting fluid like mercury)
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param contactAngle Contact angle between liquid and tube wall (radians)
 * @param density Liquid density (kg/m³)
 * @param tubeRadius Inner radius of capillary tube (m)
 * @param gravity Gravitational acceleration (m/s²)
 * @return Height of rise (positive) or depression (negative) (m)
 * @throws std::invalid_argument if density, radius, or gravity is non-positive
 */
inline double calculateCapillaryRise(double surfaceTension, double contactAngle,
                                     double density, double tubeRadius,
                                     double gravity = constants::GRAVITY) {
    if (density <= 0.0) {
        throw std::invalid_argument("Density must be positive");
    }
    if (tubeRadius <= 0.0) {
        throw std::invalid_argument("Tube radius must be positive");
    }
    if (gravity <= 0.0) {
        throw std::invalid_argument("Gravity must be positive");
    }

    return (2.0 * surfaceTension * std::cos(contactAngle)) /
           (density * gravity * tubeRadius);
}

/**
 * @brief Calculate capillary rise for water in glass tube
 *
 * Simplified formula assuming perfect wetting (θ = 0°, cos θ = 1)
 * h = 2γ/(ρgr)
 *
 * @param tubeRadius Inner radius of tube (m)
 * @param surfaceTension Surface tension (N/m), defaults to water
 * @param density Liquid density (kg/m³), defaults to water
 * @param gravity Gravitational acceleration (m/s²)
 * @return Height of capillary rise (m)
 */
inline double calculateWaterCapillaryRise(double tubeRadius,
                                          double surfaceTension = constants::WATER_SURFACE_TENSION,
                                          double density = constants::WATER_DENSITY,
                                          double gravity = constants::GRAVITY) {
    return calculateCapillaryRise(surfaceTension, 0.0, density, tubeRadius, gravity);
}

/**
 * @brief Calculate capillary depression for mercury in glass tube
 *
 * Mercury exhibits non-wetting behavior (θ ≈ 140°)
 *
 * @param tubeRadius Inner radius of tube (m)
 * @param surfaceTension Surface tension (N/m), defaults to mercury
 * @param density Liquid density (kg/m³), defaults to mercury
 * @param contactAngle Contact angle (radians), defaults to ~140°
 * @param gravity Gravitational acceleration (m/s²)
 * @return Height of capillary depression (negative value) (m)
 */
inline double calculateMercuryCapillaryDepression(
    double tubeRadius,
    double surfaceTension = constants::MERCURY_SURFACE_TENSION,
    double density = constants::MERCURY_DENSITY,
    double contactAngle = constants::MERCURY_GLASS_CONTACT_ANGLE,
    double gravity = constants::GRAVITY) {
    return calculateCapillaryRise(surfaceTension, contactAngle, density,
                                  tubeRadius, gravity);
}

/**
 * @brief Calculate tube radius from observed capillary rise
 *
 * r = (2γ cos θ)/(ρgh)
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param contactAngle Contact angle (radians)
 * @param density Liquid density (kg/m³)
 * @param height Observed height of rise (m)
 * @param gravity Gravitational acceleration (m/s²)
 * @return Tube radius (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double calculateTubeRadiusFromRise(double surfaceTension, double contactAngle,
                                          double density, double height,
                                          double gravity = constants::GRAVITY) {
    if (density <= 0.0 || gravity <= 0.0) {
        throw std::invalid_argument("Density and gravity must be positive");
    }
    if (std::abs(height) < 1e-15) {
        throw std::invalid_argument("Height must be non-zero");
    }

    return (2.0 * surfaceTension * std::cos(contactAngle)) /
           (density * gravity * height);
}

// ============================================================================
// CAPILLARY BETWEEN PARALLEL PLATES
// ============================================================================

/**
 * @brief Calculate capillary rise between two parallel plates
 *
 * h = (2γ cos θ)/(ρgd)
 *
 * where d is the separation between plates
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param contactAngle Contact angle (radians)
 * @param density Liquid density (kg/m³)
 * @param separation Distance between plates (m)
 * @param gravity Gravitational acceleration (m/s²)
 * @return Height of capillary rise (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double calculateCapillaryRiseBetweenPlates(double surfaceTension,
                                                  double contactAngle,
                                                  double density,
                                                  double separation,
                                                  double gravity = constants::GRAVITY) {
    if (density <= 0.0 || separation <= 0.0 || gravity <= 0.0) {
        throw std::invalid_argument("Density, separation, and gravity must be positive");
    }

    return (2.0 * surfaceTension * std::cos(contactAngle)) /
           (density * gravity * separation);
}

// ============================================================================
// MENISCUS AND CONTACT ANGLE
// ============================================================================

/**
 * @brief Calculate the shape parameter of a meniscus
 *
 * Capillary length a = √(γ/(ρg))
 *
 * Characteristic length scale for surface tension effects
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param density Liquid density (kg/m³)
 * @param gravity Gravitational acceleration (m/s²)
 * @return Capillary length (m)
 * @throws std::invalid_argument if density or gravity is non-positive
 */
inline double calculateCapillaryLength(double surfaceTension, double density,
                                       double gravity = constants::GRAVITY) {
    if (density <= 0.0 || gravity <= 0.0) {
        throw std::invalid_argument("Density and gravity must be positive");
    }
    return std::sqrt(surfaceTension / (density * gravity));
}

/**
 * @brief Calculate contact angle from force balance
 *
 * cos θ = (γ_SV - γ_SL)/γ_LV
 *
 * Young's equation relating interfacial tensions
 *
 * @param solidVaporTension Solid-vapor interfacial tension (N/m)
 * @param solidLiquidTension Solid-liquid interfacial tension (N/m)
 * @param liquidVaporTension Liquid-vapor surface tension (N/m)
 * @return Contact angle (radians)
 * @throws std::invalid_argument if liquid-vapor tension is non-positive
 */
inline double calculateContactAngle(double solidVaporTension,
                                    double solidLiquidTension,
                                    double liquidVaporTension) {
    if (liquidVaporTension <= 0.0) {
        throw std::invalid_argument("Liquid-vapor tension must be positive");
    }

    double cosTheta = (solidVaporTension - solidLiquidTension) / liquidVaporTension;

    // Clamp to valid range [-1, 1] for numerical stability
    if (cosTheta > 1.0) cosTheta = 1.0;
    if (cosTheta < -1.0) cosTheta = -1.0;

    return std::acos(cosTheta);
}

// ============================================================================
// DROPLET AND BUBBLE DYNAMICS
// ============================================================================

/**
 * @brief Calculate the radius of a droplet from its volume
 *
 * r = ∛(3V/(4π))
 *
 * @param volume Volume of droplet (m³)
 * @return Radius of droplet (m)
 * @throws std::invalid_argument if volume is non-positive
 */
inline double calculateDropletRadius(double volume) {
    if (volume <= 0.0) {
        throw std::invalid_argument("Volume must be positive");
    }
    return std::cbrt((3.0 * volume) / (4.0 * M_PI));
}

/**
 * @brief Calculate energy required to split a droplet into smaller droplets
 *
 * ΔE = γ × ΔA = γ × 4π(n × r_small² - r_large²)
 *
 * where one large droplet of radius r_large splits into n droplets of radius r_small
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param largeRadius Radius of original droplet (m)
 * @param smallRadius Radius of each smaller droplet (m)
 * @param numberOfDroplets Number of smaller droplets
 * @return Energy required (J)
 * @throws std::invalid_argument if radii are non-positive or count is less than 2
 */
inline double calculateDropletSplittingEnergy(double surfaceTension,
                                              double largeRadius,
                                              double smallRadius,
                                              int numberOfDroplets) {
    if (largeRadius <= 0.0 || smallRadius <= 0.0) {
        throw std::invalid_argument("Radii must be positive");
    }
    if (numberOfDroplets < 2) {
        throw std::invalid_argument("Must split into at least 2 droplets");
    }

    double originalArea = 4.0 * M_PI * largeRadius * largeRadius;
    double newTotalArea = numberOfDroplets * 4.0 * M_PI * smallRadius * smallRadius;
    double areaChange = newTotalArea - originalArea;

    return surfaceTension * areaChange;
}

/**
 * @brief Calculate radius of smaller droplets when one splits while conserving volume
 *
 * r_small = r_large / ∛n
 *
 * @param largeRadius Radius of original droplet (m)
 * @param numberOfDroplets Number of equal smaller droplets
 * @return Radius of each smaller droplet (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double calculateSplitDropletRadius(double largeRadius, int numberOfDroplets) {
    if (largeRadius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }
    if (numberOfDroplets < 2) {
        throw std::invalid_argument("Must split into at least 2 droplets");
    }

    return largeRadius / std::cbrt(static_cast<double>(numberOfDroplets));
}

/**
 * @brief Calculate excess pressure difference between two connected bubbles
 *
 * For two soap bubbles of different radii connected by a tube:
 * ΔP = 4γ(1/r₁ - 1/r₂)
 *
 * Air flows from smaller to larger bubble (higher to lower pressure)
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param radius1 Radius of first bubble (m)
 * @param radius2 Radius of second bubble (m)
 * @return Pressure difference (Pa)
 * @throws std::invalid_argument if radii are non-positive
 */
inline double calculateBubblePressureDifference(double surfaceTension,
                                                double radius1, double radius2) {
    if (radius1 <= 0.0 || radius2 <= 0.0) {
        throw std::invalid_argument("Radii must be positive");
    }

    return 4.0 * surfaceTension * (1.0 / radius1 - 1.0 / radius2);
}

// ============================================================================
// MISCELLANEOUS APPLICATIONS
// ============================================================================

/**
 * @brief Calculate maximum weight supported by a wire frame on liquid surface
 *
 * F_max = γ × L
 *
 * For a wire of length L resting on liquid surface
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param wireLength Length of wire in contact with surface (m)
 * @return Maximum supported force (N)
 * @throws std::invalid_argument if length is negative
 */
inline double calculateMaxSupportedWeight(double surfaceTension, double wireLength) {
    if (wireLength < 0.0) {
        throw std::invalid_argument("Wire length cannot be negative");
    }
    return surfaceTension * wireLength;
}

/**
 * @brief Calculate force needed to pull a ring from liquid surface
 *
 * F = 2γ × (2πR) = 4πγR
 *
 * Factor of 2 because ring has inner and outer surfaces
 *
 * @param surfaceTension Surface tension coefficient (N/m)
 * @param ringRadius Radius of ring (m)
 * @return Force required to detach ring (N)
 * @throws std::invalid_argument if radius is non-positive
 */
inline double calculateRingDetachmentForce(double surfaceTension, double ringRadius) {
    if (ringRadius <= 0.0) {
        throw std::invalid_argument("Ring radius must be positive");
    }
    return 4.0 * M_PI * surfaceTension * ringRadius;
}

/**
 * @brief Calculate surface tension from capillary rise measurement
 *
 * γ = (ρghr)/(2 cos θ)
 *
 * Common experimental method to determine surface tension
 *
 * @param density Liquid density (kg/m³)
 * @param height Observed capillary rise (m)
 * @param tubeRadius Tube radius (m)
 * @param contactAngle Contact angle (radians)
 * @param gravity Gravitational acceleration (m/s²)
 * @return Surface tension (N/m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double calculateSurfaceTensionFromRise(double density, double height,
                                              double tubeRadius, double contactAngle,
                                              double gravity = constants::GRAVITY) {
    if (density <= 0.0 || tubeRadius <= 0.0 || gravity <= 0.0) {
        throw std::invalid_argument("Density, radius, and gravity must be positive");
    }

    double cosTheta = std::cos(contactAngle);
    if (std::abs(cosTheta) < 1e-15) {
        throw std::invalid_argument("Contact angle must not be 90° (cos θ = 0)");
    }

    return (density * gravity * height * tubeRadius) / (2.0 * cosTheta);
}

} // namespace surface_tension
} // namespace physics

#endif // PHYSICS_SURFACE_TENSION_HPP
