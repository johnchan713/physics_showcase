#ifndef PHYSICS_OPTICS_HPP
#define PHYSICS_OPTICS_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file optics.hpp
 * @brief Comprehensive implementation of geometric and physical optics
 *
 * This module implements:
 * - Reflection and refraction (Snell's law)
 * - Oblique incidence (change in velocity and direction)
 * - Critical angle and total internal reflection
 * - Brewster's law (polarization)
 * - Thin lens equation and lens formula
 * - Principal focus and focal length
 * - Mirror equations
 * - Magnification
 * - Optical instruments (microscopes, telescopes)
 * - Lens combinations
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace optics {

/**
 * @namespace constants
 * @brief Refractive indices for common materials
 */
namespace constants {
    constexpr double REFRACTIVE_INDEX_VACUUM = 1.0;
    constexpr double REFRACTIVE_INDEX_AIR = 1.000293;
    constexpr double REFRACTIVE_INDEX_WATER = 1.333;
    constexpr double REFRACTIVE_INDEX_GLASS = 1.5;          // Typical crown glass
    constexpr double REFRACTIVE_INDEX_DIAMOND = 2.417;
    constexpr double REFRACTIVE_INDEX_ICE = 1.31;
    constexpr double REFRACTIVE_INDEX_ETHANOL = 1.361;
    constexpr double REFRACTIVE_INDEX_QUARTZ = 1.458;
}

// ============================================================================
// REFLECTION
// ============================================================================

/**
 * @brief Calculate angle of reflection
 *
 * θ_r = θ_i (Law of reflection)
 *
 * @param incidentAngle Angle of incidence (radians)
 * @return Angle of reflection (radians)
 */
inline double angleOfReflection(double incidentAngle) {
    return incidentAngle;
}

// ============================================================================
// REFRACTION AND SNELL'S LAW
// ============================================================================

/**
 * @brief Calculate angle of refraction using Snell's law
 *
 * n₁ sin(θ₁) = n₂ sin(θ₂)
 *
 * @param incidentAngle Incident angle from normal (radians)
 * @param n1 Refractive index of first medium
 * @param n2 Refractive index of second medium
 * @return Refracted angle from normal (radians)
 * @throws std::invalid_argument if refractive indices are non-positive or total internal reflection occurs
 */
inline double snellsLaw(double incidentAngle, double n1, double n2) {
    if (n1 <= 0.0 || n2 <= 0.0) {
        throw std::invalid_argument("Refractive indices must be positive");
    }

    double sinTheta2 = (n1 / n2) * std::sin(incidentAngle);

    // Check for total internal reflection
    if (sinTheta2 > 1.0) {
        throw std::invalid_argument("Total internal reflection occurs - no refracted ray");
    }

    return std::asin(sinTheta2);
}

/**
 * @brief Calculate refractive index from angles
 *
 * n = sin(θ₁) / sin(θ₂)
 *
 * Relative refractive index when going from air (n≈1) to medium
 *
 * @param incidentAngle Incident angle (radians)
 * @param refractedAngle Refracted angle (radians)
 * @return Relative refractive index
 * @throws std::invalid_argument if refracted angle sine is zero
 */
inline double refractiveIndexFromAngles(double incidentAngle, double refractedAngle) {
    double sinRefracted = std::sin(refractedAngle);
    if (std::abs(sinRefracted) < 1e-15) {
        throw std::invalid_argument("Refracted angle sine must be non-zero");
    }
    return std::sin(incidentAngle) / sinRefracted;
}

// ============================================================================
// OBLIQUE INCIDENCE: CHANGE IN VELOCITY AND DIRECTION
// ============================================================================

/**
 * @brief Calculate velocity of light in medium
 *
 * v = c / n
 *
 * @param speedInVacuum Speed of light in vacuum (m/s)
 * @param refractiveIndex Refractive index of medium
 * @return Speed in medium (m/s)
 * @throws std::invalid_argument if refractive index is non-positive
 */
inline double velocityInMedium(double speedInVacuum, double refractiveIndex) {
    if (refractiveIndex <= 0.0) {
        throw std::invalid_argument("Refractive index must be positive");
    }
    return speedInVacuum / refractiveIndex;
}

/**
 * @brief Calculate change in velocity at oblique incidence
 *
 * Δv = c(1/n₁ - 1/n₂)
 *
 * @param speedInVacuum Speed in vacuum (m/s)
 * @param n1 Initial medium refractive index
 * @param n2 Final medium refractive index
 * @return Change in velocity (m/s)
 * @throws std::invalid_argument if refractive indices are non-positive
 */
inline double velocityChange(double speedInVacuum, double n1, double n2) {
    if (n1 <= 0.0 || n2 <= 0.0) {
        throw std::invalid_argument("Refractive indices must be positive");
    }
    return speedInVacuum * (1.0/n1 - 1.0/n2);
}

/**
 * @brief Calculate deviation angle at interface
 *
 * δ = θ₁ - θ₂
 *
 * Positive deviation means ray bends toward normal
 *
 * @param incidentAngle Incident angle (radians)
 * @param refractedAngle Refracted angle (radians)
 * @return Deviation angle (radians)
 */
inline double deviationAngle(double incidentAngle, double refractedAngle) {
    return incidentAngle - refractedAngle;
}

// ============================================================================
// CRITICAL ANGLE AND TOTAL INTERNAL REFLECTION
// ============================================================================

/**
 * @brief Calculate critical angle for total internal reflection
 *
 * θ_c = arcsin(n₂/n₁)
 *
 * Valid only when n₁ > n₂ (going from denser to rarer medium)
 *
 * @param n1 Refractive index of denser medium
 * @param n2 Refractive index of rarer medium
 * @return Critical angle (radians)
 * @throws std::invalid_argument if n1 <= n2 or indices are non-positive
 */
inline double criticalAngle(double n1, double n2) {
    if (n1 <= 0.0 || n2 <= 0.0) {
        throw std::invalid_argument("Refractive indices must be positive");
    }
    if (n1 <= n2) {
        throw std::invalid_argument("Critical angle exists only when n1 > n2 (denser to rarer)");
    }
    return std::asin(n2 / n1);
}

/**
 * @brief Check if total internal reflection occurs
 *
 * TIR occurs when θ > θ_c and n₁ > n₂
 *
 * @param incidentAngle Incident angle (radians)
 * @param n1 Refractive index of first medium
 * @param n2 Refractive index of second medium
 * @return true if total internal reflection occurs
 */
inline bool isTotalInternalReflection(double incidentAngle, double n1, double n2) {
    if (n1 <= n2) return false;  // Can't have TIR going to denser medium

    double critAngle = criticalAngle(n1, n2);
    return incidentAngle > critAngle;
}

// ============================================================================
// BREWSTER'S LAW (POLARIZATION)
// ============================================================================

/**
 * @brief Calculate Brewster's angle (polarizing angle)
 *
 * tan(θ_B) = n₂/n₁
 * θ_B = arctan(n₂/n₁)
 *
 * At Brewster's angle, reflected light is completely polarized
 *
 * @param n1 Refractive index of first medium
 * @param n2 Refractive index of second medium
 * @return Brewster's angle (radians)
 * @throws std::invalid_argument if n1 is non-positive
 */
inline double brewstersAngle(double n1, double n2) {
    if (n1 <= 0.0) {
        throw std::invalid_argument("First medium refractive index must be positive");
    }
    return std::atan(n2 / n1);
}

/**
 * @brief Calculate refractive index from Brewster's angle
 *
 * n = tan(θ_B)
 *
 * For light going from air (n₁≈1) to medium
 *
 * @param brewsterAngle Brewster's angle (radians)
 * @return Refractive index of medium
 */
inline double refractiveIndexFromBrewster(double brewsterAngle) {
    return std::tan(brewsterAngle);
}

/**
 * @brief Verify Brewster's condition
 *
 * At Brewster's angle: θ_B + θ_r = π/2
 *
 * Reflected and refracted rays are perpendicular
 *
 * @param brewsterAngle Brewster's angle (radians)
 * @param refractedAngle Refracted angle (radians)
 * @return true if Brewster's condition is satisfied (within tolerance)
 */
inline bool verifyBrewsterCondition(double brewsterAngle, double refractedAngle) {
    constexpr double tolerance = 1e-6;
    return std::abs(brewsterAngle + refractedAngle - M_PI/2.0) < tolerance;
}

// ============================================================================
// THIN LENS EQUATION AND LENS FORMULA
// ============================================================================

/**
 * @brief Calculate focal length from object and image distances (lens formula)
 *
 * 1/f = 1/v - 1/u
 *
 * Sign convention: u is negative for real objects, v is positive for real images
 *
 * @param objectDistance Distance from lens to object (u) - negative for real object
 * @param imageDistance Distance from lens to image (v) - positive for real image
 * @return Focal length (m)
 * @throws std::invalid_argument if distances are zero
 */
inline double lensFormula(double objectDistance, double imageDistance) {
    if (std::abs(objectDistance) < 1e-15 || std::abs(imageDistance) < 1e-15) {
        throw std::invalid_argument("Object and image distances must be non-zero");
    }
    return 1.0 / (1.0/imageDistance - 1.0/objectDistance);
}

/**
 * @brief Calculate image distance from lens formula
 *
 * 1/v = 1/f + 1/u
 *
 * @param focalLength Focal length (m)
 * @param objectDistance Object distance (m)
 * @return Image distance (m)
 * @throws std::invalid_argument if calculation is invalid
 */
inline double imageDistance(double focalLength, double objectDistance) {
    if (std::abs(focalLength) < 1e-15) {
        throw std::invalid_argument("Focal length must be non-zero");
    }
    if (std::abs(objectDistance) < 1e-15) {
        throw std::invalid_argument("Object distance must be non-zero");
    }

    double reciprocal = 1.0/focalLength + 1.0/objectDistance;
    if (std::abs(reciprocal) < 1e-15) {
        throw std::invalid_argument("Image at infinity");
    }

    return 1.0 / reciprocal;
}

/**
 * @brief Calculate object distance from lens formula
 *
 * 1/u = 1/v - 1/f
 *
 * @param focalLength Focal length (m)
 * @param imageDistance Image distance (m)
 * @return Object distance (m)
 * @throws std::invalid_argument if calculation is invalid
 */
inline double objectDistance(double focalLength, double imageDistance) {
    if (std::abs(focalLength) < 1e-15) {
        throw std::invalid_argument("Focal length must be non-zero");
    }
    if (std::abs(imageDistance) < 1e-15) {
        throw std::invalid_argument("Image distance must be non-zero");
    }

    double reciprocal = 1.0/imageDistance - 1.0/focalLength;
    if (std::abs(reciprocal) < 1e-15) {
        throw std::invalid_argument("Object at infinity");
    }

    return 1.0 / reciprocal;
}

/**
 * @brief Calculate lens power (inverse of focal length)
 *
 * P = 1/f
 *
 * Power is measured in diopters (D = m⁻¹)
 *
 * @param focalLength Focal length (m)
 * @return Power (diopters)
 * @throws std::invalid_argument if focal length is zero
 */
inline double lensPower(double focalLength) {
    if (std::abs(focalLength) < 1e-15) {
        throw std::invalid_argument("Focal length must be non-zero");
    }
    return 1.0 / focalLength;
}

/**
 * @brief Calculate focal length from lens power
 *
 * f = 1/P
 *
 * @param power Lens power (diopters)
 * @return Focal length (m)
 * @throws std::invalid_argument if power is zero
 */
inline double focalLengthFromPower(double power) {
    if (std::abs(power) < 1e-15) {
        throw std::invalid_argument("Power must be non-zero");
    }
    return 1.0 / power;
}

// ============================================================================
// PRINCIPAL FOCUS AND FOCAL LENGTH
// ============================================================================

/**
 * @brief Calculate focal length from lensmaker's equation
 *
 * 1/f = (n-1)(1/R₁ - 1/R₂)
 *
 * R is positive for convex surface, negative for concave
 *
 * @param refractiveIndex Refractive index of lens material
 * @param radius1 Radius of curvature of first surface (m)
 * @param radius2 Radius of curvature of second surface (m)
 * @return Focal length (m)
 * @throws std::invalid_argument if radii are zero
 */
inline double lensmakersEquation(double refractiveIndex, double radius1, double radius2) {
    if (std::abs(radius1) < 1e-15 || std::abs(radius2) < 1e-15) {
        throw std::invalid_argument("Radii of curvature must be non-zero");
    }

    double reciprocalF = (refractiveIndex - 1.0) * (1.0/radius1 - 1.0/radius2);

    if (std::abs(reciprocalF) < 1e-15) {
        throw std::invalid_argument("Infinite focal length (plane parallel plate)");
    }

    return 1.0 / reciprocalF;
}

/**
 * @brief Calculate focal length for thin lens in medium
 *
 * 1/f = (n_lens/n_medium - 1)(1/R₁ - 1/R₂)
 *
 * @param lensIndex Refractive index of lens
 * @param mediumIndex Refractive index of surrounding medium
 * @param radius1 First surface radius (m)
 * @param radius2 Second surface radius (m)
 * @return Focal length (m)
 * @throws std::invalid_argument if parameters are invalid
 */
inline double focalLengthInMedium(double lensIndex, double mediumIndex,
                                  double radius1, double radius2) {
    if (mediumIndex <= 0.0) {
        throw std::invalid_argument("Medium refractive index must be positive");
    }
    if (std::abs(radius1) < 1e-15 || std::abs(radius2) < 1e-15) {
        throw std::invalid_argument("Radii must be non-zero");
    }

    double relativeIndex = lensIndex / mediumIndex;
    double reciprocalF = (relativeIndex - 1.0) * (1.0/radius1 - 1.0/radius2);

    if (std::abs(reciprocalF) < 1e-15) {
        throw std::invalid_argument("Infinite focal length");
    }

    return 1.0 / reciprocalF;
}

// ============================================================================
// MAGNIFICATION
// ============================================================================

/**
 * @brief Calculate linear magnification
 *
 * m = v/u = h_i/h_o
 *
 * @param imageDistance Image distance (m)
 * @param objectDistance Object distance (m)
 * @return Magnification (negative for inverted image)
 * @throws std::invalid_argument if object distance is zero
 */
inline double linearMagnification(double imageDistance, double objectDistance) {
    if (std::abs(objectDistance) < 1e-15) {
        throw std::invalid_argument("Object distance must be non-zero");
    }
    return imageDistance / objectDistance;
}

/**
 * @brief Calculate image height from magnification
 *
 * h_i = m × h_o
 *
 * @param objectHeight Object height (m)
 * @param magnification Magnification
 * @return Image height (m)
 */
inline double imageHeight(double objectHeight, double magnification) {
    return magnification * objectHeight;
}

/**
 * @brief Calculate magnification from focal length and object distance
 *
 * m = f / (f + u)
 *
 * @param focalLength Focal length (m)
 * @param objectDistance Object distance (m)
 * @return Magnification
 * @throws std::invalid_argument if denominator is zero
 */
inline double magnificationFromFocal(double focalLength, double objectDistance) {
    double denominator = focalLength + objectDistance;
    if (std::abs(denominator) < 1e-15) {
        throw std::invalid_argument("Invalid configuration");
    }
    return focalLength / denominator;
}

// ============================================================================
// MIRROR EQUATIONS
// ============================================================================

/**
 * @brief Calculate focal length of spherical mirror
 *
 * f = R/2
 *
 * where R is radius of curvature
 *
 * @param radiusOfCurvature Radius of curvature (m)
 * @return Focal length (m)
 */
inline double mirrorFocalLength(double radiusOfCurvature) {
    return radiusOfCurvature / 2.0;
}

/**
 * @brief Mirror formula (same form as lens formula)
 *
 * 1/f = 1/v + 1/u
 *
 * @param objectDistance Object distance (m)
 * @param imageDistance Image distance (m)
 * @return Focal length (m)
 * @throws std::invalid_argument if distances are zero
 */
inline double mirrorFormula(double objectDistance, double imageDistance) {
    if (std::abs(objectDistance) < 1e-15 || std::abs(imageDistance) < 1e-15) {
        throw std::invalid_argument("Distances must be non-zero");
    }
    return 1.0 / (1.0/imageDistance + 1.0/objectDistance);
}

// ============================================================================
// OPTICAL INSTRUMENTS
// ============================================================================

/**
 * @brief Calculate magnifying power of simple microscope (magnifying glass)
 *
 * M = 1 + D/f
 *
 * where D is least distance of distinct vision (typically 25 cm)
 *
 * @param focalLength Focal length of lens (m)
 * @param leastDistance Least distance of distinct vision (m), default 0.25 m
 * @return Magnifying power
 * @throws std::invalid_argument if focal length is zero
 */
inline double simpleMicroscopeMagnification(double focalLength, double leastDistance = 0.25) {
    if (std::abs(focalLength) < 1e-15) {
        throw std::invalid_argument("Focal length must be non-zero");
    }
    return 1.0 + leastDistance / focalLength;
}

/**
 * @brief Calculate magnifying power of compound microscope
 *
 * M = (v_o/u_o) × (D/f_e)
 *
 * where v_o, u_o are image and object distances for objective,
 * D is least distance, f_e is eyepiece focal length
 *
 * @param objectiveMagnification Magnification of objective
 * @param eyepieceFocal Focal length of eyepiece (m)
 * @param leastDistance Least distance of distinct vision (m)
 * @return Total magnifying power
 * @throws std::invalid_argument if eyepiece focal length is zero
 */
inline double compoundMicroscopeMagnification(double objectiveMagnification,
                                              double eyepieceFocal,
                                              double leastDistance = 0.25) {
    if (std::abs(eyepieceFocal) < 1e-15) {
        throw std::invalid_argument("Eyepiece focal length must be non-zero");
    }
    return objectiveMagnification * (leastDistance / eyepieceFocal);
}

/**
 * @brief Calculate magnifying power of astronomical telescope (normal adjustment)
 *
 * M = f_o / f_e
 *
 * where f_o is objective focal length, f_e is eyepiece focal length
 *
 * @param objectiveFocal Focal length of objective (m)
 * @param eyepieceFocal Focal length of eyepiece (m)
 * @return Magnifying power
 * @throws std::invalid_argument if eyepiece focal length is zero
 */
inline double telescopeMagnification(double objectiveFocal, double eyepieceFocal) {
    if (std::abs(eyepieceFocal) < 1e-15) {
        throw std::invalid_argument("Eyepiece focal length must be non-zero");
    }
    return objectiveFocal / eyepieceFocal;
}

/**
 * @brief Calculate length of astronomical telescope (normal adjustment)
 *
 * L = f_o + f_e
 *
 * @param objectiveFocal Focal length of objective (m)
 * @param eyepieceFocal Focal length of eyepiece (m)
 * @return Tube length (m)
 */
inline double telescopeLength(double objectiveFocal, double eyepieceFocal) {
    return objectiveFocal + eyepieceFocal;
}

/**
 * @brief Calculate resolving power of telescope
 *
 * R = D / (1.22λ)
 *
 * where D is aperture diameter
 *
 * @param apertureDiameter Diameter of objective (m)
 * @param wavelength Wavelength of light (m)
 * @return Resolving power (rad⁻¹)
 * @throws std::invalid_argument if wavelength is zero
 */
inline double telescopeResolvingPower(double apertureDiameter, double wavelength) {
    if (std::abs(wavelength) < 1e-15) {
        throw std::invalid_argument("Wavelength must be non-zero");
    }
    return apertureDiameter / (1.22 * wavelength);
}

// ============================================================================
// LENS COMBINATIONS
// ============================================================================

/**
 * @brief Calculate combined focal length of two thin lenses in contact
 *
 * 1/F = 1/f₁ + 1/f₂
 *
 * @param focal1 Focal length of first lens (m)
 * @param focal2 Focal length of second lens (m)
 * @return Combined focal length (m)
 * @throws std::invalid_argument if focal lengths are zero
 */
inline double combinedFocalLength(double focal1, double focal2) {
    if (std::abs(focal1) < 1e-15 || std::abs(focal2) < 1e-15) {
        throw std::invalid_argument("Focal lengths must be non-zero");
    }
    double reciprocalF = 1.0/focal1 + 1.0/focal2;
    if (std::abs(reciprocalF) < 1e-15) {
        throw std::invalid_argument("Infinite combined focal length");
    }
    return 1.0 / reciprocalF;
}

/**
 * @brief Calculate combined power of lenses in contact
 *
 * P = P₁ + P₂
 *
 * @param power1 Power of first lens (diopters)
 * @param power2 Power of second lens (diopters)
 * @return Combined power (diopters)
 */
inline double combinedPower(double power1, double power2) {
    return power1 + power2;
}

/**
 * @brief Calculate effective focal length of separated lenses
 *
 * 1/F = 1/f₁ + 1/f₂ - d/(f₁f₂)
 *
 * @param focal1 Focal length of first lens (m)
 * @param focal2 Focal length of second lens (m)
 * @param separation Distance between lenses (m)
 * @return Effective focal length (m)
 * @throws std::invalid_argument if focal lengths are zero
 */
inline double separatedLensesFocalLength(double focal1, double focal2, double separation) {
    if (std::abs(focal1) < 1e-15 || std::abs(focal2) < 1e-15) {
        throw std::invalid_argument("Focal lengths must be non-zero");
    }

    double reciprocalF = 1.0/focal1 + 1.0/focal2 - separation/(focal1 * focal2);

    if (std::abs(reciprocalF) < 1e-15) {
        throw std::invalid_argument("Infinite effective focal length");
    }

    return 1.0 / reciprocalF;
}

} // namespace optics
} // namespace physics

#endif // PHYSICS_OPTICS_HPP
