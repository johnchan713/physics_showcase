#ifndef PHYSICS_ADVANCED_OPTICS_HPP
#define PHYSICS_ADVANCED_OPTICS_HPP

#include <cmath>
#include <complex>
#include <vector>
#include <stdexcept>
#include <array>

/**
 * @file advanced_optics.hpp
 * @brief Comprehensive implementation of advanced optical phenomena and techniques
 *
 * This module implements:
 * - Wave optics: phase difference, spherical waves
 * - Refraction at spherical surfaces
 * - Aspheric lenses
 * - Stops and apertures
 * - Lens aberrations (spherical, chromatic, coma, astigmatism, etc.)
 * - Laser Doppler Velocimetry (LDV)
 * - Interferometry (Michelson, Mach-Zehnder, Fabry-Perot)
 * - Optical Coherence Tomography (OCT)
 * - Diffraction (single slit, multiple slits, gratings)
 * - Fourier optics
 * - Coherent and incoherent transfer functions
 * - Phase-modulated sinusoidal gratings
 * - Holography
 * - Optical triangulation
 * - Moire technique
 * - Photoelasticity
 * - Polarized light (Stokes parameters, Mueller matrices)
 * - Ellipsometry
 * - Fringe analysis
 * - Phase unwrapping
 * - Fiber optics (numerical aperture, modes, dispersion)
 *
 * All calculations use SI units unless otherwise specified.
 */

namespace physics {
namespace advanced_optics {

// ============================================================================
// CONSTANTS
// ============================================================================

namespace constants {
    constexpr double SPEED_OF_LIGHT = 299792458.0;  // m/s
    constexpr double PLANCK_CONSTANT = 6.62607015e-34;  // J·s
    constexpr double PI = 3.141592653589793;
}

// ============================================================================
// PHASE DIFFERENCE AND WAVE OPTICS
// ============================================================================

/**
 * @brief Calculate phase difference between two waves
 *
 * Δφ = (2π/λ) × Δx
 *
 * @param pathDifference Optical path difference (m)
 * @param wavelength Wavelength (m)
 * @return Phase difference (radians)
 */
inline double phaseDifference(double pathDifference, double wavelength) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return (2.0 * constants::PI * pathDifference) / wavelength;
}

/**
 * @brief Calculate optical path difference from phase difference
 *
 * @param phaseDiff Phase difference (radians)
 * @param wavelength Wavelength (m)
 * @return Optical path difference (m)
 */
inline double opticalPathDifference(double phaseDiff, double wavelength) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return (phaseDiff * wavelength) / (2.0 * constants::PI);
}

/**
 * @brief Calculate optical path length in medium
 *
 * OPL = n × d
 *
 * @param refractiveIndex Refractive index
 * @param geometricPath Geometric path length (m)
 * @return Optical path length (m)
 */
inline double opticalPathLength(double refractiveIndex, double geometricPath) {
    if (refractiveIndex <= 0.0) {
        throw std::invalid_argument("Refractive index must be positive");
    }
    return refractiveIndex * geometricPath;
}

/**
 * @brief Check for constructive interference condition
 *
 * Constructive: Δφ = 2πn, where n is integer
 *
 * @param phaseDiff Phase difference (radians)
 * @param tolerance Tolerance for comparison
 * @return true if constructive interference occurs
 */
inline bool isConstructiveInterference(double phaseDiff, double tolerance = 1e-6) {
    double n = phaseDiff / (2.0 * constants::PI);
    return std::abs(n - std::round(n)) < tolerance;
}

/**
 * @brief Check for destructive interference condition
 *
 * Destructive: Δφ = (2n+1)π, where n is integer
 *
 * @param phaseDiff Phase difference (radians)
 * @param tolerance Tolerance for comparison
 * @return true if destructive interference occurs
 */
inline bool isDestructiveInterference(double phaseDiff, double tolerance = 1e-6) {
    double n = (phaseDiff / constants::PI - 1.0) / 2.0;
    return std::abs(n - std::round(n)) < tolerance;
}

// ============================================================================
// SPHERICAL WAVES
// ============================================================================

/**
 * @brief Calculate amplitude of spherical wave at distance r
 *
 * A(r) = A₀/r
 *
 * @param initialAmplitude Initial amplitude at unit distance
 * @param distance Distance from source (m)
 * @return Amplitude at distance r
 */
inline double sphericalWaveAmplitude(double initialAmplitude, double distance) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return initialAmplitude / distance;
}

/**
 * @brief Calculate intensity of spherical wave
 *
 * I(r) = I₀/r²
 *
 * @param initialIntensity Intensity at unit distance (W/m²)
 * @param distance Distance from source (m)
 * @return Intensity at distance r (W/m²)
 */
inline double sphericalWaveIntensity(double initialIntensity, double distance) {
    if (distance <= 0.0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return initialIntensity / (distance * distance);
}

/**
 * @brief Calculate wavefront curvature radius
 *
 * For spherical wave emanating from point source
 *
 * @param distance Distance from source (m)
 * @return Radius of curvature (m)
 */
inline double wavefrontCurvature(double distance) {
    return distance;
}

/**
 * @brief Calculate Gouy phase shift for focused Gaussian beam
 *
 * ψ(z) = arctan(z/z₀)
 *
 * @param z Axial distance from waist (m)
 * @param rayleighRange Rayleigh range z₀ (m)
 * @return Gouy phase (radians)
 */
inline double gouyPhase(double z, double rayleighRange) {
    if (rayleighRange <= 0.0) {
        throw std::invalid_argument("Rayleigh range must be positive");
    }
    return std::atan(z / rayleighRange);
}

// ============================================================================
// REFRACTION AT SPHERICAL SURFACE
// ============================================================================

/**
 * @brief Refraction formula for spherical surface
 *
 * n₂/v - n₁/u = (n₂-n₁)/R
 *
 * @param n1 Refractive index of first medium
 * @param n2 Refractive index of second medium
 * @param objectDistance Object distance (m)
 * @param radiusCurvature Radius of curvature (m, positive for convex)
 * @return Image distance (m)
 */
inline double sphericalSurfaceRefraction(double n1, double n2,
                                         double objectDistance,
                                         double radiusCurvature) {
    if (n1 <= 0.0 || n2 <= 0.0) {
        throw std::invalid_argument("Refractive indices must be positive");
    }
    if (std::abs(objectDistance) < 1e-15) {
        throw std::invalid_argument("Object distance must be non-zero");
    }
    if (std::abs(radiusCurvature) < 1e-15) {
        throw std::invalid_argument("Radius must be non-zero");
    }

    double reciprocalV = (n2 - n1) / radiusCurvature + n1 / objectDistance;

    if (std::abs(reciprocalV) < 1e-15) {
        throw std::invalid_argument("Image at infinity");
    }

    return n2 / reciprocalV;
}

/**
 * @brief Calculate power of spherical refracting surface
 *
 * P = (n₂-n₁)/R
 *
 * @param n1 First medium refractive index
 * @param n2 Second medium refractive index
 * @param radiusCurvature Radius of curvature (m)
 * @return Power (diopters)
 */
inline double sphericalSurfacePower(double n1, double n2, double radiusCurvature) {
    if (std::abs(radiusCurvature) < 1e-15) {
        throw std::invalid_argument("Radius must be non-zero");
    }
    return (n2 - n1) / radiusCurvature;
}

// ============================================================================
// ASPHERIC LENSES
// ============================================================================

/**
 * @brief Calculate sag of aspheric surface
 *
 * z(r) = cr²/[1+√(1-(1+k)c²r²)] + Σ Aᵢr²ⁱ
 *
 * @param r Radial distance from optical axis (m)
 * @param curvature Curvature c = 1/R (m⁻¹)
 * @param conicConstant Conic constant k
 * @param asphericCoeffs Higher-order aspheric coefficients
 * @return Surface sag (m)
 */
inline double asphericSag(double r, double curvature, double conicConstant,
                          const std::vector<double>& asphericCoeffs = {}) {
    double r2 = r * r;
    double c2r2 = curvature * curvature * r2;

    // Base conic term
    double sqrt_term = 1.0 - (1.0 + conicConstant) * c2r2;
    if (sqrt_term < 0.0) {
        throw std::invalid_argument("Invalid conic/curvature combination");
    }

    double sag = (curvature * r2) / (1.0 + std::sqrt(sqrt_term));

    // Add higher-order aspheric terms
    double r_power = r2;
    for (size_t i = 0; i < asphericCoeffs.size(); ++i) {
        r_power *= r2;  // r⁴, r⁶, r⁸, ...
        sag += asphericCoeffs[i] * r_power;
    }

    return sag;
}

/**
 * @brief Classify conic section from conic constant
 *
 * k < -1: Hyperbola
 * k = -1: Parabola
 * -1 < k < 0: Ellipse
 * k = 0: Sphere
 * k > 0: Oblate ellipsoid
 *
 * @param conicConstant Conic constant k
 * @return Classification code: 0=sphere, 1=ellipse, 2=parabola, 3=hyperbola, 4=oblate
 */
inline int classifyConicSurface(double conicConstant) {
    constexpr double tolerance = 1e-9;

    if (std::abs(conicConstant) < tolerance) return 0;  // Sphere
    if (std::abs(conicConstant + 1.0) < tolerance) return 2;  // Parabola
    if (conicConstant < -1.0) return 3;  // Hyperbola
    if (conicConstant > 0.0) return 4;  // Oblate ellipsoid
    return 1;  // Ellipse
}

// ============================================================================
// STOPS AND APERTURES
// ============================================================================

/**
 * @brief Calculate numerical aperture (NA)
 *
 * NA = n sin(θ)
 *
 * @param refractiveIndex Refractive index of medium
 * @param halfAngle Half-angle of acceptance cone (radians)
 * @return Numerical aperture
 */
inline double numericalAperture(double refractiveIndex, double halfAngle) {
    if (refractiveIndex <= 0.0) {
        throw std::invalid_argument("Refractive index must be positive");
    }
    return refractiveIndex * std::sin(halfAngle);
}

/**
 * @brief Calculate f-number (f/#)
 *
 * f/# = f/D
 *
 * @param focalLength Focal length (m)
 * @param apertureDiameter Aperture diameter (m)
 * @return F-number
 */
inline double fNumber(double focalLength, double apertureDiameter) {
    if (apertureDiameter <= 0.0) {
        throw std::invalid_argument("Aperture diameter must be positive");
    }
    return focalLength / apertureDiameter;
}

/**
 * @brief Calculate f-number from numerical aperture
 *
 * f/# = 1/(2·NA)
 *
 * @param na Numerical aperture
 * @return F-number
 */
inline double fNumberFromNA(double na) {
    if (na <= 0.0) {
        throw std::invalid_argument("Numerical aperture must be positive");
    }
    return 1.0 / (2.0 * na);
}

/**
 * @brief Calculate depth of field
 *
 * DOF = 2·λ·(f/#)²
 *
 * @param wavelength Wavelength (m)
 * @param fNum F-number
 * @return Depth of field (m)
 */
inline double depthOfField(double wavelength, double fNum) {
    if (wavelength <= 0.0 || fNum <= 0.0) {
        throw std::invalid_argument("Wavelength and f-number must be positive");
    }
    return 2.0 * wavelength * fNum * fNum;
}

/**
 * @brief Calculate entrance pupil diameter
 *
 * @param apertureDiameter Physical aperture diameter (m)
 * @param magnification Magnification from aperture to entrance pupil
 * @return Entrance pupil diameter (m)
 */
inline double entrancePupilDiameter(double apertureDiameter, double magnification) {
    return std::abs(magnification) * apertureDiameter;
}

// ============================================================================
// LENS ABERRATIONS
// ============================================================================

/**
 * @brief Calculate longitudinal spherical aberration
 *
 * LSA ≈ -h⁴/(128·n·f³·(n-1)²)
 *
 * @param aperture Aperture height (m)
 * @param focalLength Focal length (m)
 * @param refractiveIndex Lens refractive index
 * @return Longitudinal spherical aberration (m)
 */
inline double longitudinalSphericalAberration(double aperture, double focalLength,
                                              double refractiveIndex) {
    if (focalLength <= 0.0 || refractiveIndex <= 0.0) {
        throw std::invalid_argument("Focal length and refractive index must be positive");
    }

    double h4 = aperture * aperture * aperture * aperture;
    double f3 = focalLength * focalLength * focalLength;
    double n_minus_1 = refractiveIndex - 1.0;

    return -h4 / (128.0 * refractiveIndex * f3 * n_minus_1 * n_minus_1);
}

/**
 * @brief Calculate chromatic aberration (longitudinal)
 *
 * δf = f·(n_F - n_C)/(n_D - 1) = f/ν
 *
 * where ν is Abbe number
 *
 * @param focalLength Focal length (m)
 * @param abbeNumber Abbe number (V_d)
 * @return Chromatic aberration (m)
 */
inline double chromaticAberration(double focalLength, double abbeNumber) {
    if (abbeNumber <= 0.0) {
        throw std::invalid_argument("Abbe number must be positive");
    }
    return focalLength / abbeNumber;
}

/**
 * @brief Calculate Abbe number (dispersion)
 *
 * V_d = (n_d - 1)/(n_F - n_C)
 *
 * @param n_d Refractive index at d-line (587.6 nm)
 * @param n_F Refractive index at F-line (486.1 nm)
 * @param n_C Refractive index at C-line (656.3 nm)
 * @return Abbe number
 */
inline double abbeNumber(double n_d, double n_F, double n_C) {
    double denominator = n_F - n_C;
    if (std::abs(denominator) < 1e-15) {
        throw std::invalid_argument("Dispersion must be non-zero");
    }
    return (n_d - 1.0) / denominator;
}

/**
 * @brief Calculate coma aberration
 *
 * Coma ≈ (h³/2f²)·sin(α)
 *
 * @param aperture Aperture height (m)
 * @param focalLength Focal length (m)
 * @param fieldAngle Field angle (radians)
 * @return Coma blur diameter (m)
 */
inline double comaAberration(double aperture, double focalLength, double fieldAngle) {
    if (focalLength <= 0.0) {
        throw std::invalid_argument("Focal length must be positive");
    }

    double h3 = aperture * aperture * aperture;
    double f2 = focalLength * focalLength;

    return (h3 / (2.0 * f2)) * std::sin(fieldAngle);
}

/**
 * @brief Calculate astigmatism
 *
 * Astigmatism = f·tan²(θ)
 *
 * @param focalLength Focal length (m)
 * @param fieldAngle Field angle (radians)
 * @return Astigmatic difference (m)
 */
inline double astigmatism(double focalLength, double fieldAngle) {
    double tanTheta = std::tan(fieldAngle);
    return focalLength * tanTheta * tanTheta;
}

/**
 * @brief Calculate field curvature (Petzval sum)
 *
 * 1/R_P = Σ(1/(n_i·f_i))
 *
 * @param focalLengths Vector of focal lengths (m)
 * @param refractiveIndices Vector of refractive indices
 * @return Petzval radius of curvature (m)
 */
inline double petzvalCurvature(const std::vector<double>& focalLengths,
                               const std::vector<double>& refractiveIndices) {
    if (focalLengths.size() != refractiveIndices.size()) {
        throw std::invalid_argument("Vectors must have same size");
    }

    double petzvalSum = 0.0;
    for (size_t i = 0; i < focalLengths.size(); ++i) {
        if (std::abs(focalLengths[i]) < 1e-15 || refractiveIndices[i] <= 0.0) {
            throw std::invalid_argument("Invalid focal length or refractive index");
        }
        petzvalSum += 1.0 / (refractiveIndices[i] * focalLengths[i]);
    }

    if (std::abs(petzvalSum) < 1e-15) {
        throw std::invalid_argument("Zero Petzval sum");
    }

    return 1.0 / petzvalSum;
}

/**
 * @brief Calculate distortion percentage
 *
 * Distortion = 100·(h' - h'_ideal)/h'_ideal
 *
 * @param actualHeight Actual image height (m)
 * @param idealHeight Ideal image height (m)
 * @return Distortion percentage
 */
inline double distortion(double actualHeight, double idealHeight) {
    if (std::abs(idealHeight) < 1e-15) {
        throw std::invalid_argument("Ideal height must be non-zero");
    }
    return 100.0 * (actualHeight - idealHeight) / idealHeight;
}

// ============================================================================
// LASER DOPPLER VELOCIMETRY (LDV)
// ============================================================================

/**
 * @brief Calculate Doppler frequency shift
 *
 * f_D = (2·v·sin(θ/2))/λ
 *
 * @param velocity Particle velocity (m/s)
 * @param wavelength Laser wavelength (m)
 * @param beamAngle Angle between beams (radians)
 * @return Doppler frequency (Hz)
 */
inline double dopplerFrequencyShift(double velocity, double wavelength, double beamAngle) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return (2.0 * velocity * std::sin(beamAngle / 2.0)) / wavelength;
}

/**
 * @brief Calculate fringe spacing in LDV measurement volume
 *
 * d_f = λ/(2·sin(θ/2))
 *
 * @param wavelength Laser wavelength (m)
 * @param beamAngle Angle between beams (radians)
 * @return Fringe spacing (m)
 */
inline double ldvFringeSpacing(double wavelength, double beamAngle) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    double sinHalfAngle = std::sin(beamAngle / 2.0);
    if (std::abs(sinHalfAngle) < 1e-15) {
        throw std::invalid_argument("Beam angle must be non-zero");
    }
    return wavelength / (2.0 * sinHalfAngle);
}

/**
 * @brief Calculate particle velocity from Doppler frequency
 *
 * v = (f_D·λ)/(2·sin(θ/2))
 *
 * @param dopplerFreq Doppler frequency (Hz)
 * @param wavelength Laser wavelength (m)
 * @param beamAngle Angle between beams (radians)
 * @return Particle velocity (m/s)
 */
inline double velocityFromDoppler(double dopplerFreq, double wavelength, double beamAngle) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    double sinHalfAngle = std::sin(beamAngle / 2.0);
    if (std::abs(sinHalfAngle) < 1e-15) {
        throw std::invalid_argument("Beam angle must be non-zero");
    }
    return (dopplerFreq * wavelength) / (2.0 * sinHalfAngle);
}

// ============================================================================
// INTERFEROMETRY
// ============================================================================

/**
 * @brief Calculate interference fringe visibility (contrast)
 *
 * V = (I_max - I_min)/(I_max + I_min)
 *
 * @param maxIntensity Maximum intensity
 * @param minIntensity Minimum intensity
 * @return Visibility (0 to 1)
 */
inline double fringeVisibility(double maxIntensity, double minIntensity) {
    if (maxIntensity < 0.0 || minIntensity < 0.0) {
        throw std::invalid_argument("Intensities must be non-negative");
    }
    double sum = maxIntensity + minIntensity;
    if (sum < 1e-15) {
        throw std::invalid_argument("Sum of intensities must be positive");
    }
    return (maxIntensity - minIntensity) / sum;
}

/**
 * @brief Calculate intensity in two-beam interference
 *
 * I = I₁ + I₂ + 2√(I₁I₂)cos(Δφ)
 *
 * @param intensity1 Intensity of beam 1 (W/m²)
 * @param intensity2 Intensity of beam 2 (W/m²)
 * @param phaseDiff Phase difference (radians)
 * @return Total intensity (W/m²)
 */
inline double twoBeamInterferenceIntensity(double intensity1, double intensity2,
                                           double phaseDiff) {
    if (intensity1 < 0.0 || intensity2 < 0.0) {
        throw std::invalid_argument("Intensities must be non-negative");
    }
    return intensity1 + intensity2 +
           2.0 * std::sqrt(intensity1 * intensity2) * std::cos(phaseDiff);
}

/**
 * @brief Calculate fringe spacing in double-slit interference
 *
 * β = λD/d
 *
 * @param wavelength Wavelength (m)
 * @param screenDistance Distance to screen (m)
 * @param slitSeparation Slit separation (m)
 * @return Fringe spacing (m)
 */
inline double doubleSlit FringeSpacing(double wavelength, double screenDistance,
                                       double slitSeparation) {
    if (wavelength <= 0.0 || screenDistance <= 0.0 || slitSeparation <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return (wavelength * screenDistance) / slitSeparation;
}

/**
 * @brief Calculate finesse of Fabry-Perot interferometer
 *
 * F = FSR/FWHM = π√R/(1-R)
 *
 * @param reflectivity Mirror reflectivity (0 to 1)
 * @return Finesse
 */
inline double fabryPerotFinesse(double reflectivity) {
    if (reflectivity < 0.0 || reflectivity >= 1.0) {
        throw std::invalid_argument("Reflectivity must be in [0, 1)");
    }
    return (constants::PI * std::sqrt(reflectivity)) / (1.0 - reflectivity);
}

/**
 * @brief Calculate free spectral range (FSR) of Fabry-Perot
 *
 * FSR = c/(2nL)
 *
 * @param cavityLength Cavity length (m)
 * @param refractiveIndex Refractive index
 * @return Free spectral range (Hz)
 */
inline double freeSpectralRange(double cavityLength, double refractiveIndex = 1.0) {
    if (cavityLength <= 0.0 || refractiveIndex <= 0.0) {
        throw std::invalid_argument("Cavity length and refractive index must be positive");
    }
    return constants::SPEED_OF_LIGHT / (2.0 * refractiveIndex * cavityLength);
}

/**
 * @brief Calculate transmitted intensity through Fabry-Perot
 *
 * I_t/I_i = 1/(1 + F·sin²(δ/2))
 *
 * where F = 4R/(1-R)² is coefficient of finesse
 *
 * @param incidentIntensity Incident intensity
 * @param reflectivity Mirror reflectivity
 * @param phaseShift Round-trip phase shift (radians)
 * @return Transmitted intensity
 */
inline double fabryPerotTransmission(double incidentIntensity, double reflectivity,
                                     double phaseShift) {
    if (reflectivity < 0.0 || reflectivity >= 1.0) {
        throw std::invalid_argument("Reflectivity must be in [0, 1)");
    }

    double coeffFinesse = (4.0 * reflectivity) / ((1.0 - reflectivity) * (1.0 - reflectivity));
    double sinHalf = std::sin(phaseShift / 2.0);

    return incidentIntensity / (1.0 + coeffFinesse * sinHalf * sinHalf);
}

// ============================================================================
// OPTICAL COHERENCE TOMOGRAPHY (OCT)
// ============================================================================

/**
 * @brief Calculate axial resolution in OCT
 *
 * Δz = (2·ln(2)/π)·(λ₀²/Δλ)
 *
 * @param centerWavelength Center wavelength (m)
 * @param bandwidth Spectral bandwidth (m)
 * @return Axial resolution (m)
 */
inline double octAxialResolution(double centerWavelength, double bandwidth) {
    if (centerWavelength <= 0.0 || bandwidth <= 0.0) {
        throw std::invalid_argument("Wavelength and bandwidth must be positive");
    }
    return (2.0 * std::log(2.0) / constants::PI) *
           (centerWavelength * centerWavelength / bandwidth);
}

/**
 * @brief Calculate lateral resolution in OCT
 *
 * Δx = 4λf/(πd)
 *
 * @param wavelength Wavelength (m)
 * @param focalLength Focal length (m)
 * @param beamDiameter Beam diameter (m)
 * @return Lateral resolution (m)
 */
inline double octLateralResolution(double wavelength, double focalLength,
                                   double beamDiameter) {
    if (wavelength <= 0.0 || focalLength <= 0.0 || beamDiameter <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return (4.0 * wavelength * focalLength) / (constants::PI * beamDiameter);
}

/**
 * @brief Calculate depth of focus in OCT
 *
 * DOF = πΔx²/(2λ)
 *
 * @param lateralResolution Lateral resolution (m)
 * @param wavelength Wavelength (m)
 * @return Depth of focus (m)
 */
inline double octDepthOfFocus(double lateralResolution, double wavelength) {
    if (lateralResolution <= 0.0 || wavelength <= 0.0) {
        throw std::invalid_argument("Parameters must be positive");
    }
    return (constants::PI * lateralResolution * lateralResolution) / (2.0 * wavelength);
}

/**
 * @brief Calculate sensitivity in OCT
 *
 * S = 10·log₁₀(R_sample/NEP)
 *
 * @param sampleReflectivity Sample reflectivity
 * @param noiseEquivalentPower Noise equivalent power
 * @return Sensitivity (dB)
 */
inline double octSensitivity(double sampleReflectivity, double noiseEquivalentPower) {
    if (sampleReflectivity <= 0.0 || noiseEquivalentPower <= 0.0) {
        throw std::invalid_argument("Reflectivity and NEP must be positive");
    }
    return 10.0 * std::log10(sampleReflectivity / noiseEquivalentPower);
}

// ============================================================================
// DIFFRACTION - SINGLE SLIT
// ============================================================================

/**
 * @brief Calculate intensity pattern for single slit diffraction
 *
 * I(θ) = I₀·[sin(β)/β]²
 * where β = (π·a·sin(θ))/λ
 *
 * @param slitWidth Slit width (m)
 * @param wavelength Wavelength (m)
 * @param angle Observation angle (radians)
 * @param peakIntensity Peak intensity (W/m²)
 * @return Intensity at angle θ (W/m²)
 */
inline double singleSlitIntensity(double slitWidth, double wavelength,
                                  double angle, double peakIntensity) {
    if (slitWidth <= 0.0 || wavelength <= 0.0) {
        throw std::invalid_argument("Slit width and wavelength must be positive");
    }

    double beta = (constants::PI * slitWidth * std::sin(angle)) / wavelength;

    if (std::abs(beta) < 1e-10) {
        return peakIntensity;  // Central maximum
    }

    double sinc = std::sin(beta) / beta;
    return peakIntensity * sinc * sinc;
}

/**
 * @brief Calculate position of nth minimum in single slit pattern
 *
 * a·sin(θ) = n·λ
 *
 * @param slitWidth Slit width (m)
 * @param wavelength Wavelength (m)
 * @param order Order of minimum (n = ±1, ±2, ...)
 * @return Angle to nth minimum (radians)
 */
inline double singleSlitMinimumAngle(double slitWidth, double wavelength, int order) {
    if (slitWidth <= 0.0 || wavelength <= 0.0) {
        throw std::invalid_argument("Slit width and wavelength must be positive");
    }
    if (order == 0) {
        throw std::invalid_argument("Order must be non-zero");
    }

    double ratio = (order * wavelength) / slitWidth;
    if (std::abs(ratio) > 1.0) {
        throw std::invalid_argument("Minimum does not exist for this order");
    }

    return std::asin(ratio);
}

/**
 * @brief Calculate angular width of central maximum
 *
 * Δθ = 2λ/a
 *
 * @param slitWidth Slit width (m)
 * @param wavelength Wavelength (m)
 * @return Angular width (radians)
 */
inline double singleSlitCentralMaxWidth(double slitWidth, double wavelength) {
    if (slitWidth <= 0.0 || wavelength <= 0.0) {
        throw std::invalid_argument("Slit width and wavelength must be positive");
    }
    return (2.0 * wavelength) / slitWidth;
}

// ============================================================================
// GRATING EQUATION
// ============================================================================

/**
 * @brief Calculate diffraction angle using grating equation
 *
 * d·(sin(θ_m) - sin(θ_i)) = m·λ
 *
 * @param gratingSpacing Grating spacing d (m)
 * @param wavelength Wavelength (m)
 * @param order Diffraction order m
 * @param incidentAngle Incident angle (radians)
 * @return Diffraction angle θ_m (radians)
 */
inline double gratingDiffractionAngle(double gratingSpacing, double wavelength,
                                      int order, double incidentAngle = 0.0) {
    if (gratingSpacing <= 0.0 || wavelength <= 0.0) {
        throw std::invalid_argument("Grating spacing and wavelength must be positive");
    }

    double sinTheta = (order * wavelength / gratingSpacing) + std::sin(incidentAngle);

    if (std::abs(sinTheta) > 1.0) {
        throw std::invalid_argument("Diffraction order does not exist");
    }

    return std::asin(sinTheta);
}

/**
 * @brief Calculate angular dispersion of grating
 *
 * dθ/dλ = m/(d·cos(θ))
 *
 * @param gratingSpacing Grating spacing (m)
 * @param order Diffraction order
 * @param diffractionAngle Diffraction angle (radians)
 * @return Angular dispersion (rad/m)
 */
inline double gratingAngularDispersion(double gratingSpacing, int order,
                                       double diffractionAngle) {
    if (gratingSpacing <= 0.0) {
        throw std::invalid_argument("Grating spacing must be positive");
    }

    double cosTheta = std::cos(diffractionAngle);
    if (std::abs(cosTheta) < 1e-15) {
        throw std::invalid_argument("Invalid diffraction angle");
    }

    return order / (gratingSpacing * cosTheta);
}

/**
 * @brief Calculate resolving power of grating
 *
 * R = λ/Δλ = m·N
 *
 * where N is number of illuminated rulings
 *
 * @param order Diffraction order
 * @param numberOfRulings Number of illuminated grating lines
 * @return Resolving power
 */
inline double gratingResolvingPower(int order, int numberOfRulings) {
    if (numberOfRulings <= 0) {
        throw std::invalid_argument("Number of rulings must be positive");
    }
    return std::abs(order) * numberOfRulings;
}

/**
 * @brief Calculate free spectral range of grating
 *
 * FSR = λ/m
 *
 * @param wavelength Operating wavelength (m)
 * @param order Diffraction order
 * @return Free spectral range (m)
 */
inline double gratingFSR(double wavelength, int order) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    if (order == 0) {
        throw std::invalid_argument("Order must be non-zero");
    }
    return wavelength / std::abs(order);
}

// ============================================================================
// FOURIER OPTICS
// ============================================================================

/**
 * @brief Calculate spatial frequency
 *
 * f_x = sin(θ)/λ
 *
 * @param angle Propagation angle (radians)
 * @param wavelength Wavelength (m)
 * @return Spatial frequency (m⁻¹)
 */
inline double spatialFrequency(double angle, double wavelength) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return std::sin(angle) / wavelength;
}

/**
 * @brief Calculate diffraction-limited spot size (Airy disk)
 *
 * d = 2.44·λ·f/D
 *
 * @param wavelength Wavelength (m)
 * @param focalLength Focal length (m)
 * @param apertureDiameter Aperture diameter (m)
 * @return Airy disk diameter (m)
 */
inline double airyDiskDiameter(double wavelength, double focalLength,
                               double apertureDiameter) {
    if (wavelength <= 0.0 || focalLength <= 0.0 || apertureDiameter <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return (2.44 * wavelength * focalLength) / apertureDiameter;
}

/**
 * @brief Calculate Rayleigh criterion for resolution
 *
 * θ_min = 1.22·λ/D
 *
 * @param wavelength Wavelength (m)
 * @param apertureDiameter Aperture diameter (m)
 * @return Minimum resolvable angle (radians)
 */
inline double rayleighCriterion(double wavelength, double apertureDiameter) {
    if (wavelength <= 0.0 || apertureDiameter <= 0.0) {
        throw std::invalid_argument("Wavelength and aperture must be positive");
    }
    return (1.22 * wavelength) / apertureDiameter;
}

/**
 * @brief Calculate cutoff spatial frequency for diffraction-limited system
 *
 * f_c = D/(λ·f)
 *
 * @param apertureDiameter Aperture diameter (m)
 * @param wavelength Wavelength (m)
 * @param focalLength Focal length (m)
 * @return Cutoff frequency (m⁻¹)
 */
inline double cutoffSpatialFrequency(double apertureDiameter, double wavelength,
                                     double focalLength) {
    if (wavelength <= 0.0 || focalLength <= 0.0 || apertureDiameter <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return apertureDiameter / (wavelength * focalLength);
}

// ============================================================================
// COHERENT TRANSFER FUNCTION (CTF)
// ============================================================================

/**
 * @brief Calculate coherent transfer function (CTF) for diffraction-limited system
 *
 * CTF(f) = circ(f/f_c)
 *
 * @param spatialFreq Spatial frequency (m⁻¹)
 * @param cutoffFreq Cutoff frequency (m⁻¹)
 * @return CTF value (0 or 1 for ideal case)
 */
inline double coherentTransferFunction(double spatialFreq, double cutoffFreq) {
    if (cutoffFreq <= 0.0) {
        throw std::invalid_argument("Cutoff frequency must be positive");
    }
    return (std::abs(spatialFreq) <= cutoffFreq) ? 1.0 : 0.0;
}

/**
 * @brief Calculate phase transfer function (PTF) for defocused system
 *
 * PTF(f) = -π·W₂₀·f²
 *
 * where W₂₀ is defocus coefficient
 *
 * @param spatialFreq Spatial frequency (m⁻¹)
 * @param defocusCoeff Defocus coefficient (waves)
 * @return Phase (radians)
 */
inline double phaseTransferFunction(double spatialFreq, double defocusCoeff) {
    return -constants::PI * defocusCoeff * spatialFreq * spatialFreq;
}

/**
 * @brief Calculate Strehl ratio for aberrated system
 *
 * S ≈ exp(-(2πσ/λ)²)
 *
 * where σ is RMS wavefront error
 *
 * @param rmsWavefrontError RMS wavefront error (m)
 * @param wavelength Wavelength (m)
 * @return Strehl ratio (0 to 1)
 */
inline double strehlRatio(double rmsWavefrontError, double wavelength) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    double ratio = (2.0 * constants::PI * rmsWavefrontError) / wavelength;
    return std::exp(-ratio * ratio);
}

// ============================================================================
// INCOHERENT TRANSFER FUNCTION (OTF/MTF)
// ============================================================================

/**
 * @brief Calculate optical transfer function (OTF) for diffraction-limited system
 *
 * OTF(f) = (2/π)[arccos(f/f_c) - (f/f_c)√(1-(f/f_c)²)]
 *
 * @param spatialFreq Spatial frequency (m⁻¹)
 * @param cutoffFreq Cutoff frequency (m⁻¹)
 * @return OTF value (0 to 1)
 */
inline double opticalTransferFunction(double spatialFreq, double cutoffFreq) {
    if (cutoffFreq <= 0.0) {
        throw std::invalid_argument("Cutoff frequency must be positive");
    }

    double f_norm = std::abs(spatialFreq) / cutoffFreq;

    if (f_norm >= 1.0) {
        return 0.0;
    }

    if (f_norm < 1e-10) {
        return 1.0;
    }

    double sqrt_term = std::sqrt(1.0 - f_norm * f_norm);
    return (2.0 / constants::PI) * (std::acos(f_norm) - f_norm * sqrt_term);
}

/**
 * @brief Calculate modulation transfer function (MTF)
 *
 * MTF = |OTF|
 *
 * @param spatialFreq Spatial frequency (m⁻¹)
 * @param cutoffFreq Cutoff frequency (m⁻¹)
 * @return MTF value (0 to 1)
 */
inline double modulationTransferFunction(double spatialFreq, double cutoffFreq) {
    return std::abs(opticalTransferFunction(spatialFreq, cutoffFreq));
}

/**
 * @brief Calculate contrast reduction through optical system
 *
 * C_out = MTF × C_in
 *
 * @param inputContrast Input contrast (0 to 1)
 * @param mtf MTF value at spatial frequency of interest
 * @return Output contrast (0 to 1)
 */
inline double contrastReduction(double inputContrast, double mtf) {
    if (inputContrast < 0.0 || inputContrast > 1.0) {
        throw std::invalid_argument("Input contrast must be in [0, 1]");
    }
    if (mtf < 0.0 || mtf > 1.0) {
        throw std::invalid_argument("MTF must be in [0, 1]");
    }
    return inputContrast * mtf;
}

// ============================================================================
// PHASE-MODULATED SINUSOIDAL GRATING
// ============================================================================

/**
 * @brief Calculate intensity distribution of phase grating
 *
 * I(x) = I₀·{1 + Σ J_n(φ₀)·cos(2πnx/Λ)}
 *
 * For small phase modulation φ₀:
 * I(x) ≈ I₀·{1 + φ₀·cos(2πx/Λ)}
 *
 * @param position Position x (m)
 * @param period Grating period Λ (m)
 * @param phaseModulation Maximum phase modulation φ₀ (radians)
 * @param baseIntensity Base intensity I₀
 * @return Intensity at position x
 */
inline double phaseGratingIntensity(double position, double period,
                                    double phaseModulation, double baseIntensity) {
    if (period <= 0.0 || baseIntensity < 0.0) {
        throw std::invalid_argument("Period and intensity must be positive");
    }

    // Small angle approximation
    if (std::abs(phaseModulation) < 0.1) {
        return baseIntensity * (1.0 + phaseModulation *
                               std::cos(2.0 * constants::PI * position / period));
    }

    // For larger modulation, use Bessel function approximation (first few orders)
    // J_0, J_1, J_2 contributions
    double arg = 2.0 * constants::PI * position / period;

    // Simplified using first-order Bessel functions
    // J_0(x) ≈ 1 - x²/4, J_1(x) ≈ x/2 for small x
    double j0 = 1.0 - (phaseModulation * phaseModulation) / 4.0;
    double j1 = phaseModulation / 2.0;

    return baseIntensity * (j0 + 2.0 * j1 * std::cos(arg));
}

/**
 * @brief Calculate diffraction efficiency of phase grating
 *
 * η_m = J_m²(φ₀)
 *
 * where J_m is mth order Bessel function
 *
 * For ±1 orders: η₁ = sin²(φ₀/2)
 *
 * @param order Diffraction order
 * @param phaseModulation Maximum phase modulation (radians)
 * @return Diffraction efficiency
 */
inline double phaseGratingEfficiency(int order, double phaseModulation) {
    // First-order approximation for ±1 orders
    if (std::abs(order) == 1) {
        double sinHalf = std::sin(phaseModulation / 2.0);
        return sinHalf * sinHalf;
    }

    // Zero order
    if (order == 0) {
        double cosHalf = std::cos(phaseModulation / 2.0);
        return cosHalf * cosHalf;
    }

    // Higher orders - simplified approximation
    return 0.0;  // Negligible for higher orders in simple case
}

/**
 * @brief Calculate phase modulation from refractive index modulation
 *
 * φ(x) = (2π/λ)·Δn(x)·d
 *
 * @param indexModulation Refractive index modulation amplitude
 * @param thickness Grating thickness (m)
 * @param wavelength Wavelength (m)
 * @return Phase modulation amplitude (radians)
 */
inline double phaseModulationFromIndex(double indexModulation, double thickness,
                                       double wavelength) {
    if (wavelength <= 0.0 || thickness < 0.0) {
        throw std::invalid_argument("Wavelength must be positive, thickness non-negative");
    }
    return (2.0 * constants::PI * indexModulation * thickness) / wavelength;
}

// ============================================================================
// HOLOGRAPHY
// ============================================================================

/**
 * @brief Calculate holographic recording intensity
 *
 * I = |E_object + E_reference|² = I_o + I_r + 2√(I_o·I_r)·cos(Δφ)
 *
 * @param objectIntensity Object beam intensity
 * @param referenceIntensity Reference beam intensity
 * @param phaseDiff Phase difference (radians)
 * @return Recorded intensity
 */
inline double holographicRecordingIntensity(double objectIntensity,
                                            double referenceIntensity,
                                            double phaseDiff) {
    if (objectIntensity < 0.0 || referenceIntensity < 0.0) {
        throw std::invalid_argument("Intensities must be non-negative");
    }
    return objectIntensity + referenceIntensity +
           2.0 * std::sqrt(objectIntensity * referenceIntensity) * std::cos(phaseDiff);
}

/**
 * @brief Calculate holographic fringe spacing
 *
 * Λ = λ/(2·sin(θ/2))
 *
 * @param wavelength Recording wavelength (m)
 * @param beamAngle Angle between object and reference beams (radians)
 * @return Fringe spacing (m)
 */
inline double holographicFringeSpacing(double wavelength, double beamAngle) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    double sinHalf = std::sin(beamAngle / 2.0);
    if (std::abs(sinHalf) < 1e-15) {
        throw std::invalid_argument("Beam angle must be non-zero");
    }
    return wavelength / (2.0 * sinHalf);
}

/**
 * @brief Calculate diffraction efficiency of volume hologram
 *
 * η = sin²(πΔn·d/(λ·cos(θ)))
 *
 * @param indexModulation Refractive index modulation
 * @param thickness Hologram thickness (m)
 * @param wavelength Wavelength (m)
 * @param braggAngle Bragg angle inside medium (radians)
 * @return Diffraction efficiency (0 to 1)
 */
inline double volumeHologramEfficiency(double indexModulation, double thickness,
                                       double wavelength, double braggAngle) {
    if (wavelength <= 0.0 || thickness <= 0.0) {
        throw std::invalid_argument("Wavelength and thickness must be positive");
    }

    double arg = (constants::PI * indexModulation * thickness) /
                 (wavelength * std::cos(braggAngle));
    double sinArg = std::sin(arg);

    return sinArg * sinArg;
}

/**
 * @brief Calculate Bragg angle for volume hologram
 *
 * sin(θ_B) = λ/(2Λ)
 *
 * @param wavelength Readout wavelength (m)
 * @param fringeSpacing Fringe spacing in hologram (m)
 * @return Bragg angle (radians)
 */
inline double braggAngle(double wavelength, double fringeSpacing) {
    if (wavelength <= 0.0 || fringeSpacing <= 0.0) {
        throw std::invalid_argument("Wavelength and fringe spacing must be positive");
    }

    double sinBragg = wavelength / (2.0 * fringeSpacing);
    if (std::abs(sinBragg) > 1.0) {
        throw std::invalid_argument("Bragg condition cannot be satisfied");
    }

    return std::asin(sinBragg);
}

// ============================================================================
// OPTICAL TRIANGULATION
// ============================================================================

/**
 * @brief Calculate distance using triangulation
 *
 * Z = (b·f)/d
 *
 * where b is baseline, f is focal length, d is disparity
 *
 * @param baseline Distance between two measurement points (m)
 * @param focalLength Focal length of optics (m)
 * @param disparity Measured disparity (m or pixels × pixel size)
 * @return Distance to object (m)
 */
inline double triangulationDistance(double baseline, double focalLength, double disparity) {
    if (baseline <= 0.0 || focalLength <= 0.0) {
        throw std::invalid_argument("Baseline and focal length must be positive");
    }
    if (std::abs(disparity) < 1e-15) {
        throw std::invalid_argument("Disparity must be non-zero");
    }
    return (baseline * focalLength) / disparity;
}

/**
 * @brief Calculate triangulation angle
 *
 * @param distance Distance to object (m)
 * @param baseline Distance between measurement points (m)
 * @return Triangulation angle (radians)
 */
inline double triangulationAngle(double distance, double baseline) {
    if (distance <= 0.0 || baseline <= 0.0) {
        throw std::invalid_argument("Distance and baseline must be positive");
    }
    return std::atan(baseline / distance);
}

/**
 * @brief Calculate depth resolution in triangulation
 *
 * δZ = Z²/(b·f)·δd
 *
 * @param distance Distance to object (m)
 * @param baseline Baseline (m)
 * @param focalLength Focal length (m)
 * @param disparityResolution Disparity measurement resolution (m)
 * @return Depth resolution (m)
 */
inline double triangulationDepthResolution(double distance, double baseline,
                                           double focalLength, double disparityResolution) {
    if (distance <= 0.0 || baseline <= 0.0 || focalLength <= 0.0) {
        throw std::invalid_argument("Distance, baseline, and focal length must be positive");
    }
    return (distance * distance * disparityResolution) / (baseline * focalLength);
}

// ============================================================================
// MOIRE TECHNIQUE
// ============================================================================

/**
 * @brief Calculate moire fringe spacing
 *
 * Λ_m = (p₁·p₂)/|p₁ - p₂|
 *
 * @param period1 Period of first grating (m)
 * @param period2 Period of second grating (m)
 * @return Moire fringe spacing (m)
 */
inline double moireFringeSpacing(double period1, double period2) {
    if (period1 <= 0.0 || period2 <= 0.0) {
        throw std::invalid_argument("Periods must be positive");
    }

    double diff = std::abs(period1 - period2);
    if (diff < 1e-15) {
        throw std::invalid_argument("Periods too similar - infinite moire spacing");
    }

    return (period1 * period2) / diff;
}

/**
 * @brief Calculate moire rotation angle
 *
 * θ_m = |θ₁ - θ₂|
 *
 * For small angles: Λ_m ≈ p/(2·sin(θ_m/2))
 *
 * @param angle1 Orientation of first grating (radians)
 * @param angle2 Orientation of second grating (radians)
 * @return Moire angle (radians)
 */
inline double moireRotationAngle(double angle1, double angle2) {
    return std::abs(angle1 - angle2);
}

/**
 * @brief Calculate moire magnification
 *
 * M = p₁/(p₁ - p₂)
 *
 * @param period1 Period of first (reference) grating (m)
 * @param period2 Period of second (deformed) grating (m)
 * @return Moire magnification
 */
inline double moireMagnification(double period1, double period2) {
    if (period1 <= 0.0 || period2 <= 0.0) {
        throw std::invalid_argument("Periods must be positive");
    }

    double diff = period1 - period2;
    if (std::abs(diff) < 1e-15) {
        throw std::invalid_argument("Infinite moire magnification");
    }

    return period1 / diff;
}

/**
 * @brief Calculate strain from moire pattern
 *
 * ε = Δp/p = 1/M
 *
 * @param moireMag Moire magnification
 * @return Strain (dimensionless)
 */
inline double strainFromMoire(double moireMag) {
    if (std::abs(moireMag) < 1e-15) {
        throw std::invalid_argument("Moire magnification must be non-zero");
    }
    return 1.0 / moireMag;
}

// ============================================================================
// PHOTOELASTICITY
// ============================================================================

/**
 * @brief Calculate stress-optic law
 *
 * Δn = C·σ
 *
 * where C is stress-optic coefficient
 *
 * @param stressOpticCoeff Stress-optic coefficient (Pa⁻¹)
 * @param stress Applied stress (Pa)
 * @return Birefringence (dimensionless)
 */
inline double photoelasticBirefringence(double stressOpticCoeff, double stress) {
    return stressOpticCoeff * stress;
}

/**
 * @brief Calculate retardation in photoelastic material
 *
 * δ = 2π·Δn·t/λ = 2π·C·σ·t/λ
 *
 * @param stressOpticCoeff Stress-optic coefficient (Pa⁻¹)
 * @param stress Principal stress difference (Pa)
 * @param thickness Sample thickness (m)
 * @param wavelength Wavelength (m)
 * @return Retardation (radians)
 */
inline double photoelasticRetardation(double stressOpticCoeff, double stress,
                                      double thickness, double wavelength) {
    if (wavelength <= 0.0 || thickness < 0.0) {
        throw std::invalid_argument("Wavelength must be positive, thickness non-negative");
    }

    double birefringence = photoelasticBirefringence(stressOpticCoeff, stress);
    return (2.0 * constants::PI * birefringence * thickness) / wavelength;
}

/**
 * @brief Calculate fringe order in photoelasticity
 *
 * N = δ/(2π) = C·σ·t/λ
 *
 * @param retardation Retardation (radians)
 * @return Fringe order N
 */
inline double photoelasticFringeOrder(double retardation) {
    return retardation / (2.0 * constants::PI);
}

/**
 * @brief Calculate stress from fringe order
 *
 * σ = N·λ/(C·t)
 *
 * @param fringeOrder Fringe order N
 * @param wavelength Wavelength (m)
 * @param stressOpticCoeff Stress-optic coefficient (Pa⁻¹)
 * @param thickness Sample thickness (m)
 * @return Stress (Pa)
 */
inline double stressFromFringeOrder(double fringeOrder, double wavelength,
                                    double stressOpticCoeff, double thickness) {
    if (std::abs(stressOpticCoeff) < 1e-30 || thickness <= 0.0 || wavelength <= 0.0) {
        throw std::invalid_argument("Invalid parameters");
    }
    return (fringeOrder * wavelength) / (stressOpticCoeff * thickness);
}

// ============================================================================
// POLARIZED LIGHT - STOKES PARAMETERS
// ============================================================================

/**
 * @class StokesVector
 * @brief Represents polarization state using Stokes parameters
 *
 * S = [S₀, S₁, S₂, S₃]ᵀ
 * S₀ = total intensity
 * S₁ = horizontal vs vertical linear polarization
 * S₂ = +45° vs -45° linear polarization
 * S₃ = right vs left circular polarization
 */
class StokesVector {
public:
    double S0, S1, S2, S3;

    StokesVector(double s0 = 0.0, double s1 = 0.0, double s2 = 0.0, double s3 = 0.0)
        : S0(s0), S1(s1), S2(s2), S3(s3) {}

    /**
     * @brief Calculate degree of polarization
     *
     * DOP = √(S₁² + S₂² + S₃²)/S₀
     */
    double degreeOfPolarization() const {
        if (S0 < 1e-15) {
            throw std::invalid_argument("S0 must be positive");
        }
        return std::sqrt(S1*S1 + S2*S2 + S3*S3) / S0;
    }

    /**
     * @brief Calculate degree of linear polarization
     *
     * DOLP = √(S₁² + S₂²)/S₀
     */
    double degreeOfLinearPolarization() const {
        if (S0 < 1e-15) {
            throw std::invalid_argument("S0 must be positive");
        }
        return std::sqrt(S1*S1 + S2*S2) / S0;
    }

    /**
     * @brief Calculate degree of circular polarization
     *
     * DOCP = S₃/S₀
     */
    double degreeOfCircularPolarization() const {
        if (S0 < 1e-15) {
            throw std::invalid_argument("S0 must be positive");
        }
        return S3 / S0;
    }

    /**
     * @brief Calculate polarization ellipse orientation angle
     *
     * ψ = (1/2)·arctan(S₂/S₁)
     */
    double orientationAngle() const {
        return 0.5 * std::atan2(S2, S1);
    }

    /**
     * @brief Calculate polarization ellipse ellipticity angle
     *
     * χ = (1/2)·arcsin(S₃/√(S₁² + S₂² + S₃²))
     */
    double ellipticityAngle() const {
        double magnitude = std::sqrt(S1*S1 + S2*S2 + S3*S3);
        if (magnitude < 1e-15) {
            return 0.0;
        }
        return 0.5 * std::asin(S3 / magnitude);
    }
};

/**
 * @brief Create Stokes vector for linearly polarized light
 *
 * @param intensity Total intensity
 * @param angle Polarization angle from horizontal (radians)
 * @return Stokes vector
 */
inline StokesVector linearlyPolarizedStokes(double intensity, double angle) {
    return StokesVector(
        intensity,
        intensity * std::cos(2.0 * angle),
        intensity * std::sin(2.0 * angle),
        0.0
    );
}

/**
 * @brief Create Stokes vector for circularly polarized light
 *
 * @param intensity Total intensity
 * @param rightHanded true for right-handed, false for left-handed
 * @return Stokes vector
 */
inline StokesVector circularlyPolarizedStokes(double intensity, bool rightHanded) {
    return StokesVector(
        intensity,
        0.0,
        0.0,
        rightHanded ? intensity : -intensity
    );
}

// ============================================================================
// ELLIPSOMETRY
// ============================================================================

/**
 * @brief Calculate ellipsometric parameters psi and delta
 *
 * tan(ψ)·exp(iΔ) = r_p/r_s
 *
 * @param r_p_magnitude Magnitude of p-polarization reflection coefficient
 * @param r_s_magnitude Magnitude of s-polarization reflection coefficient
 * @param phaseDifference Phase difference between p and s (radians)
 * @return std::pair<psi, delta> in radians
 */
inline std::pair<double, double> ellipsometricParameters(double r_p_magnitude,
                                                         double r_s_magnitude,
                                                         double phaseDifference) {
    if (r_s_magnitude < 1e-15) {
        throw std::invalid_argument("r_s magnitude must be positive");
    }

    double psi = std::atan(r_p_magnitude / r_s_magnitude);
    double delta = phaseDifference;

    return {psi, delta};
}

/**
 * @brief Calculate refractive index from ellipsometric data
 *
 * For simple case (transparent substrate in air):
 * n = sin(θ)·√(1 + tan²(ψ)·sin²(Δ))
 *
 * @param incidentAngle Incident angle (radians)
 * @param psi Ellipsometric psi (radians)
 * @param delta Ellipsometric delta (radians)
 * @return Refractive index
 */
inline double refractiveIndexFromEllipsometry(double incidentAngle, double psi,
                                              double delta) {
    double tanPsi = std::tan(psi);
    double sinDelta = std::sin(delta);
    double sinTheta = std::sin(incidentAngle);

    return sinTheta * std::sqrt(1.0 + tanPsi * tanPsi * sinDelta * sinDelta);
}

/**
 * @brief Calculate film thickness from ellipsometry
 *
 * For thin film: d = λ·Δ/(4π·n·cos(θ_t))
 *
 * Simplified approximation
 *
 * @param wavelength Wavelength (m)
 * @param delta Ellipsometric delta (radians)
 * @param filmIndex Film refractive index
 * @param refractedAngle Angle in film (radians)
 * @return Film thickness (m)
 */
inline double filmThicknessFromEllipsometry(double wavelength, double delta,
                                            double filmIndex, double refractedAngle) {
    if (wavelength <= 0.0 || filmIndex <= 0.0) {
        throw std::invalid_argument("Wavelength and refractive index must be positive");
    }

    return (wavelength * delta) / (4.0 * constants::PI * filmIndex * std::cos(refractedAngle));
}

// ============================================================================
// FRINGE ANALYSIS
// ============================================================================

/**
 * @brief Calculate phase from fringe intensity (phase-stepping)
 *
 * φ = arctan[(I₃ - I₁)/(I₂ - I₄)]
 *
 * For 4-step algorithm with π/2 phase shifts
 *
 * @param I1 Intensity at phase 0
 * @param I2 Intensity at phase π/2
 * @param I3 Intensity at phase π
 * @param I4 Intensity at phase 3π/2
 * @return Phase (radians, -π to π)
 */
inline double phaseFrom4Step(double I1, double I2, double I3, double I4) {
    return std::atan2(I3 - I1, I2 - I4);
}

/**
 * @brief Calculate modulation (fringe contrast) in interferogram
 *
 * γ = 2·√[(I₃-I₁)² + (I₂-I₄)²]/(I₁+I₂+I₃+I₄)
 *
 * @param I1 Intensity at phase 0
 * @param I2 Intensity at phase π/2
 * @param I3 Intensity at phase π
 * @param I4 Intensity at phase 3π/2
 * @return Modulation (0 to 1)
 */
inline double fringeModulation(double I1, double I2, double I3, double I4) {
    double sum = I1 + I2 + I3 + I4;
    if (sum < 1e-15) {
        throw std::invalid_argument("Total intensity must be positive");
    }

    double diff1 = I3 - I1;
    double diff2 = I2 - I4;

    return 2.0 * std::sqrt(diff1*diff1 + diff2*diff2) / sum;
}

/**
 * @brief Calculate average intensity (DC component)
 *
 * I_avg = (I₁ + I₂ + I₃ + I₄)/4
 *
 * @param I1 Intensity at phase 0
 * @param I2 Intensity at phase π/2
 * @param I3 Intensity at phase π
 * @param I4 Intensity at phase 3π/2
 * @return Average intensity
 */
inline double averageIntensity(double I1, double I2, double I3, double I4) {
    return (I1 + I2 + I3 + I4) / 4.0;
}

// ============================================================================
// PHASE UNWRAPPING
// ============================================================================

/**
 * @brief Unwrap phase difference (1D)
 *
 * Adds/subtracts 2π to maintain continuity
 *
 * @param wrappedPhase Wrapped phase (radians, -π to π)
 * @param previousUnwrapped Previous unwrapped phase (radians)
 * @return Unwrapped phase (radians)
 */
inline double unwrapPhase(double wrappedPhase, double previousUnwrapped) {
    double diff = wrappedPhase - std::fmod(previousUnwrapped, 2.0 * constants::PI);

    // Normalize diff to [-π, π]
    while (diff > constants::PI) diff -= 2.0 * constants::PI;
    while (diff < -constants::PI) diff += 2.0 * constants::PI;

    return previousUnwrapped + diff;
}

/**
 * @brief Calculate phase gradient quality map
 *
 * Used to guide phase unwrapping
 *
 * @param phase1 Phase at point 1 (radians)
 * @param phase2 Phase at adjacent point 2 (radians)
 * @param spacing Spatial spacing between points (m)
 * @return Phase gradient (rad/m)
 */
inline double phaseGradient(double phase1, double phase2, double spacing) {
    if (spacing <= 0.0) {
        throw std::invalid_argument("Spacing must be positive");
    }

    double diff = phase2 - phase1;

    // Wrap difference to [-π, π]
    while (diff > constants::PI) diff -= 2.0 * constants::PI;
    while (diff < -constants::PI) diff += 2.0 * constants::PI;

    return diff / spacing;
}

/**
 * @brief Calculate reliability for phase unwrapping (second difference)
 *
 * R = |∇²φ|
 *
 * Lower values indicate higher reliability
 *
 * @param phasePrev Phase at previous point (radians)
 * @param phaseCurr Phase at current point (radians)
 * @param phaseNext Phase at next point (radians)
 * @return Reliability metric (lower is better)
 */
inline double phaseUnwrappingReliability(double phasePrev, double phaseCurr,
                                         double phaseNext) {
    // Second difference approximation
    double diff1 = phaseCurr - phasePrev;
    double diff2 = phaseNext - phaseCurr;

    // Wrap to [-π, π]
    while (diff1 > constants::PI) diff1 -= 2.0 * constants::PI;
    while (diff1 < -constants::PI) diff1 += 2.0 * constants::PI;
    while (diff2 > constants::PI) diff2 -= 2.0 * constants::PI;
    while (diff2 < -constants::PI) diff2 += 2.0 * constants::PI;

    return std::abs(diff2 - diff1);
}

// ============================================================================
// FIBER OPTICS
// ============================================================================

/**
 * @brief Calculate numerical aperture of optical fiber
 *
 * NA = √(n_core² - n_cladding²)
 *
 * @param coreIndex Core refractive index
 * @param claddingIndex Cladding refractive index
 * @return Numerical aperture
 */
inline double fiberNumericalAperture(double coreIndex, double claddingIndex) {
    if (coreIndex <= 0.0 || claddingIndex <= 0.0) {
        throw std::invalid_argument("Refractive indices must be positive");
    }
    if (coreIndex < claddingIndex) {
        throw std::invalid_argument("Core index must be greater than cladding index");
    }

    return std::sqrt(coreIndex*coreIndex - claddingIndex*claddingIndex);
}

/**
 * @brief Calculate acceptance angle of fiber
 *
 * θ_max = arcsin(NA)
 *
 * @param numericalAperture Numerical aperture
 * @return Maximum acceptance angle (radians)
 */
inline double fiberAcceptanceAngle(double numericalAperture) {
    if (numericalAperture < 0.0 || numericalAperture > 1.0) {
        throw std::invalid_argument("Numerical aperture must be in [0, 1]");
    }
    return std::asin(numericalAperture);
}

/**
 * @brief Calculate number of modes in step-index fiber
 *
 * M ≈ V²/2 for large V
 * V = (2πa/λ)·NA (normalized frequency)
 *
 * @param coreRadius Core radius (m)
 * @param wavelength Wavelength (m)
 * @param numericalAperture Numerical aperture
 * @return Approximate number of modes
 */
inline double fiberNumberOfModes(double coreRadius, double wavelength,
                                 double numericalAperture) {
    if (coreRadius <= 0.0 || wavelength <= 0.0 || numericalAperture <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }

    double V = (2.0 * constants::PI * coreRadius * numericalAperture) / wavelength;

    if (V < 2.405) {
        return 1.0;  // Single mode
    }

    return (V * V) / 2.0;
}

/**
 * @brief Calculate V-number (normalized frequency) of fiber
 *
 * V = (2πa/λ)·NA
 *
 * @param coreRadius Core radius (m)
 * @param wavelength Wavelength (m)
 * @param numericalAperture Numerical aperture
 * @return V-number (dimensionless)
 */
inline double fiberVNumber(double coreRadius, double wavelength,
                           double numericalAperture) {
    if (coreRadius <= 0.0 || wavelength <= 0.0 || numericalAperture <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    return (2.0 * constants::PI * coreRadius * numericalAperture) / wavelength;
}

/**
 * @brief Calculate cutoff wavelength for single-mode operation
 *
 * λ_c = (2πa/2.405)·NA
 *
 * @param coreRadius Core radius (m)
 * @param numericalAperture Numerical aperture
 * @return Cutoff wavelength (m)
 */
inline double fiberCutoffWavelength(double coreRadius, double numericalAperture) {
    if (coreRadius <= 0.0 || numericalAperture <= 0.0) {
        throw std::invalid_argument("Core radius and NA must be positive");
    }
    return (2.0 * constants::PI * coreRadius * numericalAperture) / 2.405;
}

/**
 * @brief Calculate intermodal dispersion in multimode fiber
 *
 * Δτ = (L·NA²)/(2c·n_core)
 *
 * @param fiberLength Fiber length (m)
 * @param numericalAperture Numerical aperture
 * @param coreIndex Core refractive index
 * @return Pulse broadening (s)
 */
inline double fiberIntermodalDispersion(double fiberLength, double numericalAperture,
                                        double coreIndex) {
    if (fiberLength <= 0.0 || numericalAperture <= 0.0 || coreIndex <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }

    return (fiberLength * numericalAperture * numericalAperture) /
           (2.0 * constants::SPEED_OF_LIGHT * coreIndex);
}

/**
 * @brief Calculate chromatic dispersion parameter
 *
 * D = -(λ/c)·(d²n/dλ²)
 *
 * Typical units: ps/(nm·km)
 *
 * @param wavelength Wavelength (m)
 * @param secondDerivative Second derivative of refractive index (m⁻²)
 * @return Dispersion parameter (s/m²)
 */
inline double fiberChromaticDispersion(double wavelength, double secondDerivative) {
    if (wavelength <= 0.0) {
        throw std::invalid_argument("Wavelength must be positive");
    }
    return -(wavelength / constants::SPEED_OF_LIGHT) * secondDerivative;
}

/**
 * @brief Calculate pulse broadening due to chromatic dispersion
 *
 * Δτ = D·L·Δλ
 *
 * @param dispersionParam Dispersion parameter D (s/m²)
 * @param fiberLength Fiber length (m)
 * @param spectralWidth Spectral width Δλ (m)
 * @return Pulse broadening (s)
 */
inline double fiberPulseBroadening(double dispersionParam, double fiberLength,
                                   double spectralWidth) {
    return std::abs(dispersionParam * fiberLength * spectralWidth);
}

/**
 * @brief Calculate fiber attenuation in dB
 *
 * α_dB = -10·log₁₀(P_out/P_in)/L
 *
 * @param inputPower Input power (W)
 * @param outputPower Output power (W)
 * @param fiberLength Fiber length (m)
 * @return Attenuation coefficient (dB/m)
 */
inline double fiberAttenuation(double inputPower, double outputPower,
                               double fiberLength) {
    if (inputPower <= 0.0 || outputPower <= 0.0 || fiberLength <= 0.0) {
        throw std::invalid_argument("All parameters must be positive");
    }
    if (outputPower > inputPower) {
        throw std::invalid_argument("Output power cannot exceed input power");
    }

    return -10.0 * std::log10(outputPower / inputPower) / fiberLength;
}

/**
 * @brief Calculate output power given attenuation
 *
 * P_out = P_in·10^(-αL/10)
 *
 * @param inputPower Input power (W)
 * @param attenuation Attenuation coefficient (dB/m)
 * @param fiberLength Fiber length (m)
 * @return Output power (W)
 */
inline double fiberOutputPower(double inputPower, double attenuation,
                               double fiberLength) {
    if (inputPower <= 0.0 || fiberLength < 0.0) {
        throw std::invalid_argument("Input power must be positive, length non-negative");
    }

    return inputPower * std::pow(10.0, -attenuation * fiberLength / 10.0);
}

} // namespace advanced_optics
} // namespace physics

#endif // PHYSICS_ADVANCED_OPTICS_HPP
