#ifndef PHYSICS_ADVANCED_FLUID_DYNAMICS_VORTICITY_HPP
#define PHYSICS_ADVANCED_FLUID_DYNAMICS_VORTICITY_HPP

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <functional>
#include <stdexcept>

/**
 * @file vorticity.hpp
 * @brief Vorticity dynamics and circulation
 *
 * Implements:
 * - Vorticity vector ω = ∇×u
 * - Vorticity transport equation
 * - Circulation Γ = ∮ u·dl
 * - Kelvin's circulation theorem
 * - Vortex stretching
 * - Biot-Savart law
 */

namespace physics::advanced::fluid_dynamics {

/**
 * @class Vorticity
 * @brief Vorticity vector and related operations
 *
 * Vorticity ω = ∇×u measures local rotation rate
 */
class Vorticity {
public:
    /**
     * @brief Compute vorticity from velocity gradient
     *
     * ωᵢ = εᵢⱼₖ ∂uₖ/∂xⱼ
     *
     * or in 3D:
     * ω = [∂w/∂y - ∂v/∂z, ∂u/∂z - ∂w/∂x, ∂v/∂x - ∂u/∂y]ᵀ
     *
     * @param velocity_gradient ∂uᵢ/∂xⱼ (3×3 matrix)
     * @return Vorticity vector ω (1/s)
     */
    static Eigen::Vector3d compute(const Eigen::Matrix3d& velocity_gradient) {
        Eigen::Vector3d omega;

        // ωx = ∂w/∂y - ∂v/∂z
        omega(0) = velocity_gradient(2, 1) - velocity_gradient(1, 2);

        // ωy = ∂u/∂z - ∂w/∂x
        omega(1) = velocity_gradient(0, 2) - velocity_gradient(2, 0);

        // ωz = ∂v/∂x - ∂u/∂y
        omega(2) = velocity_gradient(1, 0) - velocity_gradient(0, 1);

        return omega;
    }

    /**
     * @brief Vorticity magnitude
     *
     * |ω| = √(ωx² + ωy² + ωz²)
     *
     * @param omega Vorticity vector (1/s)
     * @return Magnitude (1/s)
     */
    static double magnitude(const Eigen::Vector3d& omega) {
        return omega.norm();
    }

    /**
     * @brief Enstrophy (vorticity squared)
     *
     * Ω = ½ω·ω = ½|ω|²
     *
     * Enstrophy is conserved in 2D inviscid flow
     *
     * @param omega Vorticity vector (1/s)
     * @return Enstrophy (1/s²)
     */
    static double enstrophy(const Eigen::Vector3d& omega) {
        return 0.5 * omega.squaredNorm();
    }

    /**
     * @brief Check if flow is irrotational
     *
     * Irrotational: ω = 0 everywhere
     *
     * @param omega Vorticity vector
     * @param tolerance Threshold for zero
     * @return true if irrotational
     */
    static bool isIrrotational(const Eigen::Vector3d& omega,
                               double tolerance = 1e-10) {
        return omega.norm() < tolerance;
    }

    /**
     * @brief Vortex lines
     *
     * Lines tangent to vorticity vector ω at every point
     * Analogous to streamlines but for vorticity field
     */
    static Eigen::Vector3d vortexLineDirection(const Eigen::Vector3d& omega) {
        double norm = omega.norm();
        if (norm < 1e-10) {
            return Eigen::Vector3d::Zero();
        }
        return omega / norm;
    }
};

/**
 * @class VorticityTransportEquation
 * @brief Vorticity transport equation (incompressible)
 *
 * Dω/Dt = (ω·∇)u + ν∇²ω
 *
 * Terms:
 * - Dω/Dt: Material derivative of vorticity
 * - (ω·∇)u: Vortex stretching/tilting
 * - ν∇²ω: Viscous diffusion
 */
class VorticityTransportEquation {
public:
    /**
     * @brief Vortex stretching term
     *
     * (ω·∇)u = [ωⱼ ∂uᵢ/∂xⱼ]
     *
     * In 3D: Stretches and tilts vortex lines
     * In 2D: Zero (no vortex stretching in 2D)
     *
     * @param omega Vorticity vector ω (1/s)
     * @param velocity_gradient ∂uᵢ/∂xⱼ
     * @return Vortex stretching (1/s²)
     */
    static Eigen::Vector3d vortexStretching(
        const Eigen::Vector3d& omega,
        const Eigen::Matrix3d& velocity_gradient) {

        return velocity_gradient.transpose() * omega;
    }

    /**
     * @brief Viscous diffusion term
     *
     * ν∇²ω
     *
     * @param vorticity_laplacian ∇²ω (1/s·m²)
     * @param kinematic_viscosity ν (m²/s)
     * @return Viscous diffusion (1/s²)
     */
    static Eigen::Vector3d viscousDiffusion(
        const Eigen::Vector3d& vorticity_laplacian,
        double kinematic_viscosity) {

        return kinematic_viscosity * vorticity_laplacian;
    }

    /**
     * @brief Total vorticity evolution
     *
     * Dω/Dt = (ω·∇)u + ν∇²ω
     *
     * @param omega Current vorticity (1/s)
     * @param velocity_gradient ∂uᵢ/∂xⱼ
     * @param vorticity_laplacian ∇²ω (1/s·m²)
     * @param kinematic_viscosity ν (m²/s)
     * @return Dω/Dt (1/s²)
     */
    static Eigen::Vector3d totalEvolution(
        const Eigen::Vector3d& omega,
        const Eigen::Matrix3d& velocity_gradient,
        const Eigen::Vector3d& vorticity_laplacian,
        double kinematic_viscosity) {

        Eigen::Vector3d stretching = vortexStretching(omega, velocity_gradient);
        Eigen::Vector3d diffusion = viscousDiffusion(vorticity_laplacian,
                                                     kinematic_viscosity);

        return stretching + diffusion;
    }

    /**
     * @brief 2D vorticity transport equation
     *
     * In 2D: Dω/Dt = ν∇²ω (no vortex stretching)
     *
     * @param vorticity_laplacian ∇²ω (1/s·m²)
     * @param kinematic_viscosity ν (m²/s)
     * @return Dω/Dt (1/s²)
     */
    static double evolution2D(double vorticity_laplacian,
                              double kinematic_viscosity) {
        return kinematic_viscosity * vorticity_laplacian;
    }

    /**
     * @brief Vortex stretching magnitude
     *
     * Sω = ω·(ω·∇)u / |ω|²
     *
     * Positive: vortex intensification
     * Negative: vortex weakening
     *
     * @param omega Vorticity vector (1/s)
     * @param velocity_gradient ∂uᵢ/∂xⱼ
     * @return Stretching rate (1/s²)
     */
    static double stretchingRate(const Eigen::Vector3d& omega,
                                 const Eigen::Matrix3d& velocity_gradient) {
        double omega_mag_sq = omega.squaredNorm();
        if (omega_mag_sq < 1e-10) {
            return 0.0;
        }

        Eigen::Vector3d stretching = vortexStretching(omega, velocity_gradient);
        return omega.dot(stretching) / omega_mag_sq;
    }
};

/**
 * @class Circulation
 * @brief Circulation Γ = ∮ u·dl
 *
 * Line integral of velocity around closed curve
 */
class Circulation {
public:
    /**
     * @brief Compute circulation around closed path
     *
     * Γ = ∮_C u·dl
     *
     * By Stokes' theorem: Γ = ∬_S ω·n dS
     *
     * @param velocity_field Function u(position)
     * @param path_points Points defining closed curve
     * @return Circulation (m²/s)
     */
    static double compute(
        const std::function<Eigen::Vector3d(const Eigen::Vector3d&)>& velocity_field,
        const std::vector<Eigen::Vector3d>& path_points) {

        if (path_points.size() < 2) {
            throw std::invalid_argument("Need at least 2 points for circulation");
        }

        double gamma = 0.0;

        for (size_t i = 0; i < path_points.size(); ++i) {
            size_t j = (i + 1) % path_points.size();

            Eigen::Vector3d r_i = path_points[i];
            Eigen::Vector3d r_j = path_points[j];
            Eigen::Vector3d dl = r_j - r_i;

            // Use midpoint velocity
            Eigen::Vector3d r_mid = 0.5 * (r_i + r_j);
            Eigen::Vector3d u_mid = velocity_field(r_mid);

            gamma += u_mid.dot(dl);
        }

        return gamma;
    }

    /**
     * @brief Circulation from vorticity (Stokes' theorem)
     *
     * Γ = ∬_S ω·n dS
     *
     * @param omega_field Vorticity function ω(position)
     * @param surface_points Points defining surface S
     * @param normal Surface normal n
     * @return Circulation (m²/s)
     */
    static double fromVorticity(
        const std::function<Eigen::Vector3d(const Eigen::Vector3d&)>& omega_field,
        const std::vector<Eigen::Vector3d>& surface_points,
        const Eigen::Vector3d& normal) {

        double gamma = 0.0;
        double area_element = 0.0;  // Simplified: need proper surface integration

        for (const auto& point : surface_points) {
            Eigen::Vector3d omega = omega_field(point);
            gamma += omega.dot(normal);
        }

        // Simple rectangular integration (area / N points)
        area_element = 1.0 / surface_points.size();
        gamma *= area_element;

        return gamma;
    }

    /**
     * @brief Circulation for solid-body rotation
     *
     * For solid-body rotation ω = constant:
     * Γ = ωA
     *
     * where A is enclosed area
     *
     * @param omega_magnitude Rotation rate (1/s)
     * @param enclosed_area A (m²)
     * @return Circulation (m²/s)
     */
    static double solidBodyRotation(double omega_magnitude,
                                    double enclosed_area) {
        return omega_magnitude * enclosed_area;
    }

    /**
     * @brief Circulation for point vortex
     *
     * Γ = 2πr²ω for radius r
     *
     * @param vortex_strength Strength of point vortex (m²/s)
     * @return Circulation (m²/s)
     */
    static double pointVortex(double vortex_strength) {
        return vortex_strength;
    }
};

/**
 * @class KelvinCirculationTheorem
 * @brief Kelvin's circulation theorem (inviscid, barotropic)
 *
 * DΓ/Dt = 0 (circulation conserved for material loop)
 */
class KelvinCirculationTheorem {
public:
    /**
     * @brief Check if conditions for Kelvin's theorem hold
     *
     * Requirements:
     * 1. Inviscid flow (ν = 0)
     * 2. Barotropic (ρ = ρ(p) only)
     * 3. Conservative body forces (F = -∇Φ)
     *
     * @param kinematic_viscosity ν (m²/s)
     * @param tolerance Threshold for inviscid
     * @return true if conditions satisfied
     */
    static bool checkApplicability(double kinematic_viscosity,
                                   double tolerance = 1e-10) {
        return kinematic_viscosity < tolerance;
    }

    /**
     * @brief Rate of change of circulation (viscous case)
     *
     * DΓ/Dt = ∮_C ν(∇×ω)·dl + baroclinic terms
     *
     * For inviscid barotropic: DΓ/Dt = 0
     *
     * @param kinematic_viscosity ν (m²/s)
     * @param vorticity_curl ∇×ω around curve
     * @return DΓ/Dt (m²/s²)
     */
    static double circulationEvolution(
        double kinematic_viscosity,
        const std::vector<Eigen::Vector3d>& vorticity_curl,
        const std::vector<Eigen::Vector3d>& path_elements) {

        if (vorticity_curl.size() != path_elements.size()) {
            throw std::invalid_argument("Curl and path arrays must match");
        }

        double dGamma_dt = 0.0;

        for (size_t i = 0; i < vorticity_curl.size(); ++i) {
            dGamma_dt += kinematic_viscosity *
                         vorticity_curl[i].dot(path_elements[i]);
        }

        return dGamma_dt;
    }

    /**
     * @brief Verify circulation conservation
     *
     * Check |Γ(t) - Γ(0)| < tolerance
     *
     * @param circulation_history Time series of Γ
     * @param tolerance Error threshold
     * @return true if conserved
     */
    static bool verifyConservation(
        const std::vector<double>& circulation_history,
        double tolerance = 1e-6) {

        if (circulation_history.empty()) return true;

        double gamma_0 = circulation_history[0];

        for (double gamma : circulation_history) {
            if (std::abs(gamma - gamma_0) > tolerance) {
                return false;
            }
        }

        return true;
    }
};

/**
 * @class BiotSavartLaw
 * @brief Biot-Savart law for velocity induced by vorticity
 *
 * u(r) = (1/4π) ∫ [ω(r') × (r - r')] / |r - r'|³ dV'
 */
class BiotSavartLaw {
public:
    /**
     * @brief Velocity induced by vortex filament
     *
     * For infinitesimal vortex element:
     * du = (Γ/4π) [dl × (r - r')] / |r - r'|³
     *
     * @param circulation Γ (m²/s)
     * @param filament_element dl (m)
     * @param filament_position r' (m)
     * @param field_point r (m)
     * @return Induced velocity (m/s)
     */
    static Eigen::Vector3d velocityFromFilament(
        double circulation,
        const Eigen::Vector3d& filament_element,
        const Eigen::Vector3d& filament_position,
        const Eigen::Vector3d& field_point) {

        Eigen::Vector3d r = field_point - filament_position;
        double r_mag = r.norm();

        if (r_mag < 1e-10) {
            return Eigen::Vector3d::Zero();  // Singularity at vortex core
        }

        Eigen::Vector3d cross = filament_element.cross(r);
        double r_cubed = r_mag * r_mag * r_mag;

        return (circulation / (4.0 * M_PI)) * cross / r_cubed;
    }

    /**
     * @brief Velocity induced by straight vortex filament
     *
     * Infinite straight vortex along z-axis:
     * u_θ = Γ/(2πr)
     *
     * @param circulation Γ (m²/s)
     * @param radial_distance r from vortex axis (m)
     * @return Tangential velocity (m/s)
     */
    static double velocityStraightVortex(double circulation,
                                         double radial_distance) {
        if (radial_distance < 1e-10) {
            throw std::invalid_argument("Too close to vortex core");
        }
        return circulation / (2.0 * M_PI * radial_distance);
    }

    /**
     * @brief Velocity induced by circular vortex ring
     *
     * On axis of ring at distance z:
     * u_z = (ΓR²) / [2(R² + z²)^(3/2)]
     *
     * @param circulation Γ (m²/s)
     * @param ring_radius R (m)
     * @param axial_distance z (m)
     * @return Axial velocity (m/s)
     */
    static double velocityVortexRing(double circulation,
                                     double ring_radius,
                                     double axial_distance) {
        double R2 = ring_radius * ring_radius;
        double z2 = axial_distance * axial_distance;
        double denom = std::pow(R2 + z2, 1.5);

        return (circulation * R2) / (2.0 * denom);
    }

    /**
     * @brief Self-induced velocity of vortex ring
     *
     * Ring translates with velocity:
     * U = (Γ/4πR)[ln(8R/a) - 1/4]
     *
     * where a is core radius
     *
     * @param circulation Γ (m²/s)
     * @param ring_radius R (m)
     * @param core_radius a (m)
     * @return Translation velocity (m/s)
     */
    static double vortexRingSelfVelocity(double circulation,
                                         double ring_radius,
                                         double core_radius) {
        if (core_radius >= ring_radius) {
            throw std::invalid_argument("Core radius must be less than ring radius");
        }

        double ratio = 8.0 * ring_radius / core_radius;
        return (circulation / (4.0 * M_PI * ring_radius)) *
               (std::log(ratio) - 0.25);
    }

    /**
     * @brief Velocity from vortex sheet
     *
     * Discontinuity in tangential velocity:
     * Δu = γ
     *
     * where γ is sheet strength (circulation per unit length)
     *
     * @param sheet_strength γ (m/s)
     * @return Velocity jump (m/s)
     */
    static double velocityJumpAcrossSheet(double sheet_strength) {
        return sheet_strength;
    }
};

/**
 * @class VortexDynamics
 * @brief Dynamics of vortex structures
 */
class VortexDynamics {
public:
    /**
     * @brief Interaction of two point vortices
     *
     * Vortices orbit around their centroid
     * ω = (Γ₁ + Γ₂)/(2πd²)
     *
     * @param circulation1 Γ₁ (m²/s)
     * @param circulation2 Γ₂ (m²/s)
     * @param separation d (m)
     * @return Angular velocity (1/s)
     */
    static double vortexPairRotation(double circulation1,
                                     double circulation2,
                                     double separation) {
        if (separation < 1e-10) {
            throw std::invalid_argument("Vortices too close");
        }
        return (circulation1 + circulation2) /
               (2.0 * M_PI * separation * separation);
    }

    /**
     * @brief Vortex pair translation velocity
     *
     * Two vortices with opposite circulation:
     * U = Γ/(2πd)
     *
     * @param circulation Γ (m²/s)
     * @param separation d (m)
     * @return Translation speed (m/s)
     */
    static double vortexPairVelocity(double circulation,
                                     double separation) {
        if (separation < 1e-10) {
            throw std::invalid_argument("Vortices too close");
        }
        return circulation / (2.0 * M_PI * separation);
    }

    /**
     * @brief Rankine vortex velocity profile
     *
     * Solid body rotation inside core:
     * u = ωr  for r < a
     *
     * Potential flow outside:
     * u = Γ/(2πr)  for r > a
     *
     * @param radius r (m)
     * @param core_radius a (m)
     * @param circulation Γ (m²/s)
     * @return Tangential velocity (m/s)
     */
    static double rankineVortex(double radius,
                                double core_radius,
                                double circulation) {
        if (radius < core_radius) {
            // Solid body rotation
            double omega = circulation / (M_PI * core_radius * core_radius);
            return omega * radius;
        } else {
            // Potential flow
            return circulation / (2.0 * M_PI * radius);
        }
    }
};

} // namespace physics::advanced::fluid_dynamics

#endif // PHYSICS_ADVANCED_FLUID_DYNAMICS_VORTICITY_HPP
