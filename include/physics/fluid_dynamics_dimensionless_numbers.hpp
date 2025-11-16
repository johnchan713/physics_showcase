#ifndef PHYSICS_ADVANCED_FLUID_DYNAMICS_DIMENSIONLESS_NUMBERS_HPP
#define PHYSICS_ADVANCED_FLUID_DYNAMICS_DIMENSIONLESS_NUMBERS_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file dimensionless_numbers.hpp
 * @brief Dimensionless numbers for fluid dynamics similarity
 *
 * These numbers characterize the relative importance of different
 * physical effects in fluid flow.
 */

namespace physics::advanced::fluid_dynamics {

/**
 * @class DimensionlessNumbers
 * @brief Collection of important dimensionless parameters
 */
class DimensionlessNumbers {
public:
    /**
     * @brief Reynolds Number - Inertia/Viscous forces ratio
     *
     * Re = ρUL/μ = UL/ν
     *
     * Significance:
     * - Re << 1: Viscous (Stokes) flow
     * - Re ~ 1: Transition regime
     * - Re >> 1: Inviscid flow
     * - Re > 2300 (pipe): Turbulent flow
     *
     * @param density ρ (kg/m³)
     * @param velocity U (m/s)
     * @param length_scale L (m)
     * @param dynamic_viscosity μ (Pa·s)
     * @return Reynolds number (dimensionless)
     */
    static double reynoldsNumber(double density, double velocity,
                                double length_scale,
                                double dynamic_viscosity) {
        if (dynamic_viscosity <= 0.0) {
            throw std::invalid_argument("Viscosity must be positive");
        }
        return (density * velocity * length_scale) / dynamic_viscosity;
    }

    /**
     * @brief Reynolds Number from kinematic viscosity
     *
     * Re = UL/ν
     */
    static double reynoldsNumberKinematic(double velocity,
                                         double length_scale,
                                         double kinematic_viscosity) {
        if (kinematic_viscosity <= 0.0) {
            throw std::invalid_argument("Kinematic viscosity must be positive");
        }
        return (velocity * length_scale) / kinematic_viscosity;
    }

    /**
     * @brief Froude Number - Inertia/Gravity forces ratio
     *
     * Fr = U/√(gL)
     *
     * Significance:
     * - Fr < 1: Subcritical flow (gravity dominant)
     * - Fr = 1: Critical flow
     * - Fr > 1: Supercritical flow (inertia dominant)
     *
     * Important for free surface flows, hydraulic jumps
     *
     * @param velocity U (m/s)
     * @param gravity g (m/s²)
     * @param length_scale L (m)
     * @return Froude number
     */
    static double froudeNumber(double velocity, double gravity,
                              double length_scale) {
        if (gravity <= 0.0 || length_scale <= 0.0) {
            throw std::invalid_argument("Gravity and length must be positive");
        }
        return velocity / std::sqrt(gravity * length_scale);
    }

    /**
     * @brief Mach Number - Velocity/Sound speed ratio
     *
     * Ma = U/c
     *
     * Significance:
     * - Ma < 0.3: Incompressible flow
     * - 0.3 < Ma < 0.8: Subsonic compressible
     * - Ma ≈ 1: Transonic
     * - 1 < Ma < 5: Supersonic
     * - Ma > 5: Hypersonic
     *
     * @param velocity U (m/s)
     * @param sound_speed c (m/s)
     * @return Mach number
     */
    static double machNumber(double velocity, double sound_speed) {
        if (sound_speed <= 0.0) {
            throw std::invalid_argument("Sound speed must be positive");
        }
        return velocity / sound_speed;
    }

    /**
     * @brief Sound speed in ideal gas
     *
     * c = √(γRT)
     *
     * @param gamma Specific heat ratio γ = cp/cv
     * @param gas_constant R (J/(kg·K))
     * @param temperature T (K)
     * @return Sound speed (m/s)
     */
    static double soundSpeed(double gamma, double gas_constant,
                            double temperature) {
        return std::sqrt(gamma * gas_constant * temperature);
    }

    /**
     * @brief Prandtl Number - Momentum/Thermal diffusivity ratio
     *
     * Pr = ν/α = μcp/k
     *
     * Significance:
     * - Pr << 1: Thermal diffusion dominates (liquid metals)
     * - Pr ~ 1: Momentum and thermal diffusion comparable (gases)
     * - Pr >> 1: Momentum diffusion dominates (oils)
     *
     * Typical values:
     * - Air: ~0.7
     * - Water: ~7
     * - Oil: ~100-1000
     * - Liquid metals: ~0.01
     *
     * @param kinematic_viscosity ν (m²/s)
     * @param thermal_diffusivity α (m²/s)
     * @return Prandtl number
     */
    static double prandtlNumber(double kinematic_viscosity,
                               double thermal_diffusivity) {
        if (thermal_diffusivity <= 0.0) {
            throw std::invalid_argument("Thermal diffusivity must be positive");
        }
        return kinematic_viscosity / thermal_diffusivity;
    }

    /**
     * @brief Prandtl Number from physical properties
     *
     * Pr = μcp/k
     */
    static double prandtlNumberFromProperties(
        double dynamic_viscosity,
        double specific_heat,
        double thermal_conductivity) {

        if (thermal_conductivity <= 0.0) {
            throw std::invalid_argument("Thermal conductivity must be positive");
        }
        return (dynamic_viscosity * specific_heat) / thermal_conductivity;
    }

    /**
     * @brief Grashof Number - Buoyancy/Viscous forces ratio
     *
     * Gr = gβ(Ts - T∞)L³/ν²
     *
     * Significance:
     * - Gr < 10⁴: Laminar natural convection
     * - Gr > 10⁹: Turbulent natural convection
     *
     * @param gravity g (m/s²)
     * @param thermal_expansion β (1/K)
     * @param temp_diff ΔT = Ts - T∞ (K)
     * @param length_scale L (m)
     * @param kinematic_viscosity ν (m²/s)
     * @return Grashof number
     */
    static double grashofNumber(double gravity,
                               double thermal_expansion,
                               double temp_diff,
                               double length_scale,
                               double kinematic_viscosity) {
        double nu_sq = kinematic_viscosity * kinematic_viscosity;
        if (nu_sq <= 0.0) {
            throw std::invalid_argument("Kinematic viscosity must be positive");
        }

        double L_cubed = length_scale * length_scale * length_scale;
        return (gravity * thermal_expansion * temp_diff * L_cubed) / nu_sq;
    }

    /**
     * @brief Rayleigh Number - Natural convection parameter
     *
     * Ra = Gr × Pr = gβΔTL³/(να)
     *
     * Significance:
     * - Ra < 10³: Conduction dominant
     * - Ra > 10³: Convection begins
     * - Ra > 10⁹: Turbulent natural convection
     *
     * @param grashof Gr
     * @param prandtl Pr
     * @return Rayleigh number
     */
    static double rayleighNumber(double grashof, double prandtl) {
        return grashof * prandtl;
    }

    /**
     * @brief Rayleigh Number from properties
     */
    static double rayleighNumberFromProperties(
        double gravity,
        double thermal_expansion,
        double temp_diff,
        double length_scale,
        double kinematic_viscosity,
        double thermal_diffusivity) {

        double Gr = grashofNumber(gravity, thermal_expansion, temp_diff,
                                 length_scale, kinematic_viscosity);
        double Pr = prandtlNumber(kinematic_viscosity, thermal_diffusivity);
        return Gr * Pr;
    }

    /**
     * @brief Nusselt Number - Convective/Conductive heat transfer ratio
     *
     * Nu = hL/k
     *
     * where h is convective heat transfer coefficient
     *
     * Significance:
     * - Nu = 1: Pure conduction
     * - Nu > 1: Convection enhances heat transfer
     *
     * Correlations:
     * - Laminar flat plate: Nu ~ Re^0.5 × Pr^0.33
     * - Turbulent flat plate: Nu ~ Re^0.8 × Pr^0.33
     *
     * @param heat_transfer_coeff h (W/(m²·K))
     * @param length_scale L (m)
     * @param thermal_conductivity k (W/(m·K))
     * @return Nusselt number
     */
    static double nusseltNumber(double heat_transfer_coeff,
                               double length_scale,
                               double thermal_conductivity) {
        if (thermal_conductivity <= 0.0) {
            throw std::invalid_argument("Thermal conductivity must be positive");
        }
        return (heat_transfer_coeff * length_scale) / thermal_conductivity;
    }

    /**
     * @brief Nusselt Number correlation for laminar flat plate
     *
     * Nu = 0.664 × Re^0.5 × Pr^0.33
     */
    static double nusseltLaminarFlatPlate(double reynolds, double prandtl) {
        return 0.664 * std::pow(reynolds, 0.5) * std::pow(prandtl, 1.0/3.0);
    }

    /**
     * @brief Nusselt Number correlation for turbulent flat plate
     *
     * Nu = 0.037 × Re^0.8 × Pr^0.33
     */
    static double nusseltTurbulentFlatPlate(double reynolds, double prandtl) {
        return 0.037 * std::pow(reynolds, 0.8) * std::pow(prandtl, 1.0/3.0);
    }

    /**
     * @brief Peclet Number - Advection/Diffusion ratio
     *
     * Pe = UL/α = Re × Pr (thermal)
     * Pe = UL/D (mass transfer)
     *
     * Significance:
     * - Pe << 1: Diffusion dominant
     * - Pe >> 1: Advection dominant
     *
     * @param velocity U (m/s)
     * @param length_scale L (m)
     * @param diffusivity α or D (m²/s)
     * @return Peclet number
     */
    static double pecletNumber(double velocity, double length_scale,
                              double diffusivity) {
        if (diffusivity <= 0.0) {
            throw std::invalid_argument("Diffusivity must be positive");
        }
        return (velocity * length_scale) / diffusivity;
    }

    /**
     * @brief Peclet Number from Reynolds and Prandtl
     *
     * Pe = Re × Pr (for heat transfer)
     */
    static double pecletFromReynoldsPrandtl(double reynolds, double prandtl) {
        return reynolds * prandtl;
    }

    /**
     * @brief Schmidt Number - Momentum/Mass diffusivity ratio
     *
     * Sc = ν/D
     *
     * Analogous to Prandtl number for mass transfer
     *
     * @param kinematic_viscosity ν (m²/s)
     * @param mass_diffusivity D (m²/s)
     * @return Schmidt number
     */
    static double schmidtNumber(double kinematic_viscosity,
                                double mass_diffusivity) {
        if (mass_diffusivity <= 0.0) {
            throw std::invalid_argument("Mass diffusivity must be positive");
        }
        return kinematic_viscosity / mass_diffusivity;
    }

    /**
     * @brief Sherwood Number - Mass convection/diffusion ratio
     *
     * Sh = kL/D
     *
     * Analogous to Nusselt number for mass transfer
     *
     * @param mass_transfer_coeff k (m/s)
     * @param length_scale L (m)
     * @param mass_diffusivity D (m²/s)
     * @return Sherwood number
     */
    static double sherwoodNumber(double mass_transfer_coeff,
                                 double length_scale,
                                 double mass_diffusivity) {
        if (mass_diffusivity <= 0.0) {
            throw std::invalid_argument("Mass diffusivity must be positive");
        }
        return (mass_transfer_coeff * length_scale) / mass_diffusivity;
    }

    /**
     * @brief Weber Number - Inertia/Surface tension ratio
     *
     * We = ρU²L/σ
     *
     * Significance:
     * - We << 1: Surface tension dominant (droplets remain spherical)
     * - We >> 1: Inertia dominant (droplet breakup)
     *
     * @param density ρ (kg/m³)
     * @param velocity U (m/s)
     * @param length_scale L (m)
     * @param surface_tension σ (N/m)
     * @return Weber number
     */
    static double weberNumber(double density, double velocity,
                             double length_scale, double surface_tension) {
        if (surface_tension <= 0.0) {
            throw std::invalid_argument("Surface tension must be positive");
        }
        return (density * velocity * velocity * length_scale) / surface_tension;
    }

    /**
     * @brief Capillary Number - Viscous/Surface tension ratio
     *
     * Ca = μU/σ
     *
     * @param dynamic_viscosity μ (Pa·s)
     * @param velocity U (m/s)
     * @param surface_tension σ (N/m)
     * @return Capillary number
     */
    static double capillaryNumber(double dynamic_viscosity,
                                  double velocity,
                                  double surface_tension) {
        if (surface_tension <= 0.0) {
            throw std::invalid_argument("Surface tension must be positive");
        }
        return (dynamic_viscosity * velocity) / surface_tension;
    }

    /**
     * @brief Strouhal Number - Unsteady/Convective effects ratio
     *
     * St = fL/U
     *
     * where f is frequency of vortex shedding
     *
     * @param frequency f (Hz)
     * @param length_scale L (m)
     * @param velocity U (m/s)
     * @return Strouhal number
     */
    static double strouhalNumber(double frequency,
                                 double length_scale,
                                 double velocity) {
        if (velocity <= 0.0) {
            throw std::invalid_argument("Velocity must be positive");
        }
        return (frequency * length_scale) / velocity;
    }

    /**
     * @brief Classify flow regime based on Reynolds number
     */
    enum class FlowRegime {
        STOKES,        // Re < 1
        LAMINAR,       // 1 < Re < 2300 (pipe)
        TRANSITION,    // 2300 < Re < 4000
        TURBULENT      // Re > 4000
    };

    static FlowRegime classifyFlowRegime(double reynolds) {
        if (reynolds < 1.0) return FlowRegime::STOKES;
        if (reynolds < 2300.0) return FlowRegime::LAMINAR;
        if (reynolds < 4000.0) return FlowRegime::TRANSITION;
        return FlowRegime::TURBULENT;
    }

    /**
     * @brief Classify compressibility based on Mach number
     */
    enum class CompressibilityRegime {
        INCOMPRESSIBLE,  // Ma < 0.3
        SUBSONIC,        // 0.3 < Ma < 0.8
        TRANSONIC,       // 0.8 < Ma < 1.2
        SUPERSONIC,      // 1.2 < Ma < 5
        HYPERSONIC       // Ma > 5
    };

    static CompressibilityRegime classifyCompressibility(double mach) {
        if (mach < 0.3) return CompressibilityRegime::INCOMPRESSIBLE;
        if (mach < 0.8) return CompressibilityRegime::SUBSONIC;
        if (mach < 1.2) return CompressibilityRegime::TRANSONIC;
        if (mach < 5.0) return CompressibilityRegime::SUPERSONIC;
        return CompressibilityRegime::HYPERSONIC;
    }
};

} // namespace physics::advanced::fluid_dynamics

#endif // PHYSICS_ADVANCED_FLUID_DYNAMICS_DIMENSIONLESS_NUMBERS_HPP
