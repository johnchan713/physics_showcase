#ifndef PHYSICS_ADVANCED_COSMOLOGY_ENERGY_DENSITY_HPP
#define PHYSICS_ADVANCED_COSMOLOGY_ENERGY_DENSITY_HPP

#include <cmath>
#include <string>
#include <map>

/**
 * @file energy_density.hpp
 * @brief Energy density components of the universe
 *
 * Implements:
 * - Matter (baryonic + dark matter)
 * - Radiation (photons + neutrinos)
 * - Dark energy (cosmological constant)
 * - Observed density parameters
 * - Equation of state for each component
 */

namespace physics::advanced::cosmology {

/**
 * @class EnergyDensityComponents
 * @brief Components of cosmic energy density
 *
 * ρ_total = ρ_matter + ρ_radiation + ρ_dark_energy
 */
class EnergyDensityComponents {
public:
    /**
     * @brief Matter density parameter (today)
     *
     * Ω_m = ρ_m / ρ_crit ≈ 0.315 ± 0.007 (Planck 2018)
     *
     * Includes both baryonic and dark matter
     */
    static double matterDensityParameter() {
        return 0.315;
    }

    /**
     * @brief Baryonic matter density parameter
     *
     * Ω_b = ρ_baryon / ρ_crit ≈ 0.049 ± 0.001
     *
     * Ordinary matter (protons, neutrons, electrons)
     */
    static double baryonicDensityParameter() {
        return 0.049;
    }

    /**
     * @brief Dark matter density parameter
     *
     * Ω_DM = Ω_m - Ω_b ≈ 0.266
     *
     * Non-baryonic dark matter
     */
    static double darkMatterDensityParameter() {
        return matterDensityParameter() - baryonicDensityParameter();  // ~0.266
    }

    /**
     * @brief Radiation density parameter (today)
     *
     * Ω_r ≈ 9.24×10⁻⁵ (photons + neutrinos)
     *
     * Negligible today, but dominated early universe
     */
    static double radiationDensityParameter() {
        return 9.24e-5;
    }

    /**
     * @brief Photon density parameter
     *
     * Ω_γ ≈ 5.4×10⁻⁵
     *
     * CMB photons
     */
    static double photonDensityParameter() {
        return 5.4e-5;
    }

    /**
     * @brief Neutrino density parameter
     *
     * Ω_ν ≈ 3.8×10⁻⁵
     *
     * Cosmic neutrino background (CνB)
     */
    static double neutrinoDensityParameter() {
        return 3.8e-5;
    }

    /**
     * @brief Dark energy density parameter
     *
     * Ω_Λ ≈ 0.685 ± 0.007 (Planck 2018)
     *
     * Vacuum energy / cosmological constant
     */
    static double darkEnergyDensityParameter() {
        return 0.685;
    }

    /**
     * @brief Total density parameter
     *
     * Ω_total = Ω_m + Ω_r + Ω_Λ + Ω_k
     *
     * For flat universe: Ω_total = 1
     */
    static double totalDensityParameter() {
        return matterDensityParameter() +
               radiationDensityParameter() +
               darkEnergyDensityParameter();  // ≈ 1.0
    }

    /**
     * @brief Curvature density parameter
     *
     * Ω_k = 1 - Ω_total ≈ 0.001 ± 0.002
     */
    static double curvatureDensityParameter() {
        return 1.0 - totalDensityParameter();
    }
};

/**
 * @class MatterComponent
 * @brief Matter energy density (w = 0)
 *
 * Pressure-less: p = 0
 * Scales as: ρ_m ∝ a⁻³
 */
class MatterComponent {
public:
    /**
     * @brief Equation of state parameter
     *
     * w_m = p/(ρc²) = 0 (pressure-less dust)
     */
    static double equationOfState() {
        return 0.0;
    }

    /**
     * @brief Matter density at scale factor a
     *
     * ρ_m(a) = ρ_m,0 a⁻³
     *
     * Dilutes as volume ∝ a³
     */
    static double densityAtScaleFactor(double rho_m0, double scale_factor) {
        return rho_m0 * std::pow(scale_factor, -3.0);
    }

    /**
     * @brief Matter density today
     *
     * ρ_m,0 = Ω_m × ρ_crit ≈ 2.7×10⁻²⁷ kg/m³
     */
    static double densityToday() {
        double Omega_m = EnergyDensityComponents::matterDensityParameter();
        double rho_crit = 8.5e-27;  // kg/m³
        return Omega_m * rho_crit;
    }

    /**
     * @brief Number density of baryons today
     *
     * n_b ≈ Ω_b ρ_crit / m_p ≈ 0.25 protons/m³
     */
    static double baryonNumberDensity() {
        double Omega_b = EnergyDensityComponents::baryonicDensityParameter();
        double rho_crit = 8.5e-27;  // kg/m³
        double m_proton = 1.673e-27;  // kg

        return (Omega_b * rho_crit) / m_proton;  // m⁻³
    }
};

/**
 * @class RadiationComponent
 * @brief Radiation energy density (w = 1/3)
 *
 * Relativistic particles
 * Scales as: ρ_r ∝ a⁻⁴
 */
class RadiationComponent {
public:
    /**
     * @brief Equation of state parameter
     *
     * w_r = p/(ρc²) = 1/3 (relativistic)
     */
    static double equationOfState() {
        return 1.0/3.0;
    }

    /**
     * @brief Radiation density at scale factor a
     *
     * ρ_r(a) = ρ_r,0 a⁻⁴
     *
     * Dilutes as a⁻³ (volume) × a⁻¹ (redshift)
     */
    static double densityAtScaleFactor(double rho_r0, double scale_factor) {
        return rho_r0 * std::pow(scale_factor, -4.0);
    }

    /**
     * @brief Radiation density today
     *
     * ρ_r,0 = Ω_r × ρ_crit ≈ 7.8×10⁻³¹ kg/m³
     *
     * Dominated by CMB photons
     */
    static double densityToday() {
        double Omega_r = EnergyDensityComponents::radiationDensityParameter();
        double rho_crit = 8.5e-27;  // kg/m³
        return Omega_r * rho_crit;
    }

    /**
     * @brief CMB photon number density
     *
     * n_γ ≈ 411 photons/cm³ = 4.11×10⁸ m⁻³
     */
    static double photonNumberDensity() {
        return 4.11e8;  // m⁻³
    }

    /**
     * @brief Cosmic neutrino background number density
     *
     * n_ν ≈ 339 neutrinos/cm³ per flavor
     * Total: 3 × 339 ≈ 1017 neutrinos/cm³
     */
    static double neutrinoNumberDensity() {
        return 1.017e9;  // m⁻³ (all 3 flavors)
    }

    /**
     * @brief Photon-to-baryon ratio
     *
     * η = n_γ / n_b ≈ 1.6×10⁹
     *
     * Enormous asymmetry!
     */
    static double photonToBaryonRatio() {
        double n_gamma = photonNumberDensity();
        double n_baryon = MatterComponent::baryonNumberDensity();
        return n_gamma / n_baryon;  // ~1.6×10⁹
    }
};

/**
 * @class DarkEnergyComponent
 * @brief Dark energy / cosmological constant (w = -1)
 *
 * Vacuum energy with negative pressure
 * Constant density: ρ_Λ = constant
 */
class DarkEnergyComponent {
public:
    /**
     * @brief Equation of state parameter
     *
     * w_Λ = p/(ρc²) = -1 (vacuum energy)
     *
     * Negative pressure causes acceleration
     */
    static double equationOfState() {
        return -1.0;
    }

    /**
     * @brief Dark energy density (constant!)
     *
     * ρ_Λ(a) = ρ_Λ,0 = constant
     *
     * Does not dilute with expansion
     */
    static double densityAtScaleFactor(double rho_Lambda0, double scale_factor) {
        return rho_Lambda0;  // Constant!
    }

    /**
     * @brief Dark energy density today
     *
     * ρ_Λ,0 = Ω_Λ × ρ_crit ≈ 5.8×10⁻²⁷ kg/m³
     */
    static double densityToday() {
        double Omega_Lambda = EnergyDensityComponents::darkEnergyDensityParameter();
        double rho_crit = 8.5e-27;  // kg/m³
        return Omega_Lambda * rho_crit;
    }

    /**
     * @brief Cosmological constant Λ
     *
     * Λ = 8πG ρ_Λ / c²
     *
     * Λ ≈ 1.1×10⁻⁵² m⁻²
     */
    static double cosmologicalConstant() {
        double rho_Lambda = densityToday();  // kg/m³
        double G = 6.674e-11;  // m³/(kg·s²)
        double c = 2.998e8;   // m/s

        return (8.0 * M_PI * G * rho_Lambda) / (c * c);  // m⁻²
    }

    /**
     * @brief Vacuum energy density
     *
     * ρ_Λ c² ≈ 5.2×10⁻¹⁰ J/m³ ≈ 0.003 eV/cm³
     *
     * Incredibly small!
     */
    static double vacuumEnergyDensity() {
        double rho_Lambda = densityToday();
        double c = 2.998e8;
        return rho_Lambda * c * c;  // J/m³
    }

    /**
     * @brief Cosmological constant problem
     *
     * Quantum field theory predicts ρ_vacuum ~ M_Planck⁴ ~ 10¹¹⁴ J/m³
     * Observed: ρ_Λ c² ~ 10⁻¹⁰ J/m³
     *
     * Discrepancy: ~10¹²⁴ (worst prediction in physics!)
     */
    static double cosmologicalConstantProblem() {
        return 1e124;  // Factor of discrepancy
    }
};

/**
 * @class EquationOfState
 * @brief Equation of state for cosmic fluids
 *
 * p = w ρ c²
 */
class EquationOfState {
public:
    /**
     * @brief Equation of state values
     */
    static std::map<std::string, double> parameters() {
        return {
            {"matter", 0.0},        // Dust (pressure-less)
            {"radiation", 1.0/3.0}, // Relativistic
            {"dark_energy", -1.0},  // Vacuum energy
            {"stiff_matter", 1.0},  // Maximal pressure (speed of sound = c)
            {"phantom", -1.5}       // w < -1 (hypothetical)
        };
    }

    /**
     * @brief Pressure from energy density
     *
     * p = w ρ c²
     */
    static double pressure(double energy_density, double w, double c = 2.998e8) {
        return w * energy_density * c * c;
    }

    /**
     * @brief Speed of sound squared
     *
     * c_s² = dp/dρ = w c²
     *
     * For adiabatic fluid
     */
    static double soundSpeedSquared(double w) {
        return w;  // In units of c²
    }

    /**
     * @brief Check if fluid is physically reasonable
     *
     * For stability: -1 ≤ w ≤ 1
     * (except phantom energy w < -1)
     */
    static bool isPhysical(double w) {
        return (w >= -1.0 && w <= 1.0) || w == -1.5;  // Allow phantom
    }
};

/**
 * @class DensityEvolution
 * @brief Evolution of energy densities with redshift
 */
class DensityEvolution {
public:
    /**
     * @brief Total energy density at redshift z
     *
     * ρ(z) = ρ_crit,0 [Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_Λ]
     */
    static double totalDensity(double redshift) {
        double z_plus_1 = 1.0 + redshift;
        double rho_crit_0 = 8.5e-27;  // kg/m³

        double Omega_m = EnergyDensityComponents::matterDensityParameter();
        double Omega_r = EnergyDensityComponents::radiationDensityParameter();
        double Omega_Lambda = EnergyDensityComponents::darkEnergyDensityParameter();

        double matter_term = Omega_m * std::pow(z_plus_1, 3.0);
        double rad_term = Omega_r * std::pow(z_plus_1, 4.0);
        double lambda_term = Omega_Lambda;

        return rho_crit_0 * (matter_term + rad_term + lambda_term);
    }

    /**
     * @brief Matter-radiation equality redshift
     *
     * z_eq where ρ_m(z_eq) = ρ_r(z_eq)
     *
     * Ω_m(1+z_eq)³ = Ω_r(1+z_eq)⁴
     * → 1 + z_eq = Ω_m / Ω_r
     */
    static double matterRadiationEquality() {
        double Omega_m = EnergyDensityComponents::matterDensityParameter();
        double Omega_r = EnergyDensityComponents::radiationDensityParameter();

        return (Omega_m / Omega_r) - 1.0;  // z_eq ≈ 3400
    }

    /**
     * @brief Matter-dark energy equality redshift
     *
     * z_Λ where ρ_m(z_Λ) = ρ_Λ
     *
     * Recent: z_Λ ≈ 0.4 (t ≈ 5 Gyr ago)
     */
    static double matterDarkEnergyEquality() {
        double Omega_m = EnergyDensityComponents::matterDensityParameter();
        double Omega_Lambda = EnergyDensityComponents::darkEnergyDensityParameter();

        return std::pow(Omega_Lambda / Omega_m, 1.0/3.0) - 1.0;  // z_Λ ≈ 0.4
    }

    /**
     * @brief Dominant component at redshift z
     */
    static std::string dominantComponent(double redshift) {
        double z_eq = matterRadiationEquality();  // ~3400
        double z_Lambda = matterDarkEnergyEquality();  // ~0.4

        if (redshift > z_eq) {
            return "Radiation dominated";
        } else if (redshift > z_Lambda) {
            return "Matter dominated";
        } else {
            return "Dark energy dominated";
        }
    }
};

} // namespace physics::advanced::cosmology

#endif // PHYSICS_ADVANCED_COSMOLOGY_ENERGY_DENSITY_HPP
