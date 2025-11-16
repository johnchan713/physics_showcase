#ifndef PHYSICS_ADVANCED_COSMOLOGY_EXPANDING_UNIVERSE_HPP
#define PHYSICS_ADVANCED_COSMOLOGY_EXPANDING_UNIVERSE_HPP

#include <cmath>
#include <vector>
#include <string>

/**
 * @file expanding_universe.hpp
 * @brief The expanding universe: Hubble law and cosmological redshift
 *
 * Implements:
 * - Hubble's law (v = H₀d)
 * - Cosmological redshift
 * - Scale factor a(t)
 * - Olbers' paradox resolution
 * - Observable universe
 */

namespace physics::advanced::cosmology {

/**
 * @class HubbleExpansion
 * @brief Hubble's law and the expanding universe
 *
 * The universe is expanding uniformly: v = H₀d
 * where H₀ is the Hubble constant
 */
class HubbleExpansion {
public:
    /**
     * @brief Hubble constant (present day)
     *
     * H₀ = 67.4 ± 0.5 km/s/Mpc (Planck 2018)
     * H₀ = 73.0 ± 1.0 km/s/Mpc (local measurements - Hubble tension!)
     *
     * Using Planck value as default
     */
    static double hubbleConstant() {
        return 67.4;  // km/s/Mpc
    }

    /**
     * @brief Hubble constant in SI units
     *
     * H₀ in s⁻¹
     */
    static double hubbleConstantSI() {
        double H0_km_s_Mpc = hubbleConstant();
        double Mpc_to_km = 3.086e19;  // 1 Mpc = 3.086×10¹⁹ km

        return H0_km_s_Mpc / Mpc_to_km;  // s⁻¹
    }

    /**
     * @brief Hubble time (age scale)
     *
     * t_H = 1/H₀ ≈ 14.5 Gyr
     *
     * Characteristic timescale of universe
     */
    static double hubbleTime() {
        return 1.0 / hubbleConstantSI();  // seconds
    }

    /**
     * @brief Hubble time in Gyr
     */
    static double hubbleTimeGyr() {
        double t_H_seconds = hubbleTime();
        double seconds_per_year = 3.156e7;
        double Gyr = 1e9;

        return t_H_seconds / (seconds_per_year * Gyr);  // ~14.5 Gyr
    }

    /**
     * @brief Hubble's law: recession velocity
     *
     * v = H₀ × d
     *
     * @param distance d (Mpc)
     * @return Recession velocity (km/s)
     */
    static double recessionVelocity(double distance_Mpc) {
        return hubbleConstant() * distance_Mpc;  // km/s
    }

    /**
     * @brief Distance from recession velocity
     *
     * d = v / H₀
     */
    static double distanceFromVelocity(double velocity_km_s) {
        return velocity_km_s / hubbleConstant();  // Mpc
    }

    /**
     * @brief Hubble sphere (radius where v = c)
     *
     * R_H = c / H₀
     *
     * Beyond this, galaxies recede faster than light (allowed!)
     */
    static double hubbleSphere() {
        double c_km_s = 2.998e5;  // km/s
        return c_km_s / hubbleConstant();  // Mpc ~4200 Mpc
    }

    /**
     * @brief Hubble tension
     *
     * Discrepancy between CMB (67.4) and local (73.0) measurements
     */
    static double hubbleTension() {
        double H0_CMB = 67.4;
        double H0_local = 73.0;
        return std::abs(H0_local - H0_CMB);  // ~5.6 km/s/Mpc
    }
};

/**
 * @class CosmologicalRedshift
 * @brief Redshift due to expansion of space
 *
 * z = (λ_observed - λ_emitted) / λ_emitted = Δλ/λ
 */
class CosmologicalRedshift {
public:
    /**
     * @brief Redshift from wavelength shift
     *
     * z = Δλ/λ = (λ_obs - λ_emit)/λ_emit
     *
     * @param lambda_observed Observed wavelength
     * @param lambda_emitted Emitted wavelength
     * @return Redshift z
     */
    static double fromWavelength(double lambda_observed, double lambda_emitted) {
        return (lambda_observed - lambda_emitted) / lambda_emitted;
    }

    /**
     * @brief Redshift from scale factor
     *
     * 1 + z = a(t_obs) / a(t_emit) = 1/a(t_emit)  [if a(t_obs) = 1 today]
     *
     * @param scale_factor a(t_emit)
     * @return Redshift z
     */
    static double fromScaleFactor(double scale_factor) {
        return (1.0 / scale_factor) - 1.0;
    }

    /**
     * @brief Scale factor from redshift
     *
     * a = 1/(1+z)
     */
    static double scaleFactorFromRedshift(double redshift) {
        return 1.0 / (1.0 + redshift);
    }

    /**
     * @brief Velocity from redshift (non-relativistic)
     *
     * v ≈ cz  (for z << 1)
     */
    static double velocityNonRelativistic(double redshift) {
        double c = 2.998e5;  // km/s
        return c * redshift;  // km/s
    }

    /**
     * @brief Relativistic redshift formula
     *
     * 1 + z = √[(1 + β)/(1 - β)]
     *
     * where β = v/c
     */
    static double relativisticRedshift(double beta) {
        return std::sqrt((1.0 + beta) / (1.0 - beta)) - 1.0;
    }

    /**
     * @brief Cosmological distance from redshift (approximate)
     *
     * For small z: d ≈ cz/H₀
     */
    static double luminosityDistance(double redshift) {
        double c_km_s = 2.998e5;
        double H0 = HubbleExpansion::hubbleConstant();

        // Approximate for low z
        if (redshift < 0.1) {
            return c_km_s * redshift / H0;  // Mpc
        }

        // Include expansion effects
        return c_km_s * redshift * (1.0 + redshift / 2.0) / H0;  // Mpc
    }
};

/**
 * @class ScaleFactor
 * @brief Scale factor a(t) - describes expansion
 *
 * Physical distance: r(t) = a(t) × r_comoving
 */
class ScaleFactor {
public:
    /**
     * @brief Scale factor convention
     *
     * a(t_today) = 1  (normalized to present)
     * a(t_past) < 1   (smaller in past)
     * a(t_future) > 1 (larger in future)
     */
    static double today() {
        return 1.0;
    }

    /**
     * @brief Scale factor at CMB decoupling
     *
     * z_CMB ≈ 1100 → a_CMB = 1/1101 ≈ 9×10⁻⁴
     */
    static double atCMBDecoupling() {
        double z_CMB = 1100.0;
        return 1.0 / (1.0 + z_CMB);  // ≈ 9×10⁻⁴
    }

    /**
     * @brief Scale factor at matter-radiation equality
     *
     * z_eq ≈ 3400 → a_eq ≈ 2.9×10⁻⁴
     */
    static double atMatterRadiationEquality() {
        double z_eq = 3400.0;
        return 1.0 / (1.0 + z_eq);
    }

    /**
     * @brief Scale factor at nucleosynthesis
     *
     * T_BBN ~ 1 MeV → z ~ 4×10⁸ → a ~ 2.5×10⁻⁹
     */
    static double atNucleosynthesis() {
        double z_BBN = 4e8;
        return 1.0 / (1.0 + z_BBN);
    }

    /**
     * @brief Hubble parameter H(a) = ȧ/a
     *
     * Rate of expansion at scale factor a
     *
     * @param scale_factor a
     * @param a_dot Time derivative ȧ (depends on cosmology)
     * @return H(a) in s⁻¹
     */
    static double hubbleParameter(double scale_factor, double a_dot) {
        if (scale_factor <= 0.0) {
            throw std::invalid_argument("Scale factor must be positive");
        }
        return a_dot / scale_factor;
    }
};

/**
 * @class OlbersParadox
 * @brief Olbers' paradox and its resolution
 *
 * Why is the night sky dark if universe is infinite and eternal?
 */
class OlbersParadox {
public:
    /**
     * @brief Statement of paradox
     *
     * In infinite, eternal, static universe filled with stars:
     * - Every line of sight hits a star
     * - Sky should be as bright as stellar surface
     * - But night sky is dark!
     */
    static std::string paradoxStatement() {
        return "Olbers' Paradox: If universe is infinite, eternal, and static,\n"
               "every line of sight should end on a star surface.\n"
               "Sky brightness = ∫₀^∞ (L/(4πr²)) × (4πr²n) dr → ∞\n"
               "But night sky is dark!";
    }

    /**
     * @brief Resolution 1: Finite age
     *
     * Universe has finite age (~13.8 Gyr)
     * Light from distant galaxies hasn't reached us yet
     * Observable universe is finite
     */
    static std::string resolutionFiniteAge() {
        return "Resolution 1: Finite age\n"
               "Universe is only ~13.8 Gyr old\n"
               "Horizon distance: d_H = ct ~ 13.8 Gly\n"
               "Only finite number of stars visible";
    }

    /**
     * @brief Resolution 2: Expansion
     *
     * Expansion redshifts light → energy diluted
     * Photon energy: E ∝ 1/a → decreases as universe expands
     */
    static std::string resolutionExpansion() {
        return "Resolution 2: Cosmological expansion\n"
               "Expansion redshifts photons: E ∝ 1/(1+z)\n"
               "Intensity: I ∝ 1/(1+z)⁴\n"
               "Distant sources heavily dimmed";
    }

    /**
     * @brief Resolution 3: Finite stellar lifetime
     *
     * Stars don't shine forever
     * Not enough time to fill universe with light
     */
    static std::string resolutionFiniteLifetime() {
        return "Resolution 3: Finite stellar lifetime\n"
               "Stars burn for ~10¹⁰ years, then die\n"
               "Not enough integrated light to fill sky";
    }

    /**
     * @brief Calculate sky brightness in static universe
     *
     * Assuming all stars like Sun, number density n
     *
     * @param star_luminosity L (W)
     * @param number_density n (stars/m³)
     * @param horizon_distance d_H (m)
     * @return Sky brightness (W/m²)
     */
    static double skyBrightnessStatic(double star_luminosity,
                                     double number_density,
                                     double horizon_distance) {
        // Integrate brightness over sphere
        // I = ∫₀^{d_H} (L/4πr²) × 4πr²n dr = L × n × d_H
        return star_luminosity * number_density * horizon_distance;
    }

    /**
     * @brief Sky brightness with expansion
     *
     * Include (1+z)⁴ dimming factor
     */
    static double skyBrightnessExpanding(double sky_brightness_static,
                                        double average_redshift) {
        double dimming = std::pow(1.0 + average_redshift, 4.0);
        return sky_brightness_static / dimming;
    }
};

/**
 * @class ObservableUniverse
 * @brief Observable universe and horizons
 */
class ObservableUniverse {
public:
    /**
     * @brief Particle horizon (comoving distance)
     *
     * Maximum distance light could have traveled since Big Bang
     *
     * r_p = ∫₀^t c dt'/a(t')
     *
     * For flat universe: r_p ≈ 46 Gly (today)
     */
    static double particleHorizon() {
        return 46.0;  // Gly (comoving distance)
    }

    /**
     * @brief Observable universe radius
     *
     * Physical size of observable universe
     */
    static double observableRadius() {
        return particleHorizon();  // ~46 Gly
    }

    /**
     * @brief Number of observable galaxies
     *
     * ~2 trillion galaxies within observable universe
     */
    static double numberOfGalaxies() {
        return 2e12;  // ~2 trillion
    }

    /**
     * @brief Event horizon
     *
     * Maximum distance we can ever communicate with
     * (due to accelerating expansion from dark energy)
     *
     * r_e = ∫_t^∞ c dt'/a(t')
     */
    static double eventHorizon() {
        return 16.0;  // Gly (with dark energy)
    }

    /**
     * @brief Hubble radius (present)
     *
     * c/H₀ ≈ 14 Gly
     */
    static double hubbleRadius() {
        double c = 2.998e5;  // km/s
        double H0 = HubbleExpansion::hubbleConstant();  // km/s/Mpc
        double Mpc_to_Gly = 3.262e-3;  // 1 Mpc = 3.262×10⁻³ Gly

        return (c / H0) * Mpc_to_Gly;  // Gly ~4.2 Gly
    }

    /**
     * @brief Surface of last scattering (CMB)
     *
     * Distance to CMB surface (comoving)
     */
    static double cmbrSurface() {
        return particleHorizon();  // CMB from edge of observable universe
    }
};

/**
 * @class CosmicTime
 * @brief Time since Big Bang
 */
class CosmicTime {
public:
    /**
     * @brief Age of universe (current)
     *
     * t₀ = 13.787 ± 0.020 Gyr (Planck 2018)
     */
    static double ageOfUniverse() {
        return 13.787;  // Gyr
    }

    /**
     * @brief Time at given redshift (approximate)
     *
     * For matter-dominated: t ∝ a^(3/2) ∝ (1+z)^(-3/2)
     * For radiation-dominated: t ∝ a² ∝ (1+z)^(-2)
     */
    static double timeAtRedshift(double redshift) {
        double t0 = ageOfUniverse();  // Gyr
        double a = 1.0 / (1.0 + redshift);

        // Approximate for matter-dominated era
        return t0 * std::pow(a, 1.5);  // Gyr
    }

    /**
     * @brief Lookback time
     *
     * How long ago we're seeing an object at redshift z
     *
     * t_lookback = t₀ - t(z)
     */
    static double lookbackTime(double redshift) {
        double t0 = ageOfUniverse();
        double t_z = timeAtRedshift(redshift);
        return t0 - t_z;  // Gyr
    }

    /**
     * @brief Time at CMB decoupling
     *
     * z = 1100, t ≈ 380,000 years
     */
    static double atCMBDecoupling() {
        return 380000.0 / 1e9;  // Gyr ≈ 3.8×10⁻⁴ Gyr
    }

    /**
     * @brief Time at matter-radiation equality
     *
     * z_eq ≈ 3400, t ≈ 47,000 years
     */
    static double atMatterRadiationEquality() {
        return 47000.0 / 1e9;  // Gyr
    }

    /**
     * @brief Time at nucleosynthesis
     *
     * T ~ 1 MeV, t ~ 1-3 minutes
     */
    static double atNucleosynthesis() {
        return 3.0 / (60.0 * 1e9);  // 3 minutes in Gyr
    }
};

} // namespace physics::advanced::cosmology

#endif // PHYSICS_ADVANCED_COSMOLOGY_EXPANDING_UNIVERSE_HPP
