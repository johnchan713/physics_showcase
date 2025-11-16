#ifndef PHYSICS_ADVANCED_COSMOLOGY_EARLY_UNIVERSE_HPP
#define PHYSICS_ADVANCED_COSMOLOGY_EARLY_UNIVERSE_HPP

#include <cmath>
#include <map>
#include <vector>
#include <string>

/**
 * @file early_universe.hpp
 * @brief Early universe: CMB, radiation/matter eras, nucleosynthesis, baryogenesis
 *
 * Implements:
 * - Cosmic Microwave Background (CMB)
 * - Radiation and matter dominated eras
 * - Big Bang Nucleosynthesis (BBN)
 * - Baryogenesis and matter-antimatter asymmetry
 * - Thermal history of universe
 */

namespace physics::advanced::cosmology {

/**
 * @class CosmicMicrowaveBackground
 * @brief CMB radiation - relic from early universe
 *
 * Photons from recombination era (z ~ 1100, t ~ 380,000 years)
 */
class CosmicMicrowaveBackground {
public:
    /**
     * @brief CMB temperature today
     *
     * T_CMB = 2.7255 ± 0.0006 K (COBE, WMAP, Planck)
     *
     * Perfect blackbody spectrum
     */
    static double temperatureToday() {
        return 2.7255;  // K
    }

    /**
     * @brief CMB temperature at redshift z
     *
     * T(z) = T₀(1 + z)
     *
     * Temperature scales with expansion
     */
    static double temperatureAtRedshift(double redshift) {
        return temperatureToday() * (1.0 + redshift);
    }

    /**
     * @brief Temperature at decoupling
     *
     * T_dec ~ 3000 K (z ~ 1100)
     */
    static double temperatureAtDecoupling() {
        double z_dec = 1100.0;
        return temperatureAtRedshift(z_dec);  // ~3000 K
    }

    /**
     * @brief Redshift of decoupling/recombination
     *
     * z_dec ≈ 1100
     *
     * When photons last scattered off electrons
     */
    static double decouplingRedshift() {
        return 1100.0;
    }

    /**
     * @brief Time of decoupling
     *
     * t_dec ≈ 380,000 years
     */
    static double decouplingTime() {
        return 380000.0;  // years
    }

    /**
     * @brief CMB energy density today
     *
     * ρ_γ = (π²/15)(k_B T)⁴/(ℏc)³
     *
     * Stefan-Boltzmann law for photons
     */
    static double energyDensityToday() {
        double T = temperatureToday();  // K
        double k_B = 1.381e-23;  // J/K
        double hbar = 1.055e-34;  // J·s
        double c = 2.998e8;  // m/s

        double prefactor = M_PI * M_PI / 15.0;
        double T4 = std::pow(k_B * T, 4.0);
        double denominator = std::pow(hbar * c, 3.0);

        return prefactor * T4 / denominator;  // J/m³
    }

    /**
     * @brief Planck blackbody spectrum
     *
     * B_ν(T) = (2hν³/c²) / [exp(hν/k_BT) - 1]
     *
     * @param frequency ν (Hz)
     * @param temperature T (K)
     * @return Spectral radiance (W/(m²·sr·Hz))
     */
    static double planckSpectrum(double frequency, double temperature) {
        double h = 6.626e-34;  // J·s
        double k_B = 1.381e-23;  // J/K
        double c = 2.998e8;  // m/s

        double numerator = 2.0 * h * std::pow(frequency, 3.0) / (c * c);
        double exponent = (h * frequency) / (k_B * temperature);

        if (exponent > 50.0) return 0.0;  // Avoid overflow

        double denominator = std::exp(exponent) - 1.0;

        return numerator / denominator;
    }

    /**
     * @brief CMB anisotropies
     *
     * ΔT/T ~ 10⁻⁵ (temperature fluctuations)
     *
     * Seeds of structure formation
     */
    static double anisotropyLevel() {
        return 1e-5;  // ΔT/T
    }

    /**
     * @brief Acoustic peaks in CMB power spectrum
     *
     * Locations encode cosmological parameters
     */
    static std::vector<double> acousticPeakAngularScales() {
        // Angular scales of first few peaks (degrees)
        return {1.0, 0.5, 0.33};  // ℓ ~ 220, 540, 800
    }

    /**
     * @brief Silk damping scale
     *
     * Photon diffusion erases small-scale fluctuations
     *
     * λ_D ~ 1 Mpc (comoving)
     */
    static double silkDampingScale() {
        return 1.0;  // Mpc (comoving)
    }
};

/**
 * @class RadiationEra
 * @brief Radiation-dominated epoch
 *
 * Early universe: t < t_eq ~ 47,000 years
 * ρ_radiation > ρ_matter
 */
class RadiationEra {
public:
    /**
     * @brief Matter-radiation equality
     *
     * z_eq ≈ 3400, t_eq ≈ 47,000 years
     */
    static double equalityRedshift() {
        return 3400.0;
    }

    static double equalityTime() {
        return 47000.0;  // years
    }

    /**
     * @brief Scale factor evolution in radiation era
     *
     * a(t) ∝ t^(1/2)
     *
     * For radiation-dominated (w = 1/3)
     */
    static double scaleFactorEvolution(double time_ratio) {
        return std::sqrt(time_ratio);  // a ∝ t^(1/2)
    }

    /**
     * @brief Temperature evolution in radiation era
     *
     * T ∝ a⁻¹ ∝ t^(-1/2)
     *
     * @param time t (seconds)
     * @param T_ref Reference temperature (K)
     * @param t_ref Reference time (s)
     * @return Temperature at time t
     */
    static double temperatureEvolution(double time, double T_ref, double t_ref) {
        return T_ref * std::sqrt(t_ref / time);  // T ∝ t^(-1/2)
    }

    /**
     * @brief Hubble parameter in radiation era
     *
     * H(t) = 1/(2t)
     *
     * Simple relation: age ≈ 1/(2H)
     */
    static double hubbleParameter(double time) {
        if (time <= 0.0) {
            throw std::invalid_argument("Time must be positive");
        }
        return 1.0 / (2.0 * time);  // s⁻¹
    }

    /**
     * @brief Radiation energy density
     *
     * ρ_r = (π²/30) g_* (k_B T)⁴ / (ℏc)³
     *
     * where g_* is effective degrees of freedom
     */
    static double energyDensity(double temperature, double g_eff = 3.36) {
        double k_B = 1.381e-23;  // J/K
        double hbar = 1.055e-34;  // J·s
        double c = 2.998e8;  // m/s

        double prefactor = (M_PI * M_PI / 30.0) * g_eff;
        double T4 = std::pow(k_B * temperature, 4.0);
        double denominator = std::pow(hbar * c, 3.0);

        return prefactor * T4 / denominator;  // J/m³
    }

    /**
     * @brief Effective degrees of freedom
     *
     * Counts relativistic particle species
     *
     * Today: g_*,0 = 3.36 (photons + neutrinos)
     * Early universe: g_* ~ 100 (all SM particles)
     */
    static double effectiveDegreesOfFreedom(double temperature_MeV) {
        if (temperature_MeV > 300.0) {
            return 106.75;  // All SM particles relativistic
        } else if (temperature_MeV > 1.0) {
            return 10.75;  // After quark-hadron transition
        } else if (temperature_MeV > 0.5) {
            return 3.36;  // After e+e- annihilation
        } else {
            return 3.36;  // Today (photons + neutrinos)
        }
    }
};

/**
 * @class MatterEra
 * @brief Matter-dominated epoch
 *
 * t > t_eq until dark energy takes over
 * ρ_matter > ρ_radiation
 */
class MatterEra {
public:
    /**
     * @brief Start of matter era
     *
     * z_eq ≈ 3400, t_eq ≈ 47,000 years
     */
    static double startRedshift() {
        return RadiationEra::equalityRedshift();
    }

    /**
     * @brief End of matter era (dark energy takes over)
     *
     * z_Λ ≈ 0.4, t_Λ ≈ 10 Gyr
     */
    static double endRedshift() {
        return 0.4;
    }

    /**
     * @brief Scale factor evolution in matter era
     *
     * a(t) ∝ t^(2/3)
     *
     * For matter-dominated (w = 0)
     */
    static double scaleFactorEvolution(double time_ratio) {
        return std::pow(time_ratio, 2.0/3.0);  // a ∝ t^(2/3)
    }

    /**
     * @brief Hubble parameter in matter era
     *
     * H(t) = 2/(3t)
     *
     * Age ≈ 2/(3H)
     */
    static double hubbleParameter(double time) {
        if (time <= 0.0) {
            throw std::invalid_argument("Time must be positive");
        }
        return 2.0 / (3.0 * time);  // s⁻¹
    }

    /**
     * @brief Structure formation
     *
     * Matter era allows gravitational instability
     * Density perturbations grow: δρ/ρ ∝ a (linear regime)
     */
    static double densityPerturbationGrowth(double scale_factor) {
        // Linear growth in matter era
        return scale_factor;  // δ ∝ a
    }
};

/**
 * @class BigBangNucleosynthesis
 * @brief Primordial nucleosynthesis (BBN)
 *
 * Formation of light elements in first few minutes
 * T ~ 1 MeV → 0.1 MeV, t ~ 1 s → 3 minutes
 */
class BigBangNucleosynthesis {
public:
    /**
     * @brief BBN temperature range
     *
     * T ~ 1 MeV (start) → 0.1 MeV (end)
     */
    static std::pair<double, double> temperatureRange() {
        return {1.0, 0.1};  // MeV
    }

    /**
     * @brief BBN time range
     *
     * t ~ 1 second → 3 minutes
     */
    static std::pair<double, double> timeRange() {
        return {1.0, 180.0};  // seconds
    }

    /**
     * @brief BBN redshift
     *
     * z_BBN ~ 4×10⁸
     */
    static double redshift() {
        return 4e8;
    }

    /**
     * @brief Primordial abundances (by mass)
     *
     * Predicted by BBN, observed in pristine gas
     */
    static std::map<std::string, double> primordialAbundances() {
        return {
            {"H", 0.75},          // Hydrogen: 75%
            {"He-4", 0.25},       // Helium-4: ~25%
            {"He-3", 1e-5},       // Helium-3: ~10⁻⁵
            {"D", 2.5e-5},        // Deuterium: ~2.5×10⁻⁵
            {"Li-7", 4e-10}       // Lithium-7: ~4×10⁻¹⁰
        };
    }

    /**
     * @brief Deuterium bottleneck
     *
     * Deuterium photodissociation prevents He formation until T < 0.1 MeV
     *
     * D + γ → p + n (binding energy ~2.2 MeV)
     *
     * Must wait for photon energies to drop below binding energy
     */
    static double deuteriumBindingEnergy() {
        return 2.224;  // MeV
    }

    /**
     * @brief Neutron-proton ratio
     *
     * Freezes out at T ~ 1 MeV
     * n/p ≈ 1/7 at BBN (after β-decay)
     *
     * Nearly all neutrons end up in He-4
     */
    static double neutronProtonRatio() {
        return 1.0/7.0;  // At t ~ 1 minute
    }

    /**
     * @brief Helium-4 mass fraction
     *
     * Y_p = ρ(He-4) / ρ(baryons) ≈ 0.25
     *
     * From n/p ratio: Y_p ≈ 2(n/p) / [1 + (n/p)]
     */
    static double helium4MassFraction() {
        double n_over_p = neutronProtonRatio();
        return 2.0 * n_over_p / (1.0 + n_over_p);  // ~0.25
    }

    /**
     * @brief Dependence on baryon-to-photon ratio
     *
     * Deuterium abundance sensitive to η = n_b/n_γ
     *
     * Higher η → more D burned to He
     * D/H decreases with η
     */
    static double deuteriumAbundance(double eta_10) {
        // Approximate fit: D/H = 2.5×10⁻⁵ (η_10/6)^(-1.6)
        // where η_10 = η × 10¹⁰
        return 2.5e-5 * std::pow(eta_10 / 6.0, -1.6);
    }

    /**
     * @brief Number of neutrino species constraint
     *
     * BBN constrains N_ν = 2.99 ± 0.33 (consistent with 3 generations)
     *
     * More neutrinos → faster expansion → more He-4
     */
    static double neutrinoSpeciesConstraint() {
        return 2.99;  // N_ν (from BBN)
    }
};

/**
 * @class Baryogenesis
 * @brief Generation of baryon asymmetry
 *
 * Why is universe made of matter, not antimatter?
 * η = (n_b - n_b̄) / n_γ ~ 6×10⁻¹⁰
 */
class Baryogenesis {
public:
    /**
     * @brief Baryon-to-photon ratio
     *
     * η = n_b / n_γ ≈ 6×10⁻¹⁰
     *
     * Tiny excess of matter over antimatter
     */
    static double baryonToPhotonRatio() {
        return 6.1e-10;
    }

    /**
     * @brief Baryon asymmetry parameter
     *
     * η_B = (n_b - n_b̄) / s
     *
     * where s is entropy density
     */
    static double baryonAsymmetry() {
        return 8.7e-11;  // η_B
    }

    /**
     * @brief Sakharov conditions for baryogenesis
     *
     * Three necessary conditions (Sakharov 1967):
     * 1. Baryon number (B) violation
     * 2. C and CP violation
     * 3. Departure from thermal equilibrium
     */
    static std::vector<std::string> sakharovConditions() {
        return {
            "1. Baryon number violation: Processes that change B",
            "2. C and CP violation: Distinguish matter from antimatter",
            "3. Non-equilibrium: Prevent washout of asymmetry"
        };
    }

    /**
     * @brief Electroweak baryogenesis
     *
     * Generate asymmetry at electroweak phase transition
     * T ~ 100 GeV, t ~ 10⁻¹¹ s
     *
     * Requires: first-order phase transition + CP violation
     * Standard Model CP violation too weak!
     */
    static std::string electroweakBaryogenesis() {
        return "Electroweak baryogenesis:\n"
               "Temperature: T ~ 100 GeV\n"
               "Time: t ~ 10⁻¹¹ s\n"
               "Requires: Strong 1st order phase transition\n"
               "Problem: SM CP violation insufficient";
    }

    /**
     * @brief GUT baryogenesis
     *
     * Generate asymmetry at GUT scale
     * T ~ 10¹⁵ GeV, t ~ 10⁻³⁷ s
     *
     * X, Y bosons decay: X → qq vs X̄ → q̄q̄ with different rates
     */
    static std::string gutBaryogenesis() {
        return "GUT baryogenesis:\n"
               "Temperature: T ~ 10¹⁵ GeV\n"
               "Time: t ~ 10⁻³⁷ s\n"
               "Mechanism: X boson decay (Γ(X→qq) ≠ Γ(X̄→q̄q̄))\n"
               "Naturally produces η_B ~ 10⁻¹⁰";
    }

    /**
     * @brief Leptogenesis
     *
     * Generate lepton asymmetry, convert to baryon asymmetry
     * Via sphaleron processes
     *
     * Heavy right-handed neutrinos decay
     */
    static std::string leptogenesis() {
        return "Leptogenesis:\n"
               "Generate lepton asymmetry from heavy N_R decay\n"
               "Sphalerons convert (B-L) conserving, B+L violating\n"
               "Lepton asymmetry → Baryon asymmetry\n"
               "Explains neutrino masses + baryon asymmetry!";
    }

    /**
     * @brief Sphaleron processes
     *
     * Non-perturbative processes that violate B+L
     * Conserve B-L
     *
     * Active at high temperature (T > 100 GeV)
     */
    static std::string sphalerons() {
        return "Sphaleron processes:\n"
               "Violate B+L: ΔB = ΔL\n"
               "Conserve B-L\n"
               "Active for T > T_EW ~ 100 GeV\n"
               "Rate: Γ_sph ~ α_W⁵ T⁴";
    }

    /**
     * @brief Matter-antimatter annihilation era
     *
     * Most matter and antimatter annihilated
     * Small excess (η ~ 10⁻⁹) survived as all baryons today
     */
    static std::string annihilationEra() {
        return "Matter-antimatter annihilation:\n"
               "Time: t ~ 1 second, T ~ 1 MeV\n"
               "For every 10⁹ antimatter particles, 10⁹ + 1 matter\n"
               "After annihilation: 1 proton remains per 10⁹ photons\n"
               "Result: η = n_b/n_γ ~ 6×10⁻¹⁰";
    }
};

/**
 * @class ThermalHistory
 * @brief Key epochs in cosmic history
 */
class ThermalHistory {
public:
    /**
     * @brief Planck epoch
     *
     * T > 10¹⁹ GeV, t < 10⁻⁴³ s
     * Quantum gravity era
     */
    static std::pair<double, double> planckEpoch() {
        return {1e19, 1e-43};  // GeV, seconds
    }

    /**
     * @brief GUT epoch
     *
     * T ~ 10¹⁵ GeV, t ~ 10⁻³⁷ s
     * Grand unification, X/Y bosons
     */
    static std::pair<double, double> gutEpoch() {
        return {1e15, 1e-37};  // GeV, seconds
    }

    /**
     * @brief Electroweak epoch
     *
     * T ~ 100 GeV, t ~ 10⁻¹¹ s
     * Electroweak symmetry breaking, Higgs gets VEV
     */
    static std::pair<double, double> electroweakEpoch() {
        return {100.0, 1e-11};  // GeV, seconds
    }

    /**
     * @brief QCD epoch
     *
     * T ~ 200 MeV, t ~ 10⁻⁵ s
     * Quark-hadron transition, confinement
     */
    static std::pair<double, double> qcdEpoch() {
        return {0.2, 1e-5};  // GeV, seconds
    }

    /**
     * @brief Neutrino decoupling
     *
     * T ~ 1 MeV, t ~ 1 s
     * Neutrinos free-stream, form CνB
     */
    static std::pair<double, double> neutrinoDecoupling() {
        return {0.001, 1.0};  // GeV, seconds
    }

    /**
     * @brief Big Bang Nucleosynthesis
     *
     * T ~ 1-0.1 MeV, t ~ 1-180 s
     * Light element formation
     */
    static std::pair<double, double> bbnEpoch() {
        return {0.001, 180.0};  // GeV, seconds
    }

    /**
     * @brief Matter-radiation equality
     *
     * T ~ 0.75 eV, t ~ 47,000 years
     * ρ_m = ρ_r
     */
    static std::pair<double, double> matterRadiationEquality() {
        return {7.5e-10, 47000.0 * 365.25 * 24 * 3600};  // GeV, seconds
    }

    /**
     * @brief Recombination
     *
     * T ~ 0.26 eV, t ~ 380,000 years
     * Atoms form, photons decouple → CMB
     */
    static std::pair<double, double> recombination() {
        return {2.6e-10, 380000.0 * 365.25 * 24 * 3600};  // GeV, seconds
    }

    /**
     * @brief Dark ages
     *
     * t ~ 380,000 - 100 million years
     * No light sources, only CMB
     */
    static std::pair<double, double> darkAges() {
        return {380000.0, 1e8};  // years
    }

    /**
     * @brief Reionization
     *
     * z ~ 6-20, t ~ 150 million - 1 billion years
     * First stars ionize neutral hydrogen
     */
    static std::pair<double, double> reionization() {
        return {150e6, 1e9};  // years
    }
};

} // namespace physics::advanced::cosmology

#endif // PHYSICS_ADVANCED_COSMOLOGY_EARLY_UNIVERSE_HPP
