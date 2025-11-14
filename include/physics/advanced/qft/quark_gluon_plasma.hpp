#ifndef PHYSICS_ADVANCED_QFT_QUARK_GLUON_PLASMA_HPP
#define PHYSICS_ADVANCED_QFT_QUARK_GLUON_PLASMA_HPP

#include <cmath>
#include <string>

/**
 * @file quark_gluon_plasma.hpp
 * @brief Quark-gluon plasma (QGP) - deconfined quark matter
 *
 * Implements:
 * - QGP phase transition
 * - Critical temperature
 * - Equation of state
 * - Deconfinement and chiral symmetry restoration
 * - Heavy-ion collision signatures
 */

namespace physics::advanced::qft {

/**
 * @class QuarkGluonPlasma
 * @brief Properties of deconfined quark-gluon matter
 *
 * At high temperature/density, quarks and gluons are no longer
 * confined in hadrons but form a plasma phase
 */
class QuarkGluonPlasma {
public:
    /**
     * @brief Critical temperature for QGP formation
     *
     * T_c ≈ 150-170 MeV ≈ 2×10¹² K
     *
     * Above T_c: quarks and gluons deconfined (QGP)
     * Below T_c: quarks confined in hadrons
     */
    static double criticalTemperature() {
        return 0.160;  // GeV ≈ 160 MeV
    }

    /**
     * @brief Critical temperature in Kelvin
     */
    static double criticalTemperatureKelvin() {
        double T_c_GeV = criticalTemperature();
        double k_B = 8.617e-5;  // eV/K

        return T_c_GeV * 1e9 / k_B;  // ~2×10¹² K
    }

    /**
     * @brief Critical energy density
     *
     * ε_c ≈ 0.3-1.0 GeV/fm³
     */
    static double criticalEnergyDensity() {
        return 0.5;  // GeV/fm³
    }

    /**
     * @brief QCD coupling at QGP temperatures
     *
     * α_s(T) decreases at high T (asymptotic freedom)
     */
    static double couplingAtTemperature(double temperature) {
        // Simplified running
        double T_c = criticalTemperature();
        double alpha_s_Tc = 0.3;  // At T_c

        // α_s decreases logarithmically with T
        return alpha_s_Tc / std::log(temperature / T_c + 2.718);
    }
};

/**
 * @class PhaseTransition
 * @brief Hadronic matter → QGP phase transition
 */
class PhaseTransition {
public:
    /**
     * @brief Latent heat of phase transition
     *
     * ΔE ~ 1 GeV/fm³
     *
     * Energy required to deconfine quarks
     */
    static double latentHeat() {
        return 1.0;  // GeV/fm³
    }

    /**
     * @brief Order of phase transition
     *
     * Lattice QCD suggests crossover (not first-order)
     * at μ_B = 0 (vanishing baryon chemical potential)
     */
    static std::string transitionOrder() {
        return "Crossover transition at μ_B = 0 (smooth, not first-order)";
    }

    /**
     * @brief Deconfinement temperature vs baryon density
     *
     * Phase diagram: T vs μ_B
     *
     * T_c decreases with increasing μ_B
     * Possible critical endpoint
     */
    static double transitionTemperature(double baryon_chemical_potential) {
        double T_c0 = QuarkGluonPlasma::criticalTemperature();

        // Simplified: T_c(μ_B) ≈ T_c(0)[1 - (μ_B/T_c)²]
        double ratio = baryon_chemical_potential / T_c0;

        return T_c0 * (1.0 - ratio * ratio);
    }

    /**
     * @brief Chiral symmetry restoration
     *
     * In QGP, chiral symmetry is restored
     * Quarks effectively massless (constituent mass → current mass)
     */
    static std::string chiralRestoration() {
        return "Chiral symmetry restored in QGP: m_q(constituent) → m_q(current)";
    }
};

/**
 * @class EquationOfState
 * @brief Thermodynamics of QGP
 */
class EquationOfState {
public:
    /**
     * @brief Energy density (Stefan-Boltzmann)
     *
     * ε = (π²/30) g_eff T⁴
     *
     * where g_eff is effective degrees of freedom
     *
     * For QGP: g_eff ≈ 37 (quarks) + 16 (gluons) ≈ 50-60
     */
    static double energyDensity(double temperature) {
        double g_eff = 50.0;  // Effective DOF for QGP

        return (M_PI * M_PI / 30.0) * g_eff * std::pow(temperature, 4.0);  // GeV⁴
    }

    /**
     * @brief Pressure
     *
     * P = ε/3 (relativistic ideal gas)
     */
    static double pressure(double temperature) {
        return energyDensity(temperature) / 3.0;  // GeV⁴
    }

    /**
     * @brief Entropy density
     *
     * s = (2π²/45) g_eff T³
     */
    static double entropyDensity(double temperature) {
        double g_eff = 50.0;

        return (2.0 * M_PI * M_PI / 45.0) * g_eff * std::pow(temperature, 3.0);  // GeV³
    }

    /**
     * @brief Sound speed
     *
     * c_s² = ∂P/∂ε = 1/3 (ideal QGP)
     *
     * @return c_s/c (dimensionless)
     */
    static double soundSpeed() {
        return 1.0 / std::sqrt(3.0);  // ≈ 0.577c
    }

    /**
     * @brief Trace anomaly (interaction measure)
     *
     * (ε - 3P)/T⁴
     *
     * Zero for ideal gas, non-zero for interacting QGP
     */
    static double traceAnomaly(double temperature) {
        // Measure of deviation from ideal gas
        // Lattice QCD: peaks near T_c, then decreases
        double T_c = QuarkGluonPlasma::criticalTemperature();
        double ratio = temperature / T_c;

        // Simplified behavior
        if (ratio < 1.0) return 5.0;  // Hadronic phase
        if (ratio < 2.0) return 3.0 * std::exp(-(ratio - 1.0));  // Peak near T_c
        return 0.5;  // Approach ideal gas at high T
    }
};

/**
 * @class HeavyIonCollisions
 * @brief Signatures of QGP in heavy-ion experiments
 *
 * RHIC (Brookhaven) and LHC (CERN) create QGP
 */
class HeavyIonCollisions {
public:
    /**
     * @brief Bjorken energy density estimate
     *
     * ε = (dE_T/dη) / (τ_0 A_T)
     *
     * where dE_T/dη is transverse energy per rapidity
     *       τ_0 ~ 1 fm/c is formation time
     *       A_T is transverse area
     */
    static double bjorkenEnergyDensity(double dET_deta, double formation_time,
                                      double transverse_area) {
        return dET_deta / (formation_time * transverse_area);  // GeV/fm³
    }

    /**
     * @brief QGP formation time
     *
     * τ_0 ~ 0.5-1 fm/c ≈ 10⁻²⁴ s
     */
    static double formationTime() {
        return 1.0;  // fm/c
    }

    /**
     * @brief QGP lifetime
     *
     * τ_QGP ~ 5-10 fm/c ≈ 10⁻²³ s
     *
     * Before hadronization
     */
    static double lifetimeQGP() {
        return 7.0;  // fm/c
    }

    /**
     * @brief Initial temperature at RHIC
     *
     * T_initial ~ 300-400 MeV ≈ 2T_c
     */
    static double initialTemperatureRHIC() {
        return 0.350;  // GeV
    }

    /**
     * @brief Initial temperature at LHC
     *
     * T_initial ~ 500-600 MeV ≈ 3-4 T_c
     */
    static double initialTemperatureLHC() {
        return 0.550;  // GeV
    }
};

/**
 * @class QGPSignatures
 * @brief Experimental signatures of QGP formation
 */
class QGPSignatures {
public:
    /**
     * @brief Jet quenching
     *
     * High-p_T jets lose energy traversing QGP
     * R_AA < 1 (nuclear modification factor)
     */
    static std::string jetQuenching() {
        return "Jet quenching: R_AA = (dN/dp_T)^AA / <N_coll>(dN/dp_T)^pp < 1";
    }

    /**
     * @brief Elliptic flow (v₂)
     *
     * Asymmetric momentum distribution indicates collective flow
     * Hydro behavior → strongly coupled QGP
     */
    static std::string ellipticFlow() {
        return "Elliptic flow v₂: particles flow preferentially in plane\n"
               "v₂ ~ 0.1-0.2 at RHIC/LHC → near-perfect liquid behavior";
    }

    /**
     * @brief J/ψ suppression
     *
     * Charm quarkonia melted in QGP
     * Debye screening of color charge
     */
    static std::string jpsiSuppression() {
        return "J/ψ suppression: charmonium dissociates in QGP\n"
               "Color Debye screening length λ_D ~ 1/(gT)";
    }

    /**
     * @brief Strangeness enhancement
     *
     * Increased strange quark production in QGP
     * s̄s pairs easier to create than in hadronic matter
     */
    static std::string strangenessEnhancement() {
        return "Strangeness enhancement: K/π ratio increases\n"
               "Strange quarks more abundant in QGP";
    }

    /**
     * @brief Direct photon emission
     *
     * Thermal photons emitted from QGP
     * Clean probe (no strong interaction)
     */
    static std::string directPhotons() {
        return "Direct photons: γ emitted from QGP phase\n"
               "dN_γ/dy ~ T⁴ → temperature measurement";
    }

    /**
     * @brief Shear viscosity to entropy ratio
     *
     * η/s ≈ 1/(4π) (near lower bound!)
     *
     * QGP is near-perfect fluid (minimal viscosity)
     */
    static double viscosityToEntropyRatio() {
        return 1.0 / (4.0 * M_PI);  // ≈ 0.08 (in ℏ/k_B units)
    }
};

/**
 * @class ColorDebyeScreening
 * @brief Debye screening in QGP
 */
class ColorDebyeScreening {
public:
    /**
     * @brief Debye screening length
     *
     * λ_D = 1/(gT) where g ~ √(4πα_s)
     *
     * Screens color charge over distance λ_D
     */
    static double debyeLength(double temperature) {
        double alpha_s = 0.3;  // At T ~ T_c
        double g = std::sqrt(4.0 * M_PI * alpha_s);

        return 1.0 / (g * temperature);  // GeV⁻¹ = 0.197 fm
    }

    /**
     * @brief Heavy quark potential in QGP
     *
     * V(r) ~ -α_s e^(-m_D r)/r (Yukawa)
     *
     * where m_D = g T is Debye mass
     */
    static double heavyQuarkPotential(double distance, double temperature) {
        double m_D = debyeLength(temperature);
        double alpha_s = 0.3;

        return -alpha_s * std::exp(-m_D * distance) / distance;  // GeV
    }

    /**
     * @brief Quarkonium dissociation temperature
     *
     * T_diss ~ binding energy / k_B
     *
     * J/ψ: T_diss ~ 2T_c
     * Υ: T_diss ~ 4T_c (more tightly bound)
     */
    static double dissociationTemperature(double binding_energy) {
        return binding_energy;  // T ~ E_b (in natural units)
    }
};

/**
 * @class EarlyUniverse
 * @brief QGP in cosmology
 */
class EarlyUniverse {
public:
    /**
     * @brief QGP epoch in early universe
     *
     * t ~ 10 μs after Big Bang
     * T ~ 150-200 MeV
     *
     * Universe was in QGP phase
     */
    static double qgpEpochTime() {
        return 1e-5;  // seconds (10 microseconds)
    }

    /**
     * @brief Hadronization in early universe
     *
     * t ~ 20 μs: QGP → hadrons
     *
     * Quarks confined into protons, neutrons, mesons
     */
    static double hadronizationTime() {
        return 2e-5;  // seconds
    }

    /**
     * @brief Temperature evolution
     *
     * T ∝ 1/√t (radiation dominated)
     */
    static double temperatureAtTime(double time) {
        double T_1s = 1.0e-3;  // GeV at t = 1s

        return T_1s / std::sqrt(time);  // GeV
    }
};

} // namespace physics::advanced::qft

#endif // PHYSICS_ADVANCED_QFT_QUARK_GLUON_PLASMA_HPP
