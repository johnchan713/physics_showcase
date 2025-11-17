#ifndef PHYSICS_ADVANCED_QFT_CROSS_SECTIONS_HPP
#define PHYSICS_ADVANCED_QFT_CROSS_SECTIONS_HPP

#include <cmath>
#include <string>
#include <functional>

/**
 * @file cross_sections.hpp
 * @brief Interaction cross sections for particle processes
 *
 * Implements:
 * - Scattering cross sections
 * - Total and differential cross sections
 * - Rutherford scattering
 * - e⁺e⁻ annihilation
 * - Parton distribution functions
 */

namespace physics::advanced::qft {

/**
 * @class CrossSection
 * @brief Basic cross section calculations
 *
 * Cross section σ has units of area (barns, fm², m²)
 * 1 barn = 10⁻²⁸ m² = 100 fm²
 */
class CrossSection {
public:
    /**
     * @brief Unit conversions
     */
    static constexpr double barn_to_m2() { return 1e-28; }
    static constexpr double GeV2_to_barn() { return 0.3894;  }  // GeV⁻² to mb

    /**
     * @brief Differential to total cross section
     *
     * σ_total = ∫ dσ/dΩ dΩ
     *
     * Integrate over solid angle
     */
    static double totalFromDifferential(
        std::function<double(double)> dsigma_dOmega,
        double theta_min = 0.0, double theta_max = M_PI) {
        // Numerical integration (trapezoidal)
        int N = 1000;
        double dtheta = (theta_max - theta_min) / N;
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double theta = theta_min + (i + 0.5) * dtheta;
            // dΩ = sin(θ)dθdφ, integrate φ from 0 to 2π
            sum += dsigma_dOmega(theta) * std::sin(theta) * 2.0 * M_PI * dtheta;
        }

        return sum;
    }

    /**
     * @brief Luminosity and event rate
     *
     * N_events = L × σ
     *
     * where L is integrated luminosity
     */
    static double eventRate(double luminosity, double cross_section) {
        return luminosity * cross_section;  // events
    }

    /**
     * @brief LHC instantaneous luminosity
     *
     * L ~ 10³⁴ cm⁻²s⁻¹ = 10⁴² m⁻²s⁻¹
     */
    static double lhc_luminosity() {
        return 1e34;  // cm⁻²s⁻¹
    }
};

/**
 * @class RutherfordScattering
 * @brief Rutherford scattering (classical Coulomb)
 *
 * α + nucleus → α + nucleus
 */
class RutherfordScattering {
public:
    /**
     * @brief Rutherford differential cross section
     *
     * dσ/dΩ = (Z₁Z₂α ℏc / 4E)² × 1/sin⁴(θ/2)
     *
     * @param atomic_number_1 Z₁ (projectile)
     * @param atomic_number_2 Z₂ (target)
     * @param energy E (GeV)
     * @param theta Scattering angle (radians)
     * @return dσ/dΩ (barn/steradian)
     */
    static double differential(int atomic_number_1, int atomic_number_2,
                              double energy, double theta) {
        double alpha = 1.0/137.0;
        double hbar_c = 0.1973;  // GeV·fm

        double numerator = (atomic_number_1 * atomic_number_2 * alpha * hbar_c) /
                          (4.0 * energy);
        double sin_half = std::sin(theta / 2.0);
        double sin4 = sin_half * sin_half * sin_half * sin_half;

        return (numerator * numerator) / sin4;  // fm² = 0.01 barn
    }

    /**
     * @brief Mott scattering (including spin)
     *
     * Quantum correction to Rutherford formula
     */
    static double mottCorrection(double theta) {
        double cos_half = std::cos(theta / 2.0);
        return cos_half * cos_half;  // (1 - v²/c² sin²(θ/2))
    }
};

/**
 * @class QEDProcesses
 * @brief QED scattering cross sections
 */
class QEDProcesses {
public:
    /**
     * @brief e⁺e⁻ → μ⁺μ⁻ cross section
     *
     * σ = 4πα²/(3s)
     *
     * where √s is center-of-mass energy
     *
     * @param sqrt_s √s (GeV)
     * @return Cross section (barn)
     */
    static double electronMuonScattering(double sqrt_s) {
        double alpha = 1.0/137.0;
        double s = sqrt_s * sqrt_s;

        double sigma_GeV2 = (4.0 * M_PI * alpha * alpha) / (3.0 * s);

        return sigma_GeV2 * CrossSection::GeV2_to_barn();  // barn
    }

    /**
     * @brief Bhabha scattering: e⁺e⁻ → e⁺e⁻
     *
     * More complex due to t-channel and s-channel
     */
    static double bhabhaScattering(double sqrt_s) {
        // Simplified total cross section
        double alpha = 1.0/137.0;
        double s = sqrt_s * sqrt_s;

        // Contains both annihilation and scattering
        double sigma_GeV2 = (12.0 * M_PI * alpha * alpha) / s;

        return sigma_GeV2 * CrossSection::GeV2_to_barn();  // barn
    }

    /**
     * @brief Compton scattering: γe⁻ → γe⁻
     *
     * Klein-Nishina formula at high energy
     */
    static double comptonScattering(double photon_energy) {
        double r_e = 2.818e-15;  // Classical electron radius (m)
        double m_e = 0.000511;   // GeV

        // High energy limit: σ ~ (3σ_T/8)(m/E)ln(2E/m)
        double sigma_T = (8.0 * M_PI / 3.0) * r_e * r_e;  // Thomson cross section

        if (photon_energy < m_e) {
            return sigma_T;  // Low energy (Thomson)
        }

        // High energy (Klein-Nishina)
        double ratio = m_e / photon_energy;
        return (3.0 * sigma_T / 8.0) * ratio * std::log(2.0 / ratio);  // m²
    }
};

/**
 * @class HadronicCrossSections
 * @brief Strong interaction cross sections
 */
class HadronicCrossSections {
public:
    /**
     * @brief e⁺e⁻ → hadrons total cross section
     *
     * R = σ(e⁺e⁻ → hadrons) / σ(e⁺e⁻ → μ⁺μ⁻)
     *     = 3 Σ Q_q²
     *
     * where sum is over active quark flavors
     *
     * R = 3(4/9 + 1/9 + 1/9) = 2 for u,d,s
     * R = 10/3 for u,d,s,c
     * R = 11/3 for u,d,s,c,b
     */
    static double Rratio(double sqrt_s) {
        // Number of active flavors depends on energy
        if (sqrt_s < 3.0) {
            // u, d, s quarks (below charm threshold)
            return 3.0 * (4.0/9.0 + 1.0/9.0 + 1.0/9.0);  // R = 2
        } else if (sqrt_s < 10.0) {
            // u, d, s, c quarks (below bottom threshold)
            return 3.0 * (4.0/9.0 + 1.0/9.0 + 1.0/9.0 + 4.0/9.0);  // R = 10/3
        } else {
            // u, d, s, c, b quarks
            return 3.0 * (4.0/9.0 + 1.0/9.0 + 1.0/9.0 + 4.0/9.0 + 1.0/9.0);  // R = 11/3
        }
    }

    /**
     * @brief Total hadronic cross section
     */
    static double hadronicCrossSection(double sqrt_s) {
        double muon_sigma = QEDProcesses::electronMuonScattering(sqrt_s);
        return Rratio(sqrt_s) * muon_sigma;
    }

    /**
     * @brief Geometric cross section (strong)
     *
     * σ ~ πR² where R ~ 1 fm
     *
     * Typical hadronic cross sections ~ 10-100 mb
     */
    static double geometricCrossSection(double radius = 1.0e-15) {
        return M_PI * radius * radius;  // m²
    }

    /**
     * @brief Proton-proton total cross section
     *
     * σ_tot(pp) ~ 100 mb at LHC energies (√s ~ 13 TeV)
     */
    static double protonProtonTotal(double sqrt_s) {
        // Empirical fit (mb)
        double sigma_0 = 21.7;
        double B = 0.56;
        double Y = 0.308;
        double s_0 = 29.1;

        double s = sqrt_s * sqrt_s;
        return sigma_0 + B * std::pow(std::log(s/s_0), 2.0) +
               Y * std::pow(s/s_0, -0.5);  // mb
    }
};

/**
 * @class WeakProcesses
 * @brief Weak interaction cross sections
 */
class WeakProcesses {
public:
    /**
     * @brief Neutrino-nucleon cross section
     *
     * σ(νN) ~ G_F² s / π
     *
     * where G_F is Fermi constant, s is Mandelstam variable
     *
     * Linear in energy!
     */
    static double neutrinoNucleon(double neutrino_energy) {
        double G_F = 1.166e-5;  // GeV⁻²

        // σ ~ G_F² E / π (in natural units)
        double sigma_GeV2 = (G_F * G_F * neutrino_energy) / M_PI;

        return sigma_GeV2 * CrossSection::GeV2_to_barn();  // barn
    }

    /**
     * @brief Beta decay lifetime
     *
     * τ ~ 1/(G_F² m_e⁵)
     *
     * Neutron: τ ~ 880 s
     */
    static double betaDecayLifetime(double Q_value) {
        double G_F = 1.166e-5;  // GeV⁻²
        double m_e = 0.000511;  // GeV

        // Simplified: τ ~ 1/(G_F² Q⁵)
        return 1.0 / (G_F * G_F * std::pow(Q_value, 5.0));  // s
    }

    /**
     * @brief W/Z boson production at hadron colliders
     *
     * pp → W + X, pp → Z + X
     *
     * σ(W) ~ 200 nb at LHC (√s = 13 TeV)
     * σ(Z) ~ 60 nb at LHC
     */
    static double wBosonProduction() {
        return 200.0e-9 * 1e12;  // nb to mb
    }

    static double zBosonProduction() {
        return 60.0e-9 * 1e12;  // nb to mb
    }
};

/**
 * @class PartonDistributionFunctions
 * @brief PDFs for hadron structure
 *
 * f(x, Q²) = probability to find parton with momentum fraction x
 */
class PartonDistributionFunctions {
public:
    /**
     * @brief Bjorken x (momentum fraction)
     *
     * x = Q²/(2P·q)
     *
     * Range: 0 < x < 1
     */
    static double bjorkenX(double Q2, double P_q) {
        return Q2 / (2.0 * P_q);
    }

    /**
     * @brief Valence quark distribution (simple model)
     *
     * f_valence(x) ~ x^α (1-x)^β
     *
     * Peaks at intermediate x
     */
    static double valenceQuarkPDF(double x, double alpha = 0.5, double beta = 3.0) {
        if (x <= 0.0 || x >= 1.0) return 0.0;

        return std::pow(x, alpha) * std::pow(1.0 - x, beta);
    }

    /**
     * @brief Gluon distribution
     *
     * Dominates at small x
     */
    static double gluonPDF(double x) {
        if (x <= 0.0 || x >= 1.0) return 0.0;

        // Simplified: f_g(x) ~ 1/x at small x
        return 3.0 / std::pow(x, 0.3);  // Rough approximation
    }

    /**
     * @brief Momentum sum rule
     *
     * ∫₀¹ x Σf_i(x) dx = 1
     *
     * All partons carry 100% of proton momentum
     */
    static std::string momentumSumRule() {
        return "∫₀¹ x[Σq(x) + g(x)] dx = 1";
    }
};

} // namespace physics::advanced::qft

#endif // PHYSICS_ADVANCED_QFT_CROSS_SECTIONS_HPP
