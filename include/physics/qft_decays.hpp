#ifndef PHYSICS_ADVANCED_QFT_DECAYS_HPP
#define PHYSICS_ADVANCED_QFT_DECAYS_HPP

#include <cmath>
#include <string>
#include <vector>

/**
 * @file decays.hpp
 * @brief Particle decays and resonances
 *
 * Implements:
 * - Decay rates and lifetimes
 * - Branching ratios
 * - Resonances (Breit-Wigner)
 * - Golden rule calculations
 * - Phase space factors
 */

namespace physics::advanced::qft {

/**
 * @class DecayRate
 * @brief Fundamental decay calculations
 *
 * Γ (decay width) has units of energy
 * τ (lifetime) = ℏ/Γ
 */
class DecayRate {
public:
    /**
     * @brief Lifetime from decay width
     *
     * τ = ℏ/Γ
     *
     * @param decay_width Γ (GeV)
     * @param hbar ℏ (GeV·s)
     * @return Lifetime (s)
     */
    static double lifetimeFromWidth(double decay_width, double hbar = 6.582e-25) {
        if (decay_width <= 0.0) return INFINITY;  // Stable
        return hbar / decay_width;
    }

    /**
     * @brief Decay width from lifetime
     *
     * Γ = ℏ/τ
     */
    static double widthFromLifetime(double lifetime, double hbar = 6.582e-25) {
        if (lifetime == INFINITY || lifetime == 0.0) return 0.0;  // Stable
        return hbar / lifetime;
    }

    /**
     * @brief Exponential decay law
     *
     * N(t) = N₀ exp(-t/τ) = N₀ exp(-Γt/ℏ)
     *
     * @param initial_number N₀
     * @param time t (s)
     * @param lifetime τ (s)
     * @return N(t)
     */
    static double exponentialDecay(double initial_number, double time,
                                   double lifetime) {
        return initial_number * std::exp(-time / lifetime);
    }

    /**
     * @brief Half-life
     *
     * t₁/₂ = τ ln(2) ≈ 0.693 τ
     */
    static double halfLife(double lifetime) {
        return lifetime * std::log(2.0);
    }

    /**
     * @brief Mean decay length
     *
     * L = βγcτ
     *
     * Distance traveled before decaying (relativistic)
     */
    static double decayLength(double beta, double gamma, double lifetime,
                             double c = 2.998e8) {
        return beta * gamma * c * lifetime;  // meters
    }
};

/**
 * @class BranchingRatio
 * @brief Partial decay rates and branching fractions
 */
class BranchingRatio {
public:
    /**
     * @brief Branching ratio for channel i
     *
     * BR_i = Γ_i / Γ_total
     *
     * @param partial_width Γ_i (GeV)
     * @param total_width Γ_total (GeV)
     * @return Branching fraction (0 to 1)
     */
    static double branchingFraction(double partial_width, double total_width) {
        if (total_width <= 0.0) return 0.0;
        return partial_width / total_width;
    }

    /**
     * @brief Total width from partial widths
     *
     * Γ_total = Σ Γ_i
     */
    static double totalWidth(const std::vector<double>& partial_widths) {
        double total = 0.0;
        for (double width : partial_widths) {
            total += width;
        }
        return total;
    }

    /**
     * @brief Check branching ratios sum to 1
     */
    static bool checkNormalization(const std::vector<double>& branching_ratios,
                                   double tolerance = 1e-6) {
        double sum = 0.0;
        for (double br : branching_ratios) {
            sum += br;
        }
        return std::abs(sum - 1.0) < tolerance;
    }
};

/**
 * @class FermisGoldenRule
 * @brief Fermi's golden rule for decay rates
 *
 * Γ = (2π/ℏ)|M|² ρ(E)
 *
 * where |M|² is matrix element squared, ρ(E) is density of states
 */
class FermisGoldenRule {
public:
    /**
     * @brief Two-body decay rate
     *
     * Γ = (1/8πm²)|M|²|p*|
     *
     * where |p*| is momentum in rest frame
     *
     * @param matrix_element_squared |M|² (GeV²)
     * @param parent_mass m (GeV)
     * @param momentum_cm |p*| (GeV)
     * @return Decay width (GeV)
     */
    static double twoBodyRate(double matrix_element_squared, double parent_mass,
                             double momentum_cm) {
        return (1.0 / (8.0 * M_PI * parent_mass * parent_mass)) *
               matrix_element_squared * momentum_cm;
    }

    /**
     * @brief Three-body phase space suppression
     *
     * Much smaller than two-body due to phase space
     */
    static double threeBodyRate(double matrix_element_squared, double parent_mass) {
        // Simplified: Γ_3 ~ (1/(2π)³) × M² × m⁵
        return matrix_element_squared * std::pow(parent_mass, 5.0) /
               std::pow(2.0 * M_PI, 3.0);
    }

    /**
     * @brief Phase space factor for two-body decay
     *
     * |p*| = √[λ(m²,m₁²,m₂²)] / (2m)
     *
     * where λ(a,b,c) = a² + b² + c² - 2ab - 2ac - 2bc (Källén function)
     */
    static double twoBodyMomentum(double parent_mass, double daughter1_mass,
                                  double daughter2_mass) {
        double m2 = parent_mass * parent_mass;
        double m1_2 = daughter1_mass * daughter1_mass;
        double m2_2 = daughter2_mass * daughter2_mass;

        // Källén function
        double lambda = m2 * m2 + m1_2 * m1_2 + m2_2 * m2_2 -
                       2.0 * (m2 * m1_2 + m2 * m2_2 + m1_2 * m2_2);

        if (lambda < 0.0) return 0.0;  // Kinematically forbidden

        return std::sqrt(lambda) / (2.0 * parent_mass);
    }
};

/**
 * @class Resonances
 * @brief Unstable particles as resonances
 */
class Resonances {
public:
    /**
     * @brief Breit-Wigner distribution
     *
     * σ(E) ∝ Γ² / [(E - M)² + Γ²/4]
     *
     * Resonance at E = M with width Γ
     *
     * @param energy E (GeV)
     * @param resonance_mass M (GeV)
     * @param width Γ (GeV)
     * @return Relative cross section
     */
    static double breitWigner(double energy, double resonance_mass, double width) {
        double delta_E = energy - resonance_mass;
        double denominator = delta_E * delta_E + width * width / 4.0;

        return (width * width / 4.0) / denominator;
    }

    /**
     * @brief Relativistic Breit-Wigner
     *
     * More accurate for vector mesons
     */
    static double relativisticBreitWigner(double s, double mass, double width) {
        double m2 = mass * mass;
        double denominator = (s - m2) * (s - m2) + m2 * width * width;

        return (m2 * width * width) / denominator;
    }

    /**
     * @brief Full width at half maximum (FWHM)
     *
     * FWHM = Γ (directly the decay width)
     */
    static double fwhm(double width) {
        return width;
    }

    /**
     * @brief Resonance lifetime
     *
     * τ = ℏ/Γ
     */
    static double resonanceLifetime(double width, double hbar = 6.582e-25) {
        return DecayRate::lifetimeFromWidth(width, hbar);
    }

    /**
     * @brief Examples of resonances
     */
    static std::string examples() {
        return
            "Z⁰ boson: M = 91.2 GeV, Γ = 2.5 GeV → τ ~ 3×10⁻²⁵ s\n"
            "ρ meson: M = 770 MeV, Γ = 150 MeV → τ ~ 4×10⁻²⁴ s\n"
            "Δ(1232): M = 1232 MeV, Γ = 120 MeV → τ ~ 5×10⁻²⁴ s\n"
            "W± boson: M = 80.4 GeV, Γ = 2.1 GeV → τ ~ 3×10⁻²⁵ s";
    }
};

/**
 * @class WeakDecays
 * @brief Weak interaction decays
 */
class WeakDecays {
public:
    /**
     * @brief Muon decay: μ⁻ → e⁻ + ν̄ₑ + νμ
     *
     * Γ(μ) = G_F² m_μ⁵ / (192π³)
     *
     * Lifetime: τ = 2.2 μs
     */
    static double muonDecayWidth() {
        double G_F = 1.166e-5;  // GeV⁻²
        double m_mu = 0.1057;   // GeV

        return (G_F * G_F * std::pow(m_mu, 5.0)) / (192.0 * std::pow(M_PI, 3.0));
    }

    static double muonLifetime() {
        return DecayRate::lifetimeFromWidth(muonDecayWidth());  // ~2.2 μs
    }

    /**
     * @brief Neutron beta decay: n → p + e⁻ + ν̄ₑ
     *
     * τ_n ≈ 880 s
     */
    static double neutronLifetime() {
        return 880.0;  // seconds
    }

    /**
     * @brief Pion decay: π⁺ → μ⁺ + νμ
     *
     * τ_π ≈ 26 ns
     */
    static double pionLifetime() {
        return 26.0e-9;  // seconds
    }

    /**
     * @brief Kaon decay: K⁺ → μ⁺ + νμ (and other modes)
     *
     * τ_K ≈ 12 ns
     */
    static double kaonLifetime() {
        return 12.0e-9;  // seconds
    }

    /**
     * @brief CKM matrix element magnitude
     *
     * |V_ud| ≈ 0.974 (d → u transition)
     * |V_us| ≈ 0.225 (s → u transition)
     * |V_ub| ≈ 0.004 (b → u transition)
     */
    static double ckmElement(int up_type, int down_type) {
        // Simplified approximate values
        if (up_type == 1 && down_type == 1) return 0.974;  // V_ud
        if (up_type == 1 && down_type == 2) return 0.225;  // V_us
        if (up_type == 1 && down_type == 3) return 0.004;  // V_ub
        if (up_type == 2 && down_type == 1) return 0.225;  // V_cd
        if (up_type == 2 && down_type == 2) return 0.973;  // V_cs
        if (up_type == 2 && down_type == 3) return 0.041;  // V_cb
        return 0.0;
    }
};

/**
 * @class StrongDecays
 * @brief Strong interaction decays
 */
class StrongDecays {
public:
    /**
     * @brief ρ meson decay: ρ → π + π
     *
     * Γ_ρ ≈ 150 MeV → τ ~ 4×10⁻²⁴ s
     */
    static double rhoMesonWidth() {
        return 0.150;  // GeV
    }

    /**
     * @brief Δ resonance: Δ⁺⁺ → p + π⁺
     *
     * Γ_Δ ≈ 120 MeV
     */
    static double deltaResonanceWidth() {
        return 0.120;  // GeV
    }

    /**
     * @brief Strong decays are very fast
     *
     * τ ~ 10⁻²³ - 10⁻²⁴ s (hadronic timescale)
     */
    static double typicalLifetime() {
        return 1e-23;  // seconds
    }
};

/**
 * @class ElectromagneticDecays
 * @brief EM decays (slower than strong, faster than weak)
 */
class ElectromagneticDecays {
public:
    /**
     * @brief π⁰ decay: π⁰ → γ + γ
     *
     * τ_π⁰ ≈ 8.4×10⁻¹⁷ s
     */
    static double neutralPionLifetime() {
        return 8.4e-17;  // seconds
    }

    /**
     * @brief η meson decay: η → γ + γ
     *
     * τ_η ≈ 5×10⁻¹⁹ s
     */
    static double etaMesonLifetime() {
        return 5e-19;  // seconds
    }

    /**
     * @brief EM decays scale as α²
     *
     * Γ_EM ~ α² (compared to strong ~ α_s²)
     */
    static double emSuppressionFactor() {
        double alpha = 1.0/137.0;
        return alpha * alpha / 1.0;  // ~10⁻⁵ compared to strong
    }
};

} // namespace physics::advanced::qft

#endif // PHYSICS_ADVANCED_QFT_DECAYS_HPP
