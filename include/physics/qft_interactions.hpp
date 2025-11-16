#ifndef PHYSICS_ADVANCED_QFT_INTERACTIONS_HPP
#define PHYSICS_ADVANCED_QFT_INTERACTIONS_HPP

#include <cmath>
#include <string>
#include <complex>

/**
 * @file interactions.hpp
 * @brief Fundamental interactions via boson exchange
 *
 * Implements:
 * - Boson exchange mechanism (Yukawa potential)
 * - Coupling constants (α_EM, α_s, α_W)
 * - Feynman vertices
 * - Running coupling constants
 * - Force ranges
 */

namespace physics::advanced::qft {

/**
 * @class BosonExchange
 * @brief Force mediation by virtual boson exchange
 *
 * Forces arise from exchange of virtual gauge bosons:
 * - EM: photon exchange
 * - Weak: W±, Z⁰ exchange
 * - Strong: gluon exchange
 */
class BosonExchange {
public:
    /**
     * @brief Yukawa potential from massive boson exchange
     *
     * V(r) = -g²/(4π) × (e^(-m_boson × r) / r)
     *
     * For m → 0 (massless): Coulomb potential 1/r
     * For m > 0 (massive): exponential suppression at r > 1/m
     *
     * @param coupling_constant g²/(4π)
     * @param boson_mass m (GeV)
     * @param distance r (m)
     * @param hbar ℏ (GeV·s)
     * @param c speed of light (m/s)
     * @return Potential energy (GeV)
     */
    static double yukawaPotential(double coupling_constant, double boson_mass,
                                 double distance, double hbar = 6.582e-25,
                                 double c = 2.998e8) {
        // Convert to natural units
        double m_r = boson_mass * c / hbar;  // 1/m

        double exponent = -m_r * distance;
        double r_inv = 1.0 / distance;

        return -coupling_constant * std::exp(exponent) * r_inv / (4.0 * M_PI);
    }

    /**
     * @brief Force range from boson mass
     *
     * Range ~ ℏ/(m_boson × c) = λ_Compton
     *
     * Massless (photon): infinite range
     * W/Z bosons: ~10⁻¹⁸ m (nuclear scale)
     */
    static double forceRange(double boson_mass, double hbar = 6.582e-25,
                            double c = 2.998e8) {
        if (boson_mass < 1e-10) {
            return INFINITY;  // Massless boson → infinite range
        }

        return hbar / (boson_mass * c);  // meters
    }

    /**
     * @brief Virtual particle energy-time uncertainty
     *
     * ΔE × Δt ~ ℏ
     *
     * Virtual particles can violate energy conservation
     * for time Δt ~ ℏ/ΔE
     */
    static double virtualParticleLifetime(double energy_violation,
                                         double hbar = 6.582e-25) {
        return hbar / energy_violation;  // seconds
    }
};

/**
 * @class CouplingConstants
 * @brief Fundamental coupling strengths
 *
 * Determine strength of interactions
 */
class CouplingConstants {
public:
    /**
     * @brief Fine structure constant (electromagnetic)
     *
     * α = e²/(4πε₀ℏc) ≈ 1/137.036
     *
     * At low energy (Thompson limit)
     */
    static double fineStructure() {
        return 1.0 / 137.036;  // ≈ 0.00729735
    }

    /**
     * @brief Strong coupling constant
     *
     * α_s(M_Z) ≈ 0.1179
     *
     * At Z boson mass scale (~91 GeV)
     * Much stronger than EM!
     */
    static double strong(double energy_scale = 91.2) {
        // At M_Z scale
        if (std::abs(energy_scale - 91.2) < 1.0) {
            return 0.1179;
        }

        // Running coupling (simplified)
        return runningStrongCoupling(energy_scale);
    }

    /**
     * @brief Weak coupling constant
     *
     * α_W ≈ 1/30 at M_W scale
     *
     * Similar strength to EM at high energies
     */
    static double weak() {
        return 1.0 / 30.0;  // ≈ 0.033
    }

    /**
     * @brief Gravitational coupling (dimensionless)
     *
     * α_G = Gm²/(ℏc) for proton mass
     *
     * Extremely weak! α_G ~ 10⁻³⁹
     */
    static double gravitational() {
        // For proton
        double G = 6.674e-11;       // m³/(kg·s²)
        double m_p = 1.673e-27;     // kg
        double hbar = 1.055e-34;    // J·s
        double c = 2.998e8;         // m/s

        return (G * m_p * m_p) / (hbar * c);  // ~ 5.9 × 10⁻³⁹
    }

    /**
     * @brief Relative strength comparison
     *
     * Strong : EM : Weak : Gravity
     * 1 : 10⁻² : 10⁻⁶ : 10⁻³⁹
     */
    static std::string relativeStrengths() {
        return "Strong : EM : Weak : Gravity = 1 : 10⁻² : 10⁻⁶ : 10⁻³⁹";
    }

private:
    /**
     * @brief Running strong coupling (simplified)
     *
     * α_s(Q²) = α_s(μ²) / [1 + (β₀/2π)α_s(μ²)ln(Q²/μ²)]
     *
     * QCD beta function: β₀ = 11 - (2/3)n_f
     * For n_f = 5 flavors: β₀ = 23/3
     */
    static double runningStrongCoupling(double energy_scale) {
        double alpha_s_MZ = 0.1179;
        double MZ = 91.2;  // GeV
        double beta0 = 23.0/3.0;  // 5 active flavors

        double log_ratio = std::log(energy_scale * energy_scale / (MZ * MZ));
        double denominator = 1.0 + (beta0 / (2.0 * M_PI)) * alpha_s_MZ * log_ratio;

        return alpha_s_MZ / denominator;
    }
};

/**
 * @class RunningCouplings
 * @brief Energy-dependent coupling constants
 *
 * Couplings "run" with energy scale due to quantum corrections
 */
class RunningCouplings {
public:
    /**
     * @brief Running electromagnetic coupling
     *
     * α(Q²) = α(0) / [1 - (α(0)/3π)ln(Q²/m_e²)]
     *
     * Increases with energy ("antiscreening")
     */
    static double runningEM(double energy_scale) {
        double alpha_0 = CouplingConstants::fineStructure();
        double m_e = 0.000511;  // GeV

        double log_term = std::log(energy_scale * energy_scale / (m_e * m_e));
        double denominator = 1.0 - (alpha_0 / (3.0 * M_PI)) * log_term;

        return alpha_0 / denominator;
    }

    /**
     * @brief Running strong coupling
     *
     * Decreases with energy ("asymptotic freedom")
     * Increases at low energy (confinement)
     */
    static double runningStrong(double energy_scale) {
        return CouplingConstants::strong(energy_scale);
    }

    /**
     * @brief QCD Λ parameter (scale where α_s diverges)
     *
     * Λ_QCD ≈ 200 MeV
     *
     * Below this scale, perturbative QCD breaks down
     */
    static double lambdaQCD() {
        return 0.200;  // GeV
    }

    /**
     * @brief Asymptotic freedom
     *
     * α_s → 0 as Q² → ∞
     *
     * Quarks behave as free particles at high energy
     */
    static bool isAsymptoticallyFree(double energy_scale) {
        return runningStrong(energy_scale) < 0.1;  // Perturbative regime
    }

    /**
     * @brief Confinement scale
     *
     * Below ~200 MeV, α_s becomes large
     * Quarks are confined in hadrons
     */
    static bool isConfined(double energy_scale) {
        return energy_scale < lambdaQCD();
    }
};

/**
 * @class FermionBosonCoupling
 * @brief Coupling of fermions to gauge bosons
 */
class FermionBosonCoupling {
public:
    /**
     * @brief QED vertex: fermion-photon coupling
     *
     * Vertex factor: -ieγ^μ
     *
     * where e = √(4πα) is electric charge
     */
    static std::complex<double> qedVertex(double charge) {
        double alpha = CouplingConstants::fineStructure();
        double e = std::sqrt(4.0 * M_PI * alpha);

        return std::complex<double>(0, -charge * e);  // -ieQ
    }

    /**
     * @brief QCD vertex: quark-gluon coupling
     *
     * Vertex factor: -ig_s γ^μ T^a
     *
     * where g_s = √(4πα_s) and T^a are color matrices
     */
    static std::complex<double> qcdVertex(double energy_scale = 91.2) {
        double alpha_s = CouplingConstants::strong(energy_scale);
        double g_s = std::sqrt(4.0 * M_PI * alpha_s);

        return std::complex<double>(0, -g_s);  // -ig_s
    }

    /**
     * @brief Weak charged current: fermion-W± coupling
     *
     * Vertex factor: -ig_W/√2 γ^μ(1-γ⁵)
     *
     * Only couples to left-handed fermions!
     */
    static std::complex<double> weakChargedVertex() {
        double alpha_W = CouplingConstants::weak();
        double g_W = std::sqrt(4.0 * M_PI * alpha_W);

        return std::complex<double>(0, -g_W / std::sqrt(2.0));  // -ig_W/√2
    }

    /**
     * @brief Weak neutral current: fermion-Z coupling
     *
     * Vertex factor: -ig_Z γ^μ(g_V - g_A γ⁵)
     *
     * Vector (g_V) and axial-vector (g_A) couplings
     */
    static std::pair<double, double> weakNeutralCouplings(double charge,
                                                          double isospin) {
        double sin2_theta_W = 0.231;  // Weak mixing angle

        // g_V = T₃ - 2Q sin²θ_W
        double g_V = isospin - 2.0 * charge * sin2_theta_W;

        // g_A = T₃
        double g_A = isospin;

        return {g_V, g_A};
    }

    /**
     * @brief Yukawa coupling: fermion-Higgs
     *
     * y_f = √2 m_f / v
     *
     * where v = 246 GeV is Higgs VEV
     * Explains why particles have different masses
     */
    static double yukawaCoupling(double fermion_mass) {
        double v = 246.0;  // GeV (Higgs VEV)
        return std::sqrt(2.0) * fermion_mass / v;
    }

    /**
     * @brief Top quark Yukawa coupling
     *
     * y_t ≈ 1 (largest Yukawa coupling)
     *
     * m_t ≈ 173 GeV → y_t ≈ 173√2/246 ≈ 1
     */
    static double topYukawa() {
        return yukawaCoupling(173.0);  // ≈ 1.0
    }
};

/**
 * @class FeynmanRules
 * @brief Feynman rules for QFT calculations
 */
class FeynmanRules {
public:
    /**
     * @brief Photon propagator
     *
     * D^μν(q²) = -ig^μν/q²
     *
     * For massless photon
     */
    static std::string photonPropagator() {
        return "D^μν(q²) = -ig^μν/q² (Feynman gauge)";
    }

    /**
     * @brief Massive boson propagator (W, Z, Higgs)
     *
     * D^μν(q²) = -i[g^μν - q^μq^ν/M²]/(q² - M²)
     *
     * Pole at q² = M² (on-shell)
     */
    static std::string massivePropagator() {
        return "D^μν(q²) = -i[g^μν - q^μq^ν/M²]/(q² - M²)";
    }

    /**
     * @brief Fermion propagator
     *
     * S_F(p) = i(γ·p + m)/(p² - m²)
     *
     * Dirac propagator
     */
    static std::string fermionPropagator() {
        return "S_F(p) = i(γ·p + m)/(p² - m²)";
    }

    /**
     * @brief QED vertex factor
     */
    static std::string qedVertexRule() {
        return "Vertex: -ieγ^μ (charge e at photon-fermion vertex)";
    }

    /**
     * @brief Symmetry factor
     *
     * 1/n! for n identical particles in final state
     * 1/2 for identical internal lines
     */
    static double symmetryFactor(int num_identical) {
        double factor = 1.0;
        for (int i = 2; i <= num_identical; ++i) {
            factor *= i;
        }
        return 1.0 / factor;
    }
};

/**
 * @class InteractionRanges
 * @brief Characteristic ranges of forces
 */
class InteractionRanges {
public:
    /**
     * @brief Electromagnetic range
     *
     * Infinite (photon massless)
     */
    static double electromagnetic() {
        return INFINITY;  // meters
    }

    /**
     * @brief Weak interaction range
     *
     * ~ ℏ/(M_W c) ≈ 10⁻¹⁸ m
     */
    static double weak() {
        double M_W = 80.4;  // GeV
        return BosonExchange::forceRange(M_W);  // ~10⁻¹⁸ m
    }

    /**
     * @brief Strong interaction range
     *
     * ~ 1 fm = 10⁻¹⁵ m (size of nucleon)
     *
     * Gluons massless but confined
     */
    static double strong() {
        return 1.0e-15;  // 1 fermi
    }

    /**
     * @brief Gravitational range
     *
     * Infinite (graviton massless)
     */
    static double gravitational() {
        return INFINITY;  // meters
    }
};

} // namespace physics::advanced::qft

#endif // PHYSICS_ADVANCED_QFT_INTERACTIONS_HPP
