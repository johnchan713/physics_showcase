#ifndef PHYSICS_ADVANCED_GAUGE_THEORY_HELICITY_HPP
#define PHYSICS_ADVANCED_GAUGE_THEORY_HELICITY_HPP

#include <cmath>
#include <string>
#include <complex>

/**
 * @file helicity.hpp
 * @brief Helicity and chirality in particle physics
 *
 * Implements:
 * - Helicity (spin projection along momentum)
 * - Chirality (handedness, Weyl spinors)
 * - Helicity conservation in massless limit
 * - Relation between helicity and chirality
 * - V-A structure of weak interactions
 */

namespace physics::advanced::gauge_theory {

/**
 * @class Helicity
 * @brief Helicity: projection of spin along momentum
 *
 * h = S · p̂ = ±1/2 (for spin-1/2)
 *
 * Right-handed: h = +1/2 (spin parallel to momentum)
 * Left-handed: h = -1/2 (spin antiparallel to momentum)
 */
class Helicity {
public:
    /**
     * @brief Helicity definition
     *
     * h = S · p̂ / |S|
     *
     * For spin-1/2: h = ±1/2
     * For spin-1: h = ±1, 0
     */
    static std::string definition() {
        return "Helicity: h = S · p̂\n"
               "Projection of spin along momentum direction\n"
               "\n"
               "Spin-1/2: h = ±1/2 (R or L)\n"
               "Spin-1: h = ±1, 0";
    }

    /**
     * @brief Right-handed vs left-handed
     *
     * Right-handed (R): h = +1/2 (spin ↑↑ momentum)
     * Left-handed (L): h = -1/2 (spin ↑↓ momentum)
     */
    static double rightHanded() {
        return +0.5;  // h = +1/2
    }

    static double leftHanded() {
        return -0.5;  // h = -1/2
    }

    /**
     * @brief Helicity operator
     *
     * ĥ = (Σ · p̂) / 2
     *
     * where Σ are Pauli matrices
     */
    static std::string operator_form() {
        return "Helicity operator: ĥ = (Σ · p̂)/2\n"
               "Eigenvalues: ±1/2";
    }

    /**
     * @brief Helicity is frame-dependent (unless m = 0)
     *
     * For massive particles, can boost to frame where helicity flips
     * For massless particles, helicity is Lorentz invariant
     */
    static bool isLorentzInvariant(double mass) {
        return (mass == 0.0);  // Only for massless particles
    }

    static std::string frameDependendence() {
        return "Massive particles: helicity frame-dependent\n"
               "Can boost to rest frame, then boost past → helicity flips\n"
               "\n"
               "Massless particles: helicity Lorentz invariant\n"
               "Cannot boost past massless particle (v = c always)";
    }

    /**
     * @brief Helicity for photons
     *
     * Photon (spin-1, massless): h = ±1
     * h = +1: right circular polarization
     * h = -1: left circular polarization
     * h = 0: forbidden (photon is massless)
     */
    static std::string photonHelicity() {
        return "Photon helicity: h = ±1\n"
               "h = +1: right circular polarization (R)\n"
               "h = -1: left circular polarization (L)\n"
               "h = 0: forbidden (massless gauge boson)";
    }

    /**
     * @brief Helicity for gluons
     *
     * Same as photons: h = ±1
     */
    static std::string gluonHelicity() {
        return "Gluon helicity: h = ±1\n"
               "Same as photon (massless spin-1)";
    }

    /**
     * @brief Helicity for gravitons
     *
     * Graviton (spin-2, massless): h = ±2
     */
    static std::string gravitonHelicity() {
        return "Graviton helicity: h = ±2\n"
               "Massless spin-2 particle";
    }
};

/**
 * @class Chirality
 * @brief Chirality: handedness (left/right Weyl spinors)
 *
 * γ⁵ eigenvalue: ±1
 *
 * Left-chiral: ψ_L = P_L ψ, where P_L = (1 - γ⁵)/2
 * Right-chiral: ψ_R = P_R ψ, where P_R = (1 + γ⁵)/2
 */
class Chirality {
public:
    /**
     * @brief Chirality operator γ⁵
     *
     * γ⁵ = iγ⁰γ¹γ²γ³
     *
     * Eigenvalues: ±1
     * {γ⁵, γ^μ} = 0 (anticommutes with Dirac matrices)
     */
    static std::string gammaFiveOperator() {
        return "Chirality operator: γ⁵ = iγ⁰γ¹γ²γ³\n"
               "Eigenvalues: ±1\n"
               "{γ⁵, γ^μ} = 0 (anticommutes)";
    }

    /**
     * @brief Chiral projection operators
     *
     * P_L = (1 - γ⁵)/2  (left projector)
     * P_R = (1 + γ⁵)/2  (right projector)
     *
     * P_L + P_R = 1
     * P_L² = P_L, P_R² = P_R (idempotent)
     * P_L P_R = 0 (orthogonal)
     */
    static std::string projectionOperators() {
        return "Left projector: P_L = (1 - γ⁵)/2\n"
               "Right projector: P_R = (1 + γ⁵)/2\n"
               "\n"
               "Properties:\n"
               "P_L + P_R = 1\n"
               "P_L² = P_L, P_R² = P_R\n"
               "P_L P_R = 0";
    }

    /**
     * @brief Weyl spinors (chiral fermions)
     *
     * ψ_L = P_L ψ (left-chiral)
     * ψ_R = P_R ψ (right-chiral)
     *
     * ψ = ψ_L + ψ_R (Dirac spinor)
     */
    static std::string weylSpinors() {
        return "Weyl spinors (2-component):\n"
               "ψ_L = P_L ψ (left-chiral)\n"
               "ψ_R = P_R ψ (right-chiral)\n"
               "\n"
               "Dirac spinor: ψ = ψ_L + ψ_R (4-component)";
    }

    /**
     * @brief Chirality is Lorentz invariant
     *
     * Unlike helicity, chirality does NOT change under boosts
     * γ⁵ commutes with Lorentz transformations
     */
    static bool isLorentzInvariant() {
        return true;  // Always Lorentz invariant
    }

    /**
     * @brief Massless limit: chirality = helicity
     *
     * For m = 0: chirality eigenstates = helicity eigenstates
     * ψ_L ↔ h = -1/2 (left-handed)
     * ψ_R ↔ h = +1/2 (right-handed)
     */
    static std::string masslessLimit() {
        return "Massless fermions: chirality = helicity\n"
               "ψ_L ↔ h = -1/2 (left-handed)\n"
               "ψ_R ↔ h = +1/2 (right-handed)\n"
               "\n"
               "Example: neutrinos (nearly massless) are left-handed";
    }

    /**
     * @brief Massive fermions: chirality ≠ helicity
     *
     * For m ≠ 0: chirality and helicity are different
     * Mass term mixes ψ_L and ψ_R
     */
    static std::string massiveCase() {
        return "Massive fermions: chirality ≠ helicity\n"
               "Mass term: m(ψ̄_L ψ_R + ψ̄_R ψ_L) mixes chiralities\n"
               "Helicity flips in different frames";
    }

    /**
     * @brief Chiral basis
     *
     * Dirac equation in chiral basis:
     * iσ^μ ∂_μ ψ_R = m ψ_L
     * iσ̄^μ ∂_μ ψ_L = m ψ_R
     *
     * where σ^μ = (1, σⁱ), σ̄^μ = (1, -σⁱ)
     */
    static std::string chiralDiracEquation() {
        return "Dirac equation (chiral basis):\n"
               "iσ^μ ∂_μ ψ_R = m ψ_L\n"
               "iσ̄^μ ∂_μ ψ_L = m ψ_R\n"
               "\n"
               "If m = 0: ψ_L and ψ_R decouple (massless Weyl fermions)";
    }
};

/**
 * @class HelicityConservation
 * @brief Helicity conservation in massless limit
 *
 * For massless particles, helicity is conserved in all interactions
 * (except when mass is generated)
 */
class HelicityConservation {
public:
    /**
     * @brief Massless fermions: helicity conserved
     *
     * Helicity is good quantum number for m = 0
     * QED vertex conserves helicity for massless fermions
     */
    static bool conservedForMassless() {
        return true;
    }

    /**
     * @brief QED vertex for massless fermions
     *
     * Vertex: ψ̄ γ^μ ψ
     *
     * For massless: ψ̄_L γ^μ ψ_L and ψ̄_R γ^μ ψ_R
     * Helicity conserved at each vertex!
     */
    static std::string qedVertex() {
        return "QED vertex: ψ̄ γ^μ ψ\n"
               "\n"
               "Massless fermions:\n"
               "ψ̄_L γ^μ ψ_L (left stays left)\n"
               "ψ̄_R γ^μ ψ_R (right stays right)\n"
               "ψ̄_L γ^μ ψ_R = 0 (orthogonal)\n"
               "\n"
               "Helicity conserved!";
    }

    /**
     * @brief Massive fermions: helicity violated by mass
     *
     * Mass term: m ψ̄ ψ = m(ψ̄_L ψ_R + ψ̄_R ψ_L)
     *
     * Mixes left and right, violates helicity conservation
     */
    static bool conservedForMassive() {
        return false;  // Mass violates helicity
    }

    static std::string massViolation() {
        return "Mass term: m ψ̄ ψ = m(ψ̄_L ψ_R + ψ̄_R ψ_L)\n"
               "\n"
               "Mixes chiralities → violates helicity conservation\n"
               "Helicity flip rate ∝ m/E (suppressed at high energy)";
    }

    /**
     * @brief Helicity flip rate
     *
     * Probability ~ (m/E)²
     *
     * Suppressed for relativistic particles (E >> m)
     */
    static double helicityFlipRate(double mass, double energy) {
        return (mass / energy) * (mass / energy);  // ~ (m/E)²
    }

    /**
     * @brief Example: neutrinos
     *
     * Neutrinos nearly massless → nearly pure left-handed
     * m_ν ~ 0.1 eV, E_ν ~ MeV-GeV
     * Helicity flip ~ (10⁻⁷)² ~ 10⁻¹⁴ (negligible!)
     */
    static std::string neutrinoExample() {
        return "Neutrinos: m_ν ~ 0.1 eV, E_ν ~ MeV-GeV\n"
               "Helicity flip rate ~ (m/E)² ~ 10⁻¹⁴\n"
               "\n"
               "Neutrinos are effectively pure left-handed!\n"
               "Antineutrinos are pure right-handed";
    }
};

/**
 * @class WeakInteractionChirality
 * @brief V-A structure: weak interactions couple only to left-handed fermions
 *
 * W± and Z⁰ couple only to ψ_L
 * Maximal parity violation!
 */
class WeakInteractionChirality {
public:
    /**
     * @brief V-A theory (Vector - Axial)
     *
     * Weak current: J^μ = ψ̄ γ^μ (1 - γ⁵) ψ
     *             = ψ̄ γ^μ P_L ψ
     *             = 2 ψ̄_L γ^μ ψ_L
     *
     * Only left-handed fermions couple!
     */
    static std::string vMinusAStructure() {
        return "V-A structure (Feynman & Gell-Mann, 1958):\n"
               "\n"
               "Weak current: J^μ = ψ̄ γ^μ (1 - γ⁵) ψ\n"
               "             = 2 ψ̄_L γ^μ ψ_L\n"
               "\n"
               "Only left-handed fermions couple to W±, Z⁰!\n"
               "Maximal parity violation";
    }

    /**
     * @brief Charged current (W±)
     *
     * W⁺: ψ̄_L^up γ^μ ψ_L^down
     * W⁻: ψ̄_L^down γ^μ ψ_L^up
     *
     * Examples:
     * e_L⁻ → ν_eL + W⁻
     * u_L → d_L + W⁺
     */
    static std::string chargedCurrent() {
        return "Charged current (W± coupling):\n"
               "\n"
               "W⁺: u_L ↔ d_L, ν_eL ↔ e_L⁻ (left-handed only)\n"
               "W⁻: d_L ↔ u_L, e_L⁻ ↔ ν_eL\n"
               "\n"
               "Right-handed fermions do NOT couple to W±!";
    }

    /**
     * @brief Neutral current (Z⁰)
     *
     * Z⁰: ψ̄ γ^μ (c_V - c_A γ⁵) ψ
     *
     * Left: c_L = c_V - c_A
     * Right: c_R = c_V + c_A
     *
     * Different couplings for L and R (but both nonzero)
     */
    static std::string neutralCurrent() {
        return "Neutral current (Z⁰ coupling):\n"
               "\n"
               "Z⁰: ψ̄ γ^μ (c_V - c_A γ⁵) ψ\n"
               "\n"
               "Left-handed: c_L = c_V - c_A\n"
               "Right-handed: c_R = c_V + c_A\n"
               "\n"
               "Both L and R couple, but with different strengths";
    }

    /**
     * @brief Neutrino coupling
     *
     * Neutrinos: only left-handed ν_L exists
     * Antineutrinos: only right-handed ν̄_R exists
     *
     * Massless: pure chirality eigenstates
     */
    static std::string neutrinoCoupling() {
        return "Neutrinos in weak interactions:\n"
               "\n"
               "Only ν_L couples to W± (left-handed neutrino)\n"
               "Only ν̄_R couples to W± (right-handed antineutrino)\n"
               "\n"
               "ν_R does not exist in Standard Model!\n"
               "(If neutrinos have mass, tiny ν_R component possible)";
    }

    /**
     * @brief Consequences for beta decay
     *
     * n → p + e⁻ + ν̄_e
     *
     * Emitted electron is left-handed (mostly)
     * Emitted antineutrino is right-handed (always)
     */
    static std::string betaDecay() {
        return "Beta decay: n → p + e⁻ + ν̄_e\n"
               "\n"
               "W⁻ couples to (e⁻)_L and (ν̄_e)_R\n"
               "\n"
               "Electron: mostly left-handed (h = -1/2)\n"
               "Antineutrino: purely right-handed (h = +1/2)";
    }

    /**
     * @brief Parity violation from chirality
     *
     * Left-handed coupling → parity violation
     * P: ψ_L ↔ ψ_R
     * Weak interaction distinguishes L from R!
     */
    static std::string parityViolation() {
        return "Weak chirality → parity violation:\n"
               "\n"
               "Parity: P(ψ_L) = ψ_R\n"
               "Weak couples to ψ_L only\n"
               "\n"
               "P(weak interaction) ≠ weak interaction\n"
               "Maximal parity violation!";
    }
};

/**
 * @class HelicityAmplitudes
 * @brief Helicity amplitudes in scattering
 *
 * Scattering amplitudes depend on helicity states
 * Simplified in high-energy limit (massless approximation)
 */
class HelicityAmplitudes {
public:
    /**
     * @brief Spinor helicity formalism
     *
     * Massless spinors labeled by momentum and helicity
     * |p, +⟩ (positive helicity)
     * |p, -⟩ (negative helicity)
     */
    static std::string spinorHelicityFormalism() {
        return "Spinor helicity formalism:\n"
               "\n"
               "Massless spinors: |p, ±⟩\n"
               "Lorentz invariant products: ⟨p q⟩, [p q]\n"
               "\n"
               "Simplifies high-energy scattering calculations\n"
               "Used extensively in QCD and collider physics";
    }

    /**
     * @brief Conservation rules
     *
     * In massless limit:
     * - Helicity conserved at each vertex
     * - Simplifies amplitude structure
     */
    static std::string conservationRules() {
        return "Massless limit conservation:\n"
               "\n"
               "QED vertex: helicity conserved\n"
               "QCD vertex: helicity conserved (for massless quarks)\n"
               "Weak vertex: couples to definite chirality\n"
               "\n"
               "Allows powerful simplifications in amplitudes";
    }

    /**
     * @brief Example: e⁺e⁻ → μ⁺μ⁻
     *
     * At high energy (E >> m_e, m_μ):
     * Only certain helicity combinations contribute
     */
    static std::string eeToMuMuExample() {
        return "e⁺e⁻ → γ → μ⁺μ⁻ (QED)\n"
               "\n"
               "High energy (E >> m):\n"
               "e_L⁻ e_R⁺ → μ_L⁻ μ_R⁺ (allowed)\n"
               "e_R⁻ e_L⁺ → μ_R⁻ μ_L⁺ (allowed)\n"
               "e_L⁻ e_L⁺ → suppressed by (m/E)\n"
               "e_R⁻ e_R⁺ → suppressed by (m/E)\n"
               "\n"
               "Helicity conservation simplifies calculation!";
    }
};

/**
 * @class ChiralAnomalies
 * @brief Chiral anomalies: quantum violation of classical symmetry
 *
 * Triangle diagrams can violate classical chiral symmetry
 * Important for π⁰ → γγ decay, axial U(1) problem
 */
class ChiralAnomalies {
public:
    /**
     * @brief Axial current anomaly
     *
     * Classical: ∂_μ J_A^μ = 0 (chiral symmetry)
     * Quantum: ∂_μ J_A^μ = (g²/16π²) ε^μνρσ F_μν F_ρσ (anomaly!)
     *
     * Triangle diagram with gauge bosons violates chiral conservation
     */
    static std::string axialAnomaly() {
        return "Axial (chiral) anomaly:\n"
               "\n"
               "Classical chiral symmetry: ∂_μ J_A^μ = 0\n"
               "Quantum anomaly: ∂_μ J_A^μ = (g²/16π²) F·F̃\n"
               "\n"
               "Triangle diagram violates classical conservation!\n"
               "Adler-Bell-Jackiw anomaly (1969)";
    }

    /**
     * @brief π⁰ → γγ decay
     *
     * Forbidden by classical chiral symmetry
     * Allowed by anomaly!
     *
     * Γ(π⁰ → γγ) ~ (α/π)² (determined by anomaly)
     */
    static std::string pionDecay() {
        return "π⁰ → γγ decay:\n"
               "\n"
               "Classically forbidden (chiral symmetry)\n"
               "Allowed by axial anomaly\n"
               "\n"
               "Decay rate: Γ ~ (α/π)² m_π³/f_π²\n"
               "Lifetime: τ ~ 10⁻¹⁶ s\n"
               "\n"
               "Anomaly explains why neutral pion decays!";
    }

    /**
     * @brief Anomaly cancellation in SM
     *
     * Gauge anomalies must cancel (consistency)
     * SM: anomalies cancel between quarks and leptons!
     *
     * Each generation: Tr[Q³] = 0 (miraculous!)
     */
    static std::string anomalyCancellation() {
        return "Anomaly cancellation in Standard Model:\n"
               "\n"
               "Gauge anomalies must vanish (consistency)\n"
               "Quarks: Σ Q³ = 3(2³ - 1³ - 1³)/3³ = 1/3\n"
               "Leptons: Σ Q³ = -1³ = -1\n"
               "Wait... that doesn't cancel!\n"
               "\n"
               "With 3 colors for quarks:\n"
               "3 × 1/3 - 1 = 0 ✓\n"
               "\n"
               "Miraculous cancellation! Why 3 colors?";
    }

    /**
     * @brief Global vs gauge anomalies
     *
     * Global anomaly: OK (e.g., axial U(1), explains π⁰ → γγ)
     * Gauge anomaly: FORBIDDEN (inconsistency)
     */
    static std::string globalVsGauge() {
        return "Anomalies in symmetries:\n"
               "\n"
               "Global symmetry anomaly: Allowed\n"
               "  Example: Axial U(1) → π⁰ decay\n"
               "\n"
               "Gauge symmetry anomaly: FORBIDDEN\n"
               "  Would make theory inconsistent\n"
               "  Must cancel (as in SM)";
    }
};

} // namespace physics::advanced::gauge_theory

#endif // PHYSICS_ADVANCED_GAUGE_THEORY_HELICITY_HPP
