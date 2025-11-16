#ifndef PHYSICS_ADVANCED_GAUGE_THEORY_SYMMETRIES_HPP
#define PHYSICS_ADVANCED_GAUGE_THEORY_SYMMETRIES_HPP

#include <string>
#include <cmath>
#include <complex>

/**
 * @file symmetries.hpp
 * @brief Discrete symmetries: Parity (P), Charge conjugation (C), Time reversal (T), CPT theorem
 *
 * Implements:
 * - Parity transformation and violation in weak interactions
 * - Charge conjugation
 * - Time reversal symmetry
 * - CPT theorem (always conserved)
 * - CP and T violation
 */

namespace physics::advanced::gauge_theory {

/**
 * @class ParityTransformation
 * @brief Parity (P): spatial inversion r → -r
 *
 * P|ψ(r, t)⟩ = |ψ(-r, t)⟩
 *
 * Conserved: strong, electromagnetic
 * Violated: weak interactions
 */
class ParityTransformation {
public:
    /**
     * @brief Parity operator on position
     *
     * P: r → -r
     */
    static std::string positionTransformation() {
        return "P: r → -r (spatial inversion)";
    }

    /**
     * @brief Parity transformation of vector vs pseudovector
     *
     * Vector: P(v) = -v (e.g., momentum p)
     * Pseudovector: P(a) = +a (e.g., angular momentum L)
     */
    static std::string vectorTransformation() {
        return "Vector (polar): P(v) = -v (e.g., p, E)\n"
               "Pseudovector (axial): P(a) = +a (e.g., L, B)";
    }

    /**
     * @brief Intrinsic parity of particles
     *
     * Fermions: P = ±1
     * Bosons: P = ±1
     * Fermion-antifermion: P = (-1)^(L+1)
     */
    static int fermionParity() {
        return +1;  // Conventionally +1 for fermions
    }

    static int antifermionParity() {
        return -1;  // -1 for antifermions
    }

    static int photonParity() {
        return -1;  // Photon has odd parity
    }

    /**
     * @brief Parity in strong and EM interactions
     *
     * Strong: Parity conserved
     * EM: Parity conserved
     */
    static bool conservedInStrong() {
        return true;
    }

    static bool conservedInEM() {
        return true;
    }

    /**
     * @brief Parity violation in weak interactions
     *
     * Weak: Parity maximally violated!
     * Only left-handed particles couple to W±, Z⁰
     *
     * Wu experiment (1956): ⁶⁰Co β-decay shows parity violation
     */
    static bool conservedInWeak() {
        return false;  // Violated!
    }

    static std::string wuExperiment() {
        return "Wu experiment (1956): ⁶⁰Co β-decay\n"
               "Electrons emitted preferentially opposite to nuclear spin\n"
               "Proved parity violation in weak interactions";
    }

    /**
     * @brief Lee-Yang prediction (1956, Nobel Prize 1957)
     *
     * Parity is violated in weak interactions
     */
    static std::string leeYangPrediction() {
        return "Lee & Yang (1956): Suggested parity violation to explain θ-τ puzzle\n"
               "Nobel Prize 1957 for parity violation discovery";
    }
};

/**
 * @class ChargeConjugation
 * @brief Charge conjugation (C): particle ↔ antiparticle
 *
 * C|particle⟩ = |antiparticle⟩
 *
 * Swaps: charge, baryon number, lepton number
 * Keeps: mass, spin, lifetime
 */
class ChargeConjugation {
public:
    /**
     * @brief Charge conjugation transformation
     *
     * C: e⁻ ↔ e⁺, p ↔ p̄, ν_e ↔ ν̄_e
     */
    static std::string transformation() {
        return "C: particle ↔ antiparticle\n"
               "Charge: Q → -Q\n"
               "Baryon number: B → -B\n"
               "Lepton number: L → -L";
    }

    /**
     * @brief C eigenvalues for photons
     *
     * C(γ) = -1 (photon is C-odd)
     * π⁰ → γγ requires C(π⁰) = (-1)²  = +1
     */
    static int photonEigenvalue() {
        return -1;  // C(γ) = -1
    }

    static int neutralPionEigenvalue() {
        return +1;  // C(π⁰) = +1
    }

    /**
     * @brief C-parity of neutral mesons
     *
     * π⁰, η: C = +1
     * ω, φ: C = -1
     */
    static int cParity(const std::string& meson) {
        if (meson == "pi0" || meson == "eta") return +1;
        if (meson == "omega" || meson == "phi") return -1;
        return 0;  // Not C eigenstate
    }

    /**
     * @brief C violation in weak interactions
     *
     * Weak: C violated!
     * ν_L couples to W⁻, but ν̄_L does not (only ν̄_R does)
     */
    static bool conservedInWeak() {
        return false;  // Violated!
    }

    static std::string weakCViolation() {
        return "Weak interactions violate C:\n"
               "Only left-handed neutrinos exist (ν_L)\n"
               "Only right-handed antineutrinos exist (ν̄_R)\n"
               "C: ν_L → ν̄_L (does not exist!)";
    }

    /**
     * @brief C in strong and EM
     *
     * Strong: C conserved
     * EM: C conserved
     */
    static bool conservedInStrong() {
        return true;
    }

    static bool conservedInEM() {
        return true;
    }
};

/**
 * @class TimeReversal
 * @brief Time reversal (T): t → -t, motion reversal
 *
 * T|ψ(r, p, t)⟩ = |ψ(r, -p, -t)⟩
 *
 * Reverses: momentum, angular momentum, magnetic field
 */
class TimeReversal {
public:
    /**
     * @brief Time reversal transformation
     *
     * T: t → -t
     * Position: r → r (unchanged)
     * Momentum: p → -p
     * Angular momentum: L → -L
     * Spin: S → -S
     */
    static std::string transformation() {
        return "T: t → -t (motion reversal)\n"
               "r → r (position unchanged)\n"
               "p → -p (momentum reverses)\n"
               "L → -L (angular momentum reverses)\n"
               "S → -S (spin reverses)";
    }

    /**
     * @brief T is antiunitary operator
     *
     * T is antiunitary: T(iψ) = -i T(ψ)
     * Complex conjugates coefficients
     */
    static std::string antiunitarity() {
        return "T is antiunitary operator:\n"
               "T(α|ψ⟩ + β|φ⟩) = α* T|ψ⟩ + β* T|φ⟩\n"
               "Complex conjugates amplitudes";
    }

    /**
     * @brief Magnetic field transformation
     *
     * T: B → -B (magnetic field reverses)
     * Electric field: E → E (unchanged)
     */
    static std::string electricMagneticFields() {
        return "T: E → E (electric field unchanged)\n"
               "T: B → -B (magnetic field reverses)";
    }

    /**
     * @brief T violation in weak interactions
     *
     * T is violated (same as CP violation by CPT theorem)
     * Observed in neutral kaon and B meson systems
     */
    static bool conservedInWeak() {
        return false;  // Violated!
    }

    static std::string tViolation() {
        return "T violation observed in:\n"
               "- Neutral kaon system (K⁰ - K̄⁰)\n"
               "- B meson system (B⁰ - B̄⁰)\n"
               "- Equivalent to CP violation by CPT theorem";
    }

    /**
     * @brief T in strong and EM
     */
    static bool conservedInStrong() {
        return true;
    }

    static bool conservedInEM() {
        return true;
    }
};

/**
 * @class CPTTheorem
 * @brief CPT theorem: combined CPT is always conserved
 *
 * CPT is fundamental symmetry of all quantum field theories
 * Consequences:
 * - Particle and antiparticle have identical mass
 * - Particle and antiparticle have equal lifetime
 * - CP violation ⟺ T violation
 */
class CPTTheorem {
public:
    /**
     * @brief CPT theorem statement
     *
     * All Lorentz-invariant quantum field theories
     * conserve CPT symmetry
     */
    static std::string statement() {
        return "CPT Theorem (Pauli, Lüders, Schwinger, 1951-1957):\n"
               "Any local, Lorentz-invariant QFT with Hermitian Hamiltonian\n"
               "is invariant under combined CPT transformation";
    }

    /**
     * @brief CPT transformation
     *
     * CPT: particle(r, p, t) → antiparticle(-r, -p, -t)
     */
    static std::string transformation() {
        return "CPT: |particle, r, p, t⟩ → |antiparticle, -r, -p, -t⟩\n"
               "Combines parity, charge conjugation, time reversal";
    }

    /**
     * @brief Consequences of CPT invariance
     *
     * 1. m(particle) = m(antiparticle)
     * 2. τ(particle) = τ(antiparticle)
     * 3. |Q(particle)| = |Q(antiparticle)|
     * 4. g(particle) = g(antiparticle) (g-factor)
     */
    static std::string consequences() {
        return "CPT invariance implies:\n"
               "1. Identical masses: m_e = m_ē (to 1 part in 10⁸)\n"
               "2. Equal lifetimes: τ(K⁺) = τ(K⁻)\n"
               "3. Equal charges: |Q_p| = |Q_p̄| (to 1 part in 10²¹)\n"
               "4. Equal magnetic moments: |μ_p| = |μ_p̄|";
    }

    /**
     * @brief Test of CPT: electron-positron mass
     *
     * m_e = m_e⁺ to 1 part in 10⁸
     */
    static double massEqualityTest() {
        return 1e-8;  // |m_e - m_e+|/m_e < 10⁻⁸
    }

    /**
     * @brief Test of CPT: proton-antiproton charge
     *
     * Q_p + Q_p̄ = 0 to 1 part in 10²¹
     */
    static double chargeEqualityTest() {
        return 1e-21;  // |Q_p + Q_p̄|/e < 10⁻²¹
    }

    /**
     * @brief CPT ⟺ CP violation implies T violation
     *
     * If CPT is conserved but CP is violated,
     * then T must also be violated
     */
    static std::string cpViolationImpliesTViolation() {
        return "CPT theorem ⟹ CP violation ⟺ T violation\n"
               "Observed CP violation in K⁰ and B⁰\n"
               "Implies T violation (confirmed experimentally)";
    }

    /**
     * @brief Is CPT violated?
     *
     * No evidence for CPT violation
     * All tests consistent with CPT invariance
     */
    static bool isViolated() {
        return false;  // No evidence of CPT violation
    }

    static std::string experimentalStatus() {
        return "CPT invariance tested to extreme precision:\n"
               "No evidence of violation in any system\n"
               "Most precise test: K⁰ - K̄⁰ mass difference";
    }
};

/**
 * @class CPViolation
 * @brief CP violation: combined CP not conserved
 *
 * Discovered 1964 in neutral kaon system (Cronin & Fitch, Nobel 1980)
 * Essential for matter-antimatter asymmetry (baryogenesis)
 */
class CPViolation {
public:
    /**
     * @brief CP transformation
     *
     * CP: particle(r, t) → antiparticle(-r, t)
     */
    static std::string transformation() {
        return "CP: |particle, r⟩ → |antiparticle, -r⟩\n"
               "Combines parity and charge conjugation";
    }

    /**
     * @brief Discovery: K_L → π⁺π⁻ (1964)
     *
     * Cronin & Fitch observed K_L (CP = -1) decaying to ππ (CP = +1)
     * Proved CP is violated in weak interactions!
     */
    static std::string discovery() {
        return "Cronin & Fitch (1964): K_L → π⁺π⁻\n"
               "K_L has CP = -1, but ππ has CP = +1\n"
               "Direct evidence of CP violation\n"
               "Nobel Prize 1980";
    }

    /**
     * @brief Magnitude of CP violation in kaons
     *
     * ε ≈ 2.2×10⁻³ (small but nonzero!)
     */
    static double kaonCPViolation() {
        return 2.2e-3;  // ε parameter
    }

    /**
     * @brief CP violation in B mesons
     *
     * Large CP asymmetries observed (BaBar, Belle)
     * sin(2β) ≈ 0.7 from B⁰ → J/ψ K_S
     */
    static double bMesonCPViolation() {
        return 0.7;  // sin(2β)
    }

    /**
     * @brief Sources of CP violation in SM
     *
     * 1. CKM matrix complex phase (quark sector)
     * 2. Θ_QCD (strong CP problem, Θ < 10⁻¹⁰)
     * 3. Possible PMNS phase (lepton sector, not yet observed)
     */
    static std::string sourcesInSM() {
        return "CP violation sources in Standard Model:\n"
               "1. CKM phase δ_CKM ~ 68° (quark sector) - OBSERVED\n"
               "2. Strong CP: Θ_QCD < 10⁻¹⁰ (unnaturally small!)\n"
               "3. PMNS phase δ_PMNS (lepton sector) - not yet confirmed";
    }

    /**
     * @brief Importance for cosmology
     *
     * CP violation is necessary for baryogenesis
     * (Sakharov condition)
     *
     * BUT: SM CP violation is too small to explain
     * matter-antimatter asymmetry!
     */
    static std::string cosmologicalImportance() {
        return "CP violation needed for baryogenesis (Sakharov)\n"
               "Explains matter-antimatter asymmetry\n"
               "Problem: SM CP violation too weak by ~10⁹!\n"
               "Suggests new physics beyond SM";
    }

    /**
     * @brief Electric dipole moment (EDM)
     *
     * Permanent EDM violates both P and T (hence CP)
     * Current limits: d_e < 10⁻²⁹ e·cm (extremely small!)
     */
    static double electronEDMLimit() {
        return 1e-29;  // e·cm
    }

    static std::string edmSignificance() {
        return "Electric dipole moment (EDM) violates CP:\n"
               "Neutron: d_n < 10⁻²⁶ e·cm\n"
               "Electron: d_e < 10⁻²⁹ e·cm\n"
               "Constrains new CP violation sources beyond SM";
    }
};

/**
 * @class TViolation
 * @brief T violation (equivalent to CP violation)
 *
 * By CPT theorem: CP violation ⟺ T violation
 * Directly observed in neutral kaon system
 */
class TViolation {
public:
    /**
     * @brief Direct observation of T violation
     *
     * CPLEAR experiment (1998): K⁰ - K̄⁰ oscillations
     * BaBar (2012): B meson decays
     */
    static std::string directObservation() {
        return "Direct T violation observations:\n"
               "- CPLEAR (1998): K⁰ ↔ K̄⁰ oscillations\n"
               "- BaBar (2012): B⁰ → J/ψ K_L asymmetry\n"
               "- Confirms CPT theorem: CP violation = T violation";
    }

    /**
     * @brief T violation magnitude
     *
     * Similar to CP violation (~10⁻³ in kaons)
     */
    static double magnitude() {
        return 2.2e-3;  // Same as CP violation
    }

    /**
     * @brief Implications for thermodynamics
     *
     * Microscopic T violation vs macroscopic arrow of time
     * Thermodynamic arrow from statistical mechanics, not T violation
     */
    static std::string thermodynamicArrow() {
        return "Microscopic T violation ≠ thermodynamic arrow of time\n"
               "Entropy increase is statistical (2nd law)\n"
               "Not due to fundamental T violation";
    }
};

/**
 * @class CombinedSymmetries
 * @brief Summary of discrete symmetry conservation
 */
class CombinedSymmetries {
public:
    /**
     * @brief Symmetry conservation table
     *
     * Interaction  |  C  |  P  |  T  | CP  | CPT
     * ------------------------------------------------
     * Strong       |  ✓  |  ✓  |  ✓  |  ✓  |  ✓
     * EM           |  ✓  |  ✓  |  ✓  |  ✓  |  ✓
     * Weak         |  ✗  |  ✗  |  ✗  |  ✗  |  ✓
     */
    static std::string conservationTable() {
        return
            "Discrete Symmetry Conservation:\n"
            "\n"
            "Interaction  |  C  |  P  |  T  | CP  | CPT\n"
            "-------------------------------------------\n"
            "Strong       |  ✓  |  ✓  |  ✓  |  ✓  |  ✓\n"
            "EM           |  ✓  |  ✓  |  ✓  |  ✓  |  ✓\n"
            "Weak         |  ✗  |  ✗  |  ✗  |  ✗  |  ✓\n"
            "\n"
            "CPT is ALWAYS conserved (fundamental theorem)";
    }

    /**
     * @brief Historical timeline
     */
    static std::string timeline() {
        return
            "1927: Wigner introduces parity\n"
            "1956: Lee & Yang propose parity violation\n"
            "1957: Wu experiment confirms P violation\n"
            "1957: Landau proposes CP conservation\n"
            "1964: Cronin & Fitch discover CP violation (K_L → ππ)\n"
            "1980: Nobel Prize (Cronin & Fitch)\n"
            "1998: CPLEAR observes direct T violation\n"
            "2001: BaBar/Belle observe large CP violation in B mesons\n"
            "2008: Kobayashi & Maskawa Nobel Prize (CKM mechanism)";
    }

    /**
     * @brief Physical significance
     */
    static std::string significance() {
        return
            "Symmetry violations reveal fundamental physics:\n"
            "- P violation → weak force is chiral (V-A theory)\n"
            "- C violation → neutrinos are left-handed only\n"
            "- CP violation → matter-antimatter asymmetry (baryogenesis)\n"
            "- CPT invariance → particle-antiparticle equivalence\n"
            "\n"
            "Weak interaction is unique: violates C, P, T, CP\n"
            "But CPT is sacred: always conserved!";
    }
};

} // namespace physics::advanced::gauge_theory

#endif // PHYSICS_ADVANCED_GAUGE_THEORY_SYMMETRIES_HPP
