#ifndef PHYSICS_ADVANCED_GAUGE_THEORY_CP_VIOLATION_KAONS_HPP
#define PHYSICS_ADVANCED_GAUGE_THEORY_CP_VIOLATION_KAONS_HPP

#include <cmath>
#include <string>
#include <complex>
#include <vector>

/**
 * @file cp_violation_kaons.hpp
 * @brief CP violation in neutral kaon system and CKM matrix
 *
 * Implements:
 * - Neutral kaon mixing (K⁰-K̄⁰ oscillations)
 * - CP violation discovery (Cronin-Fitch 1964)
 * - Direct and indirect CP violation
 * - CKM matrix and quark mixing
 * - Unitarity triangle
 * - B meson CP violation
 */

namespace physics::advanced::gauge_theory {

/**
 * @class NeutralKaonSystem
 * @brief K⁰-K̄⁰ mixing and oscillations
 *
 * K⁰ (d̄s) and K̄⁰ (ds̄) mix via weak interactions
 * Mass eigenstates ≠ flavor eigenstates
 */
class NeutralKaonSystem {
public:
    /**
     * @brief Flavor eigenstates
     *
     * K⁰ = |d̄s⟩ (strangeness S = +1)
     * K̄⁰ = |ds̄⟩ (strangeness S = -1)
     *
     * Produced in strong interactions
     */
    static std::string flavorEigenstates() {
        return "Flavor eigenstates:\n"
               "\n"
               "K⁰ = |d̄s⟩ (S = +1, produced with Λ, π⁺)\n"
               "K̄⁰ = |ds̄⟩ (S = -1, produced with Λ̄, π⁻)\n"
               "\n"
               "Strangeness conserved in strong/EM\n"
               "Strangeness violated in weak decays";
    }

    /**
     * @brief K⁰-K̄⁰ mixing
     *
     * Weak box diagrams with W± allow transitions:
     * K⁰ ↔ K̄⁰
     *
     * Mass eigenstates are mixtures!
     */
    static std::string mixing() {
        return "K⁰-K̄⁰ mixing:\n"
               "\n"
               "Weak interaction: K⁰ ↔ K̄⁰\n"
               "Box diagrams with W± bosons\n"
               "\n"
               "Mass eigenstates ≠ flavor eigenstates\n"
               "K_S, K_L are mixtures of K⁰ and K̄⁰";
    }

    /**
     * @brief CP eigenstates (if CP conserved)
     *
     * K₁ = (K⁰ + K̄⁰)/√2 (CP = +1)
     * K₂ = (K⁰ - K̄⁰)/√2 (CP = -1)
     *
     * If CP exact: K_S = K₁, K_L = K₂
     */
    static std::string cpEigenstates() {
        return "CP eigenstates:\n"
               "\n"
               "K₁ = (K⁰ + K̄⁰)/√2 (CP = +1)\n"
               "K₂ = (K⁰ - K̄⁰)/√2 (CP = -1)\n"
               "\n"
               "If CP conserved:\n"
               "K₁ → ππ (CP = +1)\n"
               "K₂ → πππ (CP = -1)\n"
               "\n"
               "But CP is violated!";
    }

    /**
     * @brief Mass eigenstates (actual)
     *
     * K_S = (K⁰ + K̄⁰)/√2 + O(ε) (short-lived)
     * K_L = (K⁰ - K̄⁰)/√2 + O(ε) (long-lived)
     *
     * Nearly CP eigenstates, but small ε ≠ 0
     */
    static std::string massEigenstates() {
        return "Mass eigenstates:\n"
               "\n"
               "K_S: m ≈ 497.6 MeV, τ ≈ 0.09 ns (short)\n"
               "K_L: m ≈ 497.6 MeV, τ ≈ 52 ns (long)\n"
               "\n"
               "Δm = m_L - m_S ≈ 3.5×10⁻¹² MeV\n"
               "Δm/m ~ 10⁻¹⁴ (tiny!)\n"
               "\n"
               "Nearly degenerate but different lifetimes!";
    }

    /**
     * @brief K_S and K_L masses and lifetimes
     */
    static double kShortMass() {
        return 0.4976;  // GeV
    }

    static double kLongMass() {
        return 0.4976;  // GeV (same to high precision)
    }

    static double kShortLifetime() {
        return 0.89e-10;  // seconds (~0.09 ns)
    }

    static double kLongLifetime() {
        return 5.1e-8;  // seconds (~51 ns)
    }

    static double massDifference() {
        return 3.5e-15;  // GeV (Δm = m_L - m_S)
    }

    /**
     * @brief Oscillation frequency
     *
     * Δm = m_L - m_S
     * Oscillation period: T = 2π/Δm
     *
     * K⁰ ↔ K̄⁰ oscillation period ~ 10⁻¹⁰ s
     */
    static double oscillationPeriod() {
        double hbar = 6.582e-25;  // GeV·s
        double delta_m = massDifference();
        return 2.0 * M_PI * hbar / delta_m;  // seconds
    }

    static std::string oscillations() {
        return "K⁰-K̄⁰ oscillations:\n"
               "\n"
               "|K⁰(t)⟩ = a(t)|K⁰⟩ + b(t)|K̄⁰⟩\n"
               "\n"
               "Oscillation frequency: Δm/ℏ\n"
               "Period: T ~ 10⁻¹⁰ s\n"
               "\n"
               "K⁰ → K̄⁰ → K⁰ → ...";
    }
};

/**
 * @class CPViolationDiscovery
 * @brief Discovery of CP violation (Cronin & Fitch, 1964)
 *
 * Nobel Prize 1980
 */
class CPViolationDiscovery {
public:
    /**
     * @brief K_L → π⁺π⁻ observation
     *
     * If CP conserved: K_L (CP = -1) → πππ only
     * But K_L → ππ (CP = +1) observed!
     *
     * Branching ratio: BR(K_L → π⁺π⁻) ≈ 2×10⁻³
     */
    static std::string cronin FitchExperiment() {
        return "Cronin-Fitch experiment (1964):\n"
               "\n"
               "Expected: K_L (CP = -1) → πππ only\n"
               "Observed: K_L → π⁺π⁻ (CP = +1)\n"
               "\n"
               "BR(K_L → π⁺π⁻) ≈ 2×10⁻³\n"
               "\n"
               "Direct evidence of CP violation!\n"
               "Nobel Prize 1980";
    }

    /**
     * @brief Magnitude of CP violation in kaons
     *
     * ε ≈ 2.2×10⁻³ (small but nonzero!)
     *
     * K_S = (1/√(1+|ε|²)) [K₁ + ε K₂]
     * K_L = (1/√(1+|ε|²)) [K₂ + ε K₁]
     */
    static double epsilonParameter() {
        return 2.228e-3;  // |ε|
    }

    static std::string cpViolationMagnitude() {
        return "CP violation parameter ε:\n"
               "\n"
               "|ε| ≈ 2.2×10⁻³\n"
               "\n"
               "K_L not pure CP = -1 state\n"
               "Small admixture of CP = +1\n"
               "\n"
               "K_L = K₂ + ε K₁ (ε << 1)\n"
               "\n"
               "Small but crucial for matter-antimatter asymmetry!";
    }

    /**
     * @brief Theoretical shock
     *
     * Before 1964: believed C, P, CP all conserved in weak
     * After Wu: P violated, but CP conserved
     * 1964: CP violated too!
     *
     * Only CPT remains sacred
     */
    static std::string theoreticalImpact() {
        return "Impact of CP violation discovery:\n"
               "\n"
               "Before: C, P, CP thought sacred\n"
               "1957: P violated (Wu experiment)\n"
               "1957: \"But CP is conserved\" (Landau)\n"
               "1964: CP violated! (Cronin-Fitch)\n"
               "\n"
               "Only CPT survives as exact symmetry\n"
               "\n"
               "Opened door to matter-antimatter asymmetry";
    }
};

/**
 * @class DirectVsIndirectCPViolation
 * @brief Two types of CP violation in kaons
 */
class DirectVsIndirectCPViolation {
public:
    /**
     * @brief Indirect CP violation (ε)
     *
     * CP violation in mixing (K⁰-K̄⁰ mass matrix)
     * K_L not pure CP eigenstate
     *
     * ε ≈ 2.2×10⁻³
     */
    static std::string indirectCPViolation() {
        return "Indirect CP violation (ε):\n"
               "\n"
               "CP violation in K⁰-K̄⁰ mixing\n"
               "K_L, K_S not pure CP eigenstates\n"
               "\n"
               "Allows: K_L → ππ\n"
               "\n"
               "|ε| ≈ 2.2×10⁻³\n"
               "\n"
               "Discovered 1964 (Cronin-Fitch)";
    }

    /**
     * @brief Direct CP violation (ε')
     *
     * CP violation in decay amplitudes
     * Different decay rates for K⁰ vs K̄⁰
     *
     * ε'/ε ≈ 1.66×10⁻³
     */
    static std::string directCPViolation() {
        return "Direct CP violation (ε'):\n"
               "\n"
               "CP violation in decay amplitude\n"
               "Γ(K⁰ → ππ) ≠ Γ(K̄⁰ → ππ)\n"
               "\n"
               "|ε'/ε| ≈ 1.66×10⁻³\n"
               "\n"
               "Observed 1999 (NA48, KTeV)\n"
               "Confirms CP violation in decay!";
    }

    /**
     * @brief ε'/ε measurement
     *
     * Ratio of direct to indirect CP violation
     * Re(ε'/ε) = (1.66 ± 0.23) × 10⁻³
     *
     * Nonzero → direct CP violation confirmed!
     */
    static double epsilonPrimeOverEpsilon() {
        return 1.66e-3;  // Re(ε'/ε)
    }

    static std::string epsilonPrimeSignificance() {
        return "ε'/ε measurement:\n"
               "\n"
               "Re(ε'/ε) = (1.66 ± 0.23) × 10⁻³\n"
               "\n"
               "Nonzero → CP violation in decay!\n"
               "\n"
               "Required two experiments:\n"
               "- NA48 (CERN)\n"
               "- KTeV (Fermilab)\n"
               "\n"
               "Both confirmed ε' ≠ 0 in 1999";
    }
};

/**
 * @class CKMMatrix
 * @brief Cabibbo-Kobayashi-Maskawa quark mixing matrix
 *
 * 3×3 unitary matrix describing quark flavor mixing
 * Contains CP-violating complex phase
 */
class CKMMatrix {
public:
    /**
     * @brief CKM matrix definition
     *
     * |d'⟩       |V_ud  V_us  V_ub|   |d⟩
     * |s'⟩  =    |V_cd  V_cs  V_cb|   |s⟩
     * |b'⟩       |V_td  V_ts  V_tb|   |b⟩
     *
     * Weak eigenstates (d', s', b') ≠ mass eigenstates (d, s, b)
     */
    static std::string definition() {
        return "CKM matrix:\n"
               "\n"
               "        u       c       t\n"
               "   |V_ud   V_us   V_ub|  d'\n"
               "V =|V_cd   V_cs   V_cb|  s'\n"
               "   |V_td   V_ts   V_tb|  b'\n"
               "\n"
               "Relates weak eigenstates to mass eigenstates\n"
               "3×3 unitary matrix";
    }

    /**
     * @brief Historical development
     *
     * Cabibbo (1963): 2×2 mixing (u, d, s)
     * Kobayashi-Maskawa (1973): 3×3 with CP violation
     * Nobel Prize 2008
     */
    static std::string history() {
        return "CKM matrix history:\n"
               "\n"
               "1963: Cabibbo angle (2×2 matrix)\n"
               "1973: Kobayashi-Maskawa (3×3, CP violation)\n"
               "  Predicted 3rd generation before discovery!\n"
               "1977: b quark discovered\n"
               "1995: t quark discovered\n"
               "2001: BaBar/Belle confirm large CP violation in B\n"
               "2008: Nobel Prize (Kobayashi, Maskawa)";
    }

    /**
     * @brief Unitarity constraint
     *
     * V†V = I
     *
     * 9 unitarity relations:
     * Σ_i |V_ij|² = 1 (rows)
     * Σ_j |V_ij|² = 1 (columns)
     * Σ_k V_ik V*_jk = 0 (orthogonality, i ≠ j)
     */
    static std::string unitarity() {
        return "CKM unitarity:\n"
               "\n"
               "V†V = VV† = I\n"
               "\n"
               "Row normalization:\n"
               "|V_ud|² + |V_us|² + |V_ub|² = 1\n"
               "\n"
               "Orthogonality (unitarity triangle):\n"
               "V_ud V*_ub + V_cd V*_cb + V_td V*_tb = 0";
    }

    /**
     * @brief Wolfenstein parametrization
     *
     * Approximate expansion in λ = sin θ_C ≈ 0.22:
     *
     * V ≈ |  1-λ²/2      λ         Aλ³(ρ-iη)  |
     *     | -λ          1-λ²/2     Aλ²        |
     *     |  Aλ³(1-ρ-iη) -Aλ²       1         |
     *
     * 4 parameters: λ, A, ρ, η (η ≠ 0 → CP violation!)
     */
    static std::string wolfensteinParametrization() {
        return "Wolfenstein parametrization:\n"
               "\n"
               "λ = sin θ_C ≈ 0.22 (Cabibbo angle)\n"
               "A ≈ 0.81\n"
               "ρ ≈ 0.14\n"
               "η ≈ 0.35 (CP-violating phase!)\n"
               "\n"
               "Expansion in powers of λ:\n"
               "V_us ~ λ\n"
               "V_cb ~ Aλ²\n"
               "V_ub ~ Aλ³(ρ - iη)";
    }

    /**
     * @brief CKM matrix elements (magnitudes)
     *
     * From global fits (2022)
     */
    static std::vector<std::vector<double>> matrixElements() {
        return {
            {0.97446, 0.22452, 0.00365},  // |V_ud|, |V_us|, |V_ub|
            {0.22438, 0.97359, 0.04214},  // |V_cd|, |V_cs|, |V_cb|
            {0.00896, 0.04133, 0.999105}  // |V_td|, |V_ts|, |V_tb|
        };
    }

    /**
     * @brief Hierarchy
     *
     * Diagonal elements ~ 1
     * Off-diagonal suppressed by powers of λ
     *
     * |V_ub|, |V_td| ~ λ³ (very small!)
     */
    static std::string hierarchy() {
        return "CKM hierarchy:\n"
               "\n"
               "Diagonal: |V_ud|, |V_cs|, |V_tb| ~ 1\n"
               "\n"
               "Off-diagonal:\n"
               "|V_us|, |V_cd| ~ λ ~ 0.22\n"
               "|V_cb|, |V_ts| ~ λ² ~ 0.04\n"
               "|V_ub|, |V_td| ~ λ³ ~ 0.004\n"
               "\n"
               "Strong hierarchy!";
    }

    /**
     * @brief CP-violating phase
     *
     * 3×3 unitary matrix has 1 physical phase
     * (2 flavor: no phase possible)
     *
     * KM mechanism: minimum 3 generations for CP violation!
     */
    static std::string cpPhase() {
        return "CP-violating phase:\n"
               "\n"
               "3×3 unitary: 9 real parameters\n"
               "- 3 angles\n"
               "- 1 CP-violating phase δ_CKM\n"
               "- 5 unphysical phases (absorbed)\n"
               "\n"
               "δ_CKM ~ 68° (large phase!)\n"
               "\n"
               "2×2 matrix: NO CP violation possible\n"
               "→ Need ≥3 generations for CP violation!";
    }
};

/**
 * @class UnitarityTriangle
 * @brief Geometric representation of CKM unitarity
 *
 * V_ud V*_ub + V_cd V*_cb + V_td V*_tb = 0
 *
 * Forms triangle in complex plane
 */
class UnitarityTriangle {
public:
    /**
     * @brief Triangle construction
     *
     * Divide by V_cd V*_cb to normalize:
     * (V_ud V*_ub)/(V_cd V*_cb) + 1 + (V_td V*_tb)/(V_cd V*_cb) = 0
     *
     * Triangle with vertices at (0,0), (1,0), (ρ̄, η̄)
     */
    static std::string construction() {
        return "Unitarity triangle:\n"
               "\n"
               "V_ud V*_ub + V_cd V*_cb + V_td V*_tb = 0\n"
               "\n"
               "Rescaled (divide by V_cd V*_cb):\n"
               "Vertices: (0,0), (1,0), (ρ̄, η̄)\n"
               "\n"
               "where (ρ̄, η̄) = -(V_td V*_tb)/(V_cd V*_cb)\n"
               "\n"
               "Forms triangle in complex plane";
    }

    /**
     * @brief Apex coordinates
     *
     * (ρ̄, η̄) with ρ̄ ≈ 0.16, η̄ ≈ 0.35
     *
     * Non-zero area → CP violation!
     */
    static double rhoBar() {
        return 0.156;
    }

    static double etaBar() {
        return 0.347;
    }

    /**
     * @brief Triangle angles
     *
     * α = arg[-(V_td V*_tb)/(V_ud V*_ub)]
     * β = arg[-(V_cd V*_cb)/(V_td V*_tb)]
     * γ = arg[-(V_ud V*_ub)/(V_cd V*_cb)]
     *
     * α + β + γ = 180° (triangle!)
     */
    static double alpha() {
        return 88.0;  // degrees
    }

    static double beta() {
        return 22.2;  // degrees
    }

    static double gamma() {
        return 69.8;  // degrees
    }

    static std::string angles() {
        return "Unitarity triangle angles:\n"
               "\n"
               "α ≈ 88°\n"
               "β ≈ 22°\n"
               "γ ≈ 70°\n"
               "\n"
               "α + β + γ = 180° ✓\n"
               "\n"
               "Non-zero area → CP violation!";
    }

    /**
     * @brief Experimental determinations
     *
     * Different processes measure different sides/angles:
     * - β from B⁰ → J/ψ K_S (sin 2β)
     * - α from B⁰ → ππ, ρρ
     * - γ from B → DK
     * - |V_ub|/|V_cb| from semileptonic B decays
     *
     * Overconstrained system tests CKM!
     */
    static std::string experimentalTests() {
        return "Unitarity triangle measurements:\n"
               "\n"
               "sin 2β from B⁰ → J/ψ K_S: β ≈ 22°\n"
               "α from B → ππ: α ≈ 88°\n"
               "γ from B → DK: γ ≈ 70°\n"
               "|V_ub|/|V_cb| from semileptonic B: side lengths\n"
               "\n"
               "Overconstrained → test CKM consistency\n"
               "Result: Excellent agreement! (within ~10%)";
    }

    /**
     * @brief Jarlskog invariant
     *
     * J = Im[V_us V_cb V*_ub V*_cs] = A²λ⁶ ρ η
     *
     * Measure of CP violation (area of triangle)
     * J ≈ 3×10⁻⁵
     */
    static double jarlskogInvariant() {
        return 3.0e-5;
    }

    static std::string jarlskogSignificance() {
        return "Jarlskog invariant J:\n"
               "\n"
               "J = Im[V_us V_cb V*_ub V*_cs]\n"
               "  ≈ 3×10⁻⁵\n"
               "\n"
               "Rephasing-invariant measure of CP violation\n"
               "J = 2 × (triangle area)\n"
               "\n"
               "J = 0 → no CP violation\n"
               "J ≠ 0 → CP is violated!";
    }
};

/**
 * @class BMesonCPViolation
 * @brief CP violation in B meson system
 *
 * Larger than in kaons!
 * BaBar and Belle experiments (2001)
 */
class BMesonCPViolation {
public:
    /**
     * @brief B⁰-B̄⁰ mixing
     *
     * B⁰ = |d̄b⟩, B̄⁰ = |db̄⟩
     *
     * Mix like kaons but with top quark in loop
     * Δm_B/m_B ~ 10⁻¹³ (even smaller than kaons)
     */
    static std::string bMixing() {
        return "B⁰-B̄⁰ mixing:\n"
               "\n"
               "B⁰ = |d̄b⟩ ↔ B̄⁰ = |db̄⟩\n"
               "\n"
               "Box diagram with top quark\n"
               "Δm_B ≈ 0.5 ps⁻¹\n"
               "Oscillation period: T ~ 10 ps\n"
               "\n"
               "Much faster than kaons!";
    }

    /**
     * @brief Gold-plated mode: B⁰ → J/ψ K_S
     *
     * Measures sin 2β directly
     * Clean experimental signature
     *
     * sin 2β = 0.699 ± 0.017
     */
    static double sin2Beta() {
        return 0.699;
    }

    static std::string goldPlatedMode() {
        return "B⁰ → J/ψ K_S (gold-plated):\n"
               "\n"
               "Time-dependent CP asymmetry:\n"
               "S = sin 2β = 0.699 ± 0.017\n"
               "\n"
               "β ≈ 22° (unitarity triangle angle)\n"
               "\n"
               "Discovery of large CP violation in B!\n"
               "BaBar & Belle (2001)";
    }

    /**
     * @brief CP asymmetry
     *
     * A_CP(t) = [Γ(B̄⁰(t) → f) - Γ(B⁰(t) → f)] / [...]
     *         = S sin(Δm_B t) - C cos(Δm_B t)
     *
     * S, C: CP violation parameters
     */
    static std::string cpAsymmetry() {
        return "Time-dependent CP asymmetry:\n"
               "\n"
               "A_CP(t) = S sin(Δm t) - C cos(Δm t)\n"
               "\n"
               "B⁰ → J/ψ K_S:\n"
               "S = sin 2β ≈ 0.70\n"
               "C ≈ 0 (no direct CP violation)\n"
               "\n"
               "Large CP asymmetry!\n"
               "Oscillates with time";
    }

    /**
     * @brief Other B decay modes
     *
     * B → ππ, ρρ, ρπ → measure α
     * B → DK → measure γ
     * B_s → J/ψ φ → measure β_s
     */
    static std::string otherModes() {
        return "Other B CP violation modes:\n"
               "\n"
               "B → π⁺π⁻: measures α\n"
               "B → DK: measures γ\n"
               "B_s → J/ψ φ: measures β_s\n"
               "\n"
               "Complete determination of unitarity triangle!";
    }

    /**
     * @brief LHCb contributions
     *
     * LHCb at CERN: precision B physics
     * - γ angle measurement
     * - B_s mixing and CP violation
     * - Rare B decays
     */
    static std::string lhcbContributions() {
        return "LHCb experiments:\n"
               "\n"
               "Precision measurements:\n"
               "- γ = 69.8° ± 3.5°\n"
               "- CP violation in B_s system\n"
               "- Rare decays B_s → μ⁺μ⁻\n"
               "\n"
               "Confirms and refines unitarity triangle!";
    }
};

/**
 * @class CPViolationInSM
 * @brief CP violation in the Standard Model
 */
class CPViolationInSM {
public:
    /**
     * @brief Sources of CP violation in SM
     *
     * 1. CKM phase δ_CKM ~ 68° (quark sector) ✓ observed
     * 2. QCD θ-term (strong CP) < 10⁻¹⁰ (unnaturally small!)
     * 3. PMNS phase (lepton sector) - not yet observed
     */
    static std::string sourcesInSM() {
        return "CP violation in Standard Model:\n"
               "\n"
               "1. CKM matrix phase δ_CKM ~ 68°\n"
               "   Quark sector - OBSERVED\n"
               "\n"
               "2. Strong CP: θ_QCD < 10⁻¹⁰\n"
               "   Unnatural smallness (strong CP problem)\n"
               "\n"
               "3. PMNS matrix (neutrinos)\n"
               "   Leptonic CP - not yet confirmed";
    }

    /**
     * @brief Insufficiency for baryogenesis
     *
     * SM CP violation too weak by ~10⁹
     * Cannot explain matter-antimatter asymmetry!
     *
     * Requires new physics beyond SM
     */
    static std::string baryogenesisInsufficiency() {
        return "SM CP violation vs baryogenesis:\n"
               "\n"
               "Observed: η = n_B/n_γ ~ 6×10⁻¹⁰\n"
               "SM prediction: η_SM ~ 10⁻¹⁸ (too small!)\n"
               "\n"
               "SM CP violation insufficient by ~10⁹\n"
               "\n"
               "→ Requires new CP violation beyond SM\n"
               "(SUSY? Leptogenesis? New interactions?)";
    }

    /**
     * @brief Strong CP problem
     *
     * QCD allows θ-term: L ⊃ θ (g²/32π²) G·G̃
     *
     * Would violate CP, predict neutron EDM
     * Measured: θ < 10⁻¹⁰ (why so small?)
     *
     * Peccei-Quinn solution: axion
     */
    static std::string strongCPProblem() {
        return "Strong CP problem:\n"
               "\n"
               "QCD θ-term: L ⊃ θ G·G̃ violates CP\n"
               "\n"
               "Neutron EDM: d_n ~ 10⁻¹⁶ θ e·cm\n"
               "Measured: d_n < 10⁻²⁶ e·cm\n"
               "→ θ < 10⁻¹⁰\n"
               "\n"
               "Why is θ so unnaturally small?\n"
               "Peccei-Quinn solution: θ → 0 dynamically (axion)";
    }

    /**
     * @brief Future prospects
     *
     * - More precise unitarity triangle
     * - Search for new CP violation sources
     * - Electric dipole moments (EDM)
     * - Neutrino CP violation
     */
    static std::string futureProspects() {
        return "Future CP violation tests:\n"
               "\n"
               "- Unitarity triangle to 1% precision\n"
               "- Electric dipole moments (d_e, d_n)\n"
               "  Constrain new CP sources\n"
               "- Neutrino oscillations (δ_PMNS)\n"
               "- Rare decays (new physics?)\n"
               "\n"
               "Any deviation → new physics discovery!";
    }
};

} // namespace physics::advanced::gauge_theory

#endif // PHYSICS_ADVANCED_GAUGE_THEORY_CP_VIOLATION_KAONS_HPP
