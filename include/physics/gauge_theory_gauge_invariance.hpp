#ifndef PHYSICS_ADVANCED_GAUGE_THEORY_GAUGE_INVARIANCE_HPP
#define PHYSICS_ADVANCED_GAUGE_THEORY_GAUGE_INVARIANCE_HPP

#include <string>
#include <cmath>
#include <complex>
#include <vector>

/**
 * @file gauge_invariance.hpp
 * @brief Gauge transformations and gauge invariance
 *
 * Implements:
 * - Local gauge transformations U(1), SU(2), SU(3)
 * - Gauge invariance and minimal coupling
 * - Covariant derivative
 * - Gauge fields (photon, W, Z, gluons)
 * - Electroweak gauge invariance SU(2)_L × U(1)_Y
 * - QCD gauge invariance SU(3)_C
 */

namespace physics::advanced::gauge_theory {

/**
 * @class GaugeTransformation
 * @brief Local gauge transformations
 *
 * Global symmetry: same transformation everywhere
 * Local gauge symmetry: transformation varies with spacetime point x
 */
class GaugeTransformation {
public:
    /**
     * @brief Global vs local symmetry
     *
     * Global: ψ(x) → e^(iα) ψ(x) (α constant)
     * Local gauge: ψ(x) → e^(iα(x)) ψ(x) (α(x) varies)
     */
    static std::string globalVsLocal() {
        return "Global symmetry: α = constant\n"
               "ψ(x) → e^(iα) ψ(x)\n"
               "\n"
               "Local gauge symmetry: α = α(x)\n"
               "ψ(x) → e^(iα(x)) ψ(x)\n"
               "\n"
               "Gauge invariance requires introducing gauge field!";
    }

    /**
     * @brief Why gauge symmetry?
     *
     * Gauge symmetry → interaction structure
     * Forces us to introduce gauge bosons (photon, W, Z, gluons)
     *
     * "Symmetry dictates interaction"
     */
    static std::string motivation() {
        return "Gauge symmetry is not just redundancy:\n"
               "\n"
               "Local gauge invariance → gauge bosons emerge\n"
               "U(1): photon (QED)\n"
               "SU(2): W±, Z⁰ (weak)\n"
               "SU(3): 8 gluons (QCD)\n"
               "\n"
               "\"Symmetry dictates the forces of nature!\"";
    }

    /**
     * @brief Gauge redundancy
     *
     * Different gauge field configurations A_μ and A_μ'
     * describe the same physical state
     *
     * A_μ' = A_μ + ∂_μ χ (gauge transformation)
     */
    static std::string gaugeRedundancy() {
        return "Gauge transformation of A_μ:\n"
               "\n"
               "A_μ → A_μ' = A_μ + (1/e) ∂_μ α(x)\n"
               "\n"
               "A_μ and A_μ' describe same physics\n"
               "Redundancy in description (gauge freedom)\n"
               "Physical: field strength F_μν (gauge invariant)";
    }
};

/**
 * @class U1GaugeTheory
 * @brief U(1) gauge theory (Quantum Electrodynamics)
 *
 * Simplest gauge theory: QED
 * Gauge group: U(1) (phase rotations)
 * Gauge boson: photon A_μ
 */
class U1GaugeTheory {
public:
    /**
     * @brief U(1) gauge transformation
     *
     * ψ(x) → e^(ie α(x)) ψ(x)
     * A_μ(x) → A_μ(x) + ∂_μ α(x)
     *
     * where e is electric charge
     */
    static std::string gaugeTransformation() {
        return "U(1) gauge transformation:\n"
               "\n"
               "ψ(x) → ψ'(x) = e^(ie α(x)) ψ(x)\n"
               "A_μ(x) → A_μ'(x) = A_μ(x) + ∂_μ α(x)\n"
               "\n"
               "e = electric charge\n"
               "α(x) = arbitrary function (gauge parameter)";
    }

    /**
     * @brief Covariant derivative
     *
     * D_μ = ∂_μ - ie A_μ
     *
     * Transforms covariantly:
     * D_μ ψ → e^(ieα) D_μ ψ
     */
    static std::string covariantDerivative() {
        return "Covariant derivative:\n"
               "\n"
               "D_μ = ∂_μ - ie A_μ\n"
               "\n"
               "Transforms covariantly:\n"
               "D_μ ψ → e^(ieα(x)) D_μ ψ\n"
               "\n"
               "Ordinary derivative ∂_μ ψ is NOT gauge covariant!";
    }

    /**
     * @brief Minimal coupling
     *
     * Free Dirac: iγ^μ ∂_μ ψ - m ψ = 0
     * QED: iγ^μ D_μ ψ - m ψ = 0
     *
     * Substitution: ∂_μ → D_μ = ∂_μ - ie A_μ
     *
     * Introduces electromagnetic interaction!
     */
    static std::string minimalCoupling() {
        return "Minimal coupling: ∂_μ → D_μ\n"
               "\n"
               "Free: iγ^μ ∂_μ ψ - m ψ = 0\n"
               "QED: iγ^μ (∂_μ - ie A_μ) ψ - m ψ = 0\n"
               "\n"
               "Interaction term: e ψ̄ γ^μ ψ A_μ\n"
               "\n"
               "Gauge invariance → electromagnetic interaction!";
    }

    /**
     * @brief Field strength tensor
     *
     * F_μν = ∂_μ A_ν - ∂_ν A_μ
     *
     * Gauge invariant: F_μν → F_μν (unchanged)
     */
    static std::string fieldStrength() {
        return "Field strength tensor:\n"
               "\n"
               "F_μν = ∂_μ A_ν - ∂_ν A_μ\n"
               "\n"
               "Gauge invariant: F_μν' = F_μν\n"
               "\n"
               "Physical fields: E = -∇φ - ∂A/∂t\n"
               "                 B = ∇ × A";
    }

    /**
     * @brief QED Lagrangian
     *
     * L = ψ̄ (iγ^μ D_μ - m) ψ - (1/4) F_μν F^μν
     *
     * Gauge invariant!
     */
    static std::string lagrangian() {
        return "QED Lagrangian:\n"
               "\n"
               "L = ψ̄ (iγ^μ D_μ - m) ψ - (1/4) F_μν F^μν\n"
               "\n"
               "= ψ̄ (iγ^μ ∂_μ - m) ψ - (1/4) F² - e ψ̄ γ^μ ψ A_μ\n"
               "\n"
               "Fermion + Photon + Interaction\n"
               "Fully gauge invariant under U(1)!";
    }

    /**
     * @brief Charge quantization
     *
     * Electric charges are quantized: Q = n e
     * Explained by Dirac monopole condition
     * Or by GUT unification
     */
    static std::string chargeQuantization() {
        return "Charge quantization:\n"
               "\n"
               "Observed: all charges are multiples of e\n"
               "Quark: Q = ±e/3, ±2e/3\n"
               "Lepton: Q = 0, ±e\n"
               "\n"
               "Not explained by U(1) alone\n"
               "Requires: monopoles or GUT";
    }
};

/**
 * @class SU2GaugeTheory
 * @brief SU(2) gauge theory (Weak isospin)
 *
 * Non-Abelian gauge theory
 * Gauge group: SU(2) (weak isospin)
 * Gauge bosons: W₁, W₂, W₃ (3 generators)
 */
class SU2GaugeTheory {
public:
    /**
     * @brief SU(2) gauge transformation
     *
     * ψ(x) → U(x) ψ(x) = e^(ig τᵃ αᵃ(x)) ψ(x)
     *
     * where τᵃ = σᵃ/2 (Pauli matrices/2)
     *       a = 1, 2, 3 (3 generators)
     */
    static std::string gaugeTransformation() {
        return "SU(2) gauge transformation:\n"
               "\n"
               "ψ(x) → U(x) ψ(x) = exp(ig τᵃ αᵃ(x)) ψ(x)\n"
               "\n"
               "τᵃ = σᵃ/2 (Pauli matrices)\n"
               "a = 1, 2, 3 (three generators)\n"
               "αᵃ(x) = three gauge parameters\n"
               "\n"
               "Non-Abelian: generators don't commute!";
    }

    /**
     * @brief SU(2) generators and commutation relations
     *
     * [τᵃ, τᵇ] = i ε^(abc) τ^c
     *
     * Structure constants: ε^(abc) (totally antisymmetric)
     */
    static std::string lieAlgebra() {
        return "SU(2) Lie algebra:\n"
               "\n"
               "[τᵃ, τᵇ] = i ε^(abc) τ^c\n"
               "\n"
               "Structure constants: f^(abc) = ε^(abc)\n"
               "\n"
               "Non-Abelian: generators don't commute\n"
               "(Unlike U(1) where [Q, Q] = 0)";
    }

    /**
     * @brief Covariant derivative (non-Abelian)
     *
     * D_μ = ∂_μ - ig W_μ^a τ^a = ∂_μ - ig W_μ
     *
     * where W_μ = W_μ^a τ^a (three gauge fields)
     */
    static std::string covariantDerivative() {
        return "SU(2) covariant derivative:\n"
               "\n"
               "D_μ = ∂_μ - ig W_μ^a τ^a\n"
               "\n"
               "W_μ^a (a = 1,2,3): three gauge bosons\n"
               "g: SU(2) coupling constant\n"
               "\n"
               "Transforms: D_μ → U D_μ U⁻¹";
    }

    /**
     * @brief Field strength tensor (non-Abelian)
     *
     * W_μν^a = ∂_μ W_ν^a - ∂_ν W_μ^a + g ε^(abc) W_μ^b W_ν^c
     *
     * Extra term from non-Abelian structure!
     * Gauge bosons self-interact
     */
    static std::string fieldStrength() {
        return "SU(2) field strength:\n"
               "\n"
               "W_μν^a = ∂_μ W_ν^a - ∂_ν W_μ^a + g ε^(abc) W_μ^b W_ν^c\n"
               "\n"
               "Non-Abelian term: g ε^(abc) W_μ^b W_ν^c\n"
               "\n"
               "→ W bosons self-interact!\n"
               "(Unlike photons in QED)";
    }

    /**
     * @brief Yang-Mills Lagrangian (pure gauge)
     *
     * L = -(1/4) W_μν^a W^(μν)a
     *
     * Contains cubic and quartic W self-interactions!
     */
    static std::string yangMillsLagrangian() {
        return "Yang-Mills Lagrangian:\n"
               "\n"
               "L = -(1/4) Tr(W_μν W^μν)\n"
               "\n"
               "Expands to:\n"
               "L = -(1/2) ∂_μ W_ν^a ∂^μ W^(νa)\n"
               "    + g ε^(abc) (∂_μ W_ν^a) W^(μb) W^(νc)\n"
               "    - (g²/4) ε^(abe) ε^(cde) W_μ^a W_ν^b W^(μc) W^(νd)\n"
               "\n"
               "3-point and 4-point W vertices!";
    }

    /**
     * @brief Number of gauge bosons
     *
     * SU(N): N² - 1 gauge bosons
     * SU(2): 2² - 1 = 3 (W₁, W₂, W₃)
     * SU(3): 3² - 1 = 8 (8 gluons)
     */
    static int numberOfGaugeBosons(int N) {
        return N * N - 1;
    }
};

/**
 * @class ElectroweakGaugeInvariance
 * @brief SU(2)_L × U(1)_Y gauge invariance (Electroweak theory)
 *
 * Glashow-Weinberg-Salam model (Nobel Prize 1979)
 * Unifies electromagnetic and weak interactions
 */
class ElectroweakGaugeInvariance {
public:
    /**
     * @brief Gauge group: SU(2)_L × U(1)_Y
     *
     * SU(2)_L: weak isospin (couples to left-handed fermions only)
     * U(1)_Y: weak hypercharge
     *
     * NOT SU(2) × U(1)_EM !
     */
    static std::string gaugeGroup() {
        return "Electroweak gauge group:\n"
               "\n"
               "G = SU(2)_L × U(1)_Y\n"
               "\n"
               "SU(2)_L: weak isospin (L = left-handed)\n"
               "U(1)_Y: weak hypercharge\n"
               "\n"
               "4 gauge bosons: W₁, W₂, W₃, B\n"
               "After EWSB: W±, Z⁰, γ";
    }

    /**
     * @brief Gauge fields and couplings
     *
     * SU(2)_L: W_μ^a (a = 1,2,3), coupling g
     * U(1)_Y: B_μ, coupling g'
     *
     * Covariant derivative:
     * D_μ = ∂_μ - ig W_μ^a τ^a - ig' Y B_μ
     */
    static std::string gaugeFields() {
        return "Electroweak gauge fields:\n"
               "\n"
               "SU(2)_L: W_μ^(1,2,3), coupling g\n"
               "U(1)_Y: B_μ, coupling g'\n"
               "\n"
               "Covariant derivative:\n"
               "D_μ = ∂_μ - ig W_μ^a τ^a - ig' Y B_μ/2\n"
               "\n"
               "Before EWSB: all massless";
    }

    /**
     * @brief Weak hypercharge Y
     *
     * Y = 2(Q - T₃)
     *
     * where Q = electric charge
     *       T₃ = 3rd component of weak isospin
     *
     * Left doublet (e, ν_e)_L: Y = -1
     * Right singlet e_R: Y = -2
     */
    static std::string hypercharge() {
        return "Weak hypercharge: Y = 2(Q - T₃)\n"
               "\n"
               "Left doublets:\n"
               "(ν_e, e⁻)_L: T₃ = (±1/2), Q = (0, -1), Y = -1\n"
               "(u, d)_L: T₃ = (±1/2), Q = (2/3, -1/3), Y = 1/3\n"
               "\n"
               "Right singlets:\n"
               "e_R: T₃ = 0, Q = -1, Y = -2\n"
               "u_R: T₃ = 0, Q = 2/3, Y = 4/3\n"
               "d_R: T₃ = 0, Q = -1/3, Y = -2/3";
    }

    /**
     * @brief Physical gauge bosons (after EWSB)
     *
     * W± = (W₁ ∓ i W₂)/√2 (charged)
     * Z⁰ = cos(θ_W) W₃ - sin(θ_W) B (neutral, massive)
     * γ = sin(θ_W) W₃ + cos(θ_W) B (photon, massless)
     *
     * where θ_W = Weinberg angle
     */
    static std::string physicalBosons() {
        return "Physical gauge bosons:\n"
               "\n"
               "W± = (W₁ ∓ i W₂)/√2 (m_W = 80.4 GeV)\n"
               "Z⁰ = cos θ_W W₃ - sin θ_W B (m_Z = 91.2 GeV)\n"
               "γ = sin θ_W W₃ + cos θ_W B (m_γ = 0)\n"
               "\n"
               "θ_W: Weinberg (weak mixing) angle\n"
               "sin²θ_W ≈ 0.23";
    }

    /**
     * @brief Weinberg angle
     *
     * tan θ_W = g'/g
     *
     * Relates SU(2) and U(1) couplings
     * sin²θ_W ≈ 0.231 (measured precisely)
     */
    static double weinbergAngle() {
        return std::asin(std::sqrt(0.231));  // radians
    }

    static double sin2ThetaW() {
        return 0.231;  // sin²θ_W (MS-bar scheme, m_Z scale)
    }

    /**
     * @brief Electromagnetic coupling from EW
     *
     * e = g sin θ_W = g' cos θ_W
     *
     * Electric charge emerges from EW symmetry!
     */
    static std::string electromagneticCoupling() {
        return "Electromagnetic coupling:\n"
               "\n"
               "e = g sin θ_W = g' cos θ_W\n"
               "\n"
               "EM is mixture of SU(2) and U(1)!\n"
               "Photon = mixture of W₃ and B";
    }

    /**
     * @brief Custodial SU(2) symmetry
     *
     * Approximate symmetry protecting m_W/m_Z ratio
     * m_W/m_Z = cos θ_W (tree level)
     */
    static std::string custodialSymmetry() {
        return "Custodial SU(2)_V symmetry:\n"
               "\n"
               "Protects: m_W = m_Z cos θ_W\n"
               "ρ = m_W²/(m_Z² cos²θ_W) ≈ 1\n"
               "\n"
               "Measured: ρ = 1.00038 ± 0.00020\n"
               "Confirms Higgs doublet structure!";
    }
};

/**
 * @class SU3GaugeTheory
 * @brief SU(3) gauge theory (Quantum Chromodynamics)
 *
 * Gauge group: SU(3)_C (color)
 * Gauge bosons: 8 gluons
 */
class SU3GaugeTheory {
public:
    /**
     * @brief SU(3)_C gauge group
     *
     * 3 colors: red, green, blue
     * 3² - 1 = 8 gluons
     */
    static std::string gaugeGroup() {
        return "QCD gauge group: SU(3)_C\n"
               "\n"
               "3 colors: r, g, b\n"
               "8 gluons: G^a (a = 1,...,8)\n"
               "\n"
               "Quarks: color triplet (r, g, b)\n"
               "Gluons: color octet (rḡ, rḃ, gḃ, ...)";
    }

    /**
     * @brief Gell-Mann matrices (SU(3) generators)
     *
     * λ^a (a = 1,...,8) are 3×3 traceless Hermitian matrices
     * Generators: T^a = λ^a/2
     */
    static std::string gellMannMatrices() {
        return "SU(3) generators: T^a = λ^a/2\n"
               "\n"
               "Gell-Mann matrices λ^a (a = 1,...,8)\n"
               "3×3 traceless, Hermitian\n"
               "\n"
               "[T^a, T^b] = i f^(abc) T^c\n"
               "f^(abc): SU(3) structure constants";
    }

    /**
     * @brief QCD covariant derivative
     *
     * D_μ = ∂_μ - ig_s G_μ^a T^a
     *
     * where g_s is strong coupling constant
     */
    static std::string covariantDerivative() {
        return "QCD covariant derivative:\n"
               "\n"
               "D_μ = ∂_μ - ig_s G_μ^a T^a\n"
               "\n"
               "G_μ^a (a = 1,...,8): gluon fields\n"
               "g_s: strong coupling (α_s = g_s²/4π)";
    }

    /**
     * @brief Gluon field strength
     *
     * G_μν^a = ∂_μ G_ν^a - ∂_ν G_μ^a + g_s f^(abc) G_μ^b G_ν^c
     *
     * Gluon self-interaction!
     */
    static std::string fieldStrength() {
        return "Gluon field strength:\n"
               "\n"
               "G_μν^a = ∂_μ G_ν^a - ∂_ν G_μ^a + g_s f^(abc) G_μ^b G_ν^c\n"
               "\n"
               "Non-Abelian term → gluon self-coupling\n"
               "3-gluon and 4-gluon vertices!";
    }

    /**
     * @brief QCD Lagrangian
     *
     * L = ψ̄ (iγ^μ D_μ - m) ψ - (1/4) G_μν^a G^(μν)a
     *
     * Quarks + Gluons + Interactions
     */
    static std::string lagrangian() {
        return "QCD Lagrangian:\n"
               "\n"
               "L = Σ_q ψ̄_q (iγ^μ D_μ - m_q) ψ_q - (1/4) G_μν^a G^(μν)a\n"
               "\n"
               "Sum over quark flavors (u, d, s, c, b, t)\n"
               "Pure gauge term with 3-gluon, 4-gluon vertices\n"
               "\n"
               "Asymptotic freedom + confinement!";
    }

    /**
     * @brief Color confinement
     *
     * Only color-singlet states observable
     * Quarks and gluons confined in hadrons
     *
     * Unsolved problem: analytic proof of confinement
     */
    static std::string confinement() {
        return "Color confinement:\n"
               "\n"
               "Only color-singlet (colorless) states exist\n"
               "Quarks: qqq (baryons), q̄q (mesons)\n"
               "Free quarks/gluons never observed!\n"
               "\n"
               "Mechanism: strong coupling at large distances\n"
               "α_s(r → ∞) → ∞";
    }

    /**
     * @brief Asymptotic freedom
     *
     * α_s(Q²) → 0 as Q² → ∞
     *
     * Quarks nearly free at short distances!
     * Nobel Prize 2004 (Gross, Politzer, Wilczek)
     */
    static std::string asymptoticFreedom() {
        return "Asymptotic freedom (Nobel 2004):\n"
               "\n"
               "α_s(Q²) → 0 as Q² → ∞\n"
               "\n"
               "High energy: quarks nearly free\n"
               "Low energy: quarks confined\n"
               "\n"
               "Running: dα_s/d log Q² < 0 (antiscreening!)";
    }
};

/**
 * @class StandardModelGaugeGroup
 * @brief Complete Standard Model gauge symmetry
 *
 * G_SM = SU(3)_C × SU(2)_L × U(1)_Y
 */
class StandardModelGaugeGroup {
public:
    /**
     * @brief Standard Model gauge group
     *
     * SU(3)_C: color (QCD) → 8 gluons
     * SU(2)_L: weak isospin → W₁, W₂, W₃
     * U(1)_Y: weak hypercharge → B
     *
     * Total: 12 gauge bosons before EWSB
     * After EWSB: 8 gluons + W± + Z⁰ + γ
     */
    static std::string completeGaugeGroup() {
        return "Standard Model gauge group:\n"
               "\n"
               "G_SM = SU(3)_C × SU(2)_L × U(1)_Y\n"
               "\n"
               "Gauge bosons (before EWSB): 12\n"
               "- SU(3)_C: 8 gluons\n"
               "- SU(2)_L: W₁, W₂, W₃\n"
               "- U(1)_Y: B\n"
               "\n"
               "After EWSB: 8 g + W± + Z⁰ + γ (12 total)";
    }

    /**
     * @brief Gauge coupling constants
     *
     * g_s: strong (SU(3))
     * g: weak (SU(2))
     * g': hypercharge (U(1))
     *
     * At m_Z: α_s ≈ 0.118, α ≈ 1/128, α_W ≈ 1/30
     */
    static std::string couplingConstants() {
        return "SM gauge couplings (at m_Z):\n"
               "\n"
               "α_s = g_s²/4π ≈ 0.118 (strong)\n"
               "α = e²/4π ≈ 1/128 (EM)\n"
               "α_2 = g²/4π ≈ 1/30 (SU(2))\n"
               "\n"
               "Running couplings unify at M_GUT ~ 10¹⁶ GeV!";
    }

    /**
     * @brief Quantum numbers
     *
     * Particle           SU(3)  SU(2)  U(1)_Y
     * ----------------------------------------
     * (u, d)_L           3      2      1/3
     * u_R                3      1      4/3
     * d_R                3      1      -2/3
     * (ν_e, e)_L         1      2      -1
     * e_R                1      1      -2
     * Higgs doublet      1      2      1
     */
    static std::string quantumNumbers() {
        return
            "Particle quantum numbers:\n"
            "\n"
            "Particle          SU(3)  SU(2)  U(1)_Y\n"
            "---------------------------------------\n"
            "Quark doublet_L     3      2     1/3\n"
            "u_R                 3      1     4/3\n"
            "d_R                 3      1    -2/3\n"
            "Lepton doublet_L    1      2     -1\n"
            "e_R                 1      1     -2\n"
            "Higgs doublet       1      2      1\n"
            "\n"
            "Determines all interactions!";
    }

    /**
     * @brief Anomaly cancellation
     *
     * Gauge anomalies must cancel for consistency
     * SM: miraculous cancellation between quarks and leptons!
     */
    static std::string anomalyCancellation() {
        return "Anomaly cancellation in SM:\n"
               "\n"
               "Tr[T^a {T^b, T^c}] must vanish for all a,b,c\n"
               "\n"
               "Separate cancellation:\n"
               "- SU(3)³: quarks only\n"
               "- SU(2)³: quarks + leptons\n"
               "- U(1)³: quarks + leptons (requires N_colors = 3!)\n"
               "- SU(3)² U(1): quarks\n"
               "- SU(2)² U(1): quarks + leptons\n"
               "\n"
               "All cancel → consistent theory!";
    }

    /**
     * @brief Why 3 colors?
     *
     * Anomaly cancellation requires N_c = 3
     * QCD coupling unification requires N_c = 3
     * Δ⁺⁺ baryon (uuu with J=3/2) requires N_c ≥ 3
     */
    static std::string whyThreeColors() {
        return "Why 3 colors?\n"
               "\n"
               "1. Anomaly cancellation: N_c = 3\n"
               "2. Coupling unification at M_GUT: N_c = 3\n"
               "3. Δ⁺⁺ baryon (uuu, J=3/2): N_c ≥ 3\n"
               "4. π⁰ → γγ decay rate: N_c = 3\n"
               "\n"
               "Multiple independent constraints → N_c = 3!";
    }
};

} // namespace physics::advanced::gauge_theory

#endif // PHYSICS_ADVANCED_GAUGE_THEORY_GAUGE_INVARIANCE_HPP
