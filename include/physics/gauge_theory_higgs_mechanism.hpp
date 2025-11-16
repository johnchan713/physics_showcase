#ifndef PHYSICS_ADVANCED_GAUGE_THEORY_HIGGS_MECHANISM_HPP
#define PHYSICS_ADVANCED_GAUGE_THEORY_HIGGS_MECHANISM_HPP

#include <cmath>
#include <string>
#include <complex>
#include <map>

/**
 * @file higgs_mechanism.hpp
 * @brief Higgs mechanism and spontaneous symmetry breaking
 *
 * Implements:
 * - Spontaneous symmetry breaking (SSB)
 * - Goldstone's theorem
 * - Higgs mechanism (gauge bosons acquire mass)
 * - Electroweak symmetry breaking
 * - Yukawa couplings and fermion masses
 * - Higgs boson discovery and properties
 */

namespace physics::advanced::gauge_theory {

/**
 * @class SpontaneousSymmetryBreaking
 * @brief Spontaneous breaking of global or local symmetry
 *
 * Lagrangian symmetric, but vacuum (ground state) breaks symmetry
 *
 * Example: φ⁴ theory with Mexican hat potential
 */
class SpontaneousSymmetryBreaking {
public:
    /**
     * @brief SSB definition
     *
     * Lagrangian: symmetric under transformation
     * Vacuum: not symmetric (chooses direction)
     *
     * "Symmetry is hidden, not lost"
     */
    static std::string definition() {
        return "Spontaneous Symmetry Breaking (SSB):\n"
               "\n"
               "Lagrangian L is symmetric: δL = 0\n"
               "But vacuum |0⟩ is NOT: G|0⟩ ≠ |0⟩\n"
               "\n"
               "Symmetry spontaneously broken\n"
               "\"Symmetry is hidden in the ground state\"";
    }

    /**
     * @brief Mexican hat potential
     *
     * V(φ) = -μ²|φ|² + λ|φ|⁴
     *
     * μ² > 0: unstable at φ = 0
     * Minimum at |φ| = v = μ/√(2λ) (circle of minima)
     */
    static std::string mexicanHatPotential() {
        return "Mexican hat potential:\n"
               "\n"
               "V(φ) = -μ²|φ|² + λ|φ|⁴\n"
               "\n"
               "μ² > 0 → unstable at origin\n"
               "Minimum: |φ| = v = μ/√(2λ)\n"
               "\n"
               "Circle of degenerate vacua!\n"
               "Choosing one breaks U(1) symmetry";
    }

    /**
     * @brief Vacuum expectation value (VEV)
     *
     * ⟨0|φ|0⟩ = v ≠ 0
     *
     * Non-zero VEV breaks symmetry
     */
    static std::string vacuumExpectationValue() {
        return "Vacuum expectation value (VEV):\n"
               "\n"
               "⟨0|φ|0⟩ = v ≠ 0\n"
               "\n"
               "Scalar field has non-zero VEV\n"
               "Spontaneously breaks symmetry\n"
               "\n"
               "Electroweak: v ≈ 246 GeV";
    }

    /**
     * @brief Example: ferromagnetism
     *
     * Hamiltonian: rotationally symmetric
     * Ground state: magnetization points in ONE direction
     *
     * Breaks rotational symmetry spontaneously
     */
    static std::string ferromagnetismAnalogy() {
        return "Analogy: Ferromagnetism\n"
               "\n"
               "Hamiltonian: rotationally symmetric\n"
               "Ground state: M points in definite direction\n"
               "\n"
               "Spontaneous breaking of rotation symmetry\n"
               "Goldstone mode: spin waves (massless)";
    }

    /**
     * @brief Discrete vs continuous SSB
     *
     * Discrete (Z₂): two vacua (φ = ±v)
     * Continuous (U(1)): circle of vacua
     *
     * Continuous SSB → Goldstone bosons
     */
    static std::string discreteVsContinuous() {
        return "Types of SSB:\n"
               "\n"
               "Discrete (Z₂): φ → -φ\n"
               "Vacua: φ = ±v (two choices)\n"
               "No Goldstone bosons\n"
               "\n"
               "Continuous (U(1)): φ → e^(iα) φ\n"
               "Vacua: |φ| = v, any phase\n"
               "Goldstone boson emerges!";
    }
};

/**
 * @class GoldstonesTheorem
 * @brief Goldstone bosons from spontaneously broken continuous symmetries
 *
 * For each broken generator: one massless Goldstone boson
 *
 * Goldstone (1961), Nambu (Nobel 2008)
 */
class GoldstonesTheorem {
public:
    /**
     * @brief Goldstone's theorem statement
     *
     * When continuous global symmetry is spontaneously broken,
     * massless scalar bosons (Goldstone bosons) appear
     *
     * Number of Goldstone bosons = number of broken generators
     */
    static std::string statement() {
        return "Goldstone's Theorem:\n"
               "\n"
               "If continuous global symmetry G spontaneously broken to H,\n"
               "then dim(G) - dim(H) massless Goldstone bosons appear\n"
               "\n"
               "# Goldstone bosons = # broken generators\n"
               "\n"
               "Example: U(1) → nothing: 1 Goldstone boson";
    }

    /**
     * @brief Origin: flat direction in potential
     *
     * Moving along circle of vacua costs NO energy
     * → massless excitation (Goldstone mode)
     */
    static std::string physicalOrigin() {
        return "Physical origin of Goldstone bosons:\n"
               "\n"
               "Potential flat along vacuum manifold\n"
               "Moving between degenerate vacua: zero energy cost\n"
               "\n"
               "Quantum excitation along flat direction\n"
               "→ massless Goldstone boson";
    }

    /**
     * @brief Example: U(1) SSB
     *
     * Complex scalar φ = (φ₁ + iφ₂)/√2
     * Potential: V = -μ²|φ|² + λ|φ|⁴
     *
     * SSB: ⟨φ⟩ = v/√2 (choosing real direction)
     * Goldstone: phase fluctuations (massless)
     * Higgs: radial fluctuations (massive)
     */
    static std::string u1Example() {
        return "U(1) SSB example:\n"
               "\n"
               "φ = (v + h + iχ)/√2 around minimum\n"
               "\n"
               "h: radial (massive, m² = 2μ²) - Higgs mode\n"
               "χ: angular (massless) - Goldstone boson\n"
               "\n"
               "U(1) broken → 1 Goldstone boson";
    }

    /**
     * @brief Example: SU(2) SSB
     *
     * SU(2) → U(1): 3 - 1 = 2 Goldstone bosons
     * SU(2) → nothing: 3 Goldstone bosons
     */
    static std::string su2Example() {
        return "SU(2) SSB examples:\n"
               "\n"
               "SU(2) → U(1): 2 Goldstone bosons\n"
               "Example: chiral symmetry breaking in QCD (π±, π⁰)\n"
               "\n"
               "SU(2) → nothing: 3 Goldstone bosons\n"
               "Example: Electroweak (eaten by W±, Z⁰)";
    }

    /**
     * @brief Pions as pseudo-Goldstone bosons
     *
     * Chiral SU(2)_L × SU(2)_R → SU(2)_V (vector)
     * 3 broken generators → π⁺, π⁻, π⁰
     *
     * Nearly massless (m_π ~ 140 MeV << Λ_QCD)
     * Small mass from explicit symmetry breaking (quark masses)
     */
    static std::string pionsAsGoldstoneBosons() {
        return "Pions as Goldstone bosons:\n"
               "\n"
               "Chiral symmetry: SU(2)_L × SU(2)_R → SU(2)_V\n"
               "3 broken generators → π⁺, π⁻, π⁰\n"
               "\n"
               "m_π ~ 140 MeV (small compared to m_ρ ~ 770 MeV)\n"
               "Mass from quark masses (explicit breaking)\n"
               "\n"
               "Pseudo-Goldstone bosons";
    }
};

/**
 * @class HiggsMechanism
 * @brief Higgs mechanism: Goldstone bosons eaten by gauge bosons
 *
 * Local gauge symmetry: Goldstone bosons become longitudinal
 * components of massive gauge bosons
 *
 * Higgs (1964), Englert-Brout (Nobel 2013)
 */
class HiggsMechanism {
public:
    /**
     * @brief Higgs mechanism summary
     *
     * Spontaneous breaking of LOCAL gauge symmetry:
     * - Goldstone bosons are "eaten" by gauge bosons
     * - Gauge bosons acquire mass (longitudinal DOF)
     * - Massive gauge bosons + massive Higgs
     *
     * No massless Goldstone bosons remain!
     */
    static std::string mechanism() {
        return "Higgs Mechanism:\n"
               "\n"
               "Local gauge symmetry spontaneously broken:\n"
               "1. Goldstone bosons \"eaten\" by gauge bosons\n"
               "2. Gauge bosons acquire mass (3 DOF each)\n"
               "3. Physical Higgs boson remains\n"
               "\n"
               "Massless gauge bosons (2 DOF)\n"
               "+ Goldstone bosons (N DOF)\n"
               "→ Massive gauge bosons (3 DOF)\n"
               "+ Higgs boson (1 DOF)";
    }

    /**
     * @brief Degrees of freedom counting
     *
     * Before SSB: massless gauge boson (2 DOF) + complex scalar (2 DOF) = 4
     * After SSB: massive gauge boson (3 DOF) + Higgs (1 DOF) = 4
     *
     * DOF conserved!
     */
    static std::string degreesFreedom() {
        return "Degrees of freedom:\n"
               "\n"
               "Before SSB:\n"
               "Massless gauge boson A_μ: 2 DOF (transverse)\n"
               "Complex scalar φ: 2 DOF\n"
               "Total: 4 DOF\n"
               "\n"
               "After SSB (Higgs mechanism):\n"
               "Massive gauge boson: 3 DOF (2 trans + 1 long)\n"
               "Physical Higgs h: 1 DOF\n"
               "Total: 4 DOF ✓\n"
               "\n"
               "Goldstone boson eaten → longitudinal mode!";
    }

    /**
     * @brief Abelian Higgs model (U(1))
     *
     * Lagrangian:
     * L = -(1/4)F_μν F^μν + |D_μ φ|² - V(φ)
     * V(φ) = -μ²|φ|² + λ|φ|⁴
     *
     * After SSB: photon acquires mass!
     * (Not our physical photon, just pedagogical example)
     */
    static std::string abelianHiggsModel() {
        return "Abelian Higgs model (U(1)):\n"
               "\n"
               "L = -(1/4)F² + |D_μφ|² - V(φ)\n"
               "D_μ = ∂_μ - ieA_μ\n"
               "\n"
               "SSB: ⟨φ⟩ = v/√2\n"
               "\n"
               "Result:\n"
               "m_A = ev (gauge boson mass)\n"
               "m_h = √(2μ²) = √(2λ)v (Higgs mass)\n"
               "\n"
               "Goldstone eaten by photon!";
    }

    /**
     * @brief Unitary gauge
     *
     * Gauge choice that removes Goldstone bosons explicitly
     * Shows massive gauge bosons and physical Higgs
     *
     * φ(x) = (v + h(x))/√2 (real field only)
     */
    static std::string unitaryGauge() {
        return "Unitary gauge:\n"
               "\n"
               "Use gauge freedom to eliminate Goldstone fields\n"
               "φ(x) = (v + h(x))/√2\n"
               "\n"
               "Advantages:\n"
               "- Physical spectrum clear\n"
               "- No Goldstone ghosts\n"
               "\n"
               "Disadvantages:\n"
               "- Gauge boson propagator complicated\n"
               "- High-energy behavior obscured";
    }

    /**
     * @brief Gauge boson mass from Higgs VEV
     *
     * Kinetic term: |D_μ φ|²
     * With ⟨φ⟩ = v: generates (gauge boson)² mass term
     *
     * m² ~ g² v²
     */
    static std::string massGeneration() {
        return "Gauge boson mass generation:\n"
               "\n"
               "Kinetic: |D_μφ|² = |(∂_μ - igA_μ)φ|²\n"
               "\n"
               "With ⟨φ⟩ = v:\n"
               "|D_μφ|² ⊃ (gv)² A_μ A^μ / 2\n"
               "\n"
               "Gauge boson mass:\n"
               "m_A = gv\n"
               "\n"
               "Mass proportional to coupling and VEV!";
    }

    /**
     * @brief Why gauge bosons can be massive
     *
     * Massless: required by gauge invariance
     * Higgs mechanism: gauge invariance spontaneously broken
     *
     * Goldstone theorem evaded (local, not global symmetry)
     */
    static std::string whyMassiveGaugeBosons() {
        return "Why can gauge bosons be massive?\n"
               "\n"
               "Naively: m² A_μ A^μ term breaks gauge invariance\n"
               "\n"
               "Higgs mechanism:\n"
               "- Gauge invariance preserved in Lagrangian\n"
               "- Spontaneously broken in vacuum\n"
               "- Mass term emerges from ⟨φ⟩ ≠ 0\n"
               "\n"
               "Goldstone theorem: local symmetry → eaten Goldsone\n"
               "No massless scalars remain!";
    }
};

/**
 * @class ElectroweakSymmetryBreaking
 * @brief EWSB: SU(2)_L × U(1)_Y → U(1)_EM
 *
 * Higgs doublet acquires VEV
 * W±, Z⁰ become massive, photon remains massless
 */
class ElectroweakSymmetryBreaking {
public:
    /**
     * @brief EWSB pattern
     *
     * SU(2)_L × U(1)_Y → U(1)_EM
     *
     * 4 generators - 1 unbroken = 3 Goldstone bosons
     * Eaten by W± and Z⁰
     */
    static std::string breakingPattern() {
        return "Electroweak symmetry breaking:\n"
               "\n"
               "SU(2)_L × U(1)_Y → U(1)_EM\n"
               "\n"
               "4 gauge bosons: W₁, W₂, W₃, B\n"
               "→ W±, Z⁰ (massive), γ (massless)\n"
               "\n"
               "3 Goldstone bosons eaten by W±, Z⁰\n"
               "1 Higgs boson remains (h)";
    }

    /**
     * @brief Higgs doublet
     *
     * φ = (φ⁺)
     *     (φ⁰)
     *
     * SU(2) doublet, Y = +1
     * 4 real DOF
     */
    static std::string higgsDoublet() {
        return "Higgs doublet:\n"
               "\n"
               "φ = (φ⁺)  with T = 1/2, Y = +1\n"
               "    (φ⁰)\n"
               "\n"
               "4 real DOF: Re(φ⁺), Im(φ⁺), Re(φ⁰), Im(φ⁰)\n"
               "\n"
               "After EWSB:\n"
               "3 eaten (W±, Z⁰ longitudinal modes)\n"
               "1 physical Higgs h";
    }

    /**
     * @brief Higgs potential
     *
     * V(φ) = -μ²(φ†φ) + λ(φ†φ)²
     *
     * μ² > 0 → SSB
     * Minimum: φ†φ = v²/2 where v = μ/√λ
     */
    static std::string higgsPotential() {
        return "Higgs potential:\n"
               "\n"
               "V(φ) = -μ²(φ†φ) + λ(φ†φ)²\n"
               "\n"
               "Minimum: φ†φ = v²/2\n"
               "v = μ/√λ ≈ 246 GeV\n"
               "\n"
               "Higgs VEV: v ≈ 246 GeV\n"
               "Sets electroweak scale!";
    }

    /**
     * @brief Vacuum choice (unitary gauge)
     *
     * ⟨φ⟩ = (  0  )
     *        ( v/√2 )
     *
     * Breaks SU(2)_L × U(1)_Y but preserves U(1)_EM
     */
    static std::string vacuumChoice() {
        return "Higgs VEV (unitary gauge):\n"
               "\n"
               "⟨φ⟩ = (  0  )\n"
               "      ( v/√2 )\n"
               "\n"
               "Breaks: SU(2)_L (T₃ = -1/2 → 0)\n"
               "Breaks: U(1)_Y\n"
               "Preserves: U(1)_EM (Q = T₃ + Y/2)\n"
               "\n"
               "Photon remains massless!";
    }

    /**
     * @brief Gauge boson masses
     *
     * m_W = gv/2 ≈ 80.4 GeV
     * m_Z = √(g² + g'²) v/2 ≈ 91.2 GeV
     * m_γ = 0 (massless)
     *
     * where v ≈ 246 GeV
     */
    static double wBosonMass() {
        return 80.4;  // GeV
    }

    static double zBosonMass() {
        return 91.2;  // GeV
    }

    static double higgsVEV() {
        return 246.0;  // GeV
    }

    static std::string gaugeBos onMassFormulas() {
        return "Gauge boson masses from EWSB:\n"
               "\n"
               "m_W = gv/2 ≈ 80.4 GeV\n"
               "m_Z = √(g² + g'²) v/2 ≈ 91.2 GeV\n"
               "m_γ = 0\n"
               "\n"
               "Mass ratio:\n"
               "m_W/m_Z = cos θ_W ≈ 0.88\n"
               "\n"
               "From v ≈ 246 GeV (Higgs VEV)";
    }

    /**
     * @brief ρ parameter
     *
     * ρ = m_W²/(m_Z² cos²θ_W)
     *
     * Higgs doublet: ρ = 1 (tree level)
     * Measured: ρ = 1.00038 ± 0.00020
     *
     * Confirms Higgs is doublet!
     */
    static double rhoParameter() {
        return 1.00038;
    }

    static std::string rhoParameterSignificance() {
        return "ρ parameter test:\n"
               "\n"
               "ρ = m_W²/(m_Z² cos²θ_W)\n"
               "\n"
               "Higgs doublet: ρ = 1 (tree level)\n"
               "Measured: ρ = 1.00038 ± 0.00020\n"
               "\n"
               "Confirms: Higgs is SU(2) doublet\n"
               "Rules out triplet, singlet, etc.";
    }
};

/**
 * @class FermionMasses
 * @brief Yukawa couplings give fermions mass
 *
 * Fermion mass terms not allowed by gauge symmetry
 * Yukawa coupling to Higgs → mass after EWSB
 */
class FermionMasses {
public:
    /**
     * @brief Yukawa interaction
     *
     * L_Yukawa = -y ψ̄_L φ ψ_R + h.c.
     *
     * After EWSB (φ → v):
     * L_mass = -(yv/√2) ψ̄ ψ = -m ψ̄ ψ
     *
     * Fermion mass: m = yv/√2
     */
    static std::string yukawaInteraction() {
        return "Yukawa coupling:\n"
               "\n"
               "L_Y = -y ψ̄_L φ ψ_R + h.c.\n"
               "\n"
               "After EWSB (⟨φ⟩ = v/√2):\n"
               "m = yv/√2\n"
               "\n"
               "Fermion mass proportional to Yukawa coupling!";
    }

    /**
     * @brief Why direct mass term forbidden
     *
     * m ψ̄ ψ = m(ψ̄_L ψ_R + ψ̄_R ψ_L)
     *
     * ψ_L and ψ_R have different SU(2) quantum numbers
     * Violates gauge invariance!
     */
    static std::string whyNoDirectMass() {
        return "Why no direct fermion mass?\n"
               "\n"
               "Direct: m ψ̄ ψ = m(ψ̄_L ψ_R + ψ̄_R ψ_L)\n"
               "\n"
               "Problem:\n"
               "ψ_L: SU(2) doublet\n"
               "ψ_R: SU(2) singlet\n"
               "\n"
               "ψ̄_L ψ_R violates SU(2) gauge invariance!\n"
               "\n"
               "Higgs doublet allows: ψ̄_L φ ψ_R (SU(2) singlet)";
    }

    /**
     * @brief Yukawa couplings (Standard Model)
     *
     * y_e = m_e√2/v ≈ 2×10⁻⁶ (electron)
     * y_μ = m_μ√2/v ≈ 4×10⁻⁴ (muon)
     * y_τ = m_τ√2/v ≈ 7×10⁻³ (tau)
     * y_t = m_t√2/v ≈ 1.0 (top quark)
     *
     * Hierarchy problem: why y_e << y_t ?
     */
    static double yukawaCoupling(double fermion_mass) {
        double v = 246.0;  // GeV
        return fermion_mass * std::sqrt(2.0) / v;
    }

    static std::map<std::string, double> standardModelYukawas() {
        return {
            {"electron", 2.9e-6},
            {"muon", 6.1e-4},
            {"tau", 1.0e-2},
            {"up", 1.2e-5},
            {"charm", 7.3e-3},
            {"top", 0.995},       // ~ 1 !
            {"down", 2.7e-5},
            {"strange", 5.5e-4},
            {"bottom", 2.4e-2}
        };
    }

    /**
     * @brief Top quark special
     *
     * y_t ≈ 1 (O(1) Yukawa!)
     * Strong coupling to Higgs
     *
     * Dominates: Higgs production, EWSB dynamics
     */
    static std::string topQuarkSpecial() {
        return "Top quark Yukawa:\n"
               "\n"
               "y_t ≈ 1.0 (O(1) coupling!)\n"
               "m_t ≈ 173 GeV ≈ v/√2\n"
               "\n"
               "Top is special:\n"
               "- Strongest Higgs coupling\n"
               "- Dominates Higgs production at LHC\n"
               "- Vacuum stability depends on y_t\n"
               "\n"
               "Why is top so heavy? Unsolved!";
    }

    /**
     * @brief Fermion mass hierarchy problem
     *
     * Why m_e/m_t ~ 10⁻⁶ ?
     * Equivalently: why y_e/y_t ~ 10⁻⁶ ?
     *
     * Not explained by Standard Model
     */
    static std::string hierarchyProblem() {
        return "Fermion mass hierarchy:\n"
               "\n"
               "m_e ~ 0.5 MeV\n"
               "m_t ~ 173 GeV\n"
               "\n"
               "Ratio: m_e/m_t ~ 3×10⁻⁶\n"
               "\n"
               "Why such huge hierarchy?\n"
               "Yukawa couplings are free parameters in SM\n"
               "Suggests deeper structure (flavor physics)";
    }
};

/**
 * @class HiggsBosonProperties
 * @brief Properties of the physical Higgs boson
 *
 * Discovered July 4, 2012 at LHC (ATLAS & CMS)
 * Nobel Prize 2013 (Higgs, Englert)
 */
class HiggsBosonProperties {
public:
    /**
     * @brief Higgs mass
     *
     * m_h = 125.10 ± 0.14 GeV (2022 PDG)
     *
     * Not predicted by SM! Free parameter
     * Determines quartic coupling: λ = m_h²/(2v²)
     */
    static double mass() {
        return 125.10;  // GeV
    }

    static std::string massSignificance() {
        return "Higgs mass: m_h = 125.10 GeV\n"
               "\n"
               "Not predicted by SM (free parameter)\n"
               "Determines λ = m_h²/(2v²) ≈ 0.13\n"
               "\n"
               "Implications:\n"
               "- Vacuum is metastable (near critical!)\n"
               "- Close to scale invariance\n"
               "- Suggests new physics at TeV-PeV?";
    }

    /**
     * @brief Higgs decay modes
     *
     * Primary decays:
     * h → bb̄ (58%)
     * h → WW* (21%)
     * h → ττ (6.3%)
     * h → ZZ* (2.6%)
     * h → γγ (0.23%) - discovery channel!
     */
    static std::map<std::string, double> branchingRatios() {
        return {
            {"bb", 0.58},
            {"WW*", 0.21},
            {"gg", 0.08},      // Via top loop
            {"tau tau", 0.063},
            {"cc", 0.029},
            {"ZZ*", 0.026},
            {"gamma gamma", 0.0023},  // Loop-induced
            {"Z gamma", 0.0015},
            {"mu mu", 0.00022}
        };
    }

    /**
     * @brief Higgs width
     *
     * Γ_h ≈ 4.1 MeV (very narrow!)
     * τ_h ≈ 1.6×10⁻²² s
     *
     * Too short to directly measure width
     * Use off-shell production
     */
    static double width() {
        return 0.0041;  // GeV = 4.1 MeV
    }

    /**
     * @brief Higgs couplings to particles
     *
     * g_hXX ~ m_X/v (for fermions and gauge bosons)
     *
     * Heavier particles couple more strongly!
     */
    static std::string couplingPattern() {
        return "Higgs couplings:\n"
               "\n"
               "Fermions: g_hff̄ = -y_f = -m_f√2/v\n"
               "Gauge bosons: g_hVV = 2m_V²/v\n"
               "\n"
               "Coupling proportional to mass!\n"
               "\n"
               "g_htt̄/g_hbb̄ = m_t/m_b ≈ 40\n"
               "Top couples strongest to Higgs";
    }

    /**
     * @brief Discovery channels at LHC
     *
     * h → γγ (clear signature, low background)
     * h → ZZ* → 4ℓ (gold-plated channel)
     *
     * July 4, 2012: 5σ discovery!
     */
    static std::string discoveryChannels() {
        return "Higgs discovery (July 4, 2012):\n"
               "\n"
               "Two channels:\n"
               "1. h → γγ (diphoton bump at 125 GeV)\n"
               "2. h → ZZ* → 4ℓ (four lepton events)\n"
               "\n"
               "Combined significance: > 5σ\n"
               "\n"
               "ATLAS & CMS independent confirmation\n"
               "Nobel Prize 2013 (Higgs, Englert)";
    }

    /**
     * @brief Higgs production at LHC
     *
     * Main modes:
     * - Gluon fusion gg → h (87%) [via top loop]
     * - Vector boson fusion qq → qqh (7%)
     * - Vh (Higgsstrahlung) (4%)
     * - tt̄h (1%)
     */
    static std::string productionModes() {
        return "Higgs production at LHC:\n"
               "\n"
               "gg → h (87%): gluon fusion via top loop\n"
               "qq → qqh (7%): vector boson fusion (VBF)\n"
               "Vh (4%): associated production (W/Z + h)\n"
               "tt̄h (1%): top pair + Higgs\n"
               "\n"
               "Total cross section ~ 50 pb at √s = 13 TeV";
    }

    /**
     * @brief Vacuum stability
     *
     * With m_h = 125 GeV, SM vacuum is metastable
     * Lifetime >> age of universe (OK!)
     *
     * But hints at new physics?
     */
    static std::string vacuumStability() {
        return "Vacuum stability:\n"
               "\n"
               "m_h = 125 GeV → Higgs quartic coupling λ runs negative!\n"
               "λ(Λ) < 0 at Λ ~ 10¹⁰-10¹³ GeV\n"
               "\n"
               "SM vacuum is metastable\n"
               "Tunneling time: τ >> 10²⁰⁰ years (safe!)\n"
               "\n"
               "But: we live near critical stability bound\n"
               "Suggests new physics at high scales?";
    }
};

} // namespace physics::advanced::gauge_theory

#endif // PHYSICS_ADVANCED_GAUGE_THEORY_HIGGS_MECHANISM_HPP
