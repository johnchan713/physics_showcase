#ifndef PHYSICS_ADVANCED_GAUGE_THEORY_RUNNING_COUPLINGS_HPP
#define PHYSICS_ADVANCED_GAUGE_THEORY_RUNNING_COUPLINGS_HPP

#include <cmath>
#include <string>
#include <vector>

/**
 * @file running_couplings.hpp
 * @brief Running couplings and renormalization group equations
 *
 * Implements:
 * - Renormalization group (RG) equations
 * - Running of coupling constants with energy
 * - QED: α(Q²) increases (screening)
 * - QCD: α_s(Q²) decreases (asymptotic freedom)
 * - Electroweak running
 * - Coupling unification in GUTs
 * - Experimental tests
 */

namespace physics::advanced::gauge_theory {

/**
 * @class RenormalizationGroup
 * @brief Renormalization group equations (RGE)
 *
 * Coupling constants "run" with energy scale Q²
 * β-function governs running: dg/d log Q² = β(g)
 */
class RenormalizationGroup {
public:
    /**
     * @brief Beta function
     *
     * β(g) = dg/d log μ²
     *
     * where μ is renormalization scale
     *
     * β > 0: coupling increases (QED)
     * β < 0: coupling decreases (QCD - asymptotic freedom)
     */
    static std::string betaFunction() {
        return "Beta function:\n"
               "\n"
               "β(g) = dg/d log μ²\n"
               "\n"
               "β(g) > 0: coupling increases with energy (QED)\n"
               "β(g) < 0: coupling decreases with energy (QCD)\n"
               "\n"
               "Determines running of coupling constant";
    }

    /**
     * @brief One-loop beta function
     *
     * β(g) = -b₀ g³/(16π²) + O(g⁵)
     *
     * where b₀ depends on particle content
     */
    static std::string oneLoopBeta() {
        return "One-loop beta function:\n"
               "\n"
               "β(g) = -b₀ g³/(16π²)\n"
               "\n"
               "b₀ = (11C_A - 4T_F n_f)/3 (non-Abelian)\n"
               "b₀ = -4T_F n_f/3 (Abelian, e.g. QED)\n"
               "\n"
               "C_A: Casimir of adjoint (N for SU(N))\n"
               "T_F: fermion representation (1/2 for fund.)\n"
               "n_f: number of fermion flavors";
    }

    /**
     * @brief Running coupling solution (one-loop)
     *
     * α(Q²) = α(μ²) / [1 - (b₀α(μ²)/2π) log(Q²/μ²)]
     *
     * where α = g²/(4π)
     */
    static double runningCoupling(double alpha_ref, double Q2, double mu2,
                                   double b0) {
        double log_ratio = std::log(Q2 / mu2);
        return alpha_ref / (1.0 - (b0 * alpha_ref / (2.0 * M_PI)) * log_ratio);
    }

    /**
     * @brief Landau pole
     *
     * If β > 0, coupling diverges at high energy:
     * Λ² = μ² exp(2π/b₀α(μ²))
     *
     * QED: Landau pole at ~10²⁸⁶ GeV (way beyond Planck!)
     */
    static std::string landauPole() {
        return "Landau pole:\n"
               "\n"
               "If β > 0, coupling diverges at Λ:\n"
               "α(Λ²) → ∞\n"
               "\n"
               "QED: Λ_Landau ~ 10²⁸⁶ GeV\n"
               "(Much higher than Planck scale!)\n"
               "\n"
               "Suggests QED not valid to arbitrarily high energy\n"
               "But no practical problem";
    }

    /**
     * @brief Asymptotic freedom
     *
     * If β < 0, coupling → 0 as Q² → ∞
     *
     * QCD: α_s(Q²) → 0 (Nobel Prize 2004)
     * Allows perturbative calculations at high energy
     */
    static std::string asymptoticFreedom() {
        return "Asymptotic freedom (Nobel 2004):\n"
               "\n"
               "β < 0 → α_s(Q²) → 0 as Q² → ∞\n"
               "\n"
               "High energy: quarks nearly free\n"
               "Low energy: quarks confined\n"
               "\n"
               "Gross, Politzer, Wilczek (2004)\n"
               "Essential for QCD as strong force theory";
    }
};

/**
 * @class QEDRunning
 * @brief Running of electromagnetic coupling α
 *
 * α increases with energy (vacuum polarization screening)
 */
class QEDRunning {
public:
    /**
     * @brief QED beta function (one-loop)
     *
     * β(e) = e³/(12π²) n_f
     *
     * where n_f is number of charged fermions
     *
     * Positive → α increases with Q²
     */
    static double betaCoefficient(int n_fermions) {
        return -4.0 * n_fermions / 3.0;  // b₀ for QED
    }

    /**
     * @brief Running of α (one-loop)
     *
     * α(Q²) = α(m_e²) / [1 - (α(m_e)/3π) Σ_f Q_f² log(Q²/m_f²)]
     *
     * Sum over charged fermions with mass m_f < Q
     */
    static double runningAlpha(double Q2) {
        double alpha_me = 1.0 / 137.036;  // At electron mass
        double me2 = 0.511e-3 * 0.511e-3;  // (MeV)² in GeV²

        if (Q2 <= me2) return alpha_me;

        // Simple approximation (e, μ, τ)
        double delta = (alpha_me / (3.0 * M_PI)) * 3.0 * std::log(Q2 / me2);

        return alpha_me / (1.0 - delta);
    }

    /**
     * @brief α at various scales
     *
     * α(m_e) ≈ 1/137.036 (low energy)
     * α(m_Z) ≈ 1/128 (Z boson mass)
     * α(m_t) ≈ 1/127 (top mass)
     */
    static double alphaAtElectronMass() {
        return 1.0 / 137.036;
    }

    static double alphaAtZMass() {
        return 1.0 / 127.9;  // MS-bar scheme
    }

    static double alphaAtTopMass() {
        return 1.0 / 127.0;
    }

    /**
     * @brief Vacuum polarization
     *
     * Virtual e⁺e⁻ pairs screen electric charge
     * Effective charge increases at short distances
     *
     * "Charge screening" (opposite of QCD antiscreening!)
     */
    static std::string vacuumPolarization() {
        return "Vacuum polarization in QED:\n"
               "\n"
               "Virtual e⁺e⁻ pairs screen charge\n"
               "At large r: see screened charge\n"
               "At small r: see bare charge (larger!)\n"
               "\n"
               "α(Q²) increases with Q²\n"
               "\n"
               "Opposite of QCD antiscreening!";
    }

    /**
     * @brief Experimental verification
     *
     * Measured at LEP, Tevatron, LHC
     * α(m_Z) = 1/127.9 ± 0.1
     *
     * Confirms QED running
     */
    static std::string experimentalTests() {
        return "QED running measurements:\n"
               "\n"
               "α(m_e) = 1/137.036 (low energy)\n"
               "α(m_Z) = 1/127.9 ± 0.1 (LEP, Z pole)\n"
               "\n"
               "~7% increase from m_e to m_Z\n"
               "Confirms vacuum polarization effects";
    }
};

/**
 * @class QCDRunning
 * @brief Running of strong coupling α_s
 *
 * α_s decreases with energy (asymptotic freedom)
 * Increases at low energy (confinement)
 */
class QCDRunning {
public:
    /**
     * @brief QCD beta function (one-loop)
     *
     * β₀ = (33 - 2n_f) / 3
     *
     * For n_f = 6 flavors: β₀ = 7
     * β₀ > 0 → asymptotic freedom!
     */
    static double betaCoefficient(int n_flavors) {
        return (33.0 - 2.0 * n_flavors) / 3.0;
    }

    /**
     * @brief Running of α_s (one-loop)
     *
     * α_s(Q²) = α_s(μ²) / [1 + (β₀α_s(μ²)/2π) log(Q²/μ²)]
     *
     * Note: plus sign! (β < 0 for QCD)
     */
    static double runningAlphaS(double Q2, double mu2 = 91.2 * 91.2,
                                 double alpha_s_mu = 0.118, int n_f = 5) {
        double beta0 = betaCoefficient(n_f);
        double log_ratio = std::log(Q2 / mu2);

        return alpha_s_mu / (1.0 + (beta0 * alpha_s_mu / (2.0 * M_PI)) * log_ratio);
    }

    /**
     * @brief α_s at various scales
     *
     * α_s(m_Z) ≈ 0.1179 ± 0.0009 (world average 2022)
     * α_s(m_τ) ≈ 0.33 (low energy)
     * α_s(M_Pl) ~ 0.01 (Planck scale)
     */
    static double alphaSAtZMass() {
        return 0.1179;  // PDG 2022
    }

    static double alphaSAtTauMass() {
        return 0.33;
    }

    static double alphaSAtTopMass() {
        return 0.108;
    }

    /**
     * @brief QCD scale Λ_QCD
     *
     * α_s(Q²) ≈ 2π / [β₀ log(Q²/Λ²_QCD)]
     *
     * Λ_QCD ≈ 200-300 MeV (n_f = 5)
     *
     * Scale where α_s ~ 1 (strong coupling)
     */
    static double lambdaQCD() {
        return 0.213;  // GeV (n_f = 5, MS-bar)
    }

    static std::string qcdScale() {
        return "QCD scale Λ_QCD:\n"
               "\n"
               "α_s(Q²) ~ 2π/[β₀ log(Q²/Λ²_QCD)]\n"
               "\n"
               "Λ_QCD ≈ 213 MeV (n_f = 5)\n"
               "\n"
               "Scale where α_s ~ 1\n"
               "Confinement scale!";
    }

    /**
     * @brief Asymptotic freedom condition
     *
     * β₀ > 0 requires: 33 - 2n_f > 0
     * → n_f < 16.5
     *
     * Asymptotic freedom only for n_f ≤ 16 flavors
     * (SM has 6 flavors → safe!)
     */
    static std::string asymptoticFreedomCondition() {
        return "Asymptotic freedom condition:\n"
               "\n"
               "β₀ = (33 - 2n_f)/3 > 0\n"
               "Requires: n_f < 16.5\n"
               "\n"
               "Standard Model: n_f = 6 ✓\n"
               "\n"
               "Too many fermions → no asymptotic freedom!\n"
               "(Would need different theory)";
    }

    /**
     * @brief Antiscreening in QCD
     *
     * Gluon self-interactions cause antiscreening
     * Color charge appears larger at large distances!
     *
     * Opposite of QED screening
     */
    static std::string antiscreening() {
        return "QCD antiscreening:\n"
               "\n"
               "Gluon self-coupling → antiscreening\n"
               "Color charge increases at large r\n"
               "\n"
               "Short distance: weak coupling (α_s small)\n"
               "Long distance: strong coupling (confinement)\n"
               "\n"
               "Opposite of QED!";
    }

    /**
     * @brief Experimental measurements
     *
     * α_s(m_Z) measured to ~1% precision
     * Best determined from:
     * - e⁺e⁻ → hadrons (R ratio)
     * - τ decays
     * - Deep inelastic scattering
     * - Lattice QCD
     */
    static std::string experimentalMeasurements() {
        return "α_s(m_Z) measurements:\n"
               "\n"
               "World average: 0.1179 ± 0.0009\n"
               "\n"
               "Methods:\n"
               "- τ decays: 0.1171 ± 0.0014\n"
               "- e⁺e⁻ → hadrons: 0.1198 ± 0.0015\n"
               "- Deep inelastic: 0.1156 ± 0.0021\n"
               "- Lattice QCD: 0.1184 ± 0.0012\n"
               "\n"
               "All consistent with running!";
    }
};

/**
 * @class ElectroweakRunning
 * @brief Running of electroweak couplings
 *
 * g (SU(2)), g' (U(1)) both run
 * sin²θ_W runs with energy
 */
class ElectroweakRunning {
public:
    /**
     * @brief Running of sin²θ_W
     *
     * sin²θ_W(Q²) varies with scale
     *
     * sin²θ_W(m_Z) ≈ 0.231 (MS-bar)
     * sin²θ_W(0) ≈ 0.239 (Thomson limit)
     */
    static double sin2ThetaWAtZMass() {
        return 0.23122;  // MS-bar scheme
    }

    static double sin2ThetaWAtZero() {
        return 0.23867;  // Thomson limit (Q² → 0)
    }

    /**
     * @brief Running of g and g'
     *
     * Both couplings run logarithmically
     * g increases, g' increases
     *
     * But sin²θ_W = g'²/(g² + g'²) runs slowly
     */
    static std::string couplingRunning() {
        return "Electroweak coupling running:\n"
               "\n"
               "g(Q²): SU(2) coupling increases\n"
               "g'(Q²): U(1) coupling increases\n"
               "\n"
               "sin²θ_W(Q²) = g'²/(g² + g'²) runs slowly\n"
               "\n"
               "sin²θ_W(0) ≈ 0.239\n"
               "sin²θ_W(m_Z) ≈ 0.231";
    }

    /**
     * @brief Scheme dependence
     *
     * Different renormalization schemes:
     * - MS-bar: sin²θ_W(m_Z) = 0.23122
     * - On-shell: sin²θ_W = 1 - m_W²/m_Z² = 0.2229
     */
    static std::string schemeDependence() {
        return "Scheme dependence:\n"
               "\n"
               "MS-bar: sin²θ_W(m_Z) = 0.23122\n"
               "On-shell: sin²θ_W = 1 - m_W²/m_Z² = 0.2229\n"
               "\n"
               "Different definitions!\n"
               "Must specify scheme when quoting sin²θ_W";
    }

    /**
     * @brief ρ parameter running
     *
     * ρ = m_W²/(m_Z² cos²θ_W)
     *
     * Tree level: ρ = 1
     * Loop corrections: Δρ ~ (G_F m_t²)/(8π²√2)
     */
    static double rhoParameter() {
        return 1.00038;  // Including radiative corrections
    }

    static std::string rhoCorrections() {
        return "ρ parameter corrections:\n"
               "\n"
               "Tree level: ρ = 1\n"
               "\n"
               "One-loop: Δρ ~ (3G_F m_t²)/(8π²√2)\n"
               "Δρ ≈ 0.01 (from m_t = 173 GeV)\n"
               "\n"
               "Measured: ρ = 1.00038 ± 0.00020\n"
               "Consistent with top quark contribution!";
    }
};

/**
 * @class GrandUnification
 * @brief Grand Unified Theories (GUTs)
 *
 * SU(3) × SU(2) × U(1) → SU(5), SO(10), E₆, ...
 * Couplings unify at M_GUT ~ 10¹⁶ GeV
 */
class GrandUnification {
public:
    /**
     * @brief Coupling unification
     *
     * In MSSM (with SUSY):
     * α₁⁻¹, α₂⁻¹, α₃⁻¹ meet at M_GUT ~ 2×10¹⁶ GeV
     *
     * In SM: couplings miss! (Hint for SUSY)
     */
    static std::string couplingUnification() {
        return "Gauge coupling unification:\n"
               "\n"
               "Standard Model: couplings nearly meet\n"
               "Miss by ~2% at M ~ 10¹⁴ GeV\n"
               "\n"
               "MSSM (with SUSY): couplings unify!\n"
               "M_GUT ≈ 2×10¹⁶ GeV\n"
               "α_GUT ≈ 1/25\n"
               "\n"
               "Strong evidence for SUSY + GUT!";
    }

    /**
     * @brief GUT scale
     *
     * M_GUT ~ 2×10¹⁶ GeV (SUSY)
     * M_GUT ~ 10¹⁴ GeV (SM, approximate)
     */
    static double gutScale() {
        return 2.0e16;  // GeV (MSSM)
    }

    static double gutCoupling() {
        return 1.0 / 25.0;  // α_GUT
    }

    /**
     * @brief GUT groups
     *
     * Minimal: SU(5) (Georgi-Glashow, 1974)
     * Larger: SO(10), E₆
     *
     * SU(5): one generation fits in 5̄ + 10
     */
    static std::string gutGroups() {
        return "GUT gauge groups:\n"
               "\n"
               "SU(5): minimal GUT (Georgi-Glashow 1974)\n"
               "  One generation: 5̄ + 10\n"
               "  5̄ = (d̄, d̄, d̄, e⁻, νₑ)\n"
               "  10 = (u, u, u, ū, ū, ū, d̄, ē, e⁺, ν̄)\n"
               "\n"
               "SO(10): one generation in single 16\n"
               "E₆: even larger unification";
    }

    /**
     * @brief Proton decay
     *
     * GUTs predict proton decay: p → e⁺ + π⁰
     * τ_p ~ (M_GUT/m_p)⁴ M_GUT ~ 10³⁴ years
     *
     * Experimental limit: τ_p > 10³⁴ years
     * Rules out simple SU(5)!
     */
    static std::string protonDecay() {
        return "Proton decay in GUTs:\n"
               "\n"
               "GUT: B and L violated (B - L conserved)\n"
               "Allows: p → e⁺ + π⁰\n"
               "\n"
               "Minimal SU(5): τ_p ~ 10³⁰ years\n"
               "Experimental limit: τ_p > 1.6×10³⁴ years\n"
               "\n"
               "Simple SU(5) ruled out!\n"
               "SUSY GUTs predict τ_p ~ 10³⁴⁻³⁶ years (still allowed)";
    }

    /**
     * @brief Charge quantization in GUTs
     *
     * GUT: quarks and leptons in same multiplet
     * → Q_p + Q_e = 0 (automatic!)
     *
     * Explains electric charge quantization
     */
    static std::string chargeQuantization() {
        return "Charge quantization from GUT:\n"
               "\n"
               "Quarks and leptons in same multiplet\n"
               "Tr(Q) = 0 (traceless generator)\n"
               "\n"
               "→ 3Q_q + Q_e = 0\n"
               "→ Q_p = -Q_e (automatic!)\n"
               "\n"
               "Explains why |Q_p| = |Q_e| exactly";
    }

    /**
     * @brief Neutrino masses in GUTs
     *
     * Right-handed neutrinos naturally appear
     * See-saw mechanism: m_ν ~ m_D²/M_R
     *
     * M_R ~ M_GUT → m_ν ~ eV (perfect!)
     */
    static std::string neutrinoMasses() {
        return "Neutrino masses in GUTs:\n"
               "\n"
               "SO(10): right-handed neutrino in 16\n"
               "See-saw: m_ν ~ m_D²/M_R\n"
               "\n"
               "m_D ~ 100 GeV (like charged fermions)\n"
               "M_R ~ M_GUT ~ 10¹⁶ GeV\n"
               "\n"
               "→ m_ν ~ (100)²/10¹⁶ ~ 10⁻⁴ eV\n"
               "\n"
               "Naturally explains tiny neutrino masses!";
    }
};

/**
 * @class ExperimentalTests
 * @brief Precision tests of running couplings
 */
class ExperimentalTests {
public:
    /**
     * @brief LEP precision measurements
     *
     * Z pole measurements at LEP (CERN):
     * - α_s(m_Z) = 0.1196 ± 0.0030
     * - sin²θ_W(m_Z) = 0.23122 ± 0.00003
     * - ρ = 1.00037 ± 0.00023
     *
     * Test electroweak running to high precision
     */
    static std::string lepMeasurements() {
        return "LEP precision measurements:\n"
               "\n"
               "α_s(m_Z) = 0.1196 ± 0.0030\n"
               "sin²θ_W = 0.23122 ± 0.00003\n"
               "ρ = 1.00037 ± 0.00023\n"
               "\n"
               "Confirms:\n"
               "- QCD running\n"
               "- Electroweak radiative corrections\n"
               "- Top quark contribution to ρ";
    }

    /**
     * @brief Deep inelastic scattering
     *
     * HERA, SLAC: measure α_s(Q²) over wide range
     * 1 GeV² < Q² < 10,000 GeV²
     *
     * Confirms QCD running
     */
    static std::string deepInelasticScattering() {
        return "Deep inelastic scattering:\n"
               "\n"
               "Probe α_s(Q²) from Q² ~ 1 to 10⁴ GeV²\n"
               "\n"
               "HERA data:\n"
               "α_s(1.5 GeV²) ≈ 0.3\n"
               "α_s(100 GeV²) ≈ 0.12\n"
               "\n"
               "Confirms QCD running prediction!";
    }

    /**
     * @brief Comparison with theory
     *
     * α_s(m_Z) predicted from lower scales
     * Agrees with direct measurement at m_Z
     *
     * Confirms RG evolution!
     */
    static std::string theoryComparison() {
        return "Running coupling consistency:\n"
               "\n"
               "Measure α_s at low Q², evolve to m_Z\n"
               "Compare with direct measurement at Z pole\n"
               "\n"
               "Result: Perfect agreement!\n"
               "\n"
               "Confirms renormalization group equations\n"
               "Validates QCD as correct theory";
    }

    /**
     * @brief Future prospects
     *
     * - Higher precision α_s measurements
     * - Test unification at LHC/future colliders
     * - Search for deviations (new physics?)
     */
    static std::string futureProspects() {
        return "Future tests:\n"
               "\n"
               "- Precision α_s(m_Z) to 0.1% (lattice + experiment)\n"
               "- Test unification with SUSY searches\n"
               "- Look for deviations from SM running\n"
               "  (Could indicate new particles!)\n"
               "\n"
               "Any deviation → new physics discovery";
    }
};

} // namespace physics::advanced::gauge_theory

#endif // PHYSICS_ADVANCED_GAUGE_THEORY_RUNNING_COUPLINGS_HPP
