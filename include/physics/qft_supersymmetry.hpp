#ifndef PHYSICS_ADVANCED_QFT_SUPERSYMMETRY_HPP
#define PHYSICS_ADVANCED_QFT_SUPERSYMMETRY_HPP

#include <string>
#include <vector>
#include <cmath>

/**
 * @file supersymmetry.hpp
 * @brief Supersymmetry (SUSY) - beyond Standard Model
 *
 * Implements:
 * - Supersymmetric particle spectrum
 * - SUSY algebra
 * - Soft SUSY breaking
 * - R-parity
 * - Minimal Supersymmetric Standard Model (MSSM)
 */

namespace physics::advanced::qft {

/**
 * @class Supersymmetry
 * @brief Symmetry between fermions and bosons
 *
 * SUSY: Q|fermion⟩ = |boson⟩
 *       Q|boson⟩ = |fermion⟩
 *
 * Every Standard Model particle has supersymmetric partner
 * differing by spin-1/2
 */
class Supersymmetry {
public:
    /**
     * @brief SUSY algebra commutation relation
     *
     * {Q_α, Q̄_β} = 2(γ^μ)_αβ P_μ
     *
     * where Q is SUSY generator (spinor charge)
     *       P_μ is momentum operator
     */
    static std::string susyAlgebra() {
        return "{Q_α, Q̄_β} = 2(γ^μ)_αβ P_μ\n"
               "Anticommutator of SUSY charges equals momentum";
    }

    /**
     * @brief Mass relation in unbroken SUSY
     *
     * M(boson) = M(fermion)
     *
     * Superpartners have identical masses
     * (not observed → SUSY must be broken)
     */
    static bool isUnbroken(double boson_mass, double fermion_mass,
                          double tolerance = 1e-3) {
        return std::abs(boson_mass - fermion_mass) < tolerance;
    }

    /**
     * @brief Supersymmetric multiplets
     *
     * Chiral multiplet: (fermion, 2 scalars)
     * Vector multiplet: (boson, fermion)
     */
    static std::string multipletStructure() {
        return "Chiral multiplet: (ψ, φ, φ*) - fermion + 2 complex scalars\n"
               "Vector multiplet: (A_μ, λ) - gauge boson + gaugino";
    }

    /**
     * @brief Hierarchy problem solution
     *
     * SUSY stabilizes Higgs mass against quantum corrections
     * Fermionic loops cancel bosonic loops
     */
    static std::string hierarchyProblemSolution() {
        return "Fermionic corrections: Δm² ~ +Λ²\n"
               "Bosonic corrections: Δm² ~ -Λ²\n"
               "Cancel if SUSY exact → Natural Higgs mass";
    }
};

/**
 * @class SuperparticleSpectrum
 * @brief Supersymmetric partners of Standard Model particles
 */
class SuperparticleSpectrum {
public:
    /**
     * @brief Squarks (scalar partners of quarks)
     *
     * ũ, d̃, c̃, s̃, t̃, b̃
     *
     * Spin-0 (scalar)
     * Same charge and color as quarks
     */
    static std::vector<std::string> squarks() {
        return {"ũ_L", "ũ_R", "d̃_L", "d̃_R",
                "c̃_L", "c̃_R", "s̃_L", "s̃_R",
                "t̃_1", "t̃_2", "b̃_1", "b̃_2"};  // Mass eigenstates
    }

    /**
     * @brief Sleptons (scalar partners of leptons)
     *
     * ẽ, μ̃, τ̃, ν̃_e, ν̃_μ, ν̃_τ
     *
     * Spin-0 (scalar)
     */
    static std::vector<std::string> sleptons() {
        return {"ẽ_L", "ẽ_R", "μ̃_L", "μ̃_R", "τ̃_1", "τ̃_2",
                "ν̃_e", "ν̃_μ", "ν̃_τ"};
    }

    /**
     * @brief Gauginos (fermionic partners of gauge bosons)
     *
     * Photino (γ̃), Zino (Z̃), Winos (W̃±, W̃⁰)
     * Gluino (g̃)
     *
     * Spin-1/2 (fermion)
     */
    static std::vector<std::string> gauginos() {
        return {"γ̃ (photino)", "Z̃ (zino)", "W̃± (wino)", "g̃ (gluino)"};
    }

    /**
     * @brief Higgsinos (fermionic partners of Higgs bosons)
     *
     * H̃_u, H̃_d (up-type and down-type)
     *
     * Spin-1/2
     */
    static std::vector<std::string> higgsinos() {
        return {"H̃_u⁰", "H̃_u⁺", "H̃_d⁰", "H̃_d⁻"};
    }

    /**
     * @brief Mass eigenstates (after EWSB)
     *
     * Neutralinos: χ̃₁⁰, χ̃₂⁰, χ̃₃⁰, χ̃₄⁰ (mix of γ̃, Z̃, H̃⁰)
     * Charginos: χ̃₁±, χ̃₂± (mix of W̃±, H̃±)
     */
    static std::vector<std::string> massEigenstates() {
        return {"χ̃₁⁰ (LSP)", "χ̃₂⁰", "χ̃₃⁰", "χ̃₄⁰",
                "χ̃₁±", "χ̃₂±", "t̃₁", "t̃₂", "b̃₁", "b̃₂", "τ̃₁", "τ̃₂"};
    }

    /**
     * @brief Stop quark (t̃) mass eigenvalues
     *
     * m²_t̃ = m²_Q + m²_t ± Δ
     *
     * where Δ includes mixing
     */
    static std::pair<double, double> stopMasses(double mQ, double mt,
                                                double mixing_angle) {
        double m2_avg = (mQ * mQ + mt * mt) / 2.0;
        double delta = std::abs(mQ * mQ - mt * mt) / 2.0 *
                      std::cos(2.0 * mixing_angle);

        double m1 = std::sqrt(m2_avg - delta);
        double m2 = std::sqrt(m2_avg + delta);

        return {m1, m2};  // t̃₁ (lighter), t̃₂ (heavier)
    }
};

/**
 * @class RParity
 * @brief R-parity conservation
 *
 * R = (-1)^(3B + L + 2S)
 *
 * where B = baryon number, L = lepton number, S = spin
 */
class RParity {
public:
    /**
     * @brief Calculate R-parity
     *
     * R = (-1)^(3B + L + 2S)
     *
     * Standard Model particles: R = +1
     * Superpartners: R = -1
     */
    static int calculate(int baryon_number, int lepton_number, double spin) {
        int exponent = 3 * baryon_number + lepton_number +
                      static_cast<int>(2.0 * spin);

        return (exponent % 2 == 0) ? +1 : -1;
    }

    /**
     * @brief Check if R-parity conserved
     *
     * If conserved:
     * - Superpartners produced in pairs
     * - Lightest SUSY particle (LSP) is stable
     * - LSP is dark matter candidate (χ̃₁⁰)
     */
    static bool isConserved() {
        return true;  // Assumed in MSSM
    }

    /**
     * @brief Lightest supersymmetric particle (LSP)
     *
     * If R-parity conserved, LSP is stable
     * Candidate: neutralino χ̃₁⁰
     *
     * Dark matter relic density requires:
     * Ω_χ h² ≈ 0.12
     */
    static std::string lspDarkMatter() {
        return "LSP (χ̃₁⁰) is dark matter candidate\n"
               "Thermal relic: σv ~ 3×10⁻²⁶ cm³/s for Ω_DM h² ≈ 0.12";
    }

    /**
     * @brief R-parity violating interactions
     *
     * If R-parity violated:
     * - Single sparticle production possible
     * - LSP decays (no dark matter)
     * - Potential proton decay
     */
    static std::string rParityViolation() {
        return "R-parity violating terms:\n"
               "λ LLĒ (lepton number violation)\n"
               "λ' LQD̄ (lepton & baryon number violation)\n"
               "λ'' ŪD̄D̄ (baryon number violation)";
    }
};

/**
 * @class SUSYBreaking
 * @brief Soft SUSY breaking mechanisms
 *
 * SUSY must be broken to explain mass differences
 * Soft breaking: no quadratic divergences
 */
class SUSYBreaking {
public:
    /**
     * @brief Soft breaking terms
     *
     * L_soft = -m²|φ|² - (Aλφ³ + B μφ² + h.c.) - (½M_a λ_a λ_a + h.c.)
     *
     * where m² are scalar masses
     *       A are trilinear couplings
     *       M_a are gaugino masses
     */
    static std::string softBreakingLagrangian() {
        return "L_soft = -m²|φ|² - Aλφ³ - Bμφ² - ½M_a λ_a²\n"
               "Soft masses break SUSY without reintroducing hierarchy problem";
    }

    /**
     * @brief SUSY breaking scale
     *
     * M_SUSY ~ 1 TeV (TeV scale SUSY)
     *
     * Naturalness suggests m_sparticle ~ O(TeV)
     */
    static double susyBreakingScale() {
        return 1000.0;  // GeV (TeV scale)
    }

    /**
     * @brief Gaugino mass unification
     *
     * M_1 : M_2 : M_3 = α_1 : α_2 : α_3 at GUT scale
     *
     * At weak scale:
     * M_1 : M_2 : M_3 ≈ 1 : 2 : 7
     */
    static std::vector<double> gauginoMassRatios() {
        return {1.0, 2.0, 7.0};  // M_1, M_2, M_3 (bino, wino, gluino)
    }

    /**
     * @brief Gravity mediation (mSUGRA)
     *
     * SUSY breaking in hidden sector
     * Mediated by gravity
     *
     * m_soft ~ F/M_Pl where F is SUSY breaking scale
     */
    static double gravityMediatedMass(double F_breaking, double M_Planck) {
        return F_breaking / M_Planck;  // GeV
    }

    /**
     * @brief Gauge mediation (GMSB)
     *
     * SUSY breaking mediated by gauge interactions
     * Lighter grav itino (~eV-keV)
     */
    static std::string gaugeMediatedBreaking() {
        return "Gauge mediated SUSY breaking (GMSB)\n"
               "Gravitino LSP: m_G̃ ~ F/√(F_messenger)\n"
               "Typical: m_G̃ ~ keV-MeV";
    }
};

/**
 * @class MSSM
 * @brief Minimal Supersymmetric Standard Model
 */
class MSSM {
public:
    /**
     * @brief MSSM particle content
     *
     * 2 Higgs doublets (H_u, H_d) required
     * 5 physical Higgs bosons: h⁰, H⁰, A⁰, H±
     */
    static std::vector<std::string> higgsSpectrum() {
        return {"h⁰ (light CP-even)",
                "H⁰ (heavy CP-even)",
                "A⁰ (CP-odd)",
                "H± (charged)"};
    }

    /**
     * @brief MSSM parameters (at weak scale)
     *
     * ~105 free parameters after SUSY breaking!
     * mSUGRA reduces to 5: m_0, m_1/2, A_0, tan β, sign(μ)
     */
    static int numberOfParameters() {
        return 105;  // Full MSSM
    }

    /**
     * @brief tan β (ratio of Higgs VEVs)
     *
     * tan β = v_u / v_d
     *
     * Typically: 1 < tan β < 60
     */
    static double tanBeta() {
        return 10.0;  // Typical value
    }

    /**
     * @brief μ parameter (Higgsino mass)
     *
     * |μ| ~ M_SUSY (naturalness)
     */
    static double muParameter() {
        return 200.0;  // GeV (typical)
    }

    /**
     * @brief Gauge coupling unification in MSSM
     *
     * α_1 = α_2 = α_3 at M_GUT ~ 2×10¹⁶ GeV
     *
     * SUSY improves unification compared to SM!
     */
    static double gutScale() {
        return 2.0e16;  // GeV
    }

    /**
     * @brief Lightest Higgs mass prediction
     *
     * m_h < m_Z at tree level
     * Radiative corrections: m_h ~ 125 GeV (consistent with discovery!)
     *
     * Requires heavy stops (~ TeV)
     */
    static double lightestHiggsMass() {
        return 125.0;  // GeV (observed)
    }
};

/**
 * @class SUSYPhenomenology
 * @brief Collider signatures of SUSY
 */
class SUSYPhenomenology {
public:
    /**
     * @brief Missing energy signature
     *
     * pp → sparticles → ... → χ̃₁⁰ + X
     *
     * LSP escapes detector → large E_T^miss
     */
    static std::string missingEnergy() {
        return "Missing transverse energy: E_T^miss\n"
               "From LSP (χ̃₁⁰) escaping detector\n"
               "Signature: jets + leptons + E_T^miss";
    }

    /**
     * @brief Squark/gluino production at LHC
     *
     * pp → q̃q̃*, g̃g̃, q̃g̃
     *
     * Large cross sections if masses accessible
     */
    static std::string strongProduction() {
        return "Strong production: pp → q̃q̃, g̃g̃\n"
               "σ(g̃g̃) ~ pb for m_g̃ ~ 1 TeV\n"
               "σ(q̃q̃) ~ 0.1 pb for m_q̃ ~ 1 TeV";
    }

    /**
     * @brief Current mass limits (LHC Run 2)
     *
     * m_g̃ > 2 TeV
     * m_q̃ > 1-2 TeV (depending on decay mode)
     * m_χ̃₁⁰ > 100-400 GeV
     */
    static std::vector<std::pair<std::string, double>> massLimits() {
        return {
            {"gluino", 2000.0},      // GeV
            {"squark (1st gen)", 1800.0},
            {"stop", 1000.0},
            {"neutralino", 400.0},
            {"chargino", 600.0}
        };
    }

    /**
     * @brief Fine-tuning problem
     *
     * If sparticles too heavy (> few TeV), fine-tuning required
     * to keep Higgs mass at 125 GeV
     *
     * Δ ~ (m_sparticle/m_Z)²
     */
    static double fineTuning(double sparticle_mass) {
        double m_Z = 91.2;  // GeV
        return (sparticle_mass / m_Z) * (sparticle_mass / m_Z);
    }

    /**
     * @brief Natural SUSY scenario
     *
     * Light stops, higgsinos, gluino (< 1-2 TeV)
     * Heavy squarks of 1st/2nd generation (> 2-3 TeV)
     *
     * Reduces fine-tuning while evading LHC limits
     */
    static std::string naturalSUSY() {
        return "Natural SUSY: light t̃, χ̃, g̃ (< 2 TeV)\n"
               "Heavy 1st/2nd gen squarks (> 3 TeV)\n"
               "Minimal fine-tuning: Δ < 100";
    }
};

} // namespace physics::advanced::qft

#endif // PHYSICS_ADVANCED_QFT_SUPERSYMMETRY_HPP
