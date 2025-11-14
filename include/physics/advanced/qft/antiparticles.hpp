#ifndef PHYSICS_ADVANCED_QFT_ANTIPARTICLES_HPP
#define PHYSICS_ADVANCED_QFT_ANTIPARTICLES_HPP

#include "particle_physics.hpp"
#include <string>
#include <cmath>

/**
 * @file antiparticles.hpp
 * @brief Antiparticles and charge conjugation
 *
 * Implements:
 * - Antiparticle properties
 * - Charge conjugation (C)
 * - CP symmetry
 * - Matter-antimatter asymmetry
 * - Pair production and annihilation
 */

namespace physics::advanced::qft {

/**
 * @class Antiparticle
 * @brief Properties of antiparticles
 *
 * Every particle has an antiparticle with:
 * - Same mass
 * - Opposite electric charge
 * - Opposite baryon/lepton number
 * - Same spin
 * - Same lifetime
 */
class Antiparticle {
public:
    /**
     * @brief Create antiparticle from particle
     *
     * Flips: charge, baryon number, lepton number
     * Preserves: mass, spin, lifetime
     */
    static ParticleProperties createAntiparticle(const ParticleProperties& particle) {
        ParticleProperties anti = particle;

        // Change name and symbol
        if (particle.symbol.find("⁺") != std::string::npos) {
            anti.symbol.replace(anti.symbol.find("⁺"), 3, "⁻");
            anti.name = "Anti-" + particle.name;
        } else if (particle.symbol.find("⁻") != std::string::npos) {
            anti.symbol.replace(anti.symbol.find("⁻"), 3, "⁺");
            anti.name = "Anti-" + particle.name;
        } else {
            anti.name = "Anti-" + particle.name;
            anti.symbol = "̄" + particle.symbol;  // Add bar
        }

        // Flip quantum numbers
        anti.charge = -particle.charge;
        anti.baryon_number = -particle.baryon_number;
        anti.lepton_number = -particle.lepton_number;

        // Mass, spin, and lifetime remain the same
        // (CPT theorem guarantees this)

        return anti;
    }

    /**
     * @brief Check if particle is its own antiparticle
     *
     * Examples: photon, Z boson, neutral pion
     */
    static bool isSelfConjugate(const ParticleProperties& particle) {
        // Particle is self-conjugate if all additive quantum numbers are zero
        return (std::abs(particle.charge) < 1e-10) &&
               (particle.baryon_number == 0) &&
               (particle.lepton_number == 0);
    }

    /**
     * @brief Antiparticle symbol notation
     */
    static std::string antiparticleSymbol(const std::string& symbol) {
        // Add bar over symbol for antiparticles
        return "̄" + symbol;
    }
};

/**
 * @class AntiQuarks
 * @brief Antiquarks
 */
class AntiQuarks {
public:
    static ParticleProperties up_bar() {
        return Antiparticle::createAntiparticle(Quarks::up());
    }

    static ParticleProperties down_bar() {
        return Antiparticle::createAntiparticle(Quarks::down());
    }

    static ParticleProperties charm_bar() {
        return Antiparticle::createAntiparticle(Quarks::charm());
    }

    static ParticleProperties strange_bar() {
        return Antiparticle::createAntiparticle(Quarks::strange());
    }

    static ParticleProperties top_bar() {
        return Antiparticle::createAntiparticle(Quarks::top());
    }

    static ParticleProperties bottom_bar() {
        return Antiparticle::createAntiparticle(Quarks::bottom());
    }
};

/**
 * @class AntiLeptons
 * @brief Antileptons
 */
class AntiLeptons {
public:
    /**
     * @brief Positron (e⁺) - antielectron
     */
    static ParticleProperties positron() {
        auto anti_e = Antiparticle::createAntiparticle(Leptons::electron());
        anti_e.name = "Positron";
        anti_e.symbol = "e⁺";
        return anti_e;
    }

    static ParticleProperties electron_antineutrino() {
        return Antiparticle::createAntiparticle(Leptons::electron_neutrino());
    }

    static ParticleProperties antimuon() {
        auto anti_mu = Antiparticle::createAntiparticle(Leptons::muon());
        anti_mu.symbol = "μ⁺";
        return anti_mu;
    }

    static ParticleProperties muon_antineutrino() {
        return Antiparticle::createAntiparticle(Leptons::muon_neutrino());
    }

    static ParticleProperties antitau() {
        auto anti_tau = Antiparticle::createAntiparticle(Leptons::tau());
        anti_tau.symbol = "τ⁺";
        return anti_tau;
    }

    static ParticleProperties tau_antineutrino() {
        return Antiparticle::createAntiparticle(Leptons::tau_neutrino());
    }
};

/**
 * @class ChargeConjugation
 * @brief Charge conjugation symmetry (C)
 *
 * C operation: particle ↔ antiparticle
 */
class ChargeConjugation {
public:
    /**
     * @brief Apply charge conjugation operator
     *
     * C|particle⟩ = |antiparticle⟩
     */
    static ParticleProperties apply(const ParticleProperties& particle) {
        return Antiparticle::createAntiparticle(particle);
    }

    /**
     * @brief C-parity (charge conjugation eigenvalue)
     *
     * Only defined for self-conjugate particles
     * C|state⟩ = η_C|state⟩ where η_C = ±1
     *
     * Examples:
     * - π⁰: C = +1
     * - η: C = +1
     * - γ (photon): C = -1
     */
    static int cParity(const std::string& particle_name) {
        if (particle_name == "photon") return -1;
        if (particle_name == "neutral pion") return +1;
        if (particle_name == "eta meson") return +1;

        throw std::invalid_argument("C-parity only defined for self-conjugate states");
    }

    /**
     * @brief Check C symmetry conservation
     *
     * Electromagnetic and strong interactions conserve C
     * Weak interactions violate C maximally
     */
    static bool isConserved(const std::string& interaction) {
        if (interaction == "electromagnetic") return true;
        if (interaction == "strong") return true;
        if (interaction == "weak") return false;  // C violation

        return false;
    }
};

/**
 * @class CPSymmetry
 * @brief Combined charge conjugation and parity (CP)
 *
 * CP: particle → antiparticle + spatial inversion
 */
class CPSymmetry {
public:
    /**
     * @brief CP transformation
     *
     * CP|particle(x)⟩ = |antiparticle(-x)⟩
     */
    static std::string cpTransformation() {
        return "CP: (C × P) transforms particle at x → antiparticle at -x";
    }

    /**
     * @brief CP violation
     *
     * Discovered in 1964 (Cronin and Fitch)
     * Observed in K⁰ meson decay
     *
     * CP violation is source of matter-antimatter asymmetry
     */
    static double cpViolationParameter() {
        // ε parameter in K⁰ system
        return 2.228e-3;  // |ε| ≈ 0.00228
    }

    /**
     * @brief CKM matrix CP-violating phase
     *
     * δ_CKM ≈ 1.2 radians ≈ 69°
     */
    static double ckmPhase() {
        return 1.2;  // radians
    }

    /**
     * @brief Check if CP is conserved
     *
     * Weak interactions can violate CP
     */
    static bool isConserved(const std::string& interaction) {
        if (interaction == "electromagnetic") return true;
        if (interaction == "strong") return true;  // QCD conserves CP (strong CP problem)
        if (interaction == "weak") return false;   // Weak can violate CP

        return false;
    }
};

/**
 * @class PairProduction
 * @brief Particle-antiparticle pair creation
 *
 * γ → e⁺ + e⁻ (photon to electron-positron pair)
 */
class PairProduction {
public:
    /**
     * @brief Threshold energy for pair production
     *
     * E_threshold = 2mc² (minimum energy needed)
     *
     * @param particle_mass m (GeV/c²)
     * @param c Speed of light (m/s)
     * @return Threshold energy (GeV)
     */
    static double thresholdEnergy(double particle_mass, double c = 2.998e8) {
        return 2.0 * particle_mass;  // 2mc² in natural units (c=1)
    }

    /**
     * @brief Photon energy for e⁺e⁻ pair production
     *
     * E_γ ≥ 2m_e c² = 1.022 MeV
     */
    static double electronPairThreshold() {
        return 1.022e-3;  // 1.022 MeV in GeV
    }

    /**
     * @brief Check if energy sufficient for pair production
     */
    static bool canProducePair(double photon_energy, double particle_mass) {
        return photon_energy >= thresholdEnergy(particle_mass);
    }

    /**
     * @brief Pair production cross section (order of magnitude)
     *
     * σ ~ α² r_e² ~ 10⁻²⁵ cm² for e⁺e⁻
     *
     * where α ≈ 1/137 (fine structure constant)
     *       r_e ≈ 2.8 × 10⁻¹⁵ m (classical electron radius)
     */
    static double crossSection(double photon_energy, double particle_mass) {
        double threshold = thresholdEnergy(particle_mass);

        if (photon_energy < threshold) {
            return 0.0;  // Below threshold
        }

        // Simplified: actual depends on detailed QED calculation
        double alpha = 1.0/137.0;  // Fine structure constant
        double r_e = 2.818e-15;    // Classical electron radius (m)

        // σ ~ α² r_e² (order of magnitude)
        return alpha * alpha * r_e * r_e;  // m²
    }
};

/**
 * @class Annihilation
 * @brief Particle-antiparticle annihilation
 *
 * e⁺ + e⁻ → 2γ (electron-positron to photons)
 */
class Annihilation {
public:
    /**
     * @brief Energy released in annihilation
     *
     * E = 2mc² (rest mass energy)
     *
     * @param particle_mass m (GeV/c²)
     * @return Energy (GeV)
     */
    static double energyReleased(double particle_mass) {
        return 2.0 * particle_mass;  // 2mc²
    }

    /**
     * @brief Electron-positron annihilation energy
     *
     * E = 2m_e c² = 1.022 MeV
     */
    static double electronPositronEnergy() {
        return 1.022e-3;  // GeV
    }

    /**
     * @brief Photon energy from e⁺e⁻ annihilation at rest
     *
     * Each photon gets m_e c² = 0.511 MeV
     */
    static double photonEnergyAtRest() {
        return 0.511e-3;  // GeV per photon
    }

    /**
     * @brief Annihilation cross section
     *
     * For e⁺e⁻ → γγ at low energy:
     * σ = πr_e²(c/v) where v is relative velocity
     *
     * At threshold (v→0): σ diverges (1/v behavior)
     */
    static double crossSection(double relative_velocity) {
        if (relative_velocity < 1e-10) {
            return 1e10;  // Diverges at v→0 (return large value)
        }

        double r_e = 2.818e-15;    // Classical electron radius (m)
        double c = 2.998e8;        // m/s

        return M_PI * r_e * r_e * (c / relative_velocity);  // m²
    }

    /**
     * @brief Annihilation rate
     *
     * Γ = σ × v × n
     *
     * where n is antiparticle density
     */
    static double annihilationRate(double cross_section, double relative_velocity,
                                  double density) {
        return cross_section * relative_velocity * density;  // 1/s
    }
};

/**
 * @class MatterAntimatterAsymmetry
 * @brief Baryon asymmetry of the universe
 *
 * Why is the universe made of matter, not antimatter?
 */
class MatterAntimatterAsymmetry {
public:
    /**
     * @brief Baryon-to-photon ratio
     *
     * η = n_B / n_γ ≈ 6 × 10⁻¹⁰
     *
     * Tiny excess of matter over antimatter in early universe
     */
    static double baryonToPhotonRatio() {
        return 6.0e-10;
    }

    /**
     * @brief Sakharov conditions for baryogenesis
     *
     * Three conditions needed to generate baryon asymmetry:
     * 1. Baryon number violation
     * 2. C and CP violation
     * 3. Departure from thermal equilibrium
     */
    static std::vector<std::string> sakharovConditions() {
        return {
            "1. Baryon number violation (B violation)",
            "2. C and CP violation",
            "3. Departure from thermal equilibrium"
        };
    }

    /**
     * @brief Estimated matter-antimatter asymmetry
     *
     * For every ~10⁹ antimatter particles, there were
     * ~10⁹ + 1 matter particles in early universe
     *
     * After annihilation, 1 matter particle remains per 10⁹ photons
     */
    static double asymmetryParameter() {
        return 1.0e-9;  // Order of magnitude
    }

    /**
     * @brief CP violation strength needed
     *
     * Standard Model CP violation (CKM matrix) is too small
     * to explain observed asymmetry → need BSM physics
     */
    static std::string cpViolationStrength() {
        return "Standard Model CP violation insufficient; "
               "requires Beyond Standard Model (BSM) physics";
    }
};

} // namespace physics::advanced::qft

#endif // PHYSICS_ADVANCED_QFT_ANTIPARTICLES_HPP
