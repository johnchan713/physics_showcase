#ifndef PHYSICS_ADVANCED_QFT_PARTICLE_PHYSICS_HPP
#define PHYSICS_ADVANCED_QFT_PARTICLE_PHYSICS_HPP

#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <map>

/**
 * @file particle_physics.hpp
 * @brief Standard Model particles: quarks, leptons, fermions, bosons
 *
 * Implements:
 * - Fundamental fermions (quarks and leptons)
 * - Gauge bosons (photon, W, Z, gluons)
 * - Higgs boson
 * - Particle properties (mass, charge, spin, quantum numbers)
 * - Generations and families
 */

namespace physics::advanced::qft {

/**
 * @enum ParticleType
 * @brief Classification of fundamental particles
 */
enum class ParticleType {
    QUARK,
    LEPTON,
    GAUGE_BOSON,
    HIGGS_BOSON,
    UNKNOWN
};

/**
 * @enum SpinType
 * @brief Spin classification
 */
enum class SpinType {
    FERMION,  // Half-integer spin (1/2, 3/2, ...)
    BOSON     // Integer spin (0, 1, 2, ...)
};

/**
 * @struct ParticleProperties
 * @brief Properties of a fundamental particle
 */
struct ParticleProperties {
    std::string name;
    std::string symbol;
    ParticleType type;
    SpinType spin_type;

    double mass;          // GeV/c²
    double charge;        // Elementary charge units
    double spin;          // ℏ units
    int color_charge;     // 0 (colorless), 1 (colored)

    int generation;       // 1, 2, or 3 for fermions
    int baryon_number;    // 1/3 for quarks, 0 for leptons
    int lepton_number;    // 1 for leptons, 0 for quarks

    double lifetime;      // seconds (0 for stable)

    /**
     * @brief Check if particle is fermion
     */
    bool isFermion() const {
        return spin_type == SpinType::FERMION;
    }

    /**
     * @brief Check if particle is boson
     */
    bool isBoson() const {
        return spin_type == SpinType::BOSON;
    }

    /**
     * @brief Check if particle is stable
     */
    bool isStable() const {
        return lifetime == 0.0 || lifetime > 1e10;  // Effectively stable
    }
};

/**
 * @class Quarks
 * @brief Six flavors of quarks in three generations
 *
 * Generation 1: up (u), down (d)
 * Generation 2: charm (c), strange (s)
 * Generation 3: top (t), bottom (b)
 */
class Quarks {
public:
    /**
     * @brief Up quark properties
     *
     * Mass: ~2.2 MeV/c²
     * Charge: +2/3 e
     * Spin: 1/2
     */
    static ParticleProperties up() {
        return {
            "Up quark", "u", ParticleType::QUARK, SpinType::FERMION,
            0.0022,      // 2.2 MeV
            +2.0/3.0,    // +2/3 e
            0.5,         // spin 1/2
            1,           // colored
            1,           // generation 1
            1,           // baryon number 1/3 (stored as 1 for quarks)
            0,           // lepton number
            0.0          // stable (in hadrons)
        };
    }

    /**
     * @brief Down quark properties
     *
     * Mass: ~4.7 MeV/c²
     * Charge: -1/3 e
     */
    static ParticleProperties down() {
        return {
            "Down quark", "d", ParticleType::QUARK, SpinType::FERMION,
            0.0047,      // 4.7 MeV
            -1.0/3.0,    // -1/3 e
            0.5, 1, 1, 1, 0, 0.0
        };
    }

    /**
     * @brief Charm quark properties
     *
     * Mass: ~1.28 GeV/c²
     * Charge: +2/3 e
     */
    static ParticleProperties charm() {
        return {
            "Charm quark", "c", ParticleType::QUARK, SpinType::FERMION,
            1.28,        // 1.28 GeV
            +2.0/3.0,
            0.5, 1, 2, 1, 0, 0.0
        };
    }

    /**
     * @brief Strange quark properties
     *
     * Mass: ~96 MeV/c²
     * Charge: -1/3 e
     */
    static ParticleProperties strange() {
        return {
            "Strange quark", "s", ParticleType::QUARK, SpinType::FERMION,
            0.096,       // 96 MeV
            -1.0/3.0,
            0.5, 1, 2, 1, 0, 0.0
        };
    }

    /**
     * @brief Top quark properties
     *
     * Mass: ~173 GeV/c²
     * Charge: +2/3 e
     * Lifetime: ~5×10⁻²⁵ s (decays before hadronization)
     */
    static ParticleProperties top() {
        return {
            "Top quark", "t", ParticleType::QUARK, SpinType::FERMION,
            173.0,       // 173 GeV (heaviest elementary particle)
            +2.0/3.0,
            0.5, 1, 3, 1, 0, 5e-25  // Decays before hadronization
        };
    }

    /**
     * @brief Bottom quark properties
     *
     * Mass: ~4.18 GeV/c²
     * Charge: -1/3 e
     */
    static ParticleProperties bottom() {
        return {
            "Bottom quark", "b", ParticleType::QUARK, SpinType::FERMION,
            4.18,        // 4.18 GeV
            -1.0/3.0,
            0.5, 1, 3, 1, 0, 0.0
        };
    }

    /**
     * @brief Get all quark flavors
     */
    static std::vector<ParticleProperties> allFlavors() {
        return {up(), down(), charm(), strange(), top(), bottom()};
    }

    /**
     * @brief Get quarks by generation
     */
    static std::vector<ParticleProperties> generation(int gen) {
        switch(gen) {
            case 1: return {up(), down()};
            case 2: return {charm(), strange()};
            case 3: return {top(), bottom()};
            default: return {};
        }
    }
};

/**
 * @class Leptons
 * @brief Six leptons in three generations
 *
 * Charged leptons: electron (e), muon (μ), tau (τ)
 * Neutrinos: νₑ, νμ, ντ
 */
class Leptons {
public:
    /**
     * @brief Electron properties
     *
     * Mass: 0.511 MeV/c²
     * Charge: -1 e
     * Spin: 1/2
     * Stable
     */
    static ParticleProperties electron() {
        return {
            "Electron", "e⁻", ParticleType::LEPTON, SpinType::FERMION,
            0.000511,    // 0.511 MeV
            -1.0,        // -1 e
            0.5,         // spin 1/2
            0,           // colorless
            1,           // generation 1
            0,           // baryon number
            1,           // lepton number
            0.0          // stable
        };
    }

    /**
     * @brief Electron neutrino properties
     *
     * Mass: < 1 eV/c² (very small, not zero due to oscillations)
     * Charge: 0
     */
    static ParticleProperties electron_neutrino() {
        return {
            "Electron neutrino", "νₑ", ParticleType::LEPTON, SpinType::FERMION,
            1e-9,        // < 1 eV (approximate)
            0.0,         // neutral
            0.5, 0, 1, 0, 1, 0.0  // stable
        };
    }

    /**
     * @brief Muon properties
     *
     * Mass: 105.7 MeV/c²
     * Charge: -1 e
     * Lifetime: 2.2 μs
     */
    static ParticleProperties muon() {
        return {
            "Muon", "μ⁻", ParticleType::LEPTON, SpinType::FERMION,
            0.1057,      // 105.7 MeV
            -1.0,
            0.5, 0, 2, 0, 1, 2.2e-6  // τ = 2.2 μs
        };
    }

    /**
     * @brief Muon neutrino properties
     */
    static ParticleProperties muon_neutrino() {
        return {
            "Muon neutrino", "νμ", ParticleType::LEPTON, SpinType::FERMION,
            1e-9, 0.0, 0.5, 0, 2, 0, 1, 0.0
        };
    }

    /**
     * @brief Tau properties
     *
     * Mass: 1.777 GeV/c²
     * Charge: -1 e
     * Lifetime: 290 fs
     */
    static ParticleProperties tau() {
        return {
            "Tau", "τ⁻", ParticleType::LEPTON, SpinType::FERMION,
            1.777,       // 1.777 GeV
            -1.0,
            0.5, 0, 3, 0, 1, 2.9e-13  // τ = 290 fs
        };
    }

    /**
     * @brief Tau neutrino properties
     */
    static ParticleProperties tau_neutrino() {
        return {
            "Tau neutrino", "ντ", ParticleType::LEPTON, SpinType::FERMION,
            1e-9, 0.0, 0.5, 0, 3, 0, 1, 0.0
        };
    }

    /**
     * @brief Get all leptons
     */
    static std::vector<ParticleProperties> allLeptons() {
        return {
            electron(), electron_neutrino(),
            muon(), muon_neutrino(),
            tau(), tau_neutrino()
        };
    }

    /**
     * @brief Get charged leptons only
     */
    static std::vector<ParticleProperties> chargedLeptons() {
        return {electron(), muon(), tau()};
    }

    /**
     * @brief Get neutrinos only
     */
    static std::vector<ParticleProperties> neutrinos() {
        return {electron_neutrino(), muon_neutrino(), tau_neutrino()};
    }
};

/**
 * @class GaugeBosons
 * @brief Force carriers of the Standard Model
 *
 * Photon (γ): Electromagnetic force
 * W±, Z⁰: Weak force
 * Gluons (g): Strong force
 */
class GaugeBosons {
public:
    /**
     * @brief Photon properties
     *
     * Mass: 0 (exact)
     * Charge: 0
     * Spin: 1
     * Mediates: Electromagnetic force
     */
    static ParticleProperties photon() {
        return {
            "Photon", "γ", ParticleType::GAUGE_BOSON, SpinType::BOSON,
            0.0,         // massless
            0.0,         // neutral
            1.0,         // spin 1
            0,           // colorless
            0,           // no generation
            0, 0, 0.0    // stable
        };
    }

    /**
     * @brief W+ boson properties
     *
     * Mass: 80.4 GeV/c²
     * Charge: +1 e
     * Spin: 1
     */
    static ParticleProperties W_plus() {
        return {
            "W boson", "W⁺", ParticleType::GAUGE_BOSON, SpinType::BOSON,
            80.4,        // 80.4 GeV
            +1.0,        // +1 e
            1.0, 0, 0, 0, 0, 3e-25  // τ ~ 3×10⁻²⁵ s
        };
    }

    /**
     * @brief W- boson properties
     */
    static ParticleProperties W_minus() {
        return {
            "W boson", "W⁻", ParticleType::GAUGE_BOSON, SpinType::BOSON,
            80.4, -1.0, 1.0, 0, 0, 0, 0, 3e-25
        };
    }

    /**
     * @brief Z boson properties
     *
     * Mass: 91.2 GeV/c²
     * Charge: 0
     * Spin: 1
     */
    static ParticleProperties Z_boson() {
        return {
            "Z boson", "Z⁰", ParticleType::GAUGE_BOSON, SpinType::BOSON,
            91.2,        // 91.2 GeV
            0.0,         // neutral
            1.0, 0, 0, 0, 0, 3e-25  // τ ~ 3×10⁻²⁵ s
        };
    }

    /**
     * @brief Gluon properties
     *
     * Mass: 0 (exact, but confined)
     * Charge: 0
     * Spin: 1
     * Color charge: carries color (8 gluons for SU(3))
     */
    static ParticleProperties gluon() {
        return {
            "Gluon", "g", ParticleType::GAUGE_BOSON, SpinType::BOSON,
            0.0,         // massless (confined)
            0.0,         // electrically neutral
            1.0,         // spin 1
            1,           // carries color charge
            0, 0, 0, 0.0 // confined (never free)
        };
    }

    /**
     * @brief Get all gauge bosons
     */
    static std::vector<ParticleProperties> allBosons() {
        return {photon(), W_plus(), W_minus(), Z_boson(), gluon()};
    }
};

/**
 * @class HiggsBoson
 * @brief Higgs boson - gives mass to particles
 */
class HiggsBoson {
public:
    /**
     * @brief Higgs boson properties
     *
     * Mass: 125.1 GeV/c²
     * Charge: 0
     * Spin: 0 (scalar)
     * Lifetime: ~1.6×10⁻²² s
     */
    static ParticleProperties higgs() {
        return {
            "Higgs boson", "H⁰", ParticleType::HIGGS_BOSON, SpinType::BOSON,
            125.1,       // 125.1 GeV (discovered 2012)
            0.0,         // neutral
            0.0,         // spin 0 (only fundamental scalar)
            0,           // colorless
            0, 0, 0, 1.6e-22  // τ ~ 1.6×10⁻²² s
        };
    }
};

/**
 * @class FermionBosonsClassification
 * @brief Classification by spin-statistics
 */
class FermionBosonsClassification {
public:
    /**
     * @brief Get all fermions (quarks + leptons)
     */
    static std::vector<ParticleProperties> allFermions() {
        auto quarks = Quarks::allFlavors();
        auto leptons = Leptons::allLeptons();

        std::vector<ParticleProperties> fermions;
        fermions.insert(fermions.end(), quarks.begin(), quarks.end());
        fermions.insert(fermions.end(), leptons.begin(), leptons.end());

        return fermions;
    }

    /**
     * @brief Get all bosons (gauge + Higgs)
     */
    static std::vector<ParticleProperties> allBosons() {
        auto gauge = GaugeBosons::allBosons();
        auto higgs_vec = std::vector<ParticleProperties>{HiggsBoson::higgs()};

        std::vector<ParticleProperties> bosons;
        bosons.insert(bosons.end(), gauge.begin(), gauge.end());
        bosons.insert(bosons.end(), higgs_vec.begin(), higgs_vec.end());

        return bosons;
    }

    /**
     * @brief Count particles by type
     */
    static std::map<std::string, int> particleCount() {
        return {
            {"quarks", 6},
            {"leptons", 6},
            {"gauge_bosons", 5},  // γ, W+, W-, Z, g (8 gluons counted as 1)
            {"higgs", 1},
            {"total_fermions", 12},
            {"total_bosons", 6}
        };
    }
};

/**
 * @class Generations
 * @brief Three generations of matter
 */
class Generations {
public:
    /**
     * @brief Get all particles in generation n (1, 2, or 3)
     */
    static std::vector<ParticleProperties> generation(int n) {
        std::vector<ParticleProperties> gen;

        // Add quarks
        auto quarks = Quarks::generation(n);
        gen.insert(gen.end(), quarks.begin(), quarks.end());

        // Add leptons
        switch(n) {
            case 1:
                gen.push_back(Leptons::electron());
                gen.push_back(Leptons::electron_neutrino());
                break;
            case 2:
                gen.push_back(Leptons::muon());
                gen.push_back(Leptons::muon_neutrino());
                break;
            case 3:
                gen.push_back(Leptons::tau());
                gen.push_back(Leptons::tau_neutrino());
                break;
        }

        return gen;
    }

    /**
     * @brief Mass hierarchy across generations
     *
     * Each generation is heavier than the previous
     */
    static std::array<double, 3> averageMassPerGeneration() {
        std::array<double, 3> masses = {0, 0, 0};

        for (int i = 1; i <= 3; ++i) {
            auto gen = generation(i);
            double sum = 0.0;
            for (const auto& p : gen) {
                sum += p.mass;
            }
            masses[i-1] = sum / gen.size();
        }

        return masses;
    }
};

/**
 * @class QuantumNumbers
 * @brief Conservation laws and quantum numbers
 */
class QuantumNumbers {
public:
    /**
     * @brief Calculate total electric charge
     */
    static double totalCharge(const std::vector<ParticleProperties>& particles) {
        double Q = 0.0;
        for (const auto& p : particles) {
            Q += p.charge;
        }
        return Q;
    }

    /**
     * @brief Calculate total baryon number
     */
    static double totalBaryonNumber(const std::vector<ParticleProperties>& particles) {
        double B = 0.0;
        for (const auto& p : particles) {
            if (p.type == ParticleType::QUARK) {
                B += 1.0/3.0;  // Each quark contributes 1/3
            }
        }
        return B;
    }

    /**
     * @brief Calculate total lepton number
     */
    static int totalLeptonNumber(const std::vector<ParticleProperties>& particles) {
        int L = 0;
        for (const auto& p : particles) {
            L += p.lepton_number;
        }
        return L;
    }

    /**
     * @brief Check if quantum numbers are conserved
     */
    static bool checkConservation(
        const std::vector<ParticleProperties>& initial,
        const std::vector<ParticleProperties>& final,
        double tolerance = 1e-6) {

        // Charge conservation
        double Q_i = totalCharge(initial);
        double Q_f = totalCharge(final);
        if (std::abs(Q_i - Q_f) > tolerance) return false;

        // Baryon number conservation
        double B_i = totalBaryonNumber(initial);
        double B_f = totalBaryonNumber(final);
        if (std::abs(B_i - B_f) > tolerance) return false;

        // Lepton number conservation
        int L_i = totalLeptonNumber(initial);
        int L_f = totalLeptonNumber(final);
        if (L_i != L_f) return false;

        return true;
    }
};

} // namespace physics::advanced::qft

#endif // PHYSICS_ADVANCED_QFT_PARTICLE_PHYSICS_HPP
