#ifndef PHYSICS_ADVANCED_QFT_SPIN_STATISTICS_HPP
#define PHYSICS_ADVANCED_QFT_SPIN_STATISTICS_HPP

#include <cmath>
#include <stdexcept>
#include <vector>

/**
 * @file spin_statistics.hpp
 * @brief Spin-statistics theorem and Pauli exclusion principle
 *
 * Implements:
 * - Spin-statistics connection
 * - Fermi-Dirac statistics (fermions)
 * - Bose-Einstein statistics (bosons)
 * - Pauli exclusion principle
 * - Symmetry/antisymmetry of wavefunctions
 */

namespace physics::advanced::qft {

/**
 * @class SpinStatisticsTheorem
 * @brief Fundamental connection between spin and statistics
 *
 * Theorem: Particles with half-integer spin (fermions) obey Fermi-Dirac
 * statistics and have antisymmetric wavefunctions. Particles with integer
 * spin (bosons) obey Bose-Einstein statistics and have symmetric wavefunctions.
 *
 * This is a deep result from relativistic quantum field theory.
 */
class SpinStatisticsTheorem {
public:
    /**
     * @brief Check if spin is fermionic (half-integer)
     *
     * Fermions: spin = n/2 where n is odd (1/2, 3/2, 5/2, ...)
     */
    static bool isFermionicSpin(double spin) {
        // Check if 2*spin is odd integer
        double two_spin = 2.0 * spin;
        double rounded = std::round(two_spin);

        if (std::abs(two_spin - rounded) > 1e-10) {
            return false;  // Not integer/half-integer
        }

        int n = static_cast<int>(rounded);
        return (n % 2) == 1;  // Odd means half-integer
    }

    /**
     * @brief Check if spin is bosonic (integer)
     *
     * Bosons: spin = n where n is integer (0, 1, 2, ...)
     */
    static bool isBosonicSpin(double spin) {
        double rounded = std::round(spin);
        return std::abs(spin - rounded) < 1e-10;
    }

    /**
     * @brief Get statistics type from spin
     */
    static std::string statisticsFromSpin(double spin) {
        if (isFermionicSpin(spin)) {
            return "Fermi-Dirac";
        } else if (isBosonicSpin(spin)) {
            return "Bose-Einstein";
        } else {
            return "Unknown";
        }
    }

    /**
     * @brief Wavefunction symmetry from spin
     *
     * Fermions: antisymmetric (ψ(1,2) = -ψ(2,1))
     * Bosons: symmetric (ψ(1,2) = +ψ(2,1))
     */
    static std::string wavefunctionSymmetry(double spin) {
        if (isFermionicSpin(spin)) {
            return "Antisymmetric";
        } else if (isBosonicSpin(spin)) {
            return "Symmetric";
        } else {
            return "Unknown";
        }
    }

    /**
     * @brief Exchange phase factor
     *
     * Fermions: -1
     * Bosons: +1
     */
    static int exchangePhase(double spin) {
        if (isFermionicSpin(spin)) {
            return -1;
        } else if (isBosonicSpin(spin)) {
            return +1;
        } else {
            throw std::invalid_argument("Invalid spin for spin-statistics theorem");
        }
    }
};

/**
 * @class FermiDiracStatistics
 * @brief Statistics for fermions (half-integer spin)
 *
 * Obey Pauli exclusion principle: no two identical fermions
 * can occupy the same quantum state.
 */
class FermiDiracStatistics {
public:
    /**
     * @brief Fermi-Dirac distribution
     *
     * n(E) = 1 / [exp((E - μ)/(kT)) + 1]
     *
     * @param energy E (energy level)
     * @param chemical_potential μ (Fermi energy at T=0)
     * @param temperature T (Kelvin)
     * @param k_boltzmann k_B (J/K)
     * @return Occupation number (0 to 1)
     */
    static double distribution(double energy, double chemical_potential,
                              double temperature, double k_boltzmann = 1.380649e-23) {
        if (temperature <= 0.0) {
            // T = 0 limit: step function at Fermi energy
            return (energy < chemical_potential) ? 1.0 : 0.0;
        }

        double exponent = (energy - chemical_potential) / (k_boltzmann * temperature);

        // Avoid overflow for large exponents
        if (exponent > 50.0) return 0.0;
        if (exponent < -50.0) return 1.0;

        return 1.0 / (std::exp(exponent) + 1.0);
    }

    /**
     * @brief Fermi energy (chemical potential at T=0)
     *
     * E_F = (ℏ²/2m)(3π²n)^(2/3)
     *
     * @param number_density n (particles/m³)
     * @param mass m (kg)
     * @param hbar ℏ (J·s)
     * @return Fermi energy (J)
     */
    static double fermiEnergy(double number_density, double mass,
                             double hbar = 1.054571817e-34) {
        double prefactor = (hbar * hbar) / (2.0 * mass);
        double density_term = std::pow(3.0 * M_PI * M_PI * number_density, 2.0/3.0);
        return prefactor * density_term;
    }

    /**
     * @brief Fermi temperature
     *
     * T_F = E_F / k_B
     */
    static double fermiTemperature(double fermi_energy,
                                  double k_boltzmann = 1.380649e-23) {
        return fermi_energy / k_boltzmann;
    }

    /**
     * @brief Pauli exclusion principle check
     *
     * Returns true if state is available (occupation < 1)
     */
    static bool stateAvailable(double energy, double chemical_potential,
                              double temperature, double k_boltzmann = 1.380649e-23) {
        double occupation = distribution(energy, chemical_potential, temperature, k_boltzmann);
        return occupation < 1.0;  // Can always add one more (probabilistically)
    }

    /**
     * @brief Chemical potential from density (for free fermions)
     *
     * At T=0: μ = E_F
     * At T>0: Requires numerical solution, use approximation
     */
    static double chemicalPotential(double number_density, double mass,
                                   double temperature, double hbar = 1.054571817e-34,
                                   double k_boltzmann = 1.380649e-23) {
        double E_F = fermiEnergy(number_density, mass, hbar);

        if (temperature <= 0.0) {
            return E_F;
        }

        // Low temperature approximation: μ ≈ E_F[1 - (π²/12)(kT/E_F)²]
        double T_F = E_F / k_boltzmann;
        double ratio = temperature / T_F;
        return E_F * (1.0 - (M_PI * M_PI / 12.0) * ratio * ratio);
    }
};

/**
 * @class BoseEinsteinStatistics
 * @brief Statistics for bosons (integer spin)
 *
 * Multiple bosons can occupy the same quantum state.
 * Leads to phenomena like Bose-Einstein condensation.
 */
class BoseEinsteinStatistics {
public:
    /**
     * @brief Bose-Einstein distribution
     *
     * n(E) = 1 / [exp((E - μ)/(kT)) - 1]
     *
     * Note: μ < E₀ (chemical potential below ground state)
     *
     * @param energy E (energy level)
     * @param chemical_potential μ
     * @param temperature T (Kelvin)
     * @param k_boltzmann k_B (J/K)
     * @return Occupation number (0 to ∞)
     */
    static double distribution(double energy, double chemical_potential,
                              double temperature, double k_boltzmann = 1.380649e-23) {
        if (temperature <= 0.0) {
            throw std::invalid_argument("Temperature must be positive for Bose-Einstein");
        }

        if (chemical_potential >= energy) {
            throw std::invalid_argument("Chemical potential must be below energy level");
        }

        double exponent = (energy - chemical_potential) / (k_boltzmann * temperature);

        // Avoid overflow
        if (exponent > 50.0) return 0.0;

        return 1.0 / (std::exp(exponent) - 1.0);
    }

    /**
     * @brief Critical temperature for Bose-Einstein condensation
     *
     * T_c = (2πℏ²/mk_B)(n/ζ(3/2))^(2/3)
     *
     * where ζ(3/2) ≈ 2.612 is Riemann zeta function
     *
     * @param number_density n (particles/m³)
     * @param mass m (kg)
     * @return Critical temperature (K)
     */
    static double criticalTemperature(double number_density, double mass,
                                     double hbar = 1.054571817e-34,
                                     double k_boltzmann = 1.380649e-23) {
        double zeta_3_2 = 2.612;  // ζ(3/2)

        double prefactor = (2.0 * M_PI * hbar * hbar) / (mass * k_boltzmann);
        double density_term = std::pow(number_density / zeta_3_2, 2.0/3.0);

        return prefactor * density_term;
    }

    /**
     * @brief Condensate fraction below T_c
     *
     * N₀/N = 1 - (T/T_c)^(3/2)
     *
     * @param temperature T (K)
     * @param critical_temperature T_c (K)
     * @return Fraction in ground state
     */
    static double condensateFraction(double temperature, double critical_temperature) {
        if (temperature >= critical_temperature) {
            return 0.0;  // No condensate above T_c
        }

        double ratio = temperature / critical_temperature;
        return 1.0 - std::pow(ratio, 1.5);
    }

    /**
     * @brief Check if Bose-Einstein condensation occurs
     */
    static bool hasBoseCondensation(double temperature, double critical_temperature) {
        return temperature < critical_temperature;
    }

    /**
     * @brief Photon distribution (massless bosons, μ=0)
     *
     * Planck's blackbody radiation law
     *
     * n(ω) = 1 / [exp(ℏω/(kT)) - 1]
     */
    static double photonDistribution(double frequency, double temperature,
                                    double hbar = 1.054571817e-34,
                                    double k_boltzmann = 1.380649e-23) {
        double energy = hbar * frequency;
        double exponent = energy / (k_boltzmann * temperature);

        if (exponent > 50.0) return 0.0;

        return 1.0 / (std::exp(exponent) - 1.0);
    }
};

/**
 * @class PauliExclusionPrinciple
 * @brief Pauli exclusion principle for fermions
 */
class PauliExclusionPrinciple {
public:
    /**
     * @brief Maximum fermions per quantum state
     *
     * Exactly 1 fermion per quantum state (including spin)
     */
    static constexpr int maxOccupation() {
        return 1;
    }

    /**
     * @brief Check if adding fermion violates Pauli exclusion
     *
     * @param current_occupation Current number in state
     * @return true if can add fermion
     */
    static bool canAddFermion(int current_occupation) {
        return current_occupation < maxOccupation();
    }

    /**
     * @brief Degeneracy with spin
     *
     * For spin-1/2: 2 states (spin up, spin down)
     * For spin-S: (2S+1) states
     *
     * So maximum fermions per orbital = 2S+1
     */
    static int degeneracyWithSpin(double spin) {
        return static_cast<int>(2.0 * spin + 1.0);
    }

    /**
     * @brief Maximum fermions per orbital (accounting for spin)
     */
    static int maxFermionsPerOrbital(double spin) {
        return degeneracyWithSpin(spin);
    }

    /**
     * @brief Examples of Pauli exclusion
     */
    static std::string examples() {
        return
            "1. Atomic structure: Electrons fill shells (n, l, m, s)\n"
            "2. Neutron stars: Degeneracy pressure supports against collapse\n"
            "3. White dwarfs: Electron degeneracy pressure\n"
            "4. Chemistry: Periodic table structure\n"
            "5. Stability of matter: Why atoms have size";
    }
};

/**
 * @class WavefunctionSymmetry
 * @brief Symmetry properties of many-particle wavefunctions
 */
class WavefunctionSymmetry {
public:
    /**
     * @brief Symmetrization factor for identical bosons
     *
     * ψ_sym(1,2) = (1/√2)[ψ_a(1)ψ_b(2) + ψ_a(2)ψ_b(1)]
     *
     * @return Normalization factor 1/√N!
     */
    static double bosonNormalization(int num_particles) {
        double factorial = 1.0;
        for (int i = 2; i <= num_particles; ++i) {
            factorial *= i;
        }
        return 1.0 / std::sqrt(factorial);
    }

    /**
     * @brief Antisymmetrization factor for identical fermions
     *
     * ψ_antisym(1,2) = (1/√2)[ψ_a(1)ψ_b(2) - ψ_a(2)ψ_b(1)]
     *
     * This is a Slater determinant
     *
     * @return Normalization factor 1/√N!
     */
    static double fermionNormalization(int num_particles) {
        return bosonNormalization(num_particles);  // Same normalization
    }

    /**
     * @brief Exchange parity
     *
     * Fermions: (-1)^P where P is number of exchanges
     * Bosons: (+1)^P = +1 always
     */
    static int exchangeParity(int num_exchanges, bool is_fermion) {
        if (is_fermion) {
            return (num_exchanges % 2 == 0) ? +1 : -1;
        } else {
            return +1;  // Bosons always symmetric
        }
    }

    /**
     * @brief Check wavefunction symmetry
     */
    static std::string checkSymmetry(double spin) {
        if (SpinStatisticsTheorem::isFermionicSpin(spin)) {
            return "Must be antisymmetric (fermion)";
        } else if (SpinStatisticsTheorem::isBosonicSpin(spin)) {
            return "Must be symmetric (boson)";
        } else {
            return "Unknown symmetry";
        }
    }
};

} // namespace physics::advanced::qft

#endif // PHYSICS_ADVANCED_QFT_SPIN_STATISTICS_HPP
