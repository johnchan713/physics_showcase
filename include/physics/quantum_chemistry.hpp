#ifndef PHYSICS_QUANTUM_CHEMISTRY_HPP
#define PHYSICS_QUANTUM_CHEMISTRY_HPP

#include <vector>
#include <complex>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>
#include <map>

namespace physics {
namespace quantum_chemistry {

using Complex = std::complex<double>;

/**
 * @brief Physical constants for quantum chemistry
 */
namespace constants {
    constexpr double h = 6.62607015e-34;        // Planck constant (J·s)
    constexpr double hbar = 1.054571817e-34;    // Reduced Planck constant
    constexpr double c = 299792458.0;           // Speed of light (m/s)
    constexpr double e = 1.602176634e-19;       // Elementary charge (C)
    constexpr double m_e = 9.1093837015e-31;    // Electron mass (kg)
    constexpr double a_0 = 5.29177210903e-11;   // Bohr radius (m)
    constexpr double E_h = 4.359744722e-18;     // Hartree energy (J) = 27.211 eV
    constexpr double k_e = 8.9875517923e9;      // Coulomb constant (N·m²/C²)
}

/**
 * @brief Atomic and Molecular Wave Functions
 *
 * Multi-electron wave functions and antisymmetrization
 */
class AtomicWaveFunctions {
public:
    /**
     * @brief Product wave function (not antisymmetrized)
     *
     * ψ(1,2,...,N) = φ₁(1)φ₂(2)...φₙ(N)
     */
    static double product_wavefunction(
        const std::vector<std::function<double(double)>>& orbitals,
        const std::vector<double>& positions) {

        if (orbitals.size() != positions.size()) {
            throw std::invalid_argument("Orbital and position count mismatch");
        }

        double product = 1.0;
        for (size_t i = 0; i < orbitals.size(); ++i) {
            product *= orbitals[i](positions[i]);
        }

        return product;
    }

    /**
     * @brief Slater determinant for N electrons
     *
     * Ensures antisymmetry: ψ(1,2) = -ψ(2,1)
     *
     * For 2 electrons:
     * ψ = (1/√2)[φ₁(1)φ₂(2) - φ₁(2)φ₂(1)]
     */
    static double slater_determinant_2electron(
        const std::function<double(double)>& phi1,
        const std::function<double(double)>& phi2,
        double r1,
        double r2) {

        return (phi1(r1) * phi2(r2) - phi1(r2) * phi2(r1)) / std::sqrt(2.0);
    }

    /**
     * @brief Normalization integral
     *
     * ∫|ψ|² dτ = 1
     */
    static double normalization_integral(
        const std::function<double(double)>& psi,
        double r_min,
        double r_max,
        int n_points = 1000) {

        double dr = (r_max - r_min) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double r = r_min + i * dr;
            double psi_val = psi(r);
            // Include r² for spherical coordinates
            integral += psi_val * psi_val * r * r * dr;
        }

        return integral * 4.0 * M_PI;
    }

    /**
     * @brief Spin-spatial wave function
     *
     * ψ(r,s) = ψ_spatial(r) × χ_spin(s)
     */
    struct SpinSpatial {
        enum class Spin { UP, DOWN };

        static std::string spin_label(Spin s) {
            return s == Spin::UP ? "α" : "β";
        }
    };

    /**
     * @brief Exchange symmetry verification
     *
     * For fermions: ψ(1,2) = -ψ(2,1)
     * For bosons: ψ(1,2) = ψ(2,1)
     */
    static bool verify_antisymmetry(
        double psi_12,
        double psi_21,
        double tol = 1e-10) {

        return std::abs(psi_12 + psi_21) < tol;
    }
};

/**
 * @brief The Hartree-Fock Method
 *
 * Self-consistent field theory for many-electron atoms
 */
class HartreeFockMethod {
public:
    /**
     * @brief Fock operator
     *
     * F = h + Σⱼ(2Jⱼ - Kⱼ)
     *
     * h: one-electron operator
     * Jⱼ: Coulomb operator
     * Kⱼ: Exchange operator
     */
    struct FockOperator {
        double one_electron_energy;
        double coulomb_energy;
        double exchange_energy;

        double total() const {
            return one_electron_energy + 2.0 * coulomb_energy - exchange_energy;
        }
    };

    /**
     * @brief Coulomb integral
     *
     * Jᵢⱼ = ∫∫ φᵢ*(1)φᵢ(1) (e²/r₁₂) φⱼ*(2)φⱼ(2) dτ₁dτ₂
     */
    static double coulomb_integral(
        const std::function<double(double)>& phi_i,
        const std::function<double(double)>& phi_j,
        double r_min,
        double r_max,
        int n_points = 100) {

        double dr = (r_max - r_min) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            for (int j = 0; j < n_points; ++j) {
                double r1 = r_min + i * dr;
                double r2 = r_min + j * dr;

                double r12 = std::abs(r1 - r2);
                if (r12 < 1e-10) r12 = 1e-10;  // Avoid singularity

                double rho1 = phi_i(r1) * phi_i(r1) * r1 * r1;
                double rho2 = phi_j(r2) * phi_j(r2) * r2 * r2;

                integral += rho1 * rho2 * (constants::e * constants::e / r12) * dr * dr;
            }
        }

        return integral * 16.0 * M_PI * M_PI;
    }

    /**
     * @brief Exchange integral
     *
     * Kᵢⱼ = ∫∫ φᵢ*(1)φⱼ(1) (e²/r₁₂) φⱼ*(2)φᵢ(2) dτ₁dτ₂
     */
    static double exchange_integral(
        const std::function<double(double)>& phi_i,
        const std::function<double(double)>& phi_j,
        double r_min,
        double r_max,
        int n_points = 100) {

        // Similar structure to Coulomb integral but with crossed orbitals
        double dr = (r_max - r_min) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            for (int j = 0; j < n_points; ++j) {
                double r1 = r_min + i * dr;
                double r2 = r_min + j * dr;

                double r12 = std::abs(r1 - r2);
                if (r12 < 1e-10) r12 = 1e-10;

                double cross = phi_i(r1) * phi_j(r1) * phi_j(r2) * phi_i(r2);

                integral += cross * (constants::e * constants::e / r12) *
                           r1 * r1 * r2 * r2 * dr * dr;
            }
        }

        return integral * 16.0 * M_PI * M_PI;
    }

    /**
     * @brief Hartree-Fock energy
     *
     * E_HF = Σᵢ hᵢᵢ + ½ΣᵢΣⱼ(2Jᵢⱼ - Kᵢⱼ)
     */
    static double hartree_fock_energy(
        const std::vector<double>& one_electron_energies,
        const std::vector<std::vector<double>>& coulomb_matrix,
        const std::vector<std::vector<double>>& exchange_matrix) {

        int n = one_electron_energies.size();
        double energy = 0.0;

        // One-electron terms
        for (int i = 0; i < n; ++i) {
            energy += one_electron_energies[i];
        }

        // Two-electron terms
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                energy += 0.5 * (2.0 * coulomb_matrix[i][j] - exchange_matrix[i][j]);
            }
        }

        return energy;
    }

    /**
     * @brief Self-consistent field iteration
     *
     * Iterate until orbital energies converge
     */
    static bool check_scf_convergence(
        const std::vector<double>& energies_old,
        const std::vector<double>& energies_new,
        double tolerance = 1e-6) {

        for (size_t i = 0; i < energies_old.size(); ++i) {
            if (std::abs(energies_new[i] - energies_old[i]) > tolerance) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Koopmans' theorem
     *
     * Ionization energy ≈ -εᵢ (orbital energy)
     */
    static double koopmans_ionization_energy(double orbital_energy) {
        return -orbital_energy;
    }
};

/**
 * @brief Slater Orbitals
 *
 * Slater-type orbitals (STOs) for atomic calculations
 */
class SlaterOrbitals {
public:
    /**
     * @brief Slater-type orbital (STO)
     *
     * φₙₗₘ(r,θ,φ) = N rⁿ⁻¹ e^(-ζr) Yₗₘ(θ,φ)
     *
     * ζ: Slater exponent (screening parameter)
     */
    static double slater_radial(
        int n,
        double zeta,
        double r) {

        // Normalization constant
        double factorial_2n = 1.0;
        for (int i = 1; i <= 2*n; ++i) {
            factorial_2n *= i;
        }

        double norm = std::pow(2.0 * zeta, n + 0.5) / std::sqrt(factorial_2n);

        return norm * std::pow(r, n - 1) * std::exp(-zeta * r);
    }

    /**
     * @brief Slater's rules for effective nuclear charge
     *
     * Z_eff = Z - S
     *
     * S: screening constant from Slater's rules
     */
    static double slater_screening_constant(
        int Z,
        int n,
        int l,
        const std::vector<std::pair<int, int>>& electron_config) {

        double S = 0.0;

        // Electrons in same (n,l) group: contribute 0.35 each
        // Electrons in (n-1) shell: contribute 0.85 each
        // Electrons in (n-2) or lower: contribute 1.00 each

        for (const auto& [n_i, count] : electron_config) {
            if (n_i == n) {
                S += (count - 1) * 0.35;  // Same shell (excluding the electron itself)
            } else if (n_i == n - 1) {
                S += count * 0.85;
            } else if (n_i < n - 1) {
                S += count * 1.00;
            }
        }

        return S;
    }

    /**
     * @brief Effective nuclear charge
     *
     * Z_eff = Z - S
     */
    static double effective_nuclear_charge(int Z, double screening) {
        return Z - screening;
    }

    /**
     * @brief Slater exponent
     *
     * ζ = (Z - S) / n*
     *
     * n*: effective principal quantum number
     */
    static double slater_exponent(int n, double Z_eff) {
        // Effective principal quantum number
        double n_star = n;
        if (n == 1) n_star = 1.0;
        else if (n == 2) n_star = 2.0;
        else if (n == 3) n_star = 3.0;
        else if (n >= 4) n_star = 3.7;

        return Z_eff / n_star;
    }

    /**
     * @brief STO overlap integral
     *
     * Sᵢⱼ = ∫ φᵢ* φⱼ dτ
     */
    static double overlap_integral(
        int n1,
        double zeta1,
        int n2,
        double zeta2,
        double r_min,
        double r_max,
        int n_points = 1000) {

        double dr = (r_max - r_min) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double r = r_min + i * dr;
            double phi1 = slater_radial(n1, zeta1, r);
            double phi2 = slater_radial(n2, zeta2, r);

            integral += phi1 * phi2 * r * r * dr;
        }

        return integral * 4.0 * M_PI;
    }
};

/**
 * @brief Multiplet Theory
 *
 * Term symbols and L-S coupling for multi-electron atoms
 */
class MultipletTheory {
public:
    /**
     * @brief Term symbol
     *
     * ²ˢ⁺¹Lⱼ
     *
     * S: total spin
     * L: total orbital angular momentum
     * J: total angular momentum
     */
    struct TermSymbol {
        int S;  // Total spin (actually 2S for integer storage)
        int L;  // Total orbital angular momentum
        int J;  // Total angular momentum (actually 2J)

        std::string to_string() const {
            std::string L_symbol;
            if (L == 0) L_symbol = "S";
            else if (L == 1) L_symbol = "P";
            else if (L == 2) L_symbol = "D";
            else if (L == 3) L_symbol = "F";
            else if (L == 4) L_symbol = "G";
            else L_symbol = std::to_string(L);

            int multiplicity = S + 1;  // 2S+1

            return std::to_string(multiplicity) + L_symbol +
                   std::to_string(J) + "/" + "2";
        }
    };

    /**
     * @brief L-S coupling (Russell-Saunders)
     *
     * L = Σlᵢ, S = Σsᵢ, J = L + S
     */
    static std::vector<TermSymbol> ls_coupling(
        const std::vector<int>& l_values,
        const std::vector<int>& s_values) {

        std::vector<TermSymbol> terms;

        // Total L
        int L_max = std::accumulate(l_values.begin(), l_values.end(), 0);

        // Total S (stored as 2S)
        int total_2S = std::accumulate(s_values.begin(), s_values.end(), 0);

        // J ranges from |L-S| to L+S
        int S = total_2S / 2;
        for (int L = 0; L <= L_max; ++L) {
            int J_min = std::abs(L - S);
            int J_max = L + S;

            for (int J = J_min; J <= J_max; ++J) {
                terms.push_back({total_2S, L, 2*J});
            }
        }

        return terms;
    }

    /**
     * @brief Hund's rules
     *
     * 1. Maximize total spin S
     * 2. Maximize total orbital angular momentum L
     * 3. J = |L-S| if less than half-filled, J = L+S if more than half-filled
     */
    static TermSymbol hunds_rules_ground_state(
        int n_electrons,
        int subshell_capacity) {

        // Simplified for demonstration
        // Real implementation would need detailed electron configuration

        int S = 0;  // Will store 2S
        int L = 0;

        // Maximize spin first
        int half_fill = subshell_capacity / 2;
        if (n_electrons <= half_fill) {
            S = n_electrons;  // All spins parallel (2S = n for s=1/2)
        } else {
            S = subshell_capacity - n_electrons;
        }

        // Determine J
        int J;
        if (n_electrons <= half_fill) {
            J = std::abs(L - S/2);  // Less than half-filled
        } else {
            J = L + S/2;  // More than half-filled
        }

        return {S, L, 2*J};
    }

    /**
     * @brief Spectroscopic notation
     */
    static std::string spectroscopic_notation(int L) {
        static const char* symbols[] = {"S", "P", "D", "F", "G", "H", "I", "K"};
        if (L >= 0 && L < 8) {
            return symbols[L];
        }
        return std::to_string(L);
    }

    /**
     * @brief Fine structure splitting
     *
     * ΔE ∝ [J(J+1) - L(L+1) - S(S+1)]
     */
    static double fine_structure_energy(int L, int S, int J, double coupling_constant) {
        return coupling_constant *
               (J * (J + 1) - L * (L + 1) - S * (S + 1)) / 2.0;
    }
};

/**
 * @brief The Born-Oppenheimer Approximation
 *
 * Separation of electronic and nuclear motion
 */
class BornOppenheimerApproximation {
public:
    /**
     * @brief Mass ratio
     *
     * m_e/M_n << 1 justifies the approximation
     */
    static double mass_ratio(double nuclear_mass) {
        return constants::m_e / nuclear_mass;
    }

    /**
     * @brief Total wave function factorization
     *
     * Ψ(r,R) ≈ ψ_electronic(r;R) × χ_nuclear(R)
     *
     * r: electronic coordinates
     * R: nuclear coordinates
     */
    static double total_wavefunction(
        double psi_electronic,
        double chi_nuclear) {

        return psi_electronic * chi_nuclear;
    }

    /**
     * @brief Electronic Hamiltonian (at fixed nuclear positions)
     *
     * H_el = T_el + V_el-el + V_el-nuc
     */
    static double electronic_energy_parametric(
        double nuclear_separation,
        const std::function<double(double)>& potential_curve) {

        return potential_curve(nuclear_separation);
    }

    /**
     * @brief Nuclear kinetic energy operator
     *
     * T_nuc = -ℏ²/(2M) ∇_R²
     */
    static double nuclear_kinetic_energy(
        double nuclear_mass,
        double velocity) {

        return 0.5 * nuclear_mass * velocity * velocity;
    }

    /**
     * @brief Validity criterion
     *
     * ω_vib << ω_el (vibrational frequency << electronic frequency)
     */
    static bool is_approximation_valid(
        double vibrational_frequency,
        double electronic_frequency) {

        return vibrational_frequency < 0.1 * electronic_frequency;
    }

    /**
     * @brief Adiabatic vs diabatic representation
     */
    enum class Representation {
        ADIABATIC,   // Electronic states adjust instantaneously
        DIABATIC     // Fixed electronic basis
    };
};

/**
 * @brief Nuclear Motion of Diatomic Molecules
 *
 * Rotational and vibrational spectroscopy
 */
class DiatomicMolecules {
public:
    /**
     * @brief Reduced mass
     *
     * μ = m₁m₂/(m₁ + m₂)
     */
    static double reduced_mass(double m1, double m2) {
        return (m1 * m2) / (m1 + m2);
    }

    /**
     * @brief Rotational energy levels
     *
     * E_J = BJ(J+1) where B = ℏ²/(2I)
     *
     * I = μR²: moment of inertia
     */
    static double rotational_energy(
        int J,
        double moment_of_inertia) {

        double B = (constants::hbar * constants::hbar) / (2.0 * moment_of_inertia);
        return B * J * (J + 1);
    }

    /**
     * @brief Rotational constant B
     *
     * B = h/(8π²cI) in wavenumber units (cm⁻¹)
     */
    static double rotational_constant(double moment_of_inertia) {
        return constants::h / (8.0 * M_PI * M_PI * constants::c * moment_of_inertia);
    }

    /**
     * @brief Vibrational energy levels (harmonic oscillator)
     *
     * E_v = ℏω(v + 1/2)
     */
    static double vibrational_energy_harmonic(
        int v,
        double angular_frequency) {

        return constants::hbar * angular_frequency * (v + 0.5);
    }

    /**
     * @brief Anharmonic oscillator energy
     *
     * E_v = ℏω(v + 1/2) - χₑℏω(v + 1/2)²
     *
     * χₑ: anharmonicity constant
     */
    static double vibrational_energy_anharmonic(
        int v,
        double omega,
        double chi_e) {

        double v_half = v + 0.5;
        return constants::hbar * omega * v_half * (1.0 - chi_e * v_half);
    }

    /**
     * @brief Rovibrational energy
     *
     * E(v,J) = ℏω(v + 1/2) + BJ(J+1)
     */
    static double rovibrational_energy(
        int v,
        int J,
        double omega,
        double B_constant) {

        return constants::hbar * omega * (v + 0.5) + B_constant * J * (J + 1);
    }

    /**
     * @brief Centrifugal distortion correction
     *
     * E_J = BJ(J+1) - DJ²(J+1)²
     *
     * D: centrifugal distortion constant
     */
    static double rotational_energy_with_distortion(
        int J,
        double B,
        double D) {

        return B * J * (J + 1) - D * J * J * (J + 1) * (J + 1);
    }

    /**
     * @brief Selection rules
     *
     * ΔJ = ±1 (rotational)
     * Δv = ±1 (vibrational, harmonic)
     */
    struct SelectionRules {
        static bool rotational_allowed(int J_initial, int J_final) {
            return std::abs(J_final - J_initial) == 1;
        }

        static bool vibrational_allowed_harmonic(int v_initial, int v_final) {
            return std::abs(v_final - v_initial) == 1;
        }
    };

    /**
     * @brief Morse potential
     *
     * V(R) = Dₑ[1 - e^(-a(R-Rₑ))]²
     *
     * Dₑ: dissociation energy
     * a: width parameter
     * Rₑ: equilibrium bond length
     */
    static double morse_potential(
        double R,
        double D_e,
        double a,
        double R_e) {

        double exponent = -a * (R - R_e);
        return D_e * std::pow(1.0 - std::exp(exponent), 2);
    }
};

/**
 * @brief The Hydrogen Molecular Ion H₂⁺
 *
 * Simplest molecular system
 */
class HydrogenMolecularIon {
public:
    /**
     * @brief LCAO (Linear Combination of Atomic Orbitals)
     *
     * ψ = c₁φ_A ± c₂φ_B
     *
     * φ_A, φ_B: atomic orbitals on nuclei A and B
     */
    static double lcao_molecular_orbital(
        double phi_A,
        double phi_B,
        bool bonding = true) {

        if (bonding) {
            return (phi_A + phi_B) / std::sqrt(2.0);  // σ_g (bonding)
        } else {
            return (phi_A - phi_B) / std::sqrt(2.0);  // σ_u* (antibonding)
        }
    }

    /**
     * @brief Bonding and antibonding energies
     *
     * E_± = (H_AA ± H_AB) / (1 ± S_AB)
     *
     * H_AA: Coulomb integral
     * H_AB: Resonance integral
     * S_AB: Overlap integral
     */
    static std::pair<double, double> bonding_antibonding_energies(
        double H_AA,
        double H_AB,
        double S_AB) {

        double E_bonding = (H_AA + H_AB) / (1.0 + S_AB);
        double E_antibonding = (H_AA - H_AB) / (1.0 - S_AB);

        return {E_bonding, E_antibonding};
    }

    /**
     * @brief Overlap integral S_AB
     *
     * S = ∫ φ_A φ_B dτ
     */
    static double overlap_integral_1s(
        double R,
        double a_0 = constants::a_0) {

        // For 1s orbitals separated by distance R
        double rho = R / a_0;

        return std::exp(-rho) * (1.0 + rho + rho * rho / 3.0);
    }

    /**
     * @brief Bonding energy curve
     *
     * E(R) for H₂⁺ as function of internuclear distance
     */
    static double bonding_energy_curve(
        double R,
        double dissociation_energy,
        double equilibrium_distance) {

        // Simplified Morse-like potential
        double alpha = 1.0 / equilibrium_distance;
        return -dissociation_energy *
               std::exp(-2.0 * alpha * (R - equilibrium_distance));
    }

    /**
     * @brief Equilibrium bond length
     *
     * R_e ≈ 2.5 a₀ for H₂⁺
     */
    static double equilibrium_bond_length() {
        return 2.5 * constants::a_0;
    }

    /**
     * @brief Dissociation energy
     *
     * D₀ ≈ 2.8 eV for H₂⁺
     */
    static double dissociation_energy() {
        return 2.8 * constants::e;  // Convert eV to Joules
    }
};

/**
 * @brief The Hydrogen Molecule H₂
 *
 * Molecular orbital and valence bond theory
 */
class HydrogenMolecule {
public:
    /**
     * @brief Molecular orbital configuration
     *
     * (σ_g)² for ground state H₂
     */
    static std::string ground_state_configuration() {
        return "(σ_g 1s)²";
    }

    /**
     * @brief Valence bond wave function
     *
     * ψ_VB = [φ_A(1)φ_B(2) + φ_A(2)φ_B(1)] × [α(1)β(2) - α(2)β(1)]
     *
     * Covalent structure with singlet spin
     */
    static double valence_bond_spatial(
        double phi_A1,
        double phi_B1,
        double phi_A2,
        double phi_B2) {

        return (phi_A1 * phi_B2 + phi_A2 * phi_B1) / std::sqrt(2.0);
    }

    /**
     * @brief Molecular orbital wave function
     *
     * ψ_MO = σ_g(1)σ_g(2) × [α(1)β(2) - α(2)β(1)]
     */
    static double molecular_orbital_spatial(
        double sigma_g_1,
        double sigma_g_2) {

        return sigma_g_1 * sigma_g_2;
    }

    /**
     * @brief Heitler-London approximation
     *
     * E = (Q + J) / (1 + S²) for H₂
     *
     * Q: Coulomb integral
     * J: Exchange integral
     * S: Overlap
     */
    static double heitler_london_energy(
        double Q,
        double J,
        double S) {

        return (Q + J) / (1.0 + S * S);
    }

    /**
     * @brief Bond dissociation energy
     *
     * D₀ = 4.75 eV for H₂
     */
    static double bond_dissociation_energy() {
        return 4.75 * constants::e;
    }

    /**
     * @brief Equilibrium bond length
     *
     * R_e = 0.74 Å for H₂
     */
    static double equilibrium_bond_length() {
        return 0.74e-10;  // meters
    }

    /**
     * @brief Ionic-covalent resonance
     *
     * ψ = c₁ψ_covalent + c₂ψ_ionic
     *
     * ψ_covalent: [A⁻B + AB⁻]
     * ψ_ionic: [A⁺B⁻ + A⁻B⁺]
     */
    static double resonance_wavefunction(
        double psi_covalent,
        double psi_ionic,
        double c_covalent,
        double c_ionic) {

        return c_covalent * psi_covalent + c_ionic * psi_ionic;
    }
};

/**
 * @brief The Chemical Bond
 *
 * Bonding concepts and molecular stability
 */
class ChemicalBond {
public:
    /**
     * @brief Bond order
     *
     * BO = (n_bonding - n_antibonding) / 2
     */
    static double bond_order(int n_bonding, int n_antibonding) {
        return (n_bonding - n_antibonding) / 2.0;
    }

    /**
     * @brief Sigma (σ) and pi (π) bonds
     */
    enum class BondType {
        SIGMA,        // Cylindrical symmetry along bond axis
        PI,           // Nodal plane containing bond axis
        DELTA         // Two nodal planes
    };

    /**
     * @brief Hybridization
     */
    enum class Hybridization {
        SP,      // Linear (180°)
        SP2,     // Trigonal planar (120°)
        SP3,     // Tetrahedral (109.5°)
        SP3D,    // Trigonal bipyramidal
        SP3D2    // Octahedral
    };

    /**
     * @brief Electronegativity difference
     *
     * Δχ determines ionic vs covalent character
     */
    static std::string bond_character(double delta_electronegativity) {
        if (delta_electronegativity < 0.5) {
            return "Nonpolar covalent";
        } else if (delta_electronegativity < 1.7) {
            return "Polar covalent";
        } else {
            return "Ionic";
        }
    }

    /**
     * @brief Percent ionic character
     *
     * % ionic ≈ 100[1 - e^(-0.25Δχ²)]
     */
    static double percent_ionic_character(double delta_electronegativity) {
        return 100.0 * (1.0 - std::exp(-0.25 * delta_electronegativity *
                                       delta_electronegativity));
    }

    /**
     * @brief Bond length correlation with bond order
     *
     * Higher bond order → shorter bond
     */
    static double bond_length_estimate(
        double single_bond_length,
        double bond_order) {

        // Empirical relationship
        return single_bond_length / std::pow(bond_order, 0.6);
    }

    /**
     * @brief Bond energy correlation with bond order
     *
     * Higher bond order → stronger bond
     */
    static double bond_energy_estimate(
        double single_bond_energy,
        double bond_order) {

        return single_bond_energy * bond_order;
    }

    /**
     * @brief Resonance structures
     *
     * Average over multiple Lewis structures
     */
    static double resonance_hybrid_energy(
        const std::vector<double>& structure_energies,
        const std::vector<double>& weights) {

        double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
        double energy = 0.0;

        for (size_t i = 0; i < structure_energies.size(); ++i) {
            energy += weights[i] * structure_energies[i];
        }

        return energy / total_weight;
    }
};

/**
 * @brief Structures of Simple Polyatomic Molecules
 *
 * VSEPR theory and molecular geometry
 */
class PolyatomicMolecules {
public:
    /**
     * @brief VSEPR (Valence Shell Electron Pair Repulsion) geometry
     */
    struct VSEPRGeometry {
        int electron_pairs;
        int bonding_pairs;
        int lone_pairs;
        std::string geometry;
        double bond_angle;

        static VSEPRGeometry determine(int bonding, int lone) {
            int total = bonding + lone;
            VSEPRGeometry geom;
            geom.bonding_pairs = bonding;
            geom.lone_pairs = lone;
            geom.electron_pairs = total;

            if (total == 2) {
                geom.geometry = "Linear";
                geom.bond_angle = 180.0;
            } else if (total == 3) {
                if (lone == 0) {
                    geom.geometry = "Trigonal planar";
                    geom.bond_angle = 120.0;
                } else {
                    geom.geometry = "Bent";
                    geom.bond_angle = 120.0;  // Slightly less due to lone pair
                }
            } else if (total == 4) {
                if (lone == 0) {
                    geom.geometry = "Tetrahedral";
                    geom.bond_angle = 109.5;
                } else if (lone == 1) {
                    geom.geometry = "Trigonal pyramidal";
                    geom.bond_angle = 107.0;
                } else {
                    geom.geometry = "Bent";
                    geom.bond_angle = 104.5;
                }
            } else if (total == 5) {
                geom.geometry = "Trigonal bipyramidal";
                geom.bond_angle = 90.0;  // Axial
            } else if (total == 6) {
                geom.geometry = "Octahedral";
                geom.bond_angle = 90.0;
            }

            return geom;
        }
    };

    /**
     * @brief Walsh diagrams
     *
     * Correlation between molecular geometry and orbital energies
     */
    static double orbital_energy_vs_angle(
        double angle,
        const std::string& orbital_type) {

        // Simplified model
        if (orbital_type == "bonding") {
            return -10.0 + 5.0 * std::cos(angle * M_PI / 180.0);
        } else {
            return -5.0 - 3.0 * std::cos(angle * M_PI / 180.0);
        }
    }

    /**
     * @brief Examples of common geometries
     */
    struct MolecularExamples {
        static VSEPRGeometry water() {
            return VSEPRGeometry::determine(2, 2);  // H₂O: bent
        }

        static VSEPRGeometry ammonia() {
            return VSEPRGeometry::determine(3, 1);  // NH₃: trigonal pyramidal
        }

        static VSEPRGeometry methane() {
            return VSEPRGeometry::determine(4, 0);  // CH₄: tetrahedral
        }

        static VSEPRGeometry carbon_dioxide() {
            return VSEPRGeometry::determine(2, 0);  // CO₂: linear
        }
    };

    /**
     * @brief Dipole moment
     *
     * μ = Σ q_i r_i (vector sum)
     */
    static double dipole_moment_magnitude(
        const std::vector<double>& charges,
        const std::vector<std::vector<double>>& positions) {

        double mu_x = 0.0, mu_y = 0.0, mu_z = 0.0;

        for (size_t i = 0; i < charges.size(); ++i) {
            mu_x += charges[i] * positions[i][0];
            mu_y += charges[i] * positions[i][1];
            mu_z += charges[i] * positions[i][2];
        }

        return std::sqrt(mu_x * mu_x + mu_y * mu_y + mu_z * mu_z);
    }
};

/**
 * @brief The Hückel Molecular Orbital Method
 *
 * π-electron systems in conjugated molecules
 */
class HuckelMolecularOrbital {
public:
    /**
     * @brief Hückel Hamiltonian matrix elements
     *
     * H_ii = α (Coulomb integral)
     * H_ij = β (resonance integral) for adjacent atoms
     * H_ij = 0 otherwise
     */
    static std::vector<std::vector<double>> huckel_matrix(
        int n_atoms,
        const std::vector<std::pair<int, int>>& bonds,
        double alpha,
        double beta) {

        std::vector<std::vector<double>> H(n_atoms, std::vector<double>(n_atoms, 0.0));

        // Diagonal elements (Coulomb integrals)
        for (int i = 0; i < n_atoms; ++i) {
            H[i][i] = alpha;
        }

        // Off-diagonal elements (resonance integrals)
        for (const auto& [i, j] : bonds) {
            H[i][j] = beta;
            H[j][i] = beta;
        }

        return H;
    }

    /**
     * @brief Hückel (4n+2) rule for aromaticity
     *
     * Aromatic if n_pi_electrons = 4n + 2 (n = 0, 1, 2, ...)
     */
    static bool is_aromatic_huckel(int n_pi_electrons) {
        // Check if n_pi_electrons = 4n + 2 for some integer n ≥ 0
        if (n_pi_electrons < 2) return false;
        return (n_pi_electrons - 2) % 4 == 0;
    }

    /**
     * @brief Total π-electron energy
     *
     * E_π = Σᵢ nᵢ εᵢ
     *
     * nᵢ: occupation number
     * εᵢ: orbital energy
     */
    static double total_pi_energy(
        const std::vector<double>& orbital_energies,
        const std::vector<int>& occupations) {

        double E_total = 0.0;
        for (size_t i = 0; i < orbital_energies.size(); ++i) {
            E_total += occupations[i] * orbital_energies[i];
        }
        return E_total;
    }

    /**
     * @brief Delocalization energy
     *
     * E_deloc = E_π(conjugated) - E_π(localized)
     */
    static double delocalization_energy(
        double conjugated_energy,
        double localized_energy) {

        return conjugated_energy - localized_energy;
    }

    /**
     * @brief Bond order in Hückel theory
     *
     * p_ij = Σₖ n_k c_ik c_jk
     *
     * n_k: occupation of MO k
     * c_ik, c_jk: LCAO coefficients
     */
    static double huckel_bond_order(
        const std::vector<double>& occupations,
        const std::vector<std::vector<double>>& coefficients,
        int atom_i,
        int atom_j) {

        double p_ij = 0.0;
        for (size_t k = 0; k < occupations.size(); ++k) {
            p_ij += occupations[k] * coefficients[k][atom_i] * coefficients[k][atom_j];
        }
        return p_ij;
    }

    /**
     * @brief Examples of conjugated systems
     */
    struct ConjugatedExamples {
        // Benzene: 6 π-electrons, C₆H₆
        static bool benzene_aromatic() {
            return is_aromatic_huckel(6);  // 4(1) + 2 = 6 ✓
        }

        // Cyclobutadiene: 4 π-electrons, C₄H₄
        static bool cyclobutadiene_aromatic() {
            return is_aromatic_huckel(4);  // Not 4n+2 ✗
        }

        // Naphthalene: 10 π-electrons, C₁₀H₈
        static bool naphthalene_aromatic() {
            return is_aromatic_huckel(10);  // 4(2) + 2 = 10 ✓
        }

        // Cyclopentadienyl anion: 6 π-electrons, C₅H₅⁻
        static bool cyclopentadienyl_anion_aromatic() {
            return is_aromatic_huckel(6);  // 4(1) + 2 = 6 ✓
        }
    };

    /**
     * @brief Charge density at atom i
     *
     * q_i = Σₖ n_k |c_ik|²
     */
    static double charge_density(
        const std::vector<double>& occupations,
        const std::vector<std::vector<double>>& coefficients,
        int atom_i) {

        double q_i = 0.0;
        for (size_t k = 0; k < occupations.size(); ++k) {
            q_i += occupations[k] * coefficients[k][atom_i] * coefficients[k][atom_i];
        }
        return q_i;
    }

    /**
     * @brief Resonance energy (aromatic stabilization)
     *
     * Difference between actual energy and hypothetical localized structure
     */
    static double resonance_energy_benzene(double beta) {
        // Benzene: E_π = 6α + 8β
        // Three isolated double bonds: E = 6α + 6β
        // Resonance energy = 2β
        return 2.0 * std::abs(beta);
    }
};

} // namespace quantum_chemistry
} // namespace physics

#endif // PHYSICS_QUANTUM_CHEMISTRY_HPP
