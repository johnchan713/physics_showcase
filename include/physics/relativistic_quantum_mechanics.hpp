#ifndef PHYSICS_RELATIVISTIC_QUANTUM_MECHANICS_HPP
#define PHYSICS_RELATIVISTIC_QUANTUM_MECHANICS_HPP

#include <vector>
#include <complex>
#include <functional>
#include <cmath>
#include <array>
#include <algorithm>

namespace physics {
namespace relativistic_quantum {

using Complex = std::complex<double>;

/**
 * @brief Physical constants for relativistic quantum mechanics
 */
namespace constants {
    constexpr double h = 6.62607015e-34;        // Planck constant (J·s)
    constexpr double hbar = 1.054571817e-34;    // Reduced Planck constant
    constexpr double c = 299792458.0;           // Speed of light (m/s)
    constexpr double e = 1.602176634e-19;       // Elementary charge (C)
    constexpr double m_e = 9.1093837015e-31;    // Electron mass (kg)
    constexpr double a_0 = 5.29177210903e-11;   // Bohr radius (m)
    constexpr double alpha = 7.2973525693e-3;   // Fine structure constant ≈ 1/137
    constexpr double mu_B = 9.2740100783e-24;   // Bohr magneton (J/T)
    constexpr double g_e = 2.00231930436256;    // Electron g-factor
}

/**
 * @brief Degenerate Position Eigenstates
 *
 * Handling degeneracy in quantum systems
 */
class DegeneratePositionEigenstates {
public:
    /**
     * @brief Degeneracy of energy level
     *
     * For hydrogen: g_n = n² (including spin: 2n²)
     */
    static int degeneracy_hydrogen(int n, bool include_spin = false) {
        int g = n * n;
        if (include_spin) {
            g *= 2;
        }
        return g;
    }

    /**
     * @brief Degeneracy for 3D harmonic oscillator
     *
     * g_N = (N+1)(N+2)/2
     */
    static int degeneracy_3d_oscillator(int N) {
        return (N + 1) * (N + 2) / 2;
    }

    /**
     * @brief Degenerate perturbation theory
     *
     * Need to diagonalize W within degenerate subspace
     */
    struct DegenerateSubspace {
        int dimension;
        std::vector<std::vector<Complex>> perturbation_matrix;
        std::vector<double> eigenvalues;
        std::vector<std::vector<Complex>> eigenvectors;
    };

    /**
     * @brief Lifted degeneracy
     *
     * Perturbation splits degenerate levels
     */
    static std::vector<double> split_degenerate_level(
        double E_0,
        const std::vector<double>& perturbation_eigenvalues) {

        std::vector<double> split_levels;
        for (double delta_E : perturbation_eigenvalues) {
            split_levels.push_back(E_0 + delta_E);
        }
        return split_levels;
    }

    /**
     * @brief Good quantum numbers
     *
     * Quantum numbers that commute with H
     */
    static bool is_good_quantum_number(
        const std::string& operator_name,
        const std::string& hamiltonian_symmetry) {

        // If [O, H] = 0, then O is a good quantum number
        return operator_name == hamiltonian_symmetry;
    }

    /**
     * @brief Accidental degeneracy
     *
     * Degeneracy not required by symmetry (e.g., hydrogen l-degeneracy)
     */
    static bool is_accidental_degeneracy(
        const std::string& system,
        int n,
        int l1,
        int l2) {

        if (system == "hydrogen") {
            // In hydrogen, different l values are degenerate (accidental)
            return l1 != l2 && l1 < n && l2 < n;
        }
        return false;
    }
};

/**
 * @brief Spin-Half Particles
 *
 * Comprehensive treatment of spin-1/2 systems
 */
class SpinHalfParticles {
public:
    /**
     * @brief Spin angular momentum operators
     *
     * S_x, S_y, S_z with [S_i, S_j] = iℏε_ijk S_k
     */
    struct SpinOperators {
        // Pauli matrices (in units of ℏ/2)
        static std::array<std::array<Complex, 2>, 2> sigma_x() {
            return {{
                {Complex(0, 0), Complex(1, 0)},
                {Complex(1, 0), Complex(0, 0)}
            }};
        }

        static std::array<std::array<Complex, 2>, 2> sigma_y() {
            return {{
                {Complex(0, 0), Complex(0, -1)},
                {Complex(0, 1), Complex(0, 0)}
            }};
        }

        static std::array<std::array<Complex, 2>, 2> sigma_z() {
            return {{
                {Complex(1, 0), Complex(0, 0)},
                {Complex(0, 0), Complex(-1, 0)}
            }};
        }

        // S_i = (ℏ/2) σ_i
        static std::array<std::array<Complex, 2>, 2> S_x() {
            auto sigma = sigma_x();
            double factor = constants::hbar / 2.0;
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    sigma[i][j] *= factor;
                }
            }
            return sigma;
        }
    };

    /**
     * @brief Spin eigenstates
     *
     * |↑⟩ = (1, 0)ᵀ, |↓⟩ = (0, 1)ᵀ
     */
    static std::array<Complex, 2> spin_up() {
        return {Complex(1, 0), Complex(0, 0)};
    }

    static std::array<Complex, 2> spin_down() {
        return {Complex(0, 0), Complex(1, 0)};
    }

    /**
     * @brief General spin state
     *
     * |χ⟩ = cos(θ/2)|↑⟩ + e^(iφ)sin(θ/2)|↓⟩
     */
    static std::array<Complex, 2> general_spin_state(double theta, double phi) {
        return {
            Complex(std::cos(theta / 2.0), 0),
            Complex(std::sin(theta / 2.0) * std::cos(phi),
                   std::sin(theta / 2.0) * std::sin(phi))
        };
    }

    /**
     * @brief Spin expectation values
     *
     * ⟨S_x⟩, ⟨S_y⟩, ⟨S_z⟩
     */
    static std::array<double, 3> spin_expectation(
        const std::array<Complex, 2>& state) {

        // ⟨S_z⟩ = (ℏ/2)(|α|² - |β|²)
        double S_z = (constants::hbar / 2.0) *
                    (std::norm(state[0]) - std::norm(state[1]));

        // ⟨S_x⟩ = (ℏ/2)Re(α*β)
        double S_x = (constants::hbar / 2.0) *
                    (std::real(std::conj(state[0]) * state[1]) +
                     std::real(std::conj(state[1]) * state[0]));

        // ⟨S_y⟩ = (ℏ/2)Im(α*β)
        double S_y = (constants::hbar / 2.0) *
                    (std::imag(std::conj(state[0]) * state[1]) -
                     std::imag(std::conj(state[1]) * state[0]));

        return {S_x, S_y, S_z};
    }

    /**
     * @brief Spin precession in magnetic field
     *
     * ω = -γB (Larmor frequency)
     * γ = g_e μ_B/ℏ (gyromagnetic ratio)
     */
    static double larmor_frequency(double magnetic_field) {
        double gamma = constants::g_e * constants::mu_B / constants::hbar;
        return gamma * magnetic_field;
    }

    /**
     * @brief Time evolution of spin state
     *
     * |χ(t)⟩ = e^(-iH t/ℏ)|χ(0)⟩
     */
    static std::array<Complex, 2> evolve_spin_state(
        const std::array<Complex, 2>& initial_state,
        double magnetic_field,
        double time) {

        double omega = larmor_frequency(magnetic_field);
        double phase = omega * time / 2.0;

        Complex phase_up = std::exp(Complex(0, -phase));
        Complex phase_down = std::exp(Complex(0, phase));

        return {
            initial_state[0] * phase_up,
            initial_state[1] * phase_down
        };
    }

    /**
     * @brief Spin-1/2 density matrix
     *
     * ρ = (1/2)(I + r⃗·σ⃗)
     *
     * r⃗: Bloch vector
     */
    static std::array<std::array<Complex, 2>, 2> density_matrix(
        double r_x,
        double r_y,
        double r_z) {

        return {{
            {Complex(0.5 * (1.0 + r_z), 0),
             Complex(0.5 * (r_x - r_y), 0)},
            {Complex(0.5 * (r_x + r_y), 0),
             Complex(0.5 * (1.0 - r_z), 0)}
        }};
    }

    /**
     * @brief Purity of spin state
     *
     * Pure state: Tr(ρ²) = 1
     * Mixed state: Tr(ρ²) < 1
     */
    static double purity(double r_x, double r_y, double r_z) {
        return r_x * r_x + r_y * r_y + r_z * r_z;
    }
};

/**
 * @brief Spin Magnetic Moment and Stern-Gerlach Experiment
 *
 * Interaction of spin with magnetic field
 */
class SpinMagneticMoment {
public:
    /**
     * @brief Magnetic moment operator
     *
     * μ⃗ = -g_e (μ_B/ℏ) S⃗
     */
    static double magnetic_moment_z(int m_s_sign) {
        return -constants::g_e * constants::mu_B * (m_s_sign * 0.5);
    }

    /**
     * @brief Interaction energy with magnetic field
     *
     * H = -μ⃗·B⃗ = g_e μ_B m_s B
     */
    static double interaction_energy(int m_s_sign, double B_field) {
        return constants::g_e * constants::mu_B * (m_s_sign * 0.5) * B_field;
    }

    /**
     * @brief Stern-Gerlach force
     *
     * F_z = μ_z ∂B_z/∂z
     */
    static double stern_gerlach_force(int m_s_sign, double gradient) {
        double mu_z = magnetic_moment_z(m_s_sign);
        return mu_z * gradient;
    }

    /**
     * @brief Beam deflection in Stern-Gerlach apparatus
     *
     * Δz = (1/2)(F/m)t²
     */
    static double beam_deflection(
        int m_s_sign,
        double gradient,
        double apparatus_length,
        double beam_velocity) {

        double force = stern_gerlach_force(m_s_sign, gradient);
        double time = apparatus_length / beam_velocity;
        return 0.5 * (force / constants::m_e) * time * time;
    }

    /**
     * @brief Sequential Stern-Gerlach experiments
     *
     * Demonstrates quantum measurement and collapse
     */
    struct SequentialSG {
        enum class Orientation { X, Y, Z };

        static double probability_after_measurement(
            Orientation first_axis,
            Orientation second_axis,
            bool first_result_up,
            bool second_result_up) {

            // If same axis: probability is 1 if same result, 0 otherwise
            if (first_axis == second_axis) {
                return (first_result_up == second_result_up) ? 1.0 : 0.0;
            }

            // If perpendicular axes: probability is 1/2
            return 0.5;
        }
    };

    /**
     * @brief Landé g-factor
     *
     * g_J = 1 + [J(J+1) - L(L+1) + S(S+1)] / [2J(J+1)]
     */
    static double lande_g_factor(int J, int L, int S) {
        if (J == 0) return 0.0;

        double numerator = J * (J + 1) - L * (L + 1) + S * (S + 1);
        double denominator = 2.0 * J * (J + 1);

        return 1.0 + numerator / denominator;
    }

    /**
     * @brief Magnetic moment for atom
     *
     * μ_J = -g_J μ_B √(J(J+1))
     */
    static double atomic_magnetic_moment(int J, int L, int S) {
        double g_J = lande_g_factor(J, L, S);
        return g_J * constants::mu_B * std::sqrt(J * (J + 1.0));
    }
};

/**
 * @brief Spin-Orbit Coupling
 *
 * Fine structure from relativistic corrections
 */
class SpinOrbitCoupling {
public:
    /**
     * @brief Spin-orbit Hamiltonian
     *
     * H_SO = (1/(2m²c²r)) (dV/dr) L⃗·S⃗
     */
    static double spin_orbit_energy(
        int n,
        int l,
        int j,
        int Z = 1) {

        if (l == 0) return 0.0;

        // For hydrogen-like atoms
        double E_0 = 13.6 * constants::e * Z * Z;  // Rydberg energy

        // Fine structure constant
        double alpha = constants::alpha;

        // Spin-orbit energy
        double energy = (alpha * alpha * E_0) / (n * n * n) *
                       (j * (j + 1) - l * (l + 1) - 0.75) /
                       (l * (l + 0.5) * (l + 1));

        return energy;
    }

    /**
     * @brief Total angular momentum j
     *
     * j = l ± 1/2
     */
    static std::pair<int, int> possible_j_values(int l) {
        // Returns (2j) for integer storage
        if (l == 0) {
            return {1, 1};  // Only j = 1/2
        }
        return {2 * l - 1, 2 * l + 1};  // j = l - 1/2, l + 1/2
    }

    /**
     * @brief Fine structure splitting
     *
     * ΔE = E(j=l+1/2) - E(j=l-1/2)
     */
    static double fine_structure_splitting(int n, int l, int Z = 1) {
        if (l == 0) return 0.0;

        // Use actual j values (stored as 2j)
        auto [j_minus, j_plus] = possible_j_values(l);

        double E_plus = spin_orbit_energy(n, l, j_plus / 2, Z);
        double E_minus = spin_orbit_energy(n, l, j_minus / 2, Z);

        return std::abs(E_plus - E_minus);
    }

    /**
     * @brief L·S expectation value
     *
     * ⟨L⃗·S⃗⟩ = (ℏ²/2)[j(j+1) - l(l+1) - s(s+1)]
     */
    static double LS_expectation(int l, int j_times_2) {
        double j = j_times_2 / 2.0;
        double s = 0.5;

        return (constants::hbar * constants::hbar / 2.0) *
               (j * (j + 1) - l * (l + 1) - s * (s + 1));
    }

    /**
     * @brief Fine structure constant correction
     *
     * Relativistic correction scale: α² ≈ (1/137)² ≈ 5 × 10⁻⁵
     */
    static double fine_structure_scale() {
        return constants::alpha * constants::alpha;
    }

    /**
     * @brief Thomas precession factor
     *
     * Factor of 1/2 from Thomas precession in moving reference frame
     */
    static double thomas_precession_factor() {
        return 0.5;
    }
};

/**
 * @brief Zeeman Effect Revisited
 *
 * Complete treatment including anomalous Zeeman effect
 */
class ZeemanEffectRevisited {
public:
    /**
     * @brief Normal Zeeman effect (no spin)
     *
     * ΔE = μ_B m_l B
     */
    static double normal_zeeman_shift(int m_l, double B_field) {
        return constants::mu_B * m_l * B_field;
    }

    /**
     * @brief Anomalous Zeeman effect (with spin)
     *
     * ΔE = g_J μ_B m_J B
     */
    static double anomalous_zeeman_shift(
        int m_J,
        int J,
        int L,
        int S,
        double B_field) {

        double g_J = SpinMagneticMoment::lande_g_factor(J, L, S);
        return g_J * constants::mu_B * m_J * B_field;
    }

    /**
     * @brief Paschen-Back effect (strong field)
     *
     * When B is strong: decouple L and S
     * ΔE = μ_B(m_l + g_s m_s)B
     */
    static double paschen_back_shift(
        int m_l,
        int m_s,
        double B_field) {

        return constants::mu_B * (m_l + constants::g_e * m_s) * B_field;
    }

    /**
     * @brief Transition between weak and strong field regimes
     *
     * Compare Zeeman energy to fine structure
     */
    static std::string field_regime(double B_field, int n, int l) {
        double zeeman_energy = constants::mu_B * B_field;
        double fine_structure = SpinOrbitCoupling::fine_structure_splitting(n, l);

        if (zeeman_energy < fine_structure) {
            return "Weak field (anomalous Zeeman)";
        } else {
            return "Strong field (Paschen-Back)";
        }
    }

    /**
     * @brief Zeeman pattern for spectral lines
     *
     * Selection rules: ΔJ = 0, ±1 (not 0→0), Δm_J = 0, ±1
     */
    struct ZeemanPattern {
        enum class Polarization { SIGMA_PLUS, SIGMA_MINUS, PI };

        static bool transition_allowed(int J_i, int m_J_i, int J_f, int m_J_f) {
            int delta_J = J_f - J_i;
            int delta_m = m_J_f - m_J_i;

            // ΔJ = 0, ±1 (but not 0→0)
            if (std::abs(delta_J) > 1) return false;
            if (J_i == 0 && J_f == 0) return false;

            // Δm_J = 0, ±1
            if (std::abs(delta_m) > 1) return false;

            return true;
        }

        static Polarization get_polarization(int delta_m) {
            if (delta_m == +1) return Polarization::SIGMA_PLUS;
            if (delta_m == -1) return Polarization::SIGMA_MINUS;
            return Polarization::PI;
        }
    };

    /**
     * @brief Hyperfine structure (nuclear spin interaction)
     *
     * ΔE_hf = A I⃗·J⃗
     */
    static double hyperfine_splitting(
        int I,
        int J,
        int F,
        double A_constant) {

        // ⟨I⃗·J⃗⟩ = (ℏ²/2)[F(F+1) - I(I+1) - J(J+1)]
        return (A_constant / 2.0) *
               (F * (F + 1) - I * (I + 1) - J * (J + 1));
    }

    /**
     * @brief 21 cm line of hydrogen
     *
     * Hyperfine transition: F=1 → F=0
     */
    static double hydrogen_21cm_frequency() {
        return 1420405751.768;  // Hz
    }
};

// ============================================================================
// KLEIN-GORDON EQUATION: Relativistic Wave Equation for Spin-0 Particles
// ============================================================================

/**
 * @brief Notation for Relativistic Quantum Mechanics
 *
 * Standard 4-vector and tensor notation
 */
class RelativisticNotation {
public:
    /**
     * @brief 4-vector components
     *
     * x^μ = (ct, x, y, z) = (x^0, x^1, x^2, x^3)
     */
    struct FourVector {
        double x0;  // ct
        double x1;  // x
        double x2;  // y
        double x3;  // z

        static FourVector position(double t, double x, double y, double z) {
            return {constants::c * t, x, y, z};
        }

        static FourVector momentum(double E, double px, double py, double pz) {
            return {E / constants::c, px, py, pz};
        }
    };

    /**
     * @brief Metric tensor (Minkowski spacetime)
     *
     * g_μν = diag(+1, -1, -1, -1) (mostly minus signature)
     * or g_μν = diag(-1, +1, +1, +1) (mostly plus signature)
     */
    static std::array<std::array<double, 4>, 4> metric_tensor(
        bool mostly_minus = true) {

        if (mostly_minus) {
            return {{
                {1.0, 0.0, 0.0, 0.0},
                {0.0, -1.0, 0.0, 0.0},
                {0.0, 0.0, -1.0, 0.0},
                {0.0, 0.0, 0.0, -1.0}
            }};
        } else {
            return {{
                {-1.0, 0.0, 0.0, 0.0},
                {0.0, 1.0, 0.0, 0.0},
                {0.0, 0.0, 1.0, 0.0},
                {0.0, 0.0, 0.0, 1.0}
            }};
        }
    }

    /**
     * @brief D'Alembertian operator
     *
     * □ = ∂_μ ∂^μ = (1/c²)∂²/∂t² - ∇²
     */
    static double dalembertian_dispersion(
        double omega,
        double k_x,
        double k_y,
        double k_z) {

        return (omega * omega) / (constants::c * constants::c) -
               (k_x * k_x + k_y * k_y + k_z * k_z);
    }

    /**
     * @brief Lorentz invariant scalar product
     *
     * x·y = x^μ y_μ = x^0 y^0 - x⃗·y⃗
     */
    static double lorentz_product(
        const FourVector& x,
        const FourVector& y) {

        return x.x0 * y.x0 - x.x1 * y.x1 - x.x2 * y.x2 - x.x3 * y.x3;
    }

    /**
     * @brief Natural units conversion
     *
     * Set ℏ = c = 1
     */
    static double natural_units_energy(double energy_SI) {
        return energy_SI / constants::hbar;
    }

    static double natural_units_mass(double mass_SI) {
        return mass_SI * constants::c * constants::c / constants::hbar;
    }
};

/**
 * @brief The Klein-Gordon Equation
 *
 * Relativistic wave equation for spin-0 particles
 * (□ + m²c²/ℏ²)ψ = 0
 */
class KleinGordonEquation {
public:
    /**
     * @brief Klein-Gordon equation operator
     *
     * (□ + μ²)ψ = 0 where μ² = m²c²/ℏ²
     */
    static bool satisfies_klein_gordon(
        double mass,
        double energy,
        double momentum) {

        // E² = (pc)² + (mc²)²
        double lhs = energy * energy;
        double rhs = momentum * momentum * constants::c * constants::c +
                     mass * mass * constants::c * constants::c * constants::c * constants::c;

        return std::abs(lhs - rhs) < 1e-10 * rhs;
    }

    /**
     * @brief Plane wave solutions
     *
     * ψ = A e^(i(k⃗·r⃗ - ωt)) with ω² = c²k² + (mc²/ℏ)²
     */
    static Complex plane_wave(
        double x,
        double t,
        double momentum,
        double energy) {

        double phase = (momentum * x - energy * t) / constants::hbar;
        return std::exp(Complex(0, phase));
    }

    /**
     * @brief Dispersion relation
     *
     * ω² = c²k² + (mc²/ℏ)²
     */
    static double frequency_from_wavevector(double k, double mass) {
        double omega_squared = constants::c * constants::c * k * k +
                              (mass * constants::c * constants::c * mass * constants::c * constants::c) /
                              (constants::hbar * constants::hbar);
        return std::sqrt(omega_squared);
    }

    /**
     * @brief Positive and negative energy solutions
     *
     * E = ±√((pc)² + (mc²)²) = ±E_p
     */
    static std::pair<double, double> energy_solutions(
        double momentum,
        double mass) {

        double E_squared = momentum * momentum * constants::c * constants::c +
                          mass * mass * constants::c * constants::c * constants::c * constants::c;
        double E = std::sqrt(E_squared);

        return {E, -E};  // E_+ and E_-
    }

    /**
     * @brief Conserved current density
     *
     * j^μ = (iℏ/2m)(ψ*∂^μψ - ψ∂^μψ*)
     */
    static double current_density_time(
        Complex psi,
        Complex psi_t,
        double mass) {

        // j^0 = (iℏ/2mc²)(ψ*∂_t ψ - ψ ∂_t ψ*)
        Complex j0 = Complex(0, constants::hbar / (2.0 * mass * constants::c * constants::c)) *
                     (std::conj(psi) * psi_t - psi * std::conj(psi_t));
        return std::real(j0);
    }

    /**
     * @brief Probability density (not positive definite!)
     *
     * ρ = j^0/c = (iℏ/2mc²)(ψ*∂_t ψ - ψ ∂_t ψ*)
     */
    static double probability_density(
        Complex psi,
        Complex psi_t,
        double mass) {

        return current_density_time(psi, psi_t, mass) / constants::c;
    }

    /**
     * @brief Klein paradox
     *
     * For V₀ > E + 2mc², transmission coefficient can exceed 1
     */
    static double klein_paradox_transmission(
        double energy,
        double barrier_height,
        double mass,
        double barrier_width) {

        double mc2 = mass * constants::c * constants::c;

        if (barrier_height > energy + 2.0 * mc2) {
            // Pair production regime - perfect transmission
            return 1.0;
        }

        // Normal tunneling (simplified)
        if (barrier_height > energy) {
            double kappa = std::sqrt(2.0 * mass * (barrier_height - energy)) / constants::hbar;
            return std::exp(-2.0 * kappa * barrier_width);
        }

        return 1.0;  // E > V₀
    }

    /**
     * @brief Continuity equation
     *
     * ∂ρ/∂t + ∇·j⃗ = 0
     */
    static bool satisfies_continuity(
        double rho_t,
        double div_j,
        double tol = 1e-10) {

        return std::abs(rho_t + div_j) < tol;
    }
};

/**
 * @brief Nonrelativistic Limit of Klein-Gordon Equation
 *
 * Recover Schrödinger equation for v << c
 */
class KleinGordonNonrelativisticLimit {
public:
    /**
     * @brief Nonrelativistic ansatz
     *
     * ψ = φ(x,t) e^(-imc²t/ℏ)
     */
    static Complex nonrelativistic_wavefunction(
        Complex phi,
        double t,
        double mass) {

        double phase = -(mass * constants::c * constants::c * t) / constants::hbar;
        return phi * std::exp(Complex(0, phase));
    }

    /**
     * @brief Schrödinger equation emerges
     *
     * iℏ∂φ/∂t = -(ℏ²/2m)∇²φ
     */
    static double schrodinger_kinetic_energy(
        double momentum,
        double mass) {

        return (momentum * momentum) / (2.0 * mass);
    }

    /**
     * @brief Relativistic correction
     *
     * E ≈ mc² + p²/2m - p⁴/8m³c² + ...
     */
    static double energy_expansion(
        double momentum,
        double mass,
        int order = 2) {

        double mc2 = mass * constants::c * constants::c;
        double p2 = momentum * momentum;

        // Leading order: rest mass
        double E = mc2;

        // First order: nonrelativistic kinetic
        if (order >= 1) {
            E += p2 / (2.0 * mass);
        }

        // Second order: relativistic correction
        if (order >= 2) {
            E -= p2 * p2 / (8.0 * mass * mass * mass * constants::c * constants::c);
        }

        return E;
    }

    /**
     * @brief Velocity β = v/c for given momentum
     */
    static double velocity_ratio(double momentum, double mass) {
        double E = std::sqrt(momentum * momentum * constants::c * constants::c +
                           mass * mass * constants::c * constants::c * constants::c * constants::c);
        return (momentum * constants::c) / E;
    }

    /**
     * @brief Check if nonrelativistic approximation is valid
     *
     * Valid if p << mc (or v << c)
     */
    static bool is_nonrelativistic(double momentum, double mass) {
        return momentum < 0.1 * mass * constants::c;  // β < 0.1
    }
};

/**
 * @brief Free Spin-0 Particles (Klein-Gordon Solutions)
 *
 * Complete solution set for free particles
 */
class KleinGordonFreeParticles {
public:
    /**
     * @brief General solution as superposition
     *
     * ψ(x,t) = ∫[A(k)e^(i(kx-ωt)) + B(k)e^(i(kx+ωt))]dk
     */
    struct FreeSolution {
        double amplitude_positive;   // A(k)
        double amplitude_negative;   // B(k)
        double wavevector;          // k
        double frequency;           // ω
    };

    /**
     * @brief Energy eigenvalue for given momentum
     *
     * E_p = +√((pc)² + (mc²)²)
     */
    static double positive_energy(double momentum, double mass) {
        return std::sqrt(momentum * momentum * constants::c * constants::c +
                        mass * mass * constants::c * constants::c * constants::c * constants::c);
    }

    /**
     * @brief Group velocity
     *
     * v_g = dω/dk = pc²/E
     */
    static double group_velocity(double momentum, double mass) {
        double E = positive_energy(momentum, mass);
        return (momentum * constants::c * constants::c) / E;
    }

    /**
     * @brief Phase velocity
     *
     * v_p = ω/k = E/p
     */
    static double phase_velocity(double momentum, double mass) {
        double E = positive_energy(momentum, mass);
        return (E * constants::c) / (momentum * constants::c);
    }

    /**
     * @brief Wave packet for localized particle
     *
     * Gaussian wave packet: ψ(x,0) = exp(-x²/4σ²)exp(ik₀x)
     */
    static Complex gaussian_wave_packet(
        double x,
        double x0,
        double sigma,
        double k0) {

        double gaussian = std::exp(-(x - x0) * (x - x0) / (4.0 * sigma * sigma));
        double phase = k0 * x;
        return gaussian * std::exp(Complex(0, phase));
    }

    /**
     * @brief Normalization for free particle states
     *
     * ⟨p|p'⟩ = (2π)³ δ³(p⃗ - p⃗')
     */
    static double normalization_factor(double volume) {
        return 1.0 / std::sqrt(volume);
    }

    /**
     * @brief Klein-Gordon inner product
     *
     * (ψ₁, ψ₂) = i∫[ψ₁*∂_t ψ₂ - (∂_t ψ₁*)ψ₂]d³x
     */
    static Complex klein_gordon_inner_product(
        Complex psi1,
        Complex psi1_t,
        Complex psi2,
        Complex psi2_t,
        double volume) {

        return Complex(0, 1) * (std::conj(psi1) * psi2_t - std::conj(psi1_t) * psi2) * volume;
    }
};

/**
 * @brief Energy-Momentum Tensor of Klein-Gordon Field
 *
 * Stress-energy tensor T^μν for scalar field
 */
class KleinGordonEnergyMomentum {
public:
    /**
     * @brief Energy density (T^00)
     *
     * T^00 = (1/2)[(∂_t ψ)² + c²(∇ψ)² + (mc²/ℏ)²ψ²]
     */
    static double energy_density(
        Complex psi,
        Complex psi_t,
        double grad_psi_squared,
        double mass) {

        double term1 = std::norm(psi_t);
        double term2 = constants::c * constants::c * grad_psi_squared;
        double term3 = (mass * constants::c * constants::c / constants::hbar) *
                       (mass * constants::c * constants::c / constants::hbar) *
                       std::norm(psi);

        return 0.5 * (term1 + term2 + term3);
    }

    /**
     * @brief Momentum density (T^0i)
     *
     * T^0i = (∂_t ψ*)(∂_i ψ) + (∂_i ψ*)(∂_t ψ)
     */
    static double momentum_density(
        Complex psi_t,
        Complex psi_x) {

        return std::real(std::conj(psi_t) * psi_x + std::conj(psi_x) * psi_t);
    }

    /**
     * @brief Stress tensor (T^ij)
     *
     * T^ij = c²(∂_i ψ*)(∂_j ψ) + c²(∂_j ψ*)(∂_i ψ) - δ^ij L
     */
    static double stress_tensor_component(
        Complex psi_i,
        Complex psi_j,
        double lagrangian_density,
        bool is_diagonal) {

        double term = constants::c * constants::c *
                     std::real(std::conj(psi_i) * psi_j + std::conj(psi_j) * psi_i);

        if (is_diagonal) {
            term -= lagrangian_density;
        }

        return term;
    }

    /**
     * @brief Conservation of energy-momentum
     *
     * ∂_μ T^μν = 0
     */
    static bool is_conserved(
        double d_T00_dt,
        double div_T0i,
        double tol = 1e-10) {

        return std::abs(d_T00_dt + div_T0i) < tol;
    }

    /**
     * @brief Hamiltonian density
     *
     * H = T^00 (for canonical formulation)
     */
    static double hamiltonian_density(
        Complex psi,
        Complex pi,
        double grad_psi_squared,
        double mass) {

        // H = π·π* + c²|∇ψ|² + (mc²/ℏ)²|ψ|²
        return std::norm(pi) +
               constants::c * constants::c * grad_psi_squared +
               (mass * constants::c * constants::c / constants::hbar) *
               (mass * constants::c * constants::c / constants::hbar) *
               std::norm(psi);
    }
};

/**
 * @brief Klein-Gordon Equation in Schrödinger Form
 *
 * First-order in time formulation
 */
class KleinGordonSchrodingerForm {
public:
    /**
     * @brief Two-component formulation
     *
     * Ψ = (ψ, π)^T where π = ∂ψ/∂t
     */
    struct TwoComponentState {
        Complex psi;
        Complex pi;  // Conjugate momentum
    };

    /**
     * @brief Schrödinger-like equation
     *
     * iℏ ∂Ψ/∂t = H_KG Ψ
     */
    static TwoComponentState time_evolution(
        const TwoComponentState& state,
        double laplacian_psi,
        double mass,
        double dt) {

        // π̇ = ∇²ψ - (mc/ℏ)²ψ
        // ψ̇ = π

        Complex pi_dot = laplacian_psi -
                        (mass * constants::c / constants::hbar) *
                        (mass * constants::c / constants::hbar) * state.psi;

        TwoComponentState evolved;
        evolved.psi = state.psi + state.pi * dt;
        evolved.pi = state.pi + pi_dot * dt;

        return evolved;
    }

    /**
     * @brief Hamiltonian operator in matrix form
     *
     * H_KG = [[0, 1], [c²∇² - (mc²)², 0]]
     */
    static std::array<std::array<double, 2>, 2> hamiltonian_matrix(
        double k_squared,
        double mass) {

        double h11 = 0.0;
        double h12 = 1.0;
        double h21 = -constants::c * constants::c * k_squared -
                     (mass * constants::c * constants::c / constants::hbar) *
                     (mass * constants::c * constants::c / constants::hbar);
        double h22 = 0.0;

        return {{{h11, h12}, {h21, h22}}};
    }

    /**
     * @brief Positive definite norm
     *
     * ||Ψ||² = ∫[|π|² + c²|∇ψ|² + (mc²/ℏ)²|ψ|²]d³x
     */
    static double positive_norm_squared(
        const TwoComponentState& state,
        double grad_psi_squared,
        double mass) {

        return std::norm(state.pi) +
               constants::c * constants::c * grad_psi_squared +
               (mass * constants::c * constants::c / constants::hbar) *
               (mass * constants::c * constants::c / constants::hbar) *
               std::norm(state.psi);
    }
};

/**
 * @brief Charge Conjugation for Klein-Gordon Field
 *
 * Particle-antiparticle symmetry
 */
class ChargeConjugation {
public:
    /**
     * @brief Charge conjugation operator
     *
     * C: ψ → ψ* (for real scalar field ψ = ψ*)
     */
    static Complex charge_conjugate(Complex psi) {
        return std::conj(psi);
    }

    /**
     * @brief Particle and antiparticle states
     *
     * Positive frequency: particle
     * Negative frequency: antiparticle
     */
    struct ParticleAntiparticle {
        Complex particle_amplitude;      // e^(-iEt/ℏ)
        Complex antiparticle_amplitude;  // e^(+iEt/ℏ)
    };

    /**
     * @brief Charge conjugation invariance
     *
     * For neutral scalar: C|ψ⟩ = |ψ⟩
     */
    static bool is_charge_eigenstate(Complex psi, double phase_tol = 1e-10) {
        // ψ = ψ* for real scalar (self-conjugate)
        return std::abs(psi - std::conj(psi)) < phase_tol;
    }

    /**
     * @brief C-parity for neutral particles
     *
     * C = ±1 for self-conjugate states
     */
    static int c_parity(bool symmetric) {
        return symmetric ? 1 : -1;
    }

    /**
     * @brief Transformation of current under C
     *
     * j^μ → -j^μ (changes sign)
     */
    static double conjugate_current(double current) {
        return -current;
    }
};

/**
 * @brief Feshbach-Villars Representation
 *
 * Positive-definite probability density formulation
 */
class FeshbachVillarsRepresentation {
public:
    /**
     * @brief Feshbach-Villars transformation
     *
     * φ = (1/√2)(ψ + iπ/mc²)
     * χ = (1/√2)(ψ - iπ/mc²)
     */
    struct FVComponents {
        Complex phi;  // Positive energy component
        Complex chi;  // Negative energy component
    };

    /**
     * @brief Transform to FV representation
     */
    static FVComponents to_feshbach_villars(
        Complex psi,
        Complex pi,
        double mass) {

        double mc2 = mass * constants::c * constants::c;
        Complex factor = Complex(0, 1) / mc2;

        FVComponents fv;
        fv.phi = (psi + factor * pi) / std::sqrt(2.0);
        fv.chi = (psi - factor * pi) / std::sqrt(2.0);

        return fv;
    }

    /**
     * @brief Positive-definite probability density
     *
     * ρ_FV = |φ|² + |χ|² ≥ 0
     */
    static double probability_density_fv(const FVComponents& fv) {
        return std::norm(fv.phi) + std::norm(fv.chi);
    }

    /**
     * @brief FV Hamiltonian
     *
     * H_FV = βmc² + α⃗·(cp⃗)
     */
    static std::array<std::array<double, 2>, 2> fv_beta_matrix() {
        return {{{1.0, 0.0}, {0.0, -1.0}}};
    }

    /**
     * @brief Coupled equations
     *
     * iℏ∂φ/∂t = mc²φ - iℏc∇χ
     * iℏ∂χ/∂t = -mc²χ + iℏc∇φ
     */
    static FVComponents time_evolution_fv(
        const FVComponents& state,
        double grad_phi,
        double grad_chi,
        double mass,
        double dt) {

        double mc2 = mass * constants::c * constants::c;

        FVComponents evolved;
        Complex phi_dot = Complex(0, -mc2 / constants::hbar) * state.phi -
                         constants::c * grad_chi;
        Complex chi_dot = Complex(0, mc2 / constants::hbar) * state.chi +
                         constants::c * grad_phi;

        evolved.phi = state.phi + phi_dot * dt;
        evolved.chi = state.chi + chi_dot * dt;

        return evolved;
    }

    /**
     * @brief Nonrelativistic limit in FV representation
     *
     * χ → 0, φ → ψ_Schrödinger
     */
    static Complex nonrelativistic_component(const FVComponents& fv) {
        return fv.phi;  // Positive energy dominates
    }
};

/**
 * @brief Klein-Gordon Equation with Electromagnetic Field
 *
 * Minimal coupling to electromagnetism
 */
class KleinGordonElectromagneticField {
public:
    /**
     * @brief Minimal coupling (gauge covariant derivative)
     *
     * ∂_μ → D_μ = ∂_μ + (iq/ℏc)A_μ
     */
    static Complex covariant_derivative_time(
        Complex psi,
        Complex psi_t,
        double charge,
        double phi_em) {

        // D_0 ψ = (∂_t + iqφ/ℏ)ψ
        return psi_t + Complex(0, charge * phi_em / constants::hbar) * psi;
    }

    /**
     * @brief Klein-Gordon with EM field
     *
     * [(∂_μ + iqA_μ)(∂^μ + iqA^μ) + (mc/ℏ)²]ψ = 0
     */
    static bool satisfies_kg_with_field(
        double energy,
        double momentum,
        double mass,
        double charge,
        double phi_em,
        double A_vector) {

        // (E - qφ)² = (p - qA)²c² + (mc²)²
        double E_kin = energy - charge * phi_em;
        double p_kin = momentum - charge * A_vector;

        double lhs = E_kin * E_kin;
        double rhs = p_kin * p_kin * constants::c * constants::c +
                     mass * mass * constants::c * constants::c * constants::c * constants::c;

        return std::abs(lhs - rhs) < 1e-10 * rhs;
    }

    /**
     * @brief Conserved current with EM field
     *
     * j^μ = (iq/2m)[ψ*(D^μψ) - (D^μψ)*ψ]
     */
    static double current_density_with_field(
        Complex psi,
        Complex D_mu_psi,
        double charge,
        double mass) {

        return std::real(Complex(0, charge / (2.0 * mass)) *
                        (std::conj(psi) * D_mu_psi - std::conj(D_mu_psi) * psi));
    }

    /**
     * @brief Landau levels for charged scalar in B field
     *
     * E_n = √[(mc²)² + 2n|q|ℏcB]
     */
    static double landau_level_energy(
        int n,
        double charge,
        double magnetic_field,
        double mass) {

        double mc2_squared = mass * mass * constants::c * constants::c *
                            constants::c * constants::c;
        double term = 2.0 * n * std::abs(charge) * constants::hbar *
                     constants::c * magnetic_field;

        return std::sqrt(mc2_squared + term);
    }

    /**
     * @brief Cyclotron frequency
     *
     * ω_c = |q|B/(γm)
     */
    static double cyclotron_frequency(
        double charge,
        double magnetic_field,
        double mass,
        double gamma = 1.0) {

        return std::abs(charge) * magnetic_field / (gamma * mass);
    }
};

/**
 * @brief Gauge Invariance of the Coupling
 *
 * Local U(1) gauge symmetry
 */
class GaugeInvariance {
public:
    /**
     * @brief Gauge transformation of wave function
     *
     * ψ → ψ' = e^(iqΛ/ℏ)ψ
     */
    static Complex gauge_transform_wavefunction(
        Complex psi,
        double charge,
        double lambda) {

        double phase = charge * lambda / constants::hbar;
        return psi * std::exp(Complex(0, phase));
    }

    /**
     * @brief Gauge transformation of potential
     *
     * A_μ → A'_μ = A_μ - ∂_μΛ
     */
    static double gauge_transform_vector_potential(
        double A,
        double grad_lambda) {

        return A - grad_lambda;
    }

    static double gauge_transform_scalar_potential(
        double phi,
        double lambda_t) {

        return phi + lambda_t;
    }

    /**
     * @brief Gauge invariance of Klein-Gordon equation
     *
     * Equation invariant under simultaneous gauge transformation
     */
    static bool verify_gauge_invariance(
        Complex psi,
        Complex psi_prime,
        double lambda,
        double charge,
        double tol = 1e-10) {

        Complex expected = gauge_transform_wavefunction(psi, charge, lambda);
        return std::abs(psi_prime - expected) < tol;
    }

    /**
     * @brief Gauge-invariant phase
     *
     * Aharonov-Bohm phase: exp(iq/ℏ ∮A⃗·dl⃗)
     */
    static double aharonov_bohm_phase(
        double charge,
        double magnetic_flux) {

        return charge * magnetic_flux / constants::hbar;
    }

    /**
     * @brief Gauge-invariant field strength
     *
     * F_μν = ∂_μA_ν - ∂_νA_μ
     */
    static double field_strength(
        double dA_nu_dx_mu,
        double dA_mu_dx_nu) {

        return dA_nu_dx_mu - dA_mu_dx_nu;
    }
};

/**
 * @brief Nonrelativistic Limit with Electromagnetic Fields
 *
 * Recover Schrödinger equation with EM coupling
 */
class KleinGordonNonrelativisticLimitWithFields {
public:
    /**
     * @brief Pauli equation emerges
     *
     * iℏ∂ψ/∂t = [(p⃗ - qA⃗)²/2m + qφ]ψ
     */
    static double schrodinger_energy_with_field(
        double momentum,
        double charge,
        double A_vector,
        double phi_em,
        double mass) {

        double p_kin = momentum - charge * A_vector;
        return (p_kin * p_kin) / (2.0 * mass) + charge * phi_em;
    }

    /**
     * @brief Darwin term (relativistic correction)
     *
     * H_Darwin = -(ℏ²/8m²c²)∇²V
     */
    static double darwin_correction(
        double laplacian_potential,
        double mass) {

        return -(constants::hbar * constants::hbar) /
               (8.0 * mass * mass * constants::c * constants::c) *
               laplacian_potential;
    }

    /**
     * @brief Spin-orbit coupling analogue (for scalar)
     *
     * Absent for spin-0, present for spin-1/2
     */
    static double spin_orbit_energy_scalar() {
        return 0.0;  // No spin-orbit for scalar particles
    }

    /**
     * @brief Relativistic kinetic energy correction
     *
     * -(p⃗ - qA⃗)⁴/8m³c²
     */
    static double relativistic_kinetic_correction(
        double momentum,
        double charge,
        double A_vector,
        double mass) {

        double p_kin = momentum - charge * A_vector;
        return -(p_kin * p_kin * p_kin * p_kin) /
               (8.0 * mass * mass * mass * constants::c * constants::c);
    }

    /**
     * @brief Zeeman energy for scalar in magnetic field
     *
     * ΔE = -(q²ℏ²/8m²c²)B²r² (diamagnetic)
     */
    static double diamagnetic_energy(
        double charge,
        double magnetic_field,
        double r_perp,
        double mass) {

        return -(charge * charge * magnetic_field * magnetic_field * r_perp * r_perp) /
               (8.0 * mass * constants::c * constants::c);
    }
};

/**
 * @brief Interpretation of One-Particle Operators
 *
 * Physical observables in relativistic QM
 */
class OneParticleOperators {
public:
    /**
     * @brief Position operator issue
     *
     * No local position operator for relativistic particles
     */
    static std::string position_operator_interpretation() {
        return "Newton-Wigner position (non-local)";
    }

    /**
     * @brief Momentum operator
     *
     * p̂ = -iℏ∇ (well-defined)
     */
    static double momentum_eigenvalue(double k) {
        return constants::hbar * k;
    }

    /**
     * @brief Energy operator
     *
     * Ê = iℏ∂/∂t (positive and negative values)
     */
    static std::pair<double, double> energy_eigenvalues(
        double momentum,
        double mass) {

        double E = std::sqrt(momentum * momentum * constants::c * constants::c +
                           mass * mass * constants::c * constants::c * constants::c * constants::c);
        return {E, -E};
    }

    /**
     * @brief Charge density operator
     *
     * ρ̂ = iq(ψ*∂_t ψ - ψ∂_t ψ*)/2mc²
     */
    static std::string charge_density_interpretation() {
        return "Not positive-definite; needs second quantization";
    }

    /**
     * @brief Current density operator
     *
     * ĵ = (q/2m)[ψ*(-iℏ∇)ψ + c.c.]
     */
    static std::string current_interpretation() {
        return "Well-defined for Klein-Gordon field";
    }

    /**
     * @brief Angular momentum
     *
     * L̂ = r⃗ × p̂ (orbital, no spin for scalar)
     */
    static double orbital_angular_momentum_quantum_number(int l) {
        return constants::hbar * std::sqrt(l * (l + 1.0));
    }

    /**
     * @brief Second quantization interpretation
     *
     * Single-particle interpretation breaks down
     */
    static std::string second_quantization_necessity() {
        return "Negative energy states → antiparticles; need QFT";
    }

    /**
     * @brief Zitterbewegung absent
     *
     * Unlike Dirac, Klein-Gordon has no Zitterbewegung
     */
    static bool has_zitterbewegung() {
        return false;
    }

    /**
     * @brief Localization length scale
     *
     * Compton wavelength: λ_C = ℏ/(mc)
     */
    static double compton_wavelength(double mass) {
        return constants::hbar / (mass * constants::c);
    }
};

/**
 * @brief The Dirac Equation
 *
 * Relativistic equation for spin-1/2 particles
 */
class DiracEquation {
public:
    /**
     * @brief Dirac equation
     *
     * (iℏ∂/∂t)ψ = (cα⃗·p⃗ + βmc²)ψ
     */

    /**
     * @brief Dirac matrices (4×4)
     *
     * α_i and β matrices
     */
    struct DiracMatrices {
        // In standard representation (Dirac-Pauli)
        // α_i = [[0, σ_i], [σ_i, 0]]
        // β = [[I, 0], [0, -I]]

        static std::vector<std::vector<Complex>> beta() {
            return {
                {Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                {Complex(0, 0), Complex(0, 0), Complex(-1, 0), Complex(0, 0)},
                {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(-1, 0)}
            };
        }

        // Gamma matrices: γ⁰ = β, γⁱ = βα_i
        static std::vector<std::vector<Complex>> gamma0() {
            return beta();
        }
    };

    /**
     * @brief Dirac spinor (4-component)
     *
     * ψ = (ψ_1, ψ_2, ψ_3, ψ_4)ᵀ
     */
    using DiracSpinor = std::array<Complex, 4>;

    /**
     * @brief Free particle solutions (plane waves)
     *
     * u(p) e^(-iEt/ℏ + ip·x/ℏ) for positive energy
     * v(p) e^(+iEt/ℏ + ip·x/ℏ) for negative energy
     */
    static DiracSpinor positive_energy_spinor(
        double momentum,
        int spin_up) {

        double E = std::sqrt(momentum * momentum * constants::c * constants::c +
                            constants::m_e * constants::m_e * constants::c * constants::c * constants::c * constants::c);

        double N = std::sqrt((E + constants::m_e * constants::c * constants::c) / (2.0 * E));

        if (spin_up) {
            return {
                Complex(N, 0),
                Complex(0, 0),
                Complex(N * momentum * constants::c / (E + constants::m_e * constants::c * constants::c), 0),
                Complex(0, 0)
            };
        } else {
            return {
                Complex(0, 0),
                Complex(N, 0),
                Complex(0, 0),
                Complex(N * momentum * constants::c / (E + constants::m_e * constants::c * constants::c), 0)
            };
        }
    }

    /**
     * @brief Positive definite probability density
     *
     * ρ = ψ†ψ (always positive!)
     */
    static double probability_density(const DiracSpinor& psi) {
        double rho = 0.0;
        for (const auto& component : psi) {
            rho += std::norm(component);
        }
        return rho;
    }

    /**
     * @brief Current density
     *
     * j⃗ = cψ†α⃗ψ
     */
    static double current_density_x(const DiracSpinor& psi) {
        // Simplified: j_x = c Re(ψ₁*ψ₃ + ψ₂*ψ₄)
        return constants::c * std::real(std::conj(psi[0]) * psi[2] +
                                        std::conj(psi[1]) * psi[3]);
    }

    /**
     * @brief Continuity equation
     *
     * ∂ρ/∂t + ∇·j⃗ = 0
     */
    static bool satisfies_continuity() {
        // Dirac equation automatically satisfies continuity
        return true;
    }

    /**
     * @brief Non-relativistic limit
     *
     * Recovers Pauli equation with spin
     */
    static std::string nonrelativistic_limit() {
        return "Pauli equation: (p²/2m - eΦ)ψ + (eℏ/2m)σ⃗·B⃗ψ";
    }
};

/**
 * @brief Spin and the Dirac Particle
 *
 * Intrinsic spin from Dirac equation
 */
class DiracParticleSpin {
public:
    /**
     * @brief Spin operator in Dirac theory
     *
     * S⃗ = (ℏ/2)Σ⃗ where Σ_i = [[σ_i, 0], [0, σ_i]]
     */
    static double spin_magnitude() {
        return constants::hbar * std::sqrt(0.5 * 1.5);  // s=1/2
    }

    /**
     * @brief Helicity operator
     *
     * h = Σ⃗·p̂ (spin projection along momentum)
     */
    static int helicity_eigenvalue(bool right_handed) {
        return right_handed ? +1 : -1;
    }

    /**
     * @brief Gyromagnetic ratio
     *
     * Dirac equation predicts g = 2 exactly
     */
    static double dirac_g_factor() {
        return 2.0;
    }

    /**
     * @brief QED correction to g-factor
     *
     * g = 2(1 + α/2π + ...) ≈ 2.00232
     */
    static double qed_g_factor() {
        return 2.0 * (1.0 + constants::alpha / (2.0 * M_PI));
    }

    /**
     * @brief Magnetic moment
     *
     * μ⃗ = -(e/m)S⃗ (from Dirac equation)
     */
    static double magnetic_moment() {
        return constants::e * constants::hbar / (2.0 * constants::m_e);
    }

    /**
     * @brief Zitterbewegung
     *
     * Rapid oscillation with frequency 2mc²/ℏ
     */
    static double zitterbewegung_frequency() {
        return 2.0 * constants::m_e * constants::c * constants::c / constants::hbar;
    }

    /**
     * @brief Zitterbewegung amplitude
     *
     * λ_C = ℏ/(mc) (Compton wavelength)
     */
    static double zitterbewegung_amplitude() {
        return constants::hbar / (constants::m_e * constants::c);
    }

    /**
     * @brief Spin precession in electromagnetic field
     *
     * Thomas-Bargmann-Michel-Telegdi equation
     */
    static double tbmt_precession_frequency(
        double E_field,
        double B_field,
        double velocity) {

        double gamma = 1.0 / std::sqrt(1.0 - velocity * velocity / (constants::c * constants::c));
        double omega_s = -(constants::e / constants::m_e) *
                        (B_field + (gamma - 1) * E_field / constants::c);
        return omega_s;
    }
};

/**
 * @brief Spin-Orbit Coupling in the Dirac Hamiltonian
 *
 * Automatic appearance of spin-orbit coupling
 */
class DiracSpinOrbitCoupling {
public:
    /**
     * @brief Dirac Hamiltonian for central potential
     *
     * Automatic spin-orbit term from Dirac equation
     */
    static double dirac_so_coupling_strength(double r, int Z) {
        // V(r) = -Ze²/(4πε₀r)
        // Leads to L·S coupling with correct Thomas factor
        double V = -Z * constants::e * constants::e * constants::k_e / r;
        double dV_dr = -V / r;

        return dV_dr / (2.0 * constants::m_e * constants::m_e * constants::c * constants::c * r);
    }

    /**
     * @brief Fine structure from Dirac equation
     *
     * E_nj = mc²[1 + (αZ)²/(n - j - 1/2 + √((j+1/2)² - (αZ)²))²]^(-1/2) - mc²
     */
    static double dirac_energy(int n, int j_times_2, int Z) {
        double j = j_times_2 / 2.0;
        double alpha_Z = constants::alpha * Z;

        double denominator_term = n - j - 0.5 +
                                 std::sqrt((j + 0.5) * (j + 0.5) - alpha_Z * alpha_Z);

        double bracket = 1.0 + (alpha_Z * alpha_Z) / (denominator_term * denominator_term);

        return constants::m_e * constants::c * constants::c *
               (1.0 / std::sqrt(bracket) - 1.0);
    }

    /**
     * @brief Darwin term
     *
     * Additional relativistic correction for s-states
     */
    static double darwin_term(int l, int Z) {
        if (l != 0) return 0.0;

        // (πℏ²/2m²c²)|ψ(0)|²(Ze²/4πε₀)
        // Represents Zitterbewegung averaging
        return (M_PI * constants::hbar * constants::hbar * Z * constants::e * constants::e) /
               (2.0 * constants::m_e * constants::m_e * constants::c * constants::c);
    }

    /**
     * @brief Kinetic energy relativistic correction
     *
     * -(p⁴/8m³c²)
     */
    static double kinetic_correction(double momentum) {
        return -std::pow(momentum, 4) / (8.0 * std::pow(constants::m_e, 3) * constants::c * constants::c);
    }

    /**
     * @brief Total fine structure
     *
     * Combination of spin-orbit + Darwin + kinetic
     */
    static double total_fine_structure(int n, int l, int j_times_2, int Z) {
        // Use Dirac energy which includes all corrections automatically
        return dirac_energy(n, j_times_2, Z);
    }
};

/**
 * @brief The Dirac Hydrogen Atom
 *
 * Exact solution including fine structure
 */
class DiracHydrogenAtom {
public:
    /**
     * @brief Quantum number j
     *
     * j = l ± 1/2 (total angular momentum)
     */
    static int quantum_number_j(int l, bool j_plus) {
        // Returns 2j for integer storage
        return j_plus ? (2 * l + 1) : (2 * l - 1);
    }

    /**
     * @brief Dirac energy levels
     *
     * E_nj = mc²[(1 + (αZ)²/(n - j - 1/2 + κ)²)^(-1/2) - 1]
     */
    static double energy_level(int n, int j_times_2, int Z = 1) {
        return DiracSpinOrbitCoupling::dirac_energy(n, j_times_2, Z);
    }

    /**
     * @brief Fine structure constant role
     *
     * α ≈ 1/137 determines size of fine structure
     */
    static double fine_structure_constant() {
        return constants::alpha;
    }

    /**
     * @brief Lamb shift
     *
     * QED correction beyond Dirac equation
     * Mainly affects s-states
     */
    static double lamb_shift(int n, int l) {
        if (l != 0) return 0.0;  // Mainly s-states

        // Approximate formula for hydrogen
        // Δν ≈ 1057.8 MHz for 2S_1/2 - 2P_1/2
        if (n == 2) {
            return 1057.8e6 * constants::h;  // Convert to energy
        }

        return 0.0;
    }

    /**
     * @brief Degeneracy with same j
     *
     * States with same n and j but different l are degenerate in Dirac theory
     * (Lamb shift breaks this)
     */
    static bool is_degenerate(int n1, int l1, int j1, int n2, int l2, int j2) {
        return (n1 == n2) && (j1 == j2);
    }

    /**
     * @brief Spectrum of hydrogen
     *
     * Including fine structure
     */
    struct HydrogenSpectrum {
        int n;
        int l;
        int j_times_2;
        double energy;
        std::string label;

        static HydrogenSpectrum ground_state() {
            return {1, 0, 1, energy_level(1, 1), "1S_1/2"};
        }

        static std::vector<HydrogenSpectrum> n2_states() {
            return {
                {2, 0, 1, energy_level(2, 1), "2S_1/2"},
                {2, 1, 1, energy_level(2, 1), "2P_1/2"},
                {2, 1, 3, energy_level(2, 3), "2P_3/2"}
            };
        }
    };

    /**
     * @brief Radial wave functions
     *
     * Modified from hydrogen to include relativity
     */
    static double radial_wavefunction(double r, int n, int kappa) {
        // Simplified - full solution requires confluent hypergeometric functions
        double a_0 = constants::a_0;
        return std::exp(-r / (n * a_0));
    }
};

/**
 * @brief The Dirac Particle in a Magnetic Field
 *
 * Landau levels and gyromagnetic ratio
 */
class DiracParticleInMagneticField {
public:
    /**
     * @brief Minimal coupling
     *
     * p⃗ → p⃗ - eA⃗
     */
    static double canonical_momentum(double kinetic_momentum, double vector_potential) {
        return kinetic_momentum + constants::e * vector_potential;
    }

    /**
     * @brief Landau levels for Dirac particle
     *
     * E_n = ±√((mc²)² + 2n|e|ℏcB)
     */
    static std::pair<double, double> landau_level_energy(int n, double B_field) {
        double mc2 = constants::m_e * constants::c * constants::c;
        double orbital_term = 2.0 * n * constants::e * constants::hbar * constants::c * B_field;

        double E_squared = mc2 * mc2 + orbital_term;
        double E = std::sqrt(E_squared);

        return {E, -E};
    }

    /**
     * @brief Cyclotron frequency
     *
     * ω_c = eB/m
     */
    static double cyclotron_frequency(double B_field) {
        return constants::e * B_field / constants::m_e;
    }

    /**
     * @brief Magnetic length
     *
     * l_B = √(ℏ/eB)
     */
    static double magnetic_length(double B_field) {
        return std::sqrt(constants::hbar / (constants::e * B_field));
    }

    /**
     * @brief Degeneracy of Landau level
     *
     * g = eBA/(2πℏ) (per unit area A)
     */
    static double landau_degeneracy_per_area(double B_field) {
        return constants::e * B_field / (2.0 * M_PI * constants::hbar);
    }

    /**
     * @brief Pauli term from Dirac equation
     *
     * -(eℏ/2m)σ⃗·B⃗ appears automatically
     */
    static double pauli_interaction_energy(int spin_sign, double B_field) {
        return -(constants::e * constants::hbar / (2.0 * constants::m_e)) *
               spin_sign * B_field;
    }

    /**
     * @brief Anomalous magnetic moment
     *
     * QED correction: g = 2(1 + α/2π) ≈ 2.00232
     */
    static double anomalous_magnetic_moment() {
        return (constants::alpha / (2.0 * M_PI)) * constants::mu_B;
    }

    /**
     * @brief Klein paradox in magnetic field
     *
     * Pair production in strong magnetic fields
     */
    static double critical_magnetic_field() {
        // B_c = m²c³/(eℏ) ≈ 4.4 × 10¹³ G
        return (constants::m_e * constants::m_e * constants::c * constants::c * constants::c) /
               (constants::e * constants::hbar);
    }

    /**
     * @brief Synchrotron radiation
     *
     * Energy loss for relativistic particle in B field
     */
    static double synchrotron_power(double gamma, double B_field) {
        // P = (2e⁴/3m²c³)γ²B²
        double e4 = std::pow(constants::e, 4);
        return (2.0 * e4 / (3.0 * constants::m_e * constants::m_e * constants::c * constants::c * constants::c)) *
               gamma * gamma * B_field * B_field;
    }
};

// ============================================================================
// LORENTZ COVARIANCE OF THE DIRAC EQUATION
// ============================================================================

/**
 * @brief Lorentz Covariance - Form Invariance
 *
 * Dirac equation maintains form under Lorentz transformations
 */
class LorentzCovariance {
public:
    /**
     * @brief Form invariance requirement
     *
     * (iℏγ^μ∂'_μ - mc)ψ'(x') = 0 if (iℏγ^μ∂_μ - mc)ψ(x) = 0
     */
    static std::string form_invariance() {
        return "Dirac equation has same form in all Lorentz frames";
    }

    /**
     * @brief Spinor transformation
     *
     * ψ'(x') = S(Λ)ψ(x) where x'^μ = Λ^μ_ν x^ν
     */
    static std::string spinor_transformation() {
        return "ψ'(x') = S(Λ)ψ(Λ⁻¹x')";
    }

    /**
     * @brief Transformation of gamma matrices
     *
     * S(Λ)γ^μS⁻¹(Λ) = Λ^μ_ν γ^ν
     */
    static bool verify_gamma_transformation(
        const std::vector<std::vector<Complex>>& S,
        const std::vector<std::vector<Complex>>& gamma_mu,
        double Lambda_mu_nu) {

        // Verify transformation property
        return true;  // Placeholder
    }

    /**
     * @brief Adjoint spinor transformation
     *
     * ψ̄'(x') = ψ̄(x)S⁻¹(Λ)
     */
    static std::string adjoint_transformation() {
        return "ψ̄' = ψ̄S⁻¹, ensuring ψ̄ψ is scalar";
    }

    /**
     * @brief Lorentz scalar
     *
     * ψ̄ψ is Lorentz invariant
     */
    static std::string lorentz_scalar() {
        return "ψ̄ψ = ψ†γ⁰ψ is Lorentz scalar";
    }

    /**
     * @brief Lorentz vector
     *
     * j^μ = ψ̄γ^μψ is a 4-vector
     */
    static std::string lorentz_vector() {
        return "j^μ = cψ̄γ^μψ transforms as 4-vector";
    }
};

/**
 * @brief Infinitesimal Lorentz Transformations
 *
 * Construction of Ŝ operator for infinitesimal transformations
 */
class InfinitesimalLorentzTransformations {
public:
    /**
     * @brief Infinitesimal Lorentz transformation
     *
     * Λ^μ_ν = δ^μ_ν + ω^μ_ν where ω_μν = -ω_νμ
     */
    static double infinitesimal_parameter(
        double omega_munu,
        bool antisymmetric = true) {

        // ω_μν = -ω_νμ (antisymmetric)
        return omega_munu;
    }

    /**
     * @brief Generator of Lorentz transformations
     *
     * S(ω) = I + (i/2)ω_μν σ^μν where σ^μν = (i/2)[γ^μ, γ^ν]
     */
    static std::string generator() {
        return "S = I + (i/2)ω_μν σ^μν, σ^μν = (i/2)[γ^μ, γ^ν]";
    }

    /**
     * @brief Commutator of gamma matrices
     *
     * σ^μν = (i/2)[γ^μ, γ^ν]
     */
    static std::vector<std::vector<Complex>> sigma_munu(
        const std::vector<std::vector<Complex>>& gamma_mu,
        const std::vector<std::vector<Complex>>& gamma_nu) {

        // Compute (i/2)[γ^μ, γ^ν]
        std::vector<std::vector<Complex>> result(4, std::vector<Complex>(4));
        return result;  // Placeholder
    }

    /**
     * @brief Infinitesimal rotation
     *
     * ω_ij = ε_ijk θ^k (spatial rotation by angle θ)
     */
    static std::string infinitesimal_rotation() {
        return "Rotation: S = I + (i/2)θ⃗·Σ⃗ where Σ_i = [[σ_i, 0], [0, σ_i]]";
    }

    /**
     * @brief Infinitesimal boost
     *
     * ω_0i = -ω_i0 = β_i (boost with velocity v⃗)
     */
    static std::string infinitesimal_boost() {
        return "Boost: S = I - (i/2)β⃗·K⃗ where K_i = (i/2)[γ⁰, γⁱ]";
    }

    /**
     * @brief Boost generator
     *
     * K_i = (i/2)γ⁰γⁱ = (i/2)βα_i
     */
    static std::string boost_generator() {
        return "K_i = (i/2)γ⁰γⁱ generates boosts";
    }
};

/**
 * @brief Finite Proper Lorentz Transformations
 *
 * Exponentiation of infinitesimal transformations
 */
class FiniteLorentzTransformations {
public:
    /**
     * @brief Finite rotation
     *
     * S_rot(θ⃗) = exp(-(i/2)θ⃗·Σ⃗)
     */
    static std::string finite_rotation() {
        return "S_rot = exp(-(i/2)θ⃗·Σ⃗) = cos(θ/2)I - i sin(θ/2)n̂·Σ⃗";
    }

    /**
     * @brief Rotation by angle θ around axis n̂
     */
    static double rotation_angle(double theta_x, double theta_y, double theta_z) {
        return std::sqrt(theta_x * theta_x + theta_y * theta_y + theta_z * theta_z);
    }

    /**
     * @brief Finite boost
     *
     * S_boost(β⃗) = exp((1/2)η n̂·K⃗) where η = tanh⁻¹(β)
     */
    static std::string finite_boost() {
        return "S_boost = cosh(η/2)I + sinh(η/2)n̂·K⃗";
    }

    /**
     * @brief Rapidity parameter
     *
     * η = tanh⁻¹(v/c) = (1/2)ln[(1+β)/(1-β)]
     */
    static double rapidity(double beta) {
        return 0.5 * std::log((1.0 + beta) / (1.0 - beta));
    }

    /**
     * @brief Lorentz factor from rapidity
     *
     * γ = cosh(η)
     */
    static double gamma_from_rapidity(double eta) {
        return std::cosh(eta);
    }

    /**
     * @brief Proper Lorentz group SO(1,3)
     */
    static std::string proper_lorentz_group() {
        return "SO(1,3): det(Λ) = +1, Λ⁰₀ ≥ +1 (proper orthochronous)";
    }

    /**
     * @brief SL(2,C) covering group
     *
     * Double cover of SO(1,3) by SL(2,C)
     */
    static std::string covering_group() {
        return "SL(2,C) is universal covering of SO↑₊(1,3)";
    }

    /**
     * @brief Spinor representation
     *
     * Spinors transform under (1/2, 0) ⊕ (0, 1/2) representation
     */
    static std::string spinor_representation() {
        return "4-component spinor: left-handed (1/2, 0) ⊕ right-handed (0, 1/2)";
    }
};

/**
 * @brief Four-Current Density
 *
 * Conserved 4-vector current from Dirac equation
 */
class DiracFourCurrent {
public:
    using DiracSpinor = std::array<Complex, 4>;

    /**
     * @brief Four-current density
     *
     * j^μ = cψ̄γ^μψ = (cρ, j⃗)
     */
    static std::array<double, 4> four_current(const DiracSpinor& psi) {
        // j^0 = cψ†ψ (charge density × c)
        double j0 = 0.0;
        for (const auto& comp : psi) {
            j0 += std::norm(comp);
        }
        j0 *= constants::c;

        // j^i = cψ†α_iψ (spatial current)
        // Simplified for j^1 only
        double j1 = constants::c * std::real(
            std::conj(psi[0]) * psi[2] + std::conj(psi[1]) * psi[3] +
            std::conj(psi[2]) * psi[0] + std::conj(psi[3]) * psi[1]
        );

        return {j0, j1, 0.0, 0.0};
    }

    /**
     * @brief Continuity equation
     *
     * ∂_μ j^μ = 0
     */
    static std::string continuity_equation() {
        return "∂_μ j^μ = ∂ρ/∂t + ∇·j⃗ = 0";
    }

    /**
     * @brief Charge conservation
     *
     * Q = ∫ρd³x is conserved
     */
    static std::string charge_conservation() {
        return "dQ/dt = 0, Q = e∫ψ†ψ d³x";
    }

    /**
     * @brief Transformation property
     *
     * j^μ transforms as 4-vector
     */
    static std::string transformation() {
        return "j'^μ = Λ^μ_ν j^ν under Lorentz transformation";
    }

    /**
     * @brief Gordon decomposition
     *
     * Separate convection and spin magnetization currents
     */
    static std::string gordon_decomposition() {
        return "j^μ = (1/2m)[ψ̄p^μψ + ∂^μ(ψ̄ψ)] + (i/2m)∂_ν(ψ̄σ^μνψ)";
    }

    /**
     * @brief Convection current
     *
     * j^μ_conv = (ψ̄p^μψ)/2m
     */
    static std::string convection_current() {
        return "Convection: (p^μ/m)ρ";
    }

    /**
     * @brief Magnetization current
     *
     * j^μ_mag = (i/2m)∂_ν(ψ̄σ^μνψ)
     */
    static std::string magnetization_current() {
        return "Magnetization: ∇ × M⃗ where M⃗ ~ ψ̄Σ⃗ψ";
    }
};

// ============================================================================
// SOLUTIONS BY LORENTZ TRANSFORMATIONS
// ============================================================================

/**
 * @brief Plane Waves in Arbitrary Directions
 *
 * Construct solutions for arbitrary momentum directions
 */
class PlaneWavesArbitraryDirections {
public:
    using DiracSpinor = std::array<Complex, 4>;

    /**
     * @brief Spinor for momentum in arbitrary direction
     *
     * u(p⃗, s) obtained by boosting rest frame spinor
     */
    static DiracSpinor arbitrary_direction_spinor(
        double px,
        double py,
        double pz,
        double mass,
        bool spin_up) {

        double p = std::sqrt(px * px + py * py + pz * pz);
        double E = std::sqrt(p * p * constants::c * constants::c +
                           mass * mass * constants::c * constants::c * constants::c * constants::c);

        double N = std::sqrt((E + mass * constants::c * constants::c) / (2.0 * E));

        // Two-component spinor χ
        Complex chi1 = spin_up ? Complex(1, 0) : Complex(0, 0);
        Complex chi2 = spin_up ? Complex(0, 0) : Complex(1, 0);

        // Small components: (σ⃗·p⃗)/(E + mc²) χ
        double sigma_dot_p = pz;  // Simplified for p⃗ = p_z ẑ
        double factor = (sigma_dot_p * constants::c) / (E + mass * constants::c * constants::c);

        return {
            N * chi1,
            N * chi2,
            Complex(N * factor, 0) * chi1,
            Complex(N * factor, 0) * chi2
        };
    }

    /**
     * @brief Helicity eigenstates
     *
     * h = Σ⃗·p̂ with eigenvalues ±1
     */
    static DiracSpinor helicity_eigenstate(
        double px,
        double py,
        double pz,
        double mass,
        int helicity) {

        // helicity = ±1 for right/left-handed
        return arbitrary_direction_spinor(px, py, pz, mass, helicity > 0);
    }

    /**
     * @brief Boost from rest frame
     *
     * u(p⃗) = S_boost(β⃗) u(0⃗)
     */
    static std::string boost_construction() {
        return "u(p⃗, s) = S_boost(p⃗/|p⃗|, β) u_rest(s)";
    }

    /**
     * @brief Rotation to arbitrary direction
     *
     * u(p⃗) = S_rot(R) u(p ẑ) where R rotates ẑ to p̂
     */
    static std::string rotation_construction() {
        return "u(p̂, s) = S_rot(ẑ → p̂) u(|p|ẑ, s)";
    }
};

/**
 * @brief General Form of Free Solutions
 *
 * Complete characterization of free Dirac solutions
 */
class GeneralFreeSolutions {
public:
    /**
     * @brief Complete solution set
     *
     * General solution: ψ = Σ[a_s u(p,s)e^(-iEt+ip·x) + b_s v(p,s)e^(iEt+ip·x)]
     */
    static std::string general_solution() {
        return "ψ = ∫d³p Σ_s [a(p,s)u(p,s)e^(-iEt+ip·x) + b(p,s)v(p,s)e^(iEt+ip·x)]";
    }

    /**
     * @brief Orthonormality relations
     *
     * ū(p,r)u(p,s) = 2Eδ_rs
     * v̄(p,r)v(p,s) = -2Eδ_rs
     * ū(p,r)v(-p,s) = 0
     */
    static std::string orthonormality() {
        return "ū(p,r)u(p,s) = 2mc²δ_rs, v̄(p,r)v(p,s) = -2mc²δ_rs";
    }

    /**
     * @brief Completeness relations
     *
     * Projection operators sum to identity
     */
    static std::string completeness() {
        return "Σ_s[u(p,s)ū(p,s) - v(p,s)v̄(p,s)] = 2mc²I";
    }

    /**
     * @brief Negative energy interpretation
     *
     * v(p) states: antiparticles (positrons)
     */
    static std::string negative_energy() {
        return "v(p,s)e^(iEt-ip·x): positron with momentum -p⃗, energy E > 0";
    }

    /**
     * @brief Charge conjugation relation
     *
     * v(p,s) = C ū^T(p,s) where C is charge conjugation matrix
     */
    static std::string charge_conjugation_relation() {
        return "v(p,s) related to u by charge conjugation";
    }
};

/**
 * @brief Polarized Electrons in Relativistic Theory
 *
 * Spin polarization for relativistic fermions
 */
class PolarizedElectrons {
public:
    /**
     * @brief Polarization vector
     *
     * s^μ = (γs⃗·v⃗/c, s⃗ + γ²s⃗·v⃗/(c(γ+1))v⃗)
     */
    static std::array<double, 4> polarization_four_vector(
        double sx,
        double sy,
        double sz,
        double vx,
        double gamma) {

        double s_dot_v = sz * vx;  // Simplified
        double s0 = gamma * s_dot_v / constants::c;
        double s1 = sx + gamma * gamma * s_dot_v / (constants::c * (gamma + 1.0)) * vx;

        return {s0, s1, sy, sz};
    }

    /**
     * @brief Spin projection operator
     *
     * P(n̂) = (1 + γ₅n̸)/2 where n̸ = γ^μn_μ
     */
    static std::string spin_projection() {
        return "P(s) = (1 + γ₅s̸)/2 projects onto spin direction s⃗";
    }

    /**
     * @brief Helicity for massless particles
     *
     * h = ±1 is conserved (chirality = helicity when m = 0)
     */
    static std::string helicity_massless() {
        return "Massless: helicity = chirality (Weyl fermions)";
    }

    /**
     * @brief Polarization density matrix
     *
     * ρ = (1 + γ₅P̸)/2 where P⃗ is polarization vector
     */
    static std::string density_matrix() {
        return "ρ = Σ_s f_s |u(p,s)⟩⟨ū(p,s)|";
    }

    /**
     * @brief Transverse polarization
     *
     * s⃗_⊥ perpendicular to momentum
     */
    static double transverse_polarization(
        double sx,
        double sy,
        double pz) {

        // s_⊥ = s⃗ - (s⃗·p̂)p̂
        return std::sqrt(sx * sx + sy * sy);
    }

    /**
     * @brief Longitudinal polarization
     *
     * s_∥ = s⃗·p̂ (helicity component)
     */
    static double longitudinal_polarization(
        double sz,
        double pz,
        double p) {

        return sz * pz / p;
    }
};

/**
 * @brief Projection Operators for Energy and Spin
 *
 * Operators projecting onto energy and spin eigenstates
 */
class ProjectionOperators {
public:
    /**
     * @brief Energy projection operators
     *
     * Λ_± = (±γ·p + mc²)/(2E) project onto E ≷ 0
     */
    static std::string energy_projectors() {
        return "Λ_+ = (γ·p + mc)/(2E), Λ_- = (-γ·p + mc)/(2E)";
    }

    /**
     * @brief Properties of energy projectors
     *
     * Λ_+² = Λ_+, Λ_-² = Λ_-, Λ_+Λ_- = 0, Λ_+ + Λ_- = I
     */
    static std::string energy_projector_properties() {
        return "Λ_±² = Λ_±, Λ_+Λ_- = 0, Λ_+ + Λ_- = I";
    }

    /**
     * @brief Spin projection operators
     *
     * P(n̂,±) = (I ± Σ⃗·n̂)/2 for non-relativistic
     */
    static std::string spin_projectors() {
        return "P_± = (I ± Σ⃗·n̂)/2 project onto spin ±ℏ/2 along n̂";
    }

    /**
     * @brief Relativistic spin projector
     *
     * S(s^μ) = (1 + γ₅s̸)/2
     */
    static std::string relativistic_spin_projector() {
        return "S = (1 + γ₅s̸)/2 for 4-vector s^μ";
    }

    /**
     * @brief Simultaneous energy-spin projection
     *
     * Λ_± P_s projects onto definite energy and spin
     */
    static std::string simultaneous_projection() {
        return "Λ_+P_+ projects onto E > 0, spin up";
    }

    /**
     * @brief Gordon identity
     *
     * ū(p')γ^μu(p) = ū(p')[(p+p')^μ/(2m) + iσ^μν q_ν/(2m)]u(p)
     */
    static std::string gordon_identity() {
        return "Separates convection and magnetization terms";
    }
};

/**
 * @brief Wave Packets of Plane Dirac Waves
 *
 * Localized wave packets from superposition
 */
class DiracWavePackets {
public:
    using DiracSpinor = std::array<Complex, 4>;

    /**
     * @brief Wave packet construction
     *
     * ψ(x,t) = ∫d³p a(p)u(p,s)e^(i(p·x-Et)/ℏ)
     */
    static std::string wave_packet() {
        return "ψ(x,t) = ∫d³p a(p)u(p,s)exp[i(p·x-E_pt)/ℏ]";
    }

    /**
     * @brief Gaussian wave packet
     *
     * a(p) = exp[-(p-p₀)²/(2σ_p²)]
     */
    static Complex gaussian_amplitude(
        double px,
        double p0x,
        double sigma_p) {

        double exponent = -(px - p0x) * (px - p0x) / (2.0 * sigma_p * sigma_p);
        return std::exp(exponent);
    }

    /**
     * @brief Group velocity
     *
     * v_g = dE/dp = pc²/E
     */
    static double group_velocity(double p, double mass) {
        double E = std::sqrt(p * p * constants::c * constants::c +
                           mass * mass * constants::c * constants::c * constants::c * constants::c);
        return (p * constants::c * constants::c) / E;
    }

    /**
     * @brief Spreading of wave packet
     *
     * Δx(t) ≈ Δx(0) + vΔp t/ℏ
     */
    static double packet_spreading(
        double delta_x0,
        double delta_p,
        double time,
        double velocity) {

        return delta_x0 + velocity * delta_p * time / constants::hbar;
    }

    /**
     * @brief Zitterbewegung in wave packets
     *
     * Interference between E > 0 and E < 0 components
     */
    static std::string zitterbewegung() {
        return "Rapid oscillation ω_Z = 2mc²/ℏ from E+/E- interference";
    }

    /**
     * @brief Localization limit
     *
     * Cannot localize better than Compton wavelength
     */
    static double minimum_localization(double mass) {
        return constants::hbar / (mass * constants::c);  // λ_C
    }
};

// ============================================================================
// DIRAC PARTICLES IN EXTERNAL FIELDS
// ============================================================================

/**
 * @brief Dirac Equation in External Fields
 *
 * Minimal coupling to electromagnetic fields
 */
class DiracExternalFields {
public:
    /**
     * @brief Minimal coupling
     *
     * p^μ → p^μ - (e/c)A^μ
     */
    static std::string minimal_coupling() {
        return "[iℏγ^μ(∂_μ + ieA_μ/ℏc) - mc]ψ = 0";
    }

    /**
     * @brief Dirac equation in EM field
     *
     * [γ^μ(iℏ∂_μ - eA_μ) - mc]ψ = 0
     */
    static std::string field_equation() {
        return "(iℏ∂/∂t + eφ)ψ = [cα⃗·(p⃗ - eA⃗/c) + βmc²]ψ";
    }

    /**
     * @brief Conserved current with field
     *
     * j^μ = cψ̄γ^μψ, ∂_μj^μ = 0
     */
    static std::string current_conservation() {
        return "∂_μj^μ = 0 still holds with EM field";
    }

    /**
     * @brief Gauge invariance
     *
     * ψ → e^(ieΛ/ℏ)ψ, A_μ → A_μ + ∂_μΛ
     */
    static std::string gauge_invariance() {
        return "U(1) gauge: ψ' = e^(ieΛ/ℏ)ψ, A'_μ = A_μ - ∂_μΛ";
    }

    /**
     * @brief Hydrogen atom (Coulomb field)
     *
     * V(r) = -Ze²/(4πε₀r)
     */
    static std::string coulomb_problem() {
        return "Dirac hydrogen: exact solution with fine structure";
    }

    /**
     * @brief Constant magnetic field
     *
     * Landau levels for Dirac particle
     */
    static double landau_level_dirac(
        int n,
        double B_field,
        double mass) {

        double mc2 = mass * constants::c * constants::c;
        double term = 2.0 * n * constants::e * constants::hbar * constants::c * B_field;
        return std::sqrt(mc2 * mc2 + term);
    }
};

/**
 * @brief Two-Centre Dirac Equation
 *
 * Dirac particle in field of two Coulomb centers
 */
class TwoCentreDiracEquation {
public:
    /**
     * @brief Two-center potential
     *
     * V(r) = -Z₁e²/r₁ - Z₂e²/r₂
     */
    static std::string potential() {
        return "V = -Z₁e²/(4πε₀r₁) - Z₂e²/(4πε₀r₂)";
    }

    /**
     * @brief Elliptic coordinates
     *
     * ξ = (r₁ + r₂)/R, η = (r₁ - r₂)/R
     */
    static std::pair<double, double> elliptic_coordinates(
        double r1,
        double r2,
        double R) {

        double xi = (r1 + r2) / R;
        double eta = (r1 - r2) / R;
        return {xi, eta};
    }

    /**
     * @brief Separation of variables
     *
     * Separable in elliptic coordinates
     */
    static std::string separation() {
        return "Separable in elliptic coordinates (ξ, η, φ)";
    }

    /**
     * @brief Application
     *
     * Molecular ions H₂⁺, molecular bonding
     */
    static std::string application() {
        return "Models H₂⁺ molecular ion and relativistic bonding";
    }

    /**
     * @brief Critical distance
     *
     * R_c determines bonding/antibonding
     */
    static std::string critical_distance() {
        return "Bonding/antibonding depends on internuclear distance R";
    }
};

// ============================================================================
// FOLDY-WOUTHUYSEN REPRESENTATION
// ============================================================================

/**
 * @brief Foldy-Wouthuysen Transformation - Free Particles
 *
 * Separate positive and negative energy states
 */
class FoldyWouthuysenFree {
public:
    /**
     * @brief FW transformation
     *
     * Unitary transformation to decouple large/small components
     */
    static std::string transformation() {
        return "U_FW = exp(iβα⃗·p⃗ tanh⁻¹(p/E)/2p)";
    }

    /**
     * @brief FW Hamiltonian for free particle
     *
     * H_FW = βE_p (diagonal in energy)
     */
    static std::string fw_hamiltonian() {
        return "H_FW = βsqrt(p² + m²) (energy eigenvalue operator)";
    }

    /**
     * @brief Advantage of FW representation
     *
     * Positive/negative energies completely separated
     */
    static std::string advantage() {
        return "E > 0 and E < 0 completely decoupled, no Zitterbewegung";
    }

    /**
     * @brief Position operator in FW
     *
     * r_FW = r + O(1/m) corrections
     */
    static std::string position_operator() {
        return "r_FW = r - iβΣ⃗×α⃗/(2E) (Newton-Wigner position)";
    }

    /**
     * @brief Velocity operator in FW
     *
     * v_FW = cp⃗/E
     */
    static double velocity_fw(double px, double E) {
        return (constants::c * px) / E;
    }
};

/**
 * @brief Foldy-Wouthuysen with External Fields
 *
 * Systematic expansion in powers of 1/m
 */
class FoldyWouthuysenFields {
public:
    /**
     * @brief FW Hamiltonian with EM field
     *
     * Systematic 1/m expansion
     */
    static std::string fw_expansion() {
        return "H_FW = βm + p²/2m + eφ + corrections O(1/m, 1/m², ...)";
    }

    /**
     * @brief Leading order (Pauli)
     *
     * H⁽⁰⁾ = p²/2m + eφ - (e/2m)σ⃗·B⃗
     */
    static std::string leading_order() {
        return "O(1): Pauli equation with g = 2";
    }

    /**
     * @brief First order corrections
     *
     * Spin-orbit, Darwin, relativistic kinetic
     */
    static std::string first_order() {
        return "O(1/m): H_SO + H_Darwin + H_kin";
    }

    /**
     * @brief Spin-orbit term
     *
     * H_SO = (e/2m²c²)(E⃗ × p⃗)·σ⃗ = (1/2m²c²r)(dV/dr)L⃗·S⃗
     */
    static double spin_orbit_fw(
        double r,
        double dV_dr,
        double L_dot_S,
        double mass) {

        return (1.0 / (2.0 * mass * mass * constants::c * constants::c * r)) *
               dV_dr * L_dot_S;
    }

    /**
     * @brief Darwin term
     *
     * H_Darwin = (ℏ²e/8m²c²)∇·E⃗ = (πℏ²e/2m²c²)|ψ(0)|²
     */
    static double darwin_fw(double field_divergence, double mass) {
        return (constants::e * constants::hbar * constants::hbar) /
               (8.0 * mass * mass * constants::c * constants::c) * field_divergence;
    }

    /**
     * @brief Advantages of FW
     *
     * Systematic perturbation theory, clear physical interpretation
     */
    static std::string advantages() {
        return "Clear separation of relativistic corrections, no Zitterbewegung";
    }
};

// ============================================================================
// HOLE THEORY AND SYMMETRIES
// ============================================================================

/**
 * @brief Hole Theory (Dirac Sea)
 *
 * Vacuum as filled Fermi sea of negative energy states
 */
class HoleTheory {
public:
    /**
     * @brief Dirac sea vacuum
     *
     * All E < 0 states filled
     */
    static std::string vacuum_state() {
        return "|0⟩ = |all E < 0 filled⟩ (Dirac sea)";
    }

    /**
     * @brief Hole interpretation
     *
     * Absence of E < 0 electron = positron with E > 0
     */
    static std::string hole_interpretation() {
        return "Hole in sea = antiparticle (positron)";
    }

    /**
     * @brief Charge of hole
     *
     * Missing electron → positive charge
     */
    static double hole_charge() {
        return +constants::e;  // Positive
    }

    /**
     * @brief Energy of hole
     *
     * E_hole = -E_electron > 0
     */
    static double hole_energy(double E_electron_negative) {
        return -E_electron_negative;  // Positive
    }

    /**
     * @brief Pair creation
     *
     * γ → e⁺ + e⁻ if ℏω ≥ 2mc²
     */
    static double pair_threshold(double mass) {
        return 2.0 * mass * constants::c * constants::c;
    }

    /**
     * @brief Pair annihilation
     *
     * e⁺ + e⁻ → 2γ
     */
    static std::string annihilation() {
        return "Positron fills hole → 2 photons (minimum)";
    }

    /**
     * @brief Vacuum polarization
     *
     * Virtual e⁺e⁻ pairs modify Coulomb potential
     */
    static std::string vacuum_polarization() {
        return "Virtual pairs screen charge: V_eff = V(1 + α/15π...)";
    }
};

/**
 * @brief Charge Conjugation
 *
 * Particle ↔ antiparticle symmetry
 */
class DiracChargeConjugation {
public:
    using DiracSpinor = std::array<Complex, 4>;

    /**
     * @brief Charge conjugation operator
     *
     * C: ψ → ψ^C = Cψ̄^T where C = iγ²γ⁰
     */
    static std::string charge_operator() {
        return "C = iγ²γ⁰ (charge conjugation matrix)";
    }

    /**
     * @brief Properties of C
     *
     * C^† = C^{-1} = -C, C^T = -C
     */
    static std::string c_properties() {
        return "C† = C⁻¹ = -C, C^T = -C, C² = -I";
    }

    /**
     * @brief Transformation of Dirac equation
     *
     * C transforms particle → antiparticle solution
     */
    static std::string transformation() {
        return "If ψ satisfies Dirac, so does ψ^C = Cψ̄^T";
    }

    /**
     * @brief Charge conjugation of plane waves
     *
     * u(p,s) ↔ v(p,s)
     */
    static std::string plane_wave_conjugation() {
        return "C: u(p,s)e^(-iEt+ip·x) → v(p,s)e^(iEt-ip·x)";
    }

    /**
     * @brief C-parity
     *
     * Eigenvalue of C for self-conjugate states
     */
    static int c_parity(bool even) {
        return even ? +1 : -1;
    }

    /**
     * @brief Majorana fermions
     *
     * Particles that are their own antiparticles
     */
    static std::string majorana() {
        return "Majorana: ψ = ψ^C (particle = antiparticle)";
    }
};

/**
 * @brief Charge Conjugation of Bound States
 *
 * C operation on hydrogen-like atoms
 */
class ChargeConjugationBoundStates {
public:
    /**
     * @brief Hydrogen atom charge conjugation
     *
     * C: (e⁻ + p) → (e⁺ + p̄) = antihydrogen
     */
    static std::string hydrogen_conjugation() {
        return "C: hydrogen → antihydrogen";
    }

    /**
     * @brief Energy levels under C
     *
     * E_n unchanged (CPT ensures H = H̄)
     */
    static std::string energy_conservation() {
        return "E_n(H) = E_n(H̄) by CPT theorem";
    }

    /**
     * @brief Wave function transformation
     *
     * ψ_nlm → ψ^C_nlm
     */
    static std::string wavefunction_c() {
        return "Cψ_nlm = ηψ_nlm for C-eigenstate";
    }

    /**
     * @brief Selection rules
     *
     * C-parity conservation in EM transitions
     */
    static std::string selection_rules() {
        return "ΔC = 0 for EM processes (if both initial/final have C)";
    }
};

/**
 * @brief Time Reversal Symmetry
 *
 * Motion reversal transformation
 */
class TimeReversalSymmetry {
public:
    /**
     * @brief Time reversal operator
     *
     * T: t → -t, ψ → Tψ
     */
    static std::string time_operator() {
        return "T = iγ¹γ³K where K is complex conjugation";
    }

    /**
     * @brief Antiunitary nature
     *
     * T is antiunitary (not unitary)
     */
    static std::string antiunitary() {
        return "T(aψ + bφ) = a*Tψ + b*Tφ (antiunitary)";
    }

    /**
     * @brief Transformation properties
     *
     * T: p⃗ → -p⃗, S⃗ → -S⃗, E → E, t → -t
     */
    static std::string transformations() {
        return "T: (t, r⃗, p⃗, S⃗) → (-t, r⃗, -p⃗, -S⃗)";
    }

    /**
     * @brief Kramers degeneracy
     *
     * T²|ψ⟩ = -|ψ⟩ for spin-1/2 → double degeneracy
     */
    static std::string kramers() {
        return "T² = -1 for fermions → Kramers degeneracy";
    }

    /**
     * @brief T-violation
     *
     * Weak interactions violate T (observed in K⁰ system)
     */
    static std::string violation() {
        return "T violated in weak interactions (K⁰, B⁰ decays)";
    }
};

/**
 * @brief PCT Theorem
 *
 * Combined CPT symmetry fundamental theorem
 */
class PCTSymmetry {
public:
    /**
     * @brief PCT transformation
     *
     * Combined parity, charge conjugation, time reversal
     */
    static std::string pct_operation() {
        return "PCT: (t,r⃗,ψ,Q) → (-t,-r⃗,ψ^C,-Q)";
    }

    /**
     * @brief PCT theorem
     *
     * All Lorentz invariant QFTs preserve CPT
     */
    static std::string pct_theorem() {
        return "CPT is exact symmetry of all local Lorentz-invariant QFTs";
    }

    /**
     * @brief Consequences
     *
     * m_particle = m_antiparticle, τ_particle = τ_antiparticle
     */
    static std::string consequences() {
        return "m_p = m_p̄, τ_p = τ_p̄, |q_p| = |q_p̄|";
    }

    /**
     * @brief CPT and causality
     *
     * Violation would break causality/locality
     */
    static std::string causality() {
        return "CPT violation → causality breakdown";
    }

    /**
     * @brief Experimental tests
     *
     * Precise tests using K⁰, B⁰, antihydrogen
     */
    static std::string tests() {
        return "Tested to <10⁻¹⁸ in K⁰ system, <10⁻⁸ in H̄";
    }

    /**
     * @brief Connection to spin-statistics
     *
     * PCT + Lorentz invariance → spin-statistics theorem
     */
    static std::string spin_statistics() {
        return "CPT + Lorentz → integer spin bosons, half-integer fermions";
    }
};

} // namespace relativistic_quantum
} // namespace physics

#endif // PHYSICS_RELATIVISTIC_QUANTUM_MECHANICS_HPP
