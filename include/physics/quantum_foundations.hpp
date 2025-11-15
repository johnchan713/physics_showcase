#ifndef PHYSICS_QUANTUM_FOUNDATIONS_HPP
#define PHYSICS_QUANTUM_FOUNDATIONS_HPP

#include <cmath>
#include <complex>
#include <vector>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <numeric>

/**
 * @file quantum_foundations.hpp
 * @brief Foundational quantum mechanics and historical development
 *
 * Comprehensive implementation of:
 * - The Discovery of Quantum Mechanics
 * - Planck and Quantization
 * - Bohr and the Hydrogen Atom
 * - Matrix Mechanics (Heisenberg)
 * - The Uncertainty Relations
 * - Wave Mechanics (Schrödinger)
 */

namespace physics {
namespace quantum_foundations {

// Physical constants
namespace constants {
    constexpr double h = 6.62607015e-34;        // Planck constant (J·s)
    constexpr double hbar = 1.054571817e-34;    // Reduced Planck constant (J·s)
    constexpr double c = 299792458.0;           // Speed of light (m/s)
    constexpr double k_B = 1.380649e-23;        // Boltzmann constant (J/K)
    constexpr double e = 1.602176634e-19;       // Elementary charge (C)
    constexpr double m_e = 9.1093837015e-31;    // Electron mass (kg)
    constexpr double epsilon_0 = 8.8541878128e-12; // Vacuum permittivity (F/m)
    constexpr double a_0 = 5.29177210903e-11;   // Bohr radius (m)
    constexpr double R_inf = 10973731.568160;   // Rydberg constant (m^-1)
}

using Complex = std::complex<double>;

/**
 * @brief Introduction to Quantum Mechanics
 *
 * Classical physics breakdown and need for quantum theory
 */
class IntroductionToQuantumMechanics {
public:
    /**
     * @brief Classical physics predictions that fail
     */

    /**
     * @brief Ultraviolet catastrophe: Rayleigh-Jeans law diverges
     *
     * Classical: u(ν,T) = (8πν²/c³)kT → ∞ as ν → ∞
     */
    static double rayleigh_jeans_law(double frequency, double temperature) {
        // Energy density per unit frequency (J·s/m³)
        return (8.0 * M_PI * frequency * frequency /
                (constants::c * constants::c * constants::c)) *
               constants::k_B * temperature;
    }

    /**
     * @brief Wien's displacement law: λ_max T = constant
     */
    static double wien_displacement_law(double temperature) {
        // Peak wavelength (m)
        constexpr double b = 2.897771955e-3; // Wien's constant (m·K)
        return b / temperature;
    }

    /**
     * @brief Stefan-Boltzmann law: total radiated power
     *
     * P = σAT⁴
     */
    static double stefan_boltzmann_law(double area, double temperature) {
        constexpr double sigma = 5.670374419e-8; // Stefan-Boltzmann constant (W·m⁻²·K⁻⁴)
        return sigma * area * std::pow(temperature, 4);
    }

    /**
     * @brief Photoelectric effect classical prediction (fails)
     *
     * Classical: kinetic energy should depend on intensity, not frequency
     */
    static double classical_photoelectric_prediction(double intensity) {
        // Classical (wrong) prediction: KE proportional to intensity
        return intensity;  // This fails experimentally!
    }

    /**
     * @brief Stability of atoms
     *
     * Classical: orbiting electron should radiate and spiral into nucleus
     */
    static double larmor_radiation_power(double acceleration, double charge) {
        // Larmor formula: P = (q²a²)/(6πε₀c³)
        return (charge * charge * acceleration * acceleration) /
               (6.0 * M_PI * constants::epsilon_0 *
                constants::c * constants::c * constants::c);
    }

    /**
     * @brief Classical orbital decay time (catastrophic)
     */
    static double classical_orbital_decay_time(double radius, double velocity) {
        // Rough estimate of collapse time (nanoseconds)
        // Classical prediction: atom collapses in ~10^-11 seconds

        double omega = velocity / radius;
        double acceleration = velocity * velocity / radius;
        double power = larmor_radiation_power(acceleration, constants::e);

        // Energy = (1/2)mv²
        double energy = 0.5 * constants::m_e * velocity * velocity;

        return energy / power;  // Collapse time
    }

    /**
     * @brief Specific heat of solids
     *
     * Classical Dulong-Petit law: C_V = 3Nk_B (independent of T)
     * Fails at low temperatures
     */
    static double dulong_petit_law(int n_atoms) {
        // Molar heat capacity (J/(mol·K))
        return 3.0 * n_atoms * constants::k_B;
    }
};

/**
 * @brief Planck and Quantization
 *
 * Birth of quantum theory (1900)
 */
class PlanckQuantization {
public:
    /**
     * @brief Planck's blackbody radiation formula
     *
     * u(ν,T) = (8πhν³/c³) / (e^(hν/kT) - 1)
     */
    static double planck_radiation_law(double frequency, double temperature) {
        double x = (constants::h * frequency) / (constants::k_B * temperature);

        if (x > 100.0) {
            // Wien limit (high frequency)
            return (8.0 * M_PI * constants::h *
                   std::pow(frequency, 3) /
                   (constants::c * constants::c * constants::c)) *
                   std::exp(-x);
        }

        // Full Planck formula
        return (8.0 * M_PI * constants::h * std::pow(frequency, 3) /
                (constants::c * constants::c * constants::c)) /
               (std::exp(x) - 1.0);
    }

    /**
     * @brief Planck's quantum hypothesis: E = nhν
     */
    static double energy_quantum(int n, double frequency) {
        return n * constants::h * frequency;
    }

    /**
     * @brief Average energy of oscillator at temperature T
     *
     * ⟨E⟩ = hν / (e^(hν/kT) - 1)
     */
    static double average_oscillator_energy(double frequency, double temperature) {
        double x = (constants::h * frequency) / (constants::k_B * temperature);
        return constants::h * frequency / (std::exp(x) - 1.0);
    }

    /**
     * @brief Partition function for quantum harmonic oscillator
     *
     * Z = 1 / (1 - e^(-hν/kT))
     */
    static double partition_function(double frequency, double temperature) {
        double x = (constants::h * frequency) / (constants::k_B * temperature);
        return 1.0 / (1.0 - std::exp(-x));
    }

    /**
     * @brief Einstein model for specific heat
     *
     * C_V = 3Nk_B (ℏω/kT)² e^(ℏω/kT) / (e^(ℏω/kT) - 1)²
     */
    static double einstein_specific_heat(
        double frequency,
        double temperature,
        int n_atoms) {

        double x = (constants::h * frequency) / (constants::k_B * temperature);
        double exp_x = std::exp(x);

        return 3.0 * n_atoms * constants::k_B *
               (x * x * exp_x) / ((exp_x - 1.0) * (exp_x - 1.0));
    }

    /**
     * @brief Debye model for specific heat (improved)
     *
     * More accurate at low temperatures
     */
    static double debye_specific_heat(
        double temperature,
        double debye_temperature,
        int n_atoms) {

        double x = debye_temperature / temperature;

        if (temperature < debye_temperature / 10.0) {
            // Low temperature limit: C_V ∝ T³
            return (12.0 * M_PI * M_PI * M_PI * M_PI / 5.0) *
                   n_atoms * constants::k_B *
                   std::pow(temperature / debye_temperature, 3);
        } else if (temperature > debye_temperature) {
            // High temperature limit: Dulong-Petit
            return 3.0 * n_atoms * constants::k_B;
        }

        // Intermediate (full Debye integral)
        // Simplified approximation
        return 3.0 * n_atoms * constants::k_B *
               (1.0 - std::exp(-3.0 * temperature / debye_temperature));
    }

    /**
     * @brief Photoelectric effect (Einstein 1905)
     *
     * E_kinetic = hν - W (work function)
     */
    static double photoelectric_kinetic_energy(
        double frequency,
        double work_function) {

        double photon_energy = constants::h * frequency;

        if (photon_energy < work_function) {
            return 0.0;  // No photoelectron emission
        }

        return photon_energy - work_function;
    }

    /**
     * @brief Threshold frequency for photoelectric effect
     */
    static double photoelectric_threshold(double work_function) {
        return work_function / constants::h;
    }

    /**
     * @brief Stopping potential in photoelectric effect
     */
    static double stopping_potential(double frequency, double work_function) {
        double ke = photoelectric_kinetic_energy(frequency, work_function);
        return ke / constants::e;  // Convert to voltage
    }
};

/**
 * @brief Bohr Model and the Hydrogen Atom
 *
 * Bohr's quantum theory of hydrogen (1913)
 */
class BohrHydrogenAtom {
public:
    /**
     * @brief Bohr radius: a₀ = ε₀h²/(πm_e e²)
     */
    static double bohr_radius() {
        return constants::a_0;
    }

    /**
     * @brief Orbital radius for quantum number n
     *
     * r_n = n² a₀
     */
    static double orbital_radius(int n) {
        if (n < 1) {
            throw std::invalid_argument("n must be ≥ 1");
        }
        return n * n * constants::a_0;
    }

    /**
     * @brief Energy levels: E_n = -13.6 eV / n²
     */
    static double energy_level(int n) {
        if (n < 1) {
            throw std::invalid_argument("n must be ≥ 1");
        }

        constexpr double E_1 = -13.605693122994; // Ground state energy (eV)
        return E_1 / (n * n);
    }

    /**
     * @brief Energy level in Joules
     */
    static double energy_level_joules(int n) {
        return energy_level(n) * constants::e;  // Convert eV to J
    }

    /**
     * @brief Orbital velocity: v_n = e²/(2ε₀hn) = c/(137n)
     */
    static double orbital_velocity(int n) {
        if (n < 1) {
            throw std::invalid_argument("n must be ≥ 1");
        }

        constexpr double alpha = 7.2973525693e-3; // Fine structure constant
        return alpha * constants::c / n;
    }

    /**
     * @brief Angular momentum quantization: L = nℏ
     */
    static double angular_momentum(int n) {
        if (n < 1) {
            throw std::invalid_argument("n must be ≥ 1");
        }
        return n * constants::hbar;
    }

    /**
     * @brief Rydberg formula for spectral lines
     *
     * 1/λ = R_∞ (1/n₁² - 1/n₂²)
     */
    static double rydberg_wavelength(int n1, int n2) {
        if (n1 < 1 || n2 <= n1) {
            throw std::invalid_argument("Need n₂ > n₁ ≥ 1");
        }

        double inv_lambda = constants::R_inf *
                           (1.0 / (n1 * n1) - 1.0 / (n2 * n2));

        return 1.0 / inv_lambda;  // Wavelength in meters
    }

    /**
     * @brief Photon energy for transition n₂ → n₁
     */
    static double transition_energy(int n1, int n2) {
        return energy_level(n2) - energy_level(n1);
    }

    /**
     * @brief Photon frequency for transition
     */
    static double transition_frequency(int n1, int n2) {
        double energy_diff = std::abs(transition_energy(n1, n2)) * constants::e;
        return energy_diff / constants::h;
    }

    /**
     * @brief Balmer series (visible): n₁ = 2
     */
    static double balmer_wavelength(int n) {
        if (n <= 2) {
            throw std::invalid_argument("n must be > 2 for Balmer series");
        }
        return rydberg_wavelength(2, n);
    }

    /**
     * @brief Lyman series (UV): n₁ = 1
     */
    static double lyman_wavelength(int n) {
        if (n <= 1) {
            throw std::invalid_argument("n must be > 1 for Lyman series");
        }
        return rydberg_wavelength(1, n);
    }

    /**
     * @brief Paschen series (IR): n₁ = 3
     */
    static double paschen_wavelength(int n) {
        if (n <= 3) {
            throw std::invalid_argument("n must be > 3 for Paschen series");
        }
        return rydberg_wavelength(3, n);
    }

    /**
     * @brief Ionization energy from level n
     */
    static double ionization_energy(int n) {
        return -energy_level(n);  // Energy to remove electron
    }

    /**
     * @brief De Broglie wavelength of electron in orbit n
     *
     * λ = h/p = 2πr_n/n
     */
    static double de_broglie_wavelength(int n) {
        double r_n = orbital_radius(n);
        return 2.0 * M_PI * r_n / n;
    }

    /**
     * @brief Bohr correspondence principle
     *
     * For large n, quantum → classical
     */
    static double classical_frequency(int n) {
        // Classical orbital frequency
        double v = orbital_velocity(n);
        double r = orbital_radius(n);
        return v / (2.0 * M_PI * r);
    }
};

/**
 * @brief Matrix Mechanics (Heisenberg 1925)
 *
 * First formulation of quantum mechanics
 */
class MatrixMechanics {
public:
    using Matrix = std::vector<std::vector<Complex>>;

    /**
     * @brief Position matrix elements for harmonic oscillator
     *
     * x_nm = √(ℏ/2mω) (√n δ_{n,m+1} + √(n+1) δ_{n,m-1})
     */
    static Matrix position_matrix(int n_max, double mass, double omega) {
        Matrix x(n_max, std::vector<Complex>(n_max, 0.0));

        double coeff = std::sqrt(constants::hbar / (2.0 * mass * omega));

        for (int n = 0; n < n_max; ++n) {
            if (n > 0) {
                x[n][n-1] = coeff * std::sqrt(n);
            }
            if (n < n_max - 1) {
                x[n][n+1] = coeff * std::sqrt(n + 1);
            }
        }

        return x;
    }

    /**
     * @brief Momentum matrix elements for harmonic oscillator
     *
     * p_nm = i√(ℏmω/2) (√n δ_{n,m+1} - √(n+1) δ_{n,m-1})
     */
    static Matrix momentum_matrix(int n_max, double mass, double omega) {
        Matrix p(n_max, std::vector<Complex>(n_max, 0.0));

        double coeff = std::sqrt(constants::hbar * mass * omega / 2.0);
        Complex i(0.0, 1.0);

        for (int n = 0; n < n_max; ++n) {
            if (n > 0) {
                p[n][n-1] = i * coeff * std::sqrt(n);
            }
            if (n < n_max - 1) {
                p[n][n+1] = -i * coeff * std::sqrt(n + 1);
            }
        }

        return p;
    }

    /**
     * @brief Hamiltonian matrix for harmonic oscillator
     *
     * H = p²/2m + (1/2)mω²x²
     */
    static Matrix hamiltonian_harmonic_oscillator(int n_max, double mass, double omega) {
        Matrix H(n_max, std::vector<Complex>(n_max, 0.0));

        // Diagonal: E_n = ℏω(n + 1/2)
        for (int n = 0; n < n_max; ++n) {
            H[n][n] = constants::hbar * omega * (n + 0.5);
        }

        return H;
    }

    /**
     * @brief Commutator [A, B] = AB - BA
     */
    static Matrix commutator(const Matrix& A, const Matrix& B) {
        int n = A.size();
        Matrix result(n, std::vector<Complex>(n, 0.0));

        // Compute AB
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        // Subtract BA
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    result[i][j] -= B[i][k] * A[k][j];
                }
            }
        }

        return result;
    }

    /**
     * @brief Canonical commutation relation: [x, p] = iℏ
     */
    static bool verify_canonical_commutation(
        const Matrix& x,
        const Matrix& p,
        double tol = 1e-10) {

        auto comm = commutator(x, p);
        int n = x.size();

        Complex ihbar(0.0, constants::hbar);

        // Check if [x, p] = iℏ I
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Complex expected = (i == j) ? ihbar : Complex(0.0, 0.0);
                if (std::abs(comm[i][j] - expected) > tol) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * @brief Ladder operators for harmonic oscillator
     *
     * a = √(mω/2ℏ) x + i√(1/2mωℏ) p
     * a† = √(mω/2ℏ) x - i√(1/2mωℏ) p
     */
    static std::pair<Matrix, Matrix> ladder_operators(
        int n_max,
        double mass,
        double omega) {

        // Annihilation operator a
        Matrix a(n_max, std::vector<Complex>(n_max, 0.0));

        for (int n = 0; n < n_max - 1; ++n) {
            a[n][n+1] = std::sqrt(n + 1);
        }

        // Creation operator a†
        Matrix a_dagger(n_max, std::vector<Complex>(n_max, 0.0));

        for (int n = 1; n < n_max; ++n) {
            a_dagger[n][n-1] = std::sqrt(n);
        }

        return {a, a_dagger};
    }

    /**
     * @brief Number operator N = a†a
     */
    static Matrix number_operator(int n_max) {
        Matrix N(n_max, std::vector<Complex>(n_max, 0.0));

        for (int n = 0; n < n_max; ++n) {
            N[n][n] = n;
        }

        return N;
    }

    /**
     * @brief Heisenberg equation of motion: dA/dt = (i/ℏ)[H, A]
     */
    static Matrix heisenberg_equation_of_motion(
        const Matrix& H,
        const Matrix& A) {

        auto comm = commutator(H, A);
        int n = H.size();

        Matrix dA_dt(n, std::vector<Complex>(n));
        Complex factor(0.0, 1.0 / constants::hbar);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                dA_dt[i][j] = factor * comm[i][j];
            }
        }

        return dA_dt;
    }
};

/**
 * @brief The Uncertainty Relations (Heisenberg 1927)
 *
 * Fundamental limits on simultaneous measurements
 */
class UncertaintyRelations {
public:
    /**
     * @brief Heisenberg uncertainty principle: ΔxΔp ≥ ℏ/2
     */
    static double heisenberg_uncertainty_minimum() {
        return constants::hbar / 2.0;
    }

    /**
     * @brief Position uncertainty for harmonic oscillator ground state
     */
    static double position_uncertainty_ground_state(double mass, double omega) {
        return std::sqrt(constants::hbar / (2.0 * mass * omega));
    }

    /**
     * @brief Momentum uncertainty for harmonic oscillator ground state
     */
    static double momentum_uncertainty_ground_state(double mass, double omega) {
        return std::sqrt(constants::hbar * mass * omega / 2.0);
    }

    /**
     * @brief Verify minimum uncertainty product for ground state
     */
    static bool verify_minimum_uncertainty(double mass, double omega) {
        double dx = position_uncertainty_ground_state(mass, omega);
        double dp = momentum_uncertainty_ground_state(mass, omega);

        double product = dx * dp;
        double minimum = heisenberg_uncertainty_minimum();

        return std::abs(product - minimum) < 1e-10;
    }

    /**
     * @brief Energy-time uncertainty: ΔEΔt ≥ ℏ/2
     */
    static double energy_time_uncertainty_minimum() {
        return constants::hbar / 2.0;
    }

    /**
     * @brief Minimum time for energy change ΔE
     */
    static double minimum_time_for_energy_change(double delta_E) {
        return constants::hbar / (2.0 * delta_E);
    }

    /**
     * @brief Natural line width from lifetime
     *
     * ΔE = ℏ/Δt
     */
    static double natural_line_width(double lifetime) {
        return constants::hbar / lifetime;
    }

    /**
     * @brief Robertson-Schrödinger uncertainty relation
     *
     * ΔAΔ B ≥ (1/2)|⟨[A,B]⟩|
     */
    static double robertson_uncertainty_bound(
        double expectation_commutator) {

        return 0.5 * std::abs(expectation_commutator);
    }

    /**
     * @brief Gaussian wave packet uncertainties
     */
    static std::pair<double, double> gaussian_packet_uncertainties(double sigma_x) {
        // For Gaussian: Δx = σ_x
        double delta_x = sigma_x;

        // Δp = ℏ/(2σ_x) for minimum uncertainty
        double delta_p = constants::hbar / (2.0 * sigma_x);

        return {delta_x, delta_p};
    }

    /**
     * @brief Wave packet spreading
     *
     * σ_x(t) = σ_x(0)√(1 + (ℏt/2mσ_x²)²)
     */
    static double wave_packet_width(
        double sigma_x_0,
        double mass,
        double time) {

        double factor = (constants::hbar * time) / (2.0 * mass * sigma_x_0 * sigma_x_0);
        return sigma_x_0 * std::sqrt(1.0 + factor * factor);
    }

    /**
     * @brief Coherent state (minimum uncertainty state)
     */
    static bool is_coherent_state(
        double delta_x,
        double delta_p) {

        double product = delta_x * delta_p;
        double minimum = heisenberg_uncertainty_minimum();

        return std::abs(product - minimum) < 1e-10;
    }

    /**
     * @brief Squeezed state (Δx < Δx_min or Δp < Δp_min)
     */
    static bool is_squeezed_state(
        double delta_x,
        double delta_p,
        double delta_x_coherent,
        double delta_p_coherent) {

        return (delta_x < delta_x_coherent) || (delta_p < delta_p_coherent);
    }
};

/**
 * @brief Wave Mechanics (Schrödinger 1926)
 *
 * Wave function formulation of quantum mechanics
 */
class WaveMechanics {
public:
    /**
     * @brief De Broglie relation: λ = h/p
     */
    static double de_broglie_wavelength(double momentum) {
        return constants::h / momentum;
    }

    /**
     * @brief De Broglie wavelength from velocity
     */
    static double de_broglie_wavelength_velocity(double mass, double velocity) {
        double momentum = mass * velocity;
        return de_broglie_wavelength(momentum);
    }

    /**
     * @brief Wave number: k = 2π/λ = p/ℏ
     */
    static double wave_number(double momentum) {
        return momentum / constants::hbar;
    }

    /**
     * @brief Angular frequency: ω = E/ℏ
     */
    static double angular_frequency(double energy) {
        return energy / constants::hbar;
    }

    /**
     * @brief Phase velocity: v_p = ω/k = E/p
     */
    static double phase_velocity(double energy, double momentum) {
        return energy / momentum;
    }

    /**
     * @brief Group velocity: v_g = dω/dk = ∂E/∂p
     *
     * For free particle: v_g = p/m (equals classical velocity)
     */
    static double group_velocity(double momentum, double mass) {
        return momentum / mass;
    }

    /**
     * @brief Free particle wave function: ψ(x,t) = A e^(i(kx - ωt))
     */
    static Complex free_particle_wave_function(
        double x,
        double t,
        double momentum,
        double energy,
        double amplitude = 1.0) {

        double k = wave_number(momentum);
        double omega = angular_frequency(energy);

        Complex i(0.0, 1.0);
        return amplitude * std::exp(i * (k * x - omega * t));
    }

    /**
     * @brief Gaussian wave packet
     *
     * ψ(x,0) = (1/(2πσ²))^(1/4) exp(-x²/4σ² + ik₀x)
     */
    static Complex gaussian_wave_packet(
        double x,
        double sigma,
        double k0) {

        double norm = std::pow(1.0 / (2.0 * M_PI * sigma * sigma), 0.25);
        Complex i(0.0, 1.0);

        return norm * std::exp(-x * x / (4.0 * sigma * sigma) + i * k0 * x);
    }

    /**
     * @brief Time-evolved Gaussian wave packet
     */
    static Complex gaussian_wave_packet_time_evolved(
        double x,
        double t,
        double sigma,
        double k0,
        double mass) {

        Complex i(0.0, 1.0);
        double hbar_over_m = constants::hbar / mass;

        // Spreading parameter
        Complex sigma_t = sigma * std::sqrt(Complex(1.0, hbar_t * t / (mass * sigma * sigma)));

        double norm = std::pow(1.0 / (2.0 * M_PI * sigma * sigma), 0.25);

        double x_shifted = x - (constants::hbar * k0 * t / mass);

        Complex phase = i * (k0 * x - (constants::hbar * k0 * k0 * t) / (2.0 * mass));

        return norm * std::exp(-x_shifted * x_shifted / (4.0 * sigma * sigma)) * std::exp(phase);
    }

    /**
     * @brief Harmonic oscillator ground state wave function
     *
     * ψ₀(x) = (mω/πℏ)^(1/4) exp(-mωx²/2ℏ)
     */
    static double harmonic_oscillator_ground_state(
        double x,
        double mass,
        double omega) {

        double norm = std::pow(mass * omega / (M_PI * constants::hbar), 0.25);
        double exponent = -(mass * omega * x * x) / (2.0 * constants::hbar);

        return norm * std::exp(exponent);
    }

    /**
     * @brief Hermite polynomials H_n(ξ)
     */
    static double hermite_polynomial(int n, double xi) {
        if (n == 0) return 1.0;
        if (n == 1) return 2.0 * xi;

        // Recursion: H_{n+1} = 2ξH_n - 2nH_{n-1}
        double H_prev = 1.0;
        double H_curr = 2.0 * xi;

        for (int k = 1; k < n; ++k) {
            double H_next = 2.0 * xi * H_curr - 2.0 * k * H_prev;
            H_prev = H_curr;
            H_curr = H_next;
        }

        return H_curr;
    }

    /**
     * @brief Harmonic oscillator eigenstate n
     *
     * ψ_n(x) = (1/√(2ⁿn!)) (mω/πℏ)^(1/4) H_n(√(mω/ℏ)x) exp(-mωx²/2ℏ)
     */
    static double harmonic_oscillator_eigenstate(
        int n,
        double x,
        double mass,
        double omega) {

        double xi = std::sqrt(mass * omega / constants::hbar) * x;
        double H_n = hermite_polynomial(n, xi);

        // Normalization
        double factorial_n = 1.0;
        for (int k = 1; k <= n; ++k) {
            factorial_n *= k;
        }

        double norm = std::pow(mass * omega / (M_PI * constants::hbar), 0.25) /
                     std::sqrt(std::pow(2.0, n) * factorial_n);

        double gaussian = std::exp(-(mass * omega * x * x) / (2.0 * constants::hbar));

        return norm * H_n * gaussian;
    }

    /**
     * @brief Probability density |ψ(x)|²
     */
    static double probability_density(Complex psi) {
        return std::norm(psi);  // |ψ|² = ψ*ψ
    }

    /**
     * @brief Probability current density
     *
     * j = (ℏ/2mi)(ψ*∇ψ - ψ∇ψ*)
     */
    static double probability_current_1d(
        Complex psi,
        Complex psi_prime) {

        // j = (ℏ/2mi)(ψ*ψ' - ψ(ψ*)')
        Complex i(0.0, 1.0);
        Complex current = (std::conj(psi) * psi_prime - psi * std::conj(psi_prime));

        return (constants::hbar / (2.0 * constants::m_e)) * current.imag();
    }

    /**
     * @brief Continuity equation: ∂ρ/∂t + ∇·j = 0
     */
    static bool verify_continuity_equation(
        double drho_dt,
        double div_j,
        double tol = 1e-10) {

        return std::abs(drho_dt + div_j) < tol;
    }

    /**
     * @brief Normalization integral ∫|ψ|²dx = 1
     */
    static double normalization_integral_numerical(
        const std::function<Complex(double)>& psi,
        double x_min,
        double x_max,
        int n_points = 1000) {

        double dx = (x_max - x_min) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double x = x_min + i * dx;
            integral += std::norm(psi(x)) * dx;
        }

        return integral;
    }

    /**
     * @brief Expectation value ⟨A⟩ = ∫ψ*Aψ dx
     */
    static double expectation_value_position(
        const std::function<double(double)>& psi_real,
        double x_min,
        double x_max,
        int n_points = 1000) {

        double dx = (x_max - x_min) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double x = x_min + i * dx;
            double psi_val = psi_real(x);
            integral += psi_val * x * psi_val * dx;
        }

        return integral;
    }

    /**
     * @brief Born rule: probability of measurement
     *
     * P(x ∈ [a,b]) = ∫_a^b |ψ(x)|² dx
     */
    static double born_rule_probability(
        const std::function<Complex(double)>& psi,
        double a,
        double b,
        int n_points = 1000) {

        return normalization_integral_numerical(psi, a, b, n_points);
    }
};

} // namespace quantum_foundations
} // namespace physics

#endif // PHYSICS_QUANTUM_FOUNDATIONS_HPP
