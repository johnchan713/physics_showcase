/**
 * @file pde_solution_methods_demo.cpp
 * @brief Demonstration of PDE solution methods, orthogonal series, and eigenfunction expansions
 *
 * Topics covered:
 * - Series of Orthogonal Functions (Parseval's identity, completeness, convergence)
 * - Eigenfunction Expansions (Sturm-Liouville theory)
 * - Orthogonal Polynomials
 * - Fourier Series
 * - Bessel Functions
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>
#include "maths/advanced/pde/pde_solution_methods.hpp"

using namespace maths::pde;

void print_section(const std::string& title) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(80, '=') << "\n";
}

void print_subsection(const std::string& title) {
    std::cout << "\n--- " << title << " ---\n";
}

/**
 * Demonstrate orthogonal series theory
 */
void demo_orthogonal_series_theory() {
    print_section("ORTHOGONAL SERIES THEORY");

    // Test function: f(x) = x on [-π, π]
    auto f = [](double x) { return x; };
    double L = M_PI;

    print_subsection("Fourier Series Expansion");

    // Compute Fourier series coefficients
    int n_terms = 20;
    auto coeffs = FourierSeries::computeCoefficients(f, L, n_terms);

    std::cout << "Function: f(x) = x on [-π, π]\n";
    std::cout << "Number of terms: " << n_terms << "\n\n";

    // Show first few coefficients
    std::cout << "Fourier coefficients:\n";
    std::cout << "a₀ = " << coeffs.a0 << "\n";
    for (int n = 0; n < 5; ++n) {
        std::cout << "a₁₀ = " << coeffs.a_n[n]
                  << ", b₁₀ = " << coeffs.b_n[n] << "\n";
    }

    print_subsection("Parseval's Identity");

    // Compute ∥f∥² directly
    double f_norm_sq_direct = OrthogonalFunctions::innerProduct(
        f, f, [](double) { return 1.0; }, -L, L);

    std::cout << "∥f∥² (direct integration) = " << f_norm_sq_direct << "\n";

    // Compute via Parseval: ∥f∥² = L(a₀²/2 + Σ(aₙ² + bₙ²))
    double parseval_sum = L * (coeffs.a0 * coeffs.a0 / 2.0);
    for (size_t n = 0; n < coeffs.a_n.size(); ++n) {
        parseval_sum += L * (coeffs.a_n[n] * coeffs.a_n[n] +
                             coeffs.b_n[n] * coeffs.b_n[n]);
    }

    std::cout << "∥f∥² (Parseval's identity) = " << parseval_sum << "\n";
    std::cout << "Relative error: "
              << std::abs(parseval_sum - f_norm_sq_direct) / f_norm_sq_direct * 100
              << "%\n";

    print_subsection("Convergence Rate Analysis");

    // Analyze decay of Fourier coefficients
    std::vector<double> all_coeffs;
    for (size_t n = 0; n < coeffs.b_n.size(); ++n) {
        all_coeffs.push_back(std::abs(coeffs.b_n[n]));
    }

    double decay_rate = OrthogonalSeriesTheory::convergenceRate(all_coeffs);
    std::cout << "Coefficient decay rate p: " << decay_rate << "\n";
    std::cout << "Expected: p ≈ 1 for piecewise smooth functions\n";

    print_subsection("Series Approximation Error");

    // Test with different numbers of terms
    std::vector<int> N_values = {5, 10, 20, 50};

    std::cout << "\nApproximation quality:\n";
    std::cout << std::setw(10) << "N"
              << std::setw(20) << "f(π/2)"
              << std::setw(20) << "Approximation"
              << std::setw(20) << "Error\n";
    std::cout << std::string(70, '-') << "\n";

    double x_test = M_PI / 2.0;
    double f_exact = f(x_test);

    for (int N : N_values) {
        auto coeffs_N = FourierSeries::computeCoefficients(f, L, N);
        double approx = FourierSeries::evaluate(coeffs_N, x_test);
        double error = std::abs(f_exact - approx);

        std::cout << std::setw(10) << N
                  << std::setw(20) << f_exact
                  << std::setw(20) << approx
                  << std::setw(20) << error << "\n";
    }
}

/**
 * Demonstrate eigenfunction expansions and Sturm-Liouville theory
 */
void demo_eigenfunction_expansions() {
    print_section("EIGENFUNCTION EXPANSIONS - STURM-LIOUVILLE THEORY");

    print_subsection("Problem 1: Fourier Sine Series (u'' + λu = 0)");

    double L = M_PI;
    auto sl_problem = EigenfunctionExpansions::fourierSineProblem(L);

    std::cout << "Boundary value problem:\n";
    std::cout << "  u'' + λu = 0 on [0, π]\n";
    std::cout << "  u(0) = u(π) = 0\n\n";

    // Analytical eigenvalues: λₙ = n²
    std::cout << "Analytical eigenvalues: λₙ = n²\n";
    std::cout << "λ₁ = 1, λ₂ = 4, λ₃ = 9, λ₄ = 16, λ₅ = 25\n\n";

    // Compute eigenvalues numerically
    std::cout << "Computing eigenvalues numerically...\n";
    auto eigenvalues = sl_problem.computeEigenvalues(5);

    std::cout << "Numerical eigenvalues:\n";
    for (size_t i = 0; i < eigenvalues.size(); ++i) {
        double exact = (i + 1) * (i + 1);
        std::cout << "λ₁₀ = " << eigenvalues[i]
                  << " (exact: " << exact << ")"
                  << " error: " << std::abs(eigenvalues[i] - exact) << "\n";
    }

    print_subsection("Eigenfunction Expansion");

    // Expand f(x) = x(π - x) in eigenfunctions
    auto f = [L](double x) { return x * (L - x); };

    std::cout << "\nExpanding f(x) = x(π - x) in sine series:\n";

    if (eigenvalues.size() >= 3) {
        std::vector<double> first_three_eigenvalues(eigenvalues.begin(), eigenvalues.begin() + 3);
        auto expansion_coeffs = sl_problem.expandFunction(f, first_three_eigenvalues);

        std::cout << "Expansion coefficients:\n";
        for (size_t i = 0; i < expansion_coeffs.size(); ++i) {
            std::cout << "c₁₀ = " << expansion_coeffs[i] << "\n";
        }

        // Reconstruct and evaluate at test point
        double x_test = L / 2.0;
        double f_exact = f(x_test);
        double f_approx = 0.0;

        for (size_t i = 0; i < expansion_coeffs.size(); ++i) {
            auto phi = sl_problem.eigenfunction(first_three_eigenvalues[i]);
            f_approx += expansion_coeffs[i] * phi(x_test);
        }

        std::cout << "\nReconstruction at x = π/2:\n";
        std::cout << "Exact: " << f_exact << "\n";
        std::cout << "Approximation (3 terms): " << f_approx << "\n";
        std::cout << "Error: " << std::abs(f_exact - f_approx) << "\n";
    }

    print_subsection("Problem 2: Bessel Problem (Cylindrical Coordinates)");

    double R = 1.0;
    auto bessel_problem = EigenfunctionExpansions::besselProblem(R);

    std::cout << "Boundary value problem:\n";
    std::cout << "  (xu')' + λxu = 0 on [0, 1]\n";
    std::cout << "  u(0) bounded, u(1) = 0\n\n";

    std::cout << "Eigenvalues correspond to zeros of J₀(x):\n";
    auto j0_zeros = BesselFunctions::besselJZeros(0, 5);

    std::cout << "Zeros of J₀(x):\n";
    for (size_t i = 0; i < j0_zeros.size(); ++i) {
        std::cout << "j₀,₁₀ = " << j0_zeros[i]
                  << ", λ₁₀ = " << j0_zeros[i] * j0_zeros[i] << "\n";
    }
}

/**
 * Demonstrate orthogonal polynomials
 */
void demo_orthogonal_polynomials() {
    print_section("ORTHOGONAL POLYNOMIALS");

    print_subsection("Legendre Polynomials Pₙ(x)");

    std::cout << "Legendre polynomials on [-1, 1] with weight w(x) = 1\n\n";

    // Evaluate at x = 0.5
    double x = 0.5;
    std::cout << "Values at x = " << x << ":\n";
    for (int n = 0; n <= 5; ++n) {
        double Pn = OrthogonalPolynomials::legendreP(n, x);
        std::cout << "P₁₀(" << x << ") = " << Pn << "\n";
    }

    // Verify orthogonality
    print_subsection("Orthogonality Verification");

    auto P2 = [](double x) { return OrthogonalPolynomials::legendreP(2, x); };
    auto P3 = [](double x) { return OrthogonalPolynomials::legendreP(3, x); };
    auto weight = [](double) { return 1.0; };

    double inner_product = OrthogonalFunctions::innerProduct(P2, P3, weight, -1.0, 1.0);

    std::cout << "⟨P₂, P₃⟩ = " << inner_product << " (should be ≈ 0)\n";
    std::cout << "Orthogonal: " << (std::abs(inner_product) < 1e-6 ? "YES" : "NO") << "\n";

    print_subsection("Chebyshev Polynomials Tₙ(x)");

    std::cout << "\nChebyshev polynomials on [-1, 1] with weight w(x) = 1/√(1-x²)\n\n";

    x = 0.7;
    std::cout << "Values at x = " << x << ":\n";
    for (int n = 0; n <= 5; ++n) {
        double Tn = OrthogonalPolynomials::chebyshevT(n, x);
        std::cout << "T₁₀(" << x << ") = " << Tn << "\n";
    }

    print_subsection("Hermite Polynomials Hₙ(x)");

    std::cout << "\nHermite polynomials on (-∞, ∞) with weight w(x) = e^(-x²)\n\n";

    x = 1.0;
    std::cout << "Values at x = " << x << ":\n";
    for (int n = 0; n <= 5; ++n) {
        double Hn = OrthogonalPolynomials::hermiteH(n, x);
        std::cout << "H₁₀(" << x << ") = " << Hn << "\n";
    }

    print_subsection("Laguerre Polynomials Lₙ(x)");

    std::cout << "\nLaguerre polynomials on [0, ∞) with weight w(x) = e^(-x)\n\n";

    x = 2.0;
    std::cout << "Values at x = " << x << ":\n";
    for (int n = 0; n <= 5; ++n) {
        double Ln = OrthogonalPolynomials::laguerreL(n, x);
        std::cout << "L₁₀(" << x << ") = " << Ln << "\n";
    }
}

/**
 * Demonstrate Bessel functions
 */
void demo_bessel_functions() {
    print_section("BESSEL FUNCTIONS");

    print_subsection("Bessel Function of First Kind Jₙ(x)");

    double x = 3.0;
    std::cout << "Values at x = " << x << ":\n";
    for (int n = 0; n <= 5; ++n) {
        double Jn = BesselFunctions::besselJ(n, x);
        std::cout << "J₁₀(" << x << ") = " << Jn << "\n";
    }

    print_subsection("Recurrence Relation Verification");

    // Verify: J_{n+1}(x) = (2n/x)J_n(x) - J_{n-1}(x)
    int n = 3;
    double J_n_minus_1 = BesselFunctions::besselJ(n - 1, x);
    double J_n = BesselFunctions::besselJ(n, x);
    double J_n_plus_1 = BesselFunctions::besselJ(n + 1, x);

    double J_recurrence = (2.0 * n / x) * J_n - J_n_minus_1;

    std::cout << "For n = " << n << ", x = " << x << ":\n";
    std::cout << "J₄(x) (direct) = " << J_n_plus_1 << "\n";
    std::cout << "J₄(x) (recurrence) = " << J_recurrence << "\n";
    std::cout << "Error: " << std::abs(J_n_plus_1 - J_recurrence) << "\n";

    print_subsection("Zeros of J₀(x)");

    auto zeros = BesselFunctions::besselJZeros(0, 5);

    std::cout << "First 5 zeros of J₀(x):\n";
    for (size_t i = 0; i < zeros.size(); ++i) {
        double zero = zeros[i];
        double J0_at_zero = BesselFunctions::besselJ(0, zero);
        std::cout << "j₀,₁₀ = " << zero
                  << ", J₀(j₀,₁₀) = " << J0_at_zero << "\n";
    }

    print_subsection("Modified Bessel Function Iₙ(x)");

    x = 2.0;
    std::cout << "\nValues at x = " << x << ":\n";
    for (int n = 0; n <= 3; ++n) {
        double In = BesselFunctions::besselI(n, x);
        std::cout << "I₁₀(" << x << ") = " << In << "\n";
    }
}

/**
 * Main demonstration
 */
int main() {
    std::cout << std::fixed << std::setprecision(8);

    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║        PDE SOLUTION METHODS - ORTHOGONAL SERIES & EIGENFUNCTIONS          ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════════════╝\n";

    try {
        demo_orthogonal_polynomials();
        demo_bessel_functions();
        demo_orthogonal_series_theory();
        demo_eigenfunction_expansions();

        std::cout << "\n" << std::string(80, '=') << "\n";
        std::cout << "All demonstrations completed successfully!\n";
        std::cout << std::string(80, '=') << "\n\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
