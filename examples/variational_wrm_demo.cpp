/**
 * @file variational_wrm_demo.cpp
 * @brief Demonstration of Weighted Residual Methods and variational techniques
 *
 * Topics covered:
 * - Green's Identities
 * - Weighted Residual Methods framework
 * - Galerkin Method
 * - Collocation Method
 * - Subdomain Method
 * - Least Squares Method
 * - Test Function Selection
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>
#include "maths/advanced/pde/pde_variational_methods.hpp"

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
 * Demonstrate Green's Identities
 */
void demo_greens_identities() {
    print_section("GREEN'S IDENTITIES");

    // Domain: [0, 1] × [0, 1]
    double x_start = 0.0, x_end = 1.0;
    double y_start = 0.0, y_end = 1.0;

    // Test functions
    auto u = [](double x, double y) { return x * x + y * y; };
    auto v = [](double x, double y) { return std::sin(M_PI * x) * std::sin(M_PI * y); };

    // Gradients
    auto grad_u_x = [](double x, double) { return 2.0 * x; };
    auto grad_u_y = [](double, double y) { return 2.0 * y; };
    auto grad_v_x = [](double x, double y) { return M_PI * std::cos(M_PI * x) * std::sin(M_PI * y); };
    auto grad_v_y = [](double x, double y) { return M_PI * std::sin(M_PI * x) * std::cos(M_PI * y); };

    // Laplacians
    auto laplacian_u = [](double, double) { return 4.0; };
    auto laplacian_v = [](double x, double y) {
        return -2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
    };

    print_subsection("Green's First Identity");
    std::cout << "∫∫ (∇u·∇v + v∇²u) dA = ∫ v(∂u/∂n) ds\n\n";

    double first_identity = GreensIdentities::greensFirstIdentity(
        grad_u_x, grad_u_y, grad_v_x, grad_v_y, v, laplacian_u,
        x_start, x_end, y_start, y_end);

    std::cout << "Volume integral: " << first_identity << "\n";

    print_subsection("Green's Second Identity");
    std::cout << "∫∫ (v∇²u - u∇²v) dA = ∫ (v∂u/∂n - u∂v/∂n) ds\n\n";

    double second_identity = GreensIdentities::greensSecondIdentity(
        u, v, laplacian_u, laplacian_v, grad_u_x, grad_u_y, grad_v_x, grad_v_y,
        x_start, x_end, y_start, y_end);

    std::cout << "Result: " << second_identity << "\n";

    print_subsection("Integration by Parts");

    auto f = [](double x) { return x * x; };
    auto g_prime = [](double x) { return std::cos(x); };
    auto f_prime = [](double x) { return 2.0 * x; };
    auto g = [](double x) { return std::sin(x); };

    double ibp = GreensIdentities::integrationByParts(
        f, g_prime, f_prime, g, 0.0, M_PI);

    std::cout << "∫₀^π x² cos(x) dx = [x² sin(x)]₀^π - ∫₀^π 2x sin(x) dx\n";
    std::cout << "Result: " << ibp << "\n";
}

/**
 * Demonstrate Weighted Residual Methods framework
 */
void demo_weighted_residual_methods() {
    print_section("WEIGHTED RESIDUAL METHODS");

    // Example problem: -u'' = 1 on [0, 1], u(0) = u(1) = 0
    // Approximate solution: u ≈ c₁ x(1-x) + c₂ x²(1-x)

    print_subsection("Problem Setup");
    std::cout << "BVP: -u'' = 1 on [0, 1]\n";
    std::cout << "     u(0) = u(1) = 0\n";
    std::cout << "Exact solution: u(x) = x(1-x)/2\n\n";

    auto exact_solution = [](double x) { return 0.5 * x * (1.0 - x); };

    // Differential operator: L[u] = -u''
    auto L_operator = [](std::function<double(double)> u, double x, double h = 1e-5) {
        // Numerical second derivative
        return -(u(x + h) - 2.0 * u(x) + u(x - h)) / (h * h);
    };

    // Source term
    auto f = [](double) { return 1.0; };

    print_subsection("Trial Functions");

    // Trial functions (satisfy boundary conditions)
    std::vector<std::function<double(double)>> trial_functions = {
        [](double x) { return x * (1.0 - x); },           // φ₁
        [](double x) { return x * x * (1.0 - x); }        // φ₂
    };

    std::cout << "φ₁(x) = x(1-x)\n";
    std::cout << "φ₂(x) = x²(1-x)\n\n";

    // Test with u ≈ c₁φ₁ + c₂φ₂
    std::vector<double> coeffs = {0.5, 0.0};  // Start with simple approximation

    auto approx_solution = [&](double x) {
        double result = 0.0;
        for (size_t i = 0; i < coeffs.size(); ++i) {
            result += coeffs[i] * trial_functions[i](x);
        }
        return result;
    };

    print_subsection("Residual Computation");

    // Compute residual R = L[u_approx] - f
    std::vector<double> x_test = {0.25, 0.5, 0.75};
    std::cout << std::setw(10) << "x"
              << std::setw(20) << "L[u]"
              << std::setw(20) << "f(x)"
              << std::setw(20) << "Residual\n";
    std::cout << std::string(70, '-') << "\n";

    for (double x : x_test) {
        double Lu = L_operator(approx_solution, x);
        double fx = f(x);
        double residual = WeightedResidualMethods::computeResidual(
            L_operator, approx_solution, f, x);

        std::cout << std::setw(10) << x
                  << std::setw(20) << Lu
                  << std::setw(20) << fx
                  << std::setw(20) << residual << "\n";
    }

    print_subsection("Weighted Residual Equation");
    std::cout << "∫ w(x) R(x) dx = 0 for each weight function w\n\n";

    // Use trial functions as weights (Galerkin method)
    for (size_t i = 0; i < trial_functions.size(); ++i) {
        double weighted_residual = WeightedResidualMethods::weightedResidualEquation(
            L_operator, approx_solution, f, trial_functions[i], 0.0, 1.0);

        std::cout << "Weight w₁₀(x) = φ₁₀(x): ∫ w R dx = "
                  << weighted_residual << "\n";
    }
}

/**
 * Demonstrate Collocation Method
 */
void demo_collocation_method() {
    print_section("COLLOCATION METHOD");

    print_subsection("Method Description");
    std::cout << "Collocation: Enforce R(xᵢ) = 0 at collocation points\n";
    std::cout << "Weight functions: w(x) = δ(x - xᵢ)\n\n";

    // Problem: -u'' = 1, u(0) = u(1) = 0
    // Use n=3 collocation points

    print_subsection("Collocation Points");

    // Chebyshev nodes mapped to (0,1)
    auto chebyshev_nodes = OtherWeightedResidualMethods::CollocationMethod::chebyshevNodes(3);

    std::cout << "Using Chebyshev collocation points:\n";
    for (size_t i = 0; i < chebyshev_nodes.size(); ++i) {
        std::cout << "x₁₀ = " << chebyshev_nodes[i] << "\n";
    }

    print_subsection("Linear System");

    // Trial functions: φᵢ(x) = sin(iπx)
    std::vector<std::function<double(double)>> trial_funcs;
    int n_basis = 3;

    for (int i = 1; i <= n_basis; ++i) {
        trial_funcs.push_back([i](double x) {
            return std::sin(i * M_PI * x);
        });
    }

    // Differential operator: L[φᵢ] = -φᵢ''
    auto L_op = [](std::function<double(double)> phi, double x) {
        // For sin(nπx): -d²/dx² sin(nπx) = n²π² sin(nπx)
        // Numerical computation for generality
        double h = 1e-5;
        return -(phi(x + h) - 2.0 * phi(x) + phi(x - h)) / (h * h);
    };

    auto source = [](double) { return 1.0; };

    auto solution = OtherWeightedResidualMethods::CollocationMethod::solveCollocation(
        L_op, trial_funcs, source, chebyshev_nodes);

    std::cout << "\nCollocation solution coefficients:\n";
    for (size_t i = 0; i < solution.size(); ++i) {
        std::cout << "c₁₀ = " << solution[i] << "\n";
    }

    // Evaluate approximation
    print_subsection("Solution Quality");

    auto exact = [](double x) { return 0.5 * x * (1.0 - x); };
    auto approx = [&](double x) {
        double result = 0.0;
        for (size_t i = 0; i < solution.size(); ++i) {
            result += solution[i] * trial_funcs[i](x);
        }
        return result;
    };

    std::cout << std::setw(10) << "x"
              << std::setw(20) << "Exact"
              << std::setw(20) << "Approx"
              << std::setw(20) << "Error\n";
    std::cout << std::string(70, '-') << "\n";

    for (double x = 0.2; x <= 0.8; x += 0.2) {
        double u_exact = exact(x);
        double u_approx = approx(x);
        std::cout << std::setw(10) << x
                  << std::setw(20) << u_exact
                  << std::setw(20) << u_approx
                  << std::setw(20) << std::abs(u_exact - u_approx) << "\n";
    }
}

/**
 * Demonstrate Subdomain Method
 */
void demo_subdomain_method() {
    print_section("SUBDOMAIN METHOD");

    print_subsection("Method Description");
    std::cout << "Subdomain: Enforce ∫_Ωᵢ R dx = 0 on each subdomain\n";
    std::cout << "Weight functions: w(x) = 1 on Ωᵢ, 0 elsewhere\n\n";

    // Partition domain [0,1] into subdomains
    int n_subdomains = 4;
    auto subdomains = OtherWeightedResidualMethods::SubdomainMethod::partitionDomain(
        0.0, 1.0, n_subdomains);

    std::cout << "Domain [0, 1] partitioned into " << n_subdomains << " subdomains:\n";
    for (size_t i = 0; i < subdomains.size(); ++i) {
        std::cout << "Ω₁₀ = [" << subdomains[i].first
                  << ", " << subdomains[i].second << "]\n";
    }

    print_subsection("Residual Integration");

    // Example: u ≈ x(1-x) for -u'' = 1
    auto L_op = [](std::function<double(double)> u, double x, double h = 1e-5) {
        return -(u(x + h) - 2.0 * u(x) + u(x - h)) / (h * h);
    };

    auto approx_u = [](double x) { return x * (1.0 - x); };
    auto source = [](double) { return 1.0; };

    std::cout << "\nIntegrated residual on each subdomain:\n";
    for (size_t i = 0; i < subdomains.size(); ++i) {
        double residual_integral = OtherWeightedResidualMethods::SubdomainMethod::subdomainResidual(
            L_op, approx_u, source, subdomains[i].first, subdomains[i].second);

        std::cout << "Ω₁₀: ∫ R dx = " << residual_integral << "\n";
    }
}

/**
 * Demonstrate Least Squares Method
 */
void demo_least_squares_method() {
    print_section("LEAST SQUARES METHOD");

    print_subsection("Method Description");
    std::cout << "Least Squares: Minimize J = ∫ R² dx\n";
    std::cout << "Optimality condition: ∂J/∂cᵢ = 0\n\n";

    // Problem: -u'' = 1, u(0) = u(1) = 0
    auto L_op = [](std::function<double(double)> u, double x, double h = 1e-5) {
        return -(u(x + h) - 2.0 * u(x) + u(x - h)) / (h * h);
    };

    auto source = [](double) { return 1.0; };

    // Trial function: u ≈ c x(1-x)
    print_subsection("Objective Functional");

    std::vector<double> c_values = {0.3, 0.4, 0.5, 0.6, 0.7};

    std::cout << "Testing different coefficient values:\n";
    std::cout << std::setw(10) << "c"
              << std::setw(20) << "J = ∫R² dx\n";
    std::cout << std::string(30, '-') << "\n";

    double min_J = 1e10;
    double optimal_c = 0.0;

    for (double c : c_values) {
        auto u_approx = [c](double x) { return c * x * (1.0 - x); };

        double J = OtherWeightedResidualMethods::LeastSquaresMethod::objectiveFunctional(
            L_op, u_approx, source, 0.0, 1.0);

        std::cout << std::setw(10) << c
                  << std::setw(20) << J << "\n";

        if (J < min_J) {
            min_J = J;
            optimal_c = c;
        }
    }

    std::cout << "\nOptimal coefficient: c = " << optimal_c << "\n";
    std::cout << "Minimum J = " << min_J << "\n";

    print_subsection("Solution Comparison");

    auto exact = [](double x) { return 0.5 * x * (1.0 - x); };
    auto optimal_approx = [optimal_c](double x) { return optimal_c * x * (1.0 - x); };

    std::cout << "\n";
    std::cout << std::setw(10) << "x"
              << std::setw(20) << "Exact"
              << std::setw(20) << "Least Squares"
              << std::setw(20) << "Error\n";
    std::cout << std::string(70, '-') << "\n";

    for (double x = 0.2; x <= 0.8; x += 0.2) {
        double u_exact = exact(x);
        double u_approx = optimal_approx(x);
        std::cout << std::setw(10) << x
                  << std::setw(20) << u_exact
                  << std::setw(20) << u_approx
                  << std::setw(20) << std::abs(u_exact - u_approx) << "\n";
    }
}

/**
 * Demonstrate Test Function Selection
 */
void demo_test_function_selection() {
    print_section("TEST FUNCTION SELECTION");

    print_subsection("Hat Functions (Piecewise Linear)");

    auto hat_func = TestFunctionChoice::TestFunctionTypes::hatFunction(0.5, 0.2);

    std::cout << "Hat function centered at x = 0.5 with width 0.2\n";
    std::cout << std::setw(10) << "x" << std::setw(20) << "φ(x)\n";
    std::cout << std::string(30, '-') << "\n";

    for (double x = 0.2; x <= 0.8; x += 0.1) {
        std::cout << std::setw(10) << x << std::setw(20) << hat_func(x) << "\n";
    }

    print_subsection("Polynomial Test Functions");

    auto poly_test = TestFunctionChoice::TestFunctionTypes::polynomialDirichlet(3, 0.0, 1.0);

    std::cout << "\nPolynomial of degree 3 satisfying homogeneous BCs on [0,1]\n";
    std::cout << std::setw(10) << "x" << std::setw(20) << "φ(x)\n";
    std::cout << std::string(30, '-') << "\n";

    for (double x = 0.0; x <= 1.0; x += 0.2) {
        std::cout << std::setw(10) << x << std::setw(20) << poly_test(x) << "\n";
    }

    print_subsection("Trigonometric Test Functions");

    auto trig_test = TestFunctionChoice::TestFunctionTypes::trigonometric(2, 1.0);

    std::cout << "\nTrigonometric function: sin(2πx/L) with L = 1\n";
    std::cout << std::setw(10) << "x" << std::setw(20) << "φ(x)\n";
    std::cout << std::string(30, '-') << "\n";

    for (double x = 0.0; x <= 1.0; x += 0.2) {
        std::cout << std::setw(10) << x << std::setw(20) << trig_test(x) << "\n";
    }

    print_subsection("Linear Independence Check");

    std::vector<std::function<double(double)>> test_funcs = {
        [](double x) { return x; },
        [](double x) { return x * x; },
        [](double x) { return x * x * x; }
    };

    bool independent = TestFunctionChoice::SelectionCriteria::linearlyIndependent(
        test_funcs, 0.0, 1.0);

    std::cout << "\nTest functions: {x, x², x³}\n";
    std::cout << "Linearly independent: " << (independent ? "YES" : "NO") << "\n";

    // Test with linearly dependent set
    std::vector<std::function<double(double)>> dependent_funcs = {
        [](double x) { return x; },
        [](double x) { return 2.0 * x; },  // Linearly dependent!
        [](double x) { return x * x; }
    };

    independent = TestFunctionChoice::SelectionCriteria::linearlyIndependent(
        dependent_funcs, 0.0, 1.0);

    std::cout << "\nTest functions: {x, 2x, x²}\n";
    std::cout << "Linearly independent: " << (independent ? "YES" : "NO") << "\n";
}

/**
 * Main demonstration
 */
int main() {
    std::cout << std::fixed << std::setprecision(8);

    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║             WEIGHTED RESIDUAL METHODS & VARIATIONAL TECHNIQUES            ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════════════╝\n";

    try {
        demo_greens_identities();
        demo_weighted_residual_methods();
        demo_collocation_method();
        demo_subdomain_method();
        demo_least_squares_method();
        demo_test_function_selection();

        std::cout << "\n" << std::string(80, '=') << "\n";
        std::cout << "All demonstrations completed successfully!\n";
        std::cout << std::string(80, '=') << "\n\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
