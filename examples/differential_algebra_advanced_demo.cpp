/**
 * @file differential_algebra_advanced_demo.cpp
 * @brief Advanced differential algebra: intersections, orthonomic systems, and PDEs
 *
 * Demonstrates:
 * - Manifold intersections (Chapter VII)
 * - Orthonomic systems and Riquier theory (Chapter VIII)
 * - Partial differential algebra (Chapter IX)
 */

#include "maths/algebra/differential_algebra.hpp"
#include <iostream>
#include <iomanip>

using namespace maths::algebra;

void printSeparator(const std::string& title) {
    std::cout << "\n=== " << title << " ===" << std::endl;
}

// Demo 1: Manifold Intersections
void demoManifoldIntersections() {
    printSeparator("Demo 1: Manifold Intersections (Chapter VII)");

    // Create two manifolds
    // M1: y'' - y = 0 (exponential/trig solutions)
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_double_prime = DifferentialPolynomial::derivative(0, 2);
    auto eq1 = y_double_prime + y * DifferentialPolynomial::constant(-1.0);

    // M2: y' - 2y = 0 (exponential with rate 2)
    auto y_prime = DifferentialPolynomial::derivative(0, 1);
    auto eq2 = y_prime + y * DifferentialPolynomial::constant(-2.0);

    std::vector<DifferentialPolynomial> system1 = {eq1};
    std::vector<DifferentialPolynomial> system2 = {eq2};

    auto M1 = AlgebraicDifferentialManifold::fromSystem(system1);
    auto M2 = AlgebraicDifferentialManifold::fromSystem(system2);

    std::cout << "Manifold M₁: y'' - y = 0" << std::endl;
    std::cout << "  Solutions: y = C₁e^x + C₂e^(-x)" << std::endl;
    std::cout << "  Dimension: " << M1.dimension(3) << std::endl;

    std::cout << "\nManifold M₂: y' - 2y = 0" << std::endl;
    std::cout << "  Solutions: y = Ce^(2x)" << std::endl;
    std::cout << "  Dimension: " << M2.dimension(3) << std::endl;

    // Compute intersection
    auto intersection = ManifoldIntersection::intersect(M1, M2);
    std::cout << "\nIntersection M₁ ∩ M₂:" << std::endl;
    std::cout << "  Defining equations: " << intersection.defining_system.size() << std::endl;

    int dim_intersection = ManifoldIntersection::intersectionDimension(M1, M2, 3);
    std::cout << "  Dimension: " << dim_intersection << std::endl;

    int order_intersection = ManifoldIntersection::intersectionOrder(M1, M2);
    std::cout << "  Order: " << order_intersection << std::endl;

    bool do_intersect = ManifoldIntersection::doIntersect(M1, M2);
    std::cout << "  Non-trivial intersection? " << (do_intersect ? "Yes" : "No") << std::endl;

    // Kronecker-type theorem
    std::cout << "\nKronecker's theorem analogue:" << std::endl;
    std::cout << "  dim(M₁) + dim(M₂) = " << (M1.dimension(3) + M2.dimension(3)) << std::endl;
    std::cout << "  dim(M₁ ∩ M₂) = " << dim_intersection << std::endl;
}

// Demo 2: Orthonomic Systems
void demoOrthonomicSystems() {
    printSeparator("Demo 2: Orthonomic Systems (Chapter VIII - Riquier)");

    // Create a simple system: y' = y
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_prime = DifferentialPolynomial::derivative(0, 1);
    auto eq = y_prime + y * DifferentialPolynomial::constant(-1.0);

    std::vector<DifferentialPolynomial> system = {eq};

    std::cout << "System: y' - y = 0" << std::endl;

    // Construct orthonomic system
    Ranking ranking(Ranking::ORDERLY);
    auto ortho = OrthonomicSystem::construct(system, ranking);

    std::cout << "\nOrthonomic form:" << std::endl;
    std::cout << "  Number of equations: " << ortho.equations.size() << std::endl;
    std::cout << "  Principal derivatives: " << ortho.principal_derivatives.size() << std::endl;
    std::cout << "  Parametric derivatives: " << ortho.parametric_derivatives.size() << std::endl;

    bool is_ortho = ortho.isOrthonomic();
    std::cout << "\nIs orthonomic? " << (is_ortho ? "Yes" : "No") << std::endl;

    bool is_passive = ortho.isPassive();
    std::cout << "Is passive? " << (is_passive ? "Yes" : "No") << std::endl;

    bool has_solution = ortho.hasLocalSolution();
    std::cout << "Has local solution (Riquier)? " << (has_solution ? "Yes" : "No") << std::endl;

    // Taylor series dissection
    std::cout << "\nTaylor series dissection:" << std::endl;
    std::map<Derivative, double> initial_cond;
    initial_cond[Derivative(0, 0)] = 1.0;  // y(0) = 1

    auto taylor_coeffs = ortho.dissectTaylorSeries(initial_cond);
    std::cout << "  Number of Taylor coefficients: " << taylor_coeffs.size() << std::endl;
    std::cout << "  (Principal derivatives determined by equations)" << std::endl;
}

// Demo 3: Partial Differential Polynomials
void demoPartialDifferentialPolynomials() {
    printSeparator("Demo 3: Partial Differential Polynomials (Chapter IX)");

    // Create partial derivatives: u, u_x, u_y, u_xx, etc.
    auto u = PartialDifferentialPolynomial::partialDerivative(0, {0, 0});      // u
    auto u_x = PartialDifferentialPolynomial::partialDerivative(0, {1, 0});    // ∂u/∂x
    auto u_y = PartialDifferentialPolynomial::partialDerivative(0, {0, 1});    // ∂u/∂y
    auto u_xx = PartialDifferentialPolynomial::partialDerivative(0, {2, 0});   // ∂²u/∂x²
    auto u_yy = PartialDifferentialPolynomial::partialDerivative(0, {0, 2});   // ∂²u/∂y²

    std::cout << "Created partial differential polynomials:" << std::endl;
    std::cout << "  u (order 0), u_x (order 1), u_y (order 1)" << std::endl;
    std::cout << "  u_xx (order 2), u_yy (order 2)" << std::endl;

    // Laplace equation: u_xx + u_yy = 0
    auto laplace = u_xx + u_yy;
    std::cout << "\nLaplace equation: ∇²u = u_xx + u_yy = 0" << std::endl;
    std::cout << "  Total order: " << laplace.totalOrder() << std::endl;
    std::cout << "  Number of terms: " << laplace.terms.size() << std::endl;

    // Wave equation: u_tt - c²u_xx = 0 (using t as second variable)
    auto u_tt = PartialDifferentialPolynomial::partialDerivative(0, {0, 2});
    auto c_squared = PartialDifferentialPolynomial::constant(1.0);  // c² = 1
    auto wave = u_tt + u_xx * c_squared * PartialDifferentialPolynomial::constant(-1.0);

    std::cout << "\nWave equation: u_tt - c²u_xx = 0 (c = 1)" << std::endl;
    std::cout << "  Total order: " << wave.totalOrder() << std::endl;

    // Heat equation: u_t - α u_xx = 0
    auto u_t = PartialDifferentialPolynomial::partialDerivative(0, {0, 1});
    auto alpha = PartialDifferentialPolynomial::constant(1.0);  // α = 1
    auto heat = u_t + u_xx * alpha * PartialDifferentialPolynomial::constant(-1.0);

    std::cout << "\nHeat equation: u_t - αu_xx = 0 (α = 1)" << std::endl;
    std::cout << "  Total order: " << heat.totalOrder() << std::endl;
}

// Demo 4: Partial Differential Ideals
void demoPartialDifferentialIdeals() {
    printSeparator("Demo 4: Partial Differential Ideals (Chapter IX)");

    // Create Laplace equation ideal
    auto u_xx = PartialDifferentialPolynomial::partialDerivative(0, {2, 0});
    auto u_yy = PartialDifferentialPolynomial::partialDerivative(0, {0, 2});
    auto laplace = u_xx + u_yy;

    std::vector<PartialDifferentialPolynomial> generators = {laplace};
    auto ideal = PartialDifferentialIdeal::generate(generators, 2);

    std::cout << "Ideal generated by Laplace equation: ∇²u = 0" << std::endl;
    std::cout << "  Number of generators: " << ideal.generators.size() << std::endl;
    std::cout << "  Independent variables: " << ideal.num_independent_vars << std::endl;

    // Test consistency
    bool consistent = ideal.isConsistent();
    std::cout << "\nIs consistent? " << (consistent ? "Yes" : "No") << std::endl;

    // Solution dimension (number of arbitrary functions)
    int sol_dim = ideal.solutionDimension();
    std::cout << "Solution dimension: " << sol_dim << std::endl;
    std::cout << "(Number of arbitrary functions in general solution)" << std::endl;

    // Decompose into components
    auto components = ideal.decompose();
    std::cout << "\nDecomposition:" << std::endl;
    std::cout << "  Number of irreducible components: " << components.size() << std::endl;
}

// Demo 5: Characteristic Sets for PDEs
void demoPartialCharacteristicSets() {
    printSeparator("Demo 5: Characteristic Sets for PDEs (Chapter IX)");

    // Wave equation system
    auto u_tt = PartialDifferentialPolynomial::partialDerivative(0, {0, 2});
    auto u_xx = PartialDifferentialPolynomial::partialDerivative(0, {2, 0});
    auto wave = u_tt + u_xx * PartialDifferentialPolynomial::constant(-1.0);

    std::vector<PartialDifferentialPolynomial> system = {wave};
    auto ideal = PartialDifferentialIdeal::generate(system, 2);

    std::cout << "Wave equation: u_tt - u_xx = 0" << std::endl;

    // Construct characteristic set
    auto char_set = PartialCharacteristicSet::construct(ideal);

    std::cout << "\nCharacteristic set:" << std::endl;
    std::cout << "  Size: " << char_set.polynomials.size() << std::endl;
    std::cout << "  (Includes integrability conditions)" << std::endl;

    // Check for singular components
    bool has_singular = char_set.hasSingularComponents();
    std::cout << "\nHas singular components? " << (has_singular ? "Yes" : "No") << std::endl;

    // Decompose
    auto components = char_set.decompose();
    std::cout << "Number of components: " << components.size() << std::endl;
}

// Demo 6: PDE Theorems
void demoPDETheorems() {
    printSeparator("Demo 6: PDE Theorems (Chapter IX)");

    // Laplace equation
    auto u_xx = PartialDifferentialPolynomial::partialDerivative(0, {2, 0});
    auto u_yy = PartialDifferentialPolynomial::partialDerivative(0, {0, 2});
    auto laplace = u_xx + u_yy;

    std::vector<PartialDifferentialPolynomial> system = {laplace};
    auto ideal = PartialDifferentialIdeal::generate(system, 2);

    std::cout << "Laplace equation: ∇²u = 0" << std::endl;

    // Theorem of zeros
    bool has_zero = PDETheorems::hasGenericZero(ideal);
    std::cout << "\nHas generic zero (Theorem of zeros)? " << (has_zero ? "Yes" : "No") << std::endl;

    // Compatibility
    bool compatible = PDETheorems::isCompatible(system);
    std::cout << "Is compatible (mixed partials commute)? " << (compatible ? "Yes" : "No") << std::endl;

    // Cauchy-Kovalevskaya
    bool has_local = PDETheorems::hasLocalAnalyticSolution(ideal);
    std::cout << "Has local analytic solution (Cauchy-Kovalevskaya)? "
              << (has_local ? "Yes" : "No") << std::endl;

    std::cout << "\nPhysical interpretation:" << std::endl;
    std::cout << "  Laplace equation describes:" << std::endl;
    std::cout << "  - Steady-state heat distribution" << std::endl;
    std::cout << "  - Electrostatic potential" << std::endl;
    std::cout << "  - Incompressible fluid flow" << std::endl;
}

// Demo 7: Practical PDE Examples
void demoPracticalPDEs() {
    printSeparator("Demo 7: Classical PDEs in Physics");

    std::cout << "1. Laplace Equation: ∇²u = 0" << std::endl;
    auto u_xx_lap = PartialDifferentialPolynomial::partialDerivative(0, {2, 0});
    auto u_yy_lap = PartialDifferentialPolynomial::partialDerivative(0, {0, 2});
    auto laplace = u_xx_lap + u_yy_lap;
    std::cout << "   Order: " << laplace.totalOrder() << std::endl;
    std::cout << "   Type: Elliptic" << std::endl;
    std::cout << "   Applications: Electrostatics, steady heat" << std::endl;

    std::cout << "\n2. Wave Equation: u_tt - c²∇²u = 0" << std::endl;
    auto u_tt = PartialDifferentialPolynomial::partialDerivative(0, {0, 2});
    auto u_xx_wave = PartialDifferentialPolynomial::partialDerivative(0, {2, 0});
    auto wave = u_tt + u_xx_wave * PartialDifferentialPolynomial::constant(-1.0);
    std::cout << "   Order: " << wave.totalOrder() << std::endl;
    std::cout << "   Type: Hyperbolic" << std::endl;
    std::cout << "   Applications: Sound, light, vibrations" << std::endl;

    std::cout << "\n3. Heat Equation: u_t - α∇²u = 0" << std::endl;
    auto u_t = PartialDifferentialPolynomial::partialDerivative(0, {0, 1});
    auto u_xx_heat = PartialDifferentialPolynomial::partialDerivative(0, {2, 0});
    auto heat = u_t + u_xx_heat * PartialDifferentialPolynomial::constant(-1.0);
    std::cout << "   Order: " << heat.totalOrder() << std::endl;
    std::cout << "   Type: Parabolic" << std::endl;
    std::cout << "   Applications: Heat conduction, diffusion" << std::endl;

    std::cout << "\n4. Schrödinger Equation (time-independent): -ℏ²/(2m)∇²ψ + Vψ = Eψ" << std::endl;
    std::cout << "   Order: 2" << std::endl;
    std::cout << "   Type: Elliptic (eigenvalue problem)" << std::endl;
    std::cout << "   Applications: Quantum mechanics" << std::endl;

    // Analyze compatibility
    std::cout << "\nCompatibility analysis for wave equation:" << std::endl;
    std::vector<PartialDifferentialPolynomial> wave_system = {wave};
    auto wave_ideal = PartialDifferentialIdeal::generate(wave_system, 2);

    bool wave_has_solution = PDETheorems::hasLocalAnalyticSolution(wave_ideal);
    std::cout << "  Admits local solution? " << (wave_has_solution ? "Yes" : "No") << std::endl;

    int wave_sol_dim = wave_ideal.solutionDimension();
    std::cout << "  Solution space dimension: " << wave_sol_dim << std::endl;
}

int main() {
    std::cout << "============================================" << std::endl;
    std::cout << "Advanced Differential Algebra Demo" << std::endl;
    std::cout << "Chapters VII-IX: Intersections, Orthonomic Systems, PDEs" << std::endl;
    std::cout << "============================================" << std::endl;

    try {
        demoManifoldIntersections();
        demoOrthonomicSystems();
        demoPartialDifferentialPolynomials();
        demoPartialDifferentialIdeals();
        demoPartialCharacteristicSets();
        demoPDETheorems();
        demoPracticalPDEs();

        std::cout << "\n============================================" << std::endl;
        std::cout << "All advanced demos completed successfully!" << std::endl;
        std::cout << "============================================" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
