/**
 * @file variational_calculus_demo.cpp
 * @brief Comprehensive demonstration of variational calculus and Lagrangian mechanics
 *
 * Demonstrates:
 * - Contact geometry and contact structures
 * - Lagrangians and Euler-Lagrange equations
 * - Poincaré-Cartan forms and integral invariants
 * - Legendre transformations (Lagrangian ↔ Hamiltonian)
 * - Noether's theorem (symmetries → conservation laws)
 * - Second variation and Jacobi fields
 * - Bäcklund transformations
 * - Conservation laws and field theories
 */

#include "maths/geometry/variational_calculus.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace maths::geometry;

void printSeparator(const std::string& title) {
    std::cout << "\n=== " << title << " ===" << std::endl;
}

// Demo 1: Contact Geometry
void demoContactGeometry() {
    printSeparator("Demo 1: Contact Geometry and Structures");

    // Standard contact structure on J^1(R,R)
    auto contact = ContactStructure::standardContact1Jet();

    std::cout << "Contact manifold J^1(R,R) with coordinates (x, u, p)" << std::endl;
    std::cout << "  Dimension: " << contact.dimension << " (must be odd)" << std::endl;
    std::cout << "  Contact form: θ = du - p dx" << std::endl;

    // Check contact condition
    bool is_contact = contact.isContact();
    std::cout << "\nContact condition dθ^n ∧ θ ≠ 0: "
              << (is_contact ? "Satisfied ✓" : "Failed ✗") << std::endl;

    // Reeb vector field
    std::vector<double> point = {1.0, 0.5, 1.0};  // (x, u, p)
    auto R = contact.reebVectorField(point);

    std::cout << "\nReeb vector field R (satisfies θ(R) = 1, R ⌟ dθ = 0):" << std::endl;
    std::cout << "  R = (";
    for (size_t i = 0; i < R.size(); i++) {
        std::cout << R[i];
        if (i < R.size() - 1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;

    std::cout << "\nPhysical interpretation:" << std::endl;
    std::cout << "  Contact geometry describes non-holonomic constraints" << std::endl;
    std::cout << "  Examples: rolling without slipping, thermodynamics" << std::endl;
}

// Demo 2: Lagrangians and Euler-Lagrange Equations
void demoLagrangians() {
    printSeparator("Demo 2: Lagrangians and Euler-Lagrange Equations");

    // Classic harmonic oscillator: L = ½(u_x)² - ½u²
    auto L_harmonic = [](const std::vector<double>& vars) {
        double u = vars[1];
        double u_x = vars[2];
        return 0.5 * u_x * u_x - 0.5 * u * u;
    };

    Lagrangian harmonic(L_harmonic, 1);

    std::cout << "Harmonic oscillator Lagrangian: L = ½(u')² - ½u²" << std::endl;

    // Test on solution u(x) = sin(x)
    auto u = [](double x) { return std::sin(x); };
    auto u_x = [](double x) { return std::cos(x); };

    double x0 = M_PI / 4.0;
    bool satisfies = harmonic.satisfiesEulerLagrange(u, u_x, x0);

    std::cout << "\nTest solution u(x) = sin(x):" << std::endl;
    std::cout << "  Euler-Lagrange equation: u'' + u = 0" << std::endl;
    std::cout << "  Satisfied? " << (satisfies ? "Yes ✓" : "No ✗") << std::endl;

    // Evaluate Lagrangian
    std::vector<double> vars = {x0, u(x0), u_x(x0)};
    double L_val = harmonic.evaluate(vars);
    double E_L = harmonic.eulerLagrangeOperator(vars);

    std::cout << "  At x = π/4:" << std::endl;
    std::cout << "    L(x, u, u') = " << std::setprecision(4) << L_val << std::endl;
    std::cout << "    E_L(u) = " << E_L << " (should be ≈ 0)" << std::endl;

    // Geodesic Lagrangian: L = ½g_ij u^i_x u^j_x
    std::cout << "\n\nGeodesic Lagrangian: L = ½√(1 + (u')²)" << std::endl;
    std::cout << "  Extremals: straight lines (geodesics)" << std::endl;
    std::cout << "  Euler-Lagrange: u'' / (1 + u'²)^(3/2) = 0" << std::endl;
}

// Demo 3: Poincaré-Cartan Forms
void demoPoincaréCartanForms() {
    printSeparator("Demo 3: Poincaré-Cartan Forms");

    // Simple Lagrangian: L = ½(u_x)²
    auto L_simple = [](const std::vector<double>& vars) {
        double u_x = vars[2];
        return 0.5 * u_x * u_x;
    };

    Lagrangian simple(L_simple, 1);
    auto pc_form = PoincaréCartanForm::fromLagrangian(simple);

    std::cout << "Lagrangian: L = ½(u')²" << std::endl;
    std::cout << "Poincaré-Cartan form: θ_L = L dx + (∂L/∂u') du" << std::endl;
    std::cout << "                           = ½(u')² dx + u' du" << std::endl;

    // Canonical symplectic form
    auto omega = pc_form.canonicalSymplecticForm();

    std::cout << "\nCanonical symplectic form: Ω = dθ_L" << std::endl;
    std::cout << "  Ω = du ∧ dp (in phase space coordinates)" << std::endl;
    std::cout << "  Non-degenerate 2-form on phase space" << std::endl;

    // Integral invariant
    std::vector<std::vector<double>> closed_curve = {
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 1.0},
        {2.0, 0.0, -1.0},
        {1.0, -1.0, 0.0}
    };

    double invariant = pc_form.integralInvariant(closed_curve);

    std::cout << "\nCartan's integral invariant ∮_γ θ_L: " << invariant << std::endl;
    std::cout << "  Preserved under canonical transformations" << std::endl;
    std::cout << "  Generalizes Poincaré's relative integral invariants" << std::endl;
}

// Demo 4: Legendre Transformations
void demoLegendreTransform() {
    printSeparator("Demo 4: Legendre Transformations");

    // Free particle: L = ½m(u_x)²
    double m = 1.0;
    auto L_free = [m](const std::vector<double>& vars) {
        double u_x = vars[2];
        return 0.5 * m * u_x * u_x;
    };

    Lagrangian L(L_free, 1);

    std::cout << "Free particle Lagrangian: L = ½m(u')²  (m = " << m << ")" << std::endl;

    // Fiber derivative: p = ∂L/∂u_x
    std::vector<double> config_vars = {0.0, 1.0, 2.0};  // (x, u, u_x)
    auto phase_vars = LegendreTransform::fiberDerivative(L, config_vars);

    std::cout << "\nFiber derivative (Legendre map):" << std::endl;
    std::cout << "  (x, u, u') ↦ (x, u, p)" << std::endl;
    std::cout << "  p = ∂L/∂u' = m·u' = " << phase_vars[2] << std::endl;

    // Check regularity
    bool regular = LegendreTransform::isRegular(L, config_vars);
    std::cout << "  Regular (det ∂²L/∂u'² ≠ 0)? " << (regular ? "Yes ✓" : "No ✗") << std::endl;

    // Hamiltonian
    auto H = LegendreTransform::hamiltonianFromLagrangian(L);

    std::cout << "\nHamiltonian: H = p·u' - L" << std::endl;
    std::cout << "  For free particle: H = p²/(2m)" << std::endl;
    std::cout << "  H(x, u, p) = " << H.evaluate(phase_vars) << std::endl;

    // Check involution (L** = L)
    bool involution = LegendreTransform::checkInvolution(L);
    std::cout << "\nInvolution property (L** = L): "
              << (involution ? "Verified ✓" : "Failed ✗") << std::endl;

    std::cout << "\nPhysical interpretation:" << std::endl;
    std::cout << "  Legendre transform relates:" << std::endl;
    std::cout << "  - Configuration space (x, u, u') ↔ Phase space (x, u, p)" << std::endl;
    std::cout << "  - Lagrangian mechanics ↔ Hamiltonian mechanics" << std::endl;
}

// Demo 5: Noether's Theorem
void demoNoetherTheorem() {
    printSeparator("Demo 5: Noether's Theorem (Symmetries → Conservation Laws)");

    // Free particle: L = ½(u_x)²
    auto L_free = [](const std::vector<double>& vars) {
        double u_x = vars[2];
        return 0.5 * u_x * u_x;
    };

    Lagrangian L(L_free, 1);

    std::cout << "Free particle: L = ½(u')²" << std::endl;
    std::cout << "\nSymmetries and conserved quantities:" << std::endl;

    // 1. Time translation symmetry → Energy
    auto time_symm = NoetherTheorem::timeTranslation();
    std::cout << "\n1. Time translation: x → x + ε" << std::endl;
    std::cout << "   Infinitesimal generator: X = ∂/∂x" << std::endl;

    auto energy = NoetherTheorem::conservedCurrent(L, time_symm);
    std::vector<double> vars = {1.0, 0.5, 1.0};  // (x, u, u_x)
    double E = energy(vars);

    std::cout << "   Conserved current (Energy): E = " << E << std::endl;
    std::cout << "   For free particle: E = ½(u')² = const" << std::endl;

    // 2. Space translation symmetry → Momentum
    auto space_symm = NoetherTheorem::spaceTranslation();
    std::cout << "\n2. Space translation: u → u + ε" << std::endl;
    std::cout << "   Infinitesimal generator: X = ∂/∂u" << std::endl;

    auto momentum = NoetherTheorem::conservedCurrent(L, space_symm);
    double P = momentum(vars);

    std::cout << "   Conserved current (Momentum): P = " << P << std::endl;
    std::cout << "   For free particle: P = ∂L/∂u' = u' = const" << std::endl;

    // Verify conservation
    bool conserved = NoetherTheorem::verifyConservation(L, time_symm, vars);
    std::cout << "\nConservation verified: " << (conserved ? "Yes ✓" : "No ✗") << std::endl;

    std::cout << "\nNoether's Theorem:" << std::endl;
    std::cout << "  Continuous symmetry ⟹ Conservation law" << std::endl;
    std::cout << "  - Time translation → Energy conservation" << std::endl;
    std::cout << "  - Space translation → Momentum conservation" << std::endl;
    std::cout << "  - Rotation → Angular momentum conservation" << std::endl;
}

// Demo 6: Conservation Laws
void demoConservationLaws() {
    printSeparator("Demo 6: Conservation Laws for Differential Equations");

    std::cout << "Conservation law: D_t(ρ) + D_x(J) = 0" << std::endl;
    std::cout << "  ρ = conserved density" << std::endl;
    std::cout << "  J = flux" << std::endl;

    // KdV equation: u_t + u·u_x + u_xxx = 0
    std::cout << "\nKorteweg-de Vries (KdV) equation: u_t + u·u_x + u_xxx = 0" << std::endl;
    std::cout << "\nInfinitely many conservation laws:" << std::endl;

    std::cout << "  1. Mass:     ∫ u dx = const" << std::endl;
    std::cout << "     ρ₁ = u,   J₁ = ½u² + u_xx" << std::endl;

    std::cout << "\n  2. Momentum: ∫ u² dx = const" << std::endl;
    std::cout << "     ρ₂ = u²,  J₂ = ⅔u³ + 2u·u_xx - u_x²" << std::endl;

    std::cout << "\n  3. Energy:   ∫ (u³ - u_x²) dx = const" << std::endl;
    std::cout << "     ρ₃ = ½u³ - ½u_x²" << std::endl;

    std::cout << "\nPhysical significance:" << std::endl;
    std::cout << "  - Completely integrable systems have infinitely many" << std::endl;
    std::cout << "  - Conservation laws constrain dynamics" << std::endl;
    std::cout << "  - Enable exact solution via inverse scattering" << std::endl;
}

// Demo 7: Second Variation and Jacobi Fields
void demoSecondVariation() {
    printSeparator("Demo 7: Second Variation and Stability");

    // Harmonic oscillator: L = ½(u_x)² - ½u²
    auto L_harmonic = [](const std::vector<double>& vars) {
        double u = vars[1];
        double u_x = vars[2];
        return 0.5 * u_x * u_x - 0.5 * u * u;
    };

    Lagrangian L(L_harmonic, 1);
    SecondVariation sv(L);

    std::cout << "Harmonic oscillator: L = ½(u')² - ½u²" << std::endl;
    std::cout << "Extremal: u₀(x) = sin(x)" << std::endl;

    auto u0 = [](double x) { return std::sin(x); };

    // Accessory Lagrangian
    auto L2 = sv.accessoryLagrangian(u0);

    std::cout << "\nAccessory Lagrangian (second variation):" << std::endl;
    std::cout << "  L₂(η, η') = ½[(∂²L/∂u²)η² + 2(∂²L/∂u∂u')ηη' + (∂²L/∂u'²)η'²]" << std::endl;

    std::cout << "\nJacobi equation (linearization around extremal):" << std::endl;
    std::cout << "  For harmonic oscillator: η'' + η = 0" << std::endl;
    std::cout << "  Solutions: η = A sin(x) + B cos(x)" << std::endl;

    // Conjugate points
    double x0 = 0.0, x1 = 2.0 * M_PI;
    auto conjugates = sv.conjugatePoints(u0, x0, x1);

    std::cout << "\nConjugate points in [0, 2π]: " << conjugates.size() << " found" << std::endl;
    std::cout << "  (Points where Jacobi field vanishes)" << std::endl;

    // Check if minimum
    bool is_min = sv.isMinimum(u0, 0.0, M_PI / 2.0);
    std::cout << "\nIs u₀ a local minimum on [0, π/2]? "
              << (is_min ? "Yes ✓" : "No ✗") << std::endl;
    std::cout << "  (No conjugate points ⟹ local minimum)" << std::endl;
}

// Demo 8: Bäcklund Transformations
void demoBäcklundTransformations() {
    printSeparator("Demo 8: Bäcklund Transformations");

    std::cout << "Bäcklund transformation: map solution → new solution" << std::endl;

    // Sine-Gordon equation
    std::cout << "\nSine-Gordon equation: u_xt = sin(u)" << std::endl;
    std::cout << "  Describes: relativistic field theory, kinks, solitons" << std::endl;

    double a = 1.0;  // Bäcklund parameter
    auto bt = BäcklundTransformation::sineGordon(a);

    std::cout << "\nBäcklund transformation (parameter a = " << a << "):" << std::endl;
    std::cout << "  u_x + v_x = 2a sin((u+v)/2)" << std::endl;
    std::cout << "  u_t - v_t = (2/a) sin((u-v)/2)" << std::endl;

    std::cout << "\nGiven solution u, these equations determine new solution v" << std::endl;

    std::cout << "\nApplications:" << std::endl;
    std::cout << "  1. Generate new solutions from known ones" << std::endl;
    std::cout << "  2. Multi-soliton solutions via permutability theorem" << std::endl;
    std::cout << "  3. Construct conservation laws" << std::endl;
    std::cout << "  4. Discretization (integrable difference equations)" << std::endl;

    std::cout << "\nClassical equations with Bäcklund transformations:" << std::endl;
    std::cout << "  - Sine-Gordon equation" << std::endl;
    std::cout << "  - Liouville equation" << std::endl;
    std::cout << "  - KdV equation" << std::endl;
    std::cout << "  - Nonlinear Schrödinger equation" << std::endl;
}

// Demo 9: Field Theories (Multiple Fields)
void demoFieldTheories() {
    printSeparator("Demo 9: Euler-Lagrange PDE Systems (Field Theories)");

    std::cout << "Field theory: multiple dependent variables u^α(x^i)" << std::endl;

    // Scalar field: L = ½(∂_μ φ)(∂^μ φ) - ½m²φ² (Klein-Gordon)
    std::cout << "\nKlein-Gordon field:" << std::endl;
    std::cout << "  Lagrangian: L = ½(∂_μ φ)(∂^μ φ) - ½m²φ²" << std::endl;
    std::cout << "  Euler-Lagrange: □φ + m²φ = 0" << std::endl;
    std::cout << "    where □ = ∂²/∂t² - ∇²" << std::endl;

    // Electromagnetic field
    std::cout << "\nElectromagnetic field:" << std::endl;
    std::cout << "  Lagrangian: L = -¼F_μν F^μν" << std::endl;
    std::cout << "    where F_μν = ∂_μ A_ν - ∂_ν A_μ" << std::endl;
    std::cout << "  Euler-Lagrange: ∂_μ F^μν = 0 (Maxwell's equations)" << std::endl;

    // Yang-Mills
    std::cout << "\nYang-Mills gauge theory:" << std::endl;
    std::cout << "  Lagrangian: L = -¼Tr(F_μν F^μν)" << std::endl;
    std::cout << "    where F_μν = ∂_μ A_ν - ∂_ν A_μ + [A_μ, A_ν]" << std::endl;
    std::cout << "  Euler-Lagrange: D_μ F^μν = 0" << std::endl;

    std::cout << "\nDe Donder-Weyl Hamiltonian formulation:" << std::endl;
    std::cout << "  Multi-symplectic geometry" << std::endl;
    std::cout << "  Covariant formulation of field theory" << std::endl;
    std::cout << "  H = ∑_i p^α_i ∂_i u^α - L" << std::endl;
}

// Demo 10: Applications
void demoApplications() {
    printSeparator("Demo 10: Applications in Physics and Engineering");

    std::cout << "Variational calculus in physical sciences:\n" << std::endl;

    std::cout << "1. Classical Mechanics:" << std::endl;
    std::cout << "   - Hamilton's principle (action minimization)" << std::endl;
    std::cout << "   - Lagrangian and Hamiltonian mechanics" << std::endl;
    std::cout << "   - Noether's theorem (symmetries)" << std::endl;

    std::cout << "\n2. Field Theory:" << std::endl;
    std::cout << "   - Electromagnetism (Maxwell)" << std::endl;
    std::cout << "   - General relativity (Einstein-Hilbert action)" << std::endl;
    std::cout << "   - Quantum field theory (path integrals)" << std::endl;

    std::cout << "\n3. Optimal Control:" << std::endl;
    std::cout << "   - Calculus of variations" << std::endl;
    std::cout << "   - Pontryagin maximum principle" << std::endl;
    std::cout << "   - Dynamic programming" << std::endl;

    std::cout << "\n4. Geometric Mechanics:" << std::endl;
    std::cout << "   - Symplectic geometry" << std::endl;
    std::cout << "   - Contact geometry (thermodynamics)" << std::endl;
    std::cout << "   - Momentum maps and reduction" << std::endl;

    std::cout << "\n5. Integrable Systems:" << std::endl;
    std::cout << "   - Soliton equations (KdV, NLS, Sine-Gordon)" << std::endl;
    std::cout << "   - Bäcklund transformations" << std::endl;
    std::cout << "   - Inverse scattering method" << std::endl;

    std::cout << "\n6. Differential Geometry:" << std::endl;
    std::cout << "   - Geodesics and minimal surfaces" << std::endl;
    std::cout << "   - Harmonic maps" << std::endl;
    std::cout << "   - Gauge theory" << std::endl;
}

int main() {
    std::cout << "============================================" << std::endl;
    std::cout << "Variational Calculus and Lagrangian Mechanics" << std::endl;
    std::cout << "Differential Geometric Approach" << std::endl;
    std::cout << "============================================" << std::endl;

    try {
        demoContactGeometry();
        demoLagrangians();
        demoPoincaréCartanForms();
        demoLegendreTransform();
        demoNoetherTheorem();
        demoConservationLaws();
        demoSecondVariation();
        demoBäcklundTransformations();
        demoFieldTheories();
        demoApplications();

        std::cout << "\n============================================" << std::endl;
        std::cout << "All variational calculus demos completed!" << std::endl;
        std::cout << "============================================" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
