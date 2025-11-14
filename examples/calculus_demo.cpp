#include "../include/maths/optimization/metric_regularity.hpp"
#include "../include/maths/calculus/differential/frechet_calculus.hpp"
#include "../include/maths/analysis/subdifferential.hpp"
#include "../include/maths/analysis/convexity.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>

using namespace maths::optimization;
using namespace maths::calculus;
using namespace maths::analysis;

void printHeader(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(70, '=') << "\n\n";
}

void printSection(const std::string& title) {
    std::cout << "\n--- " << title << " ---\n\n";
}

// Example functions for demonstrations
double quadratic(double x) {
    return x * x;
}

double cubic(double x) {
    return x * x * x;
}

std::vector<double> gradient_quadratic(const std::vector<double>& x) {
    // ∇f(x) = 2x for f(x) = ||x||²
    std::vector<double> grad(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        grad[i] = 2.0 * x[i];
    }
    return grad;
}

double norm_squared(const std::vector<double>& x) {
    double sum = 0.0;
    for (double xi : x) {
        sum += xi * xi;
    }
    return sum;
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

    // ========================================
    // 1. METRIC REGULARITY AND STABILITY
    // ========================================
    printHeader("1. METRIC REGULARITY AND STABILITY");

    printSection("Metric Regularity");
    std::cout << MetricRegularity::definition() << "\n";

    printSection("Pseudo-Lipschitz Property (Aubin Property)");
    std::cout << MetricRegularity::pseudoLipschitz() << "\n";

    printSection("Relationship Between Properties");
    std::cout << MetricRegularity::relationships() << "\n";

    printSection("Calmness");
    std::cout << MetricRegularity::calmness() << "\n";

    printSection("Well-Posedness (Tykhonov)");
    std::cout << WellPosedness::tykhonov() << "\n";

    printSection("Levitin-Polyak Well-Posedness");
    std::cout << WellPosedness::levitinPolyak() << "\n";

    printSection("Stegall's Variational Principle");
    std::cout << WellPosedness::stegall() << "\n";

    printSection("Robinson-Ursescu Open Mapping Theorem");
    std::cout << OpenMappingProperty::robinsonUrsescu() << "\n";

    // ========================================
    // 2. DIFFERENTIAL CALCULUS
    // ========================================
    printHeader("2. DIFFERENTIAL CALCULUS IN BANACH SPACES");

    printSection("Directional Derivative");
    std::cout << DirectionalDerivative::definition() << "\n";

    printSection("One-Sided Derivatives");
    std::cout << DirectionalDerivative::oneSided() << "\n";

    printSection("Fréchet Derivative");
    std::cout << FrechetDerivative::definition() << "\n";

    // Demonstrate directional derivative for f(x) = x²
    printSection("Example: Directional Derivative of f(x) = x²");
    double x0 = 2.0;
    double direction = 1.0;
    double h = 0.0001;
    double approx_deriv = (quadratic(x0 + h * direction) - quadratic(x0)) / h;
    double exact_deriv = 2.0 * x0 * direction; // f'(x;v) = 2xv
    std::cout << "At x = " << x0 << ", direction v = " << direction << "\n";
    std::cout << "Approximate directional derivative: " << approx_deriv << "\n";
    std::cout << "Exact: 2xv = " << exact_deriv << "\n";
    std::cout << "Error: " << std::abs(approx_deriv - exact_deriv) << "\n";

    // ========================================
    // 3. CALCULUS RULES
    // ========================================
    printHeader("3. CALCULUS RULES FOR FRÉCHET DERIVATIVES");

    printSection("Calculus Rules (Sum, Product, Chain)");
    std::cout << FrechetDerivative::calculusRules() << "\n";

    printSection("Mean Value Theorem");
    std::cout << FrechetDerivative::meanValueTheorem() << "\n";

    // ========================================
    // 4. INVERSE AND IMPLICIT FUNCTIONS
    // ========================================
    printHeader("4. INVERSE AND IMPLICIT FUNCTION THEOREMS");

    printSection("Inverse Function Theorem");
    std::cout << InverseFunctionTheorem::inverseTheorem() << "\n";

    printSection("Implicit Function Theorem");
    std::cout << InverseFunctionTheorem::implicitTheorem() << "\n";

    printSection("Rank Theorem");
    std::cout << InverseFunctionTheorem::rankTheorem() << "\n";

    // ========================================
    // 5. NEWTON'S METHOD
    // ========================================
    printHeader("5. NEWTON'S METHOD");

    printSection("Newton's Method for Root Finding");
    std::cout << NewtonMethod::method() << "\n";

    printSection("Applications to Optimization");
    std::cout << NewtonMethod::optimization() << "\n";

    // Demonstrate Newton's method for f(x) = x² - 2 (find √2)
    printSection("Example: Finding √2 Using Newton's Method");
    std::cout << "Finding root of f(x) = x² - 2\n";
    std::cout << "Newton iteration: x_{n+1} = x_n - f(x_n)/f'(x_n)\n";
    std::cout << "                         = x_n - (x_n² - 2)/(2x_n)\n\n";

    double x = 1.5; // Initial guess
    std::cout << "Iteration | x_n       | f(x_n)    | Error to √2\n";
    std::cout << std::string(55, '-') << "\n";

    for (int i = 0; i < 6; ++i) {
        double fx = x * x - 2.0;
        double error = std::abs(x - std::sqrt(2.0));
        std::cout << std::setw(9) << i << " | "
                  << std::setw(9) << x << " | "
                  << std::setw(9) << fx << " | "
                  << std::setw(12) << error << "\n";

        // Newton update: x_{n+1} = x_n - f(x_n)/f'(x_n)
        double fprime_x = 2.0 * x;
        x = x - fx / fprime_x;
    }
    std::cout << "\nTrue value: √2 = " << std::sqrt(2.0) << "\n";
    std::cout << "Convergence: Quadratic (errors square at each iteration)\n";

    // ========================================
    // 6. LEGENDRE TRANSFORM
    // ========================================
    printHeader("6. LEGENDRE AND LEGENDRE-FENCHEL TRANSFORMS");

    printSection("Classical Legendre Transform");
    std::cout << LegendreTransform::classical() << "\n";

    printSection("Legendre-Fenchel Transform (Conjugate)");
    std::cout << LegendreTransform::fenchel() << "\n";

    // Example: Legendre transform of f(x) = x²/2
    printSection("Example: Legendre Transform of f(x) = x²/2");
    std::cout << "For f(x) = x²/2:\n";
    std::cout << "  f'(x) = x, so inverse gives x = p\n";
    std::cout << "  f*(p) = px - f(x) = p·p - p²/2 = p²/2\n";
    std::cout << "The Legendre transform of f(x) = x²/2 is f*(p) = p²/2\n";
    std::cout << "(The function is self-dual under Legendre transform)\n";

    // ========================================
    // 7. SUBDIFFERENTIAL CALCULUS
    // ========================================
    printHeader("7. SUBDIFFERENTIAL CALCULUS");

    printSection("Subdifferential of Convex Functions");
    std::cout << ConvexFunctions::subdifferential() << "\n";

    printSection("Sum Rule for Subdifferentials");
    std::cout << SubdifferentialCalculus::sumRule() << "\n";

    printSection("Linear Composition Rule");
    std::cout << SubdifferentialCalculus::linearComposition() << "\n";

    printSection("Maximum Rule");
    std::cout << SubdifferentialCalculus::maximumRule() << "\n";

    printSection("Marginal Functions and Danskin's Theorem");
    std::cout << SubdifferentialCalculus::marginalFunctions() << "\n";

    // Example: Subdifferential of |x|
    printSection("Example: Subdifferential of f(x) = |x|");
    std::cout << "For f(x) = |x|:\n";
    std::cout << "  ∂f(x) = {sgn(x)}     if x ≠ 0\n";
    std::cout << "  ∂f(0) = [-1, 1]       if x = 0\n\n";

    std::cout << "At x = 2:  ∂f(2) = {1}\n";
    std::cout << "At x = -1: ∂f(-1) = {-1}\n";
    std::cout << "At x = 0:  ∂f(0) = [-1, 1]  (all subgradients)\n";
    std::cout << "\nNote: The subdifferential is set-valued and multi-valued at non-smooth points.\n";

    // ========================================
    // 8. CONVEX DUALITY
    // ========================================
    printHeader("8. CONVEX DUALITY THEORY");

    printSection("Lagrangian Duality");
    std::cout << ConvexDuality::lagrangianDuality() << "\n";

    printSection("Fenchel Duality");
    std::cout << ConvexDuality::fenchelDuality() << "\n";

    // Example: Simple LP duality
    printSection("Example: Linear Programming Duality");
    std::cout << "Primal:  min  c^T x  s.t.  Ax ≥ b, x ≥ 0\n";
    std::cout << "Dual:    max  b^T y  s.t.  A^T y ≤ c, y ≥ 0\n\n";
    std::cout << "Weak Duality: For feasible x, y:  b^T y ≤ c^T x\n";
    std::cout << "Strong Duality: If both have optimal solutions,\n";
    std::cout << "                then optimal values are equal.\n";

    // ========================================
    // 9. DUAL OPTIMIZATION ALGORITHMS
    // ========================================
    printHeader("9. DUAL OPTIMIZATION ALGORITHMS");

    printSection("Dual Algorithms (ADMM, Proximal Methods)");
    std::cout << ConvexDuality::dualAlgorithms() << "\n";

    // ========================================
    // 10. SUPPORT FUNCTIONS
    // ========================================
    printHeader("10. SUPPORT FUNCTIONS AND CONJUGATE CALCULUS");

    printSection("Conjugate Calculus");
    std::cout << ConvexConjugate::conjugateCalculus() << "\n";

    printSection("Support Functions");
    std::cout << ConvexConjugate::supportFunctions() << "\n";

    // Example: Support function of unit ball
    printSection("Example: Support Function of Unit Ball");
    std::cout << "For the unit ball B = {x : ||x|| ≤ 1}:\n";
    std::cout << "  σ_B(y) = sup_{||x||≤1} ⟨y, x⟩ = ||y||\n\n";

    std::vector<double> y1 = {3.0, 4.0};
    double norm_y1 = std::sqrt(y1[0]*y1[0] + y1[1]*y1[1]);
    std::cout << "For y = [3, 4]:\n";
    std::cout << "  σ_B(y) = ||y|| = √(9 + 16) = " << norm_y1 << "\n";

    // ========================================
    // SUMMARY
    // ========================================
    printHeader("SUMMARY");

    std::cout << "This demonstration covered:\n\n";

    std::cout << "1. Metric Regularity and Stability:\n";
    std::cout << "   - Metric regularity and openness\n";
    std::cout << "   - Pseudo-Lipschitz property (Aubin property)\n";
    std::cout << "   - Calmness and stability criteria\n";
    std::cout << "   - Well-posedness (Tykhonov, Levitin-Polyak)\n";
    std::cout << "   - Stegall's principle\n";
    std::cout << "   - Robinson-Ursescu theorem\n\n";

    std::cout << "2. Differential Calculus:\n";
    std::cout << "   - Directional, Gâteaux, Fréchet derivatives\n";
    std::cout << "   - Calculus rules (sum, product, chain)\n";
    std::cout << "   - Mean Value Theorem in Banach spaces\n";
    std::cout << "   - Inverse and Implicit Function Theorems\n\n";

    std::cout << "3. Newton's Method:\n";
    std::cout << "   - Newton iteration for root finding\n";
    std::cout << "   - Kantorovich convergence theorem\n";
    std::cout << "   - Demonstrated finding √2 with quadratic convergence\n\n";

    std::cout << "4. Legendre Transform:\n";
    std::cout << "   - Classical Legendre transform\n";
    std::cout << "   - Legendre-Fenchel (conjugate) transform\n";
    std::cout << "   - Self-duality of f(x) = x²/2\n\n";

    std::cout << "5. Subdifferential Calculus:\n";
    std::cout << "   - Sum, chain, and maximum rules\n";
    std::cout << "   - Marginal functions and Danskin's theorem\n";
    std::cout << "   - Example: ∂|x| with set-valued behavior\n\n";

    std::cout << "6. Convex Duality:\n";
    std::cout << "   - Lagrangian duality and KKT conditions\n";
    std::cout << "   - Fenchel duality\n";
    std::cout << "   - Strong duality and Slater's condition\n";
    std::cout << "   - Linear programming duality\n\n";

    std::cout << "7. Dual Algorithms:\n";
    std::cout << "   - ADMM algorithm\n";
    std::cout << "   - Proximal methods\n";
    std::cout << "   - Convergence guarantees\n\n";

    std::cout << "8. Support Functions:\n";
    std::cout << "   - Definition and properties\n";
    std::cout << "   - Normal cones\n";
    std::cout << "   - Example: Unit ball support function\n\n";

    std::cout << "All implementations are header-only with zero dependencies!\n";
    std::cout << "These tools form the foundation of modern optimization theory.\n";

    return 0;
}
