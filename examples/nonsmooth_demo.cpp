#include "../include/maths/analysis/convex_calculus.hpp"
#include "../include/maths/analysis/elementary_subdifferentials.hpp"
#include "../include/maths/analysis/clarke_subdifferentials.hpp"
#include "../include/maths/analysis/convexity.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace maths::analysis;

void printHeader(const std::string& title) {
    std::cout << "\n" << std::string(75, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(75, '=') << "\n\n";
}

void printSection(const std::string& title) {
    std::cout << "\n--- " << title << " ---\n\n";
}

// Example functions for demonstrations
double absolute_value(double x) {
    return std::abs(x);
}

double max_function(double x, double y) {
    return std::max(x, y);
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

    // ========================================
    // 1. FUZZY CALCULUS RULES
    // ========================================
    printHeader("1. FUZZY (APPROXIMATE) CALCULUS FOR SUBDIFFERENTIALS");

    printSection("Fuzzy Sum Rule");
    std::cout << FuzzyCalculus::fuzzySum() << "\n";

    printSection("Fuzzy Chain Rule");
    std::cout << FuzzyCalculus::fuzzyChain() << "\n";

    printSection("Fuzzy Maximum Rule");
    std::cout << FuzzyCalculus::fuzzyMaximum() << "\n";

    printSection("Approximate Minimization Rule");
    std::cout << FuzzyCalculus::approximateMinimization() << "\n";

    // ========================================
    // 2. EXACT CALCULUS RULES
    // ========================================
    printHeader("2. EXACT CALCULUS WITH CONSTRAINT QUALIFICATIONS");

    printSection("Exact Sum Rule");
    std::cout << ExactCalculus::exactSum() << "\n";

    printSection("Exact Chain Rule");
    std::cout << ExactCalculus::exactChain() << "\n";

    printSection("Exact Maximum Rule");
    std::cout << ExactCalculus::exactMaximum() << "\n";

    // ========================================
    // 3. MEAN VALUE THEOREMS
    // ========================================
    printHeader("3. MEAN VALUE THEOREMS FOR SUBDIFFERENTIALS");

    printSection("Convex Mean Value Theorem");
    std::cout << SubdifferentialMVT::convexMVT() << "\n";

    printSection("Non-Smooth Mean Value Theorem (Lebourg)");
    std::cout << SubdifferentialMVT::nonsmoothMVT() << "\n";

    printSection("Approximate Mean Value Inequality");
    std::cout << SubdifferentialMVT::approximateMVT() << "\n";

    // Demonstrate MVT for absolute value
    printSection("Example: Mean Value for f(x) = |x|");
    std::cout << "Consider f(x) = |x| from x = -1 to x = 1\n";
    std::cout << "f(1) - f(-1) = 1 - 1 = 0\n";
    std::cout << "By Lebourg MVT: ∃z ∈ (-1,1), v* ∈ ∂_Cf(z):\n";
    std::cout << "  0 = ⟨v*, 2⟩\n";
    std::cout << "At z = 0: ∂_Cf(0) = [-1, 1]\n";
    std::cout << "Choose v* = 0 ∈ ∂_Cf(0), gives 0 = ⟨0, 2⟩ ✓\n";

    // ========================================
    // 4. SMOOTHNESS OF NORMS
    // ========================================
    printHeader("4. SMOOTHNESS PROPERTIES OF NORMS");

    printSection("Smooth and Strictly Convex Norms");
    std::cout << SmoothnessNorms::smoothNorms() << "\n";

    printSection("Favorable Classes of Banach Spaces");
    std::cout << SmoothnessNorms::favorableSpaces() << "\n";

    // ========================================
    // 5. ELEMENTARY (FRÉCHET) SUBDIFFERENTIAL
    // ========================================
    printHeader("5. ELEMENTARY (FRÉCHET) SUBDIFFERENTIAL");

    printSection("Fréchet Subdifferential Definition");
    std::cout << ElementarySubdifferentials::frechetSubdifferential() << "\n";

    printSection("Proximal Subdifferential");
    std::cout << ElementarySubdifferentials::proximalSubdifferential() << "\n";

    printSection("Limiting (Mordukhovich) Subdifferential");
    std::cout << ElementarySubdifferentials::limitingSubdifferential() << "\n";

    // Example: Subdifferentials of |x|
    printSection("Example: Subdifferentials of f(x) = |x| at x = 0");
    std::cout << "For f(x) = |x| at x = 0:\n\n";
    std::cout << "Proximal:  ∂_Pf(0) = {0}\n";
    std::cout << "  (Quadratic lower bound touches at 0)\n\n";
    std::cout << "Fréchet:   ∂_Ff(0) = ∅\n";
    std::cout << "  (No o(|x|) linear bound)\n\n";
    std::cout << "Limiting:  ∂_Lf(0) = [-1, 1]\n";
    std::cout << "  (Limits of ±1 from both sides)\n\n";
    std::cout << "Clarke:    ∂_Cf(0) = [-1, 1]\n";
    std::cout << "  (Convex hull of limiting gradients)\n\n";
    std::cout << "Hierarchy: ∂_P ⊆ ∂_F ⊆ ∂_L ⊆ ∂_C\n";

    // ========================================
    // 6. CODERIVATIVES
    // ========================================
    printHeader("6. CODERIVATIVES OF SET-VALUED MAPPINGS");

    printSection("Fréchet Coderivative");
    std::cout << Coderivatives::frechetCoderivative() << "\n";

    printSection("Limiting Coderivative");
    std::cout << Coderivatives::limitingCoderivative() << "\n";

    printSection("Normal and Tangent Cones");
    std::cout << Coderivatives::normalTangentCones() << "\n";

    // ========================================
    // 7. VISCOSITY SUBDIFFERENTIALS
    // ========================================
    printHeader("7. VISCOSITY SUBDIFFERENTIALS AND PDE SOLUTIONS");

    printSection("Viscosity Subdifferential Definition");
    std::cout << ViscositySubdifferentials::definition() << "\n";

    printSection("Viscosity Solutions of PDEs");
    std::cout << ViscositySubdifferentials::viscositySolutions() << "\n";

    // ========================================
    // 8. CLARKE SUBDIFFERENTIAL
    // ========================================
    printHeader("8. CLARKE GENERALIZED GRADIENT");

    printSection("Clarke Directional Derivative");
    std::cout << ClarkeSubdifferential::clarkeDirectionalDerivative() << "\n";

    printSection("Clarke Subdifferential Definition");
    std::cout << ClarkeSubdifferential::definition() << "\n";

    printSection("Clarke Calculus Rules");
    std::cout << ClarkeSubdifferential::calculusRules() << "\n";

    printSection("Regularity and Semi-Smoothness");
    std::cout << ClarkeSubdifferential::regularity() << "\n";

    // Example: Clarke subdifferential calculations
    printSection("Example: Clarke Subdifferential of max{x, -x, 0}");
    std::cout << "Let f(x) = max{x, -x, 0} on ℝ\n\n";
    std::cout << "At x = 0 (all three functions active):\n";
    std::cout << "  Active gradients: {1, -1, 0}\n";
    std::cout << "  ∂_Cf(0) = conv{1, -1, 0} = [-1, 1]\n\n";
    std::cout << "At x = 2 (only x active):\n";
    std::cout << "  ∂_Cf(2) = {1}\n\n";
    std::cout << "At x = -2 (only -x active):\n";
    std::cout << "  ∂_Cf(-2) = {-1}\n\n";
    std::cout << "Clarke directional derivative at x = 0:\n";
    std::cout << "  f°(0; v) = max{1·v, -1·v, 0·v} = |v|\n";

    // ========================================
    // 9. APPLICATIONS TO OPTIMIZATION
    // ========================================
    printHeader("9. APPLICATIONS TO NON-SMOOTH OPTIMIZATION");

    printSection("Non-Smooth Optimization");
    std::cout << ClarkeApplications::nonsmoothOptimization() << "\n";

    printSection("Variational Inequalities");
    std::cout << ClarkeApplications::variationalInequalities() << "\n";

    // Example optimization problem
    printSection("Example: Minimizing f(x) = |x| + x²");
    std::cout << "Consider min f(x) = |x| + x² on ℝ\n\n";
    std::cout << "Clarke subdifferential:\n";
    std::cout << "  For x > 0: ∂_Cf(x) = {1 + 2x}\n";
    std::cout << "  For x < 0: ∂_Cf(x) = {-1 + 2x}\n";
    std::cout << "  For x = 0: ∂_Cf(0) = [-1, 1] + {0} = [-1, 1]\n\n";
    std::cout << "Optimality condition: 0 ∈ ∂_Cf(x*)\n\n";
    std::cout << "Checking x* = 0:\n";
    std::cout << "  0 ∈ [-1, 1] ✓\n";
    std::cout << "So x* = 0 is a stationary point.\n\n";
    std::cout << "Second-order check:\n";
    std::cout << "  f''(0+) = 2 > 0, so x* = 0 is a local minimum.\n\n";
    std::cout << "Numerical verification:\n";
    double x_vals[] = {-0.1, 0.0, 0.1};
    for (double x : x_vals) {
        double fx = std::abs(x) + x * x;
        std::cout << "  f(" << x << ") = " << fx << "\n";
    }

    // ========================================
    // 10. PRACTICAL EXAMPLES
    // ========================================
    printHeader("10. PRACTICAL EXAMPLES");

    printSection("ℓ₁ Norm (Sum of Absolute Values)");
    std::cout << "For f(x) = ‖x‖₁ = Σ|x_i| on ℝⁿ:\n\n";
    std::cout << "Clarke subdifferential at x:\n";
    std::cout << "  ∂_Cf(x) = {v : v_i = sgn(x_i) if x_i ≠ 0,\n";
    std::cout << "                  |v_i| ≤ 1 if x_i = 0}\n\n";
    std::cout << "Example in ℝ²:\n";
    std::cout << "  At x = (1, 0): ∂_Cf(x) = {(1, v₂) : |v₂| ≤ 1}\n";
    std::cout << "  At x = (0, 0): ∂_Cf(x) = {v : ‖v‖∞ ≤ 1} (unit ℓ∞ ball)\n\n";
    std::cout << "Application: LASSO regression\n";
    std::cout << "  min (1/2)‖Ax - b‖² + λ‖x‖₁\n";
    std::cout << "Optimality: 0 ∈ A^T(Ax* - b) + λ∂‖x*‖₁\n";

    printSection("ReLU Function: max{0, x}");
    std::cout << "For f(x) = max{0, x}:\n\n";
    std::cout << "Clarke subdifferential:\n";
    std::cout << "  For x > 0: ∂_Cf(x) = {1}\n";
    std::cout << "  For x < 0: ∂_Cf(x) = {0}\n";
    std::cout << "  For x = 0: ∂_Cf(0) = [0, 1]\n\n";
    std::cout << "Neural network layer: σ(Wx + b) where σ(x) = max{0, x}\n";
    std::cout << "Backpropagation uses subdifferential:\n";
    std::cout << "  ∂σ(x) = {1 if x > 0, 0 if x < 0, [0,1] if x = 0}\n";
    std::cout << "In practice: choose ∂σ(0) = 0 or 1 (convention)\n";

    printSection("Distance Function to a Set");
    std::cout << "For d_C(x) = dist(x, C) = inf_{y∈C} ‖x - y‖:\n\n";
    std::cout << "Properties:\n";
    std::cout << "  1. d_C is 1-Lipschitz: |d_C(x) - d_C(y)| ≤ ‖x - y‖\n";
    std::cout << "  2. d_C(x) = 0 ⟺ x ∈ C̄ (closure)\n";
    std::cout << "  3. Convex C ⇒ d_C convex\n\n";
    std::cout << "Clarke subdifferential (x ∉ C):\n";
    std::cout << "  ∂_Cd_C(x) ⊆ {v* : ‖v*‖ = 1, v* ∈ N_C(p)}\n";
    std::cout << "  where p ∈ proj_C(x)\n\n";
    std::cout << "Example C = [0,∞):\n";
    std::cout << "  d_C(x) = max{0, -x}\n";
    std::cout << "  ∂_Cd_C(x) = {-1} for x < 0, {0} for x > 0, [-1,0] for x = 0\n";

    // ========================================
    // SUMMARY
    // ========================================
    printHeader("SUMMARY OF NON-SMOOTH ANALYSIS");

    std::cout << "This demonstration covered:\n\n";

    std::cout << "1. CALCULUS FRAMEWORKS:\n";
    std::cout << "   ├─ Fuzzy calculus (approximate, always valid)\n";
    std::cout << "   ├─ Exact calculus (with constraint qualifications)\n";
    std::cout << "   └─ Mean value theorems (convex, Lebourg, approximate)\n\n";

    std::cout << "2. SUBDIFFERENTIAL HIERARCHIES:\n";
    std::cout << "   Proximal ⊆ Fréchet ⊆ Limiting ⊆ Clarke\n";
    std::cout << "   ∂_P      ⊆ ∂_F     ⊆ ∂_L      ⊆ ∂_C\n\n";
    std::cout << "   Trade-offs:\n";
    std::cout << "   • Proximal: strongest, may be empty\n";
    std::cout << "   • Fréchet: good local theory, may be empty\n";
    std::cout << "   • Limiting: robust, never empty (Asplund spaces)\n";
    std::cout << "   • Clarke: best calculus, may be large\n\n";

    std::cout << "3. CODERIVATIVES:\n";
    std::cout << "   • Generalized derivatives for set-valued maps\n";
    std::cout << "   • Normal and tangent cones\n";
    std::cout << "   • Metric regularity criteria\n\n";

    std::cout << "4. VISCOSITY THEORY:\n";
    std::cout << "   • Subdifferentials via test functions\n";
    std::cout << "   • Viscosity solutions of Hamilton-Jacobi PDEs\n";
    std::cout << "   • Optimal control applications\n\n";

    std::cout << "5. CLARKE SUBDIFFERENTIAL:\n";
    std::cout << "   • Most robust calculus\n";
    std::cout << "   • Always non-empty for Lipschitz f\n";
    std::cout << "   • Exact sum and max rules\n\n";

    std::cout << "6. APPLICATIONS:\n";
    std::cout << "   • Non-smooth optimization (bundle methods)\n";
    std::cout << "   • Variational inequalities\n";
    std::cout << "   • Complementarity problems\n";
    std::cout << "   • Machine learning (ReLU, ℓ₁ regularization)\n\n";

    std::cout << "7. KEY EXAMPLES:\n";
    std::cout << "   • Absolute value: |x|\n";
    std::cout << "   • ℓ₁ norm: ‖x‖₁ (sparsity, LASSO)\n";
    std::cout << "   • Max function: max{f₁,...,f_m}\n";
    std::cout << "   • ReLU: max{0, x} (neural networks)\n";
    std::cout << "   • Distance function: d_C(x)\n\n";

    std::cout << "All implementations are header-only with zero dependencies!\n";
    std::cout << "This framework provides complete tools for modern non-smooth\n";
    std::cout << "variational analysis and optimization.\n";

    return 0;
}
