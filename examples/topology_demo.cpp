#include "../include/maths/topology/metric_spaces.hpp"
#include "../include/maths/topology/topological_spaces.hpp"
#include "../include/maths/topology/set_valued_mappings.hpp"
#include "../include/maths/analysis/convexity.hpp"
#include "../include/maths/optimization/variational_principles.hpp"
#include "../include/maths/optimization/error_bounds.hpp"

#include <iostream>
#include <iomanip>
#include <vector>

using namespace maths::topology;
using namespace maths::analysis;
using namespace maths::optimization;

void printHeader(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(60, '=') << "\n\n";
}

void printSection(const std::string& title) {
    std::cout << "\n--- " << title << " ---\n\n";
}

int main() {
    // ========================================
    // 1. METRIC SPACES
    // ========================================
    printHeader("1. METRIC SPACES");

    // Euclidean metric
    printSection("Euclidean Metric in R³");
    EuclideanSpace euclidean;
    std::vector<double> p1 = {1.0, 2.0, 3.0};
    std::vector<double> p2 = {4.0, 5.0, 6.0};

    double d_euclid = euclidean.distance(p1, p2);
    std::cout << "Distance between [1,2,3] and [4,5,6]: " << d_euclid << "\n";
    std::cout << "Expected: √27 ≈ 5.196\n";

    // Manhattan metric
    printSection("Manhattan (Taxicab) Metric");
    TaxicabSpace taxicab;
    double d_taxi = taxicab.distance(p1, p2);
    std::cout << "Manhattan distance: " << d_taxi << "\n";
    std::cout << "Expected: |4-1| + |5-2| + |6-3| = 9\n";

    // Maximum metric
    printSection("Maximum (Chebyshev) Metric");
    MaximumSpace maximum;
    double d_max = maximum.distance(p1, p2);
    std::cout << "Maximum distance: " << d_max << "\n";
    std::cout << "Expected: max(3, 3, 3) = 3\n";

    // Test convergence
    printSection("Convergence of Sequences");
    std::vector<std::vector<double>> sequence = {
        {1.0, 1.0, 1.0},
        {0.5, 0.5, 0.5},
        {0.1, 0.1, 0.1},
        {0.01, 0.01, 0.01},
        {0.001, 0.001, 0.001}
    };
    std::vector<double> limit = {0.0, 0.0, 0.0};

    bool converges = euclidean.converges(sequence, limit, 0.01);
    std::cout << "Sequence converges to [0,0,0]? "
              << (converges ? "Yes" : "No") << "\n";

    // Test Cauchy
    bool is_cauchy = euclidean.isCauchy(sequence, 0.01);
    std::cout << "Is Cauchy sequence? " << (is_cauchy ? "Yes" : "No") << "\n";

    // Open balls
    printSection("Open Balls");
    std::vector<double> center = {0.0, 0.0, 0.0};
    std::vector<double> test_point = {0.5, 0.5, 0.5};
    double radius = 1.0;

    bool in_ball = euclidean.inOpenBall(test_point, center, radius);
    std::cout << "Point [0.5,0.5,0.5] in open ball B([0,0,0], 1)? "
              << (in_ball ? "Yes" : "No") << "\n";
    std::cout << "Distance from center: " << euclidean.distance(test_point, center) << "\n";

    // ========================================
    // 2. TOPOLOGICAL SPACES
    // ========================================
    printHeader("2. TOPOLOGICAL SPACES");

    printSection("Axioms and Definitions");
    std::cout << TopologicalSpaceProperties::definition() << "\n";

    printSection("Separation Axioms");
    std::cout << TopologicalSpaceProperties::separationAxioms() << "\n";

    printSection("Baire Category Theorem");
    std::cout << BaireSpaces::categoryTheorem() << "\n";

    // ========================================
    // 3. WEAK TOPOLOGY
    // ========================================
    printHeader("3. WEAK TOPOLOGIES");

    printSection("Initial Topology");
    std::cout << WeakTopology::initialTopology() << "\n";

    printSection("Weak Topology on Banach Spaces");
    std::cout << WeakTopology::weakTopologyBanach() << "\n";

    // ========================================
    // 4. SEMICONTINUITY
    // ========================================
    printHeader("4. SEMICONTINUOUS FUNCTIONS");

    printSection("Lower Semicontinuity");
    std::cout << Semicontinuity::lowerSemicontinuous() << "\n";

    printSection("Existence Results");
    std::cout << Semicontinuity::existenceResults() << "\n";

    // ========================================
    // 5. SET-VALUED MAPPINGS
    // ========================================
    printHeader("5. SET-VALUED MAPPINGS (MULTIMAPS)");

    printSection("Definitions");
    std::cout << SetValuedMappings::definitions() << "\n";

    printSection("Continuity Concepts");
    std::cout << SetValuedMappings::continuityConcepts() << "\n";

    printSection("Kakutani Fixed-Point Theorem");
    std::cout << FixedPointTheorems::kakutani() << "\n";

    // ========================================
    // 6. CONVEX SETS AND FUNCTIONS
    // ========================================
    printHeader("6. CONVEXITY");

    printSection("Convex Sets");
    std::cout << ConvexSets::definition() << "\n";

    printSection("Extreme Points and Krein-Milman");
    std::cout << ConvexSets::extremePoints() << "\n";

    printSection("Separation Theorems");
    std::cout << SeparationTheorems::basicSeparation() << "\n";

    printSection("Hahn-Banach Theorem");
    std::cout << SeparationTheorems::hahnBanach() << "\n";

    printSection("Convex Functions");
    std::cout << ConvexFunctions::definition() << "\n";

    printSection("Subdifferential");
    std::cout << ConvexFunctions::subdifferential() << "\n";

    printSection("Conjugate Functions");
    std::cout << ConvexFunctions::conjugate() << "\n";

    // ========================================
    // 7. VARIATIONAL PRINCIPLES
    // ========================================
    printHeader("7. VARIATIONAL PRINCIPLES");

    printSection("Ekeland Variational Principle");
    std::cout << EkelandPrinciple::statement() << "\n";

    printSection("Applications of EVP");
    std::cout << EkelandPrinciple::applications() << "\n";

    printSection("Caristi Fixed-Point Theorem");
    std::cout << CaristiTheorem::statement() << "\n";

    printSection("Drop Theorem");
    std::cout << GeometricPrinciples::dropTheorem() << "\n";

    printSection("Palais-Smale Condition");
    std::cout << PalaisSmale::condition() << "\n";

    printSection("Mountain Pass Theorem");
    std::cout << PalaisSmale::mountainPass() << "\n";

    // ========================================
    // 8. ERROR BOUNDS AND PENALIZATION
    // ========================================
    printHeader("8. ERROR BOUNDS AND PENALIZATION");

    printSection("Error Bounds");
    std::cout << ErrorBounds::definition() << "\n";

    printSection("Hoffman's Error Bound");
    std::cout << ErrorBounds::hoffman() << "\n";

    printSection("Łojasiewicz Inequality");
    std::cout << ErrorBounds::lojasiewicz() << "\n";

    printSection("Armijo Rule");
    std::cout << DecreasePrinciple::armijo() << "\n";

    printSection("Exact Penalty Methods");
    std::cout << Penalization::exactPenalty() << "\n";

    printSection("Augmented Lagrangian");
    std::cout << Penalization::augmentedLagrangian() << "\n";

    printSection("Robust Optimization");
    std::cout << RobustOptimization::robust() << "\n";

    printSection("Tikhonov Regularization");
    std::cout << RobustOptimization::tikhonov() << "\n";

    // ========================================
    // SUMMARY
    // ========================================
    printHeader("SUMMARY");

    std::cout << "This demonstration covered:\n\n";
    std::cout << "1. Metric Spaces:\n";
    std::cout << "   - Euclidean, Manhattan, Maximum metrics\n";
    std::cout << "   - Convergence and Cauchy sequences\n";
    std::cout << "   - Open/closed balls\n\n";

    std::cout << "2. Topological Spaces:\n";
    std::cout << "   - General topology axioms\n";
    std::cout << "   - Separation axioms (T₀ through T₄)\n";
    std::cout << "   - Baire Category Theorem\n";
    std::cout << "   - Uniform Boundedness Theorem\n\n";

    std::cout << "3. Weak Topologies:\n";
    std::cout << "   - Initial and final topologies\n";
    std::cout << "   - Weak topology on Banach spaces\n";
    std::cout << "   - Weak* topology\n\n";

    std::cout << "4. Semicontinuity:\n";
    std::cout << "   - Lower/upper semicontinuous functions\n";
    std::cout << "   - Existence theorems\n\n";

    std::cout << "5. Set-Valued Mappings:\n";
    std::cout << "   - Upper/lower hemicontinuity\n";
    std::cout << "   - Kakutani Fixed-Point Theorem\n";
    std::cout << "   - Michael Selection Theorem\n\n";

    std::cout << "6. Convexity:\n";
    std::cout << "   - Convex sets and functions\n";
    std::cout << "   - Separation theorems\n";
    std::cout << "   - Hahn-Banach Theorem\n";
    std::cout << "   - Subdifferentials and conjugates\n\n";

    std::cout << "7. Variational Principles:\n";
    std::cout << "   - Ekeland Variational Principle\n";
    std::cout << "   - Caristi Fixed-Point Theorem\n";
    std::cout << "   - Palais-Smale condition\n";
    std::cout << "   - Mountain Pass Theorem\n\n";

    std::cout << "8. Error Bounds:\n";
    std::cout << "   - Error bounds and metric regularity\n";
    std::cout << "   - Penalization methods\n";
    std::cout << "   - Augmented Lagrangian\n";
    std::cout << "   - Robust optimization\n\n";

    std::cout << "All implementations are header-only with zero dependencies!\n";

    return 0;
}
