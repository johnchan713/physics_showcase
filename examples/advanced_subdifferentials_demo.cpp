#include "../include/maths/analysis/advanced_subdifferentials.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace maths::analysis;

void printHeader(const std::string& title) {
    std::cout << "\n" << std::string(75, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(75, '=') << "\n\n";
}

void printSection(const std::string& title) {
    std::cout << "\n--- " << title << " ---\n\n";
}

void printVector(const std::string& name, const std::vector<double>& v) {
    std::cout << name << " = [";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << std::fixed << std::setprecision(6) << v[i];
        if (i < v.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

void printSubdifferential(const std::string& name, const SubdifferentialSet& subdiff) {
    std::cout << name << " (size: " << subdiff.size() << "):\n";
    for (size_t i = 0; i < std::min(size_t(5), subdiff.size()); ++i) {
        std::cout << "  ";
        printVector("", subdiff[i]);
    }
    if (subdiff.size() > 5) {
        std::cout << "  ... (" << (subdiff.size() - 5) << " more)\n";
    }
}

int main() {
    std::srand(static_cast<unsigned>(std::time(nullptr)));
    std::cout << std::fixed << std::setprecision(6);

    // ========================================
    // 1. CLARKE JACOBIAN AND SUBDIFFERENTIAL
    // ========================================
    printHeader("1. CLARKE JACOBIAN AND SUBDIFFERENTIAL (5.1.2-5.1.3)");

    printSection("Example 1: Clarke subdifferential of |x|");
    {
        // f(x) = |x| at x = 0
        auto f_abs = [](const Point& x) -> double {
            return std::abs(x[0]);
        };

        Point x0 = {0.0};
        auto subdiff = ClarkeJacobian::computeClarkeSubdifferential(f_abs, x0);

        std::cout << "f(x) = |x| at x = 0\n";
        std::cout << "Analytical: ∂_C f(0) = [-1, 1]\n\n";
        printSubdifferential("Numerical ∂_C f(0)", subdiff);

        // Find min and max
        if (!subdiff.empty()) {
            double min_val = subdiff[0][0];
            double max_val = subdiff[0][0];
            for (const auto& sg : subdiff) {
                min_val = std::min(min_val, sg[0]);
                max_val = std::max(max_val, sg[0]);
            }
            std::cout << "\nRange: [" << min_val << ", " << max_val << "]\n";
        }
    }

    printSection("Example 2: Clarke directional derivative");
    {
        auto f = [](const Point& x) -> double {
            return std::abs(x[0]);
        };

        Point x = {0.0};
        Direction v = {1.0};

        double f_dir = ClarkeJacobian::clarkeDirectionalDerivative(f, x, v);
        std::cout << "f(x) = |x| at x = 0, direction v = 1\n";
        std::cout << "f°(0; 1) = " << f_dir << "\n";
        std::cout << "Analytical: f°(0; 1) = |1| = 1\n";
    }

    printSection("Example 3: Clarke subdifferential of max{x, -x, 0}");
    {
        auto f_max = [](const Point& x) -> double {
            return std::max({x[0], -x[0], 0.0});
        };

        Point x0 = {0.0};
        auto subdiff = ClarkeJacobian::computeClarkeSubdifferential(f_max, x0);

        std::cout << "f(x) = max{x, -x, 0} at x = 0\n";
        std::cout << "Analytical: ∂_C f(0) = conv{-1, 0, 1} = [-1, 1]\n\n";
        printSubdifferential("Numerical ∂_C f(0)", subdiff);
    }

    // ========================================
    // 2. NORMAL AND TANGENT CONES
    // ========================================
    printHeader("2. NORMAL AND TANGENT CONES (5.2)");

    printSection("Example 1: Normal cone to unit ball");
    {
        auto inBall = [](const Point& x) -> bool {
            double norm_sq = 0.0;
            for (double xi : x) {
                norm_sq += xi * xi;
            }
            return norm_sq <= 1.0 + 1e-10;
        };

        Point x = {1.0, 0.0};  // Point on boundary
        auto normals = NormalTangentCones::computeFrechetNormal(inBall, x);

        std::cout << "C = unit ball in ℝ², x = (1, 0) on boundary\n";
        std::cout << "Analytical: N_C(x) = {λ(1, 0) : λ ≥ 0}\n\n";
        printSubdifferential("Numerical N_C(x)", normals);
    }

    printSection("Example 2: Tangent cone to unit ball");
    {
        auto inBall = [](const Point& x) -> bool {
            double norm_sq = 0.0;
            for (double xi : x) {
                norm_sq += xi * xi;
            }
            return norm_sq <= 1.0 + 1e-10;
        };

        Point x = {1.0, 0.0};
        auto tangents = NormalTangentCones::computeTangentCone(inBall, x);

        std::cout << "C = unit ball in ℝ², x = (1, 0)\n";
        std::cout << "Analytical: T_C(x) = {(t₁, t₂) : t₁ ≤ 0}\n\n";
        printSubdifferential("Numerical T_C(x)", tangents);
    }

    printSection("Example 3: Normal cone to halfspace");
    {
        auto inHalfspace = [](const Point& x) -> bool {
            return x[0] >= -1e-10;  // x ≥ 0
        };

        Point x = {0.0, 0.0};
        auto normals = NormalTangentCones::computeFrechetNormal(inHalfspace, x);

        std::cout << "C = {x ∈ ℝ² : x₁ ≥ 0}, point x = (0, 0)\n";
        std::cout << "Analytical: N_C(x) = {(λ, 0) : λ ≤ 0}\n\n";
        printSubdifferential("Numerical N_C(x)", normals);
    }

    // ========================================
    // 3. LIMITING SUBDIFFERENTIALS
    // ========================================
    printHeader("3. LIMITING SUBDIFFERENTIALS (6.1)");

    printSection("Example 1: Limiting subdifferential of |x|");
    {
        auto f = [](const Point& x) -> double {
            return std::abs(x[0]);
        };

        Point x = {0.0};
        auto limiting = LimitingSubdifferential::computeLimiting(f, x);

        std::cout << "f(x) = |x| at x = 0\n";
        std::cout << "Analytical: ∂_L f(0) = [-1, 1]\n\n";
        printSubdifferential("Numerical ∂_L f(0)", limiting);
    }

    printSection("Example 2: Limiting normal to set");
    {
        auto inSet = [](const Point& x) -> bool {
            return x[0] * x[0] + x[1] * x[1] <= 1.0 + 1e-10;
        };

        Point x = {1.0, 0.0};
        auto limiting_normals = LimitingSubdifferential::computeLimitingNormal(inSet, x);

        std::cout << "C = unit ball, x = (1, 0)\n";
        std::cout << "Analytical: N_L C(x) = {λ(1, 0) : λ ≥ 0}\n\n";
        printSubdifferential("Numerical N_L C(x)", limiting_normals);
    }

    printSection("Example 3: Metric regularity check");
    {
        auto F = [](const Point& x) -> Point {
            return {x[0] * x[0], x[1] * x[1]};
        };

        Point x = {1.0, 1.0};
        Point y = F(x);

        bool is_regular = LimitingSubdifferential::isMetricallyRegular(F, x, y);
        std::cout << "F(x₁, x₂) = (x₁², x₂²) at x = (1, 1)\n";
        std::cout << "Metrically regular: " << (is_regular ? "YES" : "NO") << "\n";
        std::cout << "(At this point, Jacobian is non-singular)\n";
    }

    // ========================================
    // 4. CODERIVATIVES
    // ========================================
    printHeader("4. CODERIVATIVES (6.1.2, 6.3)");

    printSection("Example 1: Coderivative computation");
    {
        auto F = [](const Point& x) -> std::vector<Point> {
            return {{x[0] + x[1], x[0] - x[1]}};
        };

        Point x = {1.0, 2.0};
        Point y = {3.0, -1.0};
        Direction y_star = {1.0, 0.5};

        auto coderivative = Coderivatives::computeCoderivative(F, x, y, y_star);

        std::cout << "F(x₁, x₂) = (x₁ + x₂, x₁ - x₂)\n";
        printVector("x", x);
        printVector("y*", y_star);
        std::cout << "\nCoderivative D*F(x, F(x))(y*):\n";
        printSubdifferential("", coderivative);
    }

    printSection("Example 2: Normal to intersection");
    {
        auto inSet1 = [](const Point& x) -> bool {
            return x[0] * x[0] + x[1] * x[1] <= 1.0 + 1e-10;
        };

        auto inSet2 = [](const Point& x) -> bool {
            return x[0] >= -1e-10;
        };

        Point x = {1.0, 0.0};
        auto N_intersection = Coderivatives::normalToIntersection(inSet1, inSet2, x);

        std::cout << "C₁ = unit ball, C₂ = halfspace x₁ ≥ 0\n";
        std::cout << "Point x = (1, 0) on boundary of both\n\n";
        printSubdifferential("N(C₁ ∩ C₂)(x)", N_intersection);
    }

    // ========================================
    // 5. GRADED SUBDIFFERENTIALS
    // ========================================
    printHeader("5. GRADED SUBDIFFERENTIALS (7.1)");

    printSection("Example 1: ε-subdifferential");
    {
        auto f = [](const Point& x) -> double {
            return x[0] * x[0] + x[1] * x[1];
        };

        Point x = {1.0, 1.0};
        double epsilon = 0.1;

        auto eps_subdiff = GradedSubdifferential::computeEpsilonSubdifferential(
            f, x, epsilon);

        std::cout << "f(x₁, x₂) = x₁² + x₂² at x = (1, 1)\n";
        std::cout << "ε = " << epsilon << "\n";
        std::cout << "Analytical: ∂ f(x) = {∇f(x)} = {(2, 2)}\n\n";
        printSubdifferential("Numerical ∂_ε f(x)", eps_subdiff);
    }

    printSection("Example 2: Graded normal cone");
    {
        auto inSet = [](const Point& x) -> bool {
            return x[0] * x[0] + x[1] * x[1] <= 1.0 + 1e-10;
        };

        Point x = {0.707, 0.707};  // Approximately on boundary
        std::vector<double> grades = {0.5, 1.0, 2.0};

        auto graded_normals = GradedSubdifferential::computeGradedNormal(
            inSet, x, grades);

        std::cout << "C = unit ball, x ≈ (1/√2, 1/√2)\n";
        std::cout << "Computing normal cone at multiple grade levels:\n\n";

        for (size_t i = 0; i < grades.size(); ++i) {
            std::cout << "Grade " << grades[i] << ":\n";
            printSubdifferential("", graded_normals[i]);
            std::cout << "\n";
        }
    }

    // ========================================
    // 6. CALCULUS RULES
    // ========================================
    printHeader("6. SUBDIFFERENTIAL CALCULUS (5.3.3, 6.4)");

    printSection("Example 1: Sum rule");
    {
        auto f = [](const Point& x) -> double {
            return std::abs(x[0]);
        };
        auto g = [](const Point& x) -> double {
            return x[0] * x[0];
        };

        Point x = {0.0};
        auto subdiff_f = ClarkeJacobian::computeClarkeSubdifferential(f, x, 1e-6, 20);
        auto subdiff_g = ClarkeJacobian::computeClarkeSubdifferential(g, x, 1e-6, 20);

        auto subdiff_sum = SubdifferentialCalculus::sumRule(subdiff_f, subdiff_g);

        std::cout << "f(x) = |x|, g(x) = x², at x = 0\n";
        std::cout << "h(x) = f(x) + g(x) = |x| + x²\n\n";
        printSubdifferential("∂_C f(0)", subdiff_f);
        std::cout << "\n";
        printSubdifferential("∂_C g(0)", subdiff_g);
        std::cout << "\n";
        printSubdifferential("∂_C(f+g)(0) via sum rule", subdiff_sum);
    }

    printSection("Example 2: Max rule");
    {
        std::vector<std::function<double(const Point&)>> functions = {
            [](const Point& x) { return x[0]; },
            [](const Point& x) { return -x[0]; },
            [](const Point& x) { return 0.0; }
        };

        Point x = {0.0};
        std::vector<SubdifferentialSet> subdiffs;
        std::vector<double> values;

        for (const auto& func : functions) {
            subdiffs.push_back(ClarkeJacobian::computeClarkeSubdifferential(func, x));
            values.push_back(func(x));
        }

        auto max_subdiff = SubdifferentialCalculus::maxRule(subdiffs, values);

        std::cout << "f(x) = max{x, -x, 0} at x = 0\n";
        std::cout << "All three functions active at x = 0\n\n";
        printSubdifferential("∂_C max{x, -x, 0}(0)", max_subdiff);
        std::cout << "\nConvex hull should give [-1, 1]\n";
    }

    printSection("Example 3: Convex hull");
    {
        SubdifferentialSet points = {
            {-1.0},
            {1.0},
            {0.0}
        };

        auto hull = SubdifferentialCalculus::convexHull(points, 20);

        std::cout << "Computing convex hull of {-1, 0, 1}\n";
        std::cout << "Expected: interval [-1, 1]\n\n";
        printSubdifferential("Convex hull (sampled)", hull);

        if (!hull.empty()) {
            double min_val = hull[0][0];
            double max_val = hull[0][0];
            for (const auto& p : hull) {
                min_val = std::min(min_val, p[0]);
                max_val = std::max(max_val, p[0]);
            }
            std::cout << "\nRange: [" << min_val << ", " << max_val << "]\n";
        }
    }

    // ========================================
    // SUMMARY
    // ========================================
    printHeader("SUMMARY");

    std::cout << "Advanced Subdifferential Calculus - Computational Tools:\n\n";

    std::cout << "1. CLARKE CALCULUS (5.1.2-5.1.3):\n";
    std::cout << "   • Clarke Jacobian for vector functions\n";
    std::cout << "   • Clarke directional derivative f°(x; v)\n";
    std::cout << "   • Clarke subdifferential ∂_C f(x)\n\n";

    std::cout << "2. NORMAL AND TANGENT CONES (5.2):\n";
    std::cout << "   • Fréchet normal cone N̂_C(x)\n";
    std::cout << "   • Tangent cone T_C(x)\n";
    std::cout << "   • Projection onto sets\n\n";

    std::cout << "3. LIMITING SUBDIFFERENTIALS (6.1):\n";
    std::cout << "   • Limiting subdifferential ∂_L f(x)\n";
    std::cout << "   • Limiting normal cone N_L C(x)\n";
    std::cout << "   • Metric regularity criterion\n\n";

    std::cout << "4. CODERIVATIVES (6.1.2, 6.3):\n";
    std::cout << "   • Coderivative D*F(x, y)\n";
    std::cout << "   • Normal to intersection\n";
    std::cout << "   • Calculus for set-valued maps\n\n";

    std::cout << "5. GRADED SUBDIFFERENTIALS (7.1-7.2):\n";
    std::cout << "   • ε-subdifferential ∂_ε f(x)\n";
    std::cout << "   • Graded normal cones\n";
    std::cout << "   • Ioffe subdifferential\n\n";

    std::cout << "6. CALCULUS RULES (5.3.3, 6.4):\n";
    std::cout << "   • Sum rule: ∂(f+g) ⊆ ∂f + ∂g\n";
    std::cout << "   • Chain rule: ∂(g∘f) ⊆ ∂g∘Df\n";
    std::cout << "   • Max rule: ∂max ⊆ conv(∂f_i : i active)\n";
    std::cout << "   • Convex hull computation\n\n";

    std::cout << "All implementations use numerical approximations suitable\n";
    std::cout << "for practical computation in finite dimensions.\n";
    std::cout << "Header-only, zero dependencies!\n";

    return 0;
}
