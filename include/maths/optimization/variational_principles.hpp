#ifndef MATHS_OPTIMIZATION_VARIATIONAL_PRINCIPLES_HPP
#define MATHS_OPTIMIZATION_VARIATIONAL_PRINCIPLES_HPP

#include <functional>
#include <cmath>
#include <string>
#include <vector>

/**
 * @file variational_principles.hpp
 * @brief Variational principles and applications
 *
 * Implements:
 * - Ekeland Variational Principle
 * - Caristi Fixed-Point Theorem
 * - Petal Theorem
 * - Drop Theorem
 * - Palais-Smale conditions
 * - Mountain Pass Theorem
 */

namespace maths::optimization {

/**
 * @class EkelandPrinciple
 * @brief The Ekeland Variational Principle and consequences
 */
class EkelandPrinciple {
public:
    /**
     * @brief Statement of Ekeland's Principle
     */
    static std::string statement() {
        return "Ekeland Variational Principle (EVP):\n"
               "\n"
               "Let (X, d) be complete metric space, f: X → ℝ ∪ {+∞}\n"
               "lower semicontinuous, bounded below, f ≢ +∞.\n"
               "\n"
               "For any ε > 0, λ > 0, and u ∈ X with f(u) ≤ inf f + ε,\n"
               "∃v ∈ X such that:\n"
               "  1. f(v) ≤ f(u)\n"
               "  2. d(u, v) ≤ λ\n"
               "  3. f(v) < f(x) + (ε/λ)d(v, x) for all x ≠ v\n"
               "\n"
               "Interpretation:\n"
               "- v is \"almost minimizer\" of f\n"
               "- Condition 3: strict decrease impossible with slope ε/λ\n"
               "- v is minimizer of f(x) + (ε/λ)d(v, x)\n"
               "\n"
               "Key: Completeness essential!\n"
               "(ℚ with f(x) = x² + 1 on [√2, ∞) is counterexample)";
    }

    /**
     * @brief Proof sketch
     */
    static std::string proofSketch() {
        return "Proof of Ekeland's Principle:\n"
               "\n"
               "Define partial order on X:\n"
               "  x ⪯ y if f(y) + (ε/λ)d(x, y) ≤ f(x)\n"
               "\n"
               "Key observations:\n"
               "1. ⪯ is partial order\n"
               "2. Every chain has lower bound in f-value\n"
               "3. By completeness, every chain has infimum\n"
               "\n"
               "Apply Zorn's Lemma: ∃v minimal element.\n"
               "\n"
               "Minimality ⇒ v satisfies condition 3:\n"
               "  If f(x) + (ε/λ)d(v, x) < f(v) for some x,\n"
               "  then v ⪯ x, contradicting minimality of v.\n"
               "\n"
               "Alternative proof: Construct decreasing sequence\n"
               "using completeness to get limit point.";
    }

    /**
     * @brief Consequences
     */
    static std::string consequences() {
        return "Consequences of Ekeland's Principle:\n"
               "\n"
               "1. Approximate minimizers:\n"
               "   inf f attained ⟺ ∀ε > 0, ∃x: f(x) ≤ inf f + ε\n"
               "   and f(y) ≥ f(x) - ε‖y - x‖ for all y\n"
               "\n"
               "2. Existence of minimizing sequences:\n"
               "   If f bounded below, ∃xₙ: f(xₙ) → inf f\n"
               "   and ‖∂f(xₙ)‖ → 0 (approximate critical points)\n"
               "\n"
               "3. Caristi Fixed-Point Theorem (equivalent to EVP!)\n"
               "\n"
               "4. Drop Theorem, Petal Theorem, Flower Petal Theorem\n"
               "\n"
               "5. Bishop-Phelps Theorem: Dense existence of\n"
               "   support functionals\n"
               "\n"
               "6. Borwein-Preiss variational principle (smooth version)\n"
               "\n"
               "All these theorems are equivalent to completeness!";
    }

    /**
     * @brief Applications
     */
    static std::string applications() {
        return "Applications of Ekeland's Principle:\n"
               "\n"
               "1. Optimization:\n"
               "   - Existence of approximate minimizers\n"
               "   - Necessary conditions for optimality\n"
               "   - Penalization methods\n"
               "\n"
               "2. Differential Equations:\n"
               "   - Existence of solutions\n"
               "   - Variational methods\n"
               "   - Critical point theory\n"
               "\n"
               "3. Fixed-Point Theory:\n"
               "   - Caristi's theorem\n"
               "   - Equilibrium problems\n"
               "\n"
               "4. Game Theory:\n"
               "   - Existence of Nash equilibria\n"
               "   - Approximate equilibria\n"
               "\n"
               "5. Control Theory:\n"
               "   - Optimal control\n"
               "   - Calculus of variations\n"
               "\n"
               "Example: Prove continuous function on compact\n"
               "set attains minimum (Weierstrass theorem).";
    }
};

/**
 * @class CaristiTheorem
 * @brief Caristi Fixed-Point Theorem
 */
class CaristiTheorem {
public:
    /**
     * @brief Statement
     */
    static std::string statement() {
        return "Caristi Fixed-Point Theorem:\n"
               "\n"
               "Let (X, d) be complete metric space, T: X → X,\n"
               "φ: X → [0, ∞) lower semicontinuous.\n"
               "\n"
               "If ∀x ∈ X: d(x, Tx) ≤ φ(x) - φ(Tx),\n"
               "then T has a fixed point.\n"
               "\n"
               "Interpretation:\n"
               "- Moving from x to Tx decreases φ by at least d(x, Tx)\n"
               "- φ is Lyapunov function for dynamical system\n"
               "- Generalizes Banach contraction (take φ(x) = cd(x, x₀)/(1-k))\n"
               "\n"
               "Theorem: EVP ⟺ Caristi ⟺ Completeness\n"
               "\n"
               "Proof from EVP: Apply EVP to f(x) = φ(x) + d(x, x₀)\n"
               "with special choice of parameters.";
    }

    /**
     * @brief Comparison with contractions
     */
    static std::string comparisonContractions() {
        return "Caristi vs. Banach Contraction:\n"
               "\n"
               "Banach Contraction Mapping:\n"
               "- d(Tx, Ty) ≤ kd(x, y) for k < 1\n"
               "- Fixed point unique\n"
               "- Iteration converges: Tⁿx → fixed point\n"
               "- Constructive\n"
               "\n"
               "Caristi:\n"
               "- d(x, Tx) ≤ φ(x) - φ(Tx)\n"
               "- Fixed point may not be unique\n"
               "- Iteration need not converge\n"
               "- Existence only (uses Zorn's Lemma)\n"
               "\n"
               "Relation: Banach ⇒ Caristi\n"
               "(take φ(x) = cd(x, x₀)/(1-k), c > 0)\n"
               "\n"
               "But Caristi strictly more general:\n"
               "Example: T(x) = x/2 on [0,1] with\n"
               "φ(x) = x satisfies Caristi but not Banach.";
    }
};

/**
 * @class GeometricPrinciples
 * @brief Drop, Petal, and Flower Petal theorems
 */
class GeometricPrinciples {
public:
    /**
     * @brief Drop Theorem
     */
    static std::string dropTheorem() {
        return "Drop Theorem (Danes):\n"
               "\n"
               "Let C ⊆ Banach space X be closed, x₀ ∉ C, r > 0.\n"
               "Define \"drop\":\n"
               "  D(x₀, C, r) = conv({x₀} ∪ B(C, r))\n"
               "\n"
               "where B(C, r) = {x : d(x, C) ≤ r}.\n"
               "\n"
               "Then ∃z ∈ C: D(x₀, C, r) ∩ C = {z}\n"
               "(drop \"rests\" on single point of C)\n"
               "\n"
               "Geometric meaning:\n"
               "- Imagine water drop hanging from x₀ above surface C\n"
               "- Drop touches C at exactly one point\n"
               "\n"
               "Proof uses Ekeland's Principle.\n"
               "\n"
               "Consequence: Bishop-Phelps theorem on\n"
               "denseness of support points.";
    }

    /**
     * @brief Petal Theorem
     */
    static std::string petalTheorem() {
        return "Petal Theorem:\n"
               "\n"
               "For closed subset C of Banach space, point x₀,\n"
               "and function φ: C → ℝ lsc with inf φ < φ(x₀):\n"
               "\n"
               "∃z ∈ C such that \"petal\" P(z) containing x₀ satisfies:\n"
               "  P(z) ∩ C = {z}\n"
               "\n"
               "where petal is certain geometric region.\n"
               "\n"
               "Interpretation: Strengthening of Ekeland:\n"
               "- Not just one point has descent property\n"
               "- Entire geometric region excluded\n"
               "\n"
               "Flower Petal Theorem: Collection of petals\n"
               "covering region around approximate minimizer.";
    }
};

/**
 * @class PalaisSmale
 * @brief Palais-Smale condition and critical point theory
 */
class PalaisSmale {
public:
    /**
     * @brief Palais-Smale condition
     */
    static std::string condition() {
        return "Palais-Smale Condition:\n"
               "\n"
               "Functional f: X → ℝ (X Banach) satisfies (PS) if:\n"
               "Every sequence (xₙ) with:\n"
               "  1. f(xₙ) bounded\n"
               "  2. ‖df(xₙ)‖ → 0\n"
               "has a convergent subsequence.\n"
               "\n"
               "(PS)_c at level c: same condition with f(xₙ) → c\n"
               "\n"
               "Meaning:\n"
               "- Approximate critical sequences are precompact\n"
               "- Compactness condition for critical points\n"
               "- Rules out \"loss of compactness\" phenomena\n"
               "\n"
               "Examples satisfying (PS):\n"
               "- Coercive functionals (f(x) → ∞ as ‖x‖ → ∞)\n"
               "- Functionals on compact manifolds\n"
               "\n"
               "Examples failing (PS):\n"
               "- f(u) = ∫(|∇u|² - u²) on H¹(ℝⁿ) (translation invariance)\n"
               "- Lack of compactness in critical Sobolev embeddings";
    }

    /**
     * @brief Mountain Pass Theorem
     */
    static std::string mountainPass() {
        return "Mountain Pass Theorem (Ambrosetti-Rabinowitz):\n"
               "\n"
               "Let f ∈ C¹(X, ℝ) satisfy (PS), with:\n"
               "  1. f(0) = 0\n"
               "  2. ∃r, α > 0: f(x) ≥ α for ‖x‖ = r\n"
               "  3. ∃e with ‖e‖ > r: f(e) ≤ 0\n"
               "\n"
               "Then f has critical point with critical value\n"
               "  c = inf_{γ ∈ Γ} max_{t ∈ [0,1]} f(γ(t)) ≥ α\n"
               "\n"
               "where Γ = {γ ∈ C([0,1], X) : γ(0) = 0, γ(1) = e}.\n"
               "\n"
               "Geometric interpretation:\n"
               "- Valley at origin, mountain ring at ‖x‖ = r\n"
               "- Point e beyond mountain (in valley)\n"
               "- Every path from 0 to e must cross mountain\n"
               "- Critical point at \"lowest pass\"\n"
               "\n"
               "Applications: Semilinear elliptic PDEs,\n"
               "existence of solutions via variational methods.";
    }

    /**
     * @brief Deformation lemmas
     */
    static std::string deformationLemmas() {
        return "Deformation Lemmas:\n"
               "\n"
               "If f satisfies (PS)_c and c is not critical value,\n"
               "then ∃ε, δ > 0 and deformation η: [0,1] × X → X:\n"
               "\n"
               "1. η(0, x) = x\n"
               "2. f(η(1, x)) ≤ c - ε for x with f(x) ∈ [c - δ, c + δ]\n"
               "3. η(t, x) = x if |f(x) - c| ≥ δ\n"
               "\n"
               "(can deform sublevel sets past non-critical levels)\n"
               "\n"
               "Consequence: Critical points correspond to\n"
               "topological changes in sublevel sets.\n"
               "\n"
               "Morse Theory: Relate critical points to\n"
               "topology via Morse inequalities:\n"
               "  ∑ dim Hₖ(M) ≤ ∑ #{critical points of index k}";
    }
};

/**
 * @class MetricConvexity
 * @brief Metric convexity and geodesics
 */
class MetricConvexity {
public:
    /**
     * @brief Busemann convexity
     */
    static std::string busemannConvexity() {
        return "Metric Convexity (Busemann):\n"
               "\n"
               "Metric space (X, d) is (uniquely) geodesic if\n"
               "∀x, y ∈ X, ∃(!) geodesic γ: [0,1] → X with:\n"
               "  γ(0) = x, γ(1) = y\n"
               "  d(γ(s), γ(t)) = |t - s| · d(x, y)\n"
               "\n"
               "CAT(0) space: geodesic space where triangles\n"
               "are \"thinner\" than Euclidean triangles\n"
               "(non-positive curvature in Alexandrov sense)\n"
               "\n"
               "Properties:\n"
               "- Midpoints unique\n"
               "- Convex functions well-defined\n"
               "- Many geometric properties of Hilbert spaces\n"
               "\n"
               "Examples:\n"
               "- Euclidean spaces, Hilbert spaces\n"
               "- Trees, ℝ-trees\n"
               "- Simply connected Riemannian manifolds\n"
               "  with non-positive sectional curvature";
    }
};

} // namespace maths::optimization

#endif // MATHS_OPTIMIZATION_VARIATIONAL_PRINCIPLES_HPP
