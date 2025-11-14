#ifndef MATHS_OPTIMIZATION_METRIC_REGULARITY_HPP
#define MATHS_OPTIMIZATION_METRIC_REGULARITY_HPP

#include <functional>
#include <string>
#include <cmath>

/**
 * @file metric_regularity.hpp
 * @brief Advanced metric regularity, Lipschitz behavior, and openness
 *
 * Implements:
 * - Metric regularity and pseudo-Lipschitz property
 * - Calmness and stability
 * - Aubin property (Lipschitz continuity of multimaps)
 * - Well-posedness and Stegall's principle
 */

namespace maths::optimization {

/**
 * @class MetricRegularity
 * @brief Advanced metric regularity concepts
 */
class MetricRegularity {
public:
    /**
     * @brief Metric regularity definition
     */
    static std::string definition() {
        return "Metric Regularity:\n"
               "\n"
               "A multimap F: X ⇉ Y is metrically regular at (x̄, ȳ) ∈ graph(F)\n"
               "with modulus κ if ∃ neighborhoods U(x̄), V(ȳ):\n"
               "\n"
               "  d(x, F⁻¹(y)) ≤ κ · d(y, F(x))\n"
               "\n"
               "for all x ∈ U, y ∈ V.\n"
               "\n"
               "Interpretation:\n"
               "- Inverse F⁻¹ is (locally) Lipschitz continuous\n"
               "- Stability of solutions under perturbations\n"
               "- Newton's method converges linearly\n"
               "\n"
               "Special cases:\n"
               "- Single-valued: F(x) = {f(x)}\n"
               "  ⇒ Metric regularity = surjectivity of f'(x̄)\n"
               "\n"
               "- Constraint system: F(x) = {0}, x ∈ C\n"
               "  ⇒ Metric regularity = constraint qualification\n"
               "\n"
               "Applications:\n"
               "- Implicit function theorem generalizations\n"
               "- Sensitivity analysis\n"
               "- Newton method convergence\n"
               "- Constraint qualification in optimization";
    }

    /**
     * @brief Pseudo-Lipschitz property (Aubin property)
     */
    static std::string pseudoLipschitz() {
        return "Pseudo-Lipschitz Property (Aubin Property):\n"
               "\n"
               "F: X ⇉ Y has pseudo-Lipschitz property at (x̄, ȳ)\n"
               "with modulus L if ∃ neighborhoods U(x̄), V(ȳ):\n"
               "\n"
               "  F(x₁) ∩ V ⊆ F(x₂) + L·d(x₁, x₂)·B_Y\n"
               "\n"
               "for all x₁, x₂ ∈ U.\n"
               "\n"
               "Equivalently:\n"
               "  d(y₁, F(x₂)) ≤ L·d(x₁, x₂)\n"
               "\n"
               "for all x₁, x₂ ∈ U, y₁ ∈ F(x₁) ∩ V.\n"
               "\n"
               "Key insight: Graph is locally Lipschitz continuous!\n"
               "\n"
               "Relationship to metric regularity:\n"
               "  F pseudo-Lipschitz ⟺ F⁻¹ metrically regular\n"
               "\n"
               "For convex-valued F:\n"
               "- Pseudo-Lipschitz ⇒ hemicontinuity\n"
               "- Convex graph ⇒ simpler characterizations\n"
               "\n"
               "Characterizations:\n"
               "1. Limiting normal cone criterion\n"
               "2. Coderivative criterion: D*F(x̄,ȳ)(0) = {0}\n"
               "3. For optimization: LICQ, MFCQ conditions";
    }

    /**
     * @brief Calmness
     */
    static std::string calmness() {
        return "Calmness:\n"
               "\n"
               "F: X ⇉ Y is calm at (x̄, ȳ) ∈ graph(F) if:\n"
               "∃ κ, δ > 0: ∀x with d(x, x̄) < δ,\n"
               "\n"
               "  d(ȳ, F(x)) ≤ κ·d(x, x̄)\n"
               "\n"
               "Interpretation:\n"
               "- Weaker than pseudo-Lipschitz (one-sided)\n"
               "- Upper estimate on solution movement\n"
               "- Sufficient for stability analysis\n"
               "\n"
               "For solution map S: P ⇉ X of parametric problem:\n"
               "  min_{x∈C(p)} f(x, p)\n"
               "\n"
               "Calmness ⇒ Lipschitz stability of optimal value\n"
               "\n"
               "Criteria:\n"
               "1. Second-order sufficient conditions\n"
               "2. Linear independence constraint qualification\n"
               "3. Strict complementarity\n"
               "\n"
               "Applications:\n"
               "- Sensitivity of optimal solutions\n"
               "- Error bounds for algorithms\n"
               "- Stability of equilibria";
    }

    /**
     * @brief Relationships between properties
     */
    static std::string relationships() {
        return "Hierarchy of Regularity Properties:\n"
               "\n"
               "For F: X ⇉ Y at (x̄, ȳ):\n"
               "\n"
               "Metric regularity of F⁻¹\n"
               "         ⇕ (equivalent)\n"
               "Pseudo-Lipschitz property of F (Aubin)\n"
               "         ⇓ (implies)\n"
               "Upper hemicontinuity of F\n"
               "\n"
               "Also:\n"
               "Pseudo-Lipschitz ⇒ Calmness\n"
               "(but not conversely)\n"
               "\n"
               "For single-valued f: X → Y:\n"
               "- Metric regularity ⟺ surjectivity of Df(x̄)\n"
               "- Pseudo-Lipschitz ⟺ Lipschitz continuity\n"
               "\n"
               "Mordukhovich criterion:\n"
               "F pseudo-Lipschitz ⟺ D*F(x̄,ȳ)(0) = {0}\n"
               "where D*F is limiting coderivative\n"
               "\n"
               "Robinson's implicit function theorem:\n"
               "If F(x, y) = 0 has metric regularity in y,\n"
               "then solution map x ↦ y(x) is Lipschitz.";
    }
};

/**
 * @class WellPosedness
 * @brief Well-posedness and stability
 */
class WellPosedness {
public:
    /**
     * @brief Tykhonov well-posedness
     */
    static std::string tykhonov() {
        return "Tykhonov Well-Posedness:\n"
               "\n"
               "Optimization problem:\n"
               "  min_{x∈X} f(x)\n"
               "\n"
               "is Tykhonov well-posed if:\n"
               "1. Unique solution x*\n"
               "2. Every minimizing sequence converges to x*\n"
               "\n"
               "Minimizing sequence: f(xₙ) → inf f\n"
               "\n"
               "Consequences:\n"
               "- Stability under perturbations\n"
               "- Any algorithm producing minimizing sequence succeeds\n"
               "- Unique optimal solution\n"
               "\n"
               "Sufficient conditions:\n"
               "- f strictly convex + coercive\n"
               "- f strongly convex (any coercivity)\n"
               "- Compact feasible set + continuous f\n"
               "\n"
               "Generalized well-posedness:\n"
               "Allow multiple solutions, but require:\n"
               "  d(xₙ, argmin f) → 0\n"
               "when f(xₙ) → inf f";
    }

    /**
     * @brief Levitin-Polyak well-posedness
     */
    static std::string levitinPolyak() {
        return "Levitin-Polyak Well-Posedness:\n"
               "\n"
               "Problem is LP well-posed if:\n"
               "∀ε > 0, ∃δ > 0: f(x) ≤ inf f + δ ⇒ d(x, S) ≤ ε\n"
               "\n"
               "where S = argmin f.\n"
               "\n"
               "Interpretation:\n"
               "- Approximate solutions are close to exact solutions\n"
               "- Weaker than Tykhonov (allows multiple minimizers)\n"
               "- Equivalent to having error bound\n"
               "\n"
               "Equivalences:\n"
               "1. LP well-posedness\n"
               "2. Error bound: d(x, S) ≤ c(f(x) - inf f)\n"
               "3. Level sets {f ≤ α} eventually constant\n"
               "\n"
               "For convex problems:\n"
               "LP well-posed ⟺ argmin f compact and non-empty\n"
               "\n"
               "Applications:\n"
               "- Algorithm convergence analysis\n"
               "- Stopping criteria\n"
               "- Sensitivity bounds";
    }

    /**
     * @brief Stegall's variational principle
     */
    static std::string stegall() {
        return "Stegall's Variational Principle:\n"
               "\n"
               "Let X be Banach with Fréchet differentiable norm,\n"
               "f: X → ℝ∪{+∞} lsc, bounded below, x₀ ∈ dom f.\n"
               "\n"
               "Then ∃ Lipschitz convex function g: X → ℝ:\n"
               "  x₀ is unique minimizer of f + g\n"
               "\n"
               "Strengthening of Ekeland's principle:\n"
               "- Ekeland: ∃ point with small subdifferential\n"
               "- Stegall: ∃ perturbation with unique minimizer\n"
               "\n"
               "Applications:\n"
               "1. Existence of minimizers for perturbed problems\n"
               "2. Generic uniqueness of solutions\n"
               "3. Density of well-posed problems\n"
               "\n"
               "Geometric meaning:\n"
               "Most optimization problems (in Baire category sense)\n"
               "can be made well-posed by small perturbation.\n"
               "\n"
               "Requires:\n"
               "- Fréchet differentiable norm (e.g., ℓᵖ, 1<p<∞)\n"
               "- Fails in ℓ¹, ℓ∞, C[0,1]";
    }

    /**
     * @brief Well-posedness in presence of constraints
     */
    static std::string constrainedProblems() {
        return "Well-Posedness for Constrained Problems:\n"
               "\n"
               "Problem: min f(x) s.t. g(x) ≤ 0, h(x) = 0\n"
               "\n"
               "Well-posedness requires:\n"
               "1. Constraint qualification (LICQ, MFCQ, etc.)\n"
               "2. Second-order sufficient conditions\n"
               "3. Strict complementarity\n"
               "\n"
               "Under these conditions:\n"
               "- Unique local minimizer\n"
               "- Unique Lagrange multipliers\n"
               "- Lipschitz stability under perturbations\n"
               "\n"
               "Perturbation of data:\n"
               "  min f(x,p) s.t. g(x,p) ≤ 0, h(x,p) = 0\n"
               "\n"
               "Solution map p ↦ x(p) is:\n"
               "- Single-valued\n"
               "- Lipschitz continuous\n"
               "- Directionally differentiable\n"
               "\n"
               "Value function V(p) = f(x(p), p) is:\n"
               "- Lipschitz continuous\n"
               "- Directionally differentiable\n"
               "- ∇_p V(p) = ∇_p L(x(p), λ(p), p)";
    }
};

/**
 * @class OpenMappingProperty
 * @brief Openness and surjectivity
 */
class OpenMappingProperty {
public:
    /**
     * @brief Open mapping theorem connections
     */
    static std::string connections() {
        return "Openness and Metric Regularity:\n"
               "\n"
               "For F: X → Y, the following are related:\n"
               "\n"
               "1. F is open (maps open sets to open sets)\n"
               "2. F is metrically regular at all points\n"
               "3. F⁻¹ has Lipschitz behavior\n"
               "\n"
               "Banach-Schauder Open Mapping Theorem:\n"
               "If T: X → Y surjective continuous linear\n"
               "between Banach spaces, then T is open.\n"
               "\n"
               "⇒ T⁻¹: Y → X is bounded (continuous)\n"
               "\n"
               "Nonlinear version (Graves):\n"
               "If f: X → Y continuously Fréchet differentiable,\n"
               "f'(x̄) surjective, then f is metrically regular.\n"
               "\n"
               "Lyusternik-Graves theorem:\n"
               "∃ neighborhoods U(x̄), V(ȳ), ∃L > 0:\n"
               "  ∀y ∈ V, ∃x ∈ U: f(x) = y,\n"
               "  d(x, x̄) ≤ L·d(y, ȳ)\n"
               "\n"
               "Consequence: Implicit function theorem\n"
               "as special case of metric regularity.";
    }

    /**
     * @brief Robinson-Ursescu theorem
     */
    static std::string robinsonUrsescu() {
        return "Robinson-Ursescu Open Mapping Theorem:\n"
               "\n"
               "Let F: X ⇉ Y be closed convex process\n"
               "(F(λx + μy) ⊇ λF(x) + μF(y)).\n"
               "\n"
               "If 0 ∈ int(dom F) and range(F) contains 0,\n"
               "then F is metrically regular at (0,0):\n"
               "\n"
               "  0 ∈ int(F(B_X))\n"
               "\n"
               "Applications:\n"
               "1. Constraint systems: Ax ∈ K, x ∈ C\n"
               "   has solution for all b near 0\n"
               "\n"
               "2. Convex optimization: ∂f(x) ∋ 0\n"
               "   solution exists and stable\n"
               "\n"
               "3. Variational inequalities\n"
               "\n"
               "Proof technique:\n"
               "- Uses Baire category argument\n"
               "- Builds nested sequence of sets\n"
               "- Completeness essential\n"
               "\n"
               "Extensions:\n"
               "- Non-convex settings (Clarke)\n"
               "- Infinite-dimensional spaces\n"
               "- Parametric versions";
    }
};

} // namespace maths::optimization

#endif // MATHS_OPTIMIZATION_METRIC_REGULARITY_HPP
