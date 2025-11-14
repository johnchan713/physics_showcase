#ifndef MATHS_ANALYSIS_CONVEXITY_HPP
#define MATHS_ANALYSIS_CONVEXITY_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>

/**
 * @file convexity.hpp
 * @brief Convex sets, convex functions, and separation theorems
 *
 * Implements:
 * - Convex sets and their properties
 * - Convex functions and conjugates
 * - Subdifferentials
 * - Hahn-Banach and separation theorems
 * - Extreme points and Krein-Milman
 */

namespace maths::analysis {

/**
 * @class ConvexSets
 * @brief Theory of convex sets
 */
class ConvexSets {
public:
    /**
     * @brief Definition and basic properties
     */
    static std::string definition() {
        return "Convex Sets:\n"
               "\n"
               "A set C ⊆ X is convex if:\n"
               "  ∀x, y ∈ C, ∀λ ∈ [0,1]: λx + (1-λ)y ∈ C\n"
               "\n"
               "(line segment between any two points lies in C)\n"
               "\n"
               "Properties:\n"
               "- Intersection of convex sets is convex\n"
               "- Image of convex set under affine map is convex\n"
               "- Preimage of convex set under affine map is convex\n"
               "- Sum C₁ + C₂ = {x + y : x ∈ C₁, y ∈ C₂} is convex\n"
               "- Scalar multiple λC is convex\n"
               "\n"
               "Examples:\n"
               "- Halfspaces: {x : ⟨a, x⟩ ≤ b}\n"
               "- Balls: {x : ‖x - x₀‖ ≤ r}\n"
               "- Ellipsoids, polyhedra, cones";
    }

    /**
     * @brief Interior, closure, and relative interior
     */
    static std::string interiorClosure() {
        return "Interior and Closure of Convex Sets:\n"
               "\n"
               "For convex C:\n"
               "- int(C) is convex (if non-empty)\n"
               "- cl(C) is convex\n"
               "- int(cl(C)) = int(C)\n"
               "\n"
               "Relative interior ri(C):\n"
               "- Interior relative to affine hull aff(C)\n"
               "- ri(C) always non-empty for convex C\n"
               "- ri(C) ⊆ int(C) (equality if C has non-empty interior)\n"
               "\n"
               "Theorem: For convex C, D with int(C) ≠ ∅:\n"
               "  int(C + D) = int(C) + D = C + int(D) = int(C) + int(D)\n"
               "\n"
               "Relative boundary: rbd(C) = cl(C) \\ ri(C)";
    }

    /**
     * @brief Extreme points
     */
    static std::string extremePoints() {
        return "Extreme Points:\n"
               "\n"
               "x ∈ C is an extreme point if:\n"
               "  x = λy + (1-λ)z with y, z ∈ C, λ ∈ (0,1) ⇒ y = z = x\n"
               "\n"
               "(cannot be written as non-trivial convex combination)\n"
               "\n"
               "Krein-Milman Theorem:\n"
               "Every compact convex set in locally convex space\n"
               "is the closed convex hull of its extreme points.\n"
               "\n"
               "K = cl(conv(ext(K)))\n"
               "\n"
               "Examples:\n"
               "- Simplex: vertices are extreme points\n"
               "- Unit ball in ℓᵖ (1 < p < ∞): all boundary points extreme\n"
               "- Unit ball in ℓ¹: only ±eᵢ are extreme\n"
               "- Probability measures: Dirac measures are extreme\n"
               "\n"
               "Choquet's Theorem: Every point in compact convex set\n"
               "is barycenter of probability measure on extreme points.";
    }

    /**
     * @brief Supporting hyperplanes
     */
    static std::string supportingHyperplanes() {
        return "Supporting Hyperplanes:\n"
               "\n"
               "Hyperplane H = {x : ⟨a, x⟩ = b} (a ≠ 0) supports C at x₀ ∈ C if:\n"
               "  ⟨a, x₀⟩ = b and C ⊆ {x : ⟨a, x⟩ ≤ b}\n"
               "\n"
               "(H touches C at x₀, C lies on one side)\n"
               "\n"
               "Theorem: For closed convex C with non-empty interior,\n"
               "every boundary point has supporting hyperplane.\n"
               "\n"
               "Normal cone at x₀ ∈ C:\n"
               "  N_C(x₀) = {v : ⟨v, x - x₀⟩ ≤ 0 ∀x ∈ C}\n"
               "\n"
               "(vectors \"pointing outward\" from C at x₀)\n"
               "\n"
               "Properties:\n"
               "- N_C(x₀) is closed convex cone\n"
               "- N_C(x₀) = {0} for x₀ ∈ int(C)\n"
               "- v ∈ N_C(x₀) ⟺ x₀ ∈ argmin_{x ∈ C} ⟨v, x⟩";
    }
};

/**
 * @class SeparationTheorems
 * @brief Geometric and functional separation
 */
class SeparationTheorems {
public:
    /**
     * @brief Basic separation
     */
    static std::string basicSeparation() {
        return "Separation Theorems:\n"
               "\n"
               "Theorem 1 (Point and closed convex set):\n"
               "If C closed convex, x₀ ∉ C, ∃a, b:\n"
               "  ⟨a, x₀⟩ > b ≥ ⟨a, x⟩ for all x ∈ C\n"
               "\n"
               "Theorem 2 (Two disjoint closed convex sets):\n"
               "If C₁, C₂ closed convex, C₁ ∩ C₂ = ∅,\n"
               "one bounded, ∃a, b:\n"
               "  ⟨a, x⟩ < b < ⟨a, y⟩ for all x ∈ C₁, y ∈ C₂\n"
               "\n"
               "(strict separation requires boundedness)\n"
               "\n"
               "Theorem 3 (Hyperplane separation):\n"
               "If C₁, C₂ convex, int(C₁) ≠ ∅,\n"
               "int(C₁) ∩ C₂ = ∅, ∃a ≠ 0, b:\n"
               "  ⟨a, x⟩ ≤ b ≤ ⟨a, y⟩ for all x ∈ C₁, y ∈ C₂\n"
               "\n"
               "Key: Interior of at least one set needed for separation.";
    }

    /**
     * @brief Hahn-Banach Theorem
     */
    static std::string hahnBanach() {
        return "Hahn-Banach Theorem:\n"
               "\n"
               "Analytic form:\n"
               "Let p: X → ℝ sublinear (p(x+y) ≤ p(x)+p(y), p(λx) = λp(x) for λ≥0).\n"
               "If f: M → ℝ linear on subspace M with f(x) ≤ p(x) on M,\n"
               "then f extends to F: X → ℝ linear with F(x) ≤ p(x) on X.\n"
               "\n"
               "Geometric form (separation):\n"
               "If A, B disjoint convex with int(A) ≠ ∅,\n"
               "∃ hyperplane separating A and B.\n"
               "\n"
               "Consequences:\n"
               "1. X* separates points of X (X normed space)\n"
               "   (∀x ≠ 0, ∃f ∈ X*: f(x) ≠ 0)\n"
               "\n"
               "2. ‖x‖ = sup{|f(x)| : f ∈ X*, ‖f‖ ≤ 1}\n"
               "\n"
               "3. Every bounded linear functional on subspace\n"
               "   extends to whole space with same norm\n"
               "\n"
               "Proof uses Zorn's Lemma (equivalent to AC).";
    }

    /**
     * @brief Farkas Lemma and alternatives
     */
    static std::string farkasLemma() {
        return "Farkas Lemma (Theorem of Alternatives):\n"
               "\n"
               "Exactly one of the following holds:\n"
               "1. ∃x ≥ 0: Ax = b\n"
               "2. ∃y: Aᵀy ≥ 0, bᵀy < 0\n"
               "\n"
               "Geometric interpretation:\n"
               "Either b is in cone generated by columns of A,\n"
               "or there's separating hyperplane.\n"
               "\n"
               "Other forms:\n"
               "- Exactly one holds: Ax ≤ b has solution, or\n"
               "  ∃y ≥ 0: Aᵀy = 0, bᵀy < 0\n"
               "\n"
               "Applications:\n"
               "- Linear programming duality\n"
               "- Optimality conditions\n"
               "- Convex optimization\n"
               "- Game theory\n"
               "\n"
               "Generalizations: Motzkin, Gordan, Stiemke theorems.";
    }
};

/**
 * @class ConvexFunctions
 * @brief Theory of convex functions
 */
class ConvexFunctions {
public:
    /**
     * @brief Definition and properties
     */
    static std::string definition() {
        return "Convex Functions:\n"
               "\n"
               "f: C → ℝ ∪ {+∞} (C convex) is convex if:\n"
               "  f(λx + (1-λ)y) ≤ λf(x) + (1-λ)f(y)\n"
               "for all x, y ∈ C, λ ∈ [0,1].\n"
               "\n"
               "Strictly convex: inequality is strict for x ≠ y, λ ∈ (0,1)\n"
               "Strongly convex: f(·) - (μ/2)‖·‖² is convex for some μ > 0\n"
               "\n"
               "Properties:\n"
               "- Epigraph epi(f) = {(x,t) : f(x) ≤ t} is convex set\n"
               "- Sublevel sets {x : f(x) ≤ α} are convex\n"
               "- f is lower semicontinuous\n"
               "- Pointwise supremum of convex functions is convex\n"
               "- Sum of convex functions is convex\n"
               "- Composition: f convex, g: ℝ → ℝ increasing convex\n"
               "  ⇒ g ∘ f convex\n"
               "\n"
               "Examples: ‖x‖ᵖ (p ≥ 1), -log x, eˣ, max{·}, indicator functions";
    }

    /**
     * @brief Jensen's inequality
     */
    static std::string jensenInequality() {
        return "Jensen's Inequality:\n"
               "\n"
               "For convex f and probability measure μ:\n"
               "  f(∫ x dμ(x)) ≤ ∫ f(x) dμ(x)\n"
               "\n"
               "Discrete version:\n"
               "  f(∑ λᵢxᵢ) ≤ ∑ λᵢf(xᵢ)\n"
               "where λᵢ ≥ 0, ∑λᵢ = 1\n"
               "\n"
               "Special cases:\n"
               "- Arithmetic-geometric mean: (∏xᵢ)^(1/n) ≤ (∑xᵢ)/n\n"
               "- Log-sum inequality\n"
               "- Information theory inequalities\n"
               "\n"
               "Proof: Use supporting hyperplane at E[X].";
    }

    /**
     * @brief Subdifferential
     */
    static std::string subdifferential() {
        return "Subdifferential:\n"
               "\n"
               "For convex f: ℝⁿ → ℝ ∪ {+∞}, the subdifferential at x is:\n"
               "  ∂f(x) = {v ∈ ℝⁿ : f(y) ≥ f(x) + ⟨v, y-x⟩ ∀y}\n"
               "\n"
               "(set of subgradients = supporting hyperplanes to epi(f))\n"
               "\n"
               "Properties:\n"
               "- ∂f(x) is closed convex set (possibly empty)\n"
               "- ∂f(x) non-empty if x ∈ int(dom f)\n"
               "- f differentiable at x ⇒ ∂f(x) = {∇f(x)}\n"
               "- Optimality: x minimizes f ⟺ 0 ∈ ∂f(x)\n"
               "\n"
               "Calculus rules:\n"
               "- ∂(f + g)(x) ⊇ ∂f(x) + ∂g(x)\n"
               "- ∂(λf)(x) = λ∂f(x) (λ > 0)\n"
               "- Chain rule (more complex)\n"
               "\n"
               "Example: f(x) = |x|, ∂f(0) = [-1, 1]";
    }

    /**
     * @brief Conjugate functions
     */
    static std::string conjugate() {
        return "Convex Conjugate (Fenchel-Legendre Transform):\n"
               "\n"
               "For f: X → ℝ ∪ {+∞}, the conjugate is:\n"
               "  f*(y) = sup_{x ∈ X} [⟨y, x⟩ - f(x)]\n"
               "\n"
               "Properties:\n"
               "- f* always convex and lsc (even if f is not)\n"
               "- f** ≤ f (equality if f convex lsc)\n"
               "- Fenchel-Young: f(x) + f*(y) ≥ ⟨x, y⟩\n"
               "  (equality ⟺ y ∈ ∂f(x) ⟺ x ∈ ∂f*(y))\n"
               "\n"
               "Examples:\n"
               "- f(x) = (1/2)‖x‖² ⇒ f*(y) = (1/2)‖y‖²\n"
               "- f(x) = ‖x‖ ⇒ f*(y) = I_{‖y‖≤1}(y)\n"
               "- f(x) = eˣ ⇒ f*(y) = y log y - y (y>0), 0 (y=0)\n"
               "\n"
               "Applications:\n"
               "- Duality theory\n"
               "- Convex optimization\n"
               "- Statistical mechanics";
    }
};

} // namespace maths::analysis

#endif // MATHS_ANALYSIS_CONVEXITY_HPP
