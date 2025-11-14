#ifndef MATHS_ANALYSIS_SUBDIFFERENTIAL_HPP
#define MATHS_ANALYSIS_SUBDIFFERENTIAL_HPP

#include <functional>
#include <vector>
#include <string>
#include <cmath>

/**
 * @file subdifferential.hpp
 * @brief Advanced convex analysis: subdifferentials and calculus rules
 *
 * Implements:
 * - Subdifferential calculus
 * - Sum rule, chain rule, composition
 * - Marginal functions
 * - Convex duality theory
 * - Applications to optimization
 */

namespace maths::analysis {

/**
 * @class SubdifferentialCalculus
 * @brief Subdifferential and its calculus rules
 */
class SubdifferentialCalculus {
public:
    /**
     * @brief Sum rule for subdifferentials
     */
    static std::string sumRule() {
        return "Sum Rule for Subdifferentials:\n"
               "\n"
               "For convex functions f, g: X → ℝ ∪ {+∞}:\n"
               "\n"
               "  ∂(f + g)(x) ⊇ ∂f(x) + ∂g(x)\n"
               "\n"
               "Equality holds if:\n"
               "1. ri(dom f) ∩ ri(dom g) ≠ ∅ (relative interiors intersect)\n"
               "2. One function is continuous at a point in dom(f) ∩ dom(g)\n"
               "\n"
               "Counterexample (equality fails):\n"
               "  f(x) = I_{(-∞,0]}(x), g(x) = I_{[0,∞)}(x)\n"
               "  ∂f(0) = (-∞, 0], ∂g(0) = [0, ∞)\n"
               "  ∂(f+g)(0) = ℝ ≠ ∂f(0) + ∂g(0) = ℝ\n"
               "  (Actually equal in this case, but pathological)\n"
               "\n"
               "Applications:\n"
               "- Composite optimization: min f(x) + g(x)\n"
               "- Proximal algorithms\n"
               "- Optimality conditions";
    }

    /**
     * @brief Composition with linear operator
     */
    static std::string linearComposition() {
        return "Chain Rule (Linear Composition):\n"
               "\n"
               "For convex f: Y → ℝ, linear A: X → Y:\n"
               "\n"
               "  ∂(f ∘ A)(x) = A* ∂f(Ax)\n"
               "\n"
               "where A*: Y* → X* is adjoint operator.\n"
               "\n"
               "In finite dimensions:\n"
               "  A* = Aᵀ (transpose)\n"
               "  ∂(f(Ax)) = Aᵀ ∂f(Ax)\n"
               "\n"
               "Special case (linear function):\n"
               "  f(x) = ⟨c, x⟩ ⇒ ∂f(x) = {c}\n"
               "\n"
               "Example:\n"
               "  f(x) = ‖Ax - b‖²\n"
               "  ∂f(x) = 2Aᵀ(Ax - b)\n"
               "\n"
               "Applications:\n"
               "- Least squares: min ‖Ax - b‖²\n"
               "- Compressed sensing: min ‖x‖₁ s.t. Ax = b\n"
               "- Regularization: min f(x) + λ‖Ax‖";
    }

    /**
     * @brief Maximum of functions
     */
    static std::string maximumRule() {
        return "Subdifferential of Maximum:\n"
               "\n"
               "For f(x) = max{f₁(x), ..., fₘ(x)} convex:\n"
               "\n"
               "  ∂f(x) = conv{⋃_{i ∈ I(x)} ∂fᵢ(x)}\n"
               "\n"
               "where I(x) = {i : fᵢ(x) = f(x)} (active indices)\n"
               "\n"
               "Interpretation:\n"
               "- Convex hull of subdifferentials of active functions\n"
               "- Only active constraints contribute\n"
               "\n"
               "Example: f(x) = max{0, x} (ReLU)\n"
               "  ∂f(x) = {0} if x < 0\n"
               "  ∂f(x) = {1} if x > 0\n"
               "  ∂f(0) = [0, 1]\n"
               "\n"
               "Example: f(x) = ‖x‖∞ = max_i |xᵢ|\n"
               "  ∂f(x) includes sign vectors for components with |xᵢ| = ‖x‖∞\n"
               "\n"
               "Applications:\n"
               "- Piecewise linear functions\n"
               "- Non-smooth optimization\n"
               "- Support vector machines";
    }

    /**
     * @brief Infimal convolution
     */
    static std::string infimalConvolution() {
        return "Infimal Convolution and Subdifferential:\n"
               "\n"
               "Infimal convolution:\n"
               "  (f □ g)(x) = inf_{y} [f(y) + g(x - y)]\n"
               "\n"
               "Subdifferential:\n"
               "If (f □ g)(x) is attained at y = x - z:\n"
               "  ∂(f □ g)(x) = ∂f(y) ∩ ∂g(z)\n"
               "\n"
               "Conjugate duality:\n"
               "  (f □ g)* = f* + g*\n"
               "\n"
               "Regularization example:\n"
               "  Moreau envelope: e_λ(x) = inf_y [f(y) + (1/2λ)‖x - y‖²]\n"
               "  = (f □ (1/2λ)‖·‖²)(x)\n"
               "\n"
               "Proximal operator:\n"
               "  prox_λf(x) = argmin_y [f(y) + (1/2λ)‖x - y‖²]\n"
               "\n"
               "Properties:\n"
               "- e_λ is C¹ even if f non-smooth\n"
               "- ∇e_λ(x) = (x - prox_λf(x))/λ\n"
               "- e_λ(x) ↑ f(x) as λ ↓ 0\n"
               "\n"
               "Applications:\n"
               "- Smoothing non-smooth functions\n"
               "- Proximal gradient methods\n"
               "- ADMM algorithms";
    }

    /**
     * @brief Marginal functions
     */
    static std::string marginalFunctions() {
        return "Subdifferentials of Marginal Functions:\n"
               "\n"
               "Marginal function:\n"
               "  φ(x) = inf_{y ∈ C} f(x, y)\n"
               "\n"
               "If minimum attained at y(x):\n"
               "\n"
               "  ∂φ(x) = {∂_x f(x, y) : y ∈ argmin_{z∈C} f(x, z)}\n"
               "\n"
               "Danskin's theorem:\n"
               "If f(x, y) differentiable in x, C compact,\n"
               "and argmin y(x) is singleton:\n"
               "\n"
               "  ∇_x φ(x) = ∇_x f(x, y(x))\n"
               "\n"
               "Example (value function in optimization):\n"
               "  V(p) = min_{x} f(x) s.t. g(x) ≤ p\n"
               "  ∂V(p) = {-λ : λ Lagrange multiplier}\n"
               "\n"
               "Envelope theorem (economics):\n"
               "How optimal value changes with parameters\n"
               "\n"
               "Applications:\n"
               "- Sensitivity analysis\n"
               "- Duality theory\n"
               "- Bilevel optimization\n"
               "- Parametric programming";
    }
};

/**
 * @class ConvexDuality
 * @brief Convex optimization duality theory
 */
class ConvexDuality {
public:
    /**
     * @brief Lagrangian duality
     */
    static std::string lagrangianDuality() {
        return "Lagrangian Duality:\n"
               "\n"
               "Primal problem:\n"
               "  min f(x) s.t. g(x) ≤ 0, h(x) = 0\n"
               "\n"
               "Lagrangian:\n"
               "  L(x, λ, ν) = f(x) + ⟨λ, g(x)⟩ + ⟨ν, h(x)⟩\n"
               "\n"
               "Dual function:\n"
               "  d(λ, ν) = inf_x L(x, λ, ν)\n"
               "\n"
               "Dual problem:\n"
               "  max d(λ, ν) s.t. λ ≥ 0\n"
               "\n"
               "Weak duality:\n"
               "  d(λ, ν) ≤ p* (always holds)\n"
               "\n"
               "Strong duality:\n"
               "  d(λ*, ν*) = p* (holds under conditions)\n"
               "\n"
               "Sufficient conditions for strong duality:\n"
               "1. Slater condition: ∃x: g(x) < 0, h(x) = 0\n"
               "2. f, g convex, h affine\n"
               "\n"
               "KKT conditions (necessary and sufficient for convex):\n"
               "1. Stationarity: 0 ∈ ∂f(x*) + Σλᵢ*∂gᵢ(x*) + Σνᵢ*∂hᵢ(x*)\n"
               "2. Primal feasibility: g(x*) ≤ 0, h(x*) = 0\n"
               "3. Dual feasibility: λ* ≥ 0\n"
               "4. Complementarity: λᵢ*gᵢ(x*) = 0\n"
               "\n"
               "Interpretation of multipliers:\n"
               "- λᵢ* = sensitivity of optimal value to constraint i\n"
               "- ∂p*/∂bᵢ = -λᵢ* (envelope theorem)";
    }

    /**
     * @brief Fenchel duality
     */
    static std::string fenchelDuality() {
        return "Fenchel Duality:\n"
               "\n"
               "Primal problem:\n"
               "  p* = inf_x [f(x) + g(Ax)]\n"
               "\n"
               "Dual problem:\n"
               "  d* = -inf_y [f*(-A*y) + g*(y)]\n"
               "\n"
               "where f*, g* are Fenchel conjugates.\n"
               "\n"
               "Weak duality: d* ≤ p*\n"
               "\n"
               "Strong duality (Rockafellar):\n"
               "If f, g proper, convex, lsc and\n"
               "  0 ∈ int(A(dom f) - dom g),\n"
               "then d* = p* and dual attained.\n"
               "\n"
               "Optimality conditions:\n"
               "  x* primal optimal, y* dual optimal ⟺\n"
               "  -A*y* ∈ ∂f(x*) and y* ∈ ∂g(Ax*)\n"
               "\n"
               "Special case (Lagrangian duality):\n"
               "  f(x) = f₀(x), g(x) = I_C(x) (indicator)\n"
               "  Recovers Lagrangian formulation\n"
               "\n"
               "Applications:\n"
               "- Convex optimization\n"
               "- Image processing (total variation)\n"
               "- Machine learning (SVM dual)\n"
               "- Compressed sensing";
    }

    /**
     * @brief Dual optimization algorithms
     */
    static std::string dualAlgorithms() {
        return "Dual-Based Algorithms:\n"
               "\n"
               "1. Dual ascent:\n"
               "   - Maximize dual d(λ)\n"
               "   - λ^(k+1) = [λ^k + α∂d(λ^k)]₊\n"
               "   - Recover primal from dual\n"
               "\n"
               "2. Augmented Lagrangian:\n"
               "   - Add penalty: L_ρ(x,λ) = f(x) + ⟨λ,g(x)⟩ + (ρ/2)‖g(x)‖²\n"
               "   - Update: x^(k+1) = argmin L_ρ(x, λ^k)\n"
               "             λ^(k+1) = λ^k + ρg(x^(k+1))\n"
               "   - Converges with bounded ρ!\n"
               "\n"
               "3. ADMM (Alternating Direction Method of Multipliers):\n"
               "   - min f(x) + g(z) s.t. Ax + Bz = c\n"
               "   - x^(k+1) = argmin L_ρ(x, z^k, λ^k)\n"
               "   - z^(k+1) = argmin L_ρ(x^(k+1), z, λ^k)\n"
               "   - λ^(k+1) = λ^k + ρ(Ax^(k+1) + Bz^(k+1) - c)\n"
               "\n"
               "4. Proximal methods:\n"
               "   - Based on Moreau envelope\n"
               "   - Proximal gradient descent\n"
               "   - Douglas-Rachford splitting\n"
               "\n"
               "Advantages of dual methods:\n"
               "- May be simpler than primal\n"
               "- Decomposition possible\n"
               "- Parallel computation\n"
               "- Handle constraints naturally";
    }
};

/**
 * @class ConvexConjugate
 * @brief Properties of Fenchel conjugates
 */
class ConvexConjugate {
public:
    /**
     * @brief Conjugate calculus
     */
    static std::string conjugateCalculus() {
        return "Calculus of Fenchel Conjugates:\n"
               "\n"
               "Sum (no simple formula in general):\n"
               "  (f + g)* ≤ f* □ g* (infimal convolution)\n"
               "  Equality under conditions\n"
               "\n"
               "Scalar multiplication:\n"
               "  (αf)*(y) = αf*(y/α) for α > 0\n"
               "\n"
               "Translation:\n"
               "  (f(· - b))*(y) = f*(y) + ⟨y, b⟩\n"
               "\n"
               "Linear composition:\n"
               "  (f ∘ A)*(y) = f*(A*y)\n"
               "\n"
               "Pointwise supremum:\n"
               "  (sup_i fᵢ)* = inf_i fᵢ*\n"
               "\n"
               "Biconjugate:\n"
               "  f** = cl(conv(f)) (convex hull closure)\n"
               "  f** = f ⟺ f convex, lsc, proper\n"
               "\n"
               "Examples:\n"
               "- f(x) = (1/2)‖x‖² ⇒ f*(y) = (1/2)‖y‖²\n"
               "- f(x) = ‖x‖ ⇒ f*(y) = I_{‖y‖≤1}(y)\n"
               "- f(x) = I_C(x) ⇒ f*(y) = σ_C(y) (support function)\n"
               "- f(x) = e^x ⇒ f*(y) = y log y - y (y>0)\n"
               "\n"
               "Applications:\n"
               "- Duality transformations\n"
               "- Convexification\n"
               "- Regularization";
    }

    /**
     * @brief Support functions
     */
    static std::string supportFunctions() {
        return "Support Functions:\n"
               "\n"
               "For set C ⊆ X, support function:\n"
               "  σ_C(y) = sup_{x ∈ C} ⟨y, x⟩\n"
               "\n"
               "Properties:\n"
               "- σ_C is convex, lsc, positively homogeneous\n"
               "- σ_C = I_C* (conjugate of indicator)\n"
               "- (σ_C)* = I_{cl(conv(C))}\n"
               "\n"
               "Calculus:\n"
               "- σ_{C+D} = σ_C + σ_D\n"
               "- σ_{λC} = λσ_C for λ ≥ 0\n"
               "- σ_{conv(C)} = σ_C\n"
               "- σ_{cl(C)} = σ_C\n"
               "\n"
               "Subdifferential:\n"
               "  ∂σ_C(y) = argmax_{x∈C} ⟨y, x⟩\n"
               "\n"
               "Normal cone:\n"
               "  N_C(x) = {y : ⟨y, z-x⟩ ≤ 0 ∀z ∈ C}\n"
               "  = ∂I_C(x)\n"
               "\n"
               "Duality:\n"
               "  y ∈ N_C(x) ⟺ x ∈ ∂σ_C(y)\n"
               "  ⟺ σ_C(y) = ⟨y, x⟩\n"
               "\n"
               "Applications:\n"
               "- Constrained optimization\n"
               "- Convex geometry\n"
               "- Best approximation";
    }
};

} // namespace maths::analysis

#endif // MATHS_ANALYSIS_SUBDIFFERENTIAL_HPP
