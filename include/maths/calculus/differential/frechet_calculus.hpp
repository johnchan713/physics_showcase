#ifndef MATHS_CALCULUS_FRECHET_CALCULUS_HPP
#define MATHS_CALCULUS_FRECHET_CALCULUS_HPP

#include <functional>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>

/**
 * @file frechet_calculus.hpp
 * @brief Fréchet differential calculus and applications
 *
 * Implements:
 * - Directional and Fréchet derivatives
 * - Inverse and implicit function theorems
 * - Newton's method
 * - Legendre transform
 * - Applications to optimization (normal/tangent cones)
 * - Lagrange multipliers
 */

namespace maths::calculus {

/**
 * @class DirectionalDerivative
 * @brief Directional derivatives in Banach spaces
 */
class DirectionalDerivative {
public:
    /**
     * @brief Directional derivative definition
     */
    static std::string definition() {
        return "Directional Derivative:\n"
               "\n"
               "For f: X → Y, point x, direction v:\n"
               "\n"
               "  f'(x; v) = lim_{t↓0} [f(x + tv) - f(x)] / t\n"
               "\n"
               "when limit exists.\n"
               "\n"
               "Properties:\n"
               "- Positively homogeneous: f'(x; tv) = t·f'(x; v) for t > 0\n"
               "- Need not be linear in v\n"
               "- Need not be continuous in v\n"
               "\n"
               "Gâteaux derivative: if f'(x; v) is linear and continuous in v,\n"
               "  ∃ bounded linear operator Df(x): X → Y\n"
               "  f'(x; v) = Df(x)(v)\n"
               "\n"
               "Examples:\n"
               "- f(x) = ‖x‖: f'(0; v) = ‖v‖ (not linear!)\n"
               "- f(x) = ⟨a, x⟩: f'(x; v) = ⟨a, v⟩ (linear)\n"
               "\n"
               "Chain rule:\n"
               "If f: X → Y, g: Y → Z differentiable,\n"
               "then (g ∘ f)'(x; v) = g'(f(x); f'(x; v))";
    }

    /**
     * @brief One-sided directional derivatives
     */
    static std::string oneSided() {
        return "One-Sided Directional Derivatives:\n"
               "\n"
               "Right derivative: f'₊(x; v) = lim_{t↓0} [f(x+tv)-f(x)]/t\n"
               "Left derivative:  f'₋(x; v) = lim_{t↑0} [f(x+tv)-f(x)]/t\n"
               "\n"
               "For convex f:\n"
               "- f'₊(x; v) always exists\n"
               "- f'₊(x; v) is sublinear (positively homogeneous, subadditive)\n"
               "- Support function of ∂f(x):\n"
               "  f'(x; v) = max{⟨v*, v⟩ : v* ∈ ∂f(x)}\n"
               "\n"
               "Relationship to subdifferential:\n"
               "  v* ∈ ∂f(x) ⟺ f'(x; v) ≥ ⟨v*, v⟩ ∀v\n"
               "\n"
               "Applications:\n"
               "- Non-smooth optimization\n"
               "- Clarke generalized gradient\n"
               "- Variational inequalities";
    }
};

/**
 * @class FrechetDerivative
 * @brief Fréchet derivative and calculus
 */
class FrechetDerivative {
public:
    /**
     * @brief Fréchet derivative definition
     */
    static std::string definition() {
        return "Fréchet Derivative:\n"
               "\n"
               "f: X → Y is Fréchet differentiable at x if\n"
               "∃ bounded linear operator A: X → Y such that:\n"
               "\n"
               "  f(x + h) = f(x) + A(h) + o(‖h‖)\n"
               "\n"
               "i.e., ‖f(x+h) - f(x) - A(h)‖ / ‖h‖ → 0 as h → 0\n"
               "\n"
               "Notation: A = Df(x) = f'(x) (Fréchet derivative)\n"
               "\n"
               "Properties:\n"
               "- Unique if exists\n"
               "- Linear in h: Df(x)(αh + βk) = αDf(x)(h) + βDf(x)(k)\n"
               "- Continuous in h\n"
               "- Stronger than Gâteaux differentiability\n"
               "\n"
               "In ℝⁿ:\n"
               "  Df(x) = Jacobian matrix [∂fᵢ/∂xⱼ]\n"
               "\n"
               "For scalar f: ℝⁿ → ℝ:\n"
               "  Df(x) = ∇f(x) (gradient)\n"
               "  f(x+h) = f(x) + ⟨∇f(x), h⟩ + o(‖h‖)\n"
               "\n"
               "Fréchet ⇒ Gâteaux:\n"
               "  f'(x; v) = Df(x)(v)";
    }

    /**
     * @brief Calculus rules
     */
    static std::string calculusRules() {
        return "Fréchet Derivative Calculus Rules:\n"
               "\n"
               "Sum rule:\n"
               "  D(f + g)(x) = Df(x) + Dg(x)\n"
               "\n"
               "Product rule (for f,g: X → ℝ):\n"
               "  D(fg)(x) = f(x)Dg(x) + g(x)Df(x)\n"
               "\n"
               "Chain rule:\n"
               "  D(g ∘ f)(x) = Dg(f(x)) ∘ Df(x)\n"
               "\n"
               "In coordinates:\n"
               "  h = g ∘ f: ℝⁿ → ℝᵐ → ℝᵏ\n"
               "  Dh(x) = Dg(f(x)) · Df(x) (matrix multiplication)\n"
               "\n"
               "Leibniz rule (operator-valued):\n"
               "  D⟨A(·), x⟩(y) = ⟨DA(y), x⟩\n"
               "\n"
               "Higher derivatives:\n"
               "  D²f(x)(h, k) = D[Df(·)(k)](x)(h)\n"
               "  Symmetric if f ∈ C²\n"
               "\n"
               "Taylor's theorem:\n"
               "  f(x+h) = f(x) + Df(x)(h) + (1/2)D²f(x)(h,h) + o(‖h‖²)";
    }

    /**
     * @brief Mean value theorem
     */
    static std::string meanValueTheorem() {
        return "Mean Value Theorem (Banach Spaces):\n"
               "\n"
               "If f: X → ℝ Fréchet differentiable on line segment [x, y],\n"
               "then ∃θ ∈ (0, 1):\n"
               "\n"
               "  f(y) - f(x) = Df(x + θ(y-x))(y - x)\n"
               "\n"
               "Vector-valued version:\n"
               "If f: X → Y, then:\n"
               "\n"
               "  ‖f(y) - f(x)‖ ≤ sup_{t∈[0,1]} ‖Df(x + t(y-x))‖ · ‖y - x‖\n"
               "\n"
               "Consequences:\n"
               "1. Lipschitz continuity:\n"
               "   If ‖Df(x)‖ ≤ L for all x, then:\n"
               "   ‖f(y) - f(x)‖ ≤ L‖y - x‖\n"
               "\n"
               "2. Constant function:\n"
               "   If Df(x) = 0 for all x, then f is constant\n"
               "\n"
               "3. Convexity:\n"
               "   f convex ⟺ f(y) ≥ f(x) + Df(x)(y-x)\n"
               "\n"
               "Applications:\n"
               "- Contraction mapping arguments\n"
               "- Fixed point iterations\n"
               "- Newton's method convergence";
    }
};

/**
 * @class InverseFunctionTheorem
 * @brief Inverse and implicit function theorems
 */
class InverseFunctionTheorem {
public:
    /**
     * @brief Inverse function theorem
     */
    static std::string inverseTheorem() {
        return "Inverse Function Theorem:\n"
               "\n"
               "Let f: X → Y be C¹ map between Banach spaces.\n"
               "If Df(x₀): X → Y is bijective (invertible),\n"
               "then:\n"
               "\n"
               "1. ∃ neighborhoods U(x₀), V(f(x₀))\n"
               "2. f: U → V is bijective\n"
               "3. f⁻¹: V → U is C¹\n"
               "4. Df⁻¹(y) = [Df(f⁻¹(y))]⁻¹\n"
               "\n"
               "Proof outline:\n"
               "- Fixed point iteration for g(x) = x - [Df(x₀)]⁻¹[f(x) - y]\n"
               "- Contraction mapping theorem\n"
               "- Banach space completeness essential\n"
               "\n"
               "Finite dimensions (ℝⁿ → ℝⁿ):\n"
               "  Df(x₀) invertible ⟺ det(Df(x₀)) ≠ 0\n"
               "\n"
               "Applications:\n"
               "- Local coordinate changes\n"
               "- Diffeomorphisms\n"
               "- Solution uniqueness\n"
               "- Sensitivity analysis\n"
               "\n"
               "Example: Polar coordinates\n"
               "  f(r,θ) = (r cos θ, r sin θ)\n"
               "  Invertible when r ≠ 0, θ ≠ kπ";
    }

    /**
     * @brief Implicit function theorem
     */
    static std::string implicitTheorem() {
        return "Implicit Function Theorem:\n"
               "\n"
               "Let F: X × Y → Z be C¹, F(x₀, y₀) = 0.\n"
               "If D_y F(x₀, y₀): Y → Z is bijective, then:\n"
               "\n"
               "1. ∃ neighborhoods U(x₀), V(y₀)\n"
               "2. ∃! C¹ function g: U → V\n"
               "3. F(x, g(x)) = 0 for all x ∈ U\n"
               "4. Dg(x) = -[D_y F(x, g(x))]⁻¹ ∘ D_x F(x, g(x))\n"
               "\n"
               "Interpretation:\n"
               "Equation F(x, y) = 0 implicitly defines y = g(x)\n"
               "near (x₀, y₀).\n"
               "\n"
               "Proof via inverse function theorem:\n"
               "- Consider Φ(x, y) = (x, F(x, y))\n"
               "- DΦ invertible ⟺ D_y F invertible\n"
               "- Apply inverse theorem to Φ\n"
               "- Get Φ⁻¹(x, 0) = (x, g(x))\n"
               "\n"
               "Applications:\n"
               "- Constraint manifolds\n"
               "- Lagrange multipliers\n"
               "- Equilibrium conditions\n"
               "- Parametric optimization\n"
               "\n"
               "Example: Unit sphere x² + y² + z² = 1\n"
               "Near (0, 0, 1): z = √(1 - x² - y²)";
    }

    /**
     * @brief Rank theorem
     */
    static std::string rankTheorem() {
        return "Constant Rank Theorem:\n"
               "\n"
               "If f: X → Y has constant rank r in neighborhood of x₀,\n"
               "then ∃ local coordinates making f look like:\n"
               "\n"
               "  f(x₁, ..., xₙ) = (x₁, ..., x_r, 0, ..., 0)\n"
               "\n"
               "Consequences:\n"
               "1. Image is smooth manifold (locally)\n"
               "2. Level sets are smooth manifolds\n"
               "3. Regular value theorem\n"
               "\n"
               "Submersion (rank = dim Y):\n"
               "- Locally looks like projection\n"
               "- f⁻¹(y) is manifold\n"
               "\n"
               "Immersion (rank = dim X):\n"
               "- Locally looks like inclusion\n"
               "- Image is manifold\n"
               "\n"
               "Applications:\n"
               "- Constrained optimization\n"
               "- Differential geometry\n"
               "- Control theory";
    }
};

/**
 * @class NewtonMethod
 * @brief Newton's method and variations
 */
class NewtonMethod {
public:
    /**
     * @brief Newton's method for equations
     */
    static std::string method() {
        return "Newton's Method:\n"
               "\n"
               "To solve F(x) = 0 where F: X → Y:\n"
               "\n"
               "Iteration:\n"
               "  x_{k+1} = x_k - [DF(x_k)]⁻¹ F(x_k)\n"
               "\n"
               "Equivalently:\n"
               "  Solve: DF(x_k)(x_{k+1} - x_k) = -F(x_k)\n"
               "\n"
               "Convergence (Kantorovich theorem):\n"
               "If ‖[DF(x₀)]⁻¹‖ ≤ β,\n"
               "   ‖[DF(x₀)]⁻¹ F(x₀)‖ ≤ η,\n"
               "   ‖D²F‖ ≤ K,\n"
               "   h = βηK ≤ 1/2,\n"
               "\n"
               "Then:\n"
               "- Sequence converges to solution x*\n"
               "- Quadratic convergence: ‖x_{k+1} - x*‖ ≤ C‖x_k - x*‖²\n"
               "\n"
               "Conditions roughly mean:\n"
               "- Good initial guess (small η)\n"
               "- Well-conditioned Jacobian (small β)\n"
               "- Smooth nonlinearity (bounded K)\n"
               "\n"
               "Variants:\n"
               "- Quasi-Newton: approximate DF(x_k)\n"
               "- Gauss-Newton: for least squares\n"
               "- Damped Newton: x_{k+1} = x_k - α_k[DF(x_k)]⁻¹F(x_k)";
    }

    /**
     * @brief Newton for optimization
     */
    static std::string optimization() {
        return "Newton's Method for Optimization:\n"
               "\n"
               "To minimize f: ℝⁿ → ℝ:\n"
               "\n"
               "Iteration:\n"
               "  x_{k+1} = x_k - [D²f(x_k)]⁻¹ ∇f(x_k)\n"
               "\n"
               "Rationale:\n"
               "- Minimize quadratic model:\n"
               "  q(x) = f(x_k) + ⟨∇f(x_k), x-x_k⟩ + (1/2)⟨D²f(x_k)(x-x_k), x-x_k⟩\n"
               "- Critical point: ∇q = 0\n"
               "  ⇒ ∇f(x_k) + D²f(x_k)(x - x_k) = 0\n"
               "\n"
               "Convergence:\n"
               "If f ∈ C², ∇²f positive definite near x*,\n"
               "then quadratic convergence from good initial point.\n"
               "\n"
               "Advantages:\n"
               "- Quadratic convergence (very fast)\n"
               "- Affine invariant\n"
               "\n"
               "Disadvantages:\n"
               "- Requires Hessian (expensive)\n"
               "- May not descend (need line search)\n"
               "- Hessian may be singular\n"
               "\n"
               "Practical variants:\n"
               "- Trust region Newton\n"
               "- Modified Newton (regularize Hessian)\n"
               "- BFGS (quasi-Newton, approximate Hessian)";
    }
};

/**
 * @class LegendreTransform
 * @brief Legendre and Legendre-Fenchel transforms
 */
class LegendreTransform {
public:
    /**
     * @brief Classical Legendre transform
     */
    static std::string classical() {
        return "Classical Legendre Transform:\n"
               "\n"
               "For strictly convex C² function f: ℝⁿ → ℝ:\n"
               "\n"
               "Transform variables from x to p = ∇f(x):\n"
               "  g(p) = ⟨p, x(p)⟩ - f(x(p))\n"
               "\n"
               "where x(p) solves ∇f(x) = p.\n"
               "\n"
               "Properties:\n"
               "- Involution: (f*)* = f\n"
               "- Gradient duality: ∇g(p) = x\n"
               "- Hessian: D²g(p) = [D²f(x)]⁻¹\n"
               "\n"
               "Example: f(x) = (1/2)⟨Ax, x⟩ (A positive definite)\n"
               "  ∇f(x) = Ax ⇒ x = A⁻¹p\n"
               "  g(p) = (1/2)⟨A⁻¹p, p⟩\n"
               "\n"
               "Applications:\n"
               "- Classical mechanics (Lagrangian ↔ Hamiltonian)\n"
               "- Thermodynamics (energy ↔ entropy)\n"
               "- Duality in optimization\n"
               "\n"
               "Hamiltonian mechanics:\n"
               "  L(q, q̇) = kinetic - potential\n"
               "  H(q, p) = ⟨p, q̇⟩ - L(q, q̇) where p = ∂L/∂q̇";
    }

    /**
     * @brief Connection to Fenchel conjugate
     */
    static std::string fenchel() {
        return "Legendre vs Fenchel Transform:\n"
               "\n"
               "Fenchel conjugate (general):\n"
               "  f*(p) = sup_{x} [⟨p, x⟩ - f(x)]\n"
               "\n"
               "Legendre (smooth, strictly convex):\n"
               "  f*(p) = ⟨p, x*⟩ - f(x*) where ∇f(x*) = p\n"
               "\n"
               "Relationship:\n"
               "- Fenchel generalizes Legendre\n"
               "- Works for non-smooth, non-strictly convex f\n"
               "- Supremum may not have unique maximizer\n"
               "\n"
               "Key differences:\n"
               "Legendre: requires ∇f invertible\n"
               "Fenchel: always defined (possibly +∞)\n"
               "\n"
               "When f is strictly convex C¹:\n"
               "  Legendre = Fenchel\n"
               "\n"
               "Examples where they differ:\n"
               "- f(x) = |x|: no Legendre transform (not differentiable)\n"
               "  but f*(p) = I_{[-1,1]}(p)\n"
               "- f(x) = x²: Legendre exists\n"
               "  f*(p) = p²/4\n"
               "\n"
               "Both satisfy:\n"
               "  f(x) + f*(p) ≥ ⟨p, x⟩ (Fenchel-Young)\n"
               "  with equality ⟺ p ∈ ∂f(x) ⟺ x ∈ ∂f*(p)";
    }
};

} // namespace maths::calculus

#endif // MATHS_CALCULUS_FRECHET_CALCULUS_HPP
