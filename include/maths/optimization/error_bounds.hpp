#ifndef MATHS_OPTIMIZATION_ERROR_BOUNDS_HPP
#define MATHS_OPTIMIZATION_ERROR_BOUNDS_HPP

#include <functional>
#include <cmath>
#include <string>
#include <vector>

/**
 * @file error_bounds.hpp
 * @brief Error bounds, penalization, and metric estimates
 *
 * Implements:
 * - Error bounds and metric regularity
 * - Decrease principles
 * - Penalization methods
 * - Augmented Lagrangians
 * - Robust optimization
 * - Stabilized infima
 */

namespace maths::optimization {

/**
 * @class ErrorBounds
 * @brief Error bounds and metric regularity
 */
class ErrorBounds {
public:
    /**
     * @brief Definition of error bounds
     */
    static std::string definition() {
        return "Error Bounds:\n"
               "\n"
               "Function f: X → ℝ has error bound at S = {x : f(x) = 0} if:\n"
               "  ∃c, δ > 0: d(x, S) ≤ c|f(x)| for all x with d(x, S) < δ\n"
               "\n"
               "(distance to solution set controlled by residual)\n"
               "\n"
               "Global error bound: inequality holds for all x ∈ X\n"
               "\n"
               "Hölder error bound: d(x, S) ≤ c|f(x)|^α for α ∈ (0, 1]\n"
               "- α = 1: Lipschitz error bound\n"
               "- α < 1: slower convergence\n"
               "\n"
               "Applications:\n"
               "- Convergence rate analysis\n"
               "- Stopping criteria for algorithms\n"
               "- Condition number estimates\n"
               "- Perturbation analysis";
    }

    /**
     * @brief Hoffman's error bound
     */
    static std::string hoffman() {
        return "Hoffman's Error Bound:\n"
               "\n"
               "For linear system Ax = b, x ≥ 0:\n"
               "  d(x, S) ≤ c‖(Ax - b)₊‖ + c‖x₋‖\n"
               "\n"
               "where S = {x : Ax = b, x ≥ 0},\n"
               "      (·)₊ = max(·, 0), (·)₋ = max(-·, 0)\n"
               "\n"
               "Constant c depends only on A (not b)!\n"
               "\n"
               "Extensions:\n"
               "- Polyhedra: S = {x : Ax ≤ b}\n"
               "- Convex inequalities\n"
               "- Variational inequalities\n"
               "\n"
               "Key insight: Error bound holds uniformly\n"
               "for all right-hand sides b.\n"
               "\n"
               "Applications:\n"
               "- Linear programming sensitivity\n"
               "- Convergence of simplex method\n"
               "- Complementarity problems";
    }

    /**
     * @brief Łojasiewicz inequality
     */
    static std::string lojasiewicz() {
        return "Łojasiewicz Inequality:\n"
               "\n"
               "For analytic f: ℝⁿ → ℝ with critical point x*:\n"
               "  ∃c, δ, α ∈ [1/2, 1): |f(x) - f(x*)|^α ≤ c‖∇f(x)‖\n"
               "for x near x*.\n"
               "\n"
               "(gradient norm controls function value difference)\n"
               "\n"
               "Consequences:\n"
               "1. Gradient descent converges to critical point\n"
               "2. Finite length property: ∫₀^∞ ‖ẋ(t)‖ dt < ∞\n"
               "3. Convergence rate estimates\n"
               "\n"
               "KŁ (Kurdyka-Łojasiewicz) property:\n"
               "Generalization to non-smooth functions:\n"
               "  φ'(f(x) - f(x*)) · d(0, ∂f(x)) ≥ 1\n"
               "for suitable desingularizing function φ.\n"
               "\n"
               "Applications:\n"
               "- Convergence of descent methods\n"
               "- Alternating minimization\n"
               "- Proximal algorithms";
    }

    /**
     * @brief Metric regularity
     */
    static std::string metricRegularity() {
        return "Metric Regularity:\n"
               "\n"
               "Set-valued map F: X ⇉ Y is metrically regular at (x̄, ȳ) ∈ graph(F)\n"
               "if ∃c, δ > 0:\n"
               "  d(x, F⁻¹(y)) ≤ c · d(y, F(x))\n"
               "for x near x̄, y near ȳ.\n"
               "\n"
               "(inverse approximately Lipschitz)\n"
               "\n"
               "Equivalent to:\n"
               "- Aubin/Lipschitz property of F⁻¹\n"
               "- Linear convergence of Newton's method\n"
               "- Robinson's constraint qualification\n"
               "\n"
               "Ioffe's Theorem: F metrically regular ⟺\n"
               "  F⁻¹ has Aubin property\n"
               "\n"
               "Applications:\n"
               "- Implicit function theorem generalizations\n"
               "- Newton method convergence\n"
               "- Sensitivity analysis\n"
               "- Regularity of solution maps";
    }
};

/**
 * @class DecreasePrinciple
 * @brief Sufficient decrease conditions
 */
class DecreasePrinciple {
public:
    /**
     * @brief Armijo rule
     */
    static std::string armijo() {
        return "Armijo Rule (Sufficient Decrease):\n"
               "\n"
               "In line search for min f(x), accept step αₖdₖ if:\n"
               "  f(xₖ + αₖdₖ) ≤ f(xₖ) + c₁αₖ⟨∇f(xₖ), dₖ⟩\n"
               "\n"
               "with c₁ ∈ (0, 1) (typically 10⁻⁴)\n"
               "\n"
               "Interpretation:\n"
               "- Actual decrease ≥ c₁ × predicted decrease\n"
               "- Prevents too small steps\n"
               "- Ensures progress at each iteration\n"
               "\n"
               "Backtracking: Start α = 1, reduce by factor until satisfied\n"
               "\n"
               "Theorem: If ∇f Lipschitz and dₖ descent direction,\n"
               "Armijo rule satisfied for sufficiently small α.";
    }

    /**
     * @brief Wolfe conditions
     */
    static std::string wolfe() {
        return "Wolfe Conditions:\n"
               "\n"
               "Accept step if both:\n"
               "1. Armijo: f(xₖ + αdₖ) ≤ f(xₖ) + c₁α⟨∇f, dₖ⟩\n"
               "2. Curvature: ⟨∇f(xₖ + αdₖ), dₖ⟩ ≥ c₂⟨∇f(xₖ), dₖ⟩\n"
               "\n"
               "with 0 < c₁ < c₂ < 1 (typical: c₁ = 10⁻⁴, c₂ = 0.9)\n"
               "\n"
               "Strong Wolfe: Replace curvature condition with\n"
               "  |⟨∇f(xₖ + αdₖ), dₖ⟩| ≤ c₂|⟨∇f(xₖ), dₖ⟩|\n"
               "\n"
               "Ensures:\n"
               "- Sufficient decrease (condition 1)\n"
               "- Not too small step (condition 2)\n"
               "- Quasi-Newton methods converge\n"
               "\n"
               "Theorem: If ∇f Lipschitz, interval of acceptable\n"
               "α satisfying Wolfe conditions is non-empty.";
    }

    /**
     * @brief Goldstein conditions
     */
    static std::string goldstein() {
        return "Goldstein Conditions:\n"
               "\n"
               "Accept step α if:\n"
               "  f(xₖ) + (1-c)α⟨∇f, dₖ⟩ ≤ f(xₖ + αdₖ) ≤ f(xₖ) + cα⟨∇f, dₖ⟩\n"
               "\n"
               "with c ∈ (0, 1/2) (typically c = 0.25)\n"
               "\n"
               "Interpretation:\n"
               "- Lower bound: prevents too large steps\n"
               "- Upper bound: prevents too small steps\n"
               "- Sandwiches actual decrease between bounds\n"
               "\n"
               "Comparison with Wolfe:\n"
               "- Goldstein: bounds on function value\n"
               "- Wolfe: bounds on derivative\n"
               "- Wolfe more commonly used\n"
               "\n"
               "Both ensure global convergence of descent methods.";
    }
};

/**
 * @class Penalization
 * @brief Penalization and augmented Lagrangian methods
 */
class Penalization {
public:
    /**
     * @brief Exact penalties
     */
    static std::string exactPenalty() {
        return "Exact Penalty Methods:\n"
               "\n"
               "Original problem:\n"
               "  min f(x) subject to g(x) = 0, h(x) ≤ 0\n"
               "\n"
               "Penalized problem:\n"
               "  min f(x) + ρ(‖g(x)‖ + ‖h(x)₊‖)\n"
               "\n"
               "Exact penalty: ∃ρ* such that for ρ > ρ*,\n"
               "minimizer of penalized = solution of original\n"
               "\n"
               "ℓ₁ exact penalty:\n"
               "  min f(x) + ρ∑|gᵢ(x)| + ρ∑max(0, hⱼ(x))\n"
               "\n"
               "Advantages:\n"
               "- Single penalty parameter\n"
               "- Exact for finite ρ\n"
               "\n"
               "Disadvantages:\n"
               "- Non-smooth (even if f, g, h smooth)\n"
               "- May require large ρ\n"
               "\n"
               "Theorem: If constraint qualification holds and\n"
               "‖λ*‖ < ρ where λ* = Lagrange multipliers,\n"
               "then penalty is exact.";
    }

    /**
     * @brief Quadratic penalty
     */
    static std::string quadraticPenalty() {
        return "Quadratic Penalty Method:\n"
               "\n"
               "Penalized problem:\n"
               "  min f(x) + (ρ/2)(‖g(x)‖² + ‖h(x)₊‖²)\n"
               "\n"
               "Algorithm:\n"
               "1. Choose ρ₀ > 0, x₀\n"
               "2. For k = 0, 1, 2, ...:\n"
               "   a) xₖ₊₁ = argmin[f(x) + (ρₖ/2)(‖g(x)‖² + ‖h(x)₊‖²)]\n"
               "   b) ρₖ₊₁ = βρₖ (β > 1)\n"
               "\n"
               "Properties:\n"
               "- Smooth (if f, g, h smooth)\n"
               "- NOT exact: need ρₖ → ∞\n"
               "- xₖ → x* but ρₖ → ∞ (ill-conditioning)\n"
               "\n"
               "Convergence:\n"
               "- If ρₖ → ∞ sufficiently fast, xₖ → x*\n"
               "- Constraint violation ‖g(xₖ)‖ → 0\n"
               "\n"
               "Issue: Large ρ causes ill-conditioning\n"
               "(Hessian eigenvalues ~ O(ρ))";
    }

    /**
     * @brief Augmented Lagrangian
     */
    static std::string augmentedLagrangian() {
        return "Augmented Lagrangian (Method of Multipliers):\n"
               "\n"
               "For equality constraints g(x) = 0:\n"
               "  L_ρ(x, λ) = f(x) + ⟨λ, g(x)⟩ + (ρ/2)‖g(x)‖²\n"
               "\n"
               "Algorithm:\n"
               "1. Initialize x₀, λ₀, ρ₀\n"
               "2. For k = 0, 1, 2, ...:\n"
               "   a) xₖ₊₁ ≈ argmin_x L_ρₖ(x, λₖ)\n"
               "   b) λₖ₊₁ = λₖ + ρₖg(xₖ₊₁)\n"
               "   c) ρₖ₊₁ = ρₖ (or increase if needed)\n"
               "\n"
               "Advantages over quadratic penalty:\n"
               "- Converges with BOUNDED ρₖ\n"
               "- Better conditioning\n"
               "- λₖ estimates Lagrange multipliers\n"
               "\n"
               "Theorem (Rockafellar): If strong second-order\n"
               "conditions hold, method converges for fixed ρ.\n"
               "\n"
               "Extensions:\n"
               "- ADMM (Alternating Direction Method of Multipliers)\n"
               "- Douglas-Rachford splitting\n"
               "- Proximal methods";
    }

    /**
     * @brief Barrier methods
     */
    static std::string barrierMethods() {
        return "Barrier Methods (Interior Point):\n"
               "\n"
               "For h(x) ≤ 0:\n"
               "  min f(x) + (1/μ)B(x)\n"
               "\n"
               "where B(x) barrier function, e.g.:\n"
               "- Logarithmic: B(x) = -∑log(-hᵢ(x))\n"
               "- Inverse: B(x) = -∑1/hᵢ(x)\n"
               "\n"
               "Properties:\n"
               "- Keeps iterates strictly feasible\n"
               "- Barrier → ∞ at boundary\n"
               "- μ → 0: barrier influence decreases\n"
               "\n"
               "Central path: {x*(μ) : μ > 0}\n"
               "where x*(μ) minimizes f + (1/μ)B\n"
               "\n"
               "Path-following algorithms:\n"
               "- Follow central path as μ → 0\n"
               "- Polynomial-time complexity\n"
               "- Modern linear/convex programming solvers\n"
               "\n"
               "Primal-dual methods: Track both primal and dual.";
    }
};

/**
 * @class RobustOptimization
 * @brief Robust infima and stabilization
 */
class RobustOptimization {
public:
    /**
     * @brief Robust optimization
     */
    static std::string robust() {
        return "Robust Optimization:\n"
               "\n"
               "Original uncertain problem:\n"
               "  min_{x ∈ X} max_{u ∈ U} f(x, u)\n"
               "\n"
               "(minimize worst-case over uncertainty U)\n"
               "\n"
               "Robust counterpart: Reformulation that can be solved\n"
               "\n"
               "Example (uncertain linear program):\n"
               "  min cᵀx s.t. aᵀx ≤ b for all a ∈ U\n"
               "\n"
               "If U = {a : ‖a - ā‖ ≤ Γ} (ellipsoidal):\n"
               "  min cᵀx s.t. āᵀx + Γ‖x‖ ≤ b (SOCP)\n"
               "\n"
               "Approaches:\n"
               "- Worst-case optimization\n"
               "- Distributional robustness\n"
               "- Adaptive optimization\n"
               "\n"
               "Robustness-performance tradeoff:\n"
               "More robust → more conservative → higher cost";
    }

    /**
     * @brief Stabilized infima
     */
    static std::string stabilized() {
        return "Stabilized Infima:\n"
               "\n"
               "Problem: inf_x f(x) may not be attained\n"
               "\n"
               "Stabilization: Add regularization term\n"
               "  val_ε = inf_x [f(x) + εr(x)]\n"
               "\n"
               "where r(x) coercive (e.g., ‖x‖²)\n"
               "\n"
               "Properties:\n"
               "- val_ε attained for ε > 0\n"
               "- val_ε ↓ inf f as ε ↓ 0\n"
               "- x_ε bounded\n"
               "\n"
               "Tykhonov regularization:\n"
               "  min_x [‖Ax - b‖² + ε‖x‖²]\n"
               "\n"
               "Applications:\n"
               "- Ill-posed problems\n"
               "- Inverse problems\n"
               "- Machine learning (regularization)\n"
               "\n"
               "Trade-off: ε large (stable, inaccurate)\n"
               "vs ε small (unstable, accurate)";
    }

    /**
     * @brief Tikhonov regularization
     */
    static std::string tikhonov() {
        return "Tikhonov Regularization:\n"
               "\n"
               "Ill-posed problem: Ax = b (A may not have bounded inverse)\n"
               "\n"
               "Regularized problem:\n"
               "  min_x [‖Ax - b‖² + α‖x‖²]\n"
               "\n"
               "Solution: x_α = (AᵀA + αI)⁻¹Aᵀb\n"
               "\n"
               "Properties:\n"
               "- AᵀA + αI always invertible for α > 0\n"
               "- x_α continuous in b\n"
               "- As α → 0: x_α → A⁺b (pseudoinverse solution)\n"
               "\n"
               "Choosing α (regularization parameter):\n"
               "- L-curve method\n"
               "- Cross-validation\n"
               "- Discrepancy principle: ‖Ax_α - b‖ ≈ noise level\n"
               "\n"
               "Generalizations:\n"
               "- Different norms/regularizers\n"
               "- Iterative regularization\n"
               "- Bayesian interpretation (prior on x)";
    }
};

} // namespace maths::optimization

#endif // MATHS_OPTIMIZATION_ERROR_BOUNDS_HPP
