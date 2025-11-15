#ifndef MATHS_ANALYSIS_ADVANCED_SUBDIFFERENTIALS_HPP
#define MATHS_ANALYSIS_ADVANCED_SUBDIFFERENTIALS_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <algorithm>
#include <set>
#include <utility>

/**
 * @file advanced_subdifferentials.hpp
 * @brief Advanced subdifferential calculus: Clarke, limiting, graded subdifferentials
 *
 * Computational implementations of Chapters 5-7:
 * - Clarke Jacobian and subdifferential in finite dimensions
 * - Circa-normal and circa-tangent cones
 * - Limiting subdifferentials and coderivatives
 * - Graded subdifferentials
 */

namespace maths::analysis {

// Type aliases for clarity
using Point = std::vector<double>;
using Direction = std::vector<double>;
using Subgradient = std::vector<double>;
using SubdifferentialSet = std::vector<Subgradient>;
using Matrix = std::vector<std::vector<double>>;

/**
 * @class ClarkeJacobian
 * @brief Clarke Jacobian for locally Lipschitz vector functions (5.1.2-5.1.3)
 *
 * For F: ℝⁿ → ℝᵐ locally Lipschitz, the Clarke Jacobian is:
 * ∂_C F(x) = conv{ lim J_F(x_k) : x_k → x, F differentiable at x_k }
 */
class ClarkeJacobian {
public:
    /**
     * @brief Compute Clarke Jacobian numerically
     *
     * Uses finite differences to approximate generalized Jacobian
     *
     * @param F Vector function ℝⁿ → ℝᵐ
     * @param x Point of evaluation
     * @param h Step size for finite differences
     * @param n_samples Number of sampling directions
     * @return Set of Jacobian matrices (as vectors of matrices)
     */
    static std::vector<Matrix> computeClarkeJacobian(
        std::function<Point(const Point&)> F,
        const Point& x,
        double h = 1e-6,
        int n_samples = 50) {

        size_t n = x.size();
        size_t m = F(x).size();
        std::vector<Matrix> jacobians;

        // Sample Jacobians from nearby points
        for (int s = 0; s < n_samples; ++s) {
            // Generate random nearby point
            Point x_nearby = x;
            for (size_t i = 0; i < n; ++i) {
                double offset = h * (2.0 * (rand() / double(RAND_MAX)) - 1.0);
                x_nearby[i] += offset;
            }

            // Compute Jacobian at nearby point via finite differences
            Matrix J(m, std::vector<double>(n));
            Point fx = F(x_nearby);

            for (size_t j = 0; j < n; ++j) {
                Point x_perturb = x_nearby;
                x_perturb[j] += h;
                Point fx_perturb = F(x_perturb);

                for (size_t i = 0; i < m; ++i) {
                    J[i][j] = (fx_perturb[i] - fx[i]) / h;
                }
            }

            jacobians.push_back(J);
        }

        return jacobians;
    }

    /**
     * @brief Compute Clarke directional derivative for scalar function
     *
     * f°(x; v) = limsup_{y→x, t↓0} [f(y+tv) - f(y)] / t
     *
     * @param f Scalar function
     * @param x Point
     * @param v Direction
     * @param h Step size
     * @return Clarke directional derivative
     */
    static double clarkeDirectionalDerivative(
        std::function<double(const Point&)> f,
        const Point& x,
        const Direction& v,
        double h = 1e-6) {

        double max_diff = -std::numeric_limits<double>::infinity();

        // Sample nearby points and small t values
        for (int i = 0; i < 20; ++i) {
            Point y = x;
            for (size_t j = 0; j < x.size(); ++j) {
                y[j] += h * (rand() / double(RAND_MAX) - 0.5) * 2.0;
            }

            for (double t = h; t > h/10; t *= 0.5) {
                Point y_plus_tv = y;
                for (size_t j = 0; j < v.size(); ++j) {
                    y_plus_tv[j] += t * v[j];
                }

                double diff = (f(y_plus_tv) - f(y)) / t;
                max_diff = std::max(max_diff, diff);
            }
        }

        return max_diff;
    }

    /**
     * @brief Compute Clarke subdifferential ∂_C f(x) for scalar function
     *
     * ∂_C f(x) = { ξ : f°(x; v) ≥ ⟨ξ, v⟩ for all v }
     *
     * @param f Scalar function
     * @param x Point
     * @param h Step size
     * @param n_directions Number of test directions
     * @return Set of subgradients
     */
    static SubdifferentialSet computeClarkeSubdifferential(
        std::function<double(const Point&)> f,
        const Point& x,
        double h = 1e-6,
        int n_directions = 30) {

        size_t n = x.size();
        SubdifferentialSet subdiff;

        // Compute gradients at nearby points where function is "smooth"
        for (int s = 0; s < n_directions; ++s) {
            Point x_nearby = x;
            for (size_t i = 0; i < n; ++i) {
                x_nearby[i] += h * (2.0 * (rand() / double(RAND_MAX)) - 1.0);
            }

            // Compute gradient via finite differences
            Subgradient grad(n);
            double fx = f(x_nearby);

            for (size_t j = 0; j < n; ++j) {
                Point x_perturb = x_nearby;
                x_perturb[j] += h;
                grad[j] = (f(x_perturb) - fx) / h;
            }

            subdiff.push_back(grad);
        }

        return subdiff;
    }
};

/**
 * @class NormalTangentCones
 * @brief Circa-normal and circa-tangent cones (5.2)
 *
 * Normal cone: N_C(x) = { v : ⟨v, y-x⟩ ≤ 0 for all y ∈ C }
 * Tangent cone: T_C(x) = { v : ∃t_k↓0, v_k→v with x + t_k v_k ∈ C }
 */
class NormalTangentCones {
public:
    /**
     * @brief Compute Fréchet normal cone to set C at x
     *
     * N̂_C(x) = { v : limsup_{y→x, y∈C} ⟨v, y-x⟩/‖y-x‖ ≤ 0 }
     *
     * @param indicator Indicator function I_C (0 if in C, infinity otherwise)
     * @param x Point in C
     * @param h Sampling radius
     * @param n_samples Number of samples
     * @return Set of normal vectors
     */
    static SubdifferentialSet computeFrechetNormal(
        std::function<bool(const Point&)> inSet,
        const Point& x,
        double h = 1e-4,
        int n_samples = 50) {

        size_t n = x.size();
        SubdifferentialSet normals;

        // Sample directions and check normal cone condition
        for (int s = 0; s < n_samples; ++s) {
            Direction v(n);
            double norm_v = 0.0;

            for (size_t i = 0; i < n; ++i) {
                v[i] = 2.0 * (rand() / double(RAND_MAX)) - 1.0;
                norm_v += v[i] * v[i];
            }
            norm_v = std::sqrt(norm_v);

            // Normalize
            for (size_t i = 0; i < n; ++i) {
                v[i] /= norm_v;
            }

            // Check if v is in normal cone
            bool is_normal = true;
            for (int t = 0; t < 20; ++t) {
                Point y = x;
                for (size_t i = 0; i < n; ++i) {
                    y[i] += h * (2.0 * (rand() / double(RAND_MAX)) - 1.0);
                }

                if (!inSet(y)) continue;

                double inner_prod = 0.0;
                double dist = 0.0;
                for (size_t i = 0; i < n; ++i) {
                    inner_prod += v[i] * (y[i] - x[i]);
                    dist += (y[i] - x[i]) * (y[i] - x[i]);
                }
                dist = std::sqrt(dist);

                if (dist > 1e-12 && inner_prod / dist > 1e-6) {
                    is_normal = false;
                    break;
                }
            }

            if (is_normal) {
                normals.push_back(v);
            }
        }

        return normals;
    }

    /**
     * @brief Compute tangent cone to set C at x
     *
     * T_C(x) = cl{ λ(y - x) : y ∈ C, λ ≥ 0 }
     *
     * @param inSet Characteristic function of set C
     * @param x Point in C
     * @param h Sampling radius
     * @param n_samples Number of samples
     * @return Set of tangent directions
     */
    static SubdifferentialSet computeTangentCone(
        std::function<bool(const Point&)> inSet,
        const Point& x,
        double h = 1e-4,
        int n_samples = 50) {

        size_t n = x.size();
        SubdifferentialSet tangents;

        // Sample points in C near x
        for (int s = 0; s < n_samples; ++s) {
            Point y = x;
            for (size_t i = 0; i < n; ++i) {
                y[i] += h * (2.0 * (rand() / double(RAND_MAX)) - 1.0);
            }

            if (!inSet(y)) continue;

            // Compute direction y - x
            Direction v(n);
            double norm_v = 0.0;
            for (size_t i = 0; i < n; ++i) {
                v[i] = y[i] - x[i];
                norm_v += v[i] * v[i];
            }
            norm_v = std::sqrt(norm_v);

            if (norm_v < 1e-12) continue;

            // Normalize
            for (size_t i = 0; i < n; ++i) {
                v[i] /= norm_v;
            }

            tangents.push_back(v);
        }

        return tangents;
    }

    /**
     * @brief Compute projection onto convex set
     *
     * proj_C(x) = argmin_{y∈C} ‖y - x‖
     *
     * @param inSet Characteristic function
     * @param distance Distance function to set
     * @param x Point to project
     * @param lr Learning rate for gradient descent
     * @param max_iter Maximum iterations
     * @return Projection of x onto C
     */
    static Point projectOntoSet(
        std::function<bool(const Point&)> inSet,
        std::function<double(const Point&)> distance,
        const Point& x,
        double lr = 0.01,
        int max_iter = 1000) {

        Point y = x;
        double h = 1e-6;

        for (int iter = 0; iter < max_iter; ++iter) {
            // Compute gradient of distance function
            Direction grad(x.size());
            double dist_y = distance(y);

            for (size_t i = 0; i < x.size(); ++i) {
                Point y_perturb = y;
                y_perturb[i] += h;
                grad[i] = (distance(y_perturb) - dist_y) / h;
            }

            // Gradient descent step
            for (size_t i = 0; i < x.size(); ++i) {
                y[i] -= lr * grad[i];
            }

            // Check convergence
            double norm_grad = 0.0;
            for (double g : grad) {
                norm_grad += g * g;
            }
            if (norm_grad < 1e-10) break;
        }

        return y;
    }
};

/**
 * @class LimitingSubdifferential
 * @brief Limiting (Mordukhovich) subdifferential (6.1)
 *
 * ∂_L f(x) = Limsup_{x'→x} ∂̂ f(x')
 */
class LimitingSubdifferential {
public:
    /**
     * @brief Compute limiting subdifferential
     *
     * Sequential closure of Fréchet subdifferentials
     *
     * @param f Function
     * @param x Point
     * @param h Neighborhood radius
     * @param n_samples Number of nearby points
     * @return Limiting subdifferential
     */
    static SubdifferentialSet computeLimiting(
        std::function<double(const Point&)> f,
        const Point& x,
        double h = 1e-5,
        int n_samples = 40) {

        size_t n = x.size();
        SubdifferentialSet limiting;

        // Sample Fréchet subgradients from sequence x_k → x
        for (int k = 0; k < n_samples; ++k) {
            double radius = h * (1.0 - k / double(n_samples));
            Point x_k = x;

            for (size_t i = 0; i < n; ++i) {
                x_k[i] += radius * (2.0 * (rand() / double(RAND_MAX)) - 1.0);
            }

            // Compute Fréchet subgradient at x_k
            Subgradient grad(n);
            double fx_k = f(x_k);
            double eps = 1e-6;

            for (size_t j = 0; j < n; ++j) {
                Point x_perturb = x_k;
                x_perturb[j] += eps;
                grad[j] = (f(x_perturb) - fx_k) / eps;
            }

            limiting.push_back(grad);
        }

        return limiting;
    }

    /**
     * @brief Compute limiting normal cone
     *
     * N_L C(x) = Limsup_{x'→x, x'∈C} N̂_C(x')
     *
     * @param inSet Characteristic function
     * @param x Point
     * @param h Sampling radius
     * @param n_samples Number of samples
     * @return Limiting normal cone
     */
    static SubdifferentialSet computeLimitingNormal(
        std::function<bool(const Point&)> inSet,
        const Point& x,
        double h = 1e-5,
        int n_samples = 40) {

        SubdifferentialSet limiting_normals;

        // Sample Fréchet normals from x_k → x
        for (int k = 0; k < n_samples; ++k) {
            double radius = h * (1.0 - k / double(n_samples));
            Point x_k = x;

            for (size_t i = 0; i < x.size(); ++i) {
                x_k[i] += radius * (2.0 * (rand() / double(RAND_MAX)) - 1.0);
            }

            if (!inSet(x_k)) continue;

            // Compute Fréchet normal at x_k
            auto frechet_normals = NormalTangentCones::computeFrechetNormal(
                inSet, x_k, radius * 0.1, 5);

            for (const auto& normal : frechet_normals) {
                limiting_normals.push_back(normal);
            }
        }

        return limiting_normals;
    }

    /**
     * @brief Check metric regularity criterion
     *
     * F is metrically regular at (x̄, ȳ) iff
     * D*_L F(x̄, ȳ)(0) = {0}
     *
     * @param F Multifunction
     * @param x Point in domain
     * @param y Point in image
     * @param tolerance Numerical tolerance
     * @return true if metrically regular
     */
    static bool isMetricallyRegular(
        std::function<Point(const Point&)> F,
        const Point& x,
        const Point& y,
        double tolerance = 1e-6) {

        // Compute coderivative D*F(x,y)(0)
        // Simplified check: verify invertibility of Jacobian
        size_t n = x.size();
        Matrix J(n, std::vector<double>(n));
        Point fx = F(x);
        double h = 1e-6;

        for (size_t j = 0; j < n; ++j) {
            Point x_perturb = x;
            x_perturb[j] += h;
            Point fx_perturb = F(x_perturb);

            for (size_t i = 0; i < n; ++i) {
                J[i][j] = (fx_perturb[i] - fx[i]) / h;
            }
        }

        // Check if determinant is non-zero (simplified criterion)
        double det = computeDeterminant(J);
        return std::abs(det) > tolerance;
    }

private:
    static double computeDeterminant(const Matrix& M) {
        size_t n = M.size();
        if (n == 1) return M[0][0];
        if (n == 2) return M[0][0] * M[1][1] - M[0][1] * M[1][0];

        // For larger matrices, use simplified approximation
        double det = 1.0;
        for (size_t i = 0; i < n; ++i) {
            det *= M[i][i];
        }
        return det;
    }
};

/**
 * @class Coderivatives
 * @brief Coderivatives of set-valued mappings (6.1.2, 6.3)
 */
class Coderivatives {
public:
    /**
     * @brief Compute coderivative D*F(x̄, ȳ)
     *
     * For F: X ⇉ Y, the coderivative D*F(x̄, ȳ): Y* → X* is:
     * D*F(x̄, ȳ)(y*) = { x* : (x*, -y*) ∈ N_L gph(F)(x̄, ȳ) }
     *
     * @param F Set-valued map (returns vector of possible images)
     * @param x Point in domain
     * @param y Point in F(x)
     * @param y_star Covector
     * @param h Sampling radius
     * @return Set of x* in coderivative
     */
    static SubdifferentialSet computeCoderivative(
        std::function<std::vector<Point>(const Point&)> F,
        const Point& x,
        const Point& y,
        const Direction& y_star,
        double h = 1e-5) {

        size_t n = x.size();
        size_t m = y.size();
        SubdifferentialSet coderivative;

        // Sample normals to graph of F
        for (int s = 0; s < 30; ++s) {
            Point x_nearby = x;
            for (size_t i = 0; i < n; ++i) {
                x_nearby[i] += h * (2.0 * (rand() / double(RAND_MAX)) - 1.0);
            }

            // Compute Jacobian approximation
            auto Fx = F(x_nearby);
            if (Fx.empty()) continue;

            Matrix J(m, std::vector<double>(n));
            double eps = 1e-6;

            for (size_t j = 0; j < n; ++j) {
                Point x_perturb = x_nearby;
                x_perturb[j] += eps;
                auto Fx_perturb = F(x_perturb);
                if (Fx_perturb.empty()) continue;

                for (size_t i = 0; i < m; ++i) {
                    J[i][j] = (Fx_perturb[0][i] - Fx[0][i]) / eps;
                }
            }

            // Compute x* = J^T y*
            Subgradient x_star(n, 0.0);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    x_star[i] += J[j][i] * y_star[j];
                }
            }

            coderivative.push_back(x_star);
        }

        return coderivative;
    }

    /**
     * @brief Compute normal cone to intersection N(C₁ ∩ C₂)(x) (6.3.1)
     *
     * Under qualification: N(C₁ ∩ C₂)(x) = N_C₁(x) + N_C₂(x)
     *
     * @param inSet1 Characteristic function of C₁
     * @param inSet2 Characteristic function of C₂
     * @param x Point in intersection
     * @param h Sampling radius
     * @return Normal cone to intersection
     */
    static SubdifferentialSet normalToIntersection(
        std::function<bool(const Point&)> inSet1,
        std::function<bool(const Point&)> inSet2,
        const Point& x,
        double h = 1e-5) {

        auto N1 = NormalTangentCones::computeFrechetNormal(inSet1, x, h, 20);
        auto N2 = NormalTangentCones::computeFrechetNormal(inSet2, x, h, 20);

        SubdifferentialSet N_intersection;

        // Compute Minkowski sum N1 + N2
        for (const auto& n1 : N1) {
            for (const auto& n2 : N2) {
                Subgradient sum(n1.size());
                for (size_t i = 0; i < n1.size(); ++i) {
                    sum[i] = n1[i] + n2[i];
                }
                N_intersection.push_back(sum);
            }
        }

        return N_intersection;
    }
};

/**
 * @class GradedSubdifferential
 * @brief Graded subdifferentials (7.1)
 *
 * The graded subdifferential uses separable subspaces
 */
class GradedSubdifferential {
public:
    /**
     * @brief Compute ε-subdifferential
     *
     * ∂_ε f(x) = { v : f(y) ≥ f(x) + ⟨v, y-x⟩ - ε for all y }
     *
     * @param f Function
     * @param x Point
     * @param epsilon Tolerance parameter
     * @param h Sampling radius
     * @param n_samples Number of samples
     * @return ε-subdifferential
     */
    static SubdifferentialSet computeEpsilonSubdifferential(
        std::function<double(const Point&)> f,
        const Point& x,
        double epsilon,
        double h = 1e-5,
        int n_samples = 30) {

        size_t n = x.size();
        SubdifferentialSet eps_subdiff;
        double fx = f(x);

        // Sample candidate subgradients
        for (int s = 0; s < n_samples; ++s) {
            Subgradient v(n);
            for (size_t i = 0; i < n; ++i) {
                v[i] = 5.0 * (2.0 * (rand() / double(RAND_MAX)) - 1.0);
            }

            // Check ε-subgradient condition
            bool is_eps_subgrad = true;
            for (int t = 0; t < 20; ++t) {
                Point y = x;
                for (size_t i = 0; i < n; ++i) {
                    y[i] += 2.0 * h * (2.0 * (rand() / double(RAND_MAX)) - 1.0);
                }

                double fy = f(y);
                double inner_prod = 0.0;
                for (size_t i = 0; i < n; ++i) {
                    inner_prod += v[i] * (y[i] - x[i]);
                }

                if (fy < fx + inner_prod - epsilon - 1e-6) {
                    is_eps_subgrad = false;
                    break;
                }
            }

            if (is_eps_subgrad) {
                eps_subdiff.push_back(v);
            }
        }

        return eps_subdiff;
    }

    /**
     * @brief Compute graded normal cone
     *
     * Uses hierarchical decomposition of space
     *
     * @param inSet Characteristic function
     * @param x Point
     * @param grades Vector of grade parameters
     * @param h Sampling radius
     * @return Graded normal cone
     */
    static std::vector<SubdifferentialSet> computeGradedNormal(
        std::function<bool(const Point&)> inSet,
        const Point& x,
        const std::vector<double>& grades,
        double h = 1e-5) {

        std::vector<SubdifferentialSet> graded_normals;

        for (double grade : grades) {
            // Compute normals at grade level
            auto normals = NormalTangentCones::computeFrechetNormal(
                inSet, x, h * grade, 20);
            graded_normals.push_back(normals);
        }

        return graded_normals;
    }

    /**
     * @brief Compute Ioffe subdifferential
     *
     * Related to graded subdifferential with specific grade structure
     *
     * @param f Function
     * @param x Point
     * @param h Sampling radius
     * @return Ioffe subdifferential
     */
    static SubdifferentialSet computeIoffeSubdifferential(
        std::function<double(const Point&)> f,
        const Point& x,
        double h = 1e-5) {

        // Ioffe subdifferential uses special limit construction
        // Simplified implementation via limiting process
        return LimitingSubdifferential::computeLimiting(f, x, h, 30);
    }
};

/**
 * @class SubdifferentialCalculus
 * @brief Calculus rules for advanced subdifferentials (5.3.3, 6.4)
 */
class SubdifferentialCalculus {
public:
    /**
     * @brief Sum rule for subdifferentials
     *
     * Under qualification: ∂(f + g)(x) ⊆ ∂f(x) + ∂g(x)
     *
     * @param subdiff_f Subdifferential of f at x
     * @param subdiff_g Subdifferential of g at x
     * @return Subdifferential of f + g
     */
    static SubdifferentialSet sumRule(
        const SubdifferentialSet& subdiff_f,
        const SubdifferentialSet& subdiff_g) {

        SubdifferentialSet sum;

        for (const auto& sf : subdiff_f) {
            for (const auto& sg : subdiff_g) {
                Subgradient s(sf.size());
                for (size_t i = 0; i < sf.size(); ++i) {
                    s[i] = sf[i] + sg[i];
                }
                sum.push_back(s);
            }
        }

        return sum;
    }

    /**
     * @brief Chain rule for composition
     *
     * ∂(g ∘ f)(x) ⊆ ∂g(f(x)) ∘ Df(x)
     *
     * @param subdiff_g Subdifferential of g at f(x)
     * @param Df_x Jacobian of f at x
     * @return Subdifferential of composition
     */
    static SubdifferentialSet chainRule(
        const SubdifferentialSet& subdiff_g,
        const Matrix& Df_x) {

        SubdifferentialSet composition;

        for (const auto& sg : subdiff_g) {
            // Compute Df_x^T * sg
            Subgradient result(Df_x[0].size(), 0.0);
            for (size_t i = 0; i < Df_x[0].size(); ++i) {
                for (size_t j = 0; j < Df_x.size(); ++j) {
                    result[i] += Df_x[j][i] * sg[j];
                }
            }
            composition.push_back(result);
        }

        return composition;
    }

    /**
     * @brief Max rule for subdifferential of maximum
     *
     * ∂ max{f₁,...,f_m}(x) ⊆ conv{ ∂f_i(x) : i ∈ I(x) }
     * where I(x) = { i : f_i(x) = max_j f_j(x) }
     *
     * @param subdifferentials Vector of subdifferentials
     * @param function_values Values of f_i at x
     * @param tolerance Tolerance for active set
     * @return Subdifferential of maximum
     */
    static SubdifferentialSet maxRule(
        const std::vector<SubdifferentialSet>& subdifferentials,
        const std::vector<double>& function_values,
        double tolerance = 1e-10) {

        // Find active indices
        double max_val = *std::max_element(function_values.begin(), function_values.end());
        std::vector<size_t> active_indices;

        for (size_t i = 0; i < function_values.size(); ++i) {
            if (std::abs(function_values[i] - max_val) < tolerance) {
                active_indices.push_back(i);
            }
        }

        // Collect subdifferentials from active functions
        SubdifferentialSet max_subdiff;
        for (size_t idx : active_indices) {
            for (const auto& subgrad : subdifferentials[idx]) {
                max_subdiff.push_back(subgrad);
            }
        }

        return max_subdiff;
    }

    /**
     * @brief Compute convex hull of subdifferential set
     *
     * @param subdiff Set of subgradients
     * @param n_convex_combinations Number of convex combinations
     * @return Convex hull (sampled)
     */
    static SubdifferentialSet convexHull(
        const SubdifferentialSet& subdiff,
        int n_convex_combinations = 50) {

        if (subdiff.empty()) return subdiff;

        SubdifferentialSet hull;
        size_t n = subdiff[0].size();

        for (int c = 0; c < n_convex_combinations; ++c) {
            // Generate random convex combination
            std::vector<double> weights(subdiff.size());
            double sum_weights = 0.0;

            for (size_t i = 0; i < subdiff.size(); ++i) {
                weights[i] = rand() / double(RAND_MAX);
                sum_weights += weights[i];
            }

            // Normalize weights
            for (double& w : weights) {
                w /= sum_weights;
            }

            // Compute convex combination
            Subgradient combination(n, 0.0);
            for (size_t i = 0; i < subdiff.size(); ++i) {
                for (size_t j = 0; j < n; ++j) {
                    combination[j] += weights[i] * subdiff[i][j];
                }
            }

            hull.push_back(combination);
        }

        return hull;
    }
};

} // namespace maths::analysis

#endif // MATHS_ANALYSIS_ADVANCED_SUBDIFFERENTIALS_HPP
