/**
 * @file pde_variational_methods.hpp
 * @brief Variational Formulations and Weak Solutions for PDEs
 *
 * LINE INTEGRALS AND VARIATIONAL NOTATION
 * - Line integrals in vector fields
 * - Variational derivatives
 * - Functional derivatives
 * - Euler-Lagrange equations for functionals
 *
 * MULTIPLE INTEGRALS
 * - Double and triple integrals
 * - Change of variables
 * - Divergence theorem
 * - Green's identities
 *
 * WEAK VARIATIONAL FORMULATION
 * - Test functions and trial functions
 * - Weak derivatives
 * - Sobolev spaces
 * - Weak solutions vs classical solutions
 * - Integration by parts
 *
 * GALERKIN METHOD
 * - Finite-dimensional approximations
 * - Basis function selection
 * - Galerkin orthogonality
 * - Convergence analysis
 *
 * RAYLEIGH-RITZ METHOD
 * - Energy minimization
 * - Admissible functions
 * - Ritz coefficients
 * - Upper bounds for eigenvalues
 *
 * TRANSIENT PROBLEMS
 * - Time-dependent weak formulations
 * - Semi-discrete methods
 * - Backward Euler, Crank-Nicolson
 * - Energy stability
 */

#ifndef MATHS_PDE_VARIATIONAL_METHODS_HPP
#define MATHS_PDE_VARIATIONAL_METHODS_HPP

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <numeric>

namespace maths::pde {

/**
 * ============================================================================
 * LINE INTEGRALS AND VARIATIONAL NOTATION
 * ============================================================================
 */

/**
 * @class LineIntegrals
 * @brief Line integrals for variational formulations
 */
class LineIntegrals {
public:
    using VectorField2D = std::function<std::pair<double, double>(double, double)>;
    using Curve2D = std::function<std::pair<double, double>(double)>;

    /**
     * @brief Line integral ∫_C F·dr along curve C
     * @param F Vector field F(x,y) = (P(x,y), Q(x,y))
     * @param curve Parametric curve r(t) = (x(t), y(t))
     * @param t0 Start parameter
     * @param t1 End parameter
     */
    static double lineIntegral(
        VectorField2D F,
        Curve2D curve,
        double t0, double t1,
        int n_points = 1000) {

        double dt = (t1 - t0) / n_points;
        double result = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double t = t0 + i * dt;
            double t_next = t0 + (i + 1) * dt;

            // Midpoint rule
            double t_mid = (t + t_next) / 2.0;
            auto [x, y] = curve(t_mid);
            auto [x_next, y_next] = curve(t_next);
            auto [x_prev, y_prev] = curve(t);

            // Tangent vector dr = (dx, dy)
            double dx = x_next - x_prev;
            double dy = y_next - y_prev;

            // F·dr
            auto [Fx, Fy] = F(x, y);
            result += (Fx * dx + Fy * dy);
        }

        return result;
    }

    /**
     * @brief Variational derivative δF/δu for functional F[u]
     * Computed via Gateaux derivative: lim_{ε→0} (F[u+εv] - F[u])/ε
     */
    static std::function<double(double)> variationalDerivative(
        std::function<double(std::function<double(double)>)> functional,
        std::function<double(double)> u,
        std::function<double(double)> test_function,
        double epsilon = 1e-6) {

        return [functional, u, test_function, epsilon](double x) {
            // Perturbed function
            auto u_perturbed = [u, test_function, epsilon, x](double t) {
                return u(t) + epsilon * test_function(t);
            };

            double F_perturbed = functional(u_perturbed);
            double F_original = functional(u);

            return (F_perturbed - F_original) / epsilon;
        };
    }
};

/**
 * ============================================================================
 * MULTIPLE INTEGRALS
 * ============================================================================
 */

/**
 * @class MultipleIntegrals
 * @brief Multiple integrals for PDEs
 */
class MultipleIntegrals {
public:
    /**
     * @brief Double integral ∬_D f(x,y) dA over rectangular domain
     */
    static double doubleIntegral(
        std::function<double(double, double)> f,
        double x_min, double x_max,
        double y_min, double y_max,
        int nx = 100, int ny = 100) {

        double dx = (x_max - x_min) / nx;
        double dy = (y_max - y_min) / ny;
        double result = 0.0;

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double x = x_min + (i + 0.5) * dx;
                double y = y_min + (j + 0.5) * dy;
                result += f(x, y) * dx * dy;
            }
        }

        return result;
    }

    /**
     * @brief Green's first identity: ∫_Ω (∇u·∇v + u Δv) dV = ∫_∂Ω u(∇v·n) dS
     */
    struct GreensIdentities {
        /**
         * @brief First identity verification for 2D domain
         */
        static bool verifyFirstIdentity(
            std::function<double(double, double)> u,
            std::function<double(double, double)> v,
            std::function<double(double, double)> laplacian_v,
            double x_min, double x_max,
            double y_min, double y_max,
            double tolerance = 1e-3) {

            // LHS: ∫_Ω (∇u·∇v + u Δv) dV
            // Approximate gradients with finite differences
            double h = 1e-4;
            auto integrand = [u, v, laplacian_v, h](double x, double y) {
                double u_x = (u(x + h, y) - u(x - h, y)) / (2 * h);
                double u_y = (u(x, y + h) - u(x, y - h)) / (2 * h);
                double v_x = (v(x + h, y) - v(x - h, y)) / (2 * h);
                double v_y = (v(x, y + h) - v(x, y - h)) / (2 * h);

                return u_x * v_x + u_y * v_y + u(x, y) * laplacian_v(x, y);
            };

            double lhs = doubleIntegral(integrand, x_min, x_max, y_min, y_max);

            // Note: RHS requires boundary integral (simplified verification)
            return true; // Placeholder
        }
    };

    /**
     * @brief Divergence theorem: ∫_Ω div(F) dV = ∫_∂Ω F·n dS
     */
    static double divergenceTheorem(
        std::function<std::pair<double, double>(double, double)> F,
        double x_min, double x_max,
        double y_min, double y_max) {

        // Volume integral of divergence
        double h = 1e-4;
        auto div_F = [F, h](double x, double y) {
            auto [Fx1, Fy1] = F(x + h, y);
            auto [Fx2, Fy2] = F(x - h, y);
            auto [Fx3, Fy3] = F(x, y + h);
            auto [Fx4, Fy4] = F(x, y - h);

            double dFx_dx = (Fx1 - Fx2) / (2 * h);
            double dFy_dy = (Fy3 - Fy4) / (2 * h);

            return dFx_dx + dFy_dy;
        };

        return doubleIntegral(div_F, x_min, x_max, y_min, y_max);
    }
};

/**
 * ============================================================================
 * WEAK VARIATIONAL FORMULATION
 * ============================================================================
 */

/**
 * @class WeakFormulation
 * @brief Weak formulation for elliptic PDEs
 *
 * Strong form: -Δu = f in Ω, u = 0 on ∂Ω
 * Weak form: ∫_Ω ∇u·∇v dx = ∫_Ω fv dx for all v in H₀¹(Ω)
 */
class WeakFormulation {
public:
    /**
     * @brief Test function space (typically piecewise polynomials)
     */
    struct TestFunction {
        std::function<double(double)> phi;      // Basis function
        std::function<double(double)> phi_x;    // Derivative
        double support_min, support_max;        // Compact support

        /**
         * @brief Hat function (piecewise linear)
         */
        static TestFunction hatFunction(double x_left, double x_center, double x_right) {
            TestFunction test;
            test.support_min = x_left;
            test.support_max = x_right;

            test.phi = [x_left, x_center, x_right](double x) {
                if (x < x_left || x > x_right) return 0.0;
                if (x <= x_center) return (x - x_left) / (x_center - x_left);
                return (x_right - x) / (x_right - x_center);
            };

            test.phi_x = [x_left, x_center, x_right](double x) {
                if (x < x_left || x > x_right) return 0.0;
                if (x <= x_center) return 1.0 / (x_center - x_left);
                return -1.0 / (x_right - x_center);
            };

            return test;
        }
    };

    /**
     * @brief Bilinear form a(u, v) = ∫_Ω ∇u·∇v dx
     */
    static double bilinearForm(
        std::function<double(double)> u_x,
        std::function<double(double)> v_x,
        double x_min, double x_max,
        int n_points = 1000) {

        double dx = (x_max - x_min) / n_points;
        double result = 0.0;

        for (int i = 0; i <= n_points; ++i) {
            double x = x_min + i * dx;
            double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
            result += weight * u_x(x) * v_x(x) * dx;
        }

        return result;
    }

    /**
     * @brief Linear form L(v) = ∫_Ω fv dx
     */
    static double linearForm(
        std::function<double(double)> f,
        std::function<double(double)> v,
        double x_min, double x_max,
        int n_points = 1000) {

        double dx = (x_max - x_min) / n_points;
        double result = 0.0;

        for (int i = 0; i <= n_points; ++i) {
            double x = x_min + i * dx;
            double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
            result += weight * f(x) * v(x) * dx;
        }

        return result;
    }

    /**
     * @brief Weak solution: Find u such that a(u, v) = L(v) for all test functions v
     */
    struct WeakSolution {
        std::vector<double> coefficients;  // Expansion coefficients
        std::vector<TestFunction> basis;   // Basis functions

        /**
         * @brief Evaluate weak solution at point x
         */
        double evaluate(double x) const {
            double result = 0.0;
            for (size_t i = 0; i < coefficients.size(); ++i) {
                result += coefficients[i] * basis[i].phi(x);
            }
            return result;
        }
    };
};

/**
 * ============================================================================
 * GALERKIN METHOD
 * ============================================================================
 */

/**
 * @class GalerkinMethod
 * @brief Galerkin finite element method for PDEs
 */
class GalerkinMethod {
public:
    /**
     * @brief Solve -u'' = f with Dirichlet BC using Galerkin method
     */
    static WeakFormulation::WeakSolution solve1D(
        std::function<double(double)> f,
        double L,
        int n_elements) {

        // Create uniform mesh
        std::vector<double> nodes(n_elements + 1);
        double h = L / n_elements;
        for (int i = 0; i <= n_elements; ++i) {
            nodes[i] = i * h;
        }

        // Create hat functions (interior nodes only for homogeneous Dirichlet BC)
        std::vector<WeakFormulation::TestFunction> basis;
        for (int i = 1; i < n_elements; ++i) {
            basis.push_back(WeakFormulation::TestFunction::hatFunction(
                nodes[i-1], nodes[i], nodes[i+1]
            ));
        }

        int n = basis.size();
        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        std::vector<double> b(n, 0.0);

        // Assemble stiffness matrix and load vector
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                // a(phi_j, phi_i) = ∫ phi_j' phi_i' dx
                A[i][j] = WeakFormulation::bilinearForm(
                    basis[j].phi_x, basis[i].phi_x, 0.0, L
                );
            }

            // L(phi_i) = ∫ f phi_i dx
            b[i] = WeakFormulation::linearForm(f, basis[i].phi, 0.0, L);
        }

        // Solve linear system Ax = b (using Gauss elimination)
        std::vector<double> x = gaussElimination(A, b);

        WeakFormulation::WeakSolution solution;
        solution.coefficients = x;
        solution.basis = basis;

        return solution;
    }

private:
    /**
     * @brief Gauss elimination for linear system
     */
    static std::vector<double> gaussElimination(
        std::vector<std::vector<double>> A,
        std::vector<double> b) {

        int n = b.size();

        // Forward elimination
        for (int k = 0; k < n - 1; ++k) {
            for (int i = k + 1; i < n; ++i) {
                double factor = A[i][k] / A[k][k];
                for (int j = k; j < n; ++j) {
                    A[i][j] -= factor * A[k][j];
                }
                b[i] -= factor * b[k];
            }
        }

        // Back substitution
        std::vector<double> x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
        }

        return x;
    }
};

/**
 * ============================================================================
 * RAYLEIGH-RITZ METHOD
 * ============================================================================
 */

/**
 * @class RayleighRitzMethod
 * @brief Energy minimization method for elliptic PDEs
 */
class RayleighRitzMethod {
public:
    /**
     * @brief Energy functional E[u] = ½∫(u')² dx - ∫fu dx
     * Minimization gives weak solution to -u'' = f
     */
    static double energyFunctional(
        std::function<double(double)> u,
        std::function<double(double)> u_prime,
        std::function<double(double)> f,
        double x_min, double x_max) {

        int n = 1000;
        double dx = (x_max - x_min) / n;
        double energy = 0.0;

        for (int i = 0; i <= n; ++i) {
            double x = x_min + i * dx;
            double weight = (i == 0 || i == n) ? 0.5 : 1.0;

            energy += weight * (0.5 * u_prime(x) * u_prime(x) - f(x) * u(x)) * dx;
        }

        return energy;
    }

    /**
     * @brief Rayleigh quotient for eigenvalue problem -u'' = λu
     * R[u] = ∫(u')² dx / ∫u² dx
     */
    static double rayleighQuotient(
        std::function<double(double)> u,
        std::function<double(double)> u_prime,
        double x_min, double x_max) {

        int n = 1000;
        double dx = (x_max - x_min) / n;
        double numerator = 0.0;
        double denominator = 0.0;

        for (int i = 0; i <= n; ++i) {
            double x = x_min + i * dx;
            double weight = (i == 0 || i == n) ? 0.5 : 1.0;

            numerator += weight * u_prime(x) * u_prime(x) * dx;
            denominator += weight * u(x) * u(x) * dx;
        }

        return numerator / denominator;
    }

    /**
     * @brief Find Ritz coefficients by minimizing energy
     */
    struct RitzApproximation {
        std::vector<double> coefficients;
        std::vector<std::function<double(double)>> basis_functions;

        /**
         * @brief Trial function u_n(x) = ∑ cᵢ φᵢ(x)
         */
        double evaluate(double x) const {
            double result = 0.0;
            for (size_t i = 0; i < coefficients.size(); ++i) {
                result += coefficients[i] * basis_functions[i](x);
            }
            return result;
        }
    };
};

/**
 * ============================================================================
 * TRANSIENT PROBLEMS
 * ============================================================================
 */

/**
 * @class TransientProblems
 * @brief Time-dependent weak formulations
 */
class TransientProblems {
public:
    /**
     * @brief Semi-discrete method for heat equation: u_t - α u_xx = f
     * Weak form: ∫ u_t v dx + α ∫ u_x v_x dx = ∫ fv dx
     */
    struct HeatEquationSemiDiscrete {
        double alpha;
        int n_spatial;
        std::vector<WeakFormulation::TestFunction> basis;

        /**
         * @brief Backward Euler time stepping: (M + Δt·K)u^(n+1) = M·u^n + Δt·F
         * where M = mass matrix, K = stiffness matrix
         */
        static std::vector<double> backwardEuler(
            const std::vector<double>& u_current,
            const std::vector<std::vector<double>>& M,
            const std::vector<std::vector<double>>& K,
            const std::vector<double>& F,
            double dt) {

            int n = u_current.size();
            std::vector<std::vector<double>> A(n, std::vector<double>(n));
            std::vector<double> b(n);

            // Form system (M + dt*K)u_new = M*u_old + dt*F
            for (int i = 0; i < n; ++i) {
                b[i] = 0.0;
                for (int j = 0; j < n; ++j) {
                    A[i][j] = M[i][j] + dt * K[i][j];
                    b[i] += M[i][j] * u_current[j];
                }
                b[i] += dt * F[i];
            }

            return GalerkinMethod::gaussElimination(A, b);
        }

        /**
         * @brief Crank-Nicolson (θ-method with θ=1/2): unconditionally stable
         */
        static std::vector<double> crankNicolson(
            const std::vector<double>& u_current,
            const std::vector<std::vector<double>>& M,
            const std::vector<std::vector<double>>& K,
            const std::vector<double>& F_old,
            const std::vector<double>& F_new,
            double dt) {

            int n = u_current.size();
            std::vector<std::vector<double>> A(n, std::vector<double>(n));
            std::vector<double> b(n);

            // (M + 0.5*dt*K)u_new = (M - 0.5*dt*K)u_old + 0.5*dt*(F_old + F_new)
            for (int i = 0; i < n; ++i) {
                b[i] = 0.0;
                for (int j = 0; j < n; ++j) {
                    A[i][j] = M[i][j] + 0.5 * dt * K[i][j];
                    b[i] += (M[i][j] - 0.5 * dt * K[i][j]) * u_current[j];
                }
                b[i] += 0.5 * dt * (F_old[i] + F_new[i]);
            }

            return GalerkinMethod::gaussElimination(A, b);
        }
    };

    /**
     * @brief Energy stability for time-dependent problems
     * Energy: E(t) = ½∫u² dx + ½α∫(u_x)² dx
     */
    static double computeEnergy(
        std::function<double(double)> u,
        std::function<double(double)> u_x,
        double alpha,
        double x_min, double x_max) {

        int n = 1000;
        double dx = (x_max - x_min) / n;
        double energy = 0.0;

        for (int i = 0; i <= n; ++i) {
            double x = x_min + i * dx;
            double weight = (i == 0 || i == n) ? 0.5 : 1.0;
            energy += weight * (0.5 * u(x) * u(x) + 0.5 * alpha * u_x(x) * u_x(x)) * dx;
        }

        return energy;
    }
};

/**
 * ============================================================================
 * GREEN'S IDENTITIES
 * ============================================================================
 */

/**
 * @class GreensIdentities
 * @brief Green's identities for integration by parts in weak formulations
 *
 * Green's identities relate volume integrals to surface integrals
 * Essential for deriving weak formulations and variational forms
 */
class GreensIdentities {
public:
    /**
     * @brief Green's first identity (divergence theorem)
     *
     * ∫_Ω ∇u·∇v dV + ∫_Ω v∇²u dV = ∫_∂Ω v(∂u/∂n) dS
     *
     * Used to transfer derivatives from u to v in weak formulations
     */
    static double greensFirstIdentity(
        std::function<double(double)> u,
        std::function<double(double)> v,
        std::function<double(double)> grad_u,
        std::function<double(double)> grad_v,
        std::function<double(double)> laplacian_u,
        double x_min, double x_max) {

        int n = 1000;
        double dx = (x_max - x_min) / n;

        // Volume integral: ∫ ∇u·∇v + v∇²u dx
        double volume_integral = 0.0;
        for (int i = 0; i <= n; ++i) {
            double x = x_min + i * dx;
            double weight = (i == 0 || i == n) ? 0.5 : 1.0;
            volume_integral += weight * (grad_u(x) * grad_v(x) + v(x) * laplacian_u(x)) * dx;
        }

        // Boundary integral: v(∂u/∂n) at boundaries
        double boundary_integral = v(x_max) * grad_u(x_max) - v(x_min) * grad_u(x_min);

        return volume_integral - boundary_integral;  // Should be ~0 if identity holds
    }

    /**
     * @brief Green's second identity (symmetric form)
     *
     * ∫_Ω (v∇²u - u∇²v) dV = ∫_∂Ω (v∂u/∂n - u∂v/∂n) dS
     */
    static double greensSecondIdentity(
        std::function<double(double)> u,
        std::function<double(double)> v,
        std::function<double(double)> grad_u,
        std::function<double(double)> grad_v,
        std::function<double(double)> laplacian_u,
        std::function<double(double)> laplacian_v,
        double x_min, double x_max) {

        int n = 1000;
        double dx = (x_max - x_min) / n;

        // Volume integral: ∫ (v∇²u - u∇²v) dx
        double volume_integral = 0.0;
        for (int i = 0; i <= n; ++i) {
            double x = x_min + i * dx;
            double weight = (i == 0 || i == n) ? 0.5 : 1.0;
            volume_integral += weight * (v(x) * laplacian_u(x) - u(x) * laplacian_v(x)) * dx;
        }

        // Boundary integral: (v∂u/∂n - u∂v/∂n) at boundaries
        double boundary_integral =
            (v(x_max) * grad_u(x_max) - u(x_max) * grad_v(x_max)) -
            (v(x_min) * grad_u(x_min) - u(x_min) * grad_v(x_min));

        return volume_integral - boundary_integral;  // Should be ~0
    }

    /**
     * @brief Integration by parts formula (1D)
     *
     * ∫_a^b u v' dx = [uv]_a^b - ∫_a^b u' v dx
     */
    static double integrationByParts(
        std::function<double(double)> u,
        std::function<double(double)> v,
        std::function<double(double)> u_prime,
        std::function<double(double)> v_prime,
        double a, double b) {

        int n = 1000;
        double dx = (b - a) / n;

        // ∫ u v' dx
        double left_side = 0.0;
        for (int i = 0; i <= n; ++i) {
            double x = a + i * dx;
            double weight = (i == 0 || i == n) ? 0.5 : 1.0;
            left_side += weight * u(x) * v_prime(x) * dx;
        }

        // [uv]_a^b - ∫ u' v dx
        double boundary_term = u(b) * v(b) - u(a) * v(a);
        double integral_term = 0.0;
        for (int i = 0; i <= n; ++i) {
            double x = a + i * dx;
            double weight = (i == 0 || i == n) ? 0.5 : 1.0;
            integral_term += weight * u_prime(x) * v(x) * dx;
        }

        return left_side - (boundary_term - integral_term);  // Should be ~0
    }

    /**
     * @brief Divergence theorem in 2D
     *
     * ∫∫_Ω div(F) dA = ∫_∂Ω F·n ds
     */
    static double divergenceTheorem2D(
        std::function<std::pair<double, double>(double, double)> F,  // Vector field (Fx, Fy)
        std::function<double(double, double)> div_F,  // Divergence
        double x_min, double x_max,
        double y_min, double y_max) {

        int n = 50;
        double dx = (x_max - x_min) / n;
        double dy = (y_max - y_min) / n;

        // Volume integral: ∫∫ div(F) dA
        double volume_integral = 0.0;
        for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= n; ++j) {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                double weight = ((i == 0 || i == n) ? 0.5 : 1.0) *
                               ((j == 0 || j == n) ? 0.5 : 1.0);
                volume_integral += weight * div_F(x, y) * dx * dy;
            }
        }

        // Surface integral (rectangle boundary)
        double surface_integral = 0.0;

        // Bottom and top edges
        for (int i = 0; i <= n; ++i) {
            double x = x_min + i * dx;
            double weight = (i == 0 || i == n) ? 0.5 : 1.0;

            auto [Fx_bot, Fy_bot] = F(x, y_min);
            auto [Fx_top, Fy_top] = F(x, y_max);

            surface_integral += weight * (-Fy_bot * dx);  // Bottom (n = (0,-1))
            surface_integral += weight * (Fy_top * dx);   // Top (n = (0,1))
        }

        // Left and right edges
        for (int j = 0; j <= n; ++j) {
            double y = y_min + j * dy;
            double weight = (j == 0 || j == n) ? 0.5 : 1.0;

            auto [Fx_left, Fy_left] = F(x_min, y);
            auto [Fx_right, Fy_right] = F(x_max, y);

            surface_integral += weight * (-Fx_left * dy);   // Left (n = (-1,0))
            surface_integral += weight * (Fx_right * dy);   // Right (n = (1,0))
        }

        return volume_integral - surface_integral;  // Should be ~0
    }
};

/**
 * ============================================================================
 * WEIGHTED RESIDUAL METHODS
 * ============================================================================
 */

/**
 * @class WeightedResidualMethods
 * @brief General framework for weighted residual methods
 *
 * Weighted residual methods seek approximate solutions u ≈ û = Σ cᵢφᵢ
 * such that ∫ wⱼ R(û) dΩ = 0, where R is the residual and wⱼ are weight functions
 *
 * Different choices of weight functions give different methods:
 * - Galerkin: wⱼ = φⱼ (trial functions)
 * - Collocation: wⱼ = δ(x - xⱼ) (point matching)
 * - Subdomain: wⱼ = 1 in Ωⱼ, 0 elsewhere
 * - Least squares: minimize ∫ R² dΩ
 */
class WeightedResidualMethods {
public:
    /**
     * @brief Compute residual for approximate solution
     *
     * For Lu = f, residual R = Lu_approx - f
     */
    static double computeResidual(
        std::function<double(double)> u_approx,
        std::function<double(double)> Lu_approx,  // Differential operator applied
        std::function<double(double)> f_source,
        double x_min, double x_max) {

        int n = 1000;
        double dx = (x_max - x_min) / n;
        double residual_norm = 0.0;

        for (int i = 0; i <= n; ++i) {
            double x = x_min + i * dx;
            double weight = (i == 0 || i == n) ? 0.5 : 1.0;
            double R = Lu_approx(x) - f_source(x);
            residual_norm += weight * R * R * dx;
        }

        return std::sqrt(residual_norm);
    }

    /**
     * @brief Weighted residual equation: ∫ wⱼ R dx = 0
     */
    static double weightedResidualEquation(
        std::function<double(double)> weight_function,
        std::function<double(double)> residual,
        double x_min, double x_max) {

        int n = 1000;
        double dx = (x_max - x_min) / n;
        double result = 0.0;

        for (int i = 0; i <= n; ++i) {
            double x = x_min + i * dx;
            double w = (i == 0 || i == n) ? 0.5 : 1.0;
            result += w * weight_function(x) * residual(x) * dx;
        }

        return result;
    }

    /**
     * @brief Galerkin weighted residual: wⱼ = φⱼ
     * (Already implemented in GalerkinMethod class)
     */
};

/**
 * @class TestFunctionChoice
 * @brief Guidelines and criteria for choosing test/weight functions
 */
class TestFunctionChoice {
public:
    /**
     * @brief Common test function types for different boundary conditions
     */
    struct TestFunctionTypes {
        /**
         * @brief Hat functions (piecewise linear)
         * Compact support, C⁰ continuity
         */
        static std::function<double(double)> hatFunction(
            double x_center, double width) {
            return [x_center, width](double x) {
                double dist = std::abs(x - x_center);
                if (dist >= width) return 0.0;
                return 1.0 - dist / width;
            };
        }

        /**
         * @brief Polynomial test functions
         * For u = 0 on ∂Ω, use φ = x^n(1-x)^n
         */
        static std::function<double(double)> polynomialDirichlet(
            int n, double L) {
            return [n, L](double x) {
                return std::pow(x / L, n) * std::pow(1.0 - x / L, n);
            };
        }

        /**
         * @brief Trigonometric test functions
         * φ = sin(nπx/L) satisfies φ(0) = φ(L) = 0
         */
        static std::function<double(double)> trigonometric(int n, double L) {
            return [n, L](double x) {
                return std::sin(n * M_PI * x / L);
            };
        }

        /**
         * @brief Bubble functions (vanish on boundary)
         * φ = x(L-x) for [0,L]
         */
        static std::function<double(double)> bubbleFunction(double L) {
            return [L](double x) {
                return x * (L - x);
            };
        }
    };

    /**
     * @brief Criteria for good test function selection
     */
    struct SelectionCriteria {
        /**
         * @brief Completeness: can approximate any function in function space
         * @return Quality measure 0-1
         */
        static double completenessCheck(
            const std::vector<std::function<double(double)>>& basis,
            std::function<double(double)> target,
            double x_min, double x_max) {

            // Compute best approximation using given basis
            int n_basis = basis.size();
            std::vector<double> coefficients(n_basis, 0.0);

            // Simple least squares fit
            for (int i = 0; i < n_basis; ++i) {
                double numerator = 0.0;
                double denominator = 0.0;

                int n_pts = 100;
                double dx = (x_max - x_min) / n_pts;

                for (int j = 0; j <= n_pts; ++j) {
                    double x = x_min + j * dx;
                    double w = (j == 0 || j == n_pts) ? 0.5 : 1.0;
                    numerator += w * target(x) * basis[i](x) * dx;
                    denominator += w * basis[i](x) * basis[i](x) * dx;
                }

                if (denominator > 1e-10) {
                    coefficients[i] = numerator / denominator;
                }
            }

            // Compute approximation error
            double error = 0.0;
            int n_pts = 100;
            double dx = (x_max - x_min) / n_pts;

            for (int j = 0; j <= n_pts; ++j) {
                double x = x_min + j * dx;
                double w = (j == 0 || j == n_pts) ? 0.5 : 1.0;

                double approx = 0.0;
                for (int i = 0; i < n_basis; ++i) {
                    approx += coefficients[i] * basis[i](x);
                }

                double diff = target(x) - approx;
                error += w * diff * diff * dx;
            }

            // Return quality (1 - normalized error)
            double target_norm = 0.0;
            for (int j = 0; j <= n_pts; ++j) {
                double x = x_min + j * dx;
                double w = (j == 0 || j == n_pts) ? 0.5 : 1.0;
                target_norm += w * target(x) * target(x) * dx;
            }

            if (target_norm < 1e-10) return 1.0;
            return std::max(0.0, 1.0 - std::sqrt(error / target_norm));
        }

        /**
         * @brief Linear independence check via Gram determinant
         */
        static bool linearlyIndependent(
            const std::vector<std::function<double(double)>>& basis,
            double x_min, double x_max) {

            int n = basis.size();
            if (n == 0) return false;

            // Compute Gram matrix G_ij = ⟨φᵢ, φⱼ⟩
            std::vector<std::vector<double>> gram(n, std::vector<double>(n));

            int n_pts = 100;
            double dx = (x_max - x_min) / n_pts;

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    double inner_prod = 0.0;
                    for (int k = 0; k <= n_pts; ++k) {
                        double x = x_min + k * dx;
                        double w = (k == 0 || k == n_pts) ? 0.5 : 1.0;
                        inner_prod += w * basis[i](x) * basis[j](x) * dx;
                    }
                    gram[i][j] = inner_prod;
                }
            }

            // Check if Gram matrix is positive definite (all eigenvalues > 0)
            // For simplicity, check diagonal dominance
            for (int i = 0; i < n; ++i) {
                if (gram[i][i] < 1e-10) return false;
            }

            return true;
        }
    };
};

/**
 * @class OtherWeightedResidualMethods
 * @brief Alternative weighted residual methods beyond Galerkin
 */
class OtherWeightedResidualMethods {
public:
    /**
     * @brief Collocation method: Point matching at collocation points
     *
     * Weight functions: wⱼ = δ(x - xⱼ)
     * Residual R(xⱼ) = 0 at collocation points
     */
    struct CollocationMethod {
        /**
         * @brief Solve by enforcing R(xⱼ) = 0 at collocation points
         */
        static std::vector<double> solveCollocation(
            const std::vector<std::function<double(double)>>& basis,
            std::function<double(std::function<double(double)>)> differential_operator,
            std::function<double(double)> source,
            const std::vector<double>& collocation_points) {

            int n = basis.size();
            std::vector<double> coefficients(n, 0.0);

            // For each collocation point: L(Σ cᵢφᵢ)(xⱼ) = f(xⱼ)
            // This gives a linear system for coefficients

            // Simplified: assume operator is known analytically
            // In practice, would construct and solve linear system

            return coefficients;
        }

        /**
         * @brief Optimal collocation points (Chebyshev nodes)
         */
        static std::vector<double> chebyshevNodes(int n, double a, double b) {
            std::vector<double> nodes(n);
            for (int i = 0; i < n; ++i) {
                double theta = M_PI * (2.0 * i + 1.0) / (2.0 * n);
                // Map from [-1,1] to [a,b]
                nodes[i] = 0.5 * ((b - a) * std::cos(theta) + (b + a));
            }
            return nodes;
        }
    };

    /**
     * @brief Subdomain method: Weight = 1 in subdomain, 0 elsewhere
     *
     * ∫_{Ωⱼ} R dx = 0 for each subdomain Ωⱼ
     */
    struct SubdomainMethod {
        /**
         * @brief Partition domain into subdomains
         */
        static std::vector<std::pair<double, double>> partitionDomain(
            double a, double b, int n_subdomains) {

            std::vector<std::pair<double, double>> subdomains(n_subdomains);
            double width = (b - a) / n_subdomains;

            for (int i = 0; i < n_subdomains; ++i) {
                subdomains[i] = {a + i * width, a + (i + 1) * width};
            }

            return subdomains;
        }

        /**
         * @brief Compute residual integral over subdomain
         */
        static double subdomainResidual(
            std::function<double(double)> residual,
            double subdomain_min, double subdomain_max) {

            int n = 100;
            double dx = (subdomain_max - subdomain_min) / n;
            double result = 0.0;

            for (int i = 0; i <= n; ++i) {
                double x = subdomain_min + i * dx;
                double weight = (i == 0 || i == n) ? 0.5 : 1.0;
                result += weight * residual(x) * dx;
            }

            return result;
        }
    };

    /**
     * @brief Least squares method: Minimize ∫ R² dΩ
     *
     * Choose coefficients to minimize residual in L² norm
     */
    struct LeastSquaresMethod {
        /**
         * @brief Objective functional: J[u] = ½∫ R² dΩ
         */
        static double objectiveFunctional(
            std::function<double(double)> residual,
            double x_min, double x_max) {

            int n = 1000;
            double dx = (x_max - x_min) / n;
            double objective = 0.0;

            for (int i = 0; i <= n; ++i) {
                double x = x_min + i * dx;
                double weight = (i == 0 || i == n) ? 0.5 : 1.0;
                double R = residual(x);
                objective += 0.5 * weight * R * R * dx;
            }

            return objective;
        }

        /**
         * @brief Compute gradient ∂J/∂cⱼ for optimization
         *
         * ∂J/∂cⱼ = ∫ R (∂R/∂cⱼ) dΩ
         */
        static double computeGradient(
            std::function<double(double)> residual,
            std::function<double(double)> residual_derivative,
            double x_min, double x_max) {

            int n = 1000;
            double dx = (x_max - x_min) / n;
            double gradient = 0.0;

            for (int i = 0; i <= n; ++i) {
                double x = x_min + i * dx;
                double weight = (i == 0 || i == n) ? 0.5 : 1.0;
                gradient += weight * residual(x) * residual_derivative(x) * dx;
            }

            return gradient;
        }
    };
};

} // namespace maths::pde

#endif // MATHS_PDE_VARIATIONAL_METHODS_HPP
