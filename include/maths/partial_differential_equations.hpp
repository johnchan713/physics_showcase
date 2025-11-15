/**
 * @file partial_differential_equations.hpp
 * @brief Comprehensive Partial Differential Equations Theory
 *
 * PDE FUNDAMENTALS:
 * - Notation and Definitions (orders, linearity, homogeneity)
 * - Initial and Boundary Conditions (Cauchy, Dirichlet, Neumann, Robin)
 * - Classification of Second Order Equations (elliptic, parabolic, hyperbolic)
 * - Well-Known Equations (heat, wave, Laplace, Poisson, transport)
 * - Superposition Principle (linear equations)
 *
 * METHOD OF CHARACTERISTICS:
 * First Order Equations:
 * - Linear Equations (constant and variable coefficients)
 * - Quasi-linear Equations
 * - Fully Nonlinear Equations
 * - Geometrical Considerations (characteristic curves, Monge cone)
 * - Existence, uniqueness, and regularity theorems
 *
 * Second Order Equations:
 * - Linear and Quasi-linear Equations (characteristic surfaces)
 */

#ifndef MATHS_PDE_HPP
#define MATHS_PDE_HPP

#include <vector>
#include <cmath>
#include <functional>
#include <array>
#include <stdexcept>
#include <algorithm>
#include <map>
#include <string>

namespace maths::pde {

// Type aliases
using ScalarField2D = std::function<double(double, double)>;
using ScalarField3D = std::function<double(double, double, double)>;
using VectorField2D = std::function<std::array<double, 2>(double, double)>;
using VectorField3D = std::function<std::array<double, 3>(double, double, double)>;

/**
 * ============================================================================
 * PDE FUNDAMENTALS
 * ============================================================================
 */

/**
 * @class PDEClassification
 * @brief Classification and properties of PDEs
 */
class PDEClassification {
public:
    /**
     * @brief Order of PDE (highest derivative)
     */
    enum class Order {
        FIRST,
        SECOND,
        HIGHER
    };

    /**
     * @brief Linearity classification
     */
    enum class Linearity {
        LINEAR,           // Coefficients independent of u and derivatives
        QUASI_LINEAR,     // Linear in highest derivatives
        SEMI_LINEAR,      // Nonlinear in lower derivatives only
        FULLY_NONLINEAR   // Nonlinear in highest derivatives
    };

    /**
     * @brief Classification of second-order PDEs
     *
     * General form: A u_xx + 2B u_xy + C u_yy + lower order terms = 0
     * Discriminant: Δ = B² - AC
     */
    enum class SecondOrderType {
        ELLIPTIC,     // Δ < 0 (e.g., Laplace equation)
        PARABOLIC,    // Δ = 0 (e.g., heat equation)
        HYPERBOLIC    // Δ > 0 (e.g., wave equation)
    };

    /**
     * @brief Classify second-order PDE at a point (x, y)
     *
     * A u_xx + 2B u_xy + C u_yy + ... = 0
     */
    static SecondOrderType classify(double A, double B, double C) {
        double discriminant = B * B - A * C;

        if (std::abs(discriminant) < 1e-10) {
            return SecondOrderType::PARABOLIC;
        } else if (discriminant < 0) {
            return SecondOrderType::ELLIPTIC;
        } else {
            return SecondOrderType::HYPERBOLIC;
        }
    }

    /**
     * @brief Get canonical form based on classification
     */
    static std::string canonicalForm(SecondOrderType type) {
        switch (type) {
            case SecondOrderType::ELLIPTIC:
                return "u_ξξ + u_ηη = F(ξ, η, u, u_ξ, u_η)";
            case SecondOrderType::PARABOLIC:
                return "u_ηη = F(ξ, η, u, u_ξ, u_η)";
            case SecondOrderType::HYPERBOLIC:
                return "u_ξη = F(ξ, η, u, u_ξ, u_η)";
            default:
                return "Unknown";
        }
    }

    /**
     * @brief Check if PDE is homogeneous (no explicit function of independent variables)
     */
    static bool isHomogeneous(std::function<double(double, double, double, double, double)> F) {
        // Homogeneous if F(x, y, ku, ku_x, ku_y) = k F(x, y, u, u_x, u_y)
        // This is a simplified check
        double x = 1.0, y = 1.0, u = 1.0, ux = 1.0, uy = 1.0;
        double k = 2.0;

        double F1 = F(x, y, u, ux, uy);
        double F2 = F(x, y, k*u, k*ux, k*uy);

        return std::abs(F2 - k*F1) < 1e-6;
    }
};

/**
 * @class BoundaryConditions
 * @brief Initial and boundary conditions
 */
class BoundaryConditions {
public:
    /**
     * @brief Types of boundary conditions
     */
    enum class Type {
        DIRICHLET,    // u = g on boundary (fixed value)
        NEUMANN,      // ∂u/∂n = g on boundary (fixed normal derivative)
        ROBIN,        // αu + β∂u/∂n = g on boundary (mixed)
        CAUCHY        // Both u and ∂u/∂n specified (over-determined for elliptic)
    };

    /**
     * @brief Dirichlet boundary condition: u|_Γ = g
     */
    struct Dirichlet {
        ScalarField2D g;  // Boundary value

        double evaluate(double x, double y) const {
            return g(x, y);
        }

        bool satisfies(ScalarField2D u, double x, double y) const {
            return std::abs(u(x, y) - g(x, y)) < 1e-6;
        }
    };

    /**
     * @brief Neumann boundary condition: ∂u/∂n|_Γ = g
     */
    struct Neumann {
        ScalarField2D g;                    // Normal derivative value
        VectorField2D normal;               // Outward normal vector

        double evaluate(double x, double y) const {
            return g(x, y);
        }

        /**
         * @brief Check if solution satisfies Neumann condition
         */
        bool satisfies(ScalarField2D u, ScalarField2D u_x, ScalarField2D u_y,
                       double x, double y) const {
            auto n = normal(x, y);
            double du_dn = u_x(x, y) * n[0] + u_y(x, y) * n[1];
            return std::abs(du_dn - g(x, y)) < 1e-6;
        }
    };

    /**
     * @brief Robin (mixed) boundary condition: αu + β∂u/∂n = g
     */
    struct Robin {
        double alpha, beta;
        ScalarField2D g;
        VectorField2D normal;

        bool satisfies(ScalarField2D u, ScalarField2D u_x, ScalarField2D u_y,
                       double x, double y) const {
            auto n = normal(x, y);
            double du_dn = u_x(x, y) * n[0] + u_y(x, y) * n[1];
            double lhs = alpha * u(x, y) + beta * du_dn;
            return std::abs(lhs - g(x, y)) < 1e-6;
        }
    };

    /**
     * @brief Cauchy problem: specify u and ∂u/∂n on initial curve
     */
    struct Cauchy {
        ScalarField2D u0;      // Initial value
        ScalarField2D u0_n;    // Initial normal derivative

        std::pair<double, double> evaluate(double x, double y) const {
            return {u0(x, y), u0_n(x, y)};
        }
    };
};

/**
 * @class WellKnownPDEs
 * @brief Standard PDEs in mathematical physics
 */
class WellKnownPDEs {
public:
    /**
     * @brief Heat equation (parabolic): u_t = α u_xx
     *
     * Models diffusion, heat conduction
     * @param alpha Thermal diffusivity (α > 0)
     */
    struct HeatEquation {
        double alpha;

        HeatEquation(double a = 1.0) : alpha(a) {}

        /**
         * @brief Check if u(x,t) satisfies heat equation
         */
        bool verifySolution(std::function<double(double, double)> u,
                           double x, double t, double h = 1e-5) const {
            // Numerical derivatives
            double u_t = (u(x, t + h) - u(x, t - h)) / (2 * h);
            double u_xx = (u(x + h, t) - 2 * u(x, t) + u(x - h, t)) / (h * h);

            return std::abs(u_t - alpha * u_xx) < 1e-4;
        }

        /**
         * @brief Fundamental solution: G(x, t) = 1/√(4παt) exp(-x²/(4αt))
         */
        double fundamentalSolution(double x, double t) const {
            if (t <= 0) return 0.0;
            return 1.0 / std::sqrt(4 * M_PI * alpha * t) *
                   std::exp(-x * x / (4 * alpha * t));
        }
    };

    /**
     * @brief Wave equation (hyperbolic): u_tt = c² u_xx
     *
     * Models vibrations, wave propagation
     * @param c Wave speed
     */
    struct WaveEquation {
        double c;  // Wave speed

        WaveEquation(double wave_speed = 1.0) : c(wave_speed) {}

        /**
         * @brief Check if u(x,t) satisfies wave equation
         */
        bool verifySolution(std::function<double(double, double)> u,
                           double x, double t, double h = 1e-5) const {
            double u_tt = (u(x, t + h) - 2 * u(x, t) + u(x, t - h)) / (h * h);
            double u_xx = (u(x + h, t) - 2 * u(x, t) + u(x - h, t)) / (h * h);

            return std::abs(u_tt - c * c * u_xx) < 1e-4;
        }

        /**
         * @brief D'Alembert solution: u(x,t) = ½[f(x-ct) + f(x+ct)] + 1/(2c)∫_{x-ct}^{x+ct} g(s) ds
         */
        std::function<double(double, double)> dAlembertSolution(
            std::function<double(double)> f,    // Initial displacement
            std::function<double(double)> g     // Initial velocity
        ) const {
            return [this, f, g](double x, double t) {
                // First term: average of characteristic values
                double term1 = 0.5 * (f(x - c * t) + f(x + c * t));

                // Second term: integral of initial velocity
                double x_minus = x - c * t;
                double x_plus = x + c * t;
                int n_steps = 100;
                double dx = (x_plus - x_minus) / n_steps;
                double integral = 0.0;
                for (int i = 0; i < n_steps; ++i) {
                    double xi = x_minus + (i + 0.5) * dx;
                    integral += g(xi) * dx;
                }
                double term2 = integral / (2.0 * c);

                return term1 + term2;
            };
        }
    };

    /**
     * @brief Laplace equation (elliptic): Δu = u_xx + u_yy = 0
     *
     * Harmonic functions, steady-state problems
     */
    struct LaplaceEquation {
        /**
         * @brief Check if u(x,y) is harmonic (satisfies Laplace equation)
         */
        bool verifySolution(ScalarField2D u, double x, double y, double h = 1e-5) const {
            double u_xx = (u(x + h, y) - 2 * u(x, y) + u(x - h, y)) / (h * h);
            double u_yy = (u(x, y + h) - 2 * u(x, y) + u(x, y - h)) / (h * h);

            return std::abs(u_xx + u_yy) < 1e-4;
        }

        /**
         * @brief Mean value property: u(x0, y0) = 1/(2πr) ∫ u(x0 + r cos θ, y0 + r sin θ) dθ
         */
        bool verifyMeanValueProperty(ScalarField2D u, double x0, double y0, double r) const {
            int n_theta = 100;
            double dtheta = 2 * M_PI / n_theta;
            double average = 0.0;

            for (int i = 0; i < n_theta; ++i) {
                double theta = i * dtheta;
                double x = x0 + r * std::cos(theta);
                double y = y0 + r * std::sin(theta);
                average += u(x, y);
            }
            average /= n_theta;

            return std::abs(average - u(x0, y0)) < 1e-4;
        }
    };

    /**
     * @brief Poisson equation (elliptic): Δu = f
     */
    struct PoissonEquation {
        ScalarField2D f;  // Source term

        PoissonEquation(ScalarField2D source) : f(source) {}

        /**
         * @brief Check if u satisfies Poisson equation
         */
        bool verifySolution(ScalarField2D u, double x, double y, double h = 1e-5) const {
            double u_xx = (u(x + h, y) - 2 * u(x, y) + u(x - h, y)) / (h * h);
            double u_yy = (u(x, y + h) - 2 * u(x, y) + u(x, y - h)) / (h * h);

            return std::abs(u_xx + u_yy - f(x, y)) < 1e-4;
        }
    };

    /**
     * @brief Transport equation (hyperbolic): u_t + c·∇u = 0
     */
    struct TransportEquation {
        std::array<double, 2> velocity;  // Constant velocity (c_x, c_y)

        TransportEquation(double c_x, double c_y) : velocity{c_x, c_y} {}

        /**
         * @brief Solution via characteristics: u(x, y, t) = u0(x - c_x t, y - c_y t)
         */
        std::function<double(double, double, double)> solve(ScalarField2D u0) const {
            return [this, u0](double x, double y, double t) {
                return u0(x - velocity[0] * t, y - velocity[1] * t);
            };
        }
    };
};

/**
 * @class SuperpositionPrinciple
 * @brief Superposition for linear equations
 */
class SuperpositionPrinciple {
public:
    /**
     * @brief For linear homogeneous PDE: L[u] = 0
     *
     * If L[u1] = 0 and L[u2] = 0, then L[c1 u1 + c2 u2] = 0
     */
    static bool verifyLinearSuperposition(
        std::function<double(ScalarField2D)> L,
        ScalarField2D u1, ScalarField2D u2,
        double c1, double c2,
        double x, double y) {

        // Combined solution
        auto u_combined = [u1, u2, c1, c2](double x, double y) {
            return c1 * u1(x, y) + c2 * u2(x, y);
        };

        return std::abs(L(u_combined)) < 1e-6;
    }

    /**
     * @brief For inhomogeneous equation: L[u] = f
     *
     * General solution = homogeneous solution + particular solution
     */
    struct InhomogeneousSolution {
        ScalarField2D u_h;  // Homogeneous solution: L[u_h] = 0
        ScalarField2D u_p;  // Particular solution: L[u_p] = f

        ScalarField2D generalSolution() const {
            return [this](double x, double y) {
                return u_h(x, y) + u_p(x, y);
            };
        }
    };
};

/**
 * ============================================================================
 * METHOD OF CHARACTERISTICS
 * ============================================================================
 */

/**
 * @class FirstOrderLinearPDE
 * @brief First-order linear PDEs
 */
class FirstOrderLinearPDE {
public:
    /**
     * @brief Constant coefficient: a u_x + b u_y = 0
     *
     * General solution: u(x, y) = F(bx - ay) for arbitrary F
     */
    struct ConstantCoefficients {
        double a, b;

        ConstantCoefficients(double a_coeff, double b_coeff) : a(a_coeff), b(b_coeff) {}

        /**
         * @brief Characteristic curves: dy/dx = b/a
         */
        std::pair<double, double> characteristicDirection() const {
            return {a, b};
        }

        /**
         * @brief General solution u = F(ξ) where ξ = bx - ay
         */
        ScalarField2D generalSolution(std::function<double(double)> F) const {
            return [this, F](double x, double y) {
                double xi = b * x - a * y;
                return F(xi);
            };
        }

        /**
         * @brief Solve with initial condition u(x, 0) = g(x)
         */
        ScalarField2D solveIVP(std::function<double(double)> g) const {
            return [this, g](double x, double y) {
                // Along characteristic: bx - ay = bx0 - a·0 = bx0
                // So x0 = x - (a/b)y
                double x0 = x - (a / b) * y;
                return g(x0);
            };
        }
    };

    /**
     * @brief Variable coefficient: a(x,y) u_x + b(x,y) u_y = c(x,y)
     *
     * Method of characteristics:
     * dx/dt = a(x,y), dy/dt = b(x,y), du/dt = c(x,y)
     */
    struct VariableCoefficients {
        ScalarField2D a, b, c;

        VariableCoefficients(ScalarField2D a_coeff, ScalarField2D b_coeff, ScalarField2D c_coeff)
            : a(a_coeff), b(b_coeff), c(c_coeff) {}

        /**
         * @brief Characteristic equations (ODE system)
         *
         * dx/dt = a(x,y)
         * dy/dt = b(x,y)
         * du/dt = c(x,y)
         */
        struct CharacteristicSystem {
            std::vector<double> x_vals, y_vals, u_vals;
            std::vector<double> t_vals;
        };

        /**
         * @brief Solve characteristic curves using Euler method
         */
        CharacteristicSystem solveCharacteristics(
            double x0, double y0, double u0,
            double t_max, int n_steps) const {

            CharacteristicSystem solution;
            solution.x_vals.reserve(n_steps + 1);
            solution.y_vals.reserve(n_steps + 1);
            solution.u_vals.reserve(n_steps + 1);
            solution.t_vals.reserve(n_steps + 1);

            double dt = t_max / n_steps;
            double x = x0, y = y0, u = u0;
            double t = 0.0;

            for (int i = 0; i <= n_steps; ++i) {
                solution.x_vals.push_back(x);
                solution.y_vals.push_back(y);
                solution.u_vals.push_back(u);
                solution.t_vals.push_back(t);

                if (i < n_steps) {
                    // Euler step
                    x += dt * a(x, y);
                    y += dt * b(x, y);
                    u += dt * c(x, y);
                    t += dt;
                }
            }

            return solution;
        }
    };
};

/**
 * @class QuasilinearPDE
 * @brief First-order quasi-linear equations
 *
 * General form: a(x,y,u) u_x + b(x,y,u) u_y = c(x,y,u)
 */
class QuasilinearPDE {
public:
    using Coefficient = std::function<double(double, double, double)>;

    Coefficient a, b, c;

    QuasilinearPDE(Coefficient a_coeff, Coefficient b_coeff, Coefficient c_coeff)
        : a(a_coeff), b(b_coeff), c(c_coeff) {}

    /**
     * @brief Characteristic equations (Charpit's method)
     *
     * dx/dt = a(x,y,u)
     * dy/dt = b(x,y,u)
     * du/dt = c(x,y,u)
     */
    struct Characteristic {
        std::vector<double> x, y, u, t;
    };

    /**
     * @brief Solve characteristic with RK4
     */
    Characteristic solveCharacteristic(
        double x0, double y0, double u0,
        double t_max, int n_steps) const {

        Characteristic curve;
        curve.x.reserve(n_steps + 1);
        curve.y.reserve(n_steps + 1);
        curve.u.reserve(n_steps + 1);
        curve.t.reserve(n_steps + 1);

        double dt = t_max / n_steps;
        double x_curr = x0, y_curr = y0, u_curr = u0;

        for (int i = 0; i <= n_steps; ++i) {
            curve.x.push_back(x_curr);
            curve.y.push_back(y_curr);
            curve.u.push_back(u_curr);
            curve.t.push_back(i * dt);

            if (i < n_steps) {
                // RK4 step
                double k1_x = a(x_curr, y_curr, u_curr);
                double k1_y = b(x_curr, y_curr, u_curr);
                double k1_u = c(x_curr, y_curr, u_curr);

                double k2_x = a(x_curr + 0.5*dt*k1_x, y_curr + 0.5*dt*k1_y, u_curr + 0.5*dt*k1_u);
                double k2_y = b(x_curr + 0.5*dt*k1_x, y_curr + 0.5*dt*k1_y, u_curr + 0.5*dt*k1_u);
                double k2_u = c(x_curr + 0.5*dt*k1_x, y_curr + 0.5*dt*k1_y, u_curr + 0.5*dt*k1_u);

                double k3_x = a(x_curr + 0.5*dt*k2_x, y_curr + 0.5*dt*k2_y, u_curr + 0.5*dt*k2_u);
                double k3_y = b(x_curr + 0.5*dt*k2_x, y_curr + 0.5*dt*k2_y, u_curr + 0.5*dt*k2_u);
                double k3_u = c(x_curr + 0.5*dt*k2_x, y_curr + 0.5*dt*k2_y, u_curr + 0.5*dt*k2_u);

                double k4_x = a(x_curr + dt*k3_x, y_curr + dt*k3_y, u_curr + dt*k3_u);
                double k4_y = b(x_curr + dt*k3_x, y_curr + dt*k3_y, u_curr + dt*k3_u);
                double k4_u = c(x_curr + dt*k3_x, y_curr + dt*k3_y, u_curr + dt*k3_u);

                x_curr += (dt / 6.0) * (k1_x + 2*k2_x + 2*k3_x + k4_x);
                y_curr += (dt / 6.0) * (k1_y + 2*k2_y + 2*k3_y + k4_y);
                u_curr += (dt / 6.0) * (k1_u + 2*k2_u + 2*k3_u + k4_u);
            }
        }

        return curve;
    }

    /**
     * @brief Solve Cauchy problem with initial curve γ(s) = (x(s), y(s), u(s))
     */
    std::vector<Characteristic> solveFromInitialCurve(
        std::function<double(double)> x_initial,
        std::function<double(double)> y_initial,
        std::function<double(double)> u_initial,
        const std::vector<double>& s_values,
        double t_max, int n_steps_t) const {

        std::vector<Characteristic> characteristics;
        characteristics.reserve(s_values.size());

        for (double s : s_values) {
            double x0 = x_initial(s);
            double y0 = y_initial(s);
            double u0 = u_initial(s);

            characteristics.push_back(solveCharacteristic(x0, y0, u0, t_max, n_steps_t));
        }

        return characteristics;
    }

    /**
     * @brief Example: Burgers' equation u_t + u u_x = 0 (inviscid)
     */
    static QuasilinearPDE burgersEquation() {
        return QuasilinearPDE(
            [](double x, double y, double u) { return u; },    // a = u
            [](double x, double y, double u) { return 1.0; },  // b = 1 (time direction)
            [](double x, double y, double u) { return 0.0; }   // c = 0
        );
    }
};

/**
 * @class NonlinearFirstOrder
 * @brief First-order fully nonlinear PDEs
 *
 * General form: F(x, y, u, u_x, u_y) = 0
 */
class NonlinearFirstOrder {
public:
    using NonlinearFunction = std::function<double(double, double, double, double, double)>;

    /**
     * @brief Charpit's equations for F(x, y, u, p, q) = 0 where p = u_x, q = u_y
     *
     * dx/dt = F_p
     * dy/dt = F_q
     * du/dt = p F_p + q F_q
     * dp/dt = -F_x - p F_u
     * dq/dt = -F_y - q F_u
     */
    struct CharpitSystem {
        NonlinearFunction F;

        CharpitSystem(NonlinearFunction func) : F(func) {}

        /**
         * @brief Numerical partial derivatives
         */
        std::array<double, 5> gradient(double x, double y, double u, double p, double q,
                                        double h = 1e-6) const {
            double F_x = (F(x+h, y, u, p, q) - F(x-h, y, u, p, q)) / (2*h);
            double F_y = (F(x, y+h, u, p, q) - F(x, y-h, u, p, q)) / (2*h);
            double F_u = (F(x, y, u+h, p, q) - F(x, y, u-h, p, q)) / (2*h);
            double F_p = (F(x, y, u, p+h, q) - F(x, y, u, p-h, q)) / (2*h);
            double F_q = (F(x, y, u, p, q+h) - F(x, y, u, p, q-h)) / (2*h);

            return {F_x, F_y, F_u, F_p, F_q};
        }

        /**
         * @brief Solve Charpit system
         */
        struct CharacteristicStrip {
            std::vector<double> x, y, u, p, q, t;
        };

        CharacteristicStrip solveCharpit(
            double x0, double y0, double u0, double p0, double q0,
            double t_max, int n_steps) const {

            CharacteristicStrip strip;
            double dt = t_max / n_steps;

            double x = x0, y = y0, u = u0, p = p0, q = q0;

            for (int i = 0; i <= n_steps; ++i) {
                strip.x.push_back(x);
                strip.y.push_back(y);
                strip.u.push_back(u);
                strip.p.push_back(p);
                strip.q.push_back(q);
                strip.t.push_back(i * dt);

                if (i < n_steps) {
                    auto [F_x, F_y, F_u, F_p, F_q] = gradient(x, y, u, p, q);

                    // Euler step for Charpit equations
                    x += dt * F_p;
                    y += dt * F_q;
                    u += dt * (p * F_p + q * F_q);
                    p += dt * (-F_x - p * F_u);
                    q += dt * (-F_y - q * F_u);
                }
            }

            return strip;
        }
    };

    /**
     * @brief Eikonal equation: (u_x)² + (u_y)² = 1
     */
    static CharpitSystem eikonalEquation() {
        return CharpitSystem([](double x, double y, double u, double p, double q) {
            return p*p + q*q - 1.0;
        });
    }

    /**
     * @brief Hamilton-Jacobi equation: u_t + H(x, u_x) = 0
     */
    static CharpitSystem hamiltonJacobi(std::function<double(double, double)> H) {
        return CharpitSystem([H](double x, double y, double u, double p, double q) {
            // y plays role of time t
            // F(x, t, u, p, q) = q + H(x, p)
            return q + H(x, p);
        });
    }
};

/**
 * @class GeometricalInterpretation
 * @brief Geometrical considerations
 */
class GeometricalInterpretation {
public:
    /**
     * @brief Integral surface: u = u(x, y)
     *
     * Solution represented as surface in (x, y, u) space
     */
    struct IntegralSurface {
        ScalarField2D u;

        /**
         * @brief Normal vector to surface at (x, y)
         */
        std::array<double, 3> normal(double x, double y, double h = 1e-5) const {
            double u_x = (u(x+h, y) - u(x-h, y)) / (2*h);
            double u_y = (u(x, y+h) - u(x, y-h)) / (2*h);

            // Normal: (-u_x, -u_y, 1)
            double norm = std::sqrt(u_x*u_x + u_y*u_y + 1.0);
            return {-u_x/norm, -u_y/norm, 1.0/norm};
        }

        /**
         * @brief Tangent plane at (x0, y0): u = u0 + u_x(x-x0) + u_y(y-y0)
         */
        double tangentPlane(double x0, double y0, double x, double y, double h = 1e-5) const {
            double u0 = u(x0, y0);
            double u_x = (u(x0+h, y0) - u(x0-h, y0)) / (2*h);
            double u_y = (u(x0, y0+h) - u(x0, y0-h)) / (2*h);

            return u0 + u_x * (x - x0) + u_y * (y - y0);
        }
    };

    /**
     * @brief Monge cone at point (x0, y0, u0)
     *
     * For F(x, y, u, p, q) = 0, cone of admissible tangent planes
     */
    struct MongeCone {
        std::function<double(double, double, double, double, double)> F;
        double x0, y0, u0;

        /**
         * @brief Check if direction (p, q) lies on Monge cone
         */
        bool isOnCone(double p, double q) const {
            return std::abs(F(x0, y0, u0, p, q)) < 1e-6;
        }

        /**
         * @brief Get generators of Monge cone (characteristic directions)
         */
        std::vector<std::pair<double, double>> generators(int n_angles = 100) const {
            std::vector<std::pair<double, double>> dirs;

            // Sample angles and solve F(x0, y0, u0, p, q) = 0
            for (int i = 0; i < n_angles; ++i) {
                double theta = 2 * M_PI * i / n_angles;
                double p = std::cos(theta);
                double q = std::sin(theta);

                if (isOnCone(p, q)) {
                    dirs.push_back({p, q});
                }
            }

            return dirs;
        }
    };
};

/**
 * @class SecondOrderCharacteristics
 * @brief Characteristics for second-order PDEs
 */
class SecondOrderCharacteristics {
public:
    /**
     * @brief General second-order linear PDE:
     * A u_xx + 2B u_xy + C u_yy + D u_x + E u_y + F u = G
     *
     * Characteristic curves satisfy: A(dy)² - 2B dx dy + C(dx)² = 0
     */
    struct LinearSecondOrder {
        ScalarField2D A, B, C, D, E, F, G;

        LinearSecondOrder(ScalarField2D a, ScalarField2D b, ScalarField2D c,
                         ScalarField2D d, ScalarField2D e, ScalarField2D f, ScalarField2D g)
            : A(a), B(b), C(c), D(d), E(e), F(f), G(g) {}

        /**
         * @brief Characteristic equation: A(dy/dx)² - 2B(dy/dx) + C = 0
         *
         * Solve for dy/dx along characteristics
         */
        std::vector<double> characteristicSlopes(double x, double y) const {
            double a = A(x, y);
            double b = B(x, y);
            double c = C(x, y);

            double discriminant = b*b - a*c;

            if (std::abs(discriminant) < 1e-10) {
                // Parabolic: one family of characteristics
                return {b / a};
            } else if (discriminant > 0) {
                // Hyperbolic: two families
                double sqrt_disc = std::sqrt(discriminant);
                double slope1 = (b + sqrt_disc) / a;
                double slope2 = (b - sqrt_disc) / a;
                return {slope1, slope2};
            } else {
                // Elliptic: complex characteristics (no real characteristics)
                return {};
            }
        }

        /**
         * @brief Integrate characteristic curves
         */
        std::vector<std::pair<double, double>> characteristicCurve(
            double x0, double y0, int family, double t_max, int n_steps) const {

            std::vector<std::pair<double, double>> curve;
            double dt = t_max / n_steps;
            double x = x0, y = y0;

            for (int i = 0; i <= n_steps; ++i) {
                curve.push_back({x, y});

                if (i < n_steps) {
                    auto slopes = characteristicSlopes(x, y);
                    if (family < static_cast<int>(slopes.size())) {
                        double dy_dx = slopes[family];
                        // Euler step: assuming dx/dt = 1
                        x += dt;
                        y += dt * dy_dx;
                    } else {
                        break;  // No real characteristics
                    }
                }
            }

            return curve;
        }

        /**
         * @brief Classify at point
         */
        PDEClassification::SecondOrderType classify(double x, double y) const {
            return PDEClassification::classify(A(x, y), B(x, y), C(x, y));
        }
    };

    /**
     * @brief Canonical form via characteristic coordinates
     *
     * Transform to ξ, η coordinates where characteristics become coordinate curves
     */
    static std::string canonicalTransformation(PDEClassification::SecondOrderType type) {
        switch (type) {
            case PDEClassification::SecondOrderType::HYPERBOLIC:
                return "ξ = const along C₁, η = const along C₂ → u_ξη = F(...)";
            case PDEClassification::SecondOrderType::PARABOLIC:
                return "ξ = const along characteristic → u_ηη = F(...)";
            case PDEClassification::SecondOrderType::ELLIPTIC:
                return "No real characteristics → use complex transformation";
            default:
                return "Unknown";
        }
    }
};

} // namespace maths::pde

#endif // MATHS_PDE_HPP
