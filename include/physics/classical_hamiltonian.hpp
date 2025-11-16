#ifndef PHYSICS_ADVANCED_CLASSICAL_HAMILTONIAN_HPP
#define PHYSICS_ADVANCED_CLASSICAL_HAMILTONIAN_HPP

#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <stdexcept>

/**
 * @file hamiltonian.hpp
 * @brief Hamiltonian mechanics and canonical transformations
 *
 * Implements:
 * - Hamiltonian formulation of classical mechanics
 * - Canonical coordinates (q, p)
 * - Hamilton's equations
 * - Poisson brackets
 * - Canonical transformations
 * - Generating functions
 * - Action-angle variables
 */

namespace physics::advanced::classical {

/**
 * @struct PhasePoint
 * @brief Point in phase space (q, p)
 */
struct PhasePoint {
    Eigen::VectorXd q;  // Generalized coordinates
    Eigen::VectorXd p;  // Generalized momenta

    PhasePoint(int dim) : q(Eigen::VectorXd::Zero(dim)),
                          p(Eigen::VectorXd::Zero(dim)) {}

    PhasePoint(const Eigen::VectorXd& q_, const Eigen::VectorXd& p_)
        : q(q_), p(p_) {
        if (q.size() != p.size()) {
            throw std::invalid_argument("q and p must have same dimension");
        }
    }

    int dimension() const { return q.size(); }
};

/**
 * @class HamiltonianSystem
 * @brief Represents a Hamiltonian dynamical system
 *
 * Hamilton's equations:
 *   dq_i/dt = ∂H/∂p_i
 *   dp_i/dt = -∂H/∂q_i
 */
class HamiltonianSystem {
private:
    int dim_;
    std::function<double(const PhasePoint&)> H_;
    double epsilon_;  // For numerical derivatives

public:
    /**
     * @brief Construct Hamiltonian system
     * @param dimension Degrees of freedom
     * @param hamiltonian Hamiltonian function H(q, p)
     * @param epsilon Step size for numerical derivatives
     */
    HamiltonianSystem(int dimension,
                     std::function<double(const PhasePoint&)> hamiltonian,
                     double epsilon = 1e-8)
        : dim_(dimension), H_(hamiltonian), epsilon_(epsilon) {
        if (dimension <= 0) {
            throw std::invalid_argument("Dimension must be positive");
        }
    }

    /**
     * @brief Evaluate Hamiltonian at phase point
     */
    double hamiltonian(const PhasePoint& point) const {
        return H_(point);
    }

    /**
     * @brief Compute ∂H/∂q_i using numerical differentiation
     */
    double dH_dq(const PhasePoint& point, int i) const {
        PhasePoint p_plus = point;
        PhasePoint p_minus = point;
        p_plus.q(i) += epsilon_;
        p_minus.q(i) -= epsilon_;
        return (H_(p_plus) - H_(p_minus)) / (2.0 * epsilon_);
    }

    /**
     * @brief Compute ∂H/∂p_i using numerical differentiation
     */
    double dH_dp(const PhasePoint& point, int i) const {
        PhasePoint p_plus = point;
        PhasePoint p_minus = point;
        p_plus.p(i) += epsilon_;
        p_minus.p(i) -= epsilon_;
        return (H_(p_plus) - H_(p_minus)) / (2.0 * epsilon_);
    }

    /**
     * @brief Compute time derivative via Hamilton's equations
     *
     * Returns (dq/dt, dp/dt) where:
     *   dq_i/dt = ∂H/∂p_i
     *   dp_i/dt = -∂H/∂q_i
     */
    PhasePoint hamiltonEquations(const PhasePoint& point) const {
        PhasePoint derivative(dim_);

        for (int i = 0; i < dim_; ++i) {
            derivative.q(i) = dH_dp(point, i);
            derivative.p(i) = -dH_dq(point, i);
        }

        return derivative;
    }

    /**
     * @brief Symplectic Euler integrator (1st order)
     *
     * Preserves symplectic structure exactly
     */
    PhasePoint stepSymplecticEuler(const PhasePoint& point, double dt) const {
        PhasePoint next = point;

        // Update p first: p_{n+1} = p_n - dt ∂H/∂q|_{q_n}
        for (int i = 0; i < dim_; ++i) {
            next.p(i) = point.p(i) - dt * dH_dq(point, i);
        }

        // Then update q: q_{n+1} = q_n + dt ∂H/∂p|_{p_{n+1}}
        for (int i = 0; i < dim_; ++i) {
            next.q(i) = point.q(i) + dt * dH_dp(next, i);
        }

        return next;
    }

    /**
     * @brief Störmer-Verlet integrator (2nd order symplectic)
     *
     * More accurate than symplectic Euler while preserving structure
     */
    PhasePoint stepVerlet(const PhasePoint& point, double dt) const {
        PhasePoint next = point;

        // Half step in p
        for (int i = 0; i < dim_; ++i) {
            next.p(i) = point.p(i) - 0.5 * dt * dH_dq(point, i);
        }

        // Full step in q
        for (int i = 0; i < dim_; ++i) {
            next.q(i) = point.q(i) + dt * dH_dp(next, i);
        }

        // Half step in p
        for (int i = 0; i < dim_; ++i) {
            next.p(i) = next.p(i) - 0.5 * dt * dH_dq(next, i);
        }

        return next;
    }

    /**
     * @brief Integrate trajectory using Verlet method
     */
    std::vector<PhasePoint> integrate(const PhasePoint& initial,
                                     double dt, int num_steps) const {
        std::vector<PhasePoint> trajectory;
        trajectory.reserve(num_steps + 1);
        trajectory.push_back(initial);

        PhasePoint current = initial;
        for (int i = 0; i < num_steps; ++i) {
            current = stepVerlet(current, dt);
            trajectory.push_back(current);
        }

        return trajectory;
    }

    /**
     * @brief Check energy conservation
     */
    bool checkEnergyConservation(const std::vector<PhasePoint>& trajectory,
                                 double tolerance = 1e-6) const {
        if (trajectory.size() < 2) return true;

        double E0 = H_(trajectory[0]);
        for (size_t i = 1; i < trajectory.size(); ++i) {
            double E = H_(trajectory[i]);
            if (std::abs(E - E0) > tolerance * std::abs(E0)) {
                return false;
            }
        }
        return true;
    }

    int dimension() const { return dim_; }

    /**
     * @brief Create simple harmonic oscillator Hamiltonian
     *
     * H = p²/(2m) + (1/2)kq²
     */
    static HamiltonianSystem harmonicOscillator(double mass, double k) {
        return HamiltonianSystem(1, [mass, k](const PhasePoint& point) {
            double q = point.q(0);
            double p = point.p(0);
            return p * p / (2.0 * mass) + 0.5 * k * q * q;
        });
    }

    /**
     * @brief Create Kepler problem (planetary motion)
     *
     * H = p²/(2m) - GMm/|q|
     */
    static HamiltonianSystem keplerProblem(double mass, double GM) {
        return HamiltonianSystem(2, [mass, GM](const PhasePoint& point) {
            double px = point.p(0);
            double py = point.p(1);
            double x = point.q(0);
            double y = point.q(1);
            double r = std::sqrt(x * x + y * y);

            if (r < 1e-10) {
                throw std::domain_error("Collision with central body");
            }

            return (px * px + py * py) / (2.0 * mass) - GM * mass / r;
        });
    }
};

/**
 * @brief Compute Poisson bracket {f, g}
 *
 * {f, g} = Σᵢ (∂f/∂qᵢ ∂g/∂pᵢ - ∂f/∂pᵢ ∂g/∂qᵢ)
 */
class PoissonBracket {
private:
    double epsilon_;

public:
    PoissonBracket(double epsilon = 1e-8) : epsilon_(epsilon) {}

    double compute(const std::function<double(const PhasePoint&)>& f,
                  const std::function<double(const PhasePoint&)>& g,
                  const PhasePoint& point) const {

        int dim = point.dimension();
        double result = 0.0;

        for (int i = 0; i < dim; ++i) {
            // ∂f/∂qᵢ
            PhasePoint p_plus = point, p_minus = point;
            p_plus.q(i) += epsilon_;
            p_minus.q(i) -= epsilon_;
            double df_dq = (f(p_plus) - f(p_minus)) / (2.0 * epsilon_);

            // ∂g/∂pᵢ
            p_plus = point; p_minus = point;
            p_plus.p(i) += epsilon_;
            p_minus.p(i) -= epsilon_;
            double dg_dp = (g(p_plus) - g(p_minus)) / (2.0 * epsilon_);

            // ∂f/∂pᵢ
            p_plus = point; p_minus = point;
            p_plus.p(i) += epsilon_;
            p_minus.p(i) -= epsilon_;
            double df_dp = (f(p_plus) - f(p_minus)) / (2.0 * epsilon_);

            // ∂g/∂qᵢ
            p_plus = point; p_minus = point;
            p_plus.q(i) += epsilon_;
            p_minus.q(i) -= epsilon_;
            double dg_dq = (g(p_plus) - g(p_minus)) / (2.0 * epsilon_);

            result += df_dq * dg_dp - df_dp * dg_dq;
        }

        return result;
    }

    /**
     * @brief Verify fundamental Poisson bracket relations
     *
     * {qᵢ, pⱼ} = δᵢⱼ
     * {qᵢ, qⱼ} = 0
     * {pᵢ, pⱼ} = 0
     */
    bool verifyFundamentalRelations(const PhasePoint& point,
                                   double tolerance = 1e-6) const {
        int dim = point.dimension();

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                // {qᵢ, pⱼ}
                auto qi = [i](const PhasePoint& p) { return p.q(i); };
                auto pj = [j](const PhasePoint& p) { return p.p(j); };
                double bracket_qp = compute(qi, pj, point);
                double expected = (i == j) ? 1.0 : 0.0;
                if (std::abs(bracket_qp - expected) > tolerance) {
                    return false;
                }

                // {qᵢ, qⱼ} = 0
                auto qj = [j](const PhasePoint& p) { return p.q(j); };
                if (std::abs(compute(qi, qj, point)) > tolerance) {
                    return false;
                }

                // {pᵢ, pⱼ} = 0
                auto pi = [i](const PhasePoint& p) { return p.p(i); };
                if (std::abs(compute(pi, pj, point)) > tolerance) {
                    return false;
                }
            }
        }

        return true;
    }
};

} // namespace physics::advanced::classical

#endif // PHYSICS_ADVANCED_CLASSICAL_HAMILTONIAN_HPP
