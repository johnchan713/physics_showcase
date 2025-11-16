#ifndef PHYSICS_CLASSICAL_FIELD_THEORY_HPP
#define PHYSICS_CLASSICAL_FIELD_THEORY_HPP

#include <array>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>

namespace physics {
namespace classical_field_theory {

using Vector4 = std::array<double, 4>;
using Tensor2 = std::array<std::array<double, 4>, 4>;

constexpr double c = 299792458.0;
constexpr double epsilon_0 = 8.8541878128e-12;
constexpr double mu_0 = 1.25663706212e-6;

// ============================================================================
// LAGRANGIAN FIELD THEORY
// ============================================================================

class ScalarField {
protected:
    std::function<double(const Vector4&)> phi_;

public:
    ScalarField(std::function<double(const Vector4&)> phi) : phi_(phi) {}

    double operator()(const Vector4& x) const { return phi_(x); }

    double partialDerivative(int mu, const Vector4& x, double eps = 1e-6) const {
        Vector4 x_plus = x, x_minus = x;
        x_plus[mu] += eps;
        x_minus[mu] -= eps;
        return (phi_(x_plus) - phi_(x_minus)) / (2.0 * eps);
    }

    Vector4 gradient(const Vector4& x, double eps = 1e-6) const {
        Vector4 grad;
        for (int mu = 0; mu < 4; ++mu) {
            grad[mu] = partialDerivative(mu, x, eps);
        }
        return grad;
    }

    double dAlembertian(const Vector4& x, double eps = 1e-6) const {
        double box = 0.0;
        box -= partialDerivative2(0, x, eps) / (c * c);
        for (int i = 1; i < 4; ++i) {
            box += partialDerivative2(i, x, eps);
        }
        return box;
    }

private:
    double partialDerivative2(int mu, const Vector4& x, double eps) const {
        Vector4 x_plus = x, x_minus = x;
        x_plus[mu] += eps;
        x_minus[mu] -= eps;
        return (phi_(x_plus) - 2.0 * phi_(x) + phi_(x_minus)) / (eps * eps);
    }
};

class LagrangianDensity {
public:
    virtual double operator()(const Vector4& x, const ScalarField& phi,
                             const Vector4& dphi) const = 0;
    virtual ~LagrangianDensity() = default;

    double action(const ScalarField& phi,
                 const std::vector<Vector4>& spacetime_points) const {
        double S = 0.0;
        double dV = 1.0;

        for (const auto& x : spacetime_points) {
            auto dphi = phi.gradient(x);
            S += (*this)(x, phi, dphi) * dV;
        }
        return S;
    }
};

class KleinGordonLagrangian : public LagrangianDensity {
private:
    double m_;

public:
    KleinGordonLagrangian(double mass) : m_(mass) {}

    double operator()(const Vector4& x, const ScalarField& phi,
                     const Vector4& dphi) const override {
        double kinetic = 0.0;
        kinetic -= dphi[0] * dphi[0] / (c * c);
        for (int i = 1; i < 4; ++i) {
            kinetic += dphi[i] * dphi[i];
        }
        kinetic *= 0.5;

        double mass_term = -0.5 * m_ * m_ * c * c * phi(x) * phi(x);

        return kinetic + mass_term;
    }

    double mass() const { return m_; }
};

class PhiFourthLagrangian : public LagrangianDensity {
private:
    double m_;
    double lambda_;

public:
    PhiFourthLagrangian(double mass, double coupling)
        : m_(mass), lambda_(coupling) {}

    double operator()(const Vector4& x, const ScalarField& phi,
                     const Vector4& dphi) const override {
        double kinetic = 0.0;
        kinetic -= dphi[0] * dphi[0] / (c * c);
        for (int i = 1; i < 4; ++i) {
            kinetic += dphi[i] * dphi[i];
        }
        kinetic *= 0.5;

        double phi_val = phi(x);
        double potential = -0.5 * m_ * m_ * phi_val * phi_val
                          - (lambda_ / 24.0) * std::pow(phi_val, 4);

        return kinetic + potential;
    }
};

// ============================================================================
// EULER-LAGRANGE EQUATIONS
// ============================================================================

class EulerLagrangeEquations {
public:
    static double fieldEquation(
        const LagrangianDensity& L,
        const ScalarField& phi,
        const Vector4& x,
        double eps = 1e-6) {

        auto dphi = phi.gradient(x, eps);
        double dL_dphi = (L(x, phi, dphi) - L(x, ScalarField(std::function<double(const Vector4&)>([&phi, &x, eps](const Vector4&) {
            return phi(x) - eps;
        })), dphi)) / eps;

        double dL_ddphi_sum = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            Vector4 x_plus = x;
            x_plus[mu] += eps;
            auto dphi_plus = phi.gradient(x_plus, eps);

            Vector4 x_minus = x;
            x_minus[mu] -= eps;
            auto dphi_minus = phi.gradient(x_minus, eps);

            double dL_ddphi_mu = (L(x_plus, phi, dphi_plus) - L(x_minus, phi, dphi_minus))
                                / (2.0 * eps);
            dL_ddphi_sum += dL_ddphi_mu;
        }

        return dL_dphi - dL_ddphi_sum;
    }

    static bool satisfiesEquation(const LagrangianDensity& L,
                                 const ScalarField& phi,
                                 const Vector4& x,
                                 double tol = 1e-3) {
        double eq = fieldEquation(L, phi, x);
        return std::abs(eq) < tol;
    }
};

// ============================================================================
// NOETHER'S THEOREM
// ============================================================================

class NoetherTheorem {
public:
    struct Symmetry {
        std::function<Vector4(const Vector4&, double)> transformation;
        double parameter;
    };

    static Vector4 noetherCurrent(
        const LagrangianDensity& L,
        const ScalarField& phi,
        const Symmetry& symmetry,
        const Vector4& x,
        double eps = 1e-6) {

        Vector4 j_mu;
        auto dphi = phi.gradient(x, eps);

        for (int mu = 0; mu < 4; ++mu) {
            Vector4 dphi_mu;
            for (int nu = 0; nu < 4; ++nu) {
                dphi_mu[nu] = dphi[nu];
            }

            j_mu[mu] = L(x, phi, dphi_mu);
        }

        return j_mu;
    }

    static double conservedCharge(
        const std::vector<Vector4>& current_density,
        const std::vector<Vector4>& spacetime_points) {

        double Q = 0.0;
        double dV = 1.0;

        for (size_t i = 0; i < current_density.size(); ++i) {
            Q += current_density[i][0] * dV;
        }
        return Q;
    }

    static bool isConserved(
        const std::vector<double>& charges_over_time,
        double tol = 1e-3) {

        if (charges_over_time.size() < 2) return true;

        double Q0 = charges_over_time[0];
        for (double Q : charges_over_time) {
            if (std::abs(Q - Q0) > tol) {
                return false;
            }
        }
        return true;
    }
};

// ============================================================================
// STRESS-ENERGY TENSOR
// ============================================================================

class StressEnergyTensor {
private:
    const LagrangianDensity& lagrangian_;

public:
    StressEnergyTensor(const LagrangianDensity& L) : lagrangian_(L) {}

    double component(int mu, int nu,
                    const ScalarField& phi,
                    const Vector4& x,
                    double eps = 1e-6) const {

        auto dphi = phi.gradient(x, eps);

        double dL_ddphi_nu = dphi[nu];
        double T_mu_nu = dL_ddphi_nu * dphi[mu];

        double eta_mu_nu = (mu == nu) ? ((mu == 0) ? -1.0 : 1.0) : 0.0;
        T_mu_nu -= eta_mu_nu * lagrangian_(x, phi, dphi);

        return T_mu_nu;
    }

    Tensor2 allComponents(const ScalarField& phi, const Vector4& x, double eps = 1e-6) const {
        Tensor2 T;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                T[mu][nu] = component(mu, nu, phi, x, eps);
            }
        }
        return T;
    }

    double energyDensity(const ScalarField& phi, const Vector4& x, double eps = 1e-6) const {
        return component(0, 0, phi, x, eps);
    }

    std::array<double, 3> momentumDensity(const ScalarField& phi, const Vector4& x,
                                          double eps = 1e-6) const {
        return {component(0, 1, phi, x, eps),
                component(0, 2, phi, x, eps),
                component(0, 3, phi, x, eps)};
    }

    double trace(const ScalarField& phi, const Vector4& x, double eps = 1e-6) const {
        double tr = 0.0;
        tr -= component(0, 0, phi, x, eps);
        for (int i = 1; i < 4; ++i) {
            tr += component(i, i, phi, x, eps);
        }
        return tr;
    }

    bool isConserved(const std::vector<Tensor2>& T_at_points,
                    const std::vector<Vector4>& points,
                    double eps = 1e-6, double tol = 1e-3) const {

        for (int mu = 0; mu < 4; ++mu) {
            for (size_t i = 1; i < points.size() - 1; ++i) {
                double div_T_mu = 0.0;
                for (int nu = 0; nu < 4; ++nu) {
                    Vector4 x_plus = points[i];
                    x_plus[nu] += eps;
                    size_t ip = i + 1;

                    Vector4 x_minus = points[i];
                    x_minus[nu] -= eps;
                    size_t im = (i > 0) ? i - 1 : i;

                    div_T_mu += (T_at_points[ip][mu][nu] - T_at_points[im][mu][nu])
                              / (2.0 * eps);
                }

                if (std::abs(div_T_mu) > tol) {
                    return false;
                }
            }
        }
        return true;
    }
};

// ============================================================================
// ELECTROMAGNETIC FIELD TENSOR
// ============================================================================

class ElectromagneticField {
private:
    std::function<std::array<double, 3>(const Vector4&)> E_;
    std::function<std::array<double, 3>(const Vector4&)> B_;

public:
    ElectromagneticField(std::function<std::array<double, 3>(const Vector4&)> E,
                        std::function<std::array<double, 3>(const Vector4&)> B)
        : E_(E), B_(B) {}

    Tensor2 fieldTensor(const Vector4& x) const {
        auto E = E_(x);
        auto B = B_(x);

        Tensor2 F;
        for (auto& row : F) row.fill(0.0);

        F[0][1] = -E[0] / c;
        F[0][2] = -E[1] / c;
        F[0][3] = -E[2] / c;

        F[1][2] = -B[2];
        F[1][3] = B[1];
        F[2][3] = -B[0];

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = mu + 1; nu < 4; ++nu) {
                F[nu][mu] = -F[mu][nu];
            }
        }

        return F;
    }

    double firstInvariant(const Vector4& x) const {
        auto E = E_(x);
        auto B = B_(x);

        double E2 = E[0]*E[0] + E[1]*E[1] + E[2]*E[2];
        double B2 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];

        return (B2 - E2 / (c * c));
    }

    double secondInvariant(const Vector4& x) const {
        auto E = E_(x);
        auto B = B_(x);

        return (E[0]*B[0] + E[1]*B[1] + E[2]*B[2]) / c;
    }

    Tensor2 stressEnergyTensor(const Vector4& x) const {
        auto F = fieldTensor(x);
        Tensor2 T;

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                T[mu][nu] = 0.0;
                for (int rho = 0; rho < 4; ++rho) {
                    double eta_rho_rho = (rho == 0) ? -1.0 : 1.0;
                    T[mu][nu] += F[mu][rho] * F[nu][rho] * eta_rho_rho / mu_0;
                }

                double eta_mu_nu = (mu == nu) ? ((mu == 0) ? -1.0 : 1.0) : 0.0;
                double F_trace = 0.0;
                for (int alpha = 0; alpha < 4; ++alpha) {
                    for (int beta = 0; beta < 4; ++beta) {
                        double eta_alpha_beta = (alpha == beta) ? ((alpha == 0) ? -1.0 : 1.0) : 0.0;
                        F_trace += F[alpha][beta] * F[alpha][beta] * eta_alpha_beta;
                    }
                }
                T[mu][nu] -= 0.25 * eta_mu_nu * F_trace / mu_0;
            }
        }
        return T;
    }

    double energyDensity(const Vector4& x) const {
        auto E = E_(x);
        auto B = B_(x);

        double E2 = E[0]*E[0] + E[1]*E[1] + E[2]*E[2];
        double B2 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];

        return 0.5 * (epsilon_0 * E2 + B2 / mu_0);
    }

    std::array<double, 3> poyntingVector(const Vector4& x) const {
        auto E = E_(x);
        auto B = B_(x);

        return {(E[1]*B[2] - E[2]*B[1]) / mu_0,
                (E[2]*B[0] - E[0]*B[2]) / mu_0,
                (E[0]*B[1] - E[1]*B[0]) / mu_0};
    }
};

// ============================================================================
// GAUGE TRANSFORMATIONS
// ============================================================================

class GaugeTransformation {
public:
    static std::array<double, 4> vectorPotential(
        const std::function<double(const Vector4&)>& phi_gauge,
        const std::array<double, 4>& A_old,
        const Vector4& x,
        double eps = 1e-6) {

        std::array<double, 4> A_new = A_old;

        for (int mu = 0; mu < 4; ++mu) {
            Vector4 x_plus = x, x_minus = x;
            x_plus[mu] += eps;
            x_minus[mu] -= eps;
            double d_phi = (phi_gauge(x_plus) - phi_gauge(x_minus)) / (2.0 * eps);
            A_new[mu] += d_phi;
        }

        return A_new;
    }

    static bool isGaugeInvariant(
        const std::array<double, 4>& A_old,
        const std::array<double, 4>& A_new,
        const Vector4& x,
        double eps = 1e-6,
        double tol = 1e-6) {

        auto F_old = computeFieldTensor(A_old, x, eps);
        auto F_new = computeFieldTensor(A_new, x, eps);

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                if (std::abs(F_old[mu][nu] - F_new[mu][nu]) > tol) {
                    return false;
                }
            }
        }
        return true;
    }

    static Tensor2 computeFieldTensor(const std::array<double, 4>& A,
                                      const Vector4& x,
                                      double eps) {
        Tensor2 F;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                F[mu][nu] = 0.0;
            }
        }
        return F;
    }

    static std::array<double, 4> coulombGauge(const std::array<double, 4>& A,
                                              const Vector4& x,
                                              double eps = 1e-6) {
        std::array<double, 4> A_coulomb = A;

        double div_A = 0.0;
        for (int i = 1; i < 4; ++i) {
            Vector4 x_plus = x, x_minus = x;
            x_plus[i] += eps;
            x_minus[i] -= eps;
            div_A += (A[i] - A[i]) / (2.0 * eps);
        }

        return A_coulomb;
    }

    static std::array<double, 4> lorenzGauge(const std::array<double, 4>& A,
                                            const Vector4& x,
                                            double eps = 1e-6) {
        std::array<double, 4> A_lorenz = A;

        double d_mu_A_mu = 0.0;
        d_mu_A_mu -= A[0] / (c * c);
        for (int i = 1; i < 4; ++i) {
            d_mu_A_mu += A[i];
        }

        return A_lorenz;
    }
};

// ============================================================================
// SYMMETRY BREAKING
// ============================================================================

class SpontaneousSymmetryBreaking {
public:
    static double mexicanHatPotential(double phi, double mu2, double lambda) {
        return -0.5 * mu2 * phi * phi + 0.25 * lambda * phi * phi * phi * phi;
    }

    static std::pair<double, double> vacuumExpectationValues(double mu2, double lambda) {
        if (mu2 <= 0.0) {
            return {0.0, 0.0};
        }

        double v = std::sqrt(mu2 / lambda);
        return {v, -v};
    }

    static double massSquaredFluctuation(double mu2, double lambda, double vev) {
        return -mu2 + 3.0 * lambda * vev * vev;
    }

    static double goldstoneBosonMass() {
        return 0.0;
    }

    static bool hasGoldstoneBosons(double mu2) {
        return mu2 > 0.0;
    }
};

} // namespace classical_field_theory
} // namespace physics

#endif // PHYSICS_CLASSICAL_FIELD_THEORY_HPP
