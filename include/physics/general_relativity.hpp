#ifndef PHYSICS_GENERAL_RELATIVITY_HPP
#define PHYSICS_GENERAL_RELATIVITY_HPP

#include <array>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace physics {
namespace general_relativity {

using Vector4 = std::array<double, 4>;
using Tensor2 = std::array<std::array<double, 4>, 4>;
using Tensor3 = std::array<std::array<std::array<double, 4>, 4>, 4>;
using Tensor4 = std::array<std::array<std::array<std::array<double, 4>, 4>, 4>, 4>;

constexpr double c = 299792458.0;
constexpr double G = 6.67430e-11;

// ============================================================================
// METRIC TENSORS
// ============================================================================

class MetricTensor {
public:
    virtual Tensor2 covariant(const Vector4& x) const = 0;
    virtual ~MetricTensor() = default;

    Tensor2 contravariant(const Vector4& x) const {
        auto g = covariant(x);
        Tensor2 g_inv;

        double det = determinant(g);
        if (std::abs(det) < 1e-15) {
            throw std::runtime_error("Singular metric");
        }

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                g_inv[mu][nu] = cofactor(g, mu, nu) / det;
            }
        }
        return g_inv;
    }

    double determinant(const Tensor2& g) const {
        return g[0][0] * (g[1][1] * (g[2][2] * g[3][3] - g[2][3] * g[3][2])
                         - g[1][2] * (g[2][1] * g[3][3] - g[2][3] * g[3][1])
                         + g[1][3] * (g[2][1] * g[3][2] - g[2][2] * g[3][1]))
             - g[0][1] * (g[1][0] * (g[2][2] * g[3][3] - g[2][3] * g[3][2])
                         - g[1][2] * (g[2][0] * g[3][3] - g[2][3] * g[3][0])
                         + g[1][3] * (g[2][0] * g[3][2] - g[2][2] * g[3][0]))
             + g[0][2] * (g[1][0] * (g[2][1] * g[3][3] - g[2][3] * g[3][1])
                         - g[1][1] * (g[2][0] * g[3][3] - g[2][3] * g[3][0])
                         + g[1][3] * (g[2][0] * g[3][1] - g[2][1] * g[3][0]))
             - g[0][3] * (g[1][0] * (g[2][1] * g[3][2] - g[2][2] * g[3][1])
                         - g[1][1] * (g[2][0] * g[3][2] - g[2][2] * g[3][0])
                         + g[1][2] * (g[2][0] * g[3][1] - g[2][1] * g[3][0]));
    }

    double cofactor(const Tensor2& g, int row, int col) const {
        double minor[3][3];
        int r = 0;
        for (int i = 0; i < 4; ++i) {
            if (i == row) continue;
            int c = 0;
            for (int j = 0; j < 4; ++j) {
                if (j == col) continue;
                minor[r][c] = g[i][j];
                ++c;
            }
            ++r;
        }
        double det3 = minor[0][0] * (minor[1][1] * minor[2][2] - minor[1][2] * minor[2][1])
                    - minor[0][1] * (minor[1][0] * minor[2][2] - minor[1][2] * minor[2][0])
                    + minor[0][2] * (minor[1][0] * minor[2][1] - minor[1][1] * minor[2][0]);
        return ((row + col) % 2 == 0 ? 1 : -1) * det3;
    }

    double lineElement(const Vector4& x, const Vector4& dx) const {
        auto g = covariant(x);
        double ds2 = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                ds2 += g[mu][nu] * dx[mu] * dx[nu];
            }
        }
        return ds2;
    }

    double properTime(const Vector4& x, const Vector4& four_velocity) const {
        return std::sqrt(-lineElement(x, four_velocity));
    }
};

class MinkowskiMetric : public MetricTensor {
public:
    Tensor2 covariant(const Vector4& x) const override {
        Tensor2 eta = {{{-1, 0, 0, 0},
                        {0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}}};
        return eta;
    }
};

class SchwarzschildMetric : public MetricTensor {
private:
    double M_;

public:
    SchwarzschildMetric(double mass) : M_(mass) {}

    Tensor2 covariant(const Vector4& x) const override {
        double r = x[1];
        double theta = x[2];
        double r_s = 2.0 * G * M_ / (c * c);

        if (r <= r_s * 1.001) {
            r = r_s * 1.001;
        }

        double f = 1.0 - r_s / r;
        double sin_theta = std::sin(theta);

        Tensor2 g = {{
            {-f, 0, 0, 0},
            {0, 1.0 / f, 0, 0},
            {0, 0, r * r, 0},
            {0, 0, 0, r * r * sin_theta * sin_theta}
        }};
        return g;
    }

    double schwarzschildRadius() const {
        return 2.0 * G * M_ / (c * c);
    }

    double photonSphereRadius() const {
        return 1.5 * schwarzschildRadius();
    }

    double innerMostStableOrbit() const {
        return 3.0 * schwarzschildRadius();
    }
};

class KerrMetric : public MetricTensor {
private:
    double M_;
    double a_;

public:
    KerrMetric(double mass, double spin_parameter) : M_(mass), a_(spin_parameter) {}

    Tensor2 covariant(const Vector4& x) const override {
        double r = x[1];
        double theta = x[2];

        double r_s = 2.0 * G * M_ / (c * c);
        double rho2 = r * r + a_ * a_ * std::cos(theta) * std::cos(theta);
        double Delta = r * r - r_s * r + a_ * a_;
        double sin_theta = std::sin(theta);

        Tensor2 g;
        g[0][0] = -(1.0 - r_s * r / rho2);
        g[0][3] = -r_s * r * a_ * sin_theta * sin_theta / rho2;
        g[3][0] = g[0][3];
        g[1][1] = rho2 / Delta;
        g[2][2] = rho2;
        g[3][3] = (r * r + a_ * a_ + r_s * r * a_ * a_ * sin_theta * sin_theta / rho2)
                 * sin_theta * sin_theta;

        g[0][1] = g[0][2] = 0.0;
        g[1][0] = g[1][2] = g[1][3] = 0.0;
        g[2][0] = g[2][1] = g[2][3] = 0.0;
        g[3][1] = g[3][2] = 0.0;

        return g;
    }

    double ergoSphereRadius(double theta) const {
        double r_s = 2.0 * G * M_ / (c * c);
        return 0.5 * r_s * (1.0 + std::sqrt(1.0 - (a_ / (0.5 * r_s)) * (a_ / (0.5 * r_s))
                                             * std::cos(theta) * std::cos(theta)));
    }
};

class FRWMetric : public MetricTensor {
private:
    std::function<double(double)> scale_factor_;
    double k_;

public:
    FRWMetric(std::function<double(double)> a, double curvature)
        : scale_factor_(a), k_(curvature) {}

    Tensor2 covariant(const Vector4& x) const override {
        double t = x[0];
        double r = x[1];
        double theta = x[2];
        double a = scale_factor_(t);
        double sin_theta = std::sin(theta);

        double spatial_factor;
        if (std::abs(k_) < 1e-10) {
            spatial_factor = 1.0;
        } else if (k_ > 0) {
            spatial_factor = 1.0 / (1.0 - k_ * r * r);
        } else {
            spatial_factor = 1.0 / (1.0 + std::abs(k_) * r * r);
        }

        Tensor2 g = {{
            {-c * c, 0, 0, 0},
            {0, a * a * spatial_factor, 0, 0},
            {0, 0, a * a * r * r, 0},
            {0, 0, 0, a * a * r * r * sin_theta * sin_theta}
        }};
        return g;
    }
};

// ============================================================================
// CHRISTOFFEL SYMBOLS
// ============================================================================

class ChristoffelSymbols {
private:
    const MetricTensor& metric_;
    double epsilon_;

public:
    ChristoffelSymbols(const MetricTensor& g, double eps = 1e-6)
        : metric_(g), epsilon_(eps) {}

    double component(int lambda, int mu, int nu, const Vector4& x) const {
        auto g_inv = metric_.contravariant(x);
        double Gamma = 0.0;

        for (int sigma = 0; sigma < 4; ++sigma) {
            Gamma += 0.5 * g_inv[lambda][sigma] * (
                partial_derivative(sigma, mu, nu, x)
                + partial_derivative(sigma, nu, mu, x)
                - partial_derivative(mu, nu, sigma, x)
            );
        }
        return Gamma;
    }

    Tensor3 allComponents(const Vector4& x) const {
        Tensor3 Gamma;
        for (int lambda = 0; lambda < 4; ++lambda) {
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    Gamma[lambda][mu][nu] = component(lambda, mu, nu, x);
                }
            }
        }
        return Gamma;
    }

private:
    double partial_derivative(int alpha, int mu, int nu, const Vector4& x) const {
        Vector4 x_plus = x;
        Vector4 x_minus = x;
        x_plus[alpha] += epsilon_;
        x_minus[alpha] -= epsilon_;

        auto g_plus = metric_.covariant(x_plus);
        auto g_minus = metric_.covariant(x_minus);

        return (g_plus[mu][nu] - g_minus[mu][nu]) / (2.0 * epsilon_);
    }
};

// ============================================================================
// RIEMANN CURVATURE TENSOR
// ============================================================================

class RiemannTensor {
private:
    ChristoffelSymbols christoffel_;
    double epsilon_;

public:
    RiemannTensor(const ChristoffelSymbols& gamma, double eps = 1e-6)
        : christoffel_(gamma), epsilon_(eps) {}

    double component(int rho, int sigma, int mu, int nu,
                    const Vector4& x) const {

        Vector4 x_mu_plus = x, x_mu_minus = x;
        x_mu_plus[mu] += epsilon_;
        x_mu_minus[mu] -= epsilon_;
        double dGamma_mu = (christoffel_.component(rho, nu, sigma, x_mu_plus)
                          - christoffel_.component(rho, nu, sigma, x_mu_minus))
                          / (2.0 * epsilon_);

        Vector4 x_nu_plus = x, x_nu_minus = x;
        x_nu_plus[nu] += epsilon_;
        x_nu_minus[nu] -= epsilon_;
        double dGamma_nu = (christoffel_.component(rho, mu, sigma, x_nu_plus)
                          - christoffel_.component(rho, mu, sigma, x_nu_minus))
                          / (2.0 * epsilon_);

        double Gamma_Gamma1 = 0.0;
        double Gamma_Gamma2 = 0.0;
        for (int lambda = 0; lambda < 4; ++lambda) {
            Gamma_Gamma1 += christoffel_.component(rho, mu, lambda, x)
                          * christoffel_.component(lambda, nu, sigma, x);
            Gamma_Gamma2 += christoffel_.component(rho, nu, lambda, x)
                          * christoffel_.component(lambda, mu, sigma, x);
        }

        return dGamma_mu - dGamma_nu + Gamma_Gamma1 - Gamma_Gamma2;
    }

    Tensor4 allComponents(const Vector4& x) const {
        Tensor4 R;
        for (int rho = 0; rho < 4; ++rho) {
            for (int sigma = 0; sigma < 4; ++sigma) {
                for (int mu = 0; mu < 4; ++mu) {
                    for (int nu = 0; nu < 4; ++nu) {
                        R[rho][sigma][mu][nu] = component(rho, sigma, mu, nu, x);
                    }
                }
            }
        }
        return R;
    }
};

// ============================================================================
// RICCI TENSOR AND SCALAR
// ============================================================================

class RicciTensor {
private:
    RiemannTensor riemann_;

public:
    RicciTensor(const RiemannTensor& R) : riemann_(R) {}

    Tensor2 components(const Vector4& x) const {
        Tensor2 Ric;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                Ric[mu][nu] = 0.0;
                for (int lambda = 0; lambda < 4; ++lambda) {
                    Ric[mu][nu] += riemann_.component(lambda, mu, lambda, nu, x);
                }
            }
        }
        return Ric;
    }

    double scalar(const Vector4& x, const MetricTensor& g) const {
        auto Ric = components(x);
        auto g_inv = g.contravariant(x);
        double R = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                R += g_inv[mu][nu] * Ric[mu][nu];
            }
        }
        return R;
    }
};

// ============================================================================
// EINSTEIN TENSOR
// ============================================================================

class EinsteinTensor {
private:
    RicciTensor ricci_;
    const MetricTensor& metric_;

public:
    EinsteinTensor(const RicciTensor& Ric, const MetricTensor& g)
        : ricci_(Ric), metric_(g) {}

    Tensor2 components(const Vector4& x) const {
        auto Ric = ricci_.components(x);
        double R = ricci_.scalar(x, metric_);
        auto g = metric_.covariant(x);

        Tensor2 G;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                G[mu][nu] = Ric[mu][nu] - 0.5 * g[mu][nu] * R;
            }
        }
        return G;
    }

    bool verifyFieldEquations(const Tensor2& stress_energy,
                              const Vector4& x,
                              double tol = 1e-6) const {
        auto G_tensor = components(x);
        double factor = 8.0 * M_PI * G / (c * c * c * c);

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                double diff = std::abs(G_tensor[mu][nu] - factor * stress_energy[mu][nu]);
                if (diff > tol) {
                    return false;
                }
            }
        }
        return true;
    }
};

// ============================================================================
// GEODESIC EQUATIONS
// ============================================================================

class GeodesicSolver {
private:
    const MetricTensor& metric_;
    ChristoffelSymbols christoffel_;

public:
    struct State {
        Vector4 position;
        Vector4 four_velocity;
    };

    GeodesicSolver(const MetricTensor& g)
        : metric_(g), christoffel_(g) {}

    State derivative(const State& s) const {
        State ds;
        ds.position = s.four_velocity;

        for (int mu = 0; mu < 4; ++mu) {
            ds.four_velocity[mu] = 0.0;
            for (int alpha = 0; alpha < 4; ++alpha) {
                for (int beta = 0; beta < 4; ++beta) {
                    double Gamma = christoffel_.component(mu, alpha, beta, s.position);
                    ds.four_velocity[mu] -= Gamma * s.four_velocity[alpha]
                                                  * s.four_velocity[beta];
                }
            }
        }
        return ds;
    }

    std::vector<State> integrate(const State& initial,
                                 double dtau, int steps) const {
        std::vector<State> trajectory;
        trajectory.push_back(initial);

        State current = initial;
        for (int step = 0; step < steps; ++step) {
            State k1 = derivative(current);

            State temp;
            for (int i = 0; i < 4; ++i) {
                temp.position[i] = current.position[i] + 0.5 * dtau * k1.position[i];
                temp.four_velocity[i] = current.four_velocity[i]
                                       + 0.5 * dtau * k1.four_velocity[i];
            }
            State k2 = derivative(temp);

            for (int i = 0; i < 4; ++i) {
                temp.position[i] = current.position[i] + 0.5 * dtau * k2.position[i];
                temp.four_velocity[i] = current.four_velocity[i]
                                       + 0.5 * dtau * k2.four_velocity[i];
            }
            State k3 = derivative(temp);

            for (int i = 0; i < 4; ++i) {
                temp.position[i] = current.position[i] + dtau * k3.position[i];
                temp.four_velocity[i] = current.four_velocity[i]
                                       + dtau * k3.four_velocity[i];
            }
            State k4 = derivative(temp);

            for (int i = 0; i < 4; ++i) {
                current.position[i] += (dtau / 6.0) * (k1.position[i]
                                     + 2 * k2.position[i]
                                     + 2 * k3.position[i]
                                     + k4.position[i]);
                current.four_velocity[i] += (dtau / 6.0) * (k1.four_velocity[i]
                                          + 2 * k2.four_velocity[i]
                                          + 2 * k3.four_velocity[i]
                                          + k4.four_velocity[i]);
            }

            trajectory.push_back(current);
        }
        return trajectory;
    }

    double perihelionShift(double M_sun, double semi_major_axis,
                          double eccentricity, int orbits = 1) const {
        double r_s = 2.0 * G * M_sun / (c * c);
        double epsilon = 6.0 * M_PI * G * M_sun / (c * c * semi_major_axis
                                                   * (1.0 - eccentricity * eccentricity));
        return epsilon * orbits * (180.0 / M_PI) * 3600.0;
    }
};

// ============================================================================
// GRAVITATIONAL WAVES
// ============================================================================

class GravitationalWaves {
public:
    static Tensor2 linearizedMetric(const Vector4& x, double amplitude,
                                    double frequency, int polarization) {
        Tensor2 h;
        for (auto& row : h) row.fill(0.0);

        double omega = 2.0 * M_PI * frequency;
        double phase = omega * (x[0] / c - x[3] / c);

        if (polarization == 0) {
            h[1][1] = amplitude * std::cos(phase);
            h[2][2] = -amplitude * std::cos(phase);
        } else {
            h[1][2] = amplitude * std::cos(phase);
            h[2][1] = amplitude * std::cos(phase);
        }
        return h;
    }

    static double strainAmplitude(double M1, double M2, double distance,
                                 double frequency) {
        double M_chirp = std::pow(M1 * M2, 0.6) / std::pow(M1 + M2, 0.2);
        double h = (4.0 / distance) * std::pow(G * M_chirp * M_PI * frequency / c / c, 2.0 / 3.0);
        return h / c;
    }
};

// ============================================================================
// SCHWARZSCHILD SOLUTION PROPERTIES
// ============================================================================

class SchwarzschildSolution {
private:
    double M_;
    double r_s_;

public:
    SchwarzschildSolution(double mass)
        : M_(mass), r_s_(2.0 * G * mass / (c * c)) {}

    double orbitalVelocity(double r) const {
        if (r <= r_s_) return 0.0;
        return std::sqrt(G * M_ / r);
    }

    double escapVelocity(double r) const {
        if (r <= r_s_) return c;
        return std::sqrt(2.0 * G * M_ / r);
    }

    double properTimeAtRadius(double r, double coordinate_time) const {
        if (r <= r_s_) return 0.0;
        double f = std::sqrt(1.0 - r_s_ / r);
        return f * coordinate_time;
    }

    double redshift(double r_emit, double r_obs) const {
        if (r_emit <= r_s_ || r_obs <= r_s_) {
            return std::numeric_limits<double>::infinity();
        }
        double f_emit = std::sqrt(1.0 - r_s_ / r_emit);
        double f_obs = std::sqrt(1.0 - r_s_ / r_obs);
        return (f_obs / f_emit) - 1.0;
    }

    double photonOrbitRadius() const {
        return 1.5 * r_s_;
    }

    double tidalForce(double r, double mass_object, double length) const {
        if (r <= r_s_) return 0.0;
        return 2.0 * G * M_ * mass_object * length / (r * r * r);
    }
};

} // namespace general_relativity
} // namespace physics

#endif // PHYSICS_GENERAL_RELATIVITY_HPP
