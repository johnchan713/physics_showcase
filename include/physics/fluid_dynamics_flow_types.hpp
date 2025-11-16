#ifndef PHYSICS_ADVANCED_FLUID_DYNAMICS_FLOW_TYPES_HPP
#define PHYSICS_ADVANCED_FLUID_DYNAMICS_FLOW_TYPES_HPP

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <complex>

/**
 * @file flow_types.hpp
 * @brief Classical analytical flow solutions
 *
 * Implements exact solutions for fundamental flow types:
 * - Poiseuille Flow (pressure-driven pipe flow)
 * - Couette Flow (shear-driven flow between plates)
 * - Stokes Flow (creeping flow, Re << 1)
 * - Potential Flow (irrotational, inviscid)
 */

namespace physics::advanced::fluid_dynamics {

/**
 * @class PoiseuilleFlow
 * @brief Pressure-driven laminar flow in pipes and channels
 *
 * Classic solution for fully developed laminar flow
 */
class PoiseuilleFlow {
public:
    /**
     * @brief Velocity profile for circular pipe (Hagen-Poiseuille)
     *
     * u(r) = (ΔP/4μL)(R² - r²)
     *
     * Parabolic velocity profile with maximum at centerline
     *
     * @param radius r - Radial position (m), 0 ≤ r ≤ R
     * @param pipe_radius R - Pipe radius (m)
     * @param pressure_gradient dp/dx (Pa/m)
     * @param dynamic_viscosity μ (Pa·s)
     * @return Axial velocity u(r) (m/s)
     */
    static double velocityCircularPipe(double radius,
                                       double pipe_radius,
                                       double pressure_gradient,
                                       double dynamic_viscosity) {
        if (radius > pipe_radius) {
            throw std::invalid_argument("Radius exceeds pipe radius");
        }
        if (dynamic_viscosity <= 0.0) {
            throw std::invalid_argument("Viscosity must be positive");
        }

        double R2 = pipe_radius * pipe_radius;
        double r2 = radius * radius;

        return (-pressure_gradient / (4.0 * dynamic_viscosity)) * (R2 - r2);
    }

    /**
     * @brief Maximum velocity at centerline
     *
     * u_max = ΔP R²/(4μL)
     *
     * @param pipe_radius R
     * @param pressure_drop ΔP
     * @param length L
     * @param dynamic_viscosity μ
     * @return Maximum velocity (m/s)
     */
    static double maxVelocity(double pipe_radius,
                             double pressure_drop,
                             double length,
                             double dynamic_viscosity) {
        if (length <= 0.0) {
            throw std::invalid_argument("Length must be positive");
        }

        return (pressure_drop * pipe_radius * pipe_radius) /
               (4.0 * dynamic_viscosity * length);
    }

    /**
     * @brief Average (bulk) velocity
     *
     * u_avg = u_max/2 = ΔP R²/(8μL)
     *
     * @param pipe_radius R
     * @param pressure_drop ΔP
     * @param length L
     * @param dynamic_viscosity μ
     * @return Average velocity (m/s)
     */
    static double averageVelocity(double pipe_radius,
                                  double pressure_drop,
                                  double length,
                                  double dynamic_viscosity) {
        return 0.5 * maxVelocity(pipe_radius, pressure_drop,
                                length, dynamic_viscosity);
    }

    /**
     * @brief Volumetric flow rate (Hagen-Poiseuille equation)
     *
     * Q = πR⁴ΔP/(8μL)
     *
     * @param pipe_radius R
     * @param pressure_drop ΔP
     * @param length L
     * @param dynamic_viscosity μ
     * @return Flow rate Q (m³/s)
     */
    static double flowRate(double pipe_radius,
                          double pressure_drop,
                          double length,
                          double dynamic_viscosity) {
        double R4 = std::pow(pipe_radius, 4.0);
        return (M_PI * R4 * pressure_drop) / (8.0 * dynamic_viscosity * length);
    }

    /**
     * @brief Wall shear stress
     *
     * τ_w = μ(du/dr)|_wall = -ΔP R/(2L)
     *
     * @param pipe_radius R
     * @param pressure_drop ΔP
     * @param length L
     * @return Wall shear stress (Pa)
     */
    static double wallShearStress(double pipe_radius,
                                  double pressure_drop,
                                  double length) {
        return (pressure_drop * pipe_radius) / (2.0 * length);
    }

    /**
     * @brief Friction factor for laminar pipe flow
     *
     * f = 64/Re
     *
     * @param reynolds_number Re
     * @return Darcy friction factor
     */
    static double frictionFactor(double reynolds_number) {
        if (reynolds_number <= 0.0) {
            throw std::invalid_argument("Reynolds number must be positive");
        }
        return 64.0 / reynolds_number;
    }

    /**
     * @brief Velocity profile for plane Poiseuille flow (2D channel)
     *
     * u(y) = (ΔP/2μL)(h² - y²)
     *
     * where y is measured from centerline, -h ≤ y ≤ h
     *
     * @param y Distance from centerline
     * @param half_height h (half channel height)
     * @param pressure_gradient dp/dx
     * @param dynamic_viscosity μ
     * @return Velocity u(y)
     */
    static double velocityPlaneChannel(double y,
                                       double half_height,
                                       double pressure_gradient,
                                       double dynamic_viscosity) {
        if (std::abs(y) > half_height) {
            throw std::invalid_argument("y exceeds channel half-height");
        }

        double h2 = half_height * half_height;
        double y2 = y * y;

        return (-pressure_gradient / (2.0 * dynamic_viscosity)) * (h2 - y2);
    }
};

/**
 * @class CouetteFlow
 * @brief Shear-driven flow between parallel plates
 */
class CouetteFlow {
public:
    /**
     * @brief Simple Couette flow (no pressure gradient)
     *
     * u(y) = U(y/h)
     *
     * Linear velocity profile between plates
     *
     * @param y Position (0 ≤ y ≤ h)
     * @param plate_velocity U (velocity of moving plate)
     * @param gap_height h
     * @return Velocity u(y)
     */
    static double velocitySimple(double y,
                                 double plate_velocity,
                                 double gap_height) {
        if (y < 0.0 || y > gap_height) {
            throw std::invalid_argument("y must be between 0 and gap height");
        }
        return plate_velocity * (y / gap_height);
    }

    /**
     * @brief Shear stress (constant for simple Couette)
     *
     * τ = μU/h
     *
     * @param dynamic_viscosity μ
     * @param plate_velocity U
     * @param gap_height h
     * @return Shear stress (Pa)
     */
    static double shearStress(double dynamic_viscosity,
                             double plate_velocity,
                             double gap_height) {
        if (gap_height <= 0.0) {
            throw std::invalid_argument("Gap height must be positive");
        }
        return (dynamic_viscosity * plate_velocity) / gap_height;
    }

    /**
     * @brief Generalized Couette flow (with pressure gradient)
     *
     * u(y) = U(y/h) + (ΔP/2μL)(y² - hy)
     *
     * Combination of shear and pressure-driven flow
     *
     * @param y Position
     * @param plate_velocity U
     * @param gap_height h
     * @param pressure_gradient dp/dx
     * @param dynamic_viscosity μ
     * @return Velocity u(y)
     */
    static double velocityGeneralized(double y,
                                      double plate_velocity,
                                      double gap_height,
                                      double pressure_gradient,
                                      double dynamic_viscosity) {
        double linear_part = plate_velocity * (y / gap_height);
        double parabolic_part = (-pressure_gradient / (2.0 * dynamic_viscosity)) *
                               (y * y - gap_height * y);

        return linear_part + parabolic_part;
    }

    /**
     * @brief Taylor-Couette flow (rotating cylinders)
     *
     * Tangential velocity between concentric rotating cylinders
     *
     * u_θ(r) = Ar + B/r
     *
     * where A and B depend on boundary conditions
     *
     * @param radius r
     * @param inner_radius R₁
     * @param outer_radius R₂
     * @param omega1 Angular velocity of inner cylinder
     * @param omega2 Angular velocity of outer cylinder
     * @return Tangential velocity
     */
    static double taylorCouetteVelocity(double radius,
                                        double inner_radius,
                                        double outer_radius,
                                        double omega1,
                                        double omega2) {
        if (radius < inner_radius || radius > outer_radius) {
            throw std::invalid_argument("Radius out of bounds");
        }

        double R1_sq = inner_radius * inner_radius;
        double R2_sq = outer_radius * outer_radius;

        double A = (omega2 * R2_sq - omega1 * R1_sq) / (R2_sq - R1_sq);
        double B = R1_sq * R2_sq * (omega1 - omega2) / (R2_sq - R1_sq);

        return A * radius + B / radius;
    }
};

/**
 * @class StokesFlow
 * @brief Creeping flow at very low Reynolds number (Re << 1)
 *
 * Inertia is negligible compared to viscous forces
 */
class StokesFlow {
public:
    /**
     * @brief Stokes drag on a sphere
     *
     * F_D = 6πμRU
     *
     * Valid for Re < 1
     *
     * @param dynamic_viscosity μ
     * @param sphere_radius R
     * @param velocity U
     * @return Drag force (N)
     */
    static double sphereDrag(double dynamic_viscosity,
                            double sphere_radius,
                            double velocity) {
        return 6.0 * M_PI * dynamic_viscosity * sphere_radius * velocity;
    }

    /**
     * @brief Drag coefficient for sphere in Stokes flow
     *
     * C_D = 24/Re
     *
     * @param reynolds_number Re = ρUD/μ
     * @return Drag coefficient
     */
    static double sphereDragCoefficient(double reynolds_number) {
        if (reynolds_number <= 0.0) {
            throw std::invalid_argument("Reynolds number must be positive");
        }
        if (reynolds_number > 1.0) {
            // Warning: Stokes approximation may not be accurate
        }
        return 24.0 / reynolds_number;
    }

    /**
     * @brief Terminal velocity of falling sphere
     *
     * U_t = (2/9)(R²g/ν)(ρ_p - ρ_f)/ρ_f
     *
     * Balance between gravity and Stokes drag
     *
     * @param sphere_radius R
     * @param particle_density ρ_p
     * @param fluid_density ρ_f
     * @param kinematic_viscosity ν
     * @param gravity g
     * @return Terminal velocity (m/s)
     */
    static double terminalVelocity(double sphere_radius,
                                   double particle_density,
                                   double fluid_density,
                                   double kinematic_viscosity,
                                   double gravity = 9.81) {
        double density_diff = particle_density - fluid_density;
        double R2 = sphere_radius * sphere_radius;

        return (2.0 * R2 * gravity * density_diff) / (9.0 * fluid_density * kinematic_viscosity);
    }

    /**
     * @brief Oseen correction to Stokes drag (Re ~ 1)
     *
     * F_D = 6πμRU(1 + 3Re/16)
     *
     * @param dynamic_viscosity μ
     * @param sphere_radius R
     * @param velocity U
     * @param reynolds_number Re
     * @return Corrected drag force
     */
    static double oseenDrag(double dynamic_viscosity,
                           double sphere_radius,
                           double velocity,
                           double reynolds_number) {
        double stokes = sphereDrag(dynamic_viscosity, sphere_radius, velocity);
        double correction = 1.0 + (3.0 * reynolds_number / 16.0);
        return stokes * correction;
    }

    /**
     * @brief Stokes drag on cylinder (per unit length)
     *
     * F_D/L = 4πμU/ln(7.4/Re)
     *
     * For flow perpendicular to cylinder axis
     *
     * @param dynamic_viscosity μ
     * @param velocity U
     * @param reynolds_number Re
     * @return Drag force per unit length (N/m)
     */
    static double cylinderDrag(double dynamic_viscosity,
                              double velocity,
                              double reynolds_number) {
        if (reynolds_number <= 0.0 || reynolds_number > 1.0) {
            throw std::invalid_argument("Reynolds number must be 0 < Re < 1");
        }
        return (4.0 * M_PI * dynamic_viscosity * velocity) / std::log(7.4 / reynolds_number);
    }
};

/**
 * @class PotentialFlow
 * @brief Irrotational, inviscid flow (∇×u = 0, μ = 0)
 *
 * Governed by Laplace equation: ∇²φ = 0
 * where u = ∇φ (velocity potential)
 */
class PotentialFlow {
public:
    /**
     * @brief Uniform flow
     *
     * φ = Ux
     * u = U, v = 0
     *
     * @param position Position vector
     * @param velocity Uniform velocity
     * @return Velocity at position
     */
    static Eigen::Vector2d uniformFlow(const Eigen::Vector2d& position,
                                       const Eigen::Vector2d& velocity) {
        return velocity;
    }

    /**
     * @brief Point source/sink
     *
     * φ = (m/2π)ln(r)
     * u_r = m/(2πr), u_θ = 0
     *
     * @param position (x, y)
     * @param source_strength m (positive = source, negative = sink)
     * @return Velocity
     */
    static Eigen::Vector2d pointSource(const Eigen::Vector2d& position,
                                       double source_strength) {
        double r = position.norm();
        if (r < 1e-10) {
            throw std::domain_error("Singularity at origin");
        }

        double u_r = source_strength / (2.0 * M_PI * r);
        return u_r * position.normalized();
    }

    /**
     * @brief Point vortex
     *
     * φ = (Γ/2π)θ
     * u_r = 0, u_θ = Γ/(2πr)
     *
     * @param position (x, y)
     * @param circulation Γ
     * @return Velocity
     */
    static Eigen::Vector2d pointVortex(const Eigen::Vector2d& position,
                                       double circulation) {
        double r = position.norm();
        if (r < 1e-10) {
            throw std::domain_error("Singularity at origin");
        }

        double u_theta = circulation / (2.0 * M_PI * r);

        // Velocity is tangent to circles: perpendicular to radius
        Eigen::Vector2d tangent(-position.y(), position.x());
        return (u_theta / r) * tangent;
    }

    /**
     * @brief Doublet (dipole)
     *
     * φ = (K/2π)(x/(x²+y²))
     *
     * @param position (x, y)
     * @param dipole_strength K
     * @return Velocity
     */
    static Eigen::Vector2d doublet(const Eigen::Vector2d& position,
                                   double dipole_strength) {
        double x = position.x();
        double y = position.y();
        double r2 = x * x + y * y;

        if (r2 < 1e-10) {
            throw std::domain_error("Singularity at origin");
        }

        double factor = dipole_strength / (2.0 * M_PI * r2);

        double u = factor * (x * x - y * y) / r2;
        double v = factor * (2.0 * x * y) / r2;

        return Eigen::Vector2d(u, v);
    }

    /**
     * @brief Flow past cylinder (uniform + doublet)
     *
     * Velocity field for flow around circular cylinder
     *
     * @param position (x, y)
     * @param uniform_velocity U_∞
     * @param cylinder_radius R
     * @return Velocity
     */
    static Eigen::Vector2d cylinderFlow(const Eigen::Vector2d& position,
                                        double uniform_velocity,
                                        double cylinder_radius) {
        double r = position.norm();

        if (r < cylinder_radius) {
            // Inside cylinder
            return Eigen::Vector2d::Zero();
        }

        // Doublet strength K = -U R²
        double K = -uniform_velocity * cylinder_radius * cylinder_radius;

        Eigen::Vector2d uniform(uniform_velocity, 0.0);
        Eigen::Vector2d doublet_vel = doublet(position, K);

        return uniform + doublet_vel;
    }

    /**
     * @brief Lift on cylinder with circulation (Kutta-Joukowski)
     *
     * L = ρU∞Γ (per unit span)
     *
     * @param density ρ
     * @param velocity U∞
     * @param circulation Γ
     * @return Lift per unit span (N/m)
     */
    static double liftKuttaJoukowski(double density,
                                     double velocity,
                                     double circulation) {
        return density * velocity * circulation;
    }

    /**
     * @brief Stream function for 2D potential flow
     *
     * For irrotational flow: ψ exists such that
     * u = ∂ψ/∂y, v = -∂ψ/∂x
     *
     * Streamlines: ψ = constant
     */
    static double streamFunction(const Eigen::Vector2d& position,
                                 double source_strength,
                                 double circulation) {
        double x = position.x();
        double y = position.y();
        double r = position.norm();
        double theta = std::atan2(y, x);

        if (r < 1e-10) {
            return 0.0;
        }

        // Source contribution
        double psi_source = (source_strength / (2.0 * M_PI)) * theta;

        // Vortex contribution
        double psi_vortex = -(circulation / (2.0 * M_PI)) * std::log(r);

        return psi_source + psi_vortex;
    }
};

} // namespace physics::advanced::fluid_dynamics

#endif // PHYSICS_ADVANCED_FLUID_DYNAMICS_FLOW_TYPES_HPP
