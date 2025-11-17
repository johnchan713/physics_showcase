#ifndef PHYSICS_ADVANCED_FLUID_DYNAMICS_BOUNDARY_LAYER_HPP
#define PHYSICS_ADVANCED_FLUID_DYNAMICS_BOUNDARY_LAYER_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file boundary_layer.hpp
 * @brief Boundary layer theory for viscous flows
 *
 * Implements:
 * - Blasius solution (flat plate boundary layer)
 * - Boundary layer thickness definitions
 * - von Kármán momentum integral
 * - Skin friction coefficients
 */

namespace physics::advanced::fluid_dynamics {

/**
 * @class BlasiusSolution
 * @brief Similarity solution for laminar flat plate boundary layer
 *
 * Self-similar solution reduces PDE to ODE:
 * f''' + (1/2)ff'' = 0
 * where eta = ysqrt(U/(nux))
 */
class BlasiusSolution {
public:
    /**
     * @brief Similarity variable eta
     *
     * eta = ysqrt(Uinf/(nux))
     *
     * @param y Distance from wall (m)
     * @param x Distance from leading edge (m)
     * @param velocity Uinf (m/s)
     * @param kinematic_viscosity nu (m^2/s)
     * @return Similarity variable eta
     */
    static double similarityVariable(double y, double x,
                                     double velocity,
                                     double kinematic_viscosity) {
        if (x <= 0.0) {
            throw std::invalid_argument("x must be positive");
        }
        return y * std::sqrt(velocity / (kinematic_viscosity * x));
    }

    /**
     * @brief Approximate velocity profile from Blasius solution
     *
     * u/Uinf = f'(eta)
     *
     * Using polynomial approximation for f'(eta):
     * f' ≈ 2eta - 2eta³/3 + eta⁴/6  for eta ≤ 3
     * f' ≈ 1                   for eta > 3
     *
     * @param eta Similarity variable eta
     * @return Dimensionless velocity u/Uinf
     */
    static double velocityProfile(double eta) {
        if (eta < 0.0) {
            return 0.0;
        }
        if (eta > 5.0) {
            return 1.0;  // Edge of boundary layer
        }

        // Polynomial approximation (valid for small eta)
        if (eta <= 3.0) {
            return 2.0 * eta - (2.0/3.0) * std::pow(eta, 3) +
                   (1.0/6.0) * std::pow(eta, 4);
        }

        // Transition region
        return 1.0 - std::exp(-eta + 3.0);
    }

    /**
     * @brief Boundary layer thickness delta99
     *
     * delta99 = 5.0sqrt(nux/Uinf)
     *
     * Distance where u = 0.99Uinf
     *
     * @param x Distance from leading edge
     * @param velocity Uinf
     * @param kinematic_viscosity nu
     * @return Boundary layer thickness (m)
     */
    static double thickness(double x, double velocity,
                           double kinematic_viscosity) {
        return 5.0 * std::sqrt(kinematic_viscosity * x / velocity);
    }

    /**
     * @brief Displacement thickness delta*
     *
     * delta* = integral0^inf (1 - u/Uinf) dy = 1.721sqrt(nux/Uinf)
     *
     * Represents mass flux deficit
     *
     * @param x Distance from leading edge
     * @param velocity Uinf
     * @param kinematic_viscosity nu
     * @return Displacement thickness (m)
     */
    static double displacementThickness(double x, double velocity,
                                       double kinematic_viscosity) {
        return 1.721 * std::sqrt(kinematic_viscosity * x / velocity);
    }

    /**
     * @brief Momentum thickness theta
     *
     * theta = integral0^inf (u/Uinf)(1 - u/Uinf) dy = 0.664sqrt(nux/Uinf)
     *
     * Represents momentum flux deficit
     *
     * @param x Distance from leading edge
     * @param velocity Uinf
     * @param kinematic_viscosity nu
     * @return Momentum thickness (m)
     */
    static double momentumThickness(double x, double velocity,
                                   double kinematic_viscosity) {
        return 0.664 * std::sqrt(kinematic_viscosity * x / velocity);
    }

    /**
     * @brief Shape factor H
     *
     * H = delta* / theta
     *
     * For Blasius: H = 2.59
     * Increases with adverse pressure gradient
     *
     * @param displacement_thickness delta*
     * @param momentum_thickness theta
     * @return Shape factor H
     */
    static double shapeFactor(double displacement_thickness,
                             double momentum_thickness) {
        if (momentum_thickness <= 0.0) {
            throw std::invalid_argument("Momentum thickness must be positive");
        }
        return displacement_thickness / momentum_thickness;
    }

    /**
     * @brief Local skin friction coefficient
     *
     * Cf = τw/(½ρUinf^2) = 0.664/sqrtRex
     *
     * where Rex = Uinfx/nu
     *
     * @param reynolds_x Local Reynolds number Rex
     * @return Skin friction coefficient
     */
    static double skinFrictionCoefficient(double reynolds_x) {
        if (reynolds_x <= 0.0) {
            throw std::invalid_argument("Reynolds number must be positive");
        }
        return 0.664 / std::sqrt(reynolds_x);
    }

    /**
     * @brief Wall shear stress
     *
     * τw = 0.332ρUinf^2sqrt(nu/(Uinfx))
     *
     * @param density ρ
     * @param velocity Uinf
     * @param x Distance from leading edge
     * @param kinematic_viscosity nu
     * @return Wall shear stress (Pa)
     */
    static double wallShearStress(double density, double velocity,
                                 double x, double kinematic_viscosity) {
        double factor = std::sqrt(kinematic_viscosity / (velocity * x));
        return 0.332 * density * velocity * velocity * factor;
    }

    /**
     * @brief Total drag force on flat plate (both sides)
     *
     * FD = 1.328ρUinf^2Lsqrt(nu/(UinfL)) × width
     *
     * @param density ρ
     * @param velocity Uinf
     * @param length L
     * @param width W
     * @param kinematic_viscosity nu
     * @return Total drag force (N)
     */
    static double totalDrag(double density, double velocity,
                           double length, double width,
                           double kinematic_viscosity) {
        double ReL = velocity * length / kinematic_viscosity;
        double Cf_avg = 1.328 / std::sqrt(ReL);
        double area = length * width;
        return 2.0 * Cf_avg * 0.5 * density * velocity * velocity * area;
    }
};

/**
 * @class VonKarmanIntegral
 * @brief Von Kármán momentum integral equation
 *
 * Integral form of boundary layer equations:
 * τw = ρU^2(dtheta/dx + theta(2+H)/U × dU/dx)
 */
class VonKarmanIntegral {
public:
    /**
     * @brief Von Kármán momentum integral equation
     *
     * dtheta/dx = Cf/2 - theta(2+H)/U × dU/dx
     *
     * @param theta_current Current momentum thickness theta
     * @param cf_local Local skin friction coefficient
     * @param shape_factor H = delta* / theta
     * @param velocity U
     * @param velocity_gradient dU/dx
     * @return dtheta/dx
     */
    static double momentumThicknessGradient(
        double theta_current,
        double cf_local,
        double shape_factor,
        double velocity,
        double velocity_gradient) {

        double term1 = cf_local / 2.0;
        double term2 = (theta_current / velocity) *
                       (2.0 + shape_factor) * velocity_gradient;

        return term1 - term2;
    }

    /**
     * @brief Approximate Cf for laminar flow
     *
     * Cf ≈ 2(dtheta/dx)  (for zero pressure gradient)
     *
     * @param theta_gradient dtheta/dx
     * @return Skin friction coefficient
     */
    static double skinFrictionFromGradient(double theta_gradient) {
        return 2.0 * theta_gradient;
    }

    /**
     * @brief Separation criterion
     *
     * Separation occurs when τw → 0, which corresponds to
     * (dtheta/dx)(U^2/nu) → 0.0082 (laminar)
     *
     * @param theta_gradient dtheta/dx
     * @param velocity U
     * @param kinematic_viscosity nu
     * @return true if separation is imminent
     */
    static bool isSeparationImminent(double theta_gradient,
                                    double velocity,
                                    double kinematic_viscosity) {
        double param = (theta_gradient * velocity * velocity) /
                       kinematic_viscosity;
        return param < 0.01;  // Near separation
    }
};

/**
 * @class TurbulentBoundaryLayer
 * @brief Empirical relations for turbulent boundary layers
 */
class TurbulentBoundaryLayer {
public:
    /**
     * @brief Turbulent boundary layer thickness
     *
     * delta ≈ 0.37x/Re_x^(1/5)
     *
     * @param x Distance from leading edge
     * @param reynolds_x Local Reynolds number
     * @return Boundary layer thickness (m)
     */
    static double thickness(double x, double reynolds_x) {
        return 0.37 * x / std::pow(reynolds_x, 0.2);
    }

    /**
     * @brief Turbulent skin friction coefficient (flat plate)
     *
     * Cf = 0.059/Re_x^(1/5)  (for 5×10⁵ < Re < 10⁷)
     *
     * @param reynolds_x Local Reynolds number
     * @return Skin friction coefficient
     */
    static double skinFrictionCoefficient(double reynolds_x) {
        if (reynolds_x < 5e5) {
            throw std::invalid_argument("Reynolds number too low for turbulent");
        }
        return 0.059 / std::pow(reynolds_x, 0.2);
    }

    /**
     * @brief Average drag coefficient (turbulent flat plate)
     *
     * CD = 0.074/Re_L^(1/5) - 1700/Re_L
     *
     * Accounts for laminar leading edge
     *
     * @param reynolds_L Reynolds number based on length
     * @return Average drag coefficient
     */
    static double avgDragCoefficient(double reynolds_L) {
        return 0.074 / std::pow(reynolds_L, 0.2) - 1700.0 / reynolds_L;
    }

    /**
     * @brief Power law velocity profile for turbulent flow
     *
     * u/Uinf = (y/delta)^(1/n)
     *
     * Typically n = 7 for moderate Reynolds numbers
     *
     * @param y Distance from wall
     * @param delta Boundary layer thickness
     * @param n Power law exponent
     * @return Dimensionless velocity u/Uinf
     */
    static double powerLawProfile(double y, double delta, int n = 7) {
        if (y > delta) {
            return 1.0;
        }
        return std::pow(y / delta, 1.0 / n);
    }

    /**
     * @brief Log law velocity profile (inner layer)
     *
     * u⁺ = (1/κ)ln(y⁺) + B
     *
     * where κ ≈ 0.41 (von Kármán constant), B ≈ 5.0
     *
     * @param y_plus Dimensionless wall distance y⁺ = yuτ/nu
     * @param kappa von Kármán constant (default 0.41)
     * @param B Constant (default 5.0)
     * @return Dimensionless velocity u⁺ = u/uτ
     */
    static double logLawProfile(double y_plus,
                               double kappa = 0.41,
                               double B = 5.0) {
        if (y_plus < 5.0) {
            // Viscous sublayer: u⁺ = y⁺
            return y_plus;
        }
        if (y_plus < 30.0) {
            // Buffer layer: transition
            return 5.0 + (1.0/kappa) * std::log(y_plus/5.0);
        }
        // Log layer
        return (1.0 / kappa) * std::log(y_plus) + B;
    }

    /**
     * @brief Friction velocity uτ
     *
     * uτ = sqrt(τw/ρ)
     *
     * Velocity scale in wall coordinates
     *
     * @param wall_shear_stress τw
     * @param density ρ
     * @return Friction velocity (m/s)
     */
    static double frictionVelocity(double wall_shear_stress, double density) {
        return std::sqrt(wall_shear_stress / density);
    }

    /**
     * @brief Transition Reynolds number
     *
     * Re_crit ≈ 5×10⁵ for flat plate (varies with conditions)
     *
     * @return Critical Reynolds number for transition
     */
    static double transitionReynolds() {
        return 5e5;
    }
};

/**
 * @class BoundaryLayerParameters
 * @brief Comprehensive boundary layer parameter calculations
 */
class BoundaryLayerParameters {
public:
    /**
     * @brief Calculate all boundary layer thicknesses
     */
    struct Thicknesses {
        double delta;           // delta (99% thickness)
        double delta_star;      // delta* (displacement)
        double theta;           // theta (momentum)
        double H;               // Shape factor
    };

    /**
     * @brief Get all thicknesses for laminar Blasius solution
     */
    static Thicknesses blasiusThicknesses(double x, double velocity,
                                          double kinematic_viscosity) {
        Thicknesses t;
        t.delta = BlasiusSolution::thickness(x, velocity, kinematic_viscosity);
        t.delta_star = BlasiusSolution::displacementThickness(x, velocity, kinematic_viscosity);
        t.theta = BlasiusSolution::momentumThickness(x, velocity, kinematic_viscosity);
        t.H = t.delta_star / t.theta;  // Should be ≈ 2.59
        return t;
    }

    /**
     * @brief Get all thicknesses for turbulent boundary layer
     */
    static Thicknesses turbulentThicknesses(double x, double reynolds_x) {
        Thicknesses t;
        t.delta = TurbulentBoundaryLayer::thickness(x, reynolds_x);
        t.delta_star = 0.047 * x / std::pow(reynolds_x, 0.2);  // Empirical
        t.theta = 0.037 * x / std::pow(reynolds_x, 0.2);       // Empirical
        t.H = t.delta_star / t.theta;  // Should be ≈ 1.3 for turbulent
        return t;
    }
};

} // namespace physics::advanced::fluid_dynamics

#endif // PHYSICS_ADVANCED_FLUID_DYNAMICS_BOUNDARY_LAYER_HPP
