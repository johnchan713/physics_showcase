#ifndef PHYSICS_ADVANCED_FLUID_DYNAMICS_COMPRESSIBLE_FLOW_HPP
#define PHYSICS_ADVANCED_FLUID_DYNAMICS_COMPRESSIBLE_FLOW_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file compressible_flow.hpp
 * @brief Compressible flow relations
 *
 * Implements:
 * - Isentropic flow relations
 * - Normal shock relations (Rankine-Hugoniot)
 * - Oblique shock relations
 * - Prandtl-Meyer expansion
 * - Stagnation properties
 * - Critical flow conditions
 */

namespace physics::advanced::fluid_dynamics {

/**
 * @class IsentropicFlow
 * @brief Isentropic (reversible adiabatic) flow relations
 *
 * For ideal gas: prho^(-gamma) = constant
 */
class IsentropicFlow {
public:
    /**
     * @brief Temperature ratio from pressure ratio
     *
     * T2/T1 = (p2/p1)^((gamma-1)/gamma)
     *
     * @param pressure_ratio p2/p1
     * @param gamma Specific heat ratio gamma = cp/cv
     * @return Temperature ratio T2/T1
     */
    static double temperatureRatio(double pressure_ratio, double gamma) {
        double exponent = (gamma - 1.0) / gamma;
        return std::pow(pressure_ratio, exponent);
    }

    /**
     * @brief Density ratio from pressure ratio
     *
     * rho2/rho1 = (p2/p1)^(1/gamma)
     *
     * @param pressure_ratio p2/p1
     * @param gamma Specific heat ratio
     * @return Density ratio rho2/rho1
     */
    static double densityRatio(double pressure_ratio, double gamma) {
        return std::pow(pressure_ratio, 1.0 / gamma);
    }

    /**
     * @brief Isentropic process equation
     *
     * p1/rho1^gamma = p2/rho2^gamma
     *
     * @param p1 Initial pressure (Pa)
     * @param rho1 Initial density (kg/m³)
     * @param rho2 Final density (kg/m³)
     * @param gamma Specific heat ratio
     * @return Final pressure p2 (Pa)
     */
    static double pressure(double p1, double rho1, double rho2, double gamma) {
        return p1 * std::pow(rho2 / rho1, gamma);
    }

    /**
     * @brief Speed of sound
     *
     * c = sqrt(gammaRT) = sqrt(gammap/rho)
     *
     * @param gamma Specific heat ratio
     * @param pressure p (Pa)
     * @param density rho (kg/m³)
     * @return Speed of sound (m/s)
     */
    static double soundSpeed(double gamma, double pressure, double density) {
        if (density <= 0.0) {
            throw std::invalid_argument("Density must be positive");
        }
        return std::sqrt(gamma * pressure / density);
    }

    /**
     * @brief Speed of sound from temperature
     *
     * c = sqrt(gammaRT)
     *
     * @param gamma Specific heat ratio
     * @param gas_constant R (J/(kg·K))
     * @param temperature T (K)
     * @return Speed of sound (m/s)
     */
    static double soundSpeedFromTemp(double gamma, double gas_constant,
                                     double temperature) {
        return std::sqrt(gamma * gas_constant * temperature);
    }
};

/**
 * @class StagnationProperties
 * @brief Stagnation (total) properties in compressible flow
 *
 * Properties when flow is isentropically brought to rest
 */
class StagnationProperties {
public:
    /**
     * @brief Stagnation temperature
     *
     * T0 = T(1 + (gamma-1)/2 × M^2)
     *
     * @param static_temperature T (K)
     * @param mach_number M
     * @param gamma Specific heat ratio
     * @return Stagnation temperature T0 (K)
     */
    static double temperature(double static_temperature,
                             double mach_number,
                             double gamma) {
        double factor = 1.0 + 0.5 * (gamma - 1.0) * mach_number * mach_number;
        return static_temperature * factor;
    }

    /**
     * @brief Stagnation pressure
     *
     * p0 = p(1 + (gamma-1)/2 × M^2)^(gamma/(gamma-1))
     *
     * @param static_pressure p (Pa)
     * @param mach_number M
     * @param gamma Specific heat ratio
     * @return Stagnation pressure p0 (Pa)
     */
    static double pressure(double static_pressure,
                          double mach_number,
                          double gamma) {
        double factor = 1.0 + 0.5 * (gamma - 1.0) * mach_number * mach_number;
        double exponent = gamma / (gamma - 1.0);
        return static_pressure * std::pow(factor, exponent);
    }

    /**
     * @brief Stagnation density
     *
     * rho0 = rho(1 + (gamma-1)/2 × M^2)^(1/(gamma-1))
     *
     * @param static_density rho (kg/m³)
     * @param mach_number M
     * @param gamma Specific heat ratio
     * @return Stagnation density rho0 (kg/m³)
     */
    static double density(double static_density,
                         double mach_number,
                         double gamma) {
        double factor = 1.0 + 0.5 * (gamma - 1.0) * mach_number * mach_number;
        double exponent = 1.0 / (gamma - 1.0);
        return static_density * std::pow(factor, exponent);
    }

    /**
     * @brief Static temperature from stagnation
     *
     * T = T0 / (1 + (gamma-1)/2 × M^2)
     *
     * @param stagnation_temperature T0 (K)
     * @param mach_number M
     * @param gamma Specific heat ratio
     * @return Static temperature T (K)
     */
    static double staticTemperature(double stagnation_temperature,
                                   double mach_number,
                                   double gamma) {
        double factor = 1.0 + 0.5 * (gamma - 1.0) * mach_number * mach_number;
        return stagnation_temperature / factor;
    }

    /**
     * @brief Static pressure from stagnation
     *
     * p = p0 / (1 + (gamma-1)/2 × M^2)^(gamma/(gamma-1))
     *
     * @param stagnation_pressure p0 (Pa)
     * @param mach_number M
     * @param gamma Specific heat ratio
     * @return Static pressure p (Pa)
     */
    static double staticPressure(double stagnation_pressure,
                                double mach_number,
                                double gamma) {
        double factor = 1.0 + 0.5 * (gamma - 1.0) * mach_number * mach_number;
        double exponent = gamma / (gamma - 1.0);
        return stagnation_pressure / std::pow(factor, exponent);
    }
};

/**
 * @class CriticalFlow
 * @brief Critical flow conditions (M = 1)
 */
class CriticalFlow {
public:
    /**
     * @brief Critical pressure ratio
     *
     * p* / p0 = (2/(gamma+1))^(gamma/(gamma-1))
     *
     * @param gamma Specific heat ratio
     * @return Critical pressure ratio
     */
    static double criticalPressureRatio(double gamma) {
        double base = 2.0 / (gamma + 1.0);
        double exponent = gamma / (gamma - 1.0);
        return std::pow(base, exponent);
    }

    /**
     * @brief Critical temperature ratio
     *
     * T* / T0 = 2/(gamma+1)
     *
     * @param gamma Specific heat ratio
     * @return Critical temperature ratio
     */
    static double criticalTemperatureRatio(double gamma) {
        return 2.0 / (gamma + 1.0);
    }

    /**
     * @brief Critical density ratio
     *
     * rho* / rho0 = (2/(gamma+1))^(1/(gamma-1))
     *
     * @param gamma Specific heat ratio
     * @return Critical density ratio
     */
    static double criticalDensityRatio(double gamma) {
        double base = 2.0 / (gamma + 1.0);
        double exponent = 1.0 / (gamma - 1.0);
        return std::pow(base, exponent);
    }

    /**
     * @brief Critical velocity (sonic)
     *
     * c* = sqrt(2gammaRT0/(gamma+1))
     *
     * @param gamma Specific heat ratio
     * @param gas_constant R (J/(kg·K))
     * @param stagnation_temperature T0 (K)
     * @return Critical velocity (m/s)
     */
    static double criticalVelocity(double gamma, double gas_constant,
                                   double stagnation_temperature) {
        return std::sqrt(2.0 * gamma * gas_constant * stagnation_temperature /
                        (gamma + 1.0));
    }
};

/**
 * @class NormalShock
 * @brief Normal shock wave relations (Rankine-Hugoniot)
 *
 * Discontinuous jump in flow properties across shock
 */
class NormalShock {
public:
    /**
     * @brief Downstream Mach number
     *
     * M2^2 = [M1^2 + 2/(gamma-1)] / [2gammaM1^2/(gamma-1) - 1]
     *
     * @param mach_upstream M1
     * @param gamma Specific heat ratio
     * @return Downstream Mach number M2
     */
    static double downstreamMach(double mach_upstream, double gamma) {
        if (mach_upstream < 1.0) {
            throw std::invalid_argument("Upstream Mach must be supersonic (M > 1)");
        }

        double M1_sq = mach_upstream * mach_upstream;
        double gm1 = gamma - 1.0;

        double numerator = M1_sq + 2.0 / gm1;
        double denominator = 2.0 * gamma * M1_sq / gm1 - 1.0;

        return std::sqrt(numerator / denominator);
    }

    /**
     * @brief Pressure ratio across shock
     *
     * p2/p1 = 1 + 2gamma/(gamma+1) × (M1^2 - 1)
     *
     * @param mach_upstream M1
     * @param gamma Specific heat ratio
     * @return Pressure ratio p2/p1
     */
    static double pressureRatio(double mach_upstream, double gamma) {
        double M1_sq = mach_upstream * mach_upstream;
        return 1.0 + 2.0 * gamma * (M1_sq - 1.0) / (gamma + 1.0);
    }

    /**
     * @brief Density ratio across shock
     *
     * rho2/rho1 = [(gamma+1)M1^2] / [2 + (gamma-1)M1^2]
     *
     * @param mach_upstream M1
     * @param gamma Specific heat ratio
     * @return Density ratio rho2/rho1
     */
    static double densityRatio(double mach_upstream, double gamma) {
        double M1_sq = mach_upstream * mach_upstream;
        double numerator = (gamma + 1.0) * M1_sq;
        double denominator = 2.0 + (gamma - 1.0) * M1_sq;
        return numerator / denominator;
    }

    /**
     * @brief Temperature ratio across shock
     *
     * T2/T1 = (p2/p1) × (rho1/rho2)
     *
     * @param mach_upstream M1
     * @param gamma Specific heat ratio
     * @return Temperature ratio T2/T1
     */
    static double temperatureRatio(double mach_upstream, double gamma) {
        double p_ratio = pressureRatio(mach_upstream, gamma);
        double rho_ratio = densityRatio(mach_upstream, gamma);
        return p_ratio / rho_ratio;
    }

    /**
     * @brief Stagnation pressure ratio across shock
     *
     * p02/p01 < 1 (entropy increase)
     *
     * p02/p01 = [(gamma+1)M1^2/(2+(gamma-1)M1^2)]^(gamma/(gamma-1)) ×
     *           [(gamma+1)/(2gammaM1^2-(gamma-1))]^(1/(gamma-1))
     *
     * @param mach_upstream M1
     * @param gamma Specific heat ratio
     * @return Stagnation pressure ratio p02/p01
     */
    static double stagnationPressureRatio(double mach_upstream, double gamma) {
        double M1_sq = mach_upstream * mach_upstream;
        double gp1 = gamma + 1.0;
        double gm1 = gamma - 1.0;

        double term1_base = (gp1 * M1_sq) / (2.0 + gm1 * M1_sq);
        double term1 = std::pow(term1_base, gamma / gm1);

        double term2_base = gp1 / (2.0 * gamma * M1_sq - gm1);
        double term2 = std::pow(term2_base, 1.0 / gm1);

        return term1 * term2;
    }

    /**
     * @brief Entropy increase across shock
     *
     * Δs = cp ln[(T2/T1)(p1/p2)^((gamma-1)/gamma)]
     *
     * @param mach_upstream M1
     * @param gamma Specific heat ratio
     * @param specific_heat_p cp (J/(kg·K))
     * @return Entropy increase Δs (J/(kg·K))
     */
    static double entropyIncrease(double mach_upstream, double gamma,
                                 double specific_heat_p) {
        double T_ratio = temperatureRatio(mach_upstream, gamma);
        double p_ratio = pressureRatio(mach_upstream, gamma);

        double factor = T_ratio * std::pow(1.0 / p_ratio, (gamma - 1.0) / gamma);
        return specific_heat_p * std::log(factor);
    }

    /**
     * @brief Shock strength parameter
     *
     * β = (p2 - p1) / p1
     *
     * @param mach_upstream M1
     * @param gamma Specific heat ratio
     * @return Shock strength
     */
    static double shockStrength(double mach_upstream, double gamma) {
        return pressureRatio(mach_upstream, gamma) - 1.0;
    }

    /**
     * @brief Check if shock is strong
     *
     * Strong shock: M1 >> 1
     * Limiting ratios:
     * - rho2/rho1 → (gamma+1)/(gamma-1)
     * - p2/p1 → inf
     *
     * @param mach_upstream M1
     * @return true if strong shock (M1 > 3)
     */
    static bool isStrongShock(double mach_upstream) {
        return mach_upstream > 3.0;
    }
};

/**
 * @class ObliqueShock
 * @brief Oblique shock wave relations
 *
 * Shock at angle β to freestream
 */
class ObliqueShock {
public:
    /**
     * @brief Normal component of Mach number
     *
     * M1ₙ = M1 sin(β)
     *
     * where β is shock angle
     *
     * @param mach_upstream M1
     * @param shock_angle β (radians)
     * @return Normal Mach number M1ₙ
     */
    static double normalMach(double mach_upstream, double shock_angle) {
        return mach_upstream * std::sin(shock_angle);
    }

    /**
     * @brief Shock angle from deflection angle (θ-β-M relation)
     *
     * tan(θ) = 2cot(β) × (M1^2sin^2β - 1) / (M1^2(gamma+cos2β) + 2)
     *
     * Requires iterative solution
     *
     * @param mach_upstream M1
     * @param deflection_angle θ (radians)
     * @param gamma Specific heat ratio
     * @return Shock angle β (radians)
     */
    static double shockAngle(double mach_upstream,
                            double deflection_angle,
                            double gamma) {
        // Simplified: return approximate solution
        // Full implementation requires Newton-Raphson iteration

        // Initial guess (weak shock solution)
        double beta = deflection_angle + std::asin(1.0 / mach_upstream);

        // Newton-Raphson iterations
        for (int i = 0; i < 10; ++i) {
            double M1n = normalMach(mach_upstream, beta);
            double tan_theta = 2.0 / std::tan(beta) *
                (M1n * M1n - 1.0) /
                (mach_upstream * mach_upstream * (gamma + std::cos(2.0 * beta)) + 2.0);

            double f = std::atan(tan_theta) - deflection_angle;
            if (std::abs(f) < 1e-6) break;

            // Numerical derivative
            double dbeta = 1e-6;
            double M1n_p = normalMach(mach_upstream, beta + dbeta);
            double tan_theta_p = 2.0 / std::tan(beta + dbeta) *
                (M1n_p * M1n_p - 1.0) /
                (mach_upstream * mach_upstream * (gamma + std::cos(2.0 * (beta + dbeta))) + 2.0);
            double f_p = std::atan(tan_theta_p) - deflection_angle;

            double df = (f_p - f) / dbeta;
            beta -= f / df;
        }

        return beta;
    }

    /**
     * @brief Downstream Mach number
     *
     * Apply normal shock relations to normal component,
     * tangential component unchanged
     *
     * @param mach_upstream M1
     * @param shock_angle β (radians)
     * @param deflection_angle θ (radians)
     * @param gamma Specific heat ratio
     * @return Downstream Mach number M2
     */
    static double downstreamMach(double mach_upstream,
                                double shock_angle,
                                double deflection_angle,
                                double gamma) {
        double M1n = normalMach(mach_upstream, shock_angle);
        double M2n = NormalShock::downstreamMach(M1n, gamma);

        // M2 = M2ₙ / sin(β - θ)
        return M2n / std::sin(shock_angle - deflection_angle);
    }

    /**
     * @brief Maximum deflection angle
     *
     * Beyond this, shock detaches
     *
     * @param mach_upstream M1
     * @param gamma Specific heat ratio
     * @return Maximum deflection angle θmax (radians)
     */
    static double maxDeflectionAngle(double mach_upstream, double gamma) {
        // Approximate formula
        double M1_sq = mach_upstream * mach_upstream;
        double numerator = M1_sq - 1.0;
        double denominator = M1_sq * (gamma + 1.0) / 2.0;
        return std::atan(std::sqrt(numerator / denominator));
    }
};

/**
 * @class PrandtlMeyerExpansion
 * @brief Prandtl-Meyer expansion fan
 *
 * Isentropic expansion around convex corner
 */
class PrandtlMeyerExpansion {
public:
    /**
     * @brief Prandtl-Meyer function
     *
     * ν(M) = sqrt[(gamma+1)/(gamma-1)] × arctan[sqrt((gamma-1)/(gamma+1) × (M^2-1))]
     *        - arctan[sqrt(M^2-1)]
     *
     * @param mach_number M
     * @param gamma Specific heat ratio
     * @return Prandtl-Meyer angle ν (radians)
     */
    static double prandtlMeyerAngle(double mach_number, double gamma) {
        if (mach_number < 1.0) {
            throw std::invalid_argument("Mach must be supersonic (M > 1)");
        }

        double M_sq = mach_number * mach_number;
        double gp1 = gamma + 1.0;
        double gm1 = gamma - 1.0;

        double term1 = std::sqrt(gp1 / gm1) *
                      std::atan(std::sqrt((gm1 / gp1) * (M_sq - 1.0)));

        double term2 = std::atan(std::sqrt(M_sq - 1.0));

        return term1 - term2;
    }

    /**
     * @brief Downstream Mach from deflection
     *
     * ν2 = ν1 + θ
     *
     * Then solve for M2 from ν2
     *
     * @param mach_upstream M1
     * @param deflection_angle θ (radians, positive for expansion)
     * @param gamma Specific heat ratio
     * @return Downstream Mach number M2
     */
    static double downstreamMach(double mach_upstream,
                                double deflection_angle,
                                double gamma) {
        double nu1 = prandtlMeyerAngle(mach_upstream, gamma);
        double nu2 = nu1 + deflection_angle;

        // Inverse Prandtl-Meyer function (iterative)
        double M2 = mach_upstream;  // Initial guess

        for (int i = 0; i < 20; ++i) {
            double nu_M2 = prandtlMeyerAngle(M2, gamma);
            double error = nu_M2 - nu2;

            if (std::abs(error) < 1e-8) break;

            // Numerical derivative
            double dM = 0.001;
            double nu_plus = prandtlMeyerAngle(M2 + dM, gamma);
            double dnu_dM = (nu_plus - nu_M2) / dM;

            M2 -= error / dnu_dM;
        }

        return M2;
    }

    /**
     * @brief Maximum Prandtl-Meyer angle
     *
     * ν_max as M → inf
     *
     * ν_max = (π/2) × (sqrt[(gamma+1)/(gamma-1)] - 1)
     *
     * @param gamma Specific heat ratio
     * @return Maximum angle (radians)
     */
    static double maxPrandtlMeyerAngle(double gamma) {
        return (M_PI / 2.0) * (std::sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0);
    }
};

} // namespace physics::advanced::fluid_dynamics

#endif // PHYSICS_ADVANCED_FLUID_DYNAMICS_COMPRESSIBLE_FLOW_HPP
