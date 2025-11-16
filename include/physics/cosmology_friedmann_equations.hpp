#ifndef PHYSICS_ADVANCED_COSMOLOGY_FRIEDMANN_EQUATIONS_HPP
#define PHYSICS_ADVANCED_COSMOLOGY_FRIEDMANN_EQUATIONS_HPP

#include <cmath>
#include <vector>
#include <stdexcept>

/**
 * @file friedmann_equations.hpp
 * @brief Friedmann equations and FLRW cosmology
 *
 * Implements:
 * - Friedmann equation (1st)
 * - Acceleration equation (2nd Friedmann)
 * - Fluid equation (energy conservation)
 * - Flat, open, closed geometries
 * - Critical density
 */

namespace physics::advanced::cosmology {

/**
 * @class FriedmannEquations
 * @brief Friedmann-Lemaître-Robertson-Walker (FLRW) equations
 *
 * Describe evolution of homogeneous, isotropic universe
 */
class FriedmannEquations {
public:
    /**
     * @brief First Friedmann equation
     *
     * H² = (ȧ/a)² = (8πG/3)ρ - kc²/a² + Λ/3
     *
     * where H is Hubble parameter
     *       ρ is total energy density
     *       k is curvature (-1, 0, +1)
     *       Λ is cosmological constant
     *
     * @param scale_factor a
     * @param a_dot ȧ
     * @param energy_density ρ (J/m³)
     * @param curvature k (-1, 0, or +1)
     * @param cosmological_constant Λ (m⁻²)
     * @param G Gravitational constant (m³/(kg·s²))
     * @param c Speed of light (m/s)
     * @return 0 if equation satisfied (for numerical solvers)
     */
    static double firstFriedmann(double scale_factor, double a_dot,
                                double energy_density, double curvature,
                                double cosmological_constant,
                                double G = 6.674e-11, double c = 2.998e8) {
        if (scale_factor <= 0.0) {
            throw std::invalid_argument("Scale factor must be positive");
        }

        double H = a_dot / scale_factor;
        double H_squared = H * H;

        double rho_term = (8.0 * M_PI * G / 3.0) * energy_density;
        double curv_term = -curvature * c * c / (scale_factor * scale_factor);
        double lambda_term = cosmological_constant / 3.0;

        return H_squared - rho_term - curv_term - lambda_term;
    }

    /**
     * @brief Second Friedmann equation (acceleration equation)
     *
     * ä/a = -(4πG/3)(ρ + 3p/c²) + Λ/3
     *
     * @param scale_factor a
     * @param a_double_dot ä
     * @param energy_density ρ (J/m³)
     * @param pressure p (Pa)
     * @param cosmological_constant Λ (m⁻²)
     * @return 0 if equation satisfied
     */
    static double secondFriedmann(double scale_factor, double a_double_dot,
                                 double energy_density, double pressure,
                                 double cosmological_constant,
                                 double G = 6.674e-11, double c = 2.998e8) {
        if (scale_factor <= 0.0) {
            throw std::invalid_argument("Scale factor must be positive");
        }

        double a_double_dot_over_a = a_double_dot / scale_factor;

        double rho_p_term = -(4.0 * M_PI * G / 3.0) *
                           (energy_density + 3.0 * pressure / (c * c));
        double lambda_term = cosmological_constant / 3.0;

        return a_double_dot_over_a - rho_p_term - lambda_term;
    }

    /**
     * @brief Simplified first Friedmann (flat, Λ = 0)
     *
     * H² = (8πG/3)ρ
     *
     * For flat universe without dark energy
     */
    static double hubbleFromDensity(double energy_density,
                                    double G = 6.674e-11) {
        return std::sqrt((8.0 * M_PI * G / 3.0) * energy_density);
    }

    /**
     * @brief Density parameter Ω
     *
     * Ω = ρ/ρ_crit
     *
     * where ρ_crit = 3H²/(8πG)
     */
    static double densityParameter(double energy_density,
                                   double hubble_parameter,
                                   double G = 6.674e-11) {
        double rho_crit = criticalDensity(hubble_parameter, G);
        return energy_density / rho_crit;
    }

    /**
     * @brief Critical density
     *
     * ρ_crit = 3H²/(8πG)
     *
     * Density for flat universe (k=0)
     */
    static double criticalDensity(double hubble_parameter,
                                 double G = 6.674e-11) {
        return (3.0 * hubble_parameter * hubble_parameter) / (8.0 * M_PI * G);
    }

    /**
     * @brief Critical density today
     *
     * ρ_crit,0 = 3H₀²/(8πG) ≈ 8.5×10⁻²⁷ kg/m³ ≈ 5 protons/m³
     */
    static double criticalDensityToday() {
        double H0_SI = 2.2e-18;  // s⁻¹ (H₀ ~ 67 km/s/Mpc)
        return criticalDensity(H0_SI);  // kg/m³
    }
};

/**
 * @class CurvatureGeometry
 * @brief Spatial curvature of universe
 */
class CurvatureGeometry {
public:
    /**
     * @brief Curvature parameter k
     *
     * k = +1: Positive curvature (closed, sphere)
     * k = 0:  Zero curvature (flat, Euclidean)
     * k = -1: Negative curvature (open, hyperbolic)
     */
    enum class Curvature {
        CLOSED = 1,   // Positive curvature
        FLAT = 0,     // Zero curvature
        OPEN = -1     // Negative curvature
    };

    /**
     * @brief Curvature from density parameters
     *
     * Ω_total = Ω_m + Ω_r + Ω_Λ + Ω_k
     *
     * where Ω_k = -kc²/(a²H²)
     *
     * If Ω_total = 1 → k = 0 (flat)
     * If Ω_total > 1 → k = +1 (closed)
     * If Ω_total < 1 → k = -1 (open)
     */
    static Curvature fromDensityParameters(double Omega_total,
                                          double tolerance = 1e-3) {
        if (std::abs(Omega_total - 1.0) < tolerance) {
            return Curvature::FLAT;
        } else if (Omega_total > 1.0) {
            return Curvature::CLOSED;
        } else {
            return Curvature::OPEN;
        }
    }

    /**
     * @brief Curvature density parameter
     *
     * Ω_k = 1 - Ω_total = 1 - (Ω_m + Ω_r + Ω_Λ)
     */
    static double curvatureDensityParameter(double Omega_matter,
                                           double Omega_radiation,
                                           double Omega_lambda) {
        return 1.0 - (Omega_matter + Omega_radiation + Omega_lambda);
    }

    /**
     * @brief Observed curvature
     *
     * Planck 2018: Ω_k = 0.001 ± 0.002
     * Universe is very close to flat!
     */
    static double observedCurvature() {
        return 0.001;  // Very close to zero
    }

    /**
     * @brief Spatial geometry description
     */
    static std::string geometryDescription(Curvature k) {
        switch(k) {
            case Curvature::CLOSED:
                return "Closed universe (k=+1): Finite volume, positive curvature\n"
                       "Sum of angles in triangle > 180°\n"
                       "Eventually recollapses (if no dark energy)";
            case Curvature::FLAT:
                return "Flat universe (k=0): Infinite volume, zero curvature\n"
                       "Euclidean geometry\n"
                       "Sum of angles in triangle = 180°\n"
                       "Critical density: Ω = 1";
            case Curvature::OPEN:
                return "Open universe (k=-1): Infinite volume, negative curvature\n"
                       "Hyperbolic geometry\n"
                       "Sum of angles in triangle < 180°\n"
                       "Expands forever";
            default:
                return "Unknown geometry";
        }
    }
};

/**
 * @class FluidEquation
 * @brief Energy-momentum conservation
 *
 * Continuity equation for expanding universe
 */
class FluidEquation {
public:
    /**
     * @brief Fluid equation (energy conservation)
     *
     * dρ/dt + 3(ȧ/a)(ρ + p/c²) = 0
     *
     * or equivalently:
     * d(ρa³)/dt = -p d(a³)/dt / c²
     *
     * @param energy_density ρ
     * @param pressure p
     * @param scale_factor a
     * @param a_dot ȧ
     * @param c Speed of light
     * @return dρ/dt
     */
    static double energyEvolution(double energy_density, double pressure,
                                 double scale_factor, double a_dot,
                                 double c = 2.998e8) {
        if (scale_factor <= 0.0) {
            throw std::invalid_argument("Scale factor must be positive");
        }

        double H = a_dot / scale_factor;
        return -3.0 * H * (energy_density + pressure / (c * c));
    }

    /**
     * @brief Equation of state parameter w
     *
     * p = w ρ c²
     *
     * Matter: w = 0
     * Radiation: w = 1/3
     * Dark energy: w = -1
     */
    static double equationOfStateParameter(double pressure, double energy_density,
                                          double c = 2.998e8) {
        if (energy_density <= 0.0) return 0.0;
        return pressure / (energy_density * c * c);
    }

    /**
     * @brief Energy density evolution with scale factor
     *
     * For constant w:
     * ρ(a) = ρ₀ a^(-3(1+w))
     *
     * Matter (w=0): ρ ∝ a⁻³ (volume dilution)
     * Radiation (w=1/3): ρ ∝ a⁻⁴ (volume + redshift)
     * Dark energy (w=-1): ρ = constant
     *
     * @param rho_0 Initial density
     * @param scale_factor a
     * @param w Equation of state parameter
     * @return ρ(a)
     */
    static double densityScaling(double rho_0, double scale_factor, double w) {
        double exponent = -3.0 * (1.0 + w);
        return rho_0 * std::pow(scale_factor, exponent);
    }
};

/**
 * @class AccelerationEquation
 * @brief Acceleration/deceleration of expansion
 */
class AccelerationEquation {
public:
    /**
     * @brief Deceleration parameter q
     *
     * q = -ä/(aH²)
     *
     * q > 0: Deceleration (gravity wins)
     * q < 0: Acceleration (dark energy wins)
     */
    static double decelerationParameter(double scale_factor, double a_dot,
                                       double a_double_dot) {
        if (scale_factor <= 0.0 || a_dot <= 0.0) {
            throw std::invalid_argument("Invalid scale factor or derivative");
        }

        double H = a_dot / scale_factor;
        return -a_double_dot / (scale_factor * H * H);
    }

    /**
     * @brief Deceleration parameter from density parameters
     *
     * q = (1/2)Ω_m + Ω_r - Ω_Λ
     *
     * @param Omega_matter Ω_m
     * @param Omega_radiation Ω_r
     * @param Omega_lambda Ω_Λ
     * @return Deceleration parameter q
     */
    static double decelerationFromOmegas(double Omega_matter,
                                        double Omega_radiation,
                                        double Omega_lambda) {
        return 0.5 * Omega_matter + Omega_radiation - Omega_lambda;
    }

    /**
     * @brief Current deceleration parameter
     *
     * q₀ ≈ -0.55 (accelerating!)
     *
     * Universe is accelerating due to dark energy
     */
    static double current() {
        return -0.55;  // Negative → acceleration
    }

    /**
     * @brief Check if universe is accelerating
     */
    static bool isAccelerating(double q) {
        return q < 0.0;
    }

    /**
     * @brief Transition from deceleration to acceleration
     *
     * Occurs when q = 0
     * For ΛCDM: happened at z ≈ 0.7 (t ≈ 7 Gyr ago)
     */
    static double transitionRedshift() {
        return 0.7;  // z_transition
    }
};

/**
 * @class HubbleParameterEvolution
 * @brief Evolution of Hubble parameter H(z)
 */
class HubbleParameterEvolution {
public:
    /**
     * @brief Hubble parameter as function of redshift
     *
     * H(z) = H₀ √[Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_k(1+z)² + Ω_Λ]
     *
     * @param redshift z
     * @param H0 Hubble constant today
     * @param Omega_matter Ω_m
     * @param Omega_radiation Ω_r
     * @param Omega_lambda Ω_Λ
     * @param Omega_k Ω_k (curvature)
     * @return H(z)
     */
    static double hubbleAtRedshift(double redshift, double H0,
                                  double Omega_matter,
                                  double Omega_radiation,
                                  double Omega_lambda,
                                  double Omega_k = 0.0) {
        double z_plus_1 = 1.0 + redshift;

        double matter_term = Omega_matter * std::pow(z_plus_1, 3.0);
        double rad_term = Omega_radiation * std::pow(z_plus_1, 4.0);
        double curv_term = Omega_k * std::pow(z_plus_1, 2.0);
        double lambda_term = Omega_lambda;

        double sum = matter_term + rad_term + curv_term + lambda_term;
        return H0 * std::sqrt(sum);
    }

    /**
     * @brief Hubble parameter in matter era (approximate)
     *
     * H(a) = H₀ a^(-3/2) for matter-dominated
     */
    static double hubbleMatterEra(double scale_factor, double H0) {
        return H0 * std::pow(scale_factor, -1.5);
    }

    /**
     * @brief Hubble parameter in radiation era
     *
     * H(a) = H₀ a^(-2) for radiation-dominated
     */
    static double hubbleRadiationEra(double scale_factor, double H0) {
        return H0 * std::pow(scale_factor, -2.0);
    }

    /**
     * @brief Hubble parameter with dark energy
     *
     * H(a) → H₀√Ω_Λ as a → ∞ (de Sitter expansion)
     */
    static double hubbleFuture(double H0, double Omega_lambda) {
        return H0 * std::sqrt(Omega_lambda);
    }
};

/**
 * @class CosmicDynamics
 * @brief Solve cosmological evolution
 */
class CosmicDynamics {
public:
    /**
     * @brief Scale factor evolution (numerical solution)
     *
     * Solve: ȧ = aH(a)
     *
     * Using Euler method for simplicity
     */
    static std::vector<double> evolveScaleFactor(double a_initial,
                                                double t_initial,
                                                double t_final,
                                                int num_steps,
                                                double H0,
                                                double Omega_matter,
                                                double Omega_lambda) {
        std::vector<double> a_values;
        a_values.reserve(num_steps);

        double dt = (t_final - t_initial) / num_steps;
        double a = a_initial;

        for (int i = 0; i < num_steps; ++i) {
            a_values.push_back(a);

            // H(a) for flat ΛCDM
            double sqrt_term = Omega_matter * std::pow(a, -3.0) + Omega_lambda;
            double H = H0 * std::sqrt(sqrt_term);

            // da/dt = aH
            double a_dot = a * H;

            // Euler step
            a += a_dot * dt;
        }

        return a_values;
    }

    /**
     * @brief Age-redshift relation (numerical integral)
     *
     * t(z) = ∫₀^z dz' / [(1+z')H(z')]
     */
    static double ageAtRedshift(double redshift, double H0,
                               double Omega_matter, double Omega_lambda) {
        // Numerical integration (trapezoidal)
        int N = 1000;
        double dz = redshift / N;
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double z = i * dz;
            double H_z = HubbleParameterEvolution::hubbleAtRedshift(
                z, H0, Omega_matter, 0.0, Omega_lambda, 0.0
            );

            sum += 1.0 / ((1.0 + z) * H_z);
        }

        return sum * dz;  // Time in units of 1/H0
    }
};

} // namespace physics::advanced::cosmology

#endif // PHYSICS_ADVANCED_COSMOLOGY_FRIEDMANN_EQUATIONS_HPP
