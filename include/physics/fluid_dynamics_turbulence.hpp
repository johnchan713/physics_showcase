#ifndef PHYSICS_ADVANCED_FLUID_DYNAMICS_TURBULENCE_HPP
#define PHYSICS_ADVANCED_FLUID_DYNAMICS_TURBULENCE_HPP

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

/**
 * @file turbulence.hpp
 * @brief Turbulence models for Reynolds-Averaged Navier-Stokes (RANS)
 *
 * Implements:
 * - Reynolds decomposition and averaging
 * - Reynolds stress tensor
 * - k-ε turbulence model (two-equation)
 * - Mixing length theory (Prandtl)
 * - Eddy viscosity models
 */

namespace physics::advanced::fluid_dynamics {

/**
 * @class ReynoldsDecomposition
 * @brief Reynolds decomposition: u = ū + u'
 *
 * Separates flow variables into mean and fluctuating components
 */
class ReynoldsDecomposition {
public:
    /**
     * @brief Reynolds averaging operator
     *
     * For velocity: u = ū + u'
     * where ū is mean, u' is fluctuation
     *
     * Properties:
     * - ⟨ū⟩ = ū
     * - ⟨u'⟩ = 0
     * - ⟨u'v'⟩ ≠ 0 (Reynolds stress)
     */
    struct DecomposedField {
        Eigen::Vector3d mean;        // ū
        Eigen::Vector3d fluctuation; // u'
    };

    /**
     * @brief Create decomposed field
     */
    static DecomposedField decompose(const Eigen::Vector3d& instantaneous,
                                     const Eigen::Vector3d& mean) {
        return {mean, instantaneous - mean};
    }

    /**
     * @brief Verify Reynolds averaging properties
     *
     * Check: ⟨ū⟩ = ū and ⟨u'⟩ = 0
     */
    static bool verifyAveraging(
        const std::vector<Eigen::Vector3d>& fluctuations,
        double tolerance = 1e-6) {

        Eigen::Vector3d avg = Eigen::Vector3d::Zero();
        for (const auto& u_prime : fluctuations) {
            avg += u_prime;
        }
        avg /= fluctuations.size();

        return avg.norm() < tolerance;
    }
};

/**
 * @class ReynoldsStressTensor
 * @brief Reynolds stress tensor -ρ⟨u'ᵢu'ⱼ⟩
 *
 * Represents turbulent momentum transport
 */
class ReynoldsStressTensor {
public:
    /**
     * @brief Compute Reynolds stress tensor from fluctuations
     *
     * τᵢⱼᴿ = -ρ⟨u'ᵢu'ⱼ⟩
     *
     * @param density ρ (kg/m³)
     * @param velocity_fluctuations Sample of u' vectors
     * @return 3×3 Reynolds stress tensor (Pa)
     */
    static Eigen::Matrix3d compute(
        double density,
        const std::vector<Eigen::Vector3d>& velocity_fluctuations) {

        Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

        for (const auto& u_prime : velocity_fluctuations) {
            stress += u_prime * u_prime.transpose();
        }

        stress /= velocity_fluctuations.size();
        stress *= -density;

        return stress;
    }

    /**
     * @brief Turbulent kinetic energy
     *
     * k = ½⟨u'ᵢu'ᵢ⟩ = ½(⟨u'²⟩ + ⟨v'²⟩ + ⟨w'²⟩)
     *
     * @param velocity_fluctuations Sample of u' vectors
     * @return Turbulent kinetic energy (m²/s²)
     */
    static double turbulentKineticEnergy(
        const std::vector<Eigen::Vector3d>& velocity_fluctuations) {

        double k = 0.0;
        for (const auto& u_prime : velocity_fluctuations) {
            k += u_prime.squaredNorm();
        }

        return 0.5 * k / velocity_fluctuations.size();
    }

    /**
     * @brief Turbulence intensity
     *
     * I = √(⅔k) / U
     *
     * Typical values:
     * - Low turbulence: I < 1%
     * - Medium turbulence: 1% < I < 5%
     * - High turbulence: I > 5%
     *
     * @param k Turbulent kinetic energy (m²/s²)
     * @param mean_velocity Mean velocity magnitude (m/s)
     * @return Turbulence intensity (dimensionless)
     */
    static double turbulenceIntensity(double k, double mean_velocity) {
        if (mean_velocity <= 0.0) {
            throw std::invalid_argument("Mean velocity must be positive");
        }
        return std::sqrt(2.0 * k / 3.0) / mean_velocity;
    }

    /**
     * @brief Anisotropy tensor
     *
     * bᵢⱼ = ⟨u'ᵢu'ⱼ⟩/(2k) - δᵢⱼ/3
     *
     * Measures deviation from isotropic turbulence
     * For isotropic turbulence: bᵢⱼ = 0
     *
     * @param reynolds_stress Reynolds stress tensor (normalized)
     * @param k Turbulent kinetic energy
     * @return Anisotropy tensor
     */
    static Eigen::Matrix3d anisotropyTensor(const Eigen::Matrix3d& reynolds_stress,
                                            double k) {
        Eigen::Matrix3d b = reynolds_stress / (2.0 * k);
        b(0,0) -= 1.0/3.0;
        b(1,1) -= 1.0/3.0;
        b(2,2) -= 1.0/3.0;
        return b;
    }
};

/**
 * @class BoussinesqHypothesis
 * @brief Boussinesq eddy viscosity hypothesis
 *
 * Relates Reynolds stresses to mean velocity gradient:
 * -ρ⟨u'ᵢu'ⱼ⟩ = μₜ(∂ūᵢ/∂xⱼ + ∂ūⱼ/∂xᵢ) - ⅔ρkδᵢⱼ
 */
class BoussinesqHypothesis {
public:
    /**
     * @brief Compute Reynolds stress from eddy viscosity
     *
     * τᵢⱼᴿ = μₜSᵢⱼ - ⅔ρkδᵢⱼ
     *
     * where Sᵢⱼ = ∂ūᵢ/∂xⱼ + ∂ūⱼ/∂xᵢ is mean strain rate
     *
     * @param eddy_viscosity μₜ (Pa·s)
     * @param mean_strain_rate Sᵢⱼ (1/s)
     * @param density ρ (kg/m³)
     * @param turbulent_ke k (m²/s²)
     * @return Reynolds stress tensor (Pa)
     */
    static Eigen::Matrix3d reynoldsStress(
        double eddy_viscosity,
        const Eigen::Matrix3d& mean_strain_rate,
        double density,
        double turbulent_ke) {

        Eigen::Matrix3d tau_R = eddy_viscosity * mean_strain_rate;

        // Subtract ⅔ρk from diagonal
        double diag_term = (2.0/3.0) * density * turbulent_ke;
        tau_R(0,0) -= diag_term;
        tau_R(1,1) -= diag_term;
        tau_R(2,2) -= diag_term;

        return tau_R;
    }

    /**
     * @brief Mean strain rate tensor
     *
     * Sᵢⱼ = ∂ūᵢ/∂xⱼ + ∂ūⱼ/∂xᵢ
     *
     * @param velocity_gradient ∂ūᵢ/∂xⱼ
     * @return Strain rate tensor (1/s)
     */
    static Eigen::Matrix3d strainRateTensor(const Eigen::Matrix3d& velocity_gradient) {
        return velocity_gradient + velocity_gradient.transpose();
    }

    /**
     * @brief Strain rate magnitude
     *
     * S = √(2SᵢⱼSᵢⱼ)
     *
     * @param strain_rate Sᵢⱼ
     * @return Strain rate magnitude (1/s)
     */
    static double strainRateMagnitude(const Eigen::Matrix3d& strain_rate) {
        return std::sqrt(2.0 * (strain_rate.array() * strain_rate.array()).sum());
    }
};

/**
 * @class MixingLengthModel
 * @brief Prandtl mixing length theory (zero-equation model)
 *
 * Algebraic model: μₜ = ρl²|∂ū/∂y|
 */
class MixingLengthModel {
public:
    /**
     * @brief Eddy viscosity from mixing length
     *
     * μₜ = ρl²|∂ū/∂y|
     *
     * @param density ρ (kg/m³)
     * @param mixing_length l (m)
     * @param velocity_gradient ∂ū/∂y (1/s)
     * @return Eddy viscosity (Pa·s)
     */
    static double eddyViscosity(double density,
                                double mixing_length,
                                double velocity_gradient) {
        return density * mixing_length * mixing_length *
               std::abs(velocity_gradient);
    }

    /**
     * @brief Mixing length for pipe/channel flow
     *
     * l = κy(1 - y/δ)
     *
     * where κ ≈ 0.41 (von Kármán constant)
     *
     * @param y Distance from wall (m)
     * @param delta Boundary layer thickness (m)
     * @param kappa von Kármán constant (default 0.41)
     * @return Mixing length (m)
     */
    static double mixingLengthWall(double y, double delta,
                                   double kappa = 0.41) {
        return kappa * y * (1.0 - y / delta);
    }

    /**
     * @brief Mixing length for free shear layer
     *
     * l = C × δ
     *
     * where δ is shear layer thickness, C ≈ 0.07
     *
     * @param shear_layer_thickness δ (m)
     * @param constant C (default 0.07)
     * @return Mixing length (m)
     */
    static double mixingLengthShear(double shear_layer_thickness,
                                    double constant = 0.07) {
        return constant * shear_layer_thickness;
    }
};

/**
 * @class KEpsilonModel
 * @brief k-ε turbulence model (two-equation model)
 *
 * Transport equations for:
 * - k: Turbulent kinetic energy
 * - ε: Turbulent dissipation rate
 *
 * Standard k-ε model (Launder-Sharma)
 */
class KEpsilonModel {
public:
    /**
     * @brief Model constants (standard k-ε)
     */
    struct Constants {
        double C_mu = 0.09;     // Eddy viscosity constant
        double C_eps1 = 1.44;   // Production constant
        double C_eps2 = 1.92;   // Dissipation constant
        double sigma_k = 1.0;   // k diffusion constant
        double sigma_eps = 1.3; // ε diffusion constant
    };

    /**
     * @brief Eddy viscosity from k-ε
     *
     * μₜ = ρCμk²/ε
     *
     * @param density ρ (kg/m³)
     * @param k Turbulent kinetic energy (m²/s²)
     * @param epsilon Dissipation rate (m²/s³)
     * @param constants Model constants
     * @return Eddy viscosity (Pa·s)
     */
    static double eddyViscosity(double density, double k, double epsilon,
                                const Constants& constants = Constants{}) {
        if (epsilon <= 0.0) {
            throw std::invalid_argument("Epsilon must be positive");
        }
        return density * constants.C_mu * k * k / epsilon;
    }

    /**
     * @brief Turbulent dissipation rate ε
     *
     * ε = μ⟨(∂u'ᵢ/∂xⱼ)(∂u'ᵢ/∂xⱼ)⟩
     *
     * Units: m²/s³
     */
    static double dissipationRate(double eddy_viscosity, double k,
                                  double density,
                                  const Constants& constants = Constants{}) {
        // From μₜ = ρCμk²/ε, solve for ε
        return density * constants.C_mu * k * k / eddy_viscosity;
    }

    /**
     * @brief Production of turbulent kinetic energy
     *
     * Pₖ = μₜSᵢⱼ(∂ūᵢ/∂xⱼ) ≈ μₜS²
     *
     * @param eddy_viscosity μₜ (Pa·s)
     * @param strain_rate_magnitude S (1/s)
     * @return Production rate (kg/(m·s³))
     */
    static double turbulenceProduction(double eddy_viscosity,
                                       double strain_rate_magnitude) {
        return eddy_viscosity * strain_rate_magnitude * strain_rate_magnitude;
    }

    /**
     * @brief Transport equation for k
     *
     * Dk/Dt = ∂/∂xⱼ[(ν + νₜ/σₖ)∂k/∂xⱼ] + Pₖ - ε
     *
     * @param k_current Current k value
     * @param production Pₖ
     * @param epsilon ε
     * @param diffusion Diffusion term
     * @param dt Time step
     * @return New k value
     */
    static double stepKEquation(double k_current,
                                double production,
                                double epsilon,
                                double diffusion,
                                double dt) {
        double dk_dt = diffusion + production - epsilon;
        double k_new = k_current + dt * dk_dt;
        return std::max(k_new, 1e-10);  // Ensure positivity
    }

    /**
     * @brief Transport equation for ε
     *
     * Dε/Dt = ∂/∂xⱼ[(ν + νₜ/σₑ)∂ε/∂xⱼ] + C₁(ε/k)Pₖ - C₂ε²/k
     *
     * @param epsilon_current Current ε value
     * @param k Turbulent kinetic energy
     * @param production Pₖ
     * @param diffusion Diffusion term
     * @param dt Time step
     * @param constants Model constants
     * @return New ε value
     */
    static double stepEpsilonEquation(double epsilon_current,
                                      double k,
                                      double production,
                                      double diffusion,
                                      double dt,
                                      const Constants& constants = Constants{}) {
        if (k <= 0.0) {
            throw std::invalid_argument("k must be positive");
        }

        double source = constants.C_eps1 * (epsilon_current / k) * production;
        double sink = constants.C_eps2 * epsilon_current * epsilon_current / k;
        double deps_dt = diffusion + source - sink;

        double eps_new = epsilon_current + dt * deps_dt;
        return std::max(eps_new, 1e-10);  // Ensure positivity
    }

    /**
     * @brief Turbulent length scale
     *
     * lₜ = k^(3/2)/ε
     *
     * @param k Turbulent kinetic energy (m²/s²)
     * @param epsilon Dissipation rate (m²/s³)
     * @return Length scale (m)
     */
    static double turbulentLengthScale(double k, double epsilon) {
        if (epsilon <= 0.0) {
            throw std::invalid_argument("Epsilon must be positive");
        }
        return std::pow(k, 1.5) / epsilon;
    }

    /**
     * @brief Turbulent time scale
     *
     * τₜ = k/ε
     *
     * @param k Turbulent kinetic energy (m²/s²)
     * @param epsilon Dissipation rate (m²/s³)
     * @return Time scale (s)
     */
    static double turbulentTimeScale(double k, double epsilon) {
        if (epsilon <= 0.0) {
            throw std::invalid_argument("Epsilon must be positive");
        }
        return k / epsilon;
    }

    /**
     * @brief Kolmogorov length scale
     *
     * η = (ν³/ε)^(1/4)
     *
     * Smallest turbulent eddies
     *
     * @param kinematic_viscosity ν (m²/s)
     * @param epsilon ε (m²/s³)
     * @return Kolmogorov length scale (m)
     */
    static double kolmogorovLengthScale(double kinematic_viscosity,
                                        double epsilon) {
        if (epsilon <= 0.0) {
            throw std::invalid_argument("Epsilon must be positive");
        }
        double nu3 = kinematic_viscosity * kinematic_viscosity * kinematic_viscosity;
        return std::pow(nu3 / epsilon, 0.25);
    }

    /**
     * @brief Kolmogorov time scale
     *
     * τη = √(ν/ε)
     *
     * @param kinematic_viscosity ν (m²/s)
     * @param epsilon ε (m²/s³)
     * @return Kolmogorov time scale (s)
     */
    static double kolmogorovTimeScale(double kinematic_viscosity,
                                      double epsilon) {
        if (epsilon <= 0.0) {
            throw std::invalid_argument("Epsilon must be positive");
        }
        return std::sqrt(kinematic_viscosity / epsilon);
    }

    /**
     * @brief Turbulent Reynolds number
     *
     * Reₜ = k²/(νε)
     *
     * @param k Turbulent kinetic energy (m²/s²)
     * @param kinematic_viscosity ν (m²/s)
     * @param epsilon ε (m²/s³)
     * @return Turbulent Reynolds number
     */
    static double turbulentReynolds(double k,
                                    double kinematic_viscosity,
                                    double epsilon) {
        if (epsilon <= 0.0) {
            throw std::invalid_argument("Epsilon must be positive");
        }
        return k * k / (kinematic_viscosity * epsilon);
    }
};

/**
 * @class WallFunctions
 * @brief Wall functions for near-wall turbulence modeling
 */
class WallFunctions {
public:
    /**
     * @brief Dimensionless wall distance
     *
     * y⁺ = yuτ/ν = y√(τw/ρ)/ν
     *
     * @param y Distance from wall (m)
     * @param friction_velocity uτ (m/s)
     * @param kinematic_viscosity ν (m²/s)
     * @return y⁺
     */
    static double yPlus(double y, double friction_velocity,
                        double kinematic_viscosity) {
        return y * friction_velocity / kinematic_viscosity;
    }

    /**
     * @brief Wall shear stress from k-ε
     *
     * τw = ρCμ^(1/4)k^(1/2)uP
     *
     * where uP is velocity at point P
     *
     * @param density ρ (kg/m³)
     * @param k Turbulent kinetic energy at P (m²/s²)
     * @param velocity_P Velocity magnitude at P (m/s)
     * @param C_mu Model constant (default 0.09)
     * @return Wall shear stress (Pa)
     */
    static double wallShearStress(double density, double k,
                                  double velocity_P,
                                  double C_mu = 0.09) {
        return density * std::pow(C_mu, 0.25) * std::sqrt(k) * velocity_P;
    }

    /**
     * @brief Check if point is in viscous sublayer
     *
     * Viscous sublayer: y⁺ < 5
     * Buffer layer: 5 < y⁺ < 30
     * Log layer: 30 < y⁺ < 500
     *
     * @param y_plus y⁺
     * @return true if in viscous sublayer
     */
    static bool isViscousSublayer(double y_plus) {
        return y_plus < 5.0;
    }

    /**
     * @brief Check if point is in log layer
     */
    static bool isLogLayer(double y_plus) {
        return y_plus > 30.0 && y_plus < 500.0;
    }
};

} // namespace physics::advanced::fluid_dynamics

#endif // PHYSICS_ADVANCED_FLUID_DYNAMICS_TURBULENCE_HPP
