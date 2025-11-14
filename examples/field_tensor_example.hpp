#ifndef EXAMPLE_FIELD_TENSOR_HPP
#define EXAMPLE_FIELD_TENSOR_HPP

/**
 * @file field_tensor_example.hpp
 * @brief Example implementation of electromagnetic field tensor using OOP
 *
 * This demonstrates how advanced topics can be implemented beyond the header-only approach.
 *
 * Compile with:
 * g++ -std=c++17 -I/usr/include/eigen3 field_tensor_example.cpp -o field_tensor
 */

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace physics::examples {

/**
 * @class ElectromagneticFieldTensor
 * @brief Represents the electromagnetic field tensor F^μν
 *
 * The field tensor is an antisymmetric rank-2 tensor containing
 * the electric and magnetic fields:
 *
 *           ⎛  0   -Ex/c -Ey/c -Ez/c ⎞
 * F^μν  =   ⎜ Ex/c   0    -Bz    By  ⎟
 *           ⎜ Ey/c   Bz    0    -Bx  ⎟
 *           ⎝ Ez/c  -By    Bx    0   ⎠
 *
 * This encodes Maxwell's equations in covariant form.
 */
class ElectromagneticFieldTensor {
private:
    Eigen::Matrix4d F_;           // Contravariant tensor F^μν
    static constexpr double c_ = 299792458.0;  // Speed of light (m/s)

public:
    /**
     * @brief Construct from electric and magnetic field vectors
     * @param E Electric field vector (V/m)
     * @param B Magnetic field vector (T)
     */
    ElectromagneticFieldTensor(const Eigen::Vector3d& E, const Eigen::Vector3d& B) {
        F_.setZero();

        // F^0i = E_i/c (time-space components)
        F_(0, 1) =  E(0) / c_;
        F_(0, 2) =  E(1) / c_;
        F_(0, 3) =  E(2) / c_;

        // F^ij = -ε_ijk B_k (space-space components)
        F_(1, 2) = -B(2);
        F_(1, 3) =  B(1);
        F_(2, 3) = -B(0);

        // Antisymmetric: F^μν = -F^νμ
        F_(1, 0) = -F_(0, 1);
        F_(2, 0) = -F_(0, 2);
        F_(3, 0) = -F_(0, 3);
        F_(2, 1) = -F_(1, 2);
        F_(3, 1) = -F_(1, 3);
        F_(3, 2) = -F_(2, 3);
    }

    /**
     * @brief Get the raw tensor components
     */
    const Eigen::Matrix4d& tensor() const { return F_; }

    /**
     * @brief Extract electric field from tensor
     * @return Electric field vector (V/m)
     */
    Eigen::Vector3d electricField() const {
        return Eigen::Vector3d(
            F_(0, 1) * c_,
            F_(0, 2) * c_,
            F_(0, 3) * c_
        );
    }

    /**
     * @brief Extract magnetic field from tensor
     * @return Magnetic field vector (T)
     */
    Eigen::Vector3d magneticField() const {
        return Eigen::Vector3d(
            F_(2, 3),
           -F_(1, 3),
            F_(1, 2)
        );
    }

    /**
     * @brief Calculate first Lorentz invariant
     *
     * I₁ = (1/2)F_μν F^μν = B² - E²/c²
     *
     * This is invariant under Lorentz transformations
     *
     * @return First invariant (dimensionless when normalized)
     */
    double firstInvariant() const {
        Eigen::Vector3d E = electricField();
        Eigen::Vector3d B = magneticField();
        return B.squaredNorm() - E.squaredNorm() / (c_ * c_);
    }

    /**
     * @brief Calculate second Lorentz invariant
     *
     * I₂ = (1/4)F_μν *F^μν = E·B/c
     *
     * where *F is the dual tensor
     *
     * @return Second invariant
     */
    double secondInvariant() const {
        Eigen::Vector3d E = electricField();
        Eigen::Vector3d B = magneticField();
        return E.dot(B) / c_;
    }

    /**
     * @brief Apply Lorentz boost to the field tensor
     *
     * Transforms fields from one inertial frame to another moving
     * at velocity v relative to the original frame.
     *
     * E'_∥ = E_∥
     * E'_⊥ = γ(E_⊥ + v × B)
     * B'_∥ = B_∥
     * B'_⊥ = γ(B_⊥ - (v × E)/c²)
     *
     * @param velocity Boost velocity vector (m/s)
     * @return Transformed field tensor in boosted frame
     */
    ElectromagneticFieldTensor boost(const Eigen::Vector3d& velocity) const {
        double v_mag = velocity.norm();

        if (v_mag < 1e-10) {
            return *this;  // No boost
        }

        if (v_mag >= c_) {
            throw std::invalid_argument("Boost velocity must be < c");
        }

        Eigen::Vector3d v_hat = velocity / v_mag;
        double beta = v_mag / c_;
        double gamma = 1.0 / std::sqrt(1.0 - beta * beta);

        Eigen::Vector3d E = electricField();
        Eigen::Vector3d B = magneticField();

        // Decompose into parallel and perpendicular components
        Eigen::Vector3d E_parallel = E.dot(v_hat) * v_hat;
        Eigen::Vector3d E_perp = E - E_parallel;

        Eigen::Vector3d B_parallel = B.dot(v_hat) * v_hat;
        Eigen::Vector3d B_perp = B - B_parallel;

        // Transform
        Eigen::Vector3d E_prime = E_parallel + gamma * (E_perp + velocity.cross(B));
        Eigen::Vector3d B_prime = B_parallel + gamma * (B_perp - velocity.cross(E) / (c_ * c_));

        return ElectromagneticFieldTensor(E_prime, B_prime);
    }

    /**
     * @brief Calculate Poynting vector S = (E × B)/μ₀
     *
     * Energy flux density of electromagnetic field
     *
     * @param mu0 Permeability of free space (default: 4π×10⁻⁷ H/m)
     * @return Poynting vector (W/m²)
     */
    Eigen::Vector3d poyntingVector(double mu0 = 4.0 * M_PI * 1e-7) const {
        Eigen::Vector3d E = electricField();
        Eigen::Vector3d B = magneticField();
        return E.cross(B) / mu0;
    }

    /**
     * @brief Calculate electromagnetic energy density
     *
     * u = (ε₀E²/2 + B²/(2μ₀))
     *
     * @param epsilon0 Permittivity of free space (default: 8.854×10⁻¹² F/m)
     * @param mu0 Permeability of free space (default: 4π×10⁻⁷ H/m)
     * @return Energy density (J/m³)
     */
    double energyDensity(double epsilon0 = 8.854187817e-12,
                         double mu0 = 4.0 * M_PI * 1e-7) const {
        Eigen::Vector3d E = electricField();
        Eigen::Vector3d B = magneticField();
        return 0.5 * epsilon0 * E.squaredNorm() + 0.5 * B.squaredNorm() / mu0;
    }

    /**
     * @brief Check if configuration is a pure electromagnetic wave
     *
     * For a plane wave: E ⊥ B, |E| = c|B|
     *
     * @param tolerance Numerical tolerance
     * @return true if fields represent a plane wave
     */
    bool isPlaneWave(double tolerance = 1e-6) const {
        Eigen::Vector3d E = electricField();
        Eigen::Vector3d B = magneticField();

        // Check perpendicularity
        double dot_product = std::abs(E.dot(B));
        if (dot_product > tolerance * E.norm() * B.norm()) {
            return false;
        }

        // Check magnitude relation
        double ratio = E.norm() / (c_ * B.norm());
        return std::abs(ratio - 1.0) < tolerance;
    }

    /**
     * @brief Pretty print the tensor
     */
    friend std::ostream& operator<<(std::ostream& os, const ElectromagneticFieldTensor& F) {
        os << "Electromagnetic Field Tensor F^μν:\n";
        os << std::fixed << std::setprecision(6);
        os << F.F_ << "\n\n";

        Eigen::Vector3d E = F.electricField();
        Eigen::Vector3d B = F.magneticField();

        os << "E = (" << E.transpose() << ") V/m\n";
        os << "B = (" << B.transpose() << ") T\n";
        os << "First invariant: " << F.firstInvariant() << "\n";
        os << "Second invariant: " << F.secondInvariant() << "\n";

        return os;
    }

    /**
     * @brief Create field for a plane wave
     * @param E0 Electric field amplitude (V/m)
     * @param direction Propagation direction (unit vector)
     * @param polarization Polarization direction (unit vector, ⊥ to direction)
     * @return Field tensor for plane wave
     */
    static ElectromagneticFieldTensor planeWave(
        double E0,
        const Eigen::Vector3d& direction,
        const Eigen::Vector3d& polarization) {

        Eigen::Vector3d k_hat = direction.normalized();
        Eigen::Vector3d E_hat = polarization.normalized();

        // Ensure perpendicularity
        if (std::abs(k_hat.dot(E_hat)) > 1e-6) {
            throw std::invalid_argument("Polarization must be perpendicular to propagation");
        }

        Eigen::Vector3d E = E0 * E_hat;
        Eigen::Vector3d B = (E0 / c_) * k_hat.cross(E_hat);

        return ElectromagneticFieldTensor(E, B);
    }
};

/**
 * @class StressEnergyTensor
 * @brief Electromagnetic stress-energy tensor T^μν
 *
 * The stress-energy tensor describes energy and momentum distribution:
 * - T^00: Energy density
 * - T^0i: Momentum density / energy flux
 * - T^ij: Stress (momentum flux)
 */
class StressEnergyTensor {
private:
    Eigen::Matrix4d T_;
    static constexpr double c_ = 299792458.0;
    static constexpr double epsilon0_ = 8.854187817e-12;
    static constexpr double mu0_ = 4.0 * M_PI * 1e-7;

public:
    /**
     * @brief Construct electromagnetic stress-energy tensor from field tensor
     *
     * T^μν = (1/μ₀)[F^μα F^ν_α - (1/4)g^μν F^αβ F_αβ]
     *
     * @param F Electromagnetic field tensor
     */
    StressEnergyTensor(const ElectromagneticFieldTensor& F) {
        Eigen::Vector3d E = F.electricField();
        Eigen::Vector3d B = F.magneticField();

        // Energy density
        double u = 0.5 * epsilon0_ * E.squaredNorm() + 0.5 * B.squaredNorm() / mu0_;
        T_(0, 0) = u;

        // Poynting vector / c² (momentum density)
        Eigen::Vector3d S = E.cross(B) / mu0_;
        for (int i = 0; i < 3; ++i) {
            T_(0, i+1) = S(i) / (c_ * c_);
            T_(i+1, 0) = S(i) / (c_ * c_);
        }

        // Maxwell stress tensor
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double delta_ij = (i == j) ? 1.0 : 0.0;
                T_(i+1, j+1) = epsilon0_ * (E(i) * E(j) - 0.5 * delta_ij * E.squaredNorm())
                             + (1.0/mu0_) * (B(i) * B(j) - 0.5 * delta_ij * B.squaredNorm());
            }
        }
    }

    /**
     * @brief Get energy density
     */
    double energyDensity() const { return T_(0, 0); }

    /**
     * @brief Get momentum density vector
     */
    Eigen::Vector3d momentumDensity() const {
        return Eigen::Vector3d(T_(0, 1), T_(0, 2), T_(0, 3)) * c_ * c_;
    }

    /**
     * @brief Get stress tensor (3×3 spatial part)
     */
    Eigen::Matrix3d stressTensor() const {
        return T_.block<3, 3>(1, 1);
    }

    /**
     * @brief Calculate trace (should be zero for EM field)
     *
     * T^μ_μ = 0 (traceless for pure radiation)
     */
    double trace() const {
        // With Minkowski metric: T^μ_μ = T^00 - T^11 - T^22 - T^33
        return T_(0,0) - T_(1,1) - T_(2,2) - T_(3,3);
    }

    const Eigen::Matrix4d& tensor() const { return T_; }

    friend std::ostream& operator<<(std::ostream& os, const StressEnergyTensor& T) {
        os << "Stress-Energy Tensor T^μν:\n";
        os << std::fixed << std::setprecision(6);
        os << T.T_ << "\n\n";
        os << "Energy density: " << T.energyDensity() << " J/m³\n";
        os << "Momentum density: " << T.momentumDensity().transpose() << " kg/(m²·s)\n";
        os << "Trace: " << T.trace() << " (should be ~0 for EM field)\n";
        return os;
    }
};

} // namespace physics::examples

#endif // EXAMPLE_FIELD_TENSOR_HPP
