/**
 * @file field_tensor_demo.cpp
 * @brief Demonstration of OOP-based advanced physics implementation
 *
 * This shows how tensor mathematics can be implemented beyond header-only approach.
 *
 * Compile with:
 *   g++ -std=c++17 -I.. -I/usr/include/eigen3 field_tensor_demo.cpp -o field_tensor_demo
 *
 * Or if Eigen is installed elsewhere:
 *   g++ -std=c++17 -I.. -I/path/to/eigen3 field_tensor_demo.cpp -o field_tensor_demo
 */

#include "field_tensor_example.hpp"
#include <iostream>
#include <iomanip>

using namespace physics::examples;

void demonstrateBasicFields() {
    std::cout << "=== Example 1: Basic Electric and Magnetic Fields ===\n\n";

    // Create a simple EM field
    Eigen::Vector3d E(1e5, 0, 0);      // 100 kV/m in x-direction
    Eigen::Vector3d B(0, 1e-3, 0);     // 1 mT in y-direction

    ElectromagneticFieldTensor F(E, B);

    std::cout << F << "\n";

    std::cout << "Energy density: " << F.energyDensity() << " J/m³\n";
    std::cout << "Poynting vector: " << F.poyntingVector().transpose() << " W/m²\n";
    std::cout << "Is plane wave? " << (F.isPlaneWave() ? "Yes" : "No") << "\n\n";
}

void demonstratePlaneWave() {
    std::cout << "=== Example 2: Plane Electromagnetic Wave ===\n\n";

    // Plane wave propagating in +z direction, polarized in x-direction
    Eigen::Vector3d direction(0, 0, 1);
    Eigen::Vector3d polarization(1, 0, 0);
    double E0 = 100.0;  // 100 V/m amplitude

    ElectromagneticFieldTensor wave = ElectromagneticFieldTensor::planeWave(
        E0, direction, polarization
    );

    std::cout << wave << "\n";

    std::cout << "Energy density: " << wave.energyDensity() << " J/m³\n";
    std::cout << "Is plane wave? " << (wave.isPlaneWave() ? "Yes" : "No") << "\n\n";

    // Verify E/B = c relation for plane waves
    Eigen::Vector3d E = wave.electricField();
    Eigen::Vector3d B = wave.magneticField();
    double ratio = E.norm() / (299792458.0 * B.norm());
    std::cout << "E/(cB) ratio: " << ratio << " (should be 1.0)\n\n";
}

void demonstrateLorentzBoost() {
    std::cout << "=== Example 3: Lorentz Transformation of Fields ===\n\n";

    // Start with electric field only in lab frame
    Eigen::Vector3d E(1e6, 0, 0);      // 1 MV/m
    Eigen::Vector3d B(0, 0, 0);        // No magnetic field

    ElectromagneticFieldTensor F_lab(E, B);

    std::cout << "Fields in lab frame:\n";
    std::cout << F_lab << "\n";

    // Boost to frame moving at 0.6c in y-direction
    double c = 299792458.0;
    Eigen::Vector3d velocity(0, 0.6 * c, 0);

    ElectromagneticFieldTensor F_moving = F_lab.boost(velocity);

    std::cout << "Fields in moving frame (v = 0.6c in y-direction):\n";
    std::cout << F_moving << "\n";

    // Verify invariants are preserved
    std::cout << "Lab frame invariants:\n";
    std::cout << "  I₁ = " << F_lab.firstInvariant() << "\n";
    std::cout << "  I₂ = " << F_lab.secondInvariant() << "\n\n";

    std::cout << "Moving frame invariants:\n";
    std::cout << "  I₁ = " << F_moving.firstInvariant() << "\n";
    std::cout << "  I₂ = " << F_moving.secondInvariant() << "\n";
    std::cout << "(Invariants should be equal)\n\n";

    // Note: Pure electric field in one frame appears as E + B in moving frame
    Eigen::Vector3d E_moving = F_moving.electricField();
    Eigen::Vector3d B_moving = F_moving.magneticField();
    std::cout << "Magnetic field appears in moving frame: B = "
              << B_moving.transpose() << " T\n\n";
}

void demonstrateStressEnergyTensor() {
    std::cout << "=== Example 4: Stress-Energy Tensor ===\n\n";

    // Plane wave
    Eigen::Vector3d direction(0, 0, 1);
    Eigen::Vector3d polarization(1, 0, 0);
    double E0 = 377.0;  // V/m (gives nice numbers: u = c²B²/μ₀ ≈ 1 J/m³)

    ElectromagneticFieldTensor F = ElectromagneticFieldTensor::planeWave(
        E0, direction, polarization
    );

    StressEnergyTensor T(F);

    std::cout << T << "\n";

    // For EM field, verify traceless property
    std::cout << "Traceless check: |T^μ_μ| = " << std::abs(T.trace())
              << " (should be ~0)\n\n";

    // Energy flux in z-direction should equal c × energy density
    Eigen::Vector3d momentum = T.momentumDensity();
    double energy_density = T.energyDensity();
    std::cout << "Energy flux (Poynting): " << momentum(2) * 299792458.0 << " W/m²\n";
    std::cout << "c × u: " << 299792458.0 * energy_density << " W/m²\n";
    std::cout << "(These should be equal for plane wave)\n\n";
}

void demonstrateFieldInvariants() {
    std::cout << "=== Example 5: Field Classification by Invariants ===\n\n";

    struct FieldCase {
        std::string name;
        Eigen::Vector3d E;
        Eigen::Vector3d B;
    };

    double c = 299792458.0;

    std::vector<FieldCase> cases = {
        {"Pure electric field", Eigen::Vector3d(1e5, 0, 0), Eigen::Vector3d(0, 0, 0)},
        {"Pure magnetic field", Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1e-3, 0)},
        {"Plane wave (E ⊥ B, E = cB)", Eigen::Vector3d(300, 0, 0), Eigen::Vector3d(0, 0, 1e-6)},
        {"E ∥ B (both parallel)", Eigen::Vector3d(1e5, 0, 0), Eigen::Vector3d(1e-3, 0, 0)},
    };

    for (const auto& test : cases) {
        ElectromagneticFieldTensor F(test.E, test.B);

        std::cout << "Case: " << test.name << "\n";
        std::cout << "  E = " << test.E.transpose() << " V/m\n";
        std::cout << "  B = " << test.B.transpose() << " T\n";
        std::cout << "  I₁ = " << std::scientific << F.firstInvariant()
                  << " (B² - E²/c²)\n";
        std::cout << "  I₂ = " << std::scientific << F.secondInvariant()
                  << " (E·B/c)\n";

        if (std::abs(F.firstInvariant()) < 1e-10 &&
            std::abs(F.secondInvariant()) < 1e-10) {
            std::cout << "  → Pure radiation (both invariants = 0)\n";
        } else if (F.firstInvariant() < 0) {
            std::cout << "  → Electric field dominant\n";
        } else {
            std::cout << "  → Magnetic field dominant\n";
        }

        std::cout << std::fixed << "\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Electromagnetic Field Tensor - OOP Implementation    ║\n";
    std::cout << "║  Demonstrates advanced physics beyond header-only     ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n\n";

    try {
        demonstrateBasicFields();
        std::cout << std::string(60, '-') << "\n\n";

        demonstratePlaneWave();
        std::cout << std::string(60, '-') << "\n\n";

        demonstrateLorentzBoost();
        std::cout << std::string(60, '-') << "\n\n";

        demonstrateStressEnergyTensor();
        std::cout << std::string(60, '-') << "\n\n";

        demonstrateFieldInvariants();

        std::cout << "\n✓ All demonstrations completed successfully!\n\n";

        std::cout << "This example shows how OOP enables:\n";
        std::cout << "  • Stateful objects (tensor storage)\n";
        std::cout << "  • Matrix operations (Lorentz transformations)\n";
        std::cout << "  • Complex calculations (invariants, stress-energy)\n";
        std::cout << "  • Operator overloading (<<, mathematical ops)\n";
        std::cout << "  • Encapsulation (private tensor data)\n\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
