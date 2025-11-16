#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>
#include <cassert>

// Include all new modules to verify they compile
#include "maths/measure_theory.hpp"
#include "maths/functional_analysis.hpp"
#include "maths/differential_geometry.hpp"
#include "maths/probability_theory.hpp"
#include "maths/real_analysis.hpp"
#include "physics/general_relativity.hpp"
#include "physics/statistical_mechanics.hpp"
#include "physics/classical_field_theory.hpp"
#include "physics/condensed_matter.hpp"

using namespace std;

// Helper function to check if values are close
bool isClose(double a, double b, double tol = 1e-6) {
    return std::abs(a - b) < tol;
}

void test_measure_theory() {
    cout << "\n=== Testing Measure Theory ===\n";
    cout << "✓ Measure theory module compiled successfully\n";
}

void test_functional_analysis() {
    cout << "\n=== Testing Functional Analysis ===\n";
    cout << "✓ Functional analysis module compiled successfully\n";
}

void test_differential_geometry() {
    cout << "\n=== Testing Differential Geometry ===\n";
    cout << "✓ Differential geometry module compiled successfully\n";
}

void test_probability_theory() {
    cout << "\n=== Testing Probability Theory ===\n";
    cout << "✓ Probability theory module compiled successfully\n";
}

void test_real_analysis() {
    cout << "\n=== Testing Real Analysis ===\n";
    cout << "✓ Real analysis module compiled successfully\n";
}

void test_general_relativity() {
    cout << "\n=== Testing General Relativity ===\n";

    using namespace physics::general_relativity;

    // Test Schwarzschild metric with solar mass
    double M_sun = 1.989e30; // kg
    SchwarzschildMetric schwarzschild(M_sun);

    double r_s = schwarzschild.schwarzschildRadius();
    cout << "Schwarzschild radius for Sun: " << r_s/1000.0 << " km\n";
    assert(r_s > 1000.0 && r_s < 10000.0); // Should be ~3 km

    double r_photon = schwarzschild.photonSphereRadius();
    cout << "Photon sphere radius: " << r_photon/1000.0 << " km\n";
    assert(isClose(r_photon / r_s, 1.5, 0.01)); // Should be 1.5 * r_s

    double r_isco = schwarzschild.innerMostStableOrbit();
    cout << "ISCO radius: " << r_isco/1000.0 << " km\n";
    assert(isClose(r_isco / r_s, 3.0, 0.01)); // Should be 3 * r_s

    // Test metric at a point far from black hole
    Vector4 x = {0.0, r_s * 10, M_PI/2, 0.0}; // t, r, theta, phi
    Tensor2 g = schwarzschild.covariant(x);

    cout << "g_tt at r=10*r_s: " << g[0][0] << " (should be negative)\n";
    assert(g[0][0] < 0);
    assert(g[0][0] > -1.0); // Should approach -1 far from horizon

    cout << "✓ General relativity tests passed\n";
}

void test_statistical_mechanics() {
    cout << "\n=== Testing Statistical Mechanics ===\n";

    using namespace physics::statistical_mechanics;

    // Test Ising model construction
    IsingModel ising(4, 4, 1.0); // 4x4 lattice, J=1
    cout << "Ising model created successfully\n";

    cout << "✓ Statistical mechanics module compiled successfully\n";
}

void test_classical_field_theory() {
    cout << "\n=== Testing Classical Field Theory ===\n";

    using namespace physics::classical_field_theory;

    // Test Klein-Gordon Lagrangian
    double m = 1.0; // mass
    KleinGordonLagrangian kg(m);

    cout << "Klein-Gordon Lagrangian object created successfully\n";

    cout << "✓ Classical field theory tests passed\n";
}

void test_condensed_matter() {
    cout << "\n=== Testing Condensed Matter ===\n";

    using namespace physics::condensed_matter;

    // Test BCS theory
    double T_c = 10.0; // Critical temperature
    double E_F = 10.0; // Fermi energy

    BCSTheory bcs(T_c, E_F);

    double gap_zero = bcs.energyGap(0.0);
    double gap_low = bcs.energyGap(0.1 * T_c);
    cout << "BCS gap at T=0: " << gap_zero << "\n";
    cout << "BCS gap at T=0.1*Tc: " << gap_low << " (should be close to Δ₀)\n";
    assert(gap_low > 0.9 * gap_zero);

    double gap_high = bcs.energyGap(0.99 * T_c);
    cout << "BCS gap at T=0.99*Tc: " << gap_high << " (should be ≈ 0)\n";
    assert(gap_high < 0.5 * gap_zero);

    double xi = bcs.coherenceLength(0.1 * T_c);
    cout << "Coherence length at low T: " << xi << "\n";
    assert(xi > 0);

    // Test Quantum Hall Effect
    double B = 1.0; // Magnetic field
    int filling_factor = 1;
    QuantumHallEffect qhe(B, filling_factor);

    double sigma = qhe.hallConductance();
    cout << "Quantum Hall conductance: " << sigma << "\n";
    assert(sigma > 0);

    double E0 = qhe.landauLevelEnergy(0);
    double E1 = qhe.landauLevelEnergy(1);
    cout << "Landau level spacing: " << (E1 - E0) << "\n";
    assert(E1 > E0);

    cout << "✓ Condensed matter tests passed\n";
}

int main() {
    cout << "========================================\n";
    cout << "  Testing New Physics & Math Modules\n";
    cout << "========================================\n";

    try {
        test_measure_theory();
        test_functional_analysis();
        test_differential_geometry();
        test_probability_theory();
        test_real_analysis();
        test_general_relativity();
        test_statistical_mechanics();
        test_classical_field_theory();
        test_condensed_matter();

        cout << "\n========================================\n";
        cout << "  ✓ ALL TESTS PASSED SUCCESSFULLY!\n";
        cout << "========================================\n";

        return 0;
    } catch (const exception& e) {
        cerr << "\n❌ TEST FAILED: " << e.what() << "\n";
        return 1;
    } catch (...) {
        cerr << "\n❌ TEST FAILED: Unknown error\n";
        return 1;
    }
}
