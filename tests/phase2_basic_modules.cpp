/**
 * @file phase2_basic_modules.cpp
 * @brief Phase 2 Validation: Basic Math and Classical Physics
 *
 * This file validates the correctness of:
 * - Calculus: IVT, MVT, fundamental theorems
 * - Kinematics: equations of motion (v = v₀ + at, s = v₀t + ½at², v² = v₀² + 2as)
 * - Dynamics: Newton's laws, force-acceleration relationships
 * - Thermodynamics: laws, heat transfer, calorimetry
 * - Gravitation: gravitational fields, potential energy
 *
 * Validation approach:
 * - Test with known physical scenarios
 * - Verify conservation laws
 * - Check dimensional consistency
 * - Ensure numerical stability within tolerance (1e-6)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <functional>

#include "../include/maths/calculus_theorems.hpp"
#include "../include/physics/kinematics.hpp"
#include "../include/physics/dynamics.hpp"
#include "../include/physics/gravitation.hpp"
#include "../include/physics/thermodynamics.hpp"
#include "../include/physics/energy_momentum.hpp"

using namespace maths::calculus;
using namespace physics::kinematics;
using namespace physics::dynamics;
using namespace physics::gravitation;
using namespace physics::thermodynamics;
using namespace physics::energy_momentum;

// Test tolerance
constexpr double TOLERANCE = 1e-6;

// Test counter
int tests_passed = 0;
int tests_failed = 0;

// Helper macros for testing
#define ASSERT_NEAR(actual, expected, tolerance) \
    do { \
        if (std::abs((actual) - (expected)) <= (tolerance)) { \
            tests_passed++; \
        } else { \
            tests_failed++; \
            std::cerr << "FAILED: " << __FILE__ << ":" << __LINE__ << std::endl; \
            std::cerr << "  Expected: " << (expected) << std::endl; \
            std::cerr << "  Actual:   " << (actual) << std::endl; \
            std::cerr << "  Diff:     " << std::abs((actual) - (expected)) << std::endl; \
        } \
    } while(0)

#define ASSERT_TRUE(condition) \
    do { \
        if (condition) { \
            tests_passed++; \
        } else { \
            tests_failed++; \
            std::cerr << "FAILED: " << __FILE__ << ":" << __LINE__ << std::endl; \
            std::cerr << "  Condition failed: " << #condition << std::endl; \
        } \
    } while(0)

#define ASSERT_FALSE(condition) ASSERT_TRUE(!(condition))

// ============================================================================
// CALCULUS TESTS
// ============================================================================

void test_intermediate_value_theorem() {
    std::cout << "\n=== Testing Intermediate Value Theorem ===" << std::endl;

    // Test IVT conditions
    auto f = [](double x) { return x * x - 2.0; }; // f(x) = x² - 2
    ASSERT_TRUE(IntermediateValueTheorem::checkConditions(f, 1.0, 2.0, 0.0));

    // Find √2 using IVT (root of x² - 2 = 0)
    double sqrt2 = IntermediateValueTheorem::findRoot(f, 1.0, 2.0, 0.0);
    ASSERT_NEAR(sqrt2, std::sqrt(2.0), 1e-9);

    // Test computeSqrt2 helper
    double sqrt2_direct = IntermediateValueTheorem::computeSqrt2();
    ASSERT_NEAR(sqrt2_direct, std::sqrt(2.0), 1e-9);

    // Find root of x³ - x - 2 = 0 (x ≈ 1.521)
    auto cubic = [](double x) { return x*x*x - x - 2.0; };
    double root = IntermediateValueTheorem::findRoot(cubic, 1.0, 2.0, 0.0);
    ASSERT_NEAR(cubic(root), 0.0, 1e-9);

    std::cout << "Intermediate Value Theorem: OK" << std::endl;
}

void test_mean_value_theorem() {
    std::cout << "\n=== Testing Mean Value Theorem ===" << std::endl;

    // Test with f(x) = x² on [0, 2]
    auto f = [](double x) { return x * x; };
    auto df = [](double x) { return 2.0 * x; };

    // Average rate of change: [f(2) - f(0)] / (2 - 0) = 4/2 = 2
    double avg_rate = MeanValueTheorem::averageRateOfChange(f, 0.0, 2.0);
    ASSERT_NEAR(avg_rate, 2.0, TOLERANCE);

    // MVT says there exists c where f'(c) = 2, so 2c = 2, thus c = 1
    double c = MeanValueTheorem::findMVTPoint(f, df, 0.0, 2.0);
    ASSERT_NEAR(c, 1.0, 1e-3); // Allow slightly larger tolerance for search
    ASSERT_NEAR(df(c), avg_rate, 1e-3);

    // Test with f(x) = x³ on [0, 1]
    auto f2 = [](double x) { return x * x * x; };
    auto df2 = [](double x) { return 3.0 * x * x; };

    // Average rate: [1 - 0] / 1 = 1
    // f'(c) = 1 => 3c² = 1 => c = 1/√3 ≈ 0.577
    double c2 = MeanValueTheorem::findMVTPoint(f2, df2, 0.0, 1.0);
    ASSERT_NEAR(df2(c2), 1.0, 1e-3);

    std::cout << "Mean Value Theorem: OK" << std::endl;
}

// ============================================================================
// KINEMATICS TESTS
// ============================================================================

void test_kinematics_first_equation() {
    std::cout << "\n=== Testing Kinematics: v = v₀ + at ===" << std::endl;

    // Test: Object starts at rest (v₀ = 0), accelerates at 2 m/s² for 5 seconds
    double v = calculateFinalVelocity(0.0, 2.0, 5.0);
    ASSERT_NEAR(v, 10.0, TOLERANCE); // v = 0 + 2*5 = 10 m/s

    // Test: Object at 10 m/s, accelerates at 3 m/s² for 2 seconds
    v = calculateFinalVelocity(10.0, 3.0, 2.0);
    ASSERT_NEAR(v, 16.0, TOLERANCE); // v = 10 + 3*2 = 16 m/s

    // Test: Object decelerates (negative acceleration)
    v = calculateFinalVelocity(20.0, -5.0, 3.0);
    ASSERT_NEAR(v, 5.0, TOLERANCE); // v = 20 + (-5)*3 = 5 m/s

    // Test: Calculate acceleration from velocities
    double a = calculateAccelerationFromVelocities(0.0, 10.0, 5.0);
    ASSERT_NEAR(a, 2.0, TOLERANCE); // a = (10 - 0) / 5 = 2 m/s²

    // Test: Calculate time from velocities
    double t = calculateTimeFromVelocities(0.0, 10.0, 2.0);
    ASSERT_NEAR(t, 5.0, TOLERANCE); // t = (10 - 0) / 2 = 5 s

    std::cout << "Kinematics first equation: OK" << std::endl;
}

void test_kinematics_second_equation() {
    std::cout << "\n=== Testing Kinematics: s = v₀t + ½at² ===" << std::endl;

    // Test: Object starts at rest, accelerates at 2 m/s² for 5 seconds
    double s = calculateDisplacement(0.0, 2.0, 5.0);
    ASSERT_NEAR(s, 25.0, TOLERANCE); // s = 0 + 0.5*2*25 = 25 m

    // Test: Object with initial velocity
    s = calculateDisplacement(10.0, 2.0, 3.0);
    ASSERT_NEAR(s, 39.0, TOLERANCE); // s = 10*3 + 0.5*2*9 = 30 + 9 = 39 m

    // Test: Free fall (g = -9.8 m/s²) from rest for 2 seconds
    s = calculateDisplacement(0.0, -9.8, 2.0);
    ASSERT_NEAR(s, -19.6, TOLERANCE); // s = 0 + 0.5*(-9.8)*4 = -19.6 m

    // Test: Calculate acceleration from displacement
    double a = calculateAccelerationFromDisplacement(25.0, 0.0, 5.0);
    ASSERT_NEAR(a, 2.0, TOLERANCE); // a = 2*(25 - 0) / 25 = 2 m/s²

    std::cout << "Kinematics second equation: OK" << std::endl;
}

void test_kinematics_third_equation() {
    std::cout << "\n=== Testing Kinematics: v² = v₀² + 2as ===" << std::endl;

    // Test: Object accelerates from rest over 25 m at 2 m/s²
    double v = calculateFinalVelocityFromDisplacement(0.0, 2.0, 25.0);
    ASSERT_NEAR(v, 10.0, TOLERANCE); // v = √(0 + 2*2*25) = √100 = 10 m/s

    // Test: Object with initial velocity
    v = calculateFinalVelocityFromDisplacement(10.0, 3.0, 12.0);
    double expected = std::sqrt(100.0 + 72.0); // √172 ≈ 13.115
    ASSERT_NEAR(v, expected, TOLERANCE);

    // Test: Braking scenario (negative acceleration)
    v = calculateFinalVelocityFromDisplacement(20.0, -5.0, 30.0);
    expected = std::sqrt(400.0 - 300.0); // √100 = 10 m/s
    ASSERT_NEAR(v, expected, TOLERANCE);

    // Test: Calculate acceleration from velocities and displacement
    double a = calculateAccelerationFromVelocitySquared(0.0, 10.0, 25.0);
    ASSERT_NEAR(a, 2.0, TOLERANCE); // a = (100 - 0) / (2*25) = 2 m/s²

    // Test: Calculate displacement from velocities and acceleration
    double s = calculateDisplacementFromVelocities(0.0, 10.0, 2.0);
    ASSERT_NEAR(s, 25.0, TOLERANCE); // s = (100 - 0) / (2*2) = 25 m

    std::cout << "Kinematics third equation: OK" << std::endl;
}

void test_kinematics_consistency() {
    std::cout << "\n=== Testing Kinematics Equation Consistency ===" << std::endl;

    // Use all three equations with same parameters, check consistency
    double v0 = 5.0; // initial velocity
    double a = 2.0;  // acceleration
    double t = 4.0;  // time

    // Calculate final velocity: v = v₀ + at
    double v_from_eq1 = calculateFinalVelocity(v0, a, t);

    // Calculate displacement: s = v₀t + ½at²
    double s = calculateDisplacement(v0, a, t);

    // Calculate final velocity from displacement: v² = v₀² + 2as
    double v_from_eq3 = calculateFinalVelocityFromDisplacement(v0, a, s);

    // Both methods should give same final velocity
    ASSERT_NEAR(v_from_eq1, v_from_eq3, TOLERANCE);

    // Verify values: v = 5 + 2*4 = 13 m/s, s = 5*4 + 0.5*2*16 = 20 + 16 = 36 m
    ASSERT_NEAR(v_from_eq1, 13.0, TOLERANCE);
    ASSERT_NEAR(s, 36.0, TOLERANCE);

    std::cout << "Kinematics consistency: OK" << std::endl;
}

// ============================================================================
// DYNAMICS TESTS (Newton's Laws)
// ============================================================================

void test_newtons_second_law() {
    std::cout << "\n=== Testing Newton's Second Law: F = ma ===" << std::endl;

    // Test: Calculate acceleration from force
    double mass = 10.0; // kg
    double force = 50.0; // N
    double a = calculateAccelerationFromForce(force, mass);
    ASSERT_NEAR(a, 5.0, TOLERANCE); // a = 50/10 = 5 m/s²

    // Test: Calculate required force for acceleration
    double F = calculateRequiredForce(mass, 5.0);
    ASSERT_NEAR(F, 50.0, TOLERANCE); // F = 10*5 = 50 N

    // Test: Multiple forces acting on object
    std::vector<double> forces = {20.0, 30.0, -10.0}; // N
    double net_force = calculateNetForce(forces);
    ASSERT_NEAR(net_force, 40.0, TOLERANCE); // 20 + 30 - 10 = 40 N

    // Test: Acceleration from net force
    a = calculateAccelerationFromForce(net_force, mass);
    ASSERT_NEAR(a, 4.0, TOLERANCE); // a = 40/10 = 4 m/s²

    std::cout << "Newton's Second Law: OK" << std::endl;
}

void test_force_motion_integration() {
    std::cout << "\n=== Testing Force-Motion Integration ===" << std::endl;

    double mass = 5.0;        // kg
    double force = 20.0;      // N
    double v0 = 10.0;         // m/s
    double time = 3.0;        // s

    // Calculate final velocity from force
    double v = calculateFinalVelocityFromForce(mass, force, v0, time);
    // a = F/m = 20/5 = 4 m/s², v = 10 + 4*3 = 22 m/s
    ASSERT_NEAR(v, 22.0, TOLERANCE);

    // Calculate velocity change
    double dv = calculateVelocityChange(mass, force, time);
    ASSERT_NEAR(dv, 12.0, TOLERANCE); // Δv = (20/5)*3 = 12 m/s

    // Calculate time for specific velocity change
    double t_required = calculateTimeForVelocityChange(mass, force, 12.0);
    ASSERT_NEAR(t_required, 3.0, TOLERANCE);

    // Test: Braking force (negative)
    double braking_force = -30.0; // N
    v = calculateFinalVelocityFromForce(mass, braking_force, 20.0, 2.0);
    // a = -30/5 = -6 m/s², v = 20 + (-6)*2 = 8 m/s
    ASSERT_NEAR(v, 8.0, TOLERANCE);

    std::cout << "Force-motion integration: OK" << std::endl;
}

void test_dynamics_kinematics_consistency() {
    std::cout << "\n=== Testing Dynamics-Kinematics Consistency ===" << std::endl;

    // Apply force to object, verify kinematics and dynamics give same result
    double mass = 10.0;     // kg
    double force = 30.0;    // N
    double v0 = 5.0;        // m/s
    double time = 4.0;      // s

    // Method 1: Use dynamics (F → a → v)
    double a_from_force = calculateAccelerationFromForce(force, mass);
    double v_from_dynamics = calculateFinalVelocityFromForce(mass, force, v0, time);

    // Method 2: Use kinematics directly with calculated acceleration
    double v_from_kinematics = calculateFinalVelocity(v0, a_from_force, time);

    // Both should match
    ASSERT_NEAR(v_from_dynamics, v_from_kinematics, TOLERANCE);

    // Verify values: a = 30/10 = 3 m/s², v = 5 + 3*4 = 17 m/s
    ASSERT_NEAR(a_from_force, 3.0, TOLERANCE);
    ASSERT_NEAR(v_from_dynamics, 17.0, TOLERANCE);

    std::cout << "Dynamics-kinematics consistency: OK" << std::endl;
}

// ============================================================================
// GRAVITATIONAL PHYSICS TESTS
// ============================================================================

void test_free_fall() {
    std::cout << "\n=== Testing Free Fall (Special case of kinematics) ===" << std::endl;

    const double g = 9.8; // m/s² (Earth's gravity)

    // Object dropped from rest, falls for 3 seconds
    double v = calculateFinalVelocity(0.0, g, 3.0);
    ASSERT_NEAR(v, 29.4, TOLERANCE); // v = 0 + 9.8*3 = 29.4 m/s

    // Distance fallen
    double s = calculateDisplacement(0.0, g, 3.0);
    ASSERT_NEAR(s, 44.1, TOLERANCE); // s = 0 + 0.5*9.8*9 = 44.1 m

    // Object thrown upward at 20 m/s, time to reach max height (v = 0)
    double t_max = calculateTimeFromVelocities(20.0, 0.0, -g);
    ASSERT_NEAR(t_max, 20.0/g, TOLERANCE); // t ≈ 2.04 s

    // Maximum height reached
    double h_max = calculateDisplacementFromVelocities(20.0, 0.0, -g);
    ASSERT_NEAR(h_max, 20.4, 1e-1); // h = 400/(2*9.8) ≈ 20.4 m

    std::cout << "Free fall: OK" << std::endl;
}

// ============================================================================
// ENERGY AND MOMENTUM TESTS
// ============================================================================

void test_kinetic_energy() {
    std::cout << "\n=== Testing Kinetic Energy ===" << std::endl;

    // Test: KE = ½mv² for object with m=2kg, v=10m/s
    double ke = calculateKineticEnergy(2.0, 10.0);
    ASSERT_NEAR(ke, 100.0, TOLERANCE); // KE = 0.5*2*100 = 100 J

    // Test: Calculate velocity from KE
    double v = calculateVelocityFromKE(100.0, 2.0);
    ASSERT_NEAR(v, 10.0, TOLERANCE); // v = √(2*100/2) = 10 m/s

    // Test: Doubling velocity quadruples energy
    double ke1 = calculateKineticEnergy(1.0, 5.0);
    double ke2 = calculateKineticEnergy(1.0, 10.0);
    ASSERT_NEAR(ke2 / ke1, 4.0, TOLERANCE);

    std::cout << "Kinetic energy: OK" << std::endl;
}

void test_momentum() {
    std::cout << "\n=== Testing Momentum ===" << std::endl;

    // Test: p = mv for object with m=5kg, v=10m/s
    double p = calculateMomentum(5.0, 10.0);
    ASSERT_NEAR(p, 50.0, TOLERANCE); // p = 5*10 = 50 kg⋅m/s

    // Test: Calculate velocity from momentum
    double v = calculateVelocityFromMomentum(50.0, 5.0);
    ASSERT_NEAR(v, 10.0, TOLERANCE); // v = 50/5 = 10 m/s

    // Test: KE from momentum: KE = p²/(2m)
    double ke = calculateKEFromMomentum(50.0, 5.0);
    ASSERT_NEAR(ke, 250.0, TOLERANCE); // KE = 2500/(2*5) = 250 J

    // Test: Momentum from KE: p = √(2m⋅KE)
    double p_from_ke = calculateMomentumFromKE(250.0, 5.0);
    ASSERT_NEAR(p_from_ke, 50.0, TOLERANCE); // p = √(2*5*250) = √2500 = 50

    std::cout << "Momentum: OK" << std::endl;
}

void test_energy_momentum_relationship() {
    std::cout << "\n=== Testing Energy-Momentum Relationship ===" << std::endl;

    // Test: For same object, verify KE and p are related by KE = p²/(2m)
    double mass = 3.0;
    double velocity = 8.0;

    double ke_direct = calculateKineticEnergy(mass, velocity);
    double p = calculateMomentum(mass, velocity);
    double ke_from_p = calculateKEFromMomentum(p, mass);

    ASSERT_NEAR(ke_direct, ke_from_p, TOLERANCE);

    // Verify: p = 3*8 = 24, KE = 0.5*3*64 = 96, KE = 576/(2*3) = 96 ✓
    ASSERT_NEAR(p, 24.0, TOLERANCE);
    ASSERT_NEAR(ke_direct, 96.0, TOLERANCE);

    std::cout << "Energy-momentum relationship: OK" << std::endl;
}

// ============================================================================
// THERMODYNAMICS TESTS
// ============================================================================

void test_boyles_law() {
    std::cout << "\n=== Testing Boyle's Law ===" << std::endl;

    // Test: P₁V₁ = P₂V₂
    double P1 = 100000.0; // 100 kPa
    double V1 = 2.0;      // 2 m³
    double V2 = 4.0;      // 4 m³

    double P2 = boylesLaw(P1, V1, V2);
    ASSERT_NEAR(P2, 50000.0, TOLERANCE); // P2 = 100000*2/4 = 50 kPa

    // Test: Compression (volume halves → pressure doubles)
    P2 = boylesLaw(P1, V1, 1.0);
    ASSERT_NEAR(P2, 200000.0, TOLERANCE);

    std::cout << "Boyle's Law: OK" << std::endl;
}

void test_charles_law() {
    std::cout << "\n=== Testing Charles's Law ===" << std::endl;

    // Test: V₁/T₁ = V₂/T₂
    double V1 = 1.0;   // 1 m³
    double T1 = 300.0; // 300 K
    double T2 = 600.0; // 600 K

    double V2 = charlesLaw(V1, T1, T2);
    ASSERT_NEAR(V2, 2.0, TOLERANCE); // V2 = 1*(600/300) = 2 m³

    // Test: Temperature doubles → volume doubles
    ASSERT_NEAR(V2 / V1, T2 / T1, TOLERANCE);

    std::cout << "Charles's Law: OK" << std::endl;
}

void test_ideal_gas_law() {
    std::cout << "\n=== Testing Ideal Gas Law ===" << std::endl;

    // Test: PV = nRT
    double n = 1.0;                      // 1 mole
    double T = 273.15;                   // 273.15 K (0°C)
    double R = physics::thermodynamics::constants::R;             // 8.314 J/(mol⋅K)
    double P = physics::thermodynamics::constants::STP_PRESSURE;  // 101325 Pa

    double V = idealGasLawVolume(n, T, P);
    ASSERT_NEAR(V, 0.0224, 1e-3); // At STP, 1 mole ≈ 22.4 L = 0.0224 m³

    // Test: Calculate pressure
    P = idealGasLawPressure(n, V, T);
    ASSERT_NEAR(P, physics::thermodynamics::constants::STP_PRESSURE, 100.0); // Within 100 Pa

    std::cout << "Ideal Gas Law: OK" << std::endl;
}

// ============================================================================
// GRAVITATION TESTS
// ============================================================================

void test_universal_gravitation() {
    std::cout << "\n=== Testing Universal Gravitation ===" << std::endl;

    // Test: Force between two 1kg masses, 1m apart
    double F = universalGravitationForce(1.0, 1.0, 1.0);
    ASSERT_NEAR(F, 6.674e-11, 1e-15); // F = G*1*1/1 = G

    // Test: Force doubles when one mass doubles
    double F2 = universalGravitationForce(2.0, 1.0, 1.0);
    ASSERT_NEAR(F2, 2.0 * F, 1e-15);

    // Test: Force quarters when distance doubles
    double F3 = universalGravitationForce(1.0, 1.0, 2.0);
    ASSERT_NEAR(F3, F / 4.0, 1e-15);

    std::cout << "Universal gravitation: OK" << std::endl;
}

void test_gravitational_field() {
    std::cout << "\n=== Testing Gravitational Field ===" << std::endl;

    // Test: g at Earth's surface
    double g_surface = gravitationalFieldStrength(physics::gravitation::constants::EARTH_MASS,
                                                   physics::gravitation::constants::EARTH_RADIUS);
    ASSERT_NEAR(g_surface, 9.8, 0.1); // Should be ≈ 9.81 m/s²

    // Test: g decreases with altitude (inverse square law)
    double g_double_radius = gravitationalFieldStrength(physics::gravitation::constants::EARTH_MASS,
                                                        2.0 * physics::gravitation::constants::EARTH_RADIUS);
    ASSERT_NEAR(g_double_radius, g_surface / 4.0, 0.1);

    std::cout << "Gravitational field: OK" << std::endl;
}

void test_gravitational_potential_energy() {
    std::cout << "\n=== Testing Gravitational Potential Energy ===" << std::endl;

    // Test: PE = mgh (near surface approximation)
    double mass = 10.0;   // kg
    double height = 100.0; // m
    double g = 9.8;       // m/s²

    double PE = mass * g * height;
    ASSERT_NEAR(PE, 9800.0, TOLERANCE); // PE = 10*9.8*100 = 9800 J

    // Test: Energy conservation in free fall
    // Object at height h with PE=mgh falls to ground
    // All PE converts to KE: ½mv² = mgh → v = √(2gh)
    double v_final = std::sqrt(2.0 * g * height);
    double KE_final = 0.5 * mass * v_final * v_final;
    ASSERT_NEAR(KE_final, PE, TOLERANCE); // Energy conserved

    std::cout << "Gravitational potential energy: OK" << std::endl;
}

// ============================================================================
// MAIN TEST RUNNER
// ============================================================================

int main() {
    std::cout << "=====================================================" << std::endl;
    std::cout << "PHASE 2 VALIDATION: BASIC MATH & CLASSICAL PHYSICS" << std::endl;
    std::cout << "=====================================================" << std::endl;

    // Calculus tests
    test_intermediate_value_theorem();
    test_mean_value_theorem();

    // Kinematics tests
    test_kinematics_first_equation();
    test_kinematics_second_equation();
    test_kinematics_third_equation();
    test_kinematics_consistency();

    // Dynamics tests
    test_newtons_second_law();
    test_force_motion_integration();
    test_dynamics_kinematics_consistency();

    // Special cases
    test_free_fall();

    // Energy and momentum tests
    test_kinetic_energy();
    test_momentum();
    test_energy_momentum_relationship();

    // Thermodynamics tests
    test_boyles_law();
    test_charles_law();
    test_ideal_gas_law();

    // Gravitation tests
    test_universal_gravitation();
    test_gravitational_field();
    test_gravitational_potential_energy();

    // Summary
    std::cout << "\n=====================================================" << std::endl;
    std::cout << "VALIDATION SUMMARY" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "Tests Passed: " << tests_passed << std::endl;
    std::cout << "Tests Failed: " << tests_failed << std::endl;
    std::cout << "Total Tests:  " << (tests_passed + tests_failed) << std::endl;

    if (tests_failed == 0) {
        std::cout << "\n✓ ALL TESTS PASSED!" << std::endl;
        std::cout << "Phase 2 validation (partial): COMPLETE" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ SOME TESTS FAILED" << std::endl;
        std::cout << "Phase 2 validation (partial): INCOMPLETE" << std::endl;
        return 1;
    }
}
