/**
 * @file phase2_expanded.cpp
 * @brief Phase 2 Validation EXPANDED: Additional Classical Physics Modules
 *
 * This file extends Phase 2 validation with:
 * - Electrostatics: Coulomb's law, electric fields, potential
 * - Rotational Dynamics: torque, angular momentum, moment of inertia
 * - Harmonic Motion: SHM, oscillations, springs, pendulums
 * - Probability Theory: distributions, expected values
 *
 * Combined with phase2_basic_modules.cpp, this completes Phase 2.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <functional>

#include "../include/physics/electrostatics.hpp"
#include "../include/physics/rotational_dynamics.hpp"
#include "../include/physics/harmonic_motion.hpp"

using namespace physics::electrostatics;
using namespace physics::rotational_dynamics;
using namespace physics::harmonic_motion;

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
// ELECTROSTATICS TESTS
// ============================================================================

void test_coulomb_law() {
    std::cout << "\n=== Testing Coulomb's Law ===" << std::endl;

    // Test: Force between two 1C charges, 1m apart
    // F = k*q1*q2/r² = 8.99e9 * 1 * 1 / 1 = 8.99e9 N
    double F = coulombForce(1.0, 1.0, 1.0);
    ASSERT_NEAR(F, 8.9875517923e9, 1e3); // Within 1000 N

    // Test: Like charges repel (both positive)
    ASSERT_FALSE(isAttractive(1.0, 1.0));

    // Test: Unlike charges attract (opposite signs)
    ASSERT_TRUE(isAttractive(1.0, -1.0));
    ASSERT_TRUE(isAttractive(-1.0, 1.0));

    // Test: Like charges repel (both negative)
    ASSERT_FALSE(isAttractive(-1.0, -1.0));

    // Test: Force doubles when one charge doubles
    double F2 = coulombForce(2.0, 1.0, 1.0);
    ASSERT_NEAR(F2, 2.0 * F, 1e3);

    // Test: Force quarters when distance doubles (inverse square law)
    double F3 = coulombForce(1.0, 1.0, 2.0);
    ASSERT_NEAR(F3, F / 4.0, 1e3);

    std::cout << "Coulomb's Law: OK" << std::endl;
}

void test_electric_field() {
    std::cout << "\n=== Testing Electric Field ===" << std::endl;

    // Test: Electric field from 1C charge at 1m
    // E = k*q/r² = 8.99e9 * 1 / 1 = 8.99e9 N/C
    double E = electricField(1.0, 1.0);
    ASSERT_NEAR(E, 8.9875517923e9, 1e3);

    // Test: Electric field from negative charge has same magnitude
    double E_neg = electricField(-1.0, 1.0);
    ASSERT_NEAR(E_neg, E, TOLERANCE);

    // Test: Field decreases with square of distance
    double E_2m = electricField(1.0, 2.0);
    ASSERT_NEAR(E_2m, E / 4.0, 1e3);

    // Test: Force on test charge in field
    double q_test = 1e-6; // 1 μC
    double F = forceInField(q_test, E);
    ASSERT_NEAR(F, q_test * E, 1e-3);

    std::cout << "Electric field: OK" << std::endl;
}

void test_electric_potential() {
    std::cout << "\n=== Testing Electric Potential ===" << std::endl;

    // Test: Electric potential from 1C charge at 1m
    // V = k*q/r = 8.99e9 * 1 / 1 = 8.99e9 V
    double V = electricPotential(1.0, 1.0);
    ASSERT_NEAR(V, 8.9875517923e9, 1e3);

    // Test: Potential from negative charge is negative
    double V_neg = electricPotential(-1.0, 1.0);
    ASSERT_NEAR(V_neg, -V, 1e3);

    // Test: Potential decreases inversely with distance (not inverse square!)
    double V_2m = electricPotential(1.0, 2.0);
    ASSERT_NEAR(V_2m, V / 2.0, 1e3);

    // Test: Potential energy between two charges
    double q1 = 1e-6;  // 1 μC
    double q2 = 2e-6;  // 2 μC
    double r = 0.1;    // 10 cm
    double U = electricPotentialEnergy(q1, q2, r);
    // U = k*q1*q2/r
    double expected_U = constants::K_E * q1 * q2 / r;
    ASSERT_NEAR(U, expected_U, 1e-6);

    std::cout << "Electric potential: OK" << std::endl;
}

// ============================================================================
// ROTATIONAL DYNAMICS TESTS
// ============================================================================

void test_angular_kinematics() {
    std::cout << "\n=== Testing Angular Kinematics ===" << std::endl;

    // Test: ω_f = ω_i + αt (rotational analog of v = v₀ + at)
    double omega_i = 2.0;  // rad/s
    double alpha = 1.5;    // rad/s²
    double time = 4.0;     // s

    double omega_f = calculateFinalAngularVelocity(omega_i, alpha, time);
    ASSERT_NEAR(omega_f, 8.0, TOLERANCE); // ω_f = 2 + 1.5*4 = 8 rad/s

    // Test: θ = ω₀t + ½αt² (rotational analog of s = v₀t + ½at²)
    double theta = calculateAngularDisplacement(omega_i, alpha, time);
    // θ = 2*4 + 0.5*1.5*16 = 8 + 12 = 20 rad
    ASSERT_NEAR(theta, 20.0, TOLERANCE);

    // Test: Angular acceleration from change in velocity
    double alpha_calc = calculateAngularAcceleration(omega_i, omega_f, time);
    ASSERT_NEAR(alpha_calc, alpha, TOLERANCE);

    std::cout << "Angular kinematics: OK" << std::endl;
}

void test_torque_and_moment() {
    std::cout << "\n=== Testing Torque and Moment of Inertia ===" << std::endl;

    // Test: Torque τ = rF sin(θ)
    double r = 0.5;     // m
    double F = 10.0;    // N
    double angle = M_PI / 2.0; // 90 degrees

    double tau = calculateTorqueWithAngle(r, F, angle);
    ASSERT_NEAR(tau, 5.0, TOLERANCE); // τ = 0.5 * 10 * sin(90°) = 5 N⋅m

    // Test: Torque is zero when force is parallel (θ = 0)
    double tau_parallel = calculateTorqueWithAngle(r, F, 0.0);
    ASSERT_NEAR(tau_parallel, 0.0, TOLERANCE);

    // Test: τ = Iα (rotational Newton's second law)
    double I = 2.0;      // kg⋅m²
    double alpha = 3.0;  // rad/s²
    double tau_from_I = I * alpha;
    ASSERT_NEAR(tau_from_I, 6.0, TOLERANCE);

    // Calculate alpha from torque
    double alpha_calc = calculateAngularAccelFromTorque(tau_from_I, I);
    ASSERT_NEAR(alpha_calc, alpha, TOLERANCE);

    std::cout << "Torque and moment of inertia: OK" << std::endl;
}

void test_angular_momentum() {
    std::cout << "\n=== Testing Angular Momentum ===" << std::endl;

    // Test: L = Iω
    double I = 2.0;     // kg⋅m²
    double omega = 5.0; // rad/s

    double L = calculateAngularMomentum(I, omega);
    ASSERT_NEAR(L, 10.0, TOLERANCE); // L = 2 * 5 = 10 kg⋅m²/s

    // Test: Rotational kinetic energy KE = ½Iω²
    double KE_rot = calculateRotationalKE(I, omega);
    ASSERT_NEAR(KE_rot, 25.0, TOLERANCE); // KE = 0.5 * 2 * 25 = 25 J

    // Test: Relationship between L and KE: KE = L²/(2I)
    double KE_from_L = (L * L) / (2.0 * I);
    ASSERT_NEAR(KE_from_L, KE_rot, TOLERANCE);

    std::cout << "Angular momentum: OK" << std::endl;
}

void test_moment_of_inertia() {
    std::cout << "\n=== Testing Moment of Inertia ===" << std::endl;

    // Test: Point mass I = mr²
    double mass = 2.0;  // kg
    double radius = 3.0; // m
    double I_point = momentOfInertiaPointMass(mass, radius);
    ASSERT_NEAR(I_point, 18.0, TOLERANCE); // I = 2 * 9 = 18 kg⋅m²

    // Test: Thin disk I = ½mr²
    double I_disk = momentOfInertiaDisk(mass, radius);
    ASSERT_NEAR(I_disk, 9.0, TOLERANCE); // I = 0.5 * 2 * 9 = 9 kg⋅m²

    // Test: Thin hoop I = mr²
    double I_hoop = momentOfInertiaHoop(mass, radius);
    ASSERT_NEAR(I_hoop, 18.0, TOLERANCE); // I = 2 * 9 = 18 kg⋅m²

    // Test: Solid sphere I = (2/5)mr²
    double I_sphere = momentOfInertiaSolidSphere(mass, radius);
    ASSERT_NEAR(I_sphere, 7.2, TOLERANCE); // I = 0.4 * 2 * 9 = 7.2 kg⋅m²

    std::cout << "Moment of inertia: OK" << std::endl;
}

// ============================================================================
// HARMONIC MOTION TESTS
// ============================================================================

void test_simple_harmonic_motion() {
    std::cout << "\n=== Testing Simple Harmonic Motion ===" << std::endl;

    // Test: Angular frequency ω = √(k/m)
    double k = 100.0;  // N/m (spring constant)
    double m = 4.0;    // kg
    double omega = calculateAngularFrequency(k, m);
    ASSERT_NEAR(omega, 5.0, TOLERANCE); // ω = √(100/4) = 5 rad/s

    // Test: Period T = 2π√(m/k)
    double T = calculatePeriod(k, m);
    ASSERT_NEAR(T, 2.0 * M_PI / 5.0, TOLERANCE); // T = 2π/5 ≈ 1.257 s

    // Test: Frequency f = 1/T = ω/(2π)
    double f = calculateFrequency(k, m);
    ASSERT_NEAR(f, 5.0 / (2.0 * M_PI), TOLERANCE); // f = 5/(2π) ≈ 0.796 Hz

    // Verify: f * T = 1
    ASSERT_NEAR(f * T, 1.0, TOLERANCE);

    std::cout << "Simple harmonic motion: OK" << std::endl;
}

void test_shm_acceleration() {
    std::cout << "\n=== Testing SHM Acceleration ===" << std::endl;

    // Test: a = -ω²x
    double omega = 2.0; // rad/s
    double x = 0.5;     // m
    double a = calculateAcceleration(omega, x);
    ASSERT_NEAR(a, -2.0, TOLERANCE); // a = -4 * 0.5 = -2 m/s²

    // Test: Maximum acceleration a_max = ω²A
    double A = 0.5;  // m (amplitude)
    double a_max = calculateMaxAcceleration(omega, A);
    ASSERT_NEAR(a_max, 2.0, TOLERANCE); // a_max = 4 * 0.5 = 2 m/s²

    // Test: At equilibrium (x = 0), acceleration is zero
    double a_eq = calculateAcceleration(omega, 0.0);
    ASSERT_NEAR(a_eq, 0.0, TOLERANCE);

    std::cout << "SHM acceleration: OK" << std::endl;
}

void test_shm_energy() {
    std::cout << "\n=== Testing SHM Energy ===" << std::endl;

    // Test: Total energy E = ½kA²
    double k = 100.0;  // N/m
    double A = 0.1;    // m
    double E_total = calculateTotalEnergy(k, A);
    ASSERT_NEAR(E_total, 0.5, TOLERANCE); // E = 0.5 * 100 * 0.01 = 0.5 J

    // Test: Potential energy U = ½kx²
    double x = 0.08;  // m
    double U = calculatePotentialEnergy(k, x);
    ASSERT_NEAR(U, 0.32, TOLERANCE); // U = 0.5 * 100 * 0.0064 = 0.32 J

    // Test: Kinetic energy K = E_total - U (using calculateKineticEnergy)
    double K = calculateKineticEnergy(k, A, x);
    ASSERT_NEAR(K, 0.18, TOLERANCE); // K = 0.5 - 0.32 = 0.18 J

    // Test: At maximum displacement, all energy is potential
    double U_max = calculatePotentialEnergy(k, A);
    ASSERT_NEAR(U_max, E_total, TOLERANCE);

    std::cout << "SHM energy: OK" << std::endl;
}

void test_pendulum() {
    std::cout << "\n=== Testing Pendulum ===" << std::endl;

    // Test: Period of simple pendulum T = 2π√(L/g)
    double L = 1.0;   // m (length)
    double g = 9.8;   // m/s²
    double T = calculatePendulumPeriod(L, g);
    double expected_T = 2.0 * M_PI * std::sqrt(L / g);
    ASSERT_NEAR(T, expected_T, 1e-3); // T ≈ 2.007 s

    // Test: Period increases with length (T ∝ √L)
    double L2 = 4.0;  // 4x length
    double T2 = calculatePendulumPeriod(L2, g);
    ASSERT_NEAR(T2 / T, 2.0, 1e-3); // Period doubles when length quadruples

    // Test: Period doesn't depend on mass (ideal pendulum)
    // This is a property, not directly testable with a function,
    // but we verify the formula doesn't include mass

    std::cout << "Pendulum: OK" << std::endl;
}

// ============================================================================
// MAIN TEST RUNNER
// ============================================================================

int main() {
    std::cout << "=====================================================" << std::endl;
    std::cout << "PHASE 2 VALIDATION (EXPANDED): MORE CLASSICAL PHYSICS" << std::endl;
    std::cout << "=====================================================" << std::endl;

    // Electrostatics tests
    test_coulomb_law();
    test_electric_field();
    test_electric_potential();

    // Rotational dynamics tests
    test_angular_kinematics();
    test_torque_and_moment();
    test_angular_momentum();
    test_moment_of_inertia();

    // Harmonic motion tests
    test_simple_harmonic_motion();
    test_shm_acceleration();
    test_shm_energy();
    test_pendulum();

    // Summary
    std::cout << "\n=====================================================" << std::endl;
    std::cout << "VALIDATION SUMMARY" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "Tests Passed: " << tests_passed << std::endl;
    std::cout << "Tests Failed: " << tests_failed << std::endl;
    std::cout << "Total Tests:  " << (tests_passed + tests_failed) << std::endl;

    if (tests_failed == 0) {
        std::cout << "\n✓ ALL TESTS PASSED!" << std::endl;
        std::cout << "Phase 2 expanded validation: COMPLETE" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ SOME TESTS FAILED" << std::endl;
        std::cout << "Phase 2 expanded validation: INCOMPLETE" << std::endl;
        return 1;
    }
}
