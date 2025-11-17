/**
 * Phase 5 Validation: Circular Motion
 *
 * Tests the circular_motion.hpp module functions.
 *
 * Coverage:
 * - Centripetal acceleration and force
 * - Angular velocity and conversions
 * - Period and frequency calculations
 * - Kinetic energy in circular motion
 * - Banking angles and vertical circles
 * - Tension calculations
 */

#include <iostream>
#include <cmath>
#include "../include/physics/circular_motion.hpp"

const double TOLERANCE = 1e-6;
const double LOOSE_TOLERANCE = 1e-3;

#define ASSERT_NEAR(actual, expected, tolerance) \
    do { \
        if (std::abs((actual) - (expected)) > (tolerance)) { \
            std::cerr << "FAIL: " << __LINE__ << ": " << #actual \
                      << " = " << (actual) << ", expected " << (expected) \
                      << " (diff: " << std::abs((actual) - (expected)) << ")" << std::endl; \
            return false; \
        } \
    } while(0)

#define ASSERT_TRUE(condition) \
    do { \
        if (!(condition)) { \
            std::cerr << "FAIL: " << __LINE__ << ": " << #condition << std::endl; \
            return false; \
        } \
    } while(0)

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 5: Circular Motion Validation ===" << std::endl;
    std::cout << std::endl;

    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "Running: " << name << "... ";
        if (test_func()) {
            std::cout << "PASS" << std::endl;
            tests_passed++;
            return true;
        } else {
            tests_failed++;
            return false;
        }
    };

    // Centripetal Acceleration Tests
    run_test("Centripetal acceleration basic calculation", []() {
        double v = 10.0;  // m/s
        double r = 5.0;   // m
        double ac = physics::circular_motion::calculateCentripetalAcceleration(v, r);
        ASSERT_NEAR(ac, 20.0, TOLERANCE);  // v²/r = 100/5 = 20 m/s²
        return true;
    });

    run_test("Centripetal acceleration increases with velocity", []() {
        double r = 10.0;
        double ac1 = physics::circular_motion::calculateCentripetalAcceleration(5.0, r);
        double ac2 = physics::circular_motion::calculateCentripetalAcceleration(10.0, r);
        ASSERT_NEAR(ac2, 4.0 * ac1, TOLERANCE);  // Quadratic in v
        return true;
    });

    run_test("Centripetal acceleration decreases with radius", []() {
        double v = 20.0;
        double ac1 = physics::circular_motion::calculateCentripetalAcceleration(v, 5.0);
        double ac2 = physics::circular_motion::calculateCentripetalAcceleration(v, 10.0);
        ASSERT_NEAR(ac2, 0.5 * ac1, TOLERANCE);  // Inverse in r
        return true;
    });

    run_test("Centripetal force basic calculation", []() {
        double m = 2.0;   // kg
        double v = 10.0;  // m/s
        double r = 5.0;   // m
        double fc = physics::circular_motion::calculateCentripetalForce(m, v, r);
        ASSERT_NEAR(fc, 40.0, TOLERANCE);  // m*v²/r = 2*100/5 = 40 N
        return true;
    });

    run_test("Centripetal force scales with mass", []() {
        double v = 10.0, r = 5.0;
        double fc1 = physics::circular_motion::calculateCentripetalForce(1.0, v, r);
        double fc2 = physics::circular_motion::calculateCentripetalForce(3.0, v, r);
        ASSERT_NEAR(fc2, 3.0 * fc1, TOLERANCE);
        return true;
    });

    run_test("Velocity from centripetal acceleration", []() {
        double ac = 20.0;  // m/s²
        double r = 5.0;    // m
        double v = physics::circular_motion::calculateVelocityFromAcceleration(ac, r);
        ASSERT_NEAR(v, 10.0, TOLERANCE);  // √(20*5) = 10 m/s
        return true;
    });

    // Angular Velocity Tests
    run_test("Angular velocity basic calculation", []() {
        double v = 10.0;  // m/s
        double r = 5.0;   // m
        double omega = physics::circular_motion::calculateAngularVelocity(v, r);
        ASSERT_NEAR(omega, 2.0, TOLERANCE);  // v/r = 10/5 = 2 rad/s
        return true;
    });

    run_test("Tangential velocity from angular velocity", []() {
        double omega = 2.0;  // rad/s
        double r = 5.0;      // m
        double v = physics::circular_motion::calculateTangentialVelocity(omega, r);
        ASSERT_NEAR(v, 10.0, TOLERANCE);  // ωr = 2*5 = 10 m/s
        return true;
    });

    run_test("Centripetal acceleration from angular velocity", []() {
        double omega = 2.0;  // rad/s
        double r = 5.0;      // m
        double ac = physics::circular_motion::calculateCentripetalAccelFromAngular(omega, r);
        ASSERT_NEAR(ac, 20.0, TOLERANCE);  // ω²r = 4*5 = 20 m/s²
        return true;
    });

    run_test("Angular velocity consistency check", []() {
        double v = 15.0, r = 3.0;
        double omega = physics::circular_motion::calculateAngularVelocity(v, r);
        double v_back = physics::circular_motion::calculateTangentialVelocity(omega, r);
        ASSERT_NEAR(v_back, v, TOLERANCE);
        return true;
    });

    // Period and Frequency Tests
    run_test("Period basic calculation", []() {
        double v = 10.0;  // m/s
        double r = 5.0;   // m
        double T = physics::circular_motion::calculatePeriod(v, r);
        ASSERT_NEAR(T, M_PI, TOLERANCE);  // 2πr/v = 2π*5/10 = π seconds
        return true;
    });

    run_test("Period from angular velocity", []() {
        double omega = 2.0;  // rad/s
        double T = physics::circular_motion::calculatePeriodFromAngular(omega);
        ASSERT_NEAR(T, M_PI, TOLERANCE);  // 2π/ω = 2π/2 = π seconds
        return true;
    });

    run_test("Frequency basic calculation", []() {
        double v = 10.0;  // m/s
        double r = 5.0;   // m
        double f = physics::circular_motion::calculateFrequency(v, r);
        ASSERT_NEAR(f, 1.0 / M_PI, TOLERANCE);  // v/(2πr) = 10/(2π*5)
        return true;
    });

    run_test("Angular velocity from period", []() {
        double T = M_PI;  // seconds
        double omega = physics::circular_motion::calculateAngularVelocityFromPeriod(T);
        ASSERT_NEAR(omega, 2.0, TOLERANCE);  // 2π/T = 2π/π = 2 rad/s
        return true;
    });

    run_test("Angular velocity from frequency", []() {
        double f = 1.0;  // Hz
        double omega = physics::circular_motion::calculateAngularVelocityFromFrequency(f);
        ASSERT_NEAR(omega, 2.0 * M_PI, TOLERANCE);  // 2πf
        return true;
    });

    run_test("Period and frequency are inverses", []() {
        double v = 20.0, r = 10.0;
        double T = physics::circular_motion::calculatePeriod(v, r);
        double f = physics::circular_motion::calculateFrequency(v, r);
        ASSERT_NEAR(T * f, 1.0, TOLERANCE);
        return true;
    });

    // Kinetic Energy Tests
    run_test("Kinetic energy basic calculation", []() {
        double m = 2.0;   // kg
        double v = 10.0;  // m/s
        double KE = physics::circular_motion::calculateKineticEnergy(m, v);
        ASSERT_NEAR(KE, 100.0, TOLERANCE);  // 0.5*m*v² = 0.5*2*100 = 100 J
        return true;
    });

    run_test("Kinetic energy from angular velocity", []() {
        double m = 2.0;      // kg
        double omega = 2.0;  // rad/s
        double r = 5.0;      // m
        double KE = physics::circular_motion::calculateKineticEnergyAngular(m, omega, r);
        ASSERT_NEAR(KE, 100.0, TOLERANCE);  // 0.5*m*r²*ω² = 0.5*2*25*4 = 100 J
        return true;
    });

    run_test("Kinetic energy consistency between methods", []() {
        double m = 3.0, r = 4.0, omega = 1.5;
        double v = physics::circular_motion::calculateTangentialVelocity(omega, r);
        double KE1 = physics::circular_motion::calculateKineticEnergy(m, v);
        double KE2 = physics::circular_motion::calculateKineticEnergyAngular(m, omega, r);
        ASSERT_NEAR(KE1, KE2, TOLERANCE);
        return true;
    });

    // Banking Angle Tests
    run_test("Banking angle basic calculation", []() {
        double v = 20.0;  // m/s
        double r = 100.0; // m
        double g = 9.81;
        double theta = physics::circular_motion::calculateBankingAngle(v, r, g);
        double expected = std::atan(v * v / (r * g));
        ASSERT_NEAR(theta, expected, TOLERANCE);
        return true;
    });

    run_test("Banking angle increases with speed", []() {
        double r = 100.0, g = 9.81;
        double theta1 = physics::circular_motion::calculateBankingAngle(10.0, r, g);
        double theta2 = physics::circular_motion::calculateBankingAngle(20.0, r, g);
        ASSERT_TRUE(theta2 > theta1);
        return true;
    });

    run_test("Banking angle decreases with radius", []() {
        double v = 20.0, g = 9.81;
        double theta1 = physics::circular_motion::calculateBankingAngle(v, 50.0, g);
        double theta2 = physics::circular_motion::calculateBankingAngle(v, 100.0, g);
        ASSERT_TRUE(theta2 < theta1);
        return true;
    });

    // Vertical Circle Tests
    run_test("Minimum velocity at top of loop", []() {
        double r = 5.0;  // m
        double g = 9.81;
        double v_min = physics::circular_motion::calculateMinVelocityTopOfLoop(r, g);
        ASSERT_NEAR(v_min, std::sqrt(g * r), TOLERANCE);
        return true;
    });

    run_test("Min velocity increases with loop radius", []() {
        double g = 9.81;
        double v1 = physics::circular_motion::calculateMinVelocityTopOfLoop(5.0, g);
        double v2 = physics::circular_motion::calculateMinVelocityTopOfLoop(20.0, g);
        ASSERT_NEAR(v2, 2.0 * v1, TOLERANCE);  // √(20/5) = 2
        return true;
    });

    run_test("Tension at bottom of vertical circle", []() {
        double m = 2.0, v = 10.0, r = 5.0, g = 9.81;
        double T = physics::circular_motion::calculateTensionAtBottom(m, v, r, g);
        double expected = m * v * v / r + m * g;
        ASSERT_NEAR(T, expected, TOLERANCE);
        return true;
    });

    run_test("Tension at top of vertical circle", []() {
        double m = 2.0, v = 10.0, r = 5.0, g = 9.81;
        double T = physics::circular_motion::calculateTensionAtTop(m, v, r, g);
        double expected = m * v * v / r - m * g;
        ASSERT_NEAR(T, expected, TOLERANCE);
        return true;
    });

    run_test("Tension at top is zero at minimum velocity", []() {
        double m = 2.0, r = 5.0, g = 9.81;
        double v_min = physics::circular_motion::calculateMinVelocityTopOfLoop(r, g);
        double T = physics::circular_motion::calculateTensionAtTop(m, v_min, r, g);
        ASSERT_NEAR(T, 0.0, TOLERANCE);
        return true;
    });

    run_test("Tension bottom exceeds tension top", []() {
        double m = 2.0, v = 10.0, r = 5.0, g = 9.81;
        double T_bottom = physics::circular_motion::calculateTensionAtBottom(m, v, r, g);
        double T_top = physics::circular_motion::calculateTensionAtTop(m, v, r, g);
        ASSERT_TRUE(T_bottom > T_top);
        return true;
    });

    // Centripetal Force Alternative Expressions
    run_test("Centripetal force from angular velocity", []() {
        double m = 2.0, omega = 3.0, r = 4.0;
        double Fc = physics::circular_motion::calculateCentripetalForceAngular(m, omega, r);
        ASSERT_NEAR(Fc, m * omega * omega * r, TOLERANCE);
        return true;
    });

    run_test("Centripetal force from period", []() {
        double m = 2.0, r = 5.0, T = M_PI;
        double Fc = physics::circular_motion::calculateCentripetalForceFromPeriod(m, r, T);
        double expected = 4.0 * M_PI * M_PI * m * r / (T * T);
        ASSERT_NEAR(Fc, expected, TOLERANCE);
        return true;
    });

    run_test("Centripetal force from frequency", []() {
        double m = 2.0, r = 5.0, f = 2.0;
        double Fc = physics::circular_motion::calculateCentripetalForceFromFrequency(m, r, f);
        double expected = 4.0 * M_PI * M_PI * m * f * f * r;
        ASSERT_NEAR(Fc, expected, TOLERANCE);
        return true;
    });

    run_test("Centripetal force consistency across methods", []() {
        double m = 3.0, v = 12.0, r = 6.0;
        double Fc1 = physics::circular_motion::calculateCentripetalForce(m, v, r);

        double omega = physics::circular_motion::calculateAngularVelocity(v, r);
        double Fc2 = physics::circular_motion::calculateCentripetalForceAngular(m, omega, r);

        ASSERT_NEAR(Fc1, Fc2, TOLERANCE);
        return true;
    });

    // RPM Conversions
    run_test("RPM to rad per sec conversion", []() {
        double rpm = 60.0;
        double rad_per_sec = physics::circular_motion::rpmToRadPerSec(rpm);
        ASSERT_NEAR(rad_per_sec, 2.0 * M_PI, TOLERANCE);  // 60 RPM = 1 rev/s = 2π rad/s
        return true;
    });

    run_test("Rad per sec to RPM conversion", []() {
        double rad_per_sec = 2.0 * M_PI;
        double rpm = physics::circular_motion::radPerSecToRpm(rad_per_sec);
        ASSERT_NEAR(rpm, 60.0, TOLERANCE);
        return true;
    });

    run_test("RPM conversion round trip", []() {
        double rpm_orig = 120.0;
        double rad_per_sec = physics::circular_motion::rpmToRadPerSec(rpm_orig);
        double rpm_back = physics::circular_motion::radPerSecToRpm(rad_per_sec);
        ASSERT_NEAR(rpm_back, rpm_orig, TOLERANCE);
        return true;
    });

    // Comprehensive Integration Tests
    run_test("Satellite in circular orbit", []() {
        // Satellite at 400 km altitude (ISS-like orbit)
        double r = 6.771e6;  // Earth radius + altitude (m)
        double v = 7670.0;   // m/s (orbital speed)

        double ac = physics::circular_motion::calculateCentripetalAcceleration(v, r);
        ASSERT_TRUE(ac > 8.0 && ac < 9.0);  // ~8.7 m/s² at this altitude
        return true;
    });

    run_test("Car on banked curve realistic values", []() {
        double v = 25.0;   // m/s (~90 km/h)
        double r = 200.0;  // m
        double g = 9.81;
        double theta = physics::circular_motion::calculateBankingAngle(v, r, g);
        ASSERT_TRUE(theta > 0.1 && theta < 0.5);  // Reasonable banking angle
        return true;
    });

    run_test("Washing machine spin cycle", []() {
        double rpm = 1000.0;  // 1000 RPM
        double omega = physics::circular_motion::rpmToRadPerSec(rpm);
        double r = 0.25;  // 25 cm drum radius

        double v = physics::circular_motion::calculateTangentialVelocity(omega, r);
        double ac = physics::circular_motion::calculateCentripetalAccelFromAngular(omega, r);

        ASSERT_TRUE(ac > 100.0);  // Very high centripetal acceleration
        return true;
    });

    run_test("Ferris wheel period calculation", []() {
        double r = 50.0;  // m (radius of large ferris wheel)
        double v = 1.0;   // m/s (slow speed)
        double T = physics::circular_motion::calculatePeriod(v, r);
        ASSERT_NEAR(T, 100.0 * M_PI, LOOSE_TOLERANCE);  // ~5.2 minutes
        return true;
    });

    // Summary
    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 5 Results: Circular Motion" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    return (tests_failed == 0) ? 0 : 1;
}
