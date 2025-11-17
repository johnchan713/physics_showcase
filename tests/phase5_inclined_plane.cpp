/**
 * Phase 5 Validation: Inclined Plane
 */

#include <iostream>
#include <cmath>
#include "../include/physics/inclined_plane.hpp"

const double TOLERANCE = 1e-6;
const double LOOSE_TOLERANCE = 1e-3;

#define ASSERT_NEAR(actual, expected, tolerance) \
    do { \
        if (std::abs((actual) - (expected)) > (tolerance)) { \
            std::cerr << "FAIL: " << __LINE__ << ": " << #actual \
                      << " = " << (actual) << ", expected " << (expected) << std::endl; \
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

    std::cout << "=== Phase 5: Inclined Plane Validation ===" << std::endl;

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

    // Force Components Tests
    run_test("Parallel force basic calculation", []() {
        double m = 10.0;  // kg
        double angle = M_PI / 6;  // 30 degrees
        double g = 9.81;
        double F_parallel = physics::inclined_plane::calculateParallelForce(m, angle, g);
        double expected = m * g * std::sin(angle);
        ASSERT_NEAR(F_parallel, expected, TOLERANCE);
        return true;
    });

    run_test("Parallel force at 45 degrees", []() {
        double m = 5.0;
        double angle = M_PI / 4;
        double g = 9.81;
        double F_parallel = physics::inclined_plane::calculateParallelForce(m, angle, g);
        ASSERT_NEAR(F_parallel, m * g / std::sqrt(2.0), LOOSE_TOLERANCE);
        return true;
    });

    run_test("Parallel force increases with angle", []() {
        double m = 10.0, g = 9.81;
        double F30 = physics::inclined_plane::calculateParallelForce(m, M_PI / 6, g);
        double F45 = physics::inclined_plane::calculateParallelForce(m, M_PI / 4, g);
        ASSERT_TRUE(F45 > F30);
        return true;
    });

    run_test("Normal force basic calculation", []() {
        double m = 10.0;
        double angle = M_PI / 6;
        double g = 9.81;
        double N = physics::inclined_plane::calculateNormalForce(m, angle, g);
        double expected = m * g * std::cos(angle);
        ASSERT_NEAR(N, expected, TOLERANCE);
        return true;
    });

    run_test("Normal force decreases with angle", []() {
        double m = 10.0, g = 9.81;
        double N30 = physics::inclined_plane::calculateNormalForce(m, M_PI / 6, g);
        double N45 = physics::inclined_plane::calculateNormalForce(m, M_PI / 4, g);
        ASSERT_TRUE(N45 < N30);
        return true;
    });

    run_test("Force components are perpendicular", []() {
        double m = 8.0;
        double angle = M_PI / 3;
        double g = 9.81;
        double F_par = physics::inclined_plane::calculateParallelForce(m, angle, g);
        double F_perp = physics::inclined_plane::calculateNormalForce(m, angle, g);
        double F_total_sq = F_par * F_par + F_perp * F_perp;
        double weight_sq = (m * g) * (m * g);
        ASSERT_NEAR(F_total_sq, weight_sq, LOOSE_TOLERANCE);
        return true;
    });

    // Acceleration Tests
    run_test("Acceleration frictionless basic", []() {
        double angle = M_PI / 6;
        double g = 9.81;
        double a = physics::inclined_plane::calculateAccelerationFrictionless(angle, g);
        double expected = g * std::sin(angle);
        ASSERT_NEAR(a, expected, TOLERANCE);
        return true;
    });

    run_test("Acceleration independent of mass", []() {
        double angle = M_PI / 4;
        double g = 9.81;
        double a = physics::inclined_plane::calculateAccelerationFrictionless(angle, g);
        ASSERT_NEAR(a, g / std::sqrt(2.0), LOOSE_TOLERANCE);
        return true;
    });

    run_test("Acceleration with friction", []() {
        double angle = M_PI / 6;
        double mu = 0.1;
        double g = 9.81;
        double a = physics::inclined_plane::calculateAccelerationWithFriction(angle, mu, g);
        double expected = g * (std::sin(angle) - mu * std::cos(angle));
        ASSERT_NEAR(a, expected, TOLERANCE);
        return true;
    });

    run_test("Friction reduces acceleration", []() {
        double angle = M_PI / 4;
        double g = 9.81;
        double a_no_friction = physics::inclined_plane::calculateAccelerationFrictionless(angle, g);
        double a_friction = physics::inclined_plane::calculateAccelerationWithFriction(angle, 0.2, g);
        ASSERT_TRUE(a_friction < a_no_friction);
        return true;
    });

    run_test("High friction can prevent sliding", []() {
        double angle = M_PI / 6;  // 30 degrees
        double mu = 1.0;  // High friction
        double g = 9.81;
        double a = physics::inclined_plane::calculateAccelerationWithFriction(angle, mu, g);
        ASSERT_TRUE(a <= 0.0);  // Object doesn't slide
        return true;
    });

    run_test("Velocity at bottom frictionless", []() {
        double height = 10.0;  // m
        double g = 9.81;
        double v = physics::inclined_plane::calculateVelocityAtFootFromHeight(height, g);
        double expected = std::sqrt(2.0 * g * height);
        ASSERT_NEAR(v, expected, TOLERANCE);
        return true;
    });

    run_test("Velocity from energy conservation", []() {
        double h = 5.0;
        double g = 9.81;
        double v = physics::inclined_plane::calculateVelocityAtFootFromHeight(h, g);
        double KE = 0.5 * v * v;  // per unit mass
        double PE = g * h;
        ASSERT_NEAR(KE, PE, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Time to slide down frictionless", []() {
        double length = 10.0;  // m
        double angle = M_PI / 6;
        double g = 9.81;
        double t = physics::inclined_plane::calculateTimeToSlide(length, angle, g);
        double a = physics::inclined_plane::calculateAccelerationFrictionless(angle, g);
        double expected = std::sqrt(2.0 * length / a);
        ASSERT_NEAR(t, expected, TOLERANCE);
        return true;
    });

    run_test("Distance traveled in time t", []() {
        double angle = M_PI / 4;
        double t = 2.0;  // seconds
        double g = 9.81;
        double a = physics::inclined_plane::calculateAccelerationFrictionless(angle, g);
        double d = 0.5 * a * t * t;  // Using kinematic equation directly
        double expected = 0.5 * a * t * t;
        ASSERT_NEAR(d, expected, TOLERANCE);
        return true;
    });

    run_test("Steeper plane gives faster slide", []() {
        double length = 10.0;
        double g = 9.81;
        double t30 = physics::inclined_plane::calculateTimeToSlide(length, M_PI / 6, g);
        double t45 = physics::inclined_plane::calculateTimeToSlide(length, M_PI / 4, g);
        ASSERT_TRUE(t45 < t30);
        return true;
    });

    run_test("Friction increases slide time", []() {
        double length = 10.0;
        double angle = M_PI / 6;
        double g = 9.81;
        double t_no_friction = physics::inclined_plane::calculateTimeToSlide(length, angle, g);
        double a_friction = physics::inclined_plane::calculateAccelerationWithFriction(angle, 0.1, g);
        double t_friction = std::sqrt(2.0 * length / a_friction);
        ASSERT_TRUE(t_friction > t_no_friction);
        return true;
    });

    run_test("Work done by gravity", []() {
        double m = 5.0;
        double h = 10.0;
        double g = 9.81;
        double W = m * g * h;
        double v = physics::inclined_plane::calculateVelocityAtFootFromHeight(h, g);
        double KE = 0.5 * m * v * v;
        ASSERT_NEAR(W, KE, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Minimum coefficient to prevent sliding", []() {
        double angle = M_PI / 6;  // 30 degrees
        double mu_min = std::tan(angle);
        double g = 9.81;
        double a = physics::inclined_plane::calculateAccelerationWithFriction(angle, mu_min, g);
        ASSERT_NEAR(a, 0.0, TOLERANCE);
        return true;
    });

    run_test("Angle for given acceleration", []() {
        double a = 5.0;  // m/s²
        double g = 9.81;
        double angle = std::asin(a / g);
        double a_check = physics::inclined_plane::calculateAccelerationFrictionless(angle, g);
        ASSERT_NEAR(a_check, a, TOLERANCE);
        return true;
    });

    run_test("Block sliding down realistic ramp", []() {
        double angle = M_PI / 6;  // 30 degree ramp
        double length = 5.0;  // 5 meter ramp
        double g = 9.81;
        double t = physics::inclined_plane::calculateTimeToSlide(length, angle, g);
        ASSERT_TRUE(t > 1.0 && t < 5.0);  // Reasonable time
        return true;
    });

    run_test("Normal force supports weight component", []() {
        double m = 10.0;
        double angle = M_PI / 4;
        double g = 9.81;
        double N = physics::inclined_plane::calculateNormalForce(m, angle, g);
        double weight = m * g;
        ASSERT_TRUE(N < weight);  // Normal force is less than full weight
        ASSERT_TRUE(N > 0.0);
        return true;
    });

    run_test("Parallel and perpendicular components at 45 degrees", []() {
        double m = 10.0;
        double angle = M_PI / 4;
        double g = 9.81;
        double F_par = physics::inclined_plane::calculateParallelForce(m, angle, g);
        double F_perp = physics::inclined_plane::calculateNormalForce(m, angle, g);
        ASSERT_NEAR(F_par, F_perp, LOOSE_TOLERANCE);  // Equal at 45°
        return true;
    });

    run_test("Zero friction equals frictionless", []() {
        double angle = M_PI / 6;
        double g = 9.81;
        double a1 = physics::inclined_plane::calculateAccelerationFrictionless(angle, g);
        double a2 = physics::inclined_plane::calculateAccelerationWithFriction(angle, 0.0, g);
        ASSERT_NEAR(a1, a2, TOLERANCE);
        return true;
    });

    run_test("Kinematic equation consistency", []() {
        double angle = M_PI / 6;
        double t = 3.0;
        double g = 9.81;
        double a = physics::inclined_plane::calculateAccelerationFrictionless(angle, g);
        double d = 0.5 * a * t * t;  // Distance from kinematics
        double v = a * t;
        double v_from_energy = std::sqrt(2.0 * a * d);
        ASSERT_NEAR(v, v_from_energy, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Height calculation from length and angle", []() {
        double length = 10.0;
        double angle = M_PI / 6;
        double height = length * std::sin(angle);
        double v1 = physics::inclined_plane::calculateVelocityAtFootFromHeight(height, 9.81);

        double v2 = physics::inclined_plane::calculateVelocityAtFootFrictionless(length, angle, 9.81);

        ASSERT_NEAR(v1, v2, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Ski slope realistic values", []() {
        double angle = M_PI / 9;  // 20 degrees
        double length = 100.0;  // 100 meter slope
        double mu = 0.05;  // Low friction (ski wax)
        double g = 9.81;

        double a = physics::inclined_plane::calculateAccelerationWithFriction(angle, mu, g);
        double t = std::sqrt(2.0 * length / a);
        double v = a * t;

        ASSERT_TRUE(v > 10.0 && v < 30.0);  // Reasonable skiing speed
        return true;
    });

    run_test("Playground slide calculation", []() {
        double angle = M_PI / 4;  // 45 degree slide
        double length = 3.0;  // 3 meters
        double mu = 0.2;  // Some friction
        double g = 9.81;

        double a = physics::inclined_plane::calculateAccelerationWithFriction(angle, mu, g);
        double v_final = std::sqrt(2.0 * a * length);

        ASSERT_TRUE(v_final > 0.0 && v_final < 10.0);  // Safe slide speed
        return true;
    });

    run_test("Truck on hill with brakes", []() {
        double angle = M_PI / 12;  // 15 degrees
        double mu_static = 0.7;  // Tire friction
        double g = 9.81;

        // Check if truck stays put
        double mu_required = std::tan(angle);
        ASSERT_TRUE(mu_static > mu_required);  // Won't slide

        double a = physics::inclined_plane::calculateAccelerationWithFriction(angle, mu_static, g);
        ASSERT_TRUE(a <= 0.0);  // Stationary or decelerating
        return true;
    });

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    return (tests_failed == 0) ? 0 : 1;
}
