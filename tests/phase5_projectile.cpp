/**
 * Phase 5 Validation: Projectile Motion
 */

#include <iostream>
#include <cmath>
#include "../include/physics/projectile.hpp"

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

    std::cout << "=== Phase 5: Projectile Motion Validation ===" << std::endl;

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

    run_test("Horizontal velocity component", []() {
        double v0 = 50.0;  // m/s
        double angle = M_PI / 4;  // 45 degrees
        double vx = physics::projectile::calculateHorizontalVelocity(v0, angle);
        ASSERT_NEAR(vx, v0 / std::sqrt(2.0), TOLERANCE);
        return true;
    });

    run_test("Vertical velocity component", []() {
        double v0 = 50.0;
        double angle = M_PI / 4;
        double vy = physics::projectile::calculateVerticalVelocity(v0, angle);
        ASSERT_NEAR(vy, v0 / std::sqrt(2.0), TOLERANCE);
        return true;
    });

    run_test("Time of flight", []() {
        double v0 = 50.0;
        double angle = M_PI / 4;
        double g = 9.81;
        double T = physics::projectile::calculateTimeOfFlight(v0, angle, g);
        double expected = 2.0 * v0 * std::sin(angle) / g;
        ASSERT_NEAR(T, expected, TOLERANCE);
        return true;
    });

    run_test("Maximum height", []() {
        double v0 = 50.0;
        double angle = M_PI / 2;  // Straight up
        double g = 9.81;
        double h = physics::projectile::calculateMaxHeight(v0, angle, g);
        double expected = v0 * v0 / (2.0 * g);
        ASSERT_NEAR(h, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Range calculation", []() {
        double v0 = 50.0;
        double angle = M_PI / 4;  // 45 degrees for max range
        double g = 9.81;
        double R = physics::projectile::calculateRange(v0, angle, g);
        double expected = v0 * v0 * std::sin(2.0 * angle) / g;
        ASSERT_NEAR(R, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("45 degrees gives maximum range", []() {
        double v0 = 50.0;
        double g = 9.81;
        double R45 = physics::projectile::calculateRange(v0, M_PI / 4, g);
        double R30 = physics::projectile::calculateRange(v0, M_PI / 6, g);
        double R60 = physics::projectile::calculateRange(v0, M_PI / 3, g);
        ASSERT_TRUE(R45 > R30);
        ASSERT_TRUE(R45 > R60);
        return true;
    });

    run_test("Horizontal displacement at time t", []() {
        double v0 = 50.0;
        double angle = M_PI / 6;
        double t = 2.0;
        double x = physics::projectile::calculateHorizontalPosition(v0, angle, t);
        double expected = v0 * std::cos(angle) * t;
        ASSERT_NEAR(x, expected, TOLERANCE);
        return true;
    });

    run_test("Vertical displacement at time t", []() {
        double v0 = 50.0;
        double angle = M_PI / 6;
        double t = 2.0;
        double g = 9.81;
        double y = physics::projectile::calculateVerticalPosition(v0, angle, t, g);
        double expected = v0 * std::sin(angle) * t - 0.5 * g * t * t;
        ASSERT_NEAR(y, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Velocity magnitude at time t", []() {
        double v0 = 50.0;
        double angle = M_PI / 4;
        double t = 1.0;
        double g = 9.81;
        double v = physics::projectile::calculateSpeedAtTime(v0, angle, t, g);
        double vx = v0 * std::cos(angle);
        double vy = v0 * std::sin(angle) - g * t;
        double expected = std::sqrt(vx * vx + vy * vy);
        ASSERT_NEAR(v, expected, TOLERANCE);
        return true;
    });

    run_test("Speed at launch equals initial velocity", []() {
        double v0 = 50.0;
        double angle = M_PI / 4;
        double t = 0.0;  // At launch
        double v = physics::projectile::calculateSpeedAtTime(v0, angle, t);
        ASSERT_NEAR(v, v0, TOLERANCE);
        return true;
    });

    run_test("Time to reach maximum height", []() {
        double v0 = 50.0;
        double angle = M_PI / 3;
        double g = 9.81;
        double t_max = physics::projectile::calculateTimeToMaxHeight(v0, angle, g);
        double expected = v0 * std::sin(angle) / g;
        ASSERT_NEAR(t_max, expected, TOLERANCE);
        return true;
    });

    run_test("Time to max height is half total time", []() {
        double v0 = 50.0;
        double angle = M_PI / 4;
        double g = 9.81;
        double t_max = physics::projectile::calculateTimeToMaxHeight(v0, angle, g);
        double T = physics::projectile::calculateTimeOfFlight(v0, angle, g);
        ASSERT_NEAR(t_max, T / 2.0, TOLERANCE);
        return true;
    });

    run_test("Projectile symmetry - same height at t and T-t", []() {
        double v0 = 50.0;
        double angle = M_PI / 4;
        double g = 9.81;
        double T = physics::projectile::calculateTimeOfFlight(v0, angle, g);
        double t = 1.0;
        double y1 = physics::projectile::calculateVerticalPosition(v0, angle, t, g);
        double y2 = physics::projectile::calculateVerticalPosition(v0, angle, T - t, g);
        ASSERT_NEAR(y1, y2, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Range equals horizontal position at landing", []() {
        double v0 = 50.0;
        double angle = M_PI / 6;
        double g = 9.81;
        double R = physics::projectile::calculateRange(v0, angle, g);
        double T = physics::projectile::calculateTimeOfFlight(v0, angle, g);
        double x = physics::projectile::calculateHorizontalPosition(v0, angle, T);
        ASSERT_NEAR(R, x, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Vertical velocity at max height is zero", []() {
        double v0 = 50.0;
        double angle = M_PI / 3;
        double g = 9.81;
        double t_max = physics::projectile::calculateTimeToMaxHeight(v0, angle, g);
        double vy = v0 * std::sin(angle) - g * t_max;
        ASSERT_NEAR(vy, 0.0, TOLERANCE);
        return true;
    });

    run_test("Energy conservation in projectile motion", []() {
        double v0 = 50.0;
        double angle = M_PI / 4;
        double m = 2.0;  // kg
        double g = 9.81;
        double t = 2.0;

        double KE0 = 0.5 * m * v0 * v0;
        double y = physics::projectile::calculateVerticalPosition(v0, angle, t, g);
        double v = physics::projectile::calculateSpeedAtTime(v0, angle, t, g);
        double KE = 0.5 * m * v * v;
        double PE = m * g * y;

        ASSERT_NEAR(KE0, KE + PE, LOOSE_TOLERANCE * m * v0 * v0);
        return true;
    });

    run_test("Horizontal velocity remains constant", []() {
        double v0 = 50.0;
        double angle = M_PI / 6;
        double vx0 = physics::projectile::calculateHorizontalVelocity(v0, angle);

        // Check at different times
        for (double t = 0; t <= 5.0; t += 1.0) {
            double vx = v0 * std::cos(angle);  // Should be constant
            ASSERT_NEAR(vx, vx0, TOLERANCE);
        }
        return true;
    });

    run_test("Angle of launch affects time of flight", []() {
        double v0 = 50.0;
        double g = 9.81;
        double T30 = physics::projectile::calculateTimeOfFlight(v0, M_PI / 6, g);
        double T60 = physics::projectile::calculateTimeOfFlight(v0, M_PI / 3, g);
        double T90 = physics::projectile::calculateTimeOfFlight(v0, M_PI / 2, g);

        ASSERT_TRUE(T90 > T60);
        ASSERT_TRUE(T60 > T30);
        return true;
    });

    run_test("Initial velocity affects range quadratically", []() {
        double angle = M_PI / 4;
        double g = 9.81;
        double R1 = physics::projectile::calculateRange(10.0, angle, g);
        double R2 = physics::projectile::calculateRange(20.0, angle, g);
        ASSERT_NEAR(R2, 4.0 * R1, LOOSE_TOLERANCE * R1);
        return true;
    });

    run_test("Complementary angles give same range", []() {
        double v0 = 50.0;
        double g = 9.81;
        double R30 = physics::projectile::calculateRange(v0, M_PI / 6, g);
        double R60 = physics::projectile::calculateRange(v0, M_PI / 3, g);
        ASSERT_NEAR(R30, R60, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Realistic baseball throw", []() {
        double v0 = 40.0;  // m/s (~90 mph)
        double angle = M_PI / 4;
        double g = 9.81;
        double R = physics::projectile::calculateRange(v0, angle, g);
        ASSERT_TRUE(R > 150.0 && R < 170.0);  // ~163 m
        return true;
    });

    run_test("Realistic soccer kick", []() {
        double v0 = 25.0;  // m/s
        double angle = M_PI / 6;  // 30 degrees
        double g = 9.81;
        double R = physics::projectile::calculateRange(v0, angle, g);
        double h_max = physics::projectile::calculateMaxHeight(v0, angle, g);
        ASSERT_TRUE(R > 50.0 && R < 60.0);
        ASSERT_TRUE(h_max > 5.0 && h_max < 10.0);
        return true;
    });

    run_test("Cannon projectile calculation", []() {
        double v0 = 300.0;  // m/s
        double angle = M_PI / 4;
        double g = 9.81;
        double R = physics::projectile::calculateRange(v0, angle, g);
        ASSERT_TRUE(R > 9000.0);  // Over 9 km range
        return true;
    });

    run_test("Trajectory at different times", []() {
        double v0 = 50.0;
        double angle = M_PI / 4;
        double g = 9.81;
        double T = physics::projectile::calculateTimeOfFlight(v0, angle, g);

        for (double t = 0; t <= T; t += T / 10.0) {
            double y = physics::projectile::calculateVerticalPosition(v0, angle, t, g);
            ASSERT_TRUE(y >= -TOLERANCE);  // Should not go below ground
        }
        return true;
    });

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    return (tests_failed == 0) ? 0 : 1;
}
