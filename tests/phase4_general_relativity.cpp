/**
 * Phase 4 Validation: General Relativity
 *
 * Tests the general_relativity.hpp module functions.
 *
 * Coverage:
 * - Metric tensors (Minkowski, Schwarzschild)
 * - Schwarzschild solution properties (radii, velocities, time dilation, redshift)
 * - Gravitational waves (strain amplitude, polarization)
 * - Schwarzschild geometry (event horizon, photon sphere, ISCO)
 * - Basic curvature properties
 */

#include <iostream>
#include <cmath>
#include <array>
#include <limits>  // Required for std::numeric_limits used in general_relativity.hpp
#include "../include/physics/general_relativity.hpp"

// Test tolerance
const double TOLERANCE = 1e-6;
const double LOOSE_TOLERANCE = 1e-3;
const double VERY_LOOSE = 0.01;

// Test macros
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

// Physical constants from the module
using physics::general_relativity::c;
using physics::general_relativity::G;

// Type aliases
using Vector4 = physics::general_relativity::Vector4;
using Tensor2 = physics::general_relativity::Tensor2;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: General Relativity Validation ===" << std::endl;
    std::cout << std::endl;

    // Helper lambda to run tests
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

    // ========================================
    // Physical Constants Tests
    // ========================================

    run_test("Speed of light constant", []() {
        // c = 299792458 m/s (exact by definition)
        ASSERT_NEAR(c, 299792458.0, TOLERANCE);
        return true;
    });

    run_test("Gravitational constant", []() {
        // G = 6.67430×10⁻¹¹ m³/(kg·s²) (CODATA 2018)
        ASSERT_NEAR(G, 6.67430e-11, 1e-16);
        return true;
    });

    // ========================================
    // Minkowski Metric Tests
    // ========================================

    run_test("Minkowski metric diagonal components", []() {
        physics::general_relativity::MinkowskiMetric eta;
        Vector4 x = {0, 0, 0, 0};
        auto g = eta.covariant(x);

        // η = diag(-1, 1, 1, 1)
        ASSERT_NEAR(g[0][0], -1.0, TOLERANCE);
        ASSERT_NEAR(g[1][1], 1.0, TOLERANCE);
        ASSERT_NEAR(g[2][2], 1.0, TOLERANCE);
        ASSERT_NEAR(g[3][3], 1.0, TOLERANCE);
        return true;
    });

    run_test("Minkowski metric off-diagonal components", []() {
        physics::general_relativity::MinkowskiMetric eta;
        Vector4 x = {0, 0, 0, 0};
        auto g = eta.covariant(x);

        // Off-diagonal should be zero
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i != j) {
                    ASSERT_NEAR(g[i][j], 0.0, TOLERANCE);
                }
            }
        }
        return true;
    });

    run_test("Minkowski metric is position-independent", []() {
        physics::general_relativity::MinkowskiMetric eta;
        Vector4 x1 = {0, 0, 0, 0};
        Vector4 x2 = {100, 50, 25, 10};

        auto g1 = eta.covariant(x1);
        auto g2 = eta.covariant(x2);

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                ASSERT_NEAR(g1[i][j], g2[i][j], TOLERANCE);
            }
        }
        return true;
    });

    run_test("Minkowski line element (timelike)", []() {
        physics::general_relativity::MinkowskiMetric eta;
        Vector4 x = {0, 0, 0, 0};
        Vector4 dx = {1, 0, 0, 0};  // Pure time displacement

        double ds2 = eta.lineElement(x, dx);
        // ds² = -c²dt² for timelike interval
        ASSERT_NEAR(ds2, -1.0, TOLERANCE);
        return true;
    });

    run_test("Minkowski line element (spacelike)", []() {
        physics::general_relativity::MinkowskiMetric eta;
        Vector4 x = {0, 0, 0, 0};
        Vector4 dx = {0, 1, 0, 0};  // Pure spatial displacement

        double ds2 = eta.lineElement(x, dx);
        // ds² = dx² for spacelike interval
        ASSERT_NEAR(ds2, 1.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Schwarzschild Metric Tests
    // ========================================

    run_test("Schwarzschild radius formula", []() {
        double M_sun = 1.989e30;  // kg
        physics::general_relativity::SchwarzschildMetric schwarzschild(M_sun);

        double r_s = schwarzschild.schwarzschildRadius();
        // r_s = 2GM/c² ≈ 2953 m for the Sun
        double expected = 2.0 * G * M_sun / (c * c);

        ASSERT_NEAR(r_s, expected, TOLERANCE);
        return true;
    });

    run_test("Schwarzschild radius for solar mass", []() {
        double M_sun = 1.989e30;  // kg
        physics::general_relativity::SchwarzschildMetric schwarzschild(M_sun);

        double r_s = schwarzschild.schwarzschildRadius();
        // Should be approximately 2953 meters
        ASSERT_NEAR(r_s, 2953.0, 10.0);  // Within 10 meters
        return true;
    });

    run_test("Photon sphere at 1.5 r_s", []() {
        double M = 1e30;  // 1 solar mass
        physics::general_relativity::SchwarzschildMetric schwarzschild(M);

        double r_photon = schwarzschild.photonSphereRadius();
        double r_s = schwarzschild.schwarzschildRadius();

        // Photon sphere at r = 1.5 r_s
        ASSERT_NEAR(r_photon, 1.5 * r_s, TOLERANCE);
        return true;
    });

    run_test("ISCO at 3 r_s", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildMetric schwarzschild(M);

        double r_isco = schwarzschild.innerMostStableOrbit();
        double r_s = schwarzschild.schwarzschildRadius();

        // ISCO at r = 3 r_s for Schwarzschild
        ASSERT_NEAR(r_isco, 3.0 * r_s, TOLERANCE);
        return true;
    });

    run_test("Schwarzschild metric at infinity approaches Minkowski", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildMetric schwarzschild(M);

        // At r → ∞, Schwarzschild → Minkowski
        double r = 1e15;  // Very far from black hole
        Vector4 x = {0, r, M_PI/2, 0};  // (t, r, θ, φ)
        auto g = schwarzschild.covariant(x);

        // Should approach diag(-1, 1, r², r²sin²θ)
        ASSERT_NEAR(g[0][0], -1.0, LOOSE_TOLERANCE);
        ASSERT_NEAR(g[1][1], 1.0, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Schwarzschild Solution Properties Tests
    // ========================================

    run_test("Orbital velocity formula", []() {
        double M_earth = 5.972e24;  // kg
        physics::general_relativity::SchwarzschildSolution earth(M_earth);

        double r = 6.371e6 + 400e3;  // Earth radius + ISS altitude
        double v_orbital = earth.orbitalVelocity(r);

        // v = √(GM/r) ≈ 7670 m/s for ISS
        double expected = std::sqrt(G * M_earth / r);

        ASSERT_NEAR(v_orbital, expected, TOLERANCE);
        return true;
    });

    run_test("Orbital velocity at ISS altitude", []() {
        double M_earth = 5.972e24;
        physics::general_relativity::SchwarzschildSolution earth(M_earth);

        double r = 6.771e6;  // ~400 km altitude
        double v = earth.orbitalVelocity(r);

        // Should be approximately 7670 m/s
        ASSERT_NEAR(v, 7670.0, 50.0);  // Within 50 m/s
        return true;
    });

    run_test("Escape velocity formula", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r = 1e10;  // Some radius
        double v_escape = bh.escapVelocity(r);

        // v_esc = √(2GM/r)
        double expected = std::sqrt(2.0 * G * M / r);

        ASSERT_NEAR(v_escape, expected, TOLERANCE);
        return true;
    });

    run_test("Escape velocity near event horizon approaches c", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r_s = 2.0 * G * M / (c * c);
        double v_escape = bh.escapVelocity(r_s * 1.001);  // Just outside horizon

        // Near event horizon, escape velocity approaches c
        ASSERT_TRUE(v_escape > 0.99 * c);
        ASSERT_TRUE(v_escape <= c);
        return true;
    });

    run_test("Escape velocity > orbital velocity", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r = 1e10;
        double v_orbital = bh.orbitalVelocity(r);
        double v_escape = bh.escapVelocity(r);

        // v_escape = √2 × v_orbital
        ASSERT_NEAR(v_escape, std::sqrt(2.0) * v_orbital, TOLERANCE);
        return true;
    });

    // ========================================
    // Gravitational Time Dilation Tests
    // ========================================

    run_test("Proper time dilation formula", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r_s = 2.0 * G * M / (c * c);
        double r = 10 * r_s;  // Far from event horizon
        double t_coord = 100.0;  // 100 seconds coordinate time

        double tau = bh.properTimeAtRadius(r, t_coord);

        // τ = t × √(1 - r_s/r)
        double expected = t_coord * std::sqrt(1.0 - r_s / r);

        ASSERT_NEAR(tau, expected, TOLERANCE);
        return true;
    });

    run_test("Time dilation: proper time < coordinate time", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r_s = 2.0 * G * M / (c * c);
        double r = 5 * r_s;
        double t = 100.0;

        double tau = bh.properTimeAtRadius(r, t);

        // Proper time should be less than coordinate time near massive objects
        ASSERT_TRUE(tau < t);
        return true;
    });

    run_test("Time dilation approaches 1 far from source", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r = 1e15;  // Very far
        double t = 100.0;

        double tau = bh.properTimeAtRadius(r, t);

        // τ/t → 1 as r → ∞
        double ratio = tau / t;
        ASSERT_NEAR(ratio, 1.0, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Gravitational Redshift Tests
    // ========================================

    run_test("Gravitational redshift formula", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r_s = 2.0 * G * M / (c * c);
        double r_emit = 3 * r_s;
        double r_obs = 10 * r_s;

        double z = bh.redshift(r_emit, r_obs);

        // z = √(1 - r_s/r_obs) / √(1 - r_s/r_emit) - 1
        double f_emit = std::sqrt(1.0 - r_s / r_emit);
        double f_obs = std::sqrt(1.0 - r_s / r_obs);
        double expected = (f_obs / f_emit) - 1.0;

        ASSERT_NEAR(z, expected, TOLERANCE);
        return true;
    });

    run_test("Gravitational redshift is positive (light loses energy)", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r_s = 2.0 * G * M / (c * c);
        double r_emit = 2 * r_s;  // Deep in potential well
        double r_obs = 100 * r_s;  // Far away observer

        double z = bh.redshift(r_emit, r_obs);

        // Photon climbing out of gravity well is redshifted (z > 0)
        ASSERT_TRUE(z > 0);
        return true;
    });

    run_test("Gravitational redshift at same radius is zero", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r_s = 2.0 * G * M / (c * c);
        double r = 10 * r_s;

        double z = bh.redshift(r, r);

        // No redshift if emission and observation at same radius
        ASSERT_NEAR(z, 0.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Photon Orbit Tests
    // ========================================

    run_test("Photon orbit radius", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r_photon = bh.photonOrbitRadius();
        double r_s = 2.0 * G * M / (c * c);

        // Unstable photon orbit at r = 1.5 r_s
        ASSERT_NEAR(r_photon, 1.5 * r_s, TOLERANCE);
        return true;
    });

    run_test("Photon orbit is outside event horizon", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r_photon = bh.photonOrbitRadius();
        double r_s = 2.0 * G * M / (c * c);

        // Photon orbit must be outside event horizon
        ASSERT_TRUE(r_photon > r_s);
        return true;
    });

    // ========================================
    // Tidal Force Tests
    // ========================================

    run_test("Tidal force formula", []() {
        double M = 1e30;  // Black hole mass
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r = 1e10;       // Distance from center
        double m_obj = 70.0;   // 70 kg human
        double length = 1.8;   // 1.8 m tall

        double F_tidal = bh.tidalForce(r, m_obj, length);

        // F_tidal = 2GM·m·L / r³
        double expected = 2.0 * G * M * m_obj * length / (r * r * r);

        ASSERT_NEAR(F_tidal, expected, TOLERANCE);
        return true;
    });

    run_test("Tidal force increases closer to black hole", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double m = 70.0;
        double L = 1.8;

        double F_far = bh.tidalForce(1e11, m, L);
        double F_near = bh.tidalForce(1e10, m, L);

        // Tidal force ∝ 1/r³, so closer = stronger
        ASSERT_TRUE(F_near > F_far);
        return true;
    });

    run_test("Tidal force scales with mass", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r = 1e10;
        double L = 1.8;

        double F1 = bh.tidalForce(r, 70.0, L);
        double F2 = bh.tidalForce(r, 140.0, L);

        // Tidal force ∝ mass
        ASSERT_NEAR(F2, 2.0 * F1, TOLERANCE);
        return true;
    });

    // ========================================
    // Gravitational Wave Tests
    // ========================================

    run_test("Gravitational wave linearized metric", []() {
        Vector4 x = {0, 0, 0, 0};
        double amplitude = 1e-21;
        double frequency = 100.0;  // 100 Hz
        int polarization = 0;  // Plus polarization

        auto h = physics::general_relativity::GravitationalWaves::linearizedMetric(
            x, amplitude, frequency, polarization);

        // At t=0, z=0: h_+ has h_11 = A, h_22 = -A
        ASSERT_NEAR(h[1][1], amplitude, TOLERANCE);
        ASSERT_NEAR(h[2][2], -amplitude, TOLERANCE);
        return true;
    });

    run_test("Gravitational wave plus polarization", []() {
        Vector4 x = {0, 0, 0, 0};
        double A = 1e-21;
        double f = 100.0;

        auto h_plus = physics::general_relativity::GravitationalWaves::linearizedMetric(
            x, A, f, 0);

        // Plus polarization: h_11 = -h_22, h_12 = 0
        ASSERT_NEAR(h_plus[1][1], -h_plus[2][2], TOLERANCE);
        ASSERT_NEAR(h_plus[1][2], 0.0, TOLERANCE);
        return true;
    });

    run_test("Gravitational wave cross polarization", []() {
        Vector4 x = {0, 0, 0, 0};
        double A = 1e-21;
        double f = 100.0;

        auto h_cross = physics::general_relativity::GravitationalWaves::linearizedMetric(
            x, A, f, 1);

        // Cross polarization: h_11 = h_22 = 0, h_12 ≠ 0
        ASSERT_NEAR(h_cross[1][1], 0.0, TOLERANCE);
        ASSERT_NEAR(h_cross[2][2], 0.0, TOLERANCE);
        ASSERT_TRUE(std::abs(h_cross[1][2]) > 0);
        return true;
    });

    run_test("Gravitational wave strain amplitude", []() {
        double M1 = 30.0 * 1.989e30;  // 30 solar masses
        double M2 = 30.0 * 1.989e30;  // 30 solar masses
        double distance = 1e25;       // ~1 Gpc
        double frequency = 100.0;     // 100 Hz

        double h = physics::general_relativity::GravitationalWaves::strainAmplitude(
            M1, M2, distance, frequency);

        // Should be positive and physically reasonable
        // For 30+30 M_sun at ~1 Gpc and 100 Hz
        ASSERT_TRUE(h > 0);
        ASSERT_TRUE(h < 1e-10);  // Should be much smaller than 1
        return true;
    });

    run_test("GW strain decreases with distance", []() {
        double M1 = 30.0 * 1.989e30;
        double M2 = 30.0 * 1.989e30;
        double f = 100.0;

        double h_near = physics::general_relativity::GravitationalWaves::strainAmplitude(
            M1, M2, 1e24, f);
        double h_far = physics::general_relativity::GravitationalWaves::strainAmplitude(
            M1, M2, 2e24, f);

        // h ∝ 1/d, so doubling distance halves strain
        ASSERT_NEAR(h_far, 0.5 * h_near, 0.1 * h_near);
        return true;
    });

    // ========================================
    // Perihelion Shift Test
    // ========================================

    run_test("Mercury perihelion shift", []() {
        double M_sun = 1.989e30;
        physics::general_relativity::MinkowskiMetric metric;
        physics::general_relativity::GeodesicSolver geodesic(metric);

        double a = 5.791e10;    // Mercury semi-major axis (m)
        double e = 0.2056;      // Mercury eccentricity
        int orbits_per_century = 415;  // Mercury orbits in 100 years

        double shift = geodesic.perihelionShift(M_sun, a, e, orbits_per_century);

        // Mercury's perihelion shift: ~43 arcseconds per century
        // This is the famous GR prediction that matched observations
        ASSERT_NEAR(shift, 43.0, 5.0);  // Within 5 arcseconds
        return true;
    });

    // ========================================
    // Consistency Tests
    // ========================================

    run_test("Schwarzschild radius scales linearly with mass", []() {
        double M1 = 1e30;
        double M2 = 2e30;

        physics::general_relativity::SchwarzschildMetric bh1(M1);
        physics::general_relativity::SchwarzschildMetric bh2(M2);

        double r1 = bh1.schwarzschildRadius();
        double r2 = bh2.schwarzschildRadius();

        // r_s ∝ M
        ASSERT_NEAR(r2, 2.0 * r1, TOLERANCE);
        return true;
    });

    run_test("Orbital velocity decreases with distance", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double v1 = bh.orbitalVelocity(1e10);
        double v2 = bh.orbitalVelocity(2e10);

        // v ∝ 1/√r
        ASSERT_NEAR(v2, v1 / std::sqrt(2.0), TOLERANCE);
        return true;
    });

    run_test("Time dilation consistency with redshift", []() {
        double M = 1e30;
        physics::general_relativity::SchwarzschildSolution bh(M);

        double r_s = 2.0 * G * M / (c * c);
        double r = 5 * r_s;
        double r_inf = 1e15;  // Effectively infinity

        // Time dilation factor
        double tau_t = bh.properTimeAtRadius(r, 1.0);

        // Redshift from r to infinity
        double z = bh.redshift(r, r_inf);

        // (1 + z) should equal 1/√(1 - r_s/r)
        double expected_ratio = 1.0 / std::sqrt(1.0 - r_s / r);
        ASSERT_NEAR(1.0 + z, expected_ratio, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: General Relativity" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "✓ All general relativity tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Metric tensors (Minkowski, Schwarzschild)" << std::endl;
        std::cout << "  - Schwarzschild geometry (r_s, photon sphere, ISCO)" << std::endl;
        std::cout << "  - Orbital mechanics (velocities, stability)" << std::endl;
        std::cout << "  - Gravitational time dilation" << std::endl;
        std::cout << "  - Gravitational redshift" << std::endl;
        std::cout << "  - Tidal forces (spaghettification)" << std::endl;
        std::cout << "  - Gravitational waves (strain, polarization)" << std::endl;
        std::cout << "  - Perihelion precession (Mercury test)" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
