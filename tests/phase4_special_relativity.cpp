/**
 * Phase 4 Validation: Special Relativity
 *
 * Tests the special_relativity.hpp module functions.
 *
 * Coverage:
 * - Lorentz transformations
 * - Time dilation and length contraction
 * - Relativistic velocity addition
 * - Relativistic energy and momentum
 * - Energy-momentum relation
 * - Relativistic Doppler effect
 * - Four-vectors and spacetime intervals
 * - Invariant mass and rapidity
 */

#include <iostream>
#include <cmath>
#include <vector>
#include "../include/physics/special_relativity.hpp"

// Test tolerance
const double TOLERANCE = 1e-6;
const double LOOSE_TOLERANCE = 1e-3;

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

// Import functions from special_relativity namespace
using physics::special_relativity::lorentzFactor;
using physics::special_relativity::beta;
using physics::special_relativity::timeDilation;
using physics::special_relativity::lengthContraction;
using physics::special_relativity::properTime;
using physics::special_relativity::properLength;
using physics::special_relativity::velocityAddition;
using physics::special_relativity::relativisticEnergy;
using physics::special_relativity::restEnergy;
using physics::special_relativity::relativisticKineticEnergy;
using physics::special_relativity::relativisticMomentum;
using physics::special_relativity::energyFromMomentum;
using physics::special_relativity::photonEnergy;
using physics::special_relativity::relativisticDoppler;
using physics::special_relativity::cosmologicalRedshift;
using physics::special_relativity::spacetimeInterval;
using physics::special_relativity::invariantMass;
using physics::special_relativity::rapidity;
using physics::special_relativity::velocityFromRapidity;

// Use fully qualified names for constants to avoid namespace ambiguity
const double c = physics::special_relativity::constants::SPEED_OF_LIGHT;  // 299792458 m/s

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Special Relativity Validation ===" << std::endl;
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
    // Lorentz Factor Tests
    // ========================================

    run_test("Lorentz factor at rest (v=0)", []() {
        double gamma = lorentzFactor(0.0);
        ASSERT_NEAR(gamma, 1.0, TOLERANCE);
        return true;
    });

    run_test("Lorentz factor at low velocity", []() {
        double v = 1000.0;  // 1 km/s << c
        double gamma = lorentzFactor(v);
        // At v << c, gamma ≈ 1 + (1/2)(v/c)²
        double expected = 1.0 + 0.5 * (v/c) * (v/c);
        ASSERT_NEAR(gamma, expected, 1e-9);
        return true;
    });

    run_test("Lorentz factor at 0.6c", []() {
        double v = 0.6 * c;
        double gamma = lorentzFactor(v);
        // γ = 1/√(1 - 0.36) = 1/√0.64 = 1/0.8 = 1.25
        ASSERT_NEAR(gamma, 1.25, TOLERANCE);
        return true;
    });

    run_test("Lorentz factor at 0.8c", []() {
        double v = 0.8 * c;
        double gamma = lorentzFactor(v);
        // γ = 1/√(1 - 0.64) = 1/√0.36 = 1/0.6 = 5/3
        ASSERT_NEAR(gamma, 5.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Lorentz factor at 0.9c", []() {
        double v = 0.9 * c;
        double gamma = lorentzFactor(v);
        // γ = 1/√(1 - 0.81) = 1/√0.19 ≈ 2.294
        double expected = 1.0 / std::sqrt(1.0 - 0.81);
        ASSERT_NEAR(gamma, expected, TOLERANCE);
        return true;
    });

    run_test("Beta calculation", []() {
        double v = 0.5 * c;
        double b = beta(v);
        ASSERT_NEAR(b, 0.5, TOLERANCE);
        return true;
    });

    run_test("Lorentz factor always >= 1", []() {
        for (double v = 0; v < 0.99*c; v += 0.1*c) {
            double gamma = lorentzFactor(v);
            ASSERT_TRUE(gamma >= 1.0);
        }
        return true;
    });

    // ========================================
    // Time Dilation Tests
    // ========================================

    run_test("Time dilation at rest", []() {
        double properTime = 10.0;  // 10 seconds
        double dilatedTime = timeDilation(properTime, 0.0);
        ASSERT_NEAR(dilatedTime, properTime, TOLERANCE);
        return true;
    });

    run_test("Time dilation at 0.6c", []() {
        double properTime = 10.0;  // 10 seconds proper time
        double v = 0.6 * c;
        double dilatedTime = timeDilation(properTime, v);
        // γ = 1.25, so Δt = 1.25 * 10 = 12.5 s
        ASSERT_NEAR(dilatedTime, 12.5, TOLERANCE);
        return true;
    });

    run_test("Time dilation at 0.8c", []() {
        double properTime = 10.0;
        double v = 0.8 * c;
        double dilatedTime = timeDilation(properTime, v);
        // γ = 5/3, so Δt = (5/3) * 10 = 16.667 s
        ASSERT_NEAR(dilatedTime, 50.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Proper time calculation", []() {
        double dilatedTime = 12.5;  // 12.5 seconds observed
        double v = 0.6 * c;
        double tau = properTime(dilatedTime, v);
        // τ = Δt/γ = 12.5/1.25 = 10 s
        ASSERT_NEAR(tau, 10.0, TOLERANCE);
        return true;
    });

    run_test("Time dilation - proper time consistency", []() {
        double tau = 5.0;
        double v = 0.5 * c;
        double t = timeDilation(tau, v);
        double tau_back = properTime(t, v);
        ASSERT_NEAR(tau_back, tau, TOLERANCE);
        return true;
    });

    // ========================================
    // Length Contraction Tests
    // ========================================

    run_test("Length contraction at rest", []() {
        double properLength = 100.0;  // 100 m
        double contractedLength = lengthContraction(properLength, 0.0);
        ASSERT_NEAR(contractedLength, properLength, TOLERANCE);
        return true;
    });

    run_test("Length contraction at 0.6c", []() {
        double L0 = 100.0;  // 100 m proper length
        double v = 0.6 * c;
        double L = lengthContraction(L0, v);
        // L = L₀/γ = 100/1.25 = 80 m
        ASSERT_NEAR(L, 80.0, TOLERANCE);
        return true;
    });

    run_test("Length contraction at 0.8c", []() {
        double L0 = 100.0;
        double v = 0.8 * c;
        double L = lengthContraction(L0, v);
        // L = L₀/γ = 100/(5/3) = 60 m
        ASSERT_NEAR(L, 60.0, TOLERANCE);
        return true;
    });

    run_test("Proper length calculation", []() {
        double L = 80.0;  // 80 m contracted
        double v = 0.6 * c;
        double L0 = properLength(L, v);
        // L₀ = L * γ = 80 * 1.25 = 100 m
        ASSERT_NEAR(L0, 100.0, TOLERANCE);
        return true;
    });

    run_test("Length contraction - proper length consistency", []() {
        double L0 = 50.0;
        double v = 0.7 * c;
        double L = lengthContraction(L0, v);
        double L0_back = properLength(L, v);
        ASSERT_NEAR(L0_back, L0, TOLERANCE);
        return true;
    });

    // ========================================
    // Relativistic Velocity Addition Tests
    // ========================================

    run_test("Velocity addition at low speeds", []() {
        double v1 = 1000.0;  // 1 km/s
        double v2 = 2000.0;  // 2 km/s
        double v_total = velocityAddition(v1, v2);
        // At v << c, should be approximately v1 + v2
        ASSERT_NEAR(v_total, 3000.0, 1.0);  // Within 1 m/s
        return true;
    });

    run_test("Velocity addition at 0.5c + 0.5c", []() {
        double v1 = 0.5 * c;
        double v2 = 0.5 * c;
        double v_total = velocityAddition(v1, v2);
        // u = (0.5c + 0.5c)/(1 + 0.25) = c/1.25 = 0.8c
        ASSERT_NEAR(v_total, 0.8 * c, TOLERANCE);
        return true;
    });

    run_test("Velocity addition at 0.6c + 0.6c", []() {
        double v1 = 0.6 * c;
        double v2 = 0.6 * c;
        double v_total = velocityAddition(v1, v2);
        // u = (1.2c)/(1 + 0.36) = 1.2c/1.36 ≈ 0.882c
        double expected = 1.2 * c / 1.36;
        ASSERT_NEAR(v_total, expected, TOLERANCE);
        return true;
    });

    run_test("Velocity addition: c + any velocity = c", []() {
        double v1 = c;
        double v2 = 0.5 * c;
        double v_total = velocityAddition(v1, v2);
        ASSERT_NEAR(v_total, c, TOLERANCE);
        return true;
    });

    run_test("Velocity addition result always < c", []() {
        for (double v1 = 0.1*c; v1 < 0.99*c; v1 += 0.1*c) {
            for (double v2 = 0.1*c; v2 < 0.99*c; v2 += 0.1*c) {
                double v_total = velocityAddition(v1, v2);
                ASSERT_TRUE(v_total <= c);
            }
        }
        return true;
    });

    // ========================================
    // Relativistic Energy Tests
    // ========================================

    run_test("Rest energy E=mc²", []() {
        double mass = 1.0;  // 1 kg
        double E0 = restEnergy(mass);
        // E₀ = mc² = 1 kg * (3×10⁸)² = 9×10¹⁶ J
        double expected = mass * c * c;
        ASSERT_NEAR(E0, expected, TOLERANCE);
        return true;
    });

    run_test("Rest energy of electron", []() {
        double m_e = 9.10938e-31;  // kg
        double E0 = restEnergy(m_e);
        // E₀ ≈ 8.187×10⁻¹⁴ J ≈ 0.511 MeV
        double expected = m_e * c * c;
        ASSERT_NEAR(E0, expected, 1e-28);
        return true;
    });

    run_test("Relativistic energy at rest", []() {
        double mass = 1.0;
        double E = relativisticEnergy(mass, 0.0);
        double E0 = restEnergy(mass);
        ASSERT_NEAR(E, E0, TOLERANCE);
        return true;
    });

    run_test("Relativistic energy at 0.6c", []() {
        double mass = 1.0;
        double v = 0.6 * c;
        double E = relativisticEnergy(mass, v);
        // E = γmc² = 1.25 * mc²
        double E0 = restEnergy(mass);
        ASSERT_NEAR(E, 1.25 * E0, TOLERANCE);
        return true;
    });

    run_test("Relativistic energy at 0.8c", []() {
        double mass = 1.0;
        double v = 0.8 * c;
        double E = relativisticEnergy(mass, v);
        // E = γmc² = (5/3) * mc²
        double E0 = restEnergy(mass);
        ASSERT_NEAR(E, (5.0/3.0) * E0, 100.0);  // Large tolerance for c² energies
        return true;
    });

    run_test("Relativistic kinetic energy at rest", []() {
        double mass = 1.0;
        double KE = relativisticKineticEnergy(mass, 0.0);
        ASSERT_NEAR(KE, 0.0, TOLERANCE);
        return true;
    });

    run_test("Relativistic kinetic energy at low speed", []() {
        double mass = 1.0;
        double v = 1000.0;  // 1 km/s << c
        double KE_rel = relativisticKineticEnergy(mass, v);
        double KE_classical = 0.5 * mass * v * v;
        // At low speeds, relativistic ≈ classical
        ASSERT_NEAR(KE_rel, KE_classical, 50.0);  // Within 50 J
        return true;
    });

    run_test("Relativistic kinetic energy at 0.6c", []() {
        double mass = 1.0;
        double v = 0.6 * c;
        double KE = relativisticKineticEnergy(mass, v);
        // KE = (γ - 1)mc² = (1.25 - 1) * mc² = 0.25mc²
        double E0 = restEnergy(mass);
        ASSERT_NEAR(KE, 0.25 * E0, 100.0);  // Large tolerance for c² energies
        return true;
    });

    run_test("Total energy = rest energy + kinetic energy", []() {
        double mass = 1.0;
        double v = 0.8 * c;
        double E_total = relativisticEnergy(mass, v);
        double E0 = restEnergy(mass);
        double KE = relativisticKineticEnergy(mass, v);
        ASSERT_NEAR(E_total, E0 + KE, TOLERANCE);
        return true;
    });

    // ========================================
    // Relativistic Momentum Tests
    // ========================================

    run_test("Relativistic momentum at rest", []() {
        double mass = 1.0;
        double p = relativisticMomentum(mass, 0.0);
        ASSERT_NEAR(p, 0.0, TOLERANCE);
        return true;
    });

    run_test("Relativistic momentum at low speed", []() {
        double mass = 1.0;
        double v = 1000.0;  // 1 km/s << c
        double p_rel = relativisticMomentum(mass, v);
        double p_classical = mass * v;
        // At low speeds, relativistic ≈ classical
        ASSERT_NEAR(p_rel, p_classical, 0.01);
        return true;
    });

    run_test("Relativistic momentum at 0.6c", []() {
        double mass = 1.0;
        double v = 0.6 * c;
        double p = relativisticMomentum(mass, v);
        // p = γmv = 1.25 * m * 0.6c = 0.75mc
        double expected = 1.25 * mass * v;
        ASSERT_NEAR(p, expected, TOLERANCE);
        return true;
    });

    run_test("Relativistic momentum at 0.8c", []() {
        double mass = 1.0;
        double v = 0.8 * c;
        double p = relativisticMomentum(mass, v);
        // p = γmv = (5/3) * m * 0.8c = (4/3)mc
        double expected = (5.0/3.0) * mass * v;
        ASSERT_NEAR(p, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Energy-Momentum Relation Tests
    // ========================================

    run_test("Energy-momentum relation: E² = (pc)² + (mc²)²", []() {
        double mass = 1.0;
        double v = 0.6 * c;
        double E = relativisticEnergy(mass, v);
        double p = relativisticMomentum(mass, v);
        double E0 = restEnergy(mass);

        double E_squared = E * E;
        double relation = p*c*p*c + E0*E0;
        ASSERT_NEAR(E_squared, relation, 1e19);  // Very large tolerance for c⁴ terms
        return true;
    });

    run_test("Energy from momentum formula", []() {
        double mass = 1.0;
        double v = 0.8 * c;
        double p = relativisticMomentum(mass, v);
        double E = energyFromMomentum(p, mass);
        double E_direct = relativisticEnergy(mass, v);
        ASSERT_NEAR(E, E_direct, 1e10);
        return true;
    });

    run_test("Photon energy E=pc (massless particle)", []() {
        double p = 1e-27;  // kg⋅m/s
        double E = photonEnergy(p);
        double expected = p * c;
        ASSERT_NEAR(E, expected, 1e-19);
        return true;
    });

    run_test("Energy-momentum for massless particle", []() {
        double p = 1e-27;
        double E_photon = photonEnergy(p);
        double E_formula = energyFromMomentum(p, 0.0);  // m = 0
        ASSERT_NEAR(E_photon, E_formula, 1e-19);
        return true;
    });

    // ========================================
    // Relativistic Doppler Effect Tests
    // ========================================

    run_test("Doppler effect at rest", []() {
        double f0 = 1e14;  // 100 THz (optical frequency)
        double f = relativisticDoppler(f0, 0.0);
        ASSERT_NEAR(f, f0, TOLERANCE);
        return true;
    });

    run_test("Doppler redshift (source receding)", []() {
        double f0 = 1e14;  // 100 THz
        double v = 0.5 * c;  // Receding at 0.5c
        double f = relativisticDoppler(f0, -v);  // Negative for receding
        // f = f₀ * √[(1-β)/(1+β)] where β = 0.5
        // f = f₀ * √[0.5/1.5] = f₀ * √(1/3) ≈ 0.577f₀
        double expected = f0 * std::sqrt(0.5 / 1.5);
        ASSERT_NEAR(f, expected, f0 * 1e-6);
        return true;
    });

    run_test("Doppler blueshift (source approaching)", []() {
        double f0 = 1e14;
        double v = 0.5 * c;  // Approaching at 0.5c
        double f = relativisticDoppler(f0, v);  // Positive for approaching
        // f = f₀ * √[(1+β)/(1-β)] where β = 0.5
        // f = f₀ * √[1.5/0.5] = f₀ * √3 ≈ 1.732f₀
        double expected = f0 * std::sqrt(1.5 / 0.5);
        ASSERT_NEAR(f, expected, f0 * 1e-6);
        return true;
    });

    run_test("Redshift parameter z", []() {
        double v = 0.5 * c;
        double z = cosmologicalRedshift(v);
        // z = √[(1+β)/(1-β)] - 1 where β = 0.5
        double expected = std::sqrt(1.5 / 0.5) - 1.0;  // √3 - 1 ≈ 0.732
        ASSERT_NEAR(z, expected, 1e-6);
        return true;
    });

    run_test("Redshift z=0 at rest", []() {
        double z = cosmologicalRedshift(0.0);
        ASSERT_NEAR(z, 0.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Spacetime Interval Tests
    // ========================================

    run_test("Spacetime interval - timelike", []() {
        double dt = 10.0;  // 10 seconds
        double dx = 1e6;   // 1000 km (much less than c*dt)
        double s2 = spacetimeInterval(dt, dx);
        // s² = c²Δt² - Δx² > 0 for timelike
        ASSERT_TRUE(s2 > 0);
        return true;
    });

    run_test("Spacetime interval - spacelike", []() {
        double dt = 1e-6;  // 1 microsecond
        double dx = 1e6;   // 1000 km (much more than c*dt)
        double s2 = spacetimeInterval(dt, dx);
        // s² = c²Δt² - Δx² < 0 for spacelike
        ASSERT_TRUE(s2 < 0);
        return true;
    });

    run_test("Spacetime interval - lightlike", []() {
        double dt = 1.0;  // 1 second
        double dx = c;    // Light travels c meters in 1 second
        double s2 = spacetimeInterval(dt, dx);
        // s² = c²Δt² - Δx² = 0 for lightlike
        ASSERT_NEAR(s2, 0.0, 1e10);  // Large tolerance due to c²
        return true;
    });

    run_test("Spacetime interval invariance", []() {
        // In rest frame
        double dt = 5.0;
        double dx = 0.0;
        double s2_rest = spacetimeInterval(dt, dx);

        // In moving frame (proper time relates to interval)
        // For stationary object: s² = (cτ)²
        double expected = c * c * dt * dt;
        ASSERT_NEAR(s2_rest, expected, 1e10);
        return true;
    });

    // ========================================
    // Invariant Mass Tests
    // ========================================

    run_test("Invariant mass of particle at rest", []() {
        double mass = 1.0;
        double E = restEnergy(mass);
        double p = 0.0;
        double m_inv = invariantMass(E, p);
        ASSERT_NEAR(m_inv, mass, TOLERANCE);
        return true;
    });

    run_test("Invariant mass of moving particle", []() {
        double mass = 1.0;
        double v = 0.6 * c;
        double E = relativisticEnergy(mass, v);
        double p = relativisticMomentum(mass, v);
        double m_inv = invariantMass(E, p);
        // Invariant mass should equal rest mass
        ASSERT_NEAR(m_inv, mass, 1e-6);
        return true;
    });

    run_test("Invariant mass of photon is zero", []() {
        double p = 1e-27;
        double E = photonEnergy(p);
        double m_inv = invariantMass(E, p);
        ASSERT_NEAR(m_inv, 0.0, 1e-40);
        return true;
    });

    // ========================================
    // Rapidity Tests
    // ========================================

    run_test("Rapidity at rest", []() {
        double eta = rapidity(0.0);
        ASSERT_NEAR(eta, 0.0, TOLERANCE);
        return true;
    });

    run_test("Rapidity calculation", []() {
        double v = 0.6 * c;
        double eta = rapidity(v);
        // η = arctanh(v/c) = arctanh(0.6)
        double expected = std::atanh(0.6);
        ASSERT_NEAR(eta, expected, TOLERANCE);
        return true;
    });

    run_test("Rapidity to velocity conversion", []() {
        double eta = 1.0;
        double v = velocityFromRapidity(eta);
        // v = c * tanh(η) = c * tanh(1)
        double expected = c * std::tanh(1.0);
        ASSERT_NEAR(v, expected, TOLERANCE);
        return true;
    });

    run_test("Velocity to rapidity to velocity consistency", []() {
        double v_original = 0.5 * c;
        double eta = rapidity(v_original);
        double v_back = velocityFromRapidity(eta);
        ASSERT_NEAR(v_back, v_original, TOLERANCE);
        return true;
    });

    run_test("Rapidity addition is linear", []() {
        double v1 = 0.3 * c;
        double v2 = 0.4 * c;
        double eta1 = rapidity(v1);
        double eta2 = rapidity(v2);
        double eta_sum = eta1 + eta2;
        double v_sum = velocityFromRapidity(eta_sum);

        // Compare with velocity addition formula
        double v_add = velocityAddition(v1, v2);
        ASSERT_NEAR(v_sum, v_add, TOLERANCE);
        return true;
    });

    // ========================================
    // Consistency and Cross-Validation Tests
    // ========================================

    run_test("Time dilation and length contraction consistency", []() {
        // For same Lorentz factor: L/L₀ = τ/t
        double v = 0.75 * c;
        double gamma = lorentzFactor(v);

        double t = 10.0;
        double tau = properTime(t, v);
        double time_ratio = tau / t;

        double L0 = 100.0;
        double L = lengthContraction(L0, v);
        double length_ratio = L / L0;

        ASSERT_NEAR(time_ratio, length_ratio, TOLERANCE);
        ASSERT_NEAR(time_ratio, 1.0/gamma, TOLERANCE);
        return true;
    });

    run_test("E² - (pc)² = (mc²)² for various velocities", []() {
        double mass = 1.0;
        for (double v = 0.0; v < 0.9*c; v += 0.1*c) {
            double E = relativisticEnergy(mass, v);
            double p = relativisticMomentum(mass, v);
            double E0 = restEnergy(mass);

            double left = E*E - p*c*p*c;
            double right = E0 * E0;
            ASSERT_NEAR(left, right, 1e19);  // Very large tolerance for c⁴ terms
        }
        return true;
    });

    run_test("Velocity addition symmetry", []() {
        double v1 = 0.4 * c;
        double v2 = 0.5 * c;
        double v12 = velocityAddition(v1, v2);
        double v21 = velocityAddition(v2, v1);
        ASSERT_NEAR(v12, v21, TOLERANCE);
        return true;
    });

    run_test("Classical limit: Low velocity energy", []() {
        double mass = 1.0;
        double v = 100.0;  // 100 m/s << c
        double KE_rel = relativisticKineticEnergy(mass, v);
        double KE_class = 0.5 * mass * v * v;
        // At v=100 m/s, v²/c² ~ 10⁻¹³, numerical precision in c² calculations
        // can introduce small absolute differences (~8 J) even though relative error is tiny
        ASSERT_NEAR(KE_rel, KE_class, 10.0);  // Within 10 J (0.2% relative error)
        return true;
    });

    run_test("Classical limit: Low velocity momentum", []() {
        double mass = 1.0;
        double v = 100.0;  // 100 m/s << c
        double p_rel = relativisticMomentum(mass, v);
        double p_class = mass * v;
        double rel_error = std::abs(p_rel - p_class) / p_class;
        ASSERT_TRUE(rel_error < 1e-10);
        return true;
    });

    run_test("Ultra-relativistic limit: E ≈ pc for v→c", []() {
        double mass = 1e-30;  // Very small mass
        double v = 0.99999 * c;  // Nearly speed of light
        double E = relativisticEnergy(mass, v);
        double p = relativisticMomentum(mass, v);
        double pc = p * c;
        // E should be approximately equal to pc
        double rel_error = std::abs(E - pc) / E;
        ASSERT_TRUE(rel_error < 0.01);  // Within 1%
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Special Relativity" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "✓ All special relativity tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Lorentz transformations (γ, β)" << std::endl;
        std::cout << "  - Time dilation and length contraction" << std::endl;
        std::cout << "  - Relativistic velocity addition" << std::endl;
        std::cout << "  - Relativistic energy (E = γmc²)" << std::endl;
        std::cout << "  - Relativistic momentum (p = γmv)" << std::endl;
        std::cout << "  - Energy-momentum relation (E² = (pc)² + (mc²)²)" << std::endl;
        std::cout << "  - Relativistic Doppler effect" << std::endl;
        std::cout << "  - Spacetime intervals" << std::endl;
        std::cout << "  - Invariant mass" << std::endl;
        std::cout << "  - Rapidity and velocity transformations" << std::endl;
        std::cout << "  - Classical and ultra-relativistic limits" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
