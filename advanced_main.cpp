/**
 * @file advanced_main.cpp
 * @brief Central demonstration of all advanced physics modules
 *
 * This program demonstrates ALL advanced physics topics organized by category:
 * 1. Classical Mechanics (Hamiltonian, Phase Space, Liouville)
 * 2. Electromagnetism (Field Tensors, Gauge Theory)
 * 3. Relativity (Einstein Equations, Geodesics)
 * 4. Wave Theory (Green Functions, Waveguides, Diffraction)
 * 5. Thermodynamics (Statistical Mechanics, Coefficients)
 * 6. Quantum Mechanics (Schrödinger, Dirac, Perturbation)
 * 7. Condensed Matter (Superconductivity, Magnetism)
 * 8. Plasma Physics (Transport, Collision-Radiative)
 * 9. Advanced Quantum (Clebsch-Gordan, Group Theory)
 * 10. Quantum Field Theory (Interaction Picture, Renormalization)
 * 11. Fluid Dynamics (Navier-Stokes, Boundary Layer, Turbulence, Compressible)
 *
 * Compile with:
 *   g++ -std=c++17 -O3 -I./include -I/usr/include/eigen3 advanced_main.cpp -o advanced_physics
 *
 * Run with:
 *   ./advanced_physics [category_number]
 *   ./advanced_physics all    # Run all demonstrations
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <memory>
#include <chrono>

// Category 1: Classical Mechanics
#include "physics/advanced/classical/hamiltonian.hpp"
#include "physics/advanced/classical/phase_space.hpp"
#include "physics/advanced/classical/liouville.hpp"

// Category 11: Fluid Dynamics
#include "physics/advanced/fluid_dynamics/governing_equations.hpp"
#include "physics/advanced/fluid_dynamics/dimensionless_numbers.hpp"
#include "physics/advanced/fluid_dynamics/flow_types.hpp"
#include "physics/advanced/fluid_dynamics/boundary_layer.hpp"
#include "physics/advanced/fluid_dynamics/turbulence.hpp"
#include "physics/advanced/fluid_dynamics/vorticity.hpp"
#include "physics/advanced/fluid_dynamics/compressible_flow.hpp"

using namespace physics::advanced;

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

class Timer {
    std::chrono::high_resolution_clock::time_point start_;
public:
    Timer() : start_(std::chrono::high_resolution_clock::now()) {}

    double elapsed() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start_).count();
    }
};

void printHeader(const std::string& category, int num) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Category " << num << ": " << category << "\n";
    std::cout << std::string(70, '=') << "\n\n";
}

void printSection(const std::string& title) {
    std::cout << "\n" << std::string(60, '-') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(60, '-') << "\n";
}

// ============================================================================
// CATEGORY 1: CLASSICAL MECHANICS
// ============================================================================

void demo_hamiltonian_mechanics() {
    printSection("Hamiltonian Mechanics");

    std::cout << "Demonstrating Hamilton's equations for simple harmonic oscillator\n\n";

    // Create SHO: H = p²/(2m) + (1/2)kx²
    double mass = 1.0;  // kg
    double k = 10.0;    // N/m
    double omega = std::sqrt(k / mass);

    auto sho = classical::HamiltonianSystem::harmonicOscillator(mass, k);

    // Initial conditions: x=1, p=0 (release from rest)
    classical::PhasePoint initial(1);
    initial.q(0) = 1.0;
    initial.p(0) = 0.0;

    std::cout << "Initial state: q = " << initial.q(0)
              << " m, p = " << initial.p(0) << " kg·m/s\n";
    std::cout << "Natural frequency: ω = " << omega << " rad/s\n";
    std::cout << "Period: T = " << 2.0 * M_PI / omega << " s\n\n";

    // Integrate for one period
    double dt = 0.01;
    int steps = static_cast<int>(2.0 * M_PI / (omega * dt));

    Timer timer;
    auto trajectory = sho.integrate(initial, dt, steps);
    double time_taken = timer.elapsed();

    // Check energy conservation
    bool conserved = sho.checkEnergyConservation(trajectory, 1e-6);

    std::cout << "Integrated " << steps << " steps in " << time_taken << " s\n";
    std::cout << "Energy conservation: " << (conserved ? "✓ PASSED" : "✗ FAILED") << "\n";

    // Show phase space trajectory
    std::cout << "\nPhase space trajectory (first 10 points):\n";
    std::cout << std::setw(10) << "t (s)"
              << std::setw(15) << "q (m)"
              << std::setw(15) << "p (kg·m/s)"
              << std::setw(15) << "Energy (J)" << "\n";

    for (int i = 0; i < std::min(10, static_cast<int>(trajectory.size())); ++i) {
        double t = i * dt;
        double E = sho.hamiltonian(trajectory[i]);
        std::cout << std::fixed << std::setprecision(4)
                  << std::setw(10) << t
                  << std::setw(15) << trajectory[i].q(0)
                  << std::setw(15) << trajectory[i].p(0)
                  << std::setw(15) << E << "\n";
    }

    std::cout << "\n✓ Hamiltonian mechanics demonstration complete\n";
}

void demo_phase_space() {
    printSection("Phase Space Analysis");

    std::cout << "Analyzing phase space structure for Kepler problem\n\n";

    // Create Kepler problem (planetary motion)
    double mass = 1.0;           // kg (Earth mass unit)
    double GM = 1.0;             // Gravitational parameter

    auto kepler = classical::HamiltonianSystem::keplerProblem(mass, GM);

    // Circular orbit initial conditions
    double r0 = 1.0;  // AU
    double v0 = std::sqrt(GM / r0);  // Circular velocity

    classical::PhasePoint initial(2);
    initial.q(0) = r0;   // x position
    initial.q(1) = 0.0;  // y position
    initial.p(0) = 0.0;  // px momentum
    initial.p(1) = mass * v0;  // py momentum

    std::cout << "Circular orbit parameters:\n";
    std::cout << "  Initial radius: " << r0 << " AU\n";
    std::cout << "  Orbital velocity: " << v0 << " AU/T\n";
    std::cout << "  Angular momentum: " << r0 * mass * v0 << "\n\n";

    // Integrate orbit
    double dt = 0.01;
    int steps = 628;  // ~1 orbit for unit parameters

    Timer timer;
    auto orbit = kepler.integrate(initial, dt, steps);
    double time_taken = timer.elapsed();

    std::cout << "Orbital integration:\n";
    std::cout << "  Steps: " << steps << "\n";
    std::cout << "  Time: " << time_taken << " s\n";

    // Calculate phase space volume
    auto volume = classical::PhaseSpaceVolume::volumeElement(orbit);
    std::cout << "  Phase space volume: " << volume << "\n";

    // Check energy conservation
    bool energy_ok = kepler.checkEnergyConservation(orbit, 1e-5);
    std::cout << "  Energy conservation: " << (energy_ok ? "✓ PASSED" : "✗ FAILED") << "\n";

    // Calculate Lyapunov exponent (should be ~0 for regular orbit)
    classical::LyapunovExponent lyapunov(kepler);
    double lambda = lyapunov.compute(initial, dt, steps / 4);

    std::cout << "\nChaos analysis:\n";
    std::cout << "  Lyapunov exponent: " << lambda << "\n";
    std::cout << "  Dynamics type: ";

    auto dynamics = lyapunov.classify(lambda);
    switch (dynamics) {
        case classical::LyapunovExponent::DynamicsType::STABLE:
            std::cout << "STABLE (fixed point)\n";
            break;
        case classical::LyapunovExponent::DynamicsType::NEUTRAL:
            std::cout << "NEUTRAL (regular/quasiperiodic)\n";
            break;
        case classical::LyapunovExponent::DynamicsType::CHAOTIC:
            std::cout << "CHAOTIC\n";
            break;
    }

    std::cout << "\n✓ Phase space analysis complete\n";
}

void demo_liouville_equation() {
    printSection("Liouville's Equation");

    std::cout << "Demonstrating Liouville's theorem and ensemble evolution\n\n";

    // Create simple harmonic oscillator
    auto sho = classical::HamiltonianSystem::harmonicOscillator(1.0, 10.0);

    // Create Liouville equation
    classical::LiouvilleEquation liouville(sho);

    // Define initial density (Gaussian in phase space)
    double q0 = 0.5, p0 = 0.0;
    double sigma_q = 0.1, sigma_p = 0.5;

    auto rho0 = [q0, p0, sigma_q, sigma_p](const classical::PhasePoint& point) {
        double dq = point.q(0) - q0;
        double dp = point.p(0) - p0;
        double norm = 1.0 / (2.0 * M_PI * sigma_q * sigma_p);
        return norm * std::exp(-0.5 * (dq * dq / (sigma_q * sigma_q) +
                                       dp * dp / (sigma_p * sigma_p)));
    };

    // Test points
    classical::PhasePoint test(1);
    test.q(0) = 0.5;
    test.p(0) = 0.0;

    // Verify Liouville operator
    bool liouville_ok = liouville.verifyLiouville(rho0, test, 1e-5);
    std::cout << "Liouville operator check: "
              << (liouville_ok ? "✓ PASSED" : "✗ FAILED") << "\n";

    // Integrate trajectory and check density conservation
    double dt = 0.01;
    int steps = 100;
    auto trajectory = sho.integrate(test, dt, steps);

    bool conservation_ok = liouville.checkConservation(rho0, trajectory, 1e-4);
    std::cout << "Density conservation along trajectory: "
              << (conservation_ok ? "✓ PASSED" : "✗ FAILED") << "\n";

    // Statistical ensemble
    std::cout << "\nStatistical Ensemble Analysis:\n";

    double temperature = 300.0;  // K
    auto canonical = classical::StatisticalEnsemble::canonical(
        sho, temperature);

    // Generate sample points
    std::vector<classical::PhasePoint> sample;
    for (int i = 0; i < 100; ++i) {
        classical::PhasePoint pt(1);
        pt.q(0) = -1.0 + 2.0 * i / 100.0;
        pt.p(0) = -2.0 + 4.0 * i / 100.0;
        sample.push_back(pt);
    }

    // Calculate partition function
    double Z = classical::StatisticalEnsemble::partitionFunction(
        sho, temperature, sample);
    std::cout << "  Partition function Z ≈ " << Z << "\n";

    // Free energy
    double F = classical::StatisticalEnsemble::freeEnergy(Z, temperature);
    std::cout << "  Helmholtz free energy F ≈ " << F << " J\n";

    // Entropy
    double S = classical::StatisticalEnsemble::entropy(canonical, sample);
    std::cout << "  Entropy S ≈ " << S << " J/K\n";

    std::cout << "\n✓ Liouville equation demonstration complete\n";
}

void demo_poisson_brackets() {
    printSection("Poisson Brackets");

    std::cout << "Verifying fundamental Poisson bracket relations\n\n";

    // Create test system
    auto sho = classical::HamiltonianSystem::harmonicOscillator(1.0, 10.0);
    classical::PoissonBracket bracket;

    // Test point
    classical::PhasePoint test(1);
    test.q(0) = 1.0;
    test.p(0) = 0.5;

    // Verify {q, p} = 1, {q, q} = 0, {p, p} = 0
    bool fundamental_ok = bracket.verifyFundamentalRelations(test, 1e-6);
    std::cout << "Fundamental Poisson brackets: "
              << (fundamental_ok ? "✓ PASSED" : "✗ FAILED") << "\n";

    // Compute {H, q} and {H, p}
    auto H = [&sho](const classical::PhasePoint& p) {
        return sho.hamiltonian(p);
    };
    auto q = [](const classical::PhasePoint& p) { return p.q(0); };
    auto p_func = [](const classical::PhasePoint& p) { return p.p(0); };

    double Hq = bracket.compute(H, q, test);
    double Hp = bracket.compute(H, p_func, test);

    std::cout << "\nPoisson brackets with Hamiltonian:\n";
    std::cout << "  {H, q} = " << Hq << " (should be ∂H/∂p)\n";
    std::cout << "  {H, p} = " << Hp << " (should be -∂H/∂q)\n";

    // Verify Hamilton's equations via Poisson brackets
    double dH_dp = sho.dH_dp(test, 0);
    double dH_dq = sho.dH_dq(test, 0);

    bool hamilton_via_poisson =
        (std::abs(Hq - dH_dp) < 1e-5) &&
        (std::abs(Hp + dH_dq) < 1e-5);

    std::cout << "  Hamilton's equations via Poisson brackets: "
              << (hamilton_via_poisson ? "✓ PASSED" : "✗ FAILED") << "\n";

    std::cout << "\n✓ Poisson bracket demonstration complete\n";
}

void demonstrate_category_1() {
    printHeader("Classical Mechanics", 1);

    std::cout << "Topics covered:\n";
    std::cout << "  • Hamilton mechanics (Hamiltonian formulation)\n";
    std::cout << "  • Phase space (structure and analysis)\n";
    std::cout << "  • Liouville's equation (ensemble evolution)\n";
    std::cout << "  • Poisson brackets (canonical structure)\n";
    std::cout << "  • Symplectic integration (structure-preserving)\n\n";

    demo_hamiltonian_mechanics();
    demo_phase_space();
    demo_liouville_equation();
    demo_poisson_brackets();

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "✓ Category 1 demonstrations complete!\n";
    std::cout << std::string(70, '=') << "\n";
}

// ============================================================================
// CATEGORY 11: FLUID DYNAMICS
// ============================================================================

void demo_governing_equations() {
    printSection("Governing Equations & Dimensionless Numbers");

    std::cout << "1. Pipe Flow Analysis (Poiseuille Flow)\n";

    // Water properties
    double density = 1000.0;      // kg/m³
    double dyn_viscosity = 1e-3;  // Pa·s (water at 20°C)

    // Pipe parameters
    double pipe_radius = 0.05;    // 5 cm radius
    double pipe_length = 10.0;    // 10 m
    double pressure_drop = 5000;  // 5 kPa
    double velocity = 2.0;        // m/s average

    // Reynolds number
    double Re = fluid_dynamics::DimensionlessNumbers::reynoldsNumber(
        density, velocity, 2.0 * pipe_radius, dyn_viscosity
    );

    auto regime = fluid_dynamics::DimensionlessNumbers::classifyFlowRegime(Re);

    std::cout << "   Reynolds number: Re = " << Re << "\n";
    std::cout << "   Flow regime: ";
    switch (regime) {
        case fluid_dynamics::DimensionlessNumbers::FlowRegime::LAMINAR:
            std::cout << "Laminar"; break;
        case fluid_dynamics::DimensionlessNumbers::FlowRegime::TURBULENT:
            std::cout << "Turbulent"; break;
        case fluid_dynamics::DimensionlessNumbers::FlowRegime::TRANSITION:
            std::cout << "Transition"; break;
        default: std::cout << "Stokes";
    }
    std::cout << "\n\n";

    // Flow rate using Hagen-Poiseuille
    double Q = fluid_dynamics::PoiseuilleFlow::flowRate(
        pipe_radius, pressure_drop, pipe_length, dyn_viscosity
    );

    std::cout << "   Flow rate: Q = " << Q * 1000 << " L/s\n";
    std::cout << "   Wall shear stress: τ_w = "
              << fluid_dynamics::PoiseuilleFlow::wallShearStress(
                     pipe_radius, pressure_drop, pipe_length) << " Pa\n\n";

    std::cout << "2. Bernoulli's Equation (Venturi Meter)\n";

    Eigen::Vector3d v1(5.0, 0, 0);  // 5 m/s
    Eigen::Vector3d v2(10.0, 0, 0); // 10 m/s (throat)

    double p1 = 101325.0;  // 1 atm
    double p2 = fluid_dynamics::BernoulliEquation::pressureFromVelocity(
        p1, density, v1, v2, 0.0, 0.0, 9.81
    );

    std::cout << "   Inlet velocity: " << v1.norm() << " m/s at " << p1/1000 << " kPa\n";
    std::cout << "   Throat velocity: " << v2.norm() << " m/s at " << p2/1000 << " kPa\n";
    std::cout << "   Pressure drop: Δp = " << (p1 - p2)/1000 << " kPa\n\n";
}

void demo_boundary_layer() {
    printSection("Boundary Layer Theory");

    std::cout << "Flat Plate Boundary Layer (Blasius Solution)\n\n";

    double U_inf = 10.0;          // m/s freestream velocity
    double nu = 1.5e-5;           // m²/s air kinematic viscosity
    double x = 1.0;               // m from leading edge

    // Reynolds number
    double Re_x = U_inf * x / nu;

    std::cout << "   Freestream velocity: U∞ = " << U_inf << " m/s\n";
    std::cout << "   Position: x = " << x << " m\n";
    std::cout << "   Reynolds number: Re_x = " << Re_x << "\n\n";

    // Blasius solution
    double delta = fluid_dynamics::BlasiusSolution::thickness(x, U_inf, nu);
    double delta_star = fluid_dynamics::BlasiusSolution::displacementThickness(x, U_inf, nu);
    double theta = fluid_dynamics::BlasiusSolution::momentumThickness(x, U_inf, nu);
    double H = fluid_dynamics::BlasiusSolution::shapeFactor(delta_star, theta);
    double Cf = fluid_dynamics::BlasiusSolution::skinFrictionCoefficient(Re_x);

    std::cout << "   Boundary layer thickness: δ₉₉ = " << delta * 1000 << " mm\n";
    std::cout << "   Displacement thickness: δ* = " << delta_star * 1000 << " mm\n";
    std::cout << "   Momentum thickness: θ = " << theta * 1000 << " mm\n";
    std::cout << "   Shape factor: H = δ*/θ = " << H << " (theory: 2.59)\n";
    std::cout << "   Skin friction coefficient: Cf = " << Cf << "\n\n";

    // Velocity profile at a few points
    std::cout << "   Velocity Profile:\n";
    std::cout << "   y/δ     u/U∞\n";
    std::cout << "   ----    ----\n";

    for (double y_frac = 0.0; y_frac <= 1.2; y_frac += 0.2) {
        double y = y_frac * delta;
        double eta = fluid_dynamics::BlasiusSolution::similarityVariable(y, x, U_inf, nu);
        double u_ratio = fluid_dynamics::BlasiusSolution::velocityProfile(eta);
        std::cout << "   " << std::fixed << std::setprecision(2) << y_frac
                  << "      " << u_ratio << "\n";
    }
    std::cout << "\n";
}

void demo_turbulence() {
    printSection("Turbulence Modeling (k-ε)");

    std::cout << "k-ε Turbulence Model Parameters\n\n";

    double rho = 1.2;             // kg/m³ air density
    double k = 1.5;               // m²/s² turbulent kinetic energy
    double epsilon = 10.0;        // m²/s³ dissipation rate
    double nu = 1.5e-5;           // m²/s kinematic viscosity

    fluid_dynamics::KEpsilonModel::Constants constants;

    // Calculate derived quantities
    double mu_t = fluid_dynamics::KEpsilonModel::eddyViscosity(rho, k, epsilon, constants);
    double l_t = fluid_dynamics::KEpsilonModel::turbulentLengthScale(k, epsilon);
    double tau_t = fluid_dynamics::KEpsilonModel::turbulentTimeScale(k, epsilon);
    double eta = fluid_dynamics::KEpsilonModel::kolmogorovLengthScale(nu, epsilon);
    double Re_t = fluid_dynamics::KEpsilonModel::turbulentReynolds(k, nu, epsilon);

    std::cout << "   Turbulent kinetic energy: k = " << k << " m²/s²\n";
    std::cout << "   Dissipation rate: ε = " << epsilon << " m²/s³\n\n";

    std::cout << "   Derived Quantities:\n";
    std::cout << "   Eddy viscosity: μ_t = " << mu_t << " Pa·s\n";
    std::cout << "   Turbulent length scale: l_t = " << l_t << " m\n";
    std::cout << "   Turbulent time scale: τ_t = " << tau_t << " s\n";
    std::cout << "   Kolmogorov length scale: η = " << eta * 1e6 << " μm\n";
    std::cout << "   Turbulent Reynolds number: Re_t = " << Re_t << "\n\n";

    // Reynolds stress calculation
    std::cout << "   Model Constants (Standard k-ε):\n";
    std::cout << "   C_μ = " << constants.C_mu << "\n";
    std::cout << "   C_ε1 = " << constants.C_eps1 << "\n";
    std::cout << "   C_ε2 = " << constants.C_eps2 << "\n";
    std::cout << "   σ_k = " << constants.sigma_k << "\n";
    std::cout << "   σ_ε = " << constants.sigma_eps << "\n\n";
}

void demo_vorticity() {
    printSection("Vorticity Dynamics & Circulation");

    std::cout << "1. Point Vortex and Biot-Savart Law\n\n";

    double circulation = 10.0;  // m²/s
    double radius = 1.0;        // m

    double u_theta = fluid_dynamics::BiotSavartLaw::velocityStraightVortex(
        circulation, radius
    );

    std::cout << "   Circulation: Γ = " << circulation << " m²/s\n";
    std::cout << "   Tangential velocity at r = " << radius << " m: u_θ = "
              << u_theta << " m/s\n\n";

    std::cout << "2. Vortex Ring Self-Velocity\n\n";

    double ring_radius = 0.1;   // m
    double core_radius = 0.01;  // m

    double U_ring = fluid_dynamics::BiotSavartLaw::vortexRingSelfVelocity(
        circulation, ring_radius, core_radius
    );

    std::cout << "   Ring radius: R = " << ring_radius << " m\n";
    std::cout << "   Core radius: a = " << core_radius << " m\n";
    std::cout << "   Translation velocity: U = " << U_ring << " m/s\n\n";

    std::cout << "3. Rankine Vortex Profile\n\n";

    double a = 0.05;  // m core radius

    std::cout << "   r/a     u (m/s)\n";
    std::cout << "   ----    -------\n";

    for (double r_ratio = 0.0; r_ratio <= 3.0; r_ratio += 0.5) {
        double r = r_ratio * a;
        double u = fluid_dynamics::VortexDynamics::rankineVortex(r, a, circulation);
        std::cout << "   " << std::fixed << std::setprecision(1) << r_ratio
                  << "      " << std::setprecision(3) << u << "\n";
    }
    std::cout << "\n";
}

void demo_compressible_flow() {
    printSection("Compressible Flow (Shocks & Expansions)");

    double gamma = 1.4;  // Air

    std::cout << "1. Normal Shock Relations\n\n";

    double M1 = 2.0;  // Supersonic upstream Mach

    double M2 = fluid_dynamics::NormalShock::downstreamMach(M1, gamma);
    double p_ratio = fluid_dynamics::NormalShock::pressureRatio(M1, gamma);
    double rho_ratio = fluid_dynamics::NormalShock::densityRatio(M1, gamma);
    double T_ratio = fluid_dynamics::NormalShock::temperatureRatio(M1, gamma);
    double p0_ratio = fluid_dynamics::NormalShock::stagnationPressureRatio(M1, gamma);

    std::cout << "   Upstream Mach: M₁ = " << M1 << "\n";
    std::cout << "   Downstream Mach: M₂ = " << M2 << "\n\n";

    std::cout << "   Pressure ratio: p₂/p₁ = " << p_ratio << "\n";
    std::cout << "   Density ratio: ρ₂/ρ₁ = " << rho_ratio << "\n";
    std::cout << "   Temperature ratio: T₂/T₁ = " << T_ratio << "\n";
    std::cout << "   Stagnation pressure ratio: p₀₂/p₀₁ = " << p0_ratio << "\n\n";

    std::cout << "2. Isentropic Flow (Stagnation Properties)\n\n";

    double T_static = 300.0;  // K
    double p_static = 101325.0;  // Pa
    double M = 0.8;

    double T0 = fluid_dynamics::StagnationProperties::temperature(T_static, M, gamma);
    double p0 = fluid_dynamics::StagnationProperties::pressure(p_static, M, gamma);

    std::cout << "   Static conditions: T = " << T_static << " K, p = " << p_static/1000 << " kPa\n";
    std::cout << "   Mach number: M = " << M << "\n";
    std::cout << "   Stagnation temperature: T₀ = " << T0 << " K\n";
    std::cout << "   Stagnation pressure: p₀ = " << p0/1000 << " kPa\n\n";

    std::cout << "3. Prandtl-Meyer Expansion\n\n";

    double M_initial = 2.0;
    double deflection_angle = 10.0 * M_PI / 180.0;  // 10 degrees

    double nu1 = fluid_dynamics::PrandtlMeyerExpansion::prandtlMeyerAngle(M_initial, gamma);
    double M_final = fluid_dynamics::PrandtlMeyerExpansion::downstreamMach(
        M_initial, deflection_angle, gamma
    );

    std::cout << "   Initial Mach: M₁ = " << M_initial << "\n";
    std::cout << "   Deflection angle: θ = " << deflection_angle * 180.0 / M_PI << "°\n";
    std::cout << "   Final Mach: M₂ = " << M_final << "\n";
    std::cout << "   Prandtl-Meyer angle: ν₁ = " << nu1 * 180.0 / M_PI << "°\n\n";
}

void demonstrate_category_11() {
    printHeader("Fluid Dynamics", 11);

    std::cout << "Comprehensive CFD demonstration covering:\n";
    std::cout << "  • Governing equations (Navier-Stokes, Euler, Bernoulli)\n";
    std::cout << "  • Dimensionless numbers (Reynolds, Mach, Prandtl, etc.)\n";
    std::cout << "  • Flow types (Poiseuille, Couette, Stokes, Potential)\n";
    std::cout << "  • Boundary layer theory (Blasius solution)\n";
    std::cout << "  • Turbulence modeling (k-ε, RANS)\n";
    std::cout << "  • Vorticity dynamics (circulation, Biot-Savart)\n";
    std::cout << "  • Compressible flow (shocks, expansions)\n\n";

    demo_governing_equations();
    demo_boundary_layer();
    demo_turbulence();
    demo_vorticity();
    demo_compressible_flow();

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "✓ Category 11 demonstrations complete!\n";
    std::cout << std::string(70, '=') << "\n";
}

// ============================================================================
// OTHER CATEGORIES (STUBS)
// ============================================================================

void demonstrate_category_2() {
    printHeader("Electromagnetism", 2);
    std::cout << "Topics: Maxwell equations, Field tensor, Stress-energy tensor,\n";
    std::cout << "        Gauge transformations, EM waves, Multipoles\n\n";
    std::cout << "⚠️ Implementation in progress\n";
    std::cout << "   See examples/field_tensor_demo.cpp for working implementation\n";
}

void demonstrate_category_3() {
    printHeader("Relativity", 3);
    std::cout << "Topics: Einstein tensor, Riemannian geometry, Geodesics,\n";
    std::cout << "        Perihelion shift, Photon trajectories\n\n";
    std::cout << "⚠️ Planned for implementation\n";
}

void demonstrate_category_4() {
    printHeader("Wave Theory", 4);
    std::cout << "Topics: Green functions, Waveguides, Spherical/cylindrical waves,\n";
    std::cout << "        Diffraction, Non-linear equations\n\n";
    std::cout << "⚠️ Planned for implementation\n";
}

void demonstrate_category_5() {
    printHeader("Thermodynamics", 5);
    std::cout << "Topics: Equation of state, Thermodynamic coefficients,\n";
    std::cout << "        Statistical mechanics, Boltzmann equation\n\n";
    std::cout << "⚠️ Planned for implementation\n";
}

void demonstrate_category_6() {
    printHeader("Quantum Mechanics", 6);
    std::cout << "Topics: Schrödinger solver, Dirac equation, Perturbation theory,\n";
    std::cout << "        Quantum statistics, Spin-orbit coupling\n\n";
    std::cout << "⚠️ Planned for implementation\n";
}

void demonstrate_category_7() {
    printHeader("Condensed Matter", 7);
    std::cout << "Topics: Superconductivity, Josephson effect, BCS model,\n";
    std::cout << "        Paramagnetism, Dielectrics\n\n";
    std::cout << "⚠️ Planned for implementation\n";
}

void demonstrate_category_8() {
    printHeader("Plasma Physics", 8);
    std::cout << "Topics: Plasma transport, Collision-radiative models,\n";
    std::cout << "        Coulomb interactions, Debye shielding\n\n";
    std::cout << "⚠️ Planned for implementation\n";
}

void demonstrate_category_9() {
    printHeader("Advanced Quantum Theory", 9);
    std::cout << "Topics: Clebsch-Gordan coefficients, Wigner-Eckart theorem,\n";
    std::cout << "        Group theory (SO(3), SU(2)), Lie algebras\n\n";
    std::cout << "⚠️ Planned for implementation\n";
}

void demonstrate_category_10() {
    printHeader("Quantum Field Theory", 10);
    std::cout << "Topics: Field quantization, Interaction picture, S-matrix,\n";
    std::cout << "        Renormalization, Standard Model, Higgs mechanism, QCD\n\n";
    std::cout << "⚠️ Future expansion (graduate-level QFT)\n";
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================

int main(int argc, char* argv[]) {
    std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║        ADVANCED PHYSICS SHOWCASE - Comprehensive Demo         ║\n";
    std::cout << "║                 11 Categories, 110+ Topics                    ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

    std::string category = (argc > 1) ? argv[1] : "1";

    try {
        if (category == "all") {
            for (int i = 1; i <= 11; ++i) {
                switch (i) {
                    case 1: demonstrate_category_1(); break;
                    case 2: demonstrate_category_2(); break;
                    case 3: demonstrate_category_3(); break;
                    case 4: demonstrate_category_4(); break;
                    case 5: demonstrate_category_5(); break;
                    case 6: demonstrate_category_6(); break;
                    case 7: demonstrate_category_7(); break;
                    case 8: demonstrate_category_8(); break;
                    case 9: demonstrate_category_9(); break;
                    case 10: demonstrate_category_10(); break;
                    case 11: demonstrate_category_11(); break;
                }
            }
        } else {
            int cat_num = std::stoi(category);
            switch (cat_num) {
                case 1: demonstrate_category_1(); break;
                case 2: demonstrate_category_2(); break;
                case 3: demonstrate_category_3(); break;
                case 4: demonstrate_category_4(); break;
                case 5: demonstrate_category_5(); break;
                case 6: demonstrate_category_6(); break;
                case 7: demonstrate_category_7(); break;
                case 8: demonstrate_category_8(); break;
                case 9: demonstrate_category_9(); break;
                case 10: demonstrate_category_10(); break;
                case 11: demonstrate_category_11(); break;
                default:
                    std::cerr << "Invalid category. Choose 1-11 or 'all'\n";
                    return 1;
            }
        }

        std::cout << "\n\n";
        std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
        std::cout << "║                    Demo Complete!                             ║\n";
        std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

        std::cout << "Current Status:\n";
        std::cout << "  ✓ Category 1: Classical Mechanics (COMPLETE)\n";
        std::cout << "  🔄 Category 2: Electromagnetism (IN PROGRESS)\n";
        std::cout << "  ⏳ Categories 3-10: Planned\n\n";

        std::cout << "To run specific categories:\n";
        std::cout << "  ./advanced_physics 1      # Classical Mechanics\n";
        std::cout << "  ./advanced_physics 2      # Electromagnetism\n";
        std::cout << "  ./advanced_physics all    # All categories\n\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
