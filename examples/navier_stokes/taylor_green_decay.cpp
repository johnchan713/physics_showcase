/**
 * @file taylor_green_decay.cpp
 * @brief Taylor-Green Vortex Energy Decay Simulation
 *
 * DESCRIPTION:
 * Simulates the classic Taylor-Green vortex, a well-studied test case for
 * incompressible flow solvers. The vortex undergoes transition to turbulence
 * and eventual decay due to viscous dissipation.
 *
 * INITIAL CONDITION:
 *   u(x,y,z,0) = sin(x)cos(y)cos(z)
 *   v(x,y,z,0) = -cos(x)sin(y)cos(z)
 *   w(x,y,z,0) = 0
 *
 * PHYSICS:
 * - Kinetic energy decays as: E(t) ≈ E₀ exp(-2νk²t)
 * - Enstrophy: Ω(t) initially grows, then decays
 * - Remains smooth (no singularities) - reference solution for validation
 *
 * USAGE:
 *   ./taylor_green_decay [N] [Re] [T_final]
 *
 *   N       : Grid resolution (default: 64)
 *   Re      : Reynolds number (default: 100)
 *   T_final : Final simulation time (default: 10.0)
 *
 * OUTPUT:
 *   taylor_green_timeseries.dat - Energy, enstrophy evolution
 *   Console: Real-time monitoring
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "../../new_theory/navier_stokes_solver.hpp"

using namespace new_theory::navier_stokes;

void print_header() {
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║      TAYLOR-GREEN VORTEX ENERGY DECAY SIMULATION        ║\n";
    std::cout << "║                 Millennium Prize Toolkit                ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
}

void print_simulation_info(int N, double Re, double nu, double T_final) {
    std::cout << "SIMULATION PARAMETERS:\n";
    std::cout << "  Grid resolution:  " << N << "³ = " << (N*N*N) << " points\n";
    std::cout << "  Domain:           [0, 2π]³\n";
    std::cout << "  Reynolds number:  Re = " << Re << "\n";
    std::cout << "  Kinematic viscosity: ν = " << nu << "\n";
    std::cout << "  Final time:       T = " << T_final << "\n";
    std::cout << "\n";

    std::cout << "EXPECTED BEHAVIOR:\n";
    std::cout << "  • Energy decay:    E(t) ≈ E₀ exp(-2νk²t)\n";
    std::cout << "  • Enstrophy:       Ω(t) initially grows, then decays\n";
    std::cout << "  • Regularity:      Solution remains smooth (no blow-up)\n";
    std::cout << "\n";
}

void print_progress_header() {
    std::cout << "SIMULATION PROGRESS:\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << std::setw(10) << "Time"
              << std::setw(15) << "Energy"
              << std::setw(15) << "Enstrophy"
              << std::setw(15) << "||ω||_∞"
              << std::setw(12) << "dt"
              << std::setw(13) << "Status" << "\n";
    std::cout << std::string(80, '-') << "\n";
}

void print_progress(double t, double E, double Omega, double omega_max,
                   double dt, const std::string& status) {
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << t
              << std::setw(15) << std::scientific << std::setprecision(4) << E
              << std::setw(15) << Omega
              << std::setw(15) << omega_max
              << std::setw(12) << std::fixed << std::setprecision(5) << dt
              << std::setw(13) << status << "\n";
}

int main(int argc, char** argv) {
    // Parse command line arguments
    int N = 64;
    double Re = 100.0;
    double T_final = 10.0;

    if (argc > 1) N = std::atoi(argv[1]);
    if (argc > 2) Re = std::atof(argv[2]);
    if (argc > 3) T_final = std::atof(argv[3]);

    // Compute viscosity from Reynolds number
    // Re = UL/ν, with U = 1, L = 1 (characteristic velocity and length)
    double nu = 1.0 / Re;

    print_header();
    print_simulation_info(N, Re, nu, T_final);

    // Create grid and solver
    std::cout << "Initializing solver...\n";
    SpectralGrid3D grid(N, 2.0 * M_PI);
    NavierStokesSolver solver(grid, nu);
    std::cout << "\n";

    // Set Taylor-Green initial condition
    std::cout << "Setting Taylor-Green vortex initial condition...\n";
    solver.setInitialConditionTaylorGreen();

    // Initial diagnostics
    double E0 = solver.computeEnergy();
    double Omega0 = solver.computeEnstrophy();
    std::cout << "Initial state:\n";
    std::cout << "  E₀ = " << E0 << "\n";
    std::cout << "  Ω₀ = " << Omega0 << "\n";
    std::cout << "\n";

    // Simulation parameters
    double dt_output = 0.1;  // Output every 0.1 time units
    double dt_max = 0.01;    // Maximum time step
    double CFL = 0.5;        // CFL number for adaptive stepping

    print_progress_header();

    // Time integration loop
    double t_next_output = 0.0;
    int step = 0;

    while (solver.getTime() < T_final) {
        // Adaptive time step
        double dt = solver.computeCFLTimeStep(CFL);
        dt = std::min(dt, dt_max);
        dt = std::min(dt, T_final - solver.getTime());  // Don't overshoot

        // Advance one step
        solver.stepRK4(dt);
        solver.updateRegularityMonitoring();
        step++;

        // Output at regular intervals
        if (solver.getTime() >= t_next_output || solver.getTime() >= T_final - 1e-10) {
            double E = solver.computeEnergy();
            double Omega = solver.computeEnstrophy();
            double omega_max = solver.computeVorticityMax();

            auto status = solver.checkRegularity();
            std::string status_str = status.is_regular ? "REGULAR" : "BLOW-UP!";

            print_progress(solver.getTime(), E, Omega, omega_max, dt, status_str);

            // Check for blow-up
            if (!status.is_regular) {
                std::cout << "\n";
                std::cout << "╔══════════════════════════════════════════════════════════╗\n";
                std::cout << "║                  ⚠️  BLOW-UP DETECTED! ⚠️                 ║\n";
                std::cout << "╚══════════════════════════════════════════════════════════╝\n";
                std::cout << "Criterion violated: " << status.criterion_violated << "\n";
                std::cout << "Severity score: " << status.severity_score << "\n";
                std::cout << "\n";
                std::cout << "NOTE: This should NOT happen for Taylor-Green vortex!\n";
                std::cout << "This indicates a numerical issue or bug in the solver.\n";
                break;
            }

            t_next_output += dt_output;
        }
    }

    // Final summary
    std::cout << std::string(80, '=') << "\n";
    std::cout << "\n";

    std::cout << "SIMULATION COMPLETE!\n";
    std::cout << "  Total steps: " << step << "\n";
    std::cout << "  Final time:  " << solver.getTime() << "\n";
    std::cout << "\n";

    // Save time series data
    std::string output_file = "taylor_green_timeseries.dat";
    solver.saveTimeSeries(output_file);
    std::cout << "Time series data saved to: " << output_file << "\n";
    std::cout << "\n";

    // Analysis
    double E_final = solver.computeEnergy();
    double decay_rate = -std::log(E_final / E0) / T_final;
    double theoretical_decay = 2.0 * nu;  // For k=1 mode

    std::cout << "ENERGY DECAY ANALYSIS:\n";
    std::cout << "  E₀ = " << E0 << "\n";
    std::cout << "  E(T) = " << E_final << "\n";
    std::cout << "  Decay ratio: " << E_final / E0 << "\n";
    std::cout << "  Measured decay rate: " << decay_rate << "\n";
    std::cout << "  Theoretical (2νk²):  " << theoretical_decay << "\n";
    std::cout << "  Relative error: "
              << std::abs(decay_rate - theoretical_decay) / theoretical_decay * 100
              << "%\n";
    std::cout << "\n";

    std::cout << "VISUALIZATION:\n";
    std::cout << "  Plot energy vs time using:\n";
    std::cout << "    gnuplot -e \"plot '" << output_file
              << "' using 1:2 with lines; pause -1\"\n";
    std::cout << "\n";

    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║              ✓ TAYLOR-GREEN SIMULATION COMPLETE          ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";

    return 0;
}
