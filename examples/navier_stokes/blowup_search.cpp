/**
 * @file blowup_search.cpp
 * @brief Systematic Search for Navier-Stokes Blow-up Solutions
 *
 * MILLENNIUM PRIZE PROBLEM:
 * Prove or disprove the global existence and smoothness of solutions to the
 * 3D incompressible Navier-Stokes equations.
 *
 * STRATEGY:
 * This program searches for finite-time singularities by:
 * 1. Testing multiple initial conditions (ABC flow, vortex rings, random)
 * 2. Varying Reynolds numbers (low to high)
 * 3. Monitoring BKM criterion, enstrophy, vorticity growth
 * 4. Detecting potential blow-up candidates
 *
 * INITIAL CONDITIONS TESTED:
 * - ABC (Arnold-Beltrami-Childress) flow: Chaotic, potential instabilities
 * - Vortex rings: Concentrated vorticity, vortex stretching
 * - Random perturbations: Generic search
 * - Kolmogorov flow: Shear instabilities
 *
 * USAGE:
 *   ./blowup_search [N] [Re_min] [Re_max] [num_Re]
 *
 *   N       : Grid resolution (default: 64)
 *   Re_min  : Minimum Reynolds number (default: 100)
 *   Re_max  : Maximum Reynolds number (default: 5000)
 *   num_Re  : Number of Re values to test (default: 5)
 *
 * OUTPUT:
 *   blowup_search_results.dat - Summary of all runs
 *   Console: Real-time monitoring
 *
 * NOTE: A detected blow-up could be worth $1,000,000!
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "../../new_theory/navier_stokes_solver.hpp"

using namespace new_theory::navier_stokes;

struct SimulationResult {
    std::string ic_name;
    double Re;
    double nu;
    double T_final;
    bool completed;
    bool blowup_detected;
    std::string criterion_violated;
    double max_vorticity;
    double final_enstrophy;
};

void print_header() {
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║     NAVIER-STOKES BLOW-UP SEARCH - MILLENNIUM PRIZE     ║\n";
    std::cout << "║        Systematic Search for Finite-Time Singularities  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    std::cout << "Prize: $1,000,000 for proof OR counterexample\n";
    std::cout << "Goal:  Find finite-time blow-up or verify regularity\n";
    std::cout << "\n";
}

void run_simulation(
    NavierStokesSolver& solver,
    const std::string& ic_name,
    double T_final,
    SimulationResult& result)
{
    std::cout << "  Running: " << std::setw(20) << std::left << ic_name
              << " Re = " << std::setw(8) << result.Re << "... " << std::flush;

    result.ic_name = ic_name;
    result.T_final = T_final;
    result.completed = false;
    result.blowup_detected = false;
    result.max_vorticity = 0.0;

    // Run simulation
    auto status = solver.runSimulation(
        T_final,
        0.05,      // dt_output
        0.01,      // dt_max
        true       // use_adaptive
    );

    // Record results
    result.completed = true;
    result.blowup_detected = !status.is_regular;
    result.criterion_violated = status.criterion_violated;
    result.max_vorticity = solver.computeVorticityMax();
    result.final_enstrophy = solver.computeEnstrophy();

    if (result.blowup_detected) {
        std::cout << "🔥 BLOW-UP DETECTED! 🔥\n";
        std::cout << "      Criterion: " << status.criterion_violated << "\n";
        std::cout << "      ||ω||_∞ = " << result.max_vorticity << "\n";
    } else {
        std::cout << "✓ Regular (||ω||_∞ = " << std::scientific
                  << std::setprecision(2) << result.max_vorticity << ")\n";
    }
}

int main(int argc, char** argv) {
    // Parse command line arguments
    int N = 64;
    double Re_min = 100.0;
    double Re_max = 5000.0;
    int num_Re = 5;

    if (argc > 1) N = std::atoi(argv[1]);
    if (argc > 2) Re_min = std::atof(argv[2]);
    if (argc > 3) Re_max = std::atof(argv[3]);
    if (argc > 4) num_Re = std::atoi(argv[4]);

    print_header();

    std::cout << "SEARCH PARAMETERS:\n";
    std::cout << "  Grid resolution:   " << N << "³\n";
    std::cout << "  Reynolds range:    " << Re_min << " - " << Re_max << "\n";
    std::cout << "  Number of Re:      " << num_Re << "\n";
    std::cout << "  Simulation time:   T = 5.0\n";
    std::cout << "\n";

    // Generate Reynolds numbers (logarithmic spacing)
    std::vector<double> Re_values;
    if (num_Re == 1) {
        Re_values.push_back(Re_min);
    } else {
        double log_min = std::log(Re_min);
        double log_max = std::log(Re_max);
        for (int i = 0; i < num_Re; ++i) {
            double log_Re = log_min + i * (log_max - log_min) / (num_Re - 1);
            Re_values.push_back(std::exp(log_Re));
        }
    }

    // Initial conditions to test
    struct ICConfig {
        std::string name;
        std::function<void(NavierStokesSolver&)> setup;
    };

    std::vector<ICConfig> initial_conditions = {
        {"Taylor-Green", [](NavierStokesSolver& s) { s.setInitialConditionTaylorGreen(); }},
        {"ABC (A=B=C=1)", [](NavierStokesSolver& s) { s.setInitialConditionABC(1.0, 1.0, 1.0); }},
        {"ABC (A=√3,B=√2,C=1)", [](NavierStokesSolver& s) {
            s.setInitialConditionABC(std::sqrt(3), std::sqrt(2), 1.0);
        }},
        {"Vortex Ring", [](NavierStokesSolver& s) {
            s.setInitialConditionVortexRing(1.0, 0.2, 1.0);
        }},
        {"Kolmogorov n=2", [](NavierStokesSolver& s) { s.setInitialConditionKolmogorov(2); }},
        {"Random seed=42", [](NavierStokesSolver& s) { s.setInitialConditionRandom(1.0, 42); }},
    };

    std::vector<SimulationResult> results;
    double T_final = 5.0;

    std::cout << "STARTING SYSTEMATIC SEARCH:\n";
    std::cout << std::string(80, '=') << "\n";

    int total_runs = Re_values.size() * initial_conditions.size();
    int current_run = 0;

    for (double Re : Re_values) {
        double nu = 1.0 / Re;

        for (const auto& ic : initial_conditions) {
            current_run++;
            std::cout << "[" << current_run << "/" << total_runs << "] ";

            // Create solver
            SpectralGrid3D grid(N, 2.0 * M_PI);
            NavierStokesSolver solver(grid, nu);

            // Set initial condition
            ic.setup(solver);

            // Run simulation
            SimulationResult result;
            result.Re = Re;
            result.nu = nu;

            run_simulation(solver, ic.name, T_final, result);
            results.push_back(result);

            // If blow-up detected, provide detailed info
            if (result.blowup_detected) {
                std::cout << "\n";
                std::cout << "╔══════════════════════════════════════════════════════════╗\n";
                std::cout << "║           🎆 POTENTIAL BLOW-UP CANDIDATE! 🎆             ║\n";
                std::cout << "╚══════════════════════════════════════════════════════════╝\n";
                std::cout << "  Initial condition: " << ic.name << "\n";
                std::cout << "  Reynolds number:   Re = " << Re << "\n";
                std::cout << "  Viscosity:         ν = " << nu << "\n";
                std::cout << "  Violated criterion: " << result.criterion_violated << "\n";
                std::cout << "  Max vorticity:     ||ω||_∞ = " << result.max_vorticity << "\n";
                std::cout << "\n";
                std::cout << "ACTION REQUIRED:\n";
                std::cout << "  1. Verify with higher resolution (N=" << 2*N << ")\n";
                std::cout << "  2. Check if blow-up persists under grid refinement\n";
                std::cout << "  3. Analyze vorticity structure\n";
                std::cout << "  4. If confirmed → PUBLISH AND CLAIM PRIZE!\n";
                std::cout << "\n";
            }
        }
    }

    std::cout << std::string(80, '=') << "\n";
    std::cout << "\n";

    // Summary
    std::cout << "SEARCH COMPLETE!\n";
    std::cout << "  Total simulations:  " << results.size() << "\n";

    int num_blowup = 0;
    for (const auto& r : results) {
        if (r.blowup_detected) num_blowup++;
    }

    std::cout << "  Blow-ups detected:  " << num_blowup << "\n";
    std::cout << "  Regular solutions:  " << (results.size() - num_blowup) << "\n";
    std::cout << "\n";

    // Save results
    std::string output_file = "blowup_search_results.dat";
    std::ofstream file(output_file);
    file << "# Navier-Stokes Blow-up Search Results\n";
    file << "# Grid: " << N << "³, T_final = " << T_final << "\n";
    file << "# Columns: IC_name Re nu Blowup? Criterion ||ω||_∞ Ω_final\n";
    file << "#\n";

    for (const auto& r : results) {
        file << std::setw(25) << r.ic_name << " "
             << std::setw(10) << r.Re << " "
             << std::scientific << std::setprecision(6)
             << std::setw(14) << r.nu << " "
             << std::setw(8) << (r.blowup_detected ? "BLOWUP" : "REGULAR") << " "
             << std::setw(20) << r.criterion_violated << " "
             << std::setw(14) << r.max_vorticity << " "
             << std::setw(14) << r.final_enstrophy << "\n";
    }
    file.close();

    std::cout << "Results saved to: " << output_file << "\n";
    std::cout << "\n";

    if (num_blowup > 0) {
        std::cout << "╔══════════════════════════════════════════════════════════╗\n";
        std::cout << "║  ⚠️  BLOW-UP CANDIDATES FOUND - FURTHER VERIFICATION  ⚠️  ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════╝\n";
        std::cout << "\n";
        std::cout << "Next steps:\n";
        std::cout << "1. Re-run suspect cases with higher resolution\n";
        std::cout << "2. Check grid convergence\n";
        std::cout << "3. Analyze physical vs numerical blow-up\n";
        std::cout << "4. If physical → mathematical proof required for prize\n";
    } else {
        std::cout << "╔══════════════════════════════════════════════════════════╗\n";
        std::cout << "║     ✓ NO BLOW-UP DETECTED IN TESTED CONFIGURATIONS      ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════╝\n";
        std::cout << "\n";
        std::cout << "This supports regularity, but more testing needed:\n";
        std::cout << "- Higher Reynolds numbers (Re > " << Re_max << ")\n";
        std::cout << "- Longer simulation times (T > " << T_final << ")\n";
        std::cout << "- More exotic initial conditions\n";
        std::cout << "- Higher resolution simulations\n";
    }
    std::cout << "\n";

    return 0;
}
