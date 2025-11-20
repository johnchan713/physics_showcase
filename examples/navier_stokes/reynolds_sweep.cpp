/**
 * @file reynolds_sweep.cpp
 * @brief Reynolds Number Parameter Sweep for ABC Flow
 *
 * DESCRIPTION:
 * Systematically varies the Reynolds number to study the transition from
 * laminar to turbulent flow. Uses ABC (Arnold-Beltrami-Childress) flow
 * as a test case known for chaotic behavior at high Re.
 *
 * PHYSICS:
 * - Low Re (< 100):  Viscosity dominates, flow remains smooth
 * - Mid Re (100-1000): Transition regime, increasingly complex
 * - High Re (> 1000): Turbulent, potential for numerical challenges
 *
 * MEASURED QUANTITIES:
 * - Energy dissipation rate ε
 * - Enstrophy evolution
 * - Max vorticity ||ω||_∞
 * - Effective Reynolds number Re_eff based on flow statistics
 *
 * USAGE:
 *   ./reynolds_sweep [N] [Re_min] [Re_max] [num_points]
 *
 *   N          : Grid resolution (default: 64)
 *   Re_min     : Minimum Reynolds number (default: 10)
 *   Re_max     : Maximum Reynolds number (default: 1000)
 *   num_points : Number of Re values (default: 20)
 *
 * OUTPUT:
 *   reynolds_sweep.dat - Re vs E, Ω, ε, ||ω||_∞
 *   Console: Progress monitoring
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include "../../new_theory/navier_stokes_solver.hpp"

using namespace new_theory::navier_stokes;

void print_header() {
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║        REYNOLDS NUMBER PARAMETER SWEEP - ABC FLOW        ║\n";
    std::cout << "║           Transition to Turbulence Study                ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
}

struct SweepResult {
    double Re;
    double nu;
    double E_initial;
    double E_final;
    double Omega_initial;
    double Omega_final;
    double omega_Linfty_max;
    double dissipation_rate_avg;
    bool regular;
};

int main(int argc, char** argv) {
    // Parse arguments
    int N = 64;
    double Re_min = 10.0;
    double Re_max = 1000.0;
    int num_points = 20;

    if (argc > 1) N = std::atoi(argv[1]);
    if (argc > 2) Re_min = std::atof(argv[2]);
    if (argc > 3) Re_max = std::atof(argv[3]);
    if (argc > 4) num_points = std::atoi(argv[4]);

    print_header();

    std::cout << "SWEEP PARAMETERS:\n";
    std::cout << "  Grid resolution: " << N << "³ = " << (N*N*N) << " points\n";
    std::cout << "  Reynolds range:  " << Re_min << " ≤ Re ≤ " << Re_max << "\n";
    std::cout << "  Number of points: " << num_points << "\n";
    std::cout << "  Initial condition: ABC flow (A=B=C=1)\n";
    std::cout << "  Simulation time:   T = 5.0\n";
    std::cout << "\n";

    // Generate Reynolds numbers (logarithmic spacing for wide range)
    std::vector<double> Re_values;
    if (num_points == 1) {
        Re_values.push_back(Re_min);
    } else {
        double log_min = std::log(Re_min);
        double log_max = std::log(Re_max);
        for (int i = 0; i < num_points; ++i) {
            double log_Re = log_min + i * (log_max - log_min) / (num_points - 1);
            Re_values.push_back(std::exp(log_Re));
        }
    }

    std::vector<SweepResult> results;
    double T_final = 5.0;

    // Print progress header
    std::cout << "RUNNING SWEEP:\n";
    std::cout << std::string(90, '=') << "\n";
    std::cout << std::setw(5) << "#"
              << std::setw(12) << "Re"
              << std::setw(12) << "ν"
              << std::setw(12) << "E(T)"
              << std::setw(12) << "Ω(T)"
              << std::setw(15) << "||ω||_∞"
              << std::setw(12) << "ε_avg"
              << std::setw(10) << "Status" << "\n";
    std::cout << std::string(90, '-') << "\n";

    // Run simulations
    for (size_t i = 0; i < Re_values.size(); ++i) {
        double Re = Re_values[i];
        double nu = 1.0 / Re;

        // Create solver and set IC
        SpectralGrid3D grid(N, 2.0 * M_PI);
        NavierStokesSolver solver(grid, nu);
        solver.setInitialConditionABC(1.0, 1.0, 1.0);

        // Record initial state
        SweepResult result;
        result.Re = Re;
        result.nu = nu;
        result.E_initial = solver.computeEnergy();
        result.Omega_initial = solver.computeEnstrophy();

        // Run simulation
        auto status = solver.runSimulation(T_final, 0.1, 0.01, true);

        // Record final state
        result.E_final = solver.computeEnergy();
        result.Omega_final = solver.computeEnstrophy();
        result.omega_Linfty_max = solver.computeVorticityMax();
        result.dissipation_rate_avg = solver.computeDissipationRate();
        result.regular = status.is_regular;

        results.push_back(result);

        // Print progress
        std::cout << std::setw(5) << (i + 1)
                  << std::setw(12) << std::fixed << std::setprecision(2) << Re
                  << std::setw(12) << std::scientific << std::setprecision(3) << nu
                  << std::setw(12) << result.E_final
                  << std::setw(12) << result.Omega_final
                  << std::setw(15) << result.omega_Linfty_max
                  << std::setw(12) << result.dissipation_rate_avg
                  << std::setw(10) << (result.regular ? "✓" : "BLOW-UP") << "\n";
    }

    std::cout << std::string(90, '=') << "\n";
    std::cout << "\n";

    // Save results
    std::string output_file = "reynolds_sweep.dat";
    std::ofstream file(output_file);
    file << "# Reynolds Number Parameter Sweep - ABC Flow\n";
    file << "# Grid: " << N << "³, T_final = " << T_final << "\n";
    file << "# Columns: Re  ν  E₀  E(T)  Ω₀  Ω(T)  ||ω||_∞  ε_avg  Regular?\n";
    file << "#\n";

    for (const auto& r : results) {
        file << std::scientific << std::setprecision(8)
             << std::setw(16) << r.Re
             << std::setw(16) << r.nu
             << std::setw(16) << r.E_initial
             << std::setw(16) << r.E_final
             << std::setw(16) << r.Omega_initial
             << std::setw(16) << r.Omega_final
             << std::setw(16) << r.omega_Linfty_max
             << std::setw(16) << r.dissipation_rate_avg
             << std::setw(10) << (r.regular ? "1" : "0") << "\n";
    }
    file.close();

    std::cout << "SWEEP COMPLETE!\n";
    std::cout << "  Data saved to: " << output_file << "\n";
    std::cout << "\n";

    // Analysis
    std::cout << "ANALYSIS:\n";

    // Find transition region (where vorticity growth accelerates)
    double omega_threshold = 0.0;
    for (const auto& r : results) {
        omega_threshold = std::max(omega_threshold, r.omega_Linfty_max);
    }
    omega_threshold *= 0.5;  // 50% of max

    int transition_idx = -1;
    for (size_t i = 0; i < results.size() - 1; ++i) {
        if (results[i].omega_Linfty_max < omega_threshold &&
            results[i+1].omega_Linfty_max >= omega_threshold) {
            transition_idx = i;
            break;
        }
    }

    if (transition_idx >= 0) {
        std::cout << "  Transition to complex dynamics:\n";
        std::cout << "    Re ≈ " << results[transition_idx].Re << " - "
                  << results[transition_idx + 1].Re << "\n";
    }

    // Check for any blow-ups
    int num_blowup = 0;
    for (const auto& r : results) {
        if (!r.regular) num_blowup++;
    }

    if (num_blowup > 0) {
        std::cout << "  ⚠️  Blow-ups detected: " << num_blowup << " cases\n";
        std::cout << "  → Requires verification with higher resolution\n";
    } else {
        std::cout << "  ✓ All cases remained regular\n";
    }

    // Scaling analysis
    std::cout << "\n";
    std::cout << "  Dissipation rate scaling:\n";
    if (results.size() >= 2) {
        double eps_low = results[0].dissipation_rate_avg;
        double eps_high = results.back().dissipation_rate_avg;
        double Re_low = results[0].Re;
        double Re_high = results.back().Re;

        if (eps_low > 1e-12 && eps_high > 1e-12) {
            double scaling = std::log(eps_high / eps_low) / std::log(Re_high / Re_low);
            std::cout << "    ε(Re) ∝ Re^" << std::fixed << std::setprecision(2) << scaling << "\n";
            std::cout << "    (Turbulent theory predicts ε ≈ U³/L independent of Re)\n";
        }
    }

    std::cout << "\n";
    std::cout << "VISUALIZATION:\n";
    std::cout << "  Plot Reynolds scaling using gnuplot:\n";
    std::cout << "    gnuplot> set logscale xy\n";
    std::cout << "    gnuplot> plot '" << output_file << "' u 1:7 w lp title '||ω||_∞'\n";
    std::cout << "    gnuplot> replot '" << output_file << "' u 1:8 w lp title 'ε'\n";
    std::cout << "\n";

    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║           ✓ REYNOLDS SWEEP COMPLETE                     ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";

    return 0;
}
