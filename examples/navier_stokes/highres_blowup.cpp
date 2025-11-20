/**
 * @file highres_blowup.cpp
 * @brief High-Resolution Blow-Up Search for Millennium Prize
 *
 * SERIOUS MILLENNIUM PRIZE ATTEMPT
 *
 * This program performs production-quality high-resolution simulations
 * targeting the discovery of finite-time singularities in 3D Navier-Stokes.
 *
 * FEATURES:
 * - High-resolution grids (128³, 256³, 512³)
 * - Checkpoint/restart for long runs
 * - VTK output for 3D visualization
 * - Automatic memory management
 * - Progress monitoring and blow-up detection
 *
 * RESOLUTION CAPABILITIES:
 * - 64³  = 262,144 points   (~60 MB RAM)
 * - 128³ = 2,097,152 points  (~480 MB RAM)
 * - 256³ = 16,777,216 points (~3.8 GB RAM)
 * - 512³ = 134,217,728 points (~30 GB RAM) - requires HPC
 *
 * USAGE:
 *   ./highres_blowup <N> <Re> <IC_type> [checkpoint_file]
 *
 *   N            : Grid resolution (64, 128, 256, 512)
 *   Re           : Reynolds number
 *   IC_type      : Initial condition (abc, vortex, random, kolmogorov)
 *   checkpoint   : (Optional) Resume from checkpoint
 *
 * EXAMPLES:
 *   ./highres_blowup 128 1000 abc
 *   ./highres_blowup 256 5000 vortex
 *   ./highres_blowup 128 1000 abc checkpoint.dat
 *
 * OUTPUT:
 *   - checkpoint_*.dat   : Restart files (every 1.0 time units)
 *   - vtk/field_*.vtk    : 3D visualization files
 *   - highres_results.dat: Time series data
 *
 * NOTE: For 256³+ simulations, compile with FFTW for performance!
 *       cmake .. -DUSE_FFTW=ON -DCMAKE_BUILD_TYPE=Release
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <sys/stat.h>
#include "../../new_theory/navier_stokes_solver.hpp"

using namespace new_theory::navier_stokes;

void print_header() {
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║   HIGH-RESOLUTION NAVIER-STOKES BLOW-UP SEARCH          ║\n";
    std::cout << "║        MILLENNIUM PRIZE PROBLEM ($1,000,000)            ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
}

void create_directory(const std::string& dir) {
    mkdir(dir.c_str(), 0755);  // Create if doesn't exist
}

void set_initial_condition(NavierStokesSolver& solver, const std::string& ic_type) {
    std::cout << "Setting initial condition: " << ic_type << "...\n";

    if (ic_type == "abc") {
        solver.setInitialConditionABC(std::sqrt(3), std::sqrt(2), 1.0);
    } else if (ic_type == "vortex") {
        solver.setInitialConditionVortexRing(1.0, 0.15, 2.0);
    } else if (ic_type == "random") {
        solver.setInitialConditionRandom(1.0, 42);
    } else if (ic_type == "kolmogorov") {
        solver.setInitialConditionKolmogorov(4);
    } else if (ic_type == "taylor-green") {
        solver.setInitialConditionTaylorGreen();
    } else {
        std::cerr << "Unknown IC type: " << ic_type << "\n";
        std::cerr << "Valid types: abc, vortex, random, kolmogorov, taylor-green\n";
        exit(1);
    }

    std::cout << "Initial condition set.\n";
}

void print_status(double t, double E, double Omega, double omega_max,
                 double dt, int step, bool regular) {
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << t
              << std::setw(15) << std::scientific << std::setprecision(4) << E
              << std::setw(15) << Omega
              << std::setw(15) << omega_max
              << std::setw(12) << std::fixed << std::setprecision(5) << dt
              << std::setw(10) << step
              << std::setw(12) << (regular ? "✓" : "BLOW-UP!") << "\n";
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <N> <Re> <IC_type> [checkpoint]\n";
        std::cout << "\n";
        std::cout << "Arguments:\n";
        std::cout << "  N         : Grid resolution (64, 128, 256, 512)\n";
        std::cout << "  Re        : Reynolds number\n";
        std::cout << "  IC_type   : Initial condition type\n";
        std::cout << "              (abc, vortex, random, kolmogorov, taylor-green)\n";
        std::cout << "  checkpoint: (Optional) Checkpoint file to resume from\n";
        std::cout << "\n";
        std::cout << "Examples:\n";
        std::cout << "  " << argv[0] << " 128 1000 abc\n";
        std::cout << "  " << argv[0] << " 256 5000 vortex\n";
        std::cout << "  " << argv[0] << " 128 1000 abc checkpoint_t5.dat\n";
        return 1;
    }

    // Parse arguments
    int N = std::atoi(argv[1]);
    double Re = std::atof(argv[2]);
    std::string ic_type = argv[3];
    std::string checkpoint_file = (argc > 4) ? argv[4] : "";
    bool resume = !checkpoint_file.empty();

    double nu = 1.0 / Re;

    print_header();

    // Print configuration
    std::cout << "SIMULATION CONFIGURATION:\n";
    std::cout << "  Grid resolution:    " << N << "³ = " << (N*N*N) << " points\n";
    std::cout << "  Reynolds number:    Re = " << Re << "\n";
    std::cout << "  Kinematic viscosity: ν = " << nu << "\n";
    std::cout << "  Initial condition:  " << ic_type << "\n";
    if (resume) {
        std::cout << "  Restart from:       " << checkpoint_file << "\n";
    }
    std::cout << "\n";

    // Memory estimate
    double mem_estimate = NavierStokesSolver::estimateMemoryUsage(N);
    std::cout << "MEMORY REQUIREMENTS:\n";
    std::cout << "  Estimated memory: " << std::fixed << std::setprecision(2)
              << mem_estimate << " MB";
    if (mem_estimate > 1024) {
        std::cout << " (" << mem_estimate/1024.0 << " GB)";
    }
    std::cout << "\n";

    if (mem_estimate > 4000 && mem_estimate <= 32000) {
        std::cout << "  ⚠️  WARNING: This simulation requires significant RAM!\n";
    } else if (mem_estimate > 32000) {
        std::cout << "  ⚠️  WARNING: This simulation requires HPC resources!\n";
        std::cout << "  Consider running on a cluster or reducing resolution.\n";
    }
    std::cout << "\n";

    // Create output directories
    create_directory("vtk");
    create_directory("checkpoints");

    // Create solver
    std::cout << "Initializing solver...\n";
    SpectralGrid3D grid(N, 2.0 * M_PI);
    NavierStokesSolver solver(grid, nu);
    solver.printMemoryUsage();
    std::cout << "\n";

    // Setup simulation
    double start_time = 0.0;

    if (resume) {
        // Load checkpoint
        std::cout << "Loading checkpoint...\n";
        if (!solver.loadCheckpoint(checkpoint_file)) {
            std::cerr << "Failed to load checkpoint!\n";
            return 1;
        }
        start_time = solver.getTime();
    } else {
        // Set initial condition
        set_initial_condition(solver, ic_type);
    }

    // Initial diagnostics
    double E0 = solver.computeEnergy();
    double Omega0 = solver.computeEnstrophy();
    std::cout << "Initial state:\n";
    std::cout << "  E₀ = " << std::scientific << E0 << "\n";
    std::cout << "  Ω₀ = " << Omega0 << "\n";
    std::cout << "\n";

    // Simulation parameters
    double T_final = 50.0;           // Long-time integration
    double dt_output = 0.1;          // Monitor every 0.1 time units
    double dt_checkpoint = 1.0;      // Checkpoint every 1.0 time units
    double dt_vtk = 2.0;             // VTK output every 2.0 time units
    double dt_max = 0.005;           // Maximum time step (smaller for stability)
    double CFL = 0.3;                // CFL number (conservative for high-res)

    // Progress tracking
    std::cout << "SIMULATION PROGRESS:\n";
    std::cout << std::string(95, '=') << "\n";
    std::cout << std::setw(10) << "Time"
              << std::setw(15) << "Energy"
              << std::setw(15) << "Enstrophy"
              << std::setw(15) << "||ω||_∞"
              << std::setw(12) << "dt"
              << std::setw(10) << "Steps"
              << std::setw(12) << "Status" << "\n";
    std::cout << std::string(95, '-') << "\n";

    double t_next_output = start_time + dt_output;
    double t_next_checkpoint = start_time + dt_checkpoint;
    double t_next_vtk = start_time + dt_vtk;
    int step = 0;
    int checkpoint_count = 0;
    int vtk_count = 0;

    // Main simulation loop
    while (solver.getTime() < T_final) {
        // Adaptive time step
        double dt = solver.computeCFLTimeStep(CFL);
        dt = std::min(dt, dt_max);
        dt = std::min(dt, T_final - solver.getTime());

        // Advance
        solver.stepRK4(dt);
        solver.updateRegularityMonitoring();
        step++;

        // Output monitoring
        if (solver.getTime() >= t_next_output || solver.getTime() >= T_final - 1e-10) {
            double E = solver.computeEnergy();
            double Omega = solver.computeEnstrophy();
            double omega_max = solver.computeVorticityMax();

            auto status = solver.checkRegularity();

            print_status(solver.getTime(), E, Omega, omega_max, dt, step, status.is_regular);

            // Check for blow-up
            if (!status.is_regular) {
                std::cout << "\n";
                std::cout << "╔══════════════════════════════════════════════════════════╗\n";
                std::cout << "║           🎆🎆🎆 BLOW-UP DETECTED! 🎆🎆🎆                 ║\n";
                std::cout << "╚══════════════════════════════════════════════════════════╝\n";
                std::cout << "Time:               t = " << solver.getTime() << "\n";
                std::cout << "Criterion violated: " << status.criterion_violated << "\n";
                std::cout << "Severity score:     " << status.severity_score << "\n";
                std::cout << "Max vorticity:      ||ω||_∞ = " << omega_max << "\n";
                std::cout << "\n";
                std::cout << "IMMEDIATE ACTIONS:\n";
                std::cout << "1. Save emergency checkpoint\n";
                std::cout << "2. Export VTK for visualization\n";
                std::cout << "3. Re-run with 2× higher resolution to confirm\n";
                std::cout << "4. If confirmed → PUBLISH AND CLAIM $1,000,000!\n";
                std::cout << "\n";

                // Emergency save
                std::string emergency_chk = "checkpoints/BLOWUP_emergency.dat";
                solver.saveCheckpoint(emergency_chk);

                std::string emergency_vtk = "vtk/BLOWUP_field";
                solver.saveVTK(emergency_vtk, true);

                break;
            }

            t_next_output += dt_output;
        }

        // Checkpoint saving
        if (solver.getTime() >= t_next_checkpoint) {
            std::string chk_file = "checkpoints/checkpoint_t" +
                                  std::to_string((int)solver.getTime()) + ".dat";
            solver.saveCheckpoint(chk_file);
            checkpoint_count++;
            t_next_checkpoint += dt_checkpoint;
        }

        // VTK output
        if (solver.getTime() >= t_next_vtk) {
            std::string vtk_file = "vtk/field_t" +
                                  std::to_string((int)solver.getTime());
            solver.saveVTK(vtk_file, true);
            vtk_count++;
            t_next_vtk += dt_vtk;
        }
    }

    std::cout << std::string(95, '=') << "\n";
    std::cout << "\n";

    // Final summary
    std::cout << "SIMULATION COMPLETE!\n";
    std::cout << "  Total steps:        " << step << "\n";
    std::cout << "  Final time:         " << solver.getTime() << "\n";
    std::cout << "  Checkpoints saved:  " << checkpoint_count << "\n";
    std::cout << "  VTK files created:  " << vtk_count << "\n";
    std::cout << "\n";

    // Save final data
    solver.saveTimeSeries("highres_results.dat");
    solver.saveCheckpoint("checkpoints/final_state.dat");
    solver.saveVTK("vtk/final_field", true);

    std::cout << "Output files:\n";
    std::cout << "  Time series:    highres_results.dat\n";
    std::cout << "  Final state:    checkpoints/final_state.dat\n";
    std::cout << "  Final VTK:      vtk/final_field.vtk\n";
    std::cout << "\n";

    // Analysis
    double E_final = solver.computeEnergy();
    std::cout << "ENERGY ANALYSIS:\n";
    std::cout << "  E₀ = " << std::scientific << E0 << "\n";
    std::cout << "  E(T) = " << E_final << "\n";
    std::cout << "  Decay ratio: " << E_final / E0 << "\n";
    std::cout << "\n";

    auto final_status = solver.checkRegularity();
    if (final_status.is_regular) {
        std::cout << "╔══════════════════════════════════════════════════════════╗\n";
        std::cout << "║      ✓ SOLUTION REMAINED REGULAR UP TO t=" << std::setw(6)
                  << std::fixed << std::setprecision(1) << solver.getTime() << "       ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════╝\n";
        std::cout << "\n";
        std::cout << "This supports regularity conjecture for:\n";
        std::cout << "  Resolution: " << N << "³\n";
        std::cout << "  Re = " << Re << "\n";
        std::cout << "  IC: " << ic_type << "\n";
        std::cout << "\n";
        std::cout << "Consider:\n";
        std::cout << "- Higher Re for more challenging regime\n";
        std::cout << "- Longer integration time\n";
        std::cout << "- Different initial conditions\n";
    }

    std::cout << "\n";
    std::cout << "VISUALIZATION:\n";
    std::cout << "  Open VTK files in ParaView:\n";
    std::cout << "    paraview vtk/final_field.vtk\n";
    std::cout << "\n";
    std::cout << "  Create animation of time series:\n";
    std::cout << "    paraview vtk/field_t*.vtk\n";
    std::cout << "\n";

    return 0;
}
