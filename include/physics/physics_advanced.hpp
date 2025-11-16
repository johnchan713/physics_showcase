#ifndef PHYSICS_ADVANCED_HPP
#define PHYSICS_ADVANCED_HPP

/**
 * @file physics_advanced.hpp
 * @brief Central header for all advanced physics modules
 *
 * This is the main entry point for the advanced physics library.
 * Include this file to access all advanced features.
 *
 * Categories:
 * 1. Classical Mechanics (Hamiltonian, Phase Space, Liouville)
 * 2. Electromagnetism (Tensors, Gauge Theory)
 * 3. Relativity (General Relativity, Geodesics)
 * 4. Wave Theory (Green Functions, Waveguides, Diffraction)
 * 5. Thermodynamics (Statistical, Coefficients)
 * 6. Quantum Mechanics (Schrödinger, Dirac, Perturbation)
 * 7. Condensed Matter (Superconductivity, Magnetism)
 * 8. Plasma Physics (Transport, Collision-Radiative)
 * 9. Advanced Quantum Theory (Group Theory, Tensor Operators)
 *
 * Requirements:
 * - C++17 or later
 * - Eigen3 library
 * - Boost (optional, for special functions)
 */

// Category 1: Classical Mechanics
#include "classical/hamiltonian.hpp"
#include "classical/phase_space.hpp"
#include "classical/liouville.hpp"

// Category 2: Electromagnetism
#include "electromagnetism/field_tensor.hpp"
#include "electromagnetism/stress_energy_tensor.hpp"
#include "electromagnetism/gauge_theory.hpp"

// Category 3: Relativity
#include "relativity/metric_tensor.hpp"
#include "relativity/christoffel.hpp"
#include "relativity/riemann_tensor.hpp"
#include "relativity/einstein_tensor.hpp"
#include "relativity/geodesics.hpp"

// Category 4: Wave Theory
#include "waves/green_functions.hpp"
#include "waves/waveguides.hpp"
#include "waves/spherical_waves.hpp"
#include "waves/diffraction.hpp"

// Category 5: Thermodynamics
#include "thermodynamics/coefficients.hpp"
#include "thermodynamics/equation_of_state.hpp"

// Category 6: Quantum Mechanics
#include "quantum/schrodinger_solver.hpp"
#include "quantum/dirac_equation.hpp"
#include "quantum/perturbation_theory.hpp"
#include "quantum/spin_matrices.hpp"

// Category 7: Condensed Matter
#include "condensed_matter/superconductivity.hpp"
#include "condensed_matter/magnetism.hpp"

// Category 8: Plasma Physics
#include "plasma/plasma_state.hpp"
#include "plasma/transport.hpp"

// Category 9: Advanced Quantum
#include "advanced_quantum/clebsch_gordan.hpp"
#include "advanced_quantum/wigner_eckart.hpp"
#include "advanced_quantum/group_theory.hpp"

namespace physics::advanced {

/**
 * @brief Version information
 */
struct Version {
    static constexpr int MAJOR = 2;
    static constexpr int MINOR = 0;
    static constexpr int PATCH = 0;

    static std::string string() {
        return std::to_string(MAJOR) + "." +
               std::to_string(MINOR) + "." +
               std::to_string(PATCH);
    }
};

/**
 * @brief Initialize advanced physics library
 *
 * Call this before using advanced features.
 * Sets up numerical precision, threading, etc.
 */
inline void initialize() {
    // Set Eigen to use all available threads
    Eigen::initParallel();

    // Set numerical precision defaults
    // (Could customize based on application needs)
}

} // namespace physics::advanced

#endif // PHYSICS_ADVANCED_HPP
