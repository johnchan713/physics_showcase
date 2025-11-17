#include <iostream>

// This comprehensive test verifies that all modules compile successfully
// Functional testing for specific modules is in test_new_modules.cpp

// Mathematics modules
#include "maths/matrices.hpp"
#include "maths/vectors.hpp"
#include "maths/complex_analysis.hpp"
#include "maths/fourier_analysis.hpp"
#include "maths/topology.hpp"
#include "maths/number_theory.hpp"
#include "maths/group_theory_lie_groups.hpp"
#include "maths/differential_algebra.hpp"
#include "maths/calculus_theorems.hpp"
#include "maths/trigonometry_identities.hpp"
#include "maths/polar_transforms.hpp"
#include "maths/distributions.hpp"
#include "maths/monte_carlo.hpp"
#include "maths/black_scholes.hpp"
#include "maths/variational_calculus.hpp"
#include "maths/partial_differential_equations.hpp"
#include "maths/ode_dynamical_systems.hpp"
#include "maths/stochastic_differential_equations.hpp"
#include "maths/pde_numerical_methods.hpp"
#include "maths/pde_solution_methods.hpp"
#include "maths/pde_transform_methods.hpp"
#include "maths/pde_variational_methods.hpp"
#include "maths/pde_classification_solutions.hpp"
#include "maths/econometrics_regression.hpp"
#include "maths/actuarial_life_tables.hpp"
#include "maths/advanced_subdifferentials.hpp"
#include "maths/nonsmooth_algorithms.hpp"

// Newly added math modules
#include "maths/measure_theory.hpp"
#include "maths/functional_analysis.hpp"
#include "maths/differential_geometry.hpp"
#include "maths/probability_theory.hpp"
#include "maths/real_analysis.hpp"

// Physics modules
#include "physics/units.hpp"
#include "physics/kinematics.hpp"
#include "physics/dynamics.hpp"
#include "physics/newton_laws.hpp"
#include "physics/energy_momentum.hpp"
#include "physics/circular_motion.hpp"
#include "physics/harmonic_motion.hpp"
#include "physics/projectile.hpp"
#include "physics/rotational_dynamics.hpp"
#include "physics/gravitation.hpp"
#include "physics/orbital.hpp"
#include "physics/elasticity.hpp"
#include "physics/fluid_mechanics.hpp"
#include "physics/inclined_plane.hpp"
#include "physics/surface_tension.hpp"
#include "physics/thermal_expansion.hpp"
#include "physics/calorimetry.hpp"
#include "physics/heat_transfer.hpp"
#include "physics/thermodynamics.hpp"
#include "physics/electrostatics.hpp"
#include "physics/magnetism.hpp"
#include "physics/electric_circuits.hpp"
#include "physics/electromagnetic_induction.hpp"
#include "physics/electromagnetic_waves.hpp"
#include "physics/maxwell_equations.hpp"
#include "physics/optics.hpp"
#include "physics/advanced_optics.hpp"
#include "physics/wave_mechanics.hpp"
#include "physics/oscillations.hpp"
#include "physics/special_relativity.hpp"
#include "physics/quantum_basics.hpp"
#include "physics/advanced_quantum_mechanics.hpp"
#include "physics/quantum_chemistry.hpp"
#include "physics/quantum_foundations.hpp"
#include "physics/relativistic_quantum_mechanics.hpp"
#include "physics/nuclear_physics.hpp"
#include "physics/advanced_mechanics.hpp"
#include "physics/statistical_models.hpp"

// Advanced physics modules
#include "physics/cosmology_friedmann_equations.hpp"
#include "physics/cosmology_expanding_universe.hpp"
#include "physics/cosmology_energy_density.hpp"
#include "physics/cosmology_early_universe.hpp"

#include "physics/gauge_theory_gauge_invariance.hpp"
#include "physics/gauge_theory_symmetries.hpp"
#include "physics/gauge_theory_higgs_mechanism.hpp"
#include "physics/gauge_theory_running_couplings.hpp"
#include "physics/gauge_theory_helicity.hpp"
#include "physics/gauge_theory_cp_violation_kaons.hpp"

#include "physics/qft_interactions.hpp"
#include "physics/qft_particle_physics.hpp"
#include "physics/qft_antiparticles.hpp"
#include "physics/qft_decays.hpp"
#include "physics/qft_cross_sections.hpp"
#include "physics/qft_spin_statistics.hpp"
#include "physics/qft_supersymmetry.hpp"
#include "physics/qft_quark_gluon_plasma.hpp"

#include "physics/loop_quantum_gravity.hpp"
#include "physics/operator_algebras.hpp"

// Newly added physics modules
#include "physics/general_relativity.hpp"
#include "physics/statistical_mechanics.hpp"
#include "physics/classical_field_theory.hpp"
#include "physics/condensed_matter.hpp"

// Fluid dynamics modules (that don't require Eigen)
#include "physics/fluid_dynamics_boundary_layer.hpp"
#include "physics/fluid_dynamics_compressible_flow.hpp"
#include "physics/fluid_dynamics_dimensionless_numbers.hpp"

// Classical mechanics modules (now using our own Vector/Matrix instead of Eigen!)
#include "physics/classical_hamiltonian.hpp"

// More fluid dynamics modules (Eigen replaced with our Vector class!)
#include "physics/fluid_dynamics_flow_types.hpp"
#include "physics/fluid_dynamics_turbulence.hpp"
#include "physics/fluid_dynamics_vorticity.hpp"

using namespace std;

int main() {
    cout << "========================================\n";
    cout << "  Compilation Test for ALL Modules\n";
    cout << "========================================\n";

    cout << "\nTesting module inclusions...\n\n";

    int math_count = 32;  // All 32 math modules compile!
    int physics_count = 69;  // 69 physics modules compile!
    int total_count = math_count + physics_count;

    cout << "Mathematics modules:      " << math_count << " ✓\n";
    cout << "Physics modules:          " << physics_count << " ✓\n";
    cout << "Total modules compiled:   " << total_count << "\n\n";

    cout << "Excluded modules (with bugs):\n";
    cout << "  None - ALL BUGS FIXED!\n\n";

    cout << "Excluded modules (require external Eigen library):\n";
    cout << "  Fluid dynamics (1 module - heavy Eigen::VectorXd/MatrixXd usage):\n";
    cout << "    - fluid_dynamics_governing_equations.hpp\n";
    cout << "  Classical mechanics (3 modules - complex matrix operations):\n";
    cout << "    - classical_liouville.hpp (uses advanced matrix features)\n";
    cout << "    - classical_phase_space.hpp (uses .block(), .norm())\n";
    cout << "    - physics_advanced.hpp (aggregator, includes above modules)\n\n";
    cout << "SUCCESS: Replaced Eigen in 4 modules with our own Vector/Matrix classes!\n";
    cout << "  - classical_hamiltonian.hpp (Hamiltonian mechanics)\n";
    cout << "  - fluid_dynamics_flow_types.hpp (Potential flow, Stokes flow, etc.)\n";
    cout << "  - fluid_dynamics_turbulence.hpp (RANS turbulence models)\n";
    cout << "  - fluid_dynamics_vorticity.hpp (Vorticity dynamics, circulation)\n\n";

    cout << "========================================\n";
    cout << "  ✓ ALL " << total_count << " MODULES COMPILED SUCCESSFULLY!\n";
    cout << "========================================\n";

    return 0;
}
