# Physics Showcase Architecture - Advanced Topics Implementation

## Overview

This document describes how to extend the physics showcase from a header-only library to a full-featured C++ physics simulation framework supporting advanced topics.

## Current Architecture

**Header-Only Library** (31 modules)
- ✅ Pure functions with no state
- ✅ Compile-time evaluation possible
- ✅ Easy integration (just include headers)
- ❌ Limited to algebraic formulas
- ❌ Cannot handle matrix operations efficiently
- ❌ No numerical solvers
- ❌ No complex state management

## Proposed Extended Architecture

### Directory Structure

```
physics_showcase/
├── include/
│   └── physics/
│       ├── [existing 31 header files]
│       ├── tensors/              # NEW: Tensor mathematics
│       │   ├── tensor.hpp
│       │   ├── stress_energy_tensor.hpp
│       │   └── field_tensor.hpp
│       ├── quantum/               # NEW: Advanced quantum mechanics
│       │   ├── dirac_operator.hpp
│       │   ├── spin_matrices.hpp
│       │   ├── perturbation.hpp
│       │   └── quantum_state.hpp
│       ├── relativity/            # NEW: General relativity
│       │   ├── metric_tensor.hpp
│       │   ├── christoffel.hpp
│       │   ├── riemann_tensor.hpp
│       │   └── einstein_equations.hpp
│       ├── waves/                 # NEW: Advanced wave theory
│       │   ├── bessel_functions.hpp
│       │   ├── green_functions.hpp
│       │   └── waveguide.hpp
│       ├── statistical/           # NEW: Statistical mechanics
│       │   ├── phase_space.hpp
│       │   ├── liouville.hpp
│       │   └── boltzmann.hpp
│       └── plasma/                # NEW: Plasma physics
│           ├── plasma_state.hpp
│           └── transport.hpp
├── src/                           # NEW: Implementation files
│   ├── tensors/
│   ├── quantum/
│   ├── relativity/
│   ├── waves/
│   ├── statistical/
│   └── plasma/
├── external/                      # NEW: Third-party libraries
│   ├── eigen/                     # Linear algebra
│   ├── boost/                     # Special functions
│   └── gsl/                       # GNU Scientific Library
├── tests/                         # NEW: Unit tests
└── benchmarks/                    # NEW: Performance tests
```

## Implementation Approaches for Advanced Topics

### 1. Tensor Mathematics (Stress-Energy, Field Tensors)

**Approach:** Use Eigen library + OOP

```cpp
// include/physics/tensors/tensor.hpp
#pragma once
#include <Eigen/Dense>
#include <array>

namespace physics::tensors {

template<int Rank, int Dim>
class Tensor {
    // Store as flat array, provide index mapping
    std::array<int, Rank> dimensions_;
    std::vector<double> data_;
public:
    double& operator()(const std::array<int, Rank>& indices);
    Tensor<Rank, Dim> contract(int index1, int index2) const;
};

// Rank-2 tensor (matrix-like)
using Tensor2D = Eigen::Matrix<double, 4, 4>;

class StressEnergyTensor {
    Tensor2D T_;  // T^μν
public:
    StressEnergyTensor(double energy_density,
                       const Eigen::Vector3d& momentum_density,
                       const Eigen::Matrix3d& stress);

    double energyDensity() const { return T_(0, 0); }
    Eigen::Vector3d momentumDensity() const;
    double trace() const { return T_.trace(); }

    // Contract with metric
    Tensor2D lower(const Tensor2D& metric) const;
};

class ElectromagneticFieldTensor {
    Tensor2D F_;  // F^μν (antisymmetric)
public:
    ElectromagneticFieldTensor(const Eigen::Vector3d& E,
                               const Eigen::Vector3d& B);

    Eigen::Vector3d electricField() const;
    Eigen::Vector3d magneticField() const;

    // Lorentz transformation
    ElectromagneticFieldTensor boost(const Eigen::Vector3d& velocity) const;

    // Invariants
    double firstInvariant() const;  // E²-B²
    double secondInvariant() const; // E·B
};

} // namespace physics::tensors
```

```cpp
// src/tensors/field_tensor.cpp
#include "physics/tensors/field_tensor.hpp"

namespace physics::tensors {

ElectromagneticFieldTensor::ElectromagneticFieldTensor(
    const Eigen::Vector3d& E, const Eigen::Vector3d& B) {

    F_.setZero();

    // F^0i = E_i/c
    F_(0, 1) = E(0);
    F_(0, 2) = E(1);
    F_(0, 3) = E(2);

    // F^ij = -ε_ijk B_k
    F_(1, 2) = -B(2);
    F_(1, 3) =  B(1);
    F_(2, 3) = -B(0);

    // Antisymmetric
    F_(1, 0) = -F_(0, 1);
    F_(2, 0) = -F_(0, 2);
    F_(3, 0) = -F_(0, 3);
    F_(2, 1) = -F_(1, 2);
    F_(3, 1) = -F_(1, 3);
    F_(3, 2) = -F_(2, 3);
}

double ElectromagneticFieldTensor::firstInvariant() const {
    // I₁ = (1/2)F_μν F^μν = B² - E²/c²
    Eigen::Vector3d E = electricField();
    Eigen::Vector3d B = magneticField();
    return B.squaredNorm() - E.squaredNorm();
}

} // namespace physics::tensors
```

### 2. Dirac Formalism & Spin Matrices

**Approach:** Matrix classes with operator overloading

```cpp
// include/physics/quantum/dirac_operator.hpp
#pragma once
#include <Eigen/Dense>
#include <complex>

namespace physics::quantum {

using Complex = std::complex<double>;
using Matrix2cd = Eigen::Matrix<Complex, 2, 2>;
using Matrix4cd = Eigen::Matrix<Complex, 4, 4>;
using Vector4cd = Eigen::Vector<Complex, 4>;

// Pauli spin matrices
class PauliMatrices {
public:
    static Matrix2cd sigma_x();
    static Matrix2cd sigma_y();
    static Matrix2cd sigma_z();
    static std::array<Matrix2cd, 3> all();
};

// Dirac gamma matrices
class DiracMatrices {
    std::array<Matrix4cd, 4> gamma_;  // γ^μ
public:
    DiracMatrices();  // Standard representation

    const Matrix4cd& gamma(int mu) const { return gamma_[mu]; }
    Matrix4cd gamma5() const;  // γ⁵ = iγ⁰γ¹γ²γ³

    // Anticommutation: {γ^μ, γ^ν} = 2g^μν
    bool verifyAnticommutation() const;
};

// Dirac spinor (4-component)
class DiracSpinor {
    Vector4cd components_;
public:
    DiracSpinor(const Vector4cd& psi) : components_(psi) {}

    // Adjoint: ψ̄ = ψ†γ⁰
    Eigen::RowVector<Complex, 4> adjoint(const DiracMatrices& gamma) const;

    // Normalization
    Complex norm() const { return components_.squaredNorm(); }
    void normalize();

    // Projection operators
    DiracSpinor projectPositiveEnergy(const DiracMatrices& gamma) const;
    DiracSpinor projectNegativeEnergy(const DiracMatrices& gamma) const;
};

// Dirac equation solver
class DiracEquation {
    DiracMatrices gamma_;
    double mass_;
    double hbar_;
    double c_;
public:
    DiracEquation(double mass);

    // Free particle solutions
    DiracSpinor plane_wave(const Eigen::Vector3d& momentum,
                           double energy, int spin) const;

    // Hamiltonian
    Matrix4cd hamiltonian(const Eigen::Vector3d& momentum) const;

    // Spin-orbit coupling
    Matrix4cd spinOrbitCoupling(const Eigen::Vector3d& position,
                                const std::function<double(Eigen::Vector3d)>& potential) const;
};

} // namespace physics::quantum
```

### 3. Schrödinger Equation (Numerical Solver)

**Approach:** Finite difference method with OOP

```cpp
// include/physics/quantum/schrodinger_solver.hpp
#pragma once
#include <Eigen/Dense>
#include <functional>

namespace physics::quantum {

class SchrodingerSolver1D {
public:
    struct Parameters {
        double x_min, x_max;
        int num_points;
        double mass;
        double hbar;
    };

private:
    Parameters params_;
    Eigen::VectorXd x_grid_;
    Eigen::VectorXd potential_;

    // Hamiltonian matrix (sparse)
    Eigen::MatrixXd H_;

    // Eigenstates
    Eigen::VectorXd energies_;
    Eigen::MatrixXd wavefunctions_;

public:
    SchrodingerSolver1D(const Parameters& params,
                        const std::function<double(double)>& V);

    // Solve time-independent Schrödinger equation
    void solveEigenvalueProblem(int num_states);

    // Get results
    double getEnergy(int n) const { return energies_(n); }
    Eigen::VectorXd getWavefunction(int n) const { return wavefunctions_.col(n); }

    // Time evolution
    Eigen::VectorXcd evolve(const Eigen::VectorXcd& psi0, double dt, int steps) const;

    // Expectation values
    double expectation_position(const Eigen::VectorXcd& psi) const;
    double expectation_momentum(const Eigen::VectorXcd& psi) const;
    double expectation_energy(const Eigen::VectorXcd& psi) const;
};

// 2D/3D versions
class SchrodingerSolver2D { /* Similar structure */ };
class SchrodingerSolver3D { /* Similar structure */ };

} // namespace physics::quantum
```

```cpp
// src/quantum/schrodinger_solver.cpp
#include "physics/quantum/schrodinger_solver.hpp"
#include <Eigen/Eigenvalues>

namespace physics::quantum {

SchrodingerSolver1D::SchrodingerSolver1D(
    const Parameters& params,
    const std::function<double(double)>& V)
    : params_(params) {

    // Create spatial grid
    int N = params_.num_points;
    x_grid_ = Eigen::VectorXd::LinSpaced(N, params_.x_min, params_.x_max);
    double dx = x_grid_(1) - x_grid_(0);

    // Evaluate potential on grid
    potential_.resize(N);
    for (int i = 0; i < N; ++i) {
        potential_(i) = V(x_grid_(i));
    }

    // Build Hamiltonian using finite differences
    // H = T + V where T is kinetic energy operator
    H_.resize(N, N);
    H_.setZero();

    double coeff = -params_.hbar * params_.hbar / (2.0 * params_.mass * dx * dx);

    for (int i = 0; i < N; ++i) {
        H_(i, i) = -2.0 * coeff + potential_(i);  // Diagonal
        if (i > 0) H_(i, i-1) = coeff;            // Off-diagonal
        if (i < N-1) H_(i, i+1) = coeff;
    }
}

void SchrodingerSolver1D::solveEigenvalueProblem(int num_states) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_);

    energies_ = solver.eigenvalues().head(num_states);
    wavefunctions_ = solver.eigenvectors().leftCols(num_states);

    // Normalize wavefunctions
    double dx = x_grid_(1) - x_grid_(0);
    for (int n = 0; n < num_states; ++n) {
        double norm = std::sqrt((wavefunctions_.col(n).array().square() * dx).sum());
        wavefunctions_.col(n) /= norm;
    }
}

} // namespace physics::quantum
```

### 4. Green Functions & Wave Equations

**Approach:** Numerical methods + special functions

```cpp
// include/physics/waves/green_functions.hpp
#pragma once
#include <functional>
#include <complex>

namespace physics::waves {

using Complex = std::complex<double>;

// 1D Green's function for wave equation
class GreenFunction1D {
    double wave_speed_;
public:
    GreenFunction1D(double c) : wave_speed_(c) {}

    // Retarded Green's function: G_ret(x,t;x',t')
    double retarded(double x, double t, double x_prime, double t_prime) const;

    // Advanced Green's function
    double advanced(double x, double t, double x_prime, double t_prime) const;

    // Solve inhomogeneous wave equation
    std::function<double(double,double)> solve(
        const std::function<double(double,double)>& source) const;
};

// 3D Green's function
class GreenFunction3D {
    double wave_speed_;
public:
    // Retarded: G_ret(r,t;r',t') = δ(t-t'-|r-r'|/c)/(4π|r-r'|)
    double retarded(const Eigen::Vector3d& r, double t,
                   const Eigen::Vector3d& r_prime, double t_prime) const;
};

// Helmholtz equation Green's function
class HelmholtzGreen {
    double k_;  // Wave number
public:
    HelmholtzGreen(double k) : k_(k) {}

    // Free space: G(r) = -e^(ik|r-r'|)/(4π|r-r'|)
    Complex freeSpace(const Eigen::Vector3d& r,
                     const Eigen::Vector3d& r_prime) const;
};

} // namespace physics::waves
```

### 5. Waveguides & Resonating Cavities

**Approach:** Boundary value problem solver

```cpp
// include/physics/waves/waveguide.hpp
#pragma once
#include <vector>
#include <Eigen/Dense>

namespace physics::waves {

// Rectangular waveguide
class RectangularWaveguide {
    double a_, b_;  // Dimensions
    double epsilon_, mu_;  // Material properties

public:
    RectangularWaveguide(double a, double b, double epsilon, double mu);

    // TE modes: E_z = 0
    struct TEMode {
        int m, n;  // Mode indices
        double cutoff_frequency;
        double propagation_constant(double frequency) const;
        double phase_velocity(double frequency) const;
        double group_velocity(double frequency) const;
    };

    std::vector<TEMode> getTEModes(int max_m, int max_n) const;

    // TM modes: H_z = 0
    struct TMMode {
        int m, n;
        double cutoff_frequency;
        double propagation_constant(double frequency) const;
    };

    std::vector<TMMode> getTMModes(int max_m, int max_n) const;

    // Field patterns
    Eigen::MatrixXd electricFieldPattern(const TEMode& mode, double z, double t) const;
    Eigen::MatrixXd magneticFieldPattern(const TEMode& mode, double z, double t) const;
};

// Cylindrical waveguide
class CylindricalWaveguide {
    double radius_;

public:
    // Uses Bessel functions for field patterns
    struct TEMode {
        int m, n;  // Azimuthal and radial mode numbers
        double cutoff_frequency;
    };

    std::vector<TEMode> getTEModes() const;
};

// Resonating cavity
class RectangularCavity {
    double a_, b_, d_;  // Dimensions

public:
    struct ResonantMode {
        int m, n, p;  // Mode numbers
        double resonant_frequency;
        double quality_factor;
    };

    std::vector<ResonantMode> getResonantModes(double max_frequency) const;

    // Energy stored in mode
    double modalEnergy(const ResonantMode& mode, double amplitude) const;
};

} // namespace physics::waves
```

### 6. Phase Space & Liouville's Equation

**Approach:** State-space representation + OOP

```cpp
// include/physics/statistical/phase_space.hpp
#pragma once
#include <Eigen/Dense>
#include <functional>

namespace physics::statistical {

// Point in phase space (q, p)
struct PhasePoint {
    Eigen::VectorXd position;    // q
    Eigen::VectorXd momentum;    // p

    int dimension() const { return position.size(); }
};

// Phase space distribution function
class PhaseSpaceDistribution {
    int dim_;
    std::function<double(const PhasePoint&)> rho_;

public:
    PhaseSpaceDistribution(int dim,
                          const std::function<double(const PhasePoint&)>& rho)
        : dim_(dim), rho_(rho) {}

    double density(const PhasePoint& point) const { return rho_(point); }

    // Ensemble averages
    double averageEnergy(const std::function<double(PhasePoint)>& H) const;
    Eigen::VectorXd averageMomentum() const;
};

// Liouville's equation: dρ/dt = {H, ρ}
class LiouvilleEquation {
    int dim_;
    std::function<double(const PhasePoint&)> hamiltonian_;

public:
    LiouvilleEquation(int dim,
                     const std::function<double(const PhasePoint&)>& H)
        : dim_(dim), hamiltonian_(H) {}

    // Poisson bracket: {f, g} = Σ(∂f/∂q_i ∂g/∂p_i - ∂f/∂p_i ∂g/∂q_i)
    double poissonBracket(
        const std::function<double(PhasePoint)>& f,
        const std::function<double(PhasePoint)>& g,
        const PhasePoint& point) const;

    // Time evolution of distribution
    PhaseSpaceDistribution evolve(
        const PhaseSpaceDistribution& rho0, double dt, int steps) const;

    // Verify Liouville's theorem: phase space volume preservation
    bool verifyVolumePreservation(const PhasePoint& p0, double dt) const;
};

// Hamiltonian dynamics
class HamiltonianSystem {
    int dim_;
    std::function<double(const PhasePoint&)> H_;

public:
    HamiltonianSystem(int dim,
                     const std::function<double(const PhasePoint&)>& H)
        : dim_(dim), H_(H) {}

    // Hamilton's equations: dq/dt = ∂H/∂p, dp/dt = -∂H/∂q
    PhasePoint timeDerivative(const PhasePoint& point) const;

    // Integrate trajectory (symplectic integrator)
    std::vector<PhasePoint> trajectory(
        const PhasePoint& initial, double dt, int steps) const;

    // Energy conservation check
    bool checkEnergyConservation(const std::vector<PhasePoint>& traj) const;
};

} // namespace physics::statistical
```

### 7. Riemannian Geometry & General Relativity

**Approach:** Tensor calculus library

```cpp
// include/physics/relativity/metric_tensor.hpp
#pragma once
#include <Eigen/Dense>
#include <functional>

namespace physics::relativity {

using Tensor2 = Eigen::Matrix4d;
using Vector4 = Eigen::Vector4d;

// Metric tensor g_μν
class MetricTensor {
    std::function<Tensor2(const Vector4&)> g_;

public:
    MetricTensor(const std::function<Tensor2(const Vector4&)>& g) : g_(g) {}

    // Minkowski metric (flat spacetime)
    static MetricTensor minkowski();

    // Schwarzschild metric (black hole)
    static MetricTensor schwarzschild(double mass);

    // Kerr metric (rotating black hole)
    static MetricTensor kerr(double mass, double angular_momentum);

    // FRW metric (cosmology)
    static MetricTensor friedmann(double scale_factor, double curvature);

    // Metric components
    Tensor2 covariant(const Vector4& x) const { return g_(x); }
    Tensor2 contravariant(const Vector4& x) const { return g_(x).inverse(); }

    // Line element: ds² = g_μν dx^μ dx^ν
    double lineElement(const Vector4& x, const Vector4& dx) const;

    // Proper time
    double properTime(const Vector4& x, const Vector4& four_velocity) const;
};

// Christoffel symbols Γ^λ_μν
class ChristoffelSymbols {
    MetricTensor metric_;
    double epsilon_;  // For numerical derivatives

public:
    ChristoffelSymbols(const MetricTensor& g, double eps = 1e-6);

    // Γ^λ_μν = (1/2)g^λσ (∂_μ g_νσ + ∂_ν g_μσ - ∂_σ g_μν)
    double component(int lambda, int mu, int nu, const Vector4& x) const;

    // All components at point
    std::array<Tensor2, 4> all(const Vector4& x) const;
};

// Riemann curvature tensor R^ρ_σμν
class RiemannTensor {
    ChristoffelSymbols christoffel_;

public:
    RiemannTensor(const ChristoffelSymbols& gamma);

    // R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
    double component(int rho, int sigma, int mu, int nu, const Vector4& x) const;
};

// Ricci tensor R_μν = R^λ_μλν
class RicciTensor {
    RiemannTensor riemann_;

public:
    RicciTensor(const RiemannTensor& R);

    Tensor2 components(const Vector4& x) const;
    double trace(const Vector4& x, const MetricTensor& g) const;  // Ricci scalar
};

// Einstein tensor G_μν = R_μν - (1/2)g_μν R
class EinsteinTensor {
    RicciTensor ricci_;
    MetricTensor metric_;

public:
    EinsteinTensor(const RicciTensor& R, const MetricTensor& g);

    Tensor2 components(const Vector4& x) const;

    // Einstein field equations: G_μν = (8πG/c⁴) T_μν
    bool verifySolution(const Tensor2& stress_energy, const Vector4& x,
                       double G_newton, double c) const;
};

// Geodesic equation solver
class GeodesicSolver {
    MetricTensor metric_;
    ChristoffelSymbols christoffel_;

public:
    GeodesicSolver(const MetricTensor& g);

    // d²x^μ/dτ² + Γ^μ_αβ (dx^α/dτ)(dx^β/dτ) = 0
    struct GeodesicState {
        Vector4 position;
        Vector4 four_velocity;
    };

    GeodesicState derivative(const GeodesicState& state) const;

    // Integrate geodesic
    std::vector<GeodesicState> trajectory(
        const GeodesicState& initial, double dtau, int steps) const;

    // Planetary orbit (Schwarzschild)
    std::vector<GeodesicState> planetaryOrbit(
        double mass, double semi_major_axis, double eccentricity) const;

    // Perihelion shift
    double perihelionShift(double mass, double semi_major_axis,
                          double eccentricity, int orbits) const;
};

} // namespace physics::relativity
```

## Build System

### CMakeLists.txt (Extended)

```cmake
cmake_minimum_required(VERSION 3.15)
project(PhysicsShowcase VERSION 2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(BUILD_ADVANCED "Build advanced physics modules" ON)
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)

# Find dependencies
if(BUILD_ADVANCED)
    find_package(Eigen3 3.3 REQUIRED)
    find_package(Boost 1.70 COMPONENTS math REQUIRED)
    # find_package(GSL REQUIRED)  # Optional
endif()

# Header-only library (existing)
add_library(physics_basic INTERFACE)
target_include_directories(physics_basic INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Advanced library (new)
if(BUILD_ADVANCED)
    add_library(physics_advanced
        src/tensors/field_tensor.cpp
        src/tensors/stress_energy_tensor.cpp
        src/quantum/dirac_operator.cpp
        src/quantum/schrodinger_solver.cpp
        src/quantum/perturbation.cpp
        src/relativity/metric_tensor.cpp
        src/relativity/christoffel.cpp
        src/relativity/riemann_tensor.cpp
        src/relativity/geodesic_solver.cpp
        src/waves/green_functions.cpp
        src/waves/waveguide.cpp
        src/statistical/phase_space.cpp
        src/statistical/liouville.cpp
    )

    target_include_directories(physics_advanced PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )

    target_link_libraries(physics_advanced PUBLIC
        Eigen3::Eigen
        Boost::math
        physics_basic
    )

    # Compiler optimizations
    target_compile_options(physics_advanced PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-O3 -march=native -fopenmp>
        $<$<CXX_COMPILER_ID:Clang>:-O3 -march=native>
        $<$<CXX_COMPILER_ID:MSVC>:/O2 /openmp>
    )
endif()

# Tests
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Examples
add_subdirectory(examples)
```

## Dependencies

### Required Libraries

1. **Eigen** - Linear algebra, matrix operations
   ```bash
   sudo apt install libeigen3-dev  # Ubuntu/Debian
   brew install eigen              # macOS
   ```

2. **Boost** - Special functions (Bessel, etc.)
   ```bash
   sudo apt install libboost-all-dev
   ```

3. **GSL** (Optional) - GNU Scientific Library for numerical methods
   ```bash
   sudo apt install libgsl-dev
   ```

## Performance Considerations

### Optimization Strategies

1. **Template Metaprogramming** - Compile-time tensor operations
2. **Expression Templates** - Avoid temporary objects (Eigen does this)
3. **SIMD Vectorization** - Use Eigen's vectorization
4. **Parallel Algorithms** - OpenMP for large matrix operations
5. **GPU Acceleration** - CUDA/ROCm for intensive calculations (future)

### Memory Management

```cpp
// Use move semantics for large objects
ElectromagneticFieldTensor createField(const Eigen::Vector3d& E,
                                       const Eigen::Vector3d& B) {
    return ElectromagneticFieldTensor(E, B);  // RVO optimization
}

// Use const references to avoid copies
double computeInvariant(const ElectromagneticFieldTensor& F) {
    return F.firstInvariant();
}
```

## Summary: What Can Be Implemented

| Topic | Approach | Complexity | Dependencies |
|-------|----------|------------|--------------|
| Stress-Energy Tensor | OOP + Eigen | Medium | Eigen |
| Field Tensor | OOP + Eigen | Medium | Eigen |
| Dirac Formalism | Matrix classes | Medium | Eigen |
| Schrödinger Solver | Numerical PDE | Hard | Eigen, sparse solvers |
| Green Functions | Numerical methods | Hard | Boost (special functions) |
| Waveguides | Boundary value | Medium | Eigen, Boost |
| Phase Space | State space + ODE | Medium | Eigen |
| Liouville Equation | Numerical evolution | Hard | Eigen |
| Riemann Tensor | Symbolic/numerical | Hard | Eigen, automatic differentiation |
| Einstein Equations | Nonlinear PDE | Very Hard | Eigen, PETSc (optional) |
| Perturbation Theory | Iterative methods | Medium | Eigen |
| Boltzmann Transport | Monte Carlo | Hard | Random number generators |

## Next Steps

1. **Phase 1**: Implement tensor mathematics (stress-energy, field tensor)
2. **Phase 2**: Quantum operators (Dirac, spin matrices)
3. **Phase 3**: Numerical solvers (Schrödinger, wave equations)
4. **Phase 4**: General relativity (metric, geodesics)
5. **Phase 5**: Statistical mechanics (phase space, transport)

Each phase can be developed independently and added incrementally to the library.
