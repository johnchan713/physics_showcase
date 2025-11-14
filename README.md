# Mathematics & Physics Showcase - C++ Implementation

A comprehensive C++ library implementing fundamental mathematics and physics concepts with well-documented standalone functions. This showcase demonstrates mathematical foundations and physical principles through clean, reusable code organized in header-only libraries.

## Table of Contents

- [Mathematics](#mathematics)
  - [Calculus](#calculus)
  - [Trigonometry](#trigonometry)
  - [Linear Algebra](#linear-algebra)
  - [Probability & Statistics](#probability--statistics)
  - [Transforms](#transforms)
- [Physics](#physics)
  - [Basic Mechanics](#basic-mechanics)
  - [Advanced Physics](#advanced-physics)
- [Building and Running](#building-and-running)

---

# Mathematics

## Calculus

**Location**: `include/maths/calculus/`

### Fundamental Theorems (`theorems.hpp`)

#### Intermediate Value Theorem (IVT)
- `IntermediateValueTheorem::statement()` - Formal statement
- `checkConditions()` - Verify IVT applicability
- `findRoot()` - Find c such that f(c) = k (bisection method)
- `proofSqrt2Exists()` - Demonstrate √2 existence

#### Mean Value Theorem (MVT)
- `MeanValueTheorem::statement()` - Formal statement
- `averageRateOfChange()` - Calculate (f(b) - f(a))/(b - a)
- `findMVTPoint()` - Find c where f'(c) = average rate
- `geometricInterpretation()` - Tangent parallel to secant

#### Rolle's Theorem
- `RollesTheorem::statement()` - Special case of MVT
- `checkConditions()` - Verify f(a) = f(b)
- `findCriticalPoint()` - Find c where f'(c) = 0
- `relationshipToMVT()` - How Rolle's relates to MVT

#### Extreme Value Theorem (EVT)
- `ExtremeValueTheorem::statement()` - Guarantees max/min exist
- `findExtremes()` - Numerically find extreme values
- `importanceOfClosedInterval()` - Why [a,b] is essential

#### Fundamental Theorem of Calculus
- `firstFundamentalTheorem()` - d/dx[∫ₐˣ f(t)dt] = f(x)
- `secondFundamentalTheorem()` - ∫ₐᵇ f(x)dx = F(b) - F(a)
- `integrate()` - Numerical integration (trapezoidal rule)

#### L'Hôpital's Rule
- `LHopitalsRule::statement()` - For indeterminate forms
- `indeterminateForms()` - 0/0, ∞/∞, etc.
- Example applications

---

## Trigonometry

**Location**: `include/maths/trigonometry/`

### Trigonometric Functions (`identities.hpp`)

#### Basic Functions
- `sine()`, `cosine()`, `tangent()` - Primary trig functions
- `cotangent()`, `secant()`, `cosecant()` - Reciprocal functions
- `degreesToRadians()` / `radiansToDegrees()` - Angle conversions
- `unitCircleDefinitions()` - Definitions from unit circle
- `specialAngles()` - Exact values for 0°, 30°, 45°, 60°, 90°

#### Pythagorean Identities
- `fundamentalIdentity()` - sin²θ + cos²θ = 1
- `tangentIdentity()` - 1 + tan²θ = sec²θ
- `cotangentIdentity()` - 1 + cot²θ = csc²θ
- `verifyFundamental()` - Numerical verification

#### Angle Formulas
- `sineSum()` / `sineDifference()` - sin(α ± β)
- `cosineSum()` / `cosineDifference()` - cos(α ± β)
- `tangentSum()` / `tangentDifference()` - tan(α ± β)
- `doubleAngleFormulas()` - sin(2α), cos(2α), tan(2α)
- `halfAngleFormulas()` - sin(α/2), cos(α/2), tan(α/2)
- `tripleAngleFormulas()` - sin(3α), cos(3α), tan(3α)

#### Product-Sum Formulas
- `productToSum()` - Convert sin·sin, cos·cos to sums
- `sumToProduct()` - Convert sin+sin to products
- `powerReduction()` - sin²α, cos²α in terms of cos(2α)

#### Inverse Trigonometric Functions
- `arcsin()`, `arccos()`, `arctan()` - Inverse functions
- `arctan2()` - Two-argument arctangent
- `domainsAndRanges()` - Valid inputs/outputs
- `inverseIdentities()` - Relationships between inverses
- `derivatives()` - d/dx of inverse trig functions

#### Laws of Triangles
- `lawOfSines()` - a/sin(A) = b/sin(B) = c/sin(C)
- `lawOfCosines()` - c² = a² + b² - 2ab·cos(C)
- `solveWithSines()` - Solve triangle using law of sines
- `solveWithCosines()` - Solve triangle using law of cosines
- `areaFormulas()` - Triangle area (Heron's formula, etc.)

---

## Linear Algebra

**Location**: `include/maths/linear_algebra/`

### Vectors (`vectors.hpp`)

#### Vector Class
- `Vector(components)` - Construct vector
- `dimension()` - Get vector dimension
- `operator[]` - Access components
- `operator+`, `operator-` - Vector addition/subtraction
- `operator*` - Scalar multiplication

#### Vector Operations
- `dot()` - Dot product v·w
- `norm()` - Euclidean norm ‖v‖
- `normalize()` - Unit vector v̂ = v/‖v‖
- `distance()` - Distance between vectors
- `isZero()` - Check if zero vector

#### Cross Product (3D)
- `CrossProduct::compute()` - v × w
- `properties()` - Anticommutative, perpendicular
- `geometricMeaning()` - Parallelogram area

#### Projections
- `VectorProjection::scalarProjection()` - comp_w(v)
- `vectorProjection()` - proj_w(v)
- `orthogonalComponent()` - perp_w(v) = v - proj_w(v)
- `verifyDecomposition()` - v = proj + perp

#### Orthogonality
- `Orthogonality::areOrthogonal()` - Check v ⊥ w
- `isUnitVector()` - Check ‖v‖ = 1
- `areOrthonormal()` - Check pairwise orthogonal + unit
- `gramSchmidt()` - Orthogonalization process
- `gramSchmidtOrthonormal()` - Orthonormal basis
- `standardBasis()` - e₁, e₂, ..., eₙ

#### Linear Independence
- `LinearIndependence::areIndependent()` - Check linear independence
- `spanDimension()` - Dimension of span
- `extractIndependentSet()` - Find independent subset

### Matrices (`matrices.hpp`)

#### Matrix Class
- `Matrix(rows, cols)` - Construct m×n matrix
- `operator()` - Access elements
- `operator+`, `operator-` - Matrix addition/subtraction
- `operator*` (scalar) - Scalar multiplication
- `operator*` (matrix) - Matrix multiplication
- `operator*` (vector) - Matrix-vector product

#### Matrix Operations
- `transpose()` - Aᵀ
- `trace()` - Sum of diagonal elements
- `determinant()` - det(A) (Laplace expansion)
- `minor()` - Submatrix (remove row/column)
- `identity()` - Identity matrix I

#### Matrix Properties
- `isSquare()`, `isSymmetric()`, `isDiagonal()`, `isIdentity()`
- `determinantProperties()` - det(AB) = det(A)det(B), etc.
- `traceProperties()` - tr(AB) = tr(BA), etc.
- `transposeProperties()` - (AB)ᵀ = BᵀAᵀ, etc.
- `invertibilityConditions()` - When is A invertible?

#### Special Matrices
- `SpecialMatrices::diagonal()` - Diagonal matrix from vector
- `symmetricMatrices()` - A = Aᵀ properties
- `orthogonalMatrices()` - QᵀQ = I properties
- `positiveDefiniteMatrices()` - xᵀAx > 0 conditions
- `triangularMatrices()` - Upper/lower triangular

#### Matrix Decompositions
- `luDecomposition()` - A = LU
- `qrDecomposition()` - A = QR
- `svd()` - A = UΣVᵀ
- `eigenvalueDecomposition()` - A = PDP⁻¹
- `choleskyDecomposition()` - A = LLᵀ

#### Linear Transformations
- `LinearTransformations::definition()` - T(x) = Ax
- `kernel()` - Null space {x : Ax = 0}
- `range()` - Column space {Ax : x ∈ ℝⁿ}
- `rankNullityTheorem()` - rank + nullity = n

---

## Probability & Statistics

**Location**: `include/maths/probability/`

### Probability Distributions (`distributions.hpp`)

#### Discrete Distributions
- `BernoulliDistribution::pmf()` - P(X=x), x ∈ {0,1}
- `BinomialDistribution::pmf()` - C(n,k)p^k(1-p)^(n-k)
- `PoissonDistribution::pmf()` - λ^k e^(-λ) / k!
- `binomialCoefficient()` - C(n,k) = n!/(k!(n-k)!)

#### Continuous Distributions
- `UniformDistribution::pdf()` - 1/(b-a) on [a,b]
- `UniformDistribution::cdf()` - (x-a)/(b-a)
- `NormalDistribution::pdf()` - Gaussian distribution
- `NormalDistribution::cdf()` - Using error function
- `standardPDF()`, `standardCDF()` - Standard normal (μ=0, σ=1)
- `ExponentialDistribution::pdf()` - λe^(-λx)
- `ExponentialDistribution::cdf()` - 1 - e^(-λx)

#### Distribution Properties
- `mean()`, `variance()`, `stddev()` - For all distributions
- `empiricalRule()` - 68-95-99.7 rule (normal)
- `memorylessProperty()` - Exponential distribution

#### Statistical Functions
- `StatisticalFunctions::mean()` - Sample mean
- `variance()` - Sample variance (unbiased)
- `stddev()` - Sample standard deviation
- `median()` - Middle value
- `covariance()` - Cov(X,Y)
- `correlation()` - Pearson's r
- `zScore()` - Standardize value

#### Fundamental Theorems
- `centralLimitTheorem()` - X̄ → N(μ, σ²/n)
- `lawOfLargeNumbers()` - X̄ₙ → μ

#### Bayesian Inference
- `BayesianInference::bayesTheorem()` - P(A|B) = P(B|A)P(A)/P(B)
- `statement()` - Posterior, likelihood, prior
- `totalProbability()` - Law of total probability

---

## Transforms

**Location**: `include/maths/transforms/`

### Polar & Coordinate Transformations (`polar.hpp`)

#### 2D Polar Coordinates
- `PolarCoordinates::fromCartesian()` - (x,y) → (r,θ)
- `toCartesian()` - (r,θ) → (x,y)
- `distance()` - Distance between polar points
- `areaElement()` - dA = r dr dθ
- `jacobian()` - |∂(x,y)/∂(r,θ)| = r

#### Cylindrical Coordinates (3D)
- `CylindricalCoordinates::fromCartesian()` - (x,y,z) → (r,θ,z)
- `toCartesian()` - (r,θ,z) → (x,y,z)
- `volumeElement()` - dV = r dr dθ dz
- `jacobian()` - r

#### Spherical Coordinates (3D)
- `SphericalCoordinates::fromCartesian()` - (x,y,z) → (ρ,θ,φ)
- `toCartesian()` - (ρ,θ,φ) → (x,y,z)
- `volumeElement()` - dV = ρ² sin(φ) dρ dθ dφ
- `jacobian()` - ρ² sin(φ)
- `convention()` - Physics vs mathematics notation

#### Complex Numbers in Polar Form
- `ComplexPolarForm::toPolar()` - z → (r, θ)
- `fromPolar()` - (r, θ) → z = re^(iθ)
- `eulersFormula()` - e^(iθ) = cos θ + i sin θ
- `multiplicationRule()` - Multiply moduli, add arguments
- `deMoivresTheorem()` - (e^(iθ))^n = e^(inθ)
- `nthRoots()` - All nth roots of complex number
- `rootsOfUnity()` - nth roots of 1

#### Polar Curves
- `circle()`, `cardioid()`, `rose()`, `spiral()`, `lemniscate()`
- `areaFormula()` - A = ½∫r²dθ
- `arcLength()` - s = ∫√(r² + (dr/dθ)²) dθ

### Fourier Transforms (`fourier.hpp`)

#### Fourier Series
- `FourierSeries::definition()` - f(x) = a₀/2 + Σ[aₙcos(nx) + bₙsin(nx)]
- `complexForm()` - f(x) = Σ cₙe^(inx)
- `computeA0()`, `computeAn()`, `computeBn()` - Coefficients
- `dirichletConditions()` - Convergence conditions
- `parsevalTheorem()` - Energy conservation

#### Fourier Transform
- `FourierTransform::definition()` - F(ω) = ∫f(t)e^(-iωt)dt
- `properties()` - Linearity, shift, scaling, convolution
- `commonPairs()` - δ(t) ↔ 1, e^(-at) ↔ 1/(a+iω), etc.
- `uncertaintyPrinciple()` - Δt·Δω ≥ 1/2
- `convolutionTheorem()` - f*g ↔ F·G

#### Discrete Fourier Transform (DFT)
- `DiscreteFourierTransform::definition()` - X[k] = Σx[n]e^(-i2πkn/N)
- `compute()` - Naive O(N²) DFT
- `inverse()` - Inverse DFT
- `magnitudeSpectrum()`, `phaseSpectrum()`
- `nyquistTheorem()` - f_sample ≥ 2f_max

#### Fast Fourier Transform (FFT)
- `FastFourierTransform::algorithm()` - O(N log N) Cooley-Tukey
- `compute()` - Recursive FFT (power-of-2 sizes)
- `applications()` - Audio, images, communications

#### Window Functions
- `WindowFunctions::rectangular()`, `hann()`, `hamming()`
- `purpose()` - Reduce spectral leakage

---

# Physics

## Basic Mechanics

**Location**: `include/physics/`

The physics library contains 31 comprehensive modules covering classical mechanics, thermodynamics, electromagnetism, and modern physics. All functions use SI units and include detailed documentation.

### Core Modules

1. **Newton's Laws** (`newton_laws.hpp`)
   - First Law (inertia), Second Law (F=ma), Third Law (action-reaction)

2. **Kinematics** (`kinematics.hpp`)
   - Motion with constant acceleration, all kinematic equations

3. **Dynamics** (`dynamics.hpp`)
   - Force-motion integration, friction, work, power

4. **Energy & Momentum** (`energy_momentum.hpp`)
   - KE, PE, momentum, impulse, conservation laws

5. **Projectile Motion** (`projectile.hpp`)
   - 2D parabolic trajectories

6. **Circular Motion** (`circular_motion.hpp`)
   - Centripetal force, angular velocity, period

7. **Simple Harmonic Motion** (`harmonic_motion.hpp`)
   - SHM, pendulums, vibrating masses, energy

8. **Rotational Dynamics** (`rotational_dynamics.hpp`)
   - Torque, angular momentum, moment of inertia

9. **Gravitation** (`gravitation.hpp`)
   - Universal gravitation, Kepler's laws

10. **Orbital Mechanics** (`orbital.hpp`)
    - Satellite motion, escape velocity

### Thermodynamics

11. **Thermodynamics** (`thermodynamics.hpp`)
    - Gas laws, kinetic theory, RMS velocity

12. **Thermal Expansion** (`thermal_expansion.hpp`)
    - Linear, area, volume expansion

13. **Calorimetry** (`calorimetry.hpp`)
    - Heat capacity, phase changes

14. **Heat Transfer** (`heat_transfer.hpp`)
    - Conduction, radiation, heat engines, Carnot cycle

### Fluids and Materials

15. **Fluid Mechanics** (`fluid_mechanics.hpp`)
    - Continuity, Bernoulli, pipe friction

16. **Elasticity** (`elasticity.hpp`)
    - Young's modulus, beam bending

17. **Surface Tension** (`surface_tension.hpp`)
    - Capillary rise, droplet/bubble pressure

### Waves and Optics

18. **Wave Mechanics** (`wave_mechanics.hpp`)
    - Sound, Doppler effect, string vibrations

19. **Electromagnetic Waves** (`electromagnetic_waves.hpp`)
    - EM wave propagation, radiation pressure

20. **Optics** (`optics.hpp`)
    - Refraction, lenses, mirrors, telescopes

### Electricity and Magnetism

21. **Electrostatics** (`electrostatics.hpp`)
    - Coulomb's law, electric field, capacitance

22. **Magnetism** (`magnetism.hpp`)
    - Magnetic forces and fields

23. **Electric Circuits** (`electric_circuits.hpp`)
    - Ohm's law, resistance, power

24. **Electromagnetic Induction** (`electromagnetic_induction.hpp`)
    - Faraday's law, motors, generators, transformers

25. **Maxwell Equations** (`maxwell_equations.hpp`)
    - Maxwell's four equations, Poynting vector

### Advanced Topics

26. **Advanced Mechanics** (`advanced_mechanics.hpp`)
    - Lagrangian, Hamiltonian mechanics

27. **Oscillations** (`oscillations.hpp`)
    - Damped, forced, coupled oscillators, resonance

28. **Special Relativity** (`special_relativity.hpp`)
    - Lorentz transformations, relativistic mechanics

29. **Quantum Basics** (`quantum_basics.hpp`)
    - de Broglie, Compton, uncertainty, Bohr model

30. **Unit Conversions** (`units.hpp`)
    - CGS, SI, Imperial conversions

31. **Inclined Plane** (`inclined_plane.hpp`)
    - Forces and motion on slopes

## Advanced Physics

**Location**: `include/physics/advanced/`

### Particle Physics & QFT (`qft/`)
- Standard Model particles (quarks, leptons, bosons)
- Spin-statistics theorem
- Supersymmetry (SUSY, MSSM)
- Antiparticles and CP violation
- Quark-gluon plasma
- Interaction cross sections
- Decay rates and resonances

### Cosmology (`cosmology/`)
- Hubble expansion and redshift
- Friedmann equations
- Energy density components (matter, dark energy, radiation)
- Early universe (CMB, BBN, baryogenesis)
- Thermal history

### Gauge Theory (`gauge_theory/`)
- Discrete symmetries (P, C, T, CPT)
- Helicity and chirality
- Gauge transformations (U(1), SU(2), SU(3))
- Higgs mechanism
- Running couplings
- CP violation in kaons and B mesons

### Fluid Dynamics (`fluid_dynamics/`)
- Governing equations (continuity, Euler, Navier-Stokes)
- Dimensionless numbers (Reynolds, Mach, Froude)
- Vorticity and circulation
- Turbulence models

---

## Project Structure

```
maths_physics_showcase/
├── include/
│   ├── maths/
│   │   ├── calculus/
│   │   │   └── theorems.hpp
│   │   ├── trigonometry/
│   │   │   └── identities.hpp
│   │   ├── linear_algebra/
│   │   │   ├── vectors.hpp
│   │   │   └── matrices.hpp
│   │   ├── probability/
│   │   │   └── distributions.hpp
│   │   └── transforms/
│   │       ├── polar.hpp
│   │       └── fourier.hpp
│   │
│   └── physics/
│       ├── [31 basic physics modules]
│       └── advanced/
│           ├── qft/          # Particle physics
│           ├── cosmology/    # Cosmology
│           ├── gauge_theory/ # Gauge theories
│           └── fluid_dynamics/ # Advanced fluids
│
├── examples/
│   ├── main.cpp
│   ├── advanced_demo.cpp
│   └── scientific_demo.cpp
│
├── Makefile
├── CMakeLists.txt
└── README.md
```

---

## Building and Running

### Prerequisites
- C++ compiler with C++11 support or later (g++, clang++, MSVC)
- Make (optional)
- CMake 3.10+ (optional)

### Option 1: Direct Compilation

```bash
# Compile with math and physics headers
g++ -std=c++11 -I./include your_program.cpp -o program
./program
```

### Option 2: Using Make

```bash
make
make run
```

### Option 3: Using CMake

```bash
mkdir build && cd build
cmake ..
make
```

---

## Usage Examples

### Mathematics Examples

#### Calculus - Mean Value Theorem
```cpp
#include "maths/calculus/theorems.hpp"

auto f = [](double x) { return x * x; };  // f(x) = x²
auto df = [](double x) { return 2.0 * x; }; // f'(x) = 2x

double c = maths::calculus::MeanValueTheorem::findMVTPoint(f, df, 0.0, 2.0);
// Result: c = 1.0, where f'(c) = average slope
```

#### Linear Algebra - Vector Operations
```cpp
#include "maths/linear_algebra/vectors.hpp"

using namespace maths::linear_algebra;

Vector v({3.0, 4.0});
Vector w({1.0, 2.0});

double dot_product = v.dot(w);  // Result: 11.0
double v_norm = v.norm();        // Result: 5.0
Vector v_hat = v.normalize();    // Unit vector

Vector cross = CrossProduct::compute(v, w);  // 3D only
```

#### Fourier Transform - FFT
```cpp
#include "maths/transforms/fourier.hpp"

std::vector<std::complex<double>> signal = {...};
auto spectrum = maths::transforms::FastFourierTransform::compute(signal);
auto magnitude = maths::transforms::DiscreteFourierTransform::magnitudeSpectrum(spectrum);
```

### Physics Examples

#### Newton's Second Law
```cpp
#include "physics/newton_laws.hpp"

double force = physics::newton::calculateForce(10.0, 5.0);
// F = ma = 10 kg × 5 m/s² = 50 N
```

#### Projectile Motion
```cpp
#include "physics/projectile.hpp"

double range = physics::projectile::calculateRange(20.0, 45.0 * M_PI/180.0, 9.81);
// 45° gives maximum range
```

---

## Design Principles

1. **Header-Only**: All implementations in headers for easy integration
2. **Namespace Organization**: `maths::` and `physics::` with subnamespaces
3. **Comprehensive Documentation**: Doxygen-style comments
4. **Error Handling**: Input validation with exceptions
5. **SI Units**: Consistent unit system (physics)
6. **Pure Functions**: No side effects
7. **Modern C++**: C++11/14/17 features where appropriate

---

## License

This project is licensed under the terms specified in the LICENSE file.

---

## References

**Mathematics:**
- Calculus (Stewart, Spivak)
- Linear Algebra (Strang, Axler)
- Probability Theory (Ross, Feller)
- Fourier Analysis (Bracewell, Oppenheim & Schafer)

**Physics:**
- Classical Mechanics (Goldstein, Marion & Thornton)
- Electromagnetism (Griffiths, Jackson)
- Quantum Mechanics (Griffiths, Sakurai)
- Particle Physics (Griffiths, Peskin & Schroeder)
- Cosmology (Weinberg, Dodelson)

---

## Statistics

- **Total Files**: 60+ header files
- **Total Lines**: ~31,000+ lines of code
- **Mathematics Modules**: 7 modules
- **Physics Modules**: 31 basic + 4 advanced categories
- **Zero Dependencies**: Header-only, no external libraries required

---

*This showcase demonstrates both the mathematical foundations and physical applications, providing a comprehensive toolkit for scientific computing in C++.*
