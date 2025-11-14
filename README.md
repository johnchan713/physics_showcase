# Mathematics & Physics Computational Showcase

Comprehensive C++17 header-only library implementing computational algorithms from mathematics, physics, and probability theory. All implementations follow a **computational pattern**: concrete parameters, numerical results, no educational strings.

## 📚 Table of Contents

- [Project Structure](#project-structure)
- [Mathematics Modules](#mathematics-modules)
- [Physics Modules](#physics-modules)
- [Demo Programs](#demo-programs)
- [Building and Running](#building-and-running)
- [Features](#features)

## 🏗️ Project Structure

```
physics_showcase/
├── include/
│   ├── maths/
│   │   ├── basic/                 # Foundational mathematics
│   │   │   ├── calculus/          # Calculus theorems and methods
│   │   │   ├── linear_algebra/    # Matrix operations & vectors
│   │   │   ├── trigonometry/      # Trig functions and identities
│   │   │   └── transforms/        # Basic transforms (Fourier, polar)
│   │   ├── advanced/              # Advanced mathematical topics
│   │   │   ├── algebra/           # Differential Algebra (Ritt-Kolchin)
│   │   │   ├── analysis/          # Fourier analysis & subdifferentials
│   │   │   ├── geometry/          # Variational calculus & Lagrangians
│   │   │   ├── stochastic/        # Monte Carlo & MCMC methods
│   │   │   ├── probability/       # Probability distributions & statistics
│   │   │   ├── dynamical_systems/ # ODEs, chaos theory, bifurcations
│   │   │   └── pde/              # Partial differential equations
│   │   ├── finance/               # Financial mathematics
│   │   ├── actuarial/             # Actuarial science
│   │   └── econometrics/          # Econometric models
│   └── physics/
│       ├── (basic modules)        # Classical mechanics, waves, etc.
│       └── advanced/              # Advanced physics topics
│           ├── classical/         # Hamiltonian, Liouville, phase space
│           ├── cosmology/         # Friedmann equations, early universe
│           ├── fluid_dynamics/    # Turbulence, compressible flow
│           ├── gauge_theory/      # Gauge invariance, Higgs mechanism
│           └── qft/              # Quantum field theory
├── examples/                      # Demonstration programs
└── README.md
```

## 🧮 Mathematics Modules

### Advanced Mathematics

### Differential Algebra (`maths/advanced/algebra/differential_algebra.hpp`)
**Chapters I-IX from Ritt's "Differential Algebra"**

- **Differential Polynomials & Fields (Ch. I)**
  - Polynomial rings with derivation operator D
  - Leibniz rule: D(fg) = D(f)g + fD(g)
  - Order and degree computation

- **Characteristic Sets & Reduction (Ch. I, V)**
  - Ritt's algorithm for reduction
  - Orderly and eliminative rankings
  - Leader and initial extraction
  - Triangular form maintenance

- **Differential Ideals (Ch. I)**
  - Ideal generation with closure under differentiation
  - Basis computation
  - Membership testing via reduction

- **Algebraic Differential Manifolds (Ch. II)**
  - Solution manifolds of differential systems
  - Dimension = number of arbitrary constants
  - Irreducible component decomposition
  - Generic zeros and theorem of zeros

- **Manifold Intersections (Ch. VII)**
  - Intersection computation M₁ ∩ M₂
  - Kronecker's theorem analogue for dimensions
  - General solution intersection analysis

- **Orthonomic Systems (Ch. VIII - Riquier)**
  - Orthonomic form construction
  - Passive systems (closed under differentiation)
  - Taylor series dissection
  - Riquier's existence theorem

- **Partial Differential Algebra (Ch. IX)**
  - PDEs with multiple independent variables
  - Laplace, wave, heat equations
  - Cauchy-Kovalevskaya criterion
  - Characteristic sets for PDE systems
  - Low power theorem for singular solutions

**Example Applications:**
- Harmonic oscillator: y'' + y = 0
- Exponential growth: y' - y = 0
- Pendulum: θ'' + (g/L)θ = 0
- Laplace equation: ∇²u = 0
- Wave equation: u_tt - c²u_xx = 0

### Fourier Analysis (`maths/advanced/analysis/fourier_analysis.hpp`)
**Discrete & continuous Fourier theory**

- **Discrete & Fast Fourier Transform**
  - DFT: O(N²) naive algorithm
  - FFT: O(N log N) Cooley-Tukey radix-2
  - 2D FFT for image processing
  - Power spectrum, magnitude, phase

- **Circulant Matrices**
  - Diagonalization by Fourier basis
  - Eigenvalues via DFT
  - Fast multiplication O(N log N)

- **Convolution**
  - Convolution theorem: f * g = IFFT(FFT(f) · FFT(g))
  - Circular and linear convolution
  - Cross-correlation, auto-correlation

- **Wavelet Transforms**
  - Haar wavelet (simplest orthogonal)
  - Daubechies-4 with D4 coefficients
  - Perfect reconstruction

- **Fourier Series & Multipliers**
  - Coefficients for periodic functions on S¹
  - Spectral differentiation via (ik)
  - Fractional Laplacian (-Δ)^s

- **Time-Frequency Analysis**
  - Short-Time Fourier Transform (STFT)
  - Spectrograms |STFT|²
  - Gabor transform with Gaussian windows
  - Chirp signal detection

### Advanced Subdifferentials (`maths/advanced/analysis/advanced_subdifferentials.hpp`)
**Nonsmooth analysis & variational calculus**

- **Clarke Subdifferential** ∂_C f(x)
- **Limiting (Mordukhovich) Subdifferential** ∂_L f(x)
- **Fréchet, Proximal, Graded Subdifferentials**
- **Normal & Tangent Cones**: N_C(x), T_C(x)
- **Coderivatives** D*F for set-valued mappings
- **Calculus Rules**: sum, chain, max rules
- **Metric Regularity Criterion**

### Nonsmooth Algorithms (`maths/advanced/analysis/nonsmooth_algorithms.hpp`)
**Optimization algorithms**

- **Proximal Operators**
  - Soft thresholding for L1
  - Box constraints
  - Quadratic penalties

- **Subgradient Methods**
  - Subgradient descent
  - Diminishing step sizes

- **ISTA & FISTA**
  - Iterative Shrinkage-Thresholding (O(1/k))
  - Fast ISTA with Nesterov acceleration (O(1/k²))

- **ADMM**
  - Alternating Direction Method of Multipliers
  - Consensus optimization
  - Dual variable updates

### Stochastic Methods (`maths/advanced/stochastic/monte_carlo.hpp`)
**Monte Carlo methods and stochastic algorithms**

- **Monte Carlo Integration**
  - Basic Monte Carlo integration
  - Importance sampling and rejection sampling
  - Stratified sampling for variance reduction

- **Markov Chains**
  - Discrete-time Markov chains with transition matrices
  - Stationary distributions via power iteration
  - n-step transitions and irreducibility checking

- **MCMC Sampling**
  - Metropolis-Hastings algorithm (general and symmetric)
  - Gibbs sampling for multivariate distributions
  - Convergence diagnostics (ESS, Gelman-Rubin R̂, acceptance rate)

- **Stochastic Processes**
  - Standard Brownian motion (Wiener process)
  - Geometric Brownian motion (stock price model)
  - Ornstein-Uhlenbeck process (mean-reverting)

- **Boltzmann Equation and Kinetic Theory**
  - Maxwell-Boltzmann velocity distribution
  - Direct Simulation Monte Carlo (DSMC)
  - H-theorem for entropy evolution

- **Hamiltonian Monte Carlo**
  - HMC with leapfrog integration
  - Efficient sampling without random walk
  - Low autocorrelation chains

**Applications:** Finance (option pricing), physics (molecular dynamics), Bayesian inference, machine learning

### Variational Calculus (`maths/advanced/geometry/variational_calculus.hpp`)
**Lagrangians, Poincaré-Cartan Forms, and variational principles**

- **Contact Geometry**
  - Contact structures on jet bundles J^1(R,R)
  - Contact forms θ = du - p dx
  - Reeb vector fields and Legendre submanifolds

- **Lagrangians and Euler-Lagrange Equations**
  - Lagrangian functions L(x, u, u_x, ...)
  - Euler-Lagrange operator: ∂L/∂u - D(∂L/∂u_x) + ...
  - Total derivatives and jet prolongation

- **Poincaré-Cartan Forms**
  - Poincaré-Cartan form θ_L for Lagrangians
  - Canonical symplectic forms Ω = dθ_L
  - Cartan's integral invariants

- **Legendre Transformations**
  - Fiber derivative: (x, u, u_x) ↦ (x, u, p)
  - Hamiltonian: H = p·u_x - L
  - Lagrangian ↔ Hamiltonian mechanics

- **Noether's Theorem**
  - Symmetries → Conservation laws
  - Time translation → Energy conservation
  - Space translation → Momentum conservation
  - Infinitesimal symmetries and prolongations

- **Advanced Topics**
  - Second variation and Jacobi equations
  - Conjugate points and optimality
  - Bäcklund transformations (Sine-Gordon, KdV)
  - Conservation laws for PDEs
  - Field theories (Klein-Gordon, Maxwell, Yang-Mills)
  - De Donder-Weyl Hamiltonian formulation

**Applications:** Classical mechanics, field theory, optimal control, integrable systems

### Differential Equations and Dynamical Systems (`maths/advanced/dynamical_systems/ode_dynamical_systems.hpp`)
**Comprehensive ODE theory and chaos**

**Classical ODE Theory:**
- **Newton's Equations**: Second-order to first-order conversion, autonomous equations, equilibria
- **Initial Value Problems**: Euler, Heun, RK4 methods, Picard iteration, Lipschitz continuity
- **Linear Systems**: Matrix exponential, fundamental matrices, Floquet theory for periodic systems
- **Complex Domain**: Frobenius method, indicial equations, Bessel's equation
- **Boundary Value Problems**: Sturm-Liouville theory foundation

**Dynamical Systems:**
- **Flows and Trajectories**: Flows φ_t(x), fixed points, Liapunov functions, stability analysis
- **Local Behavior**: Jacobian analysis, eigenvalues, classification (nodes, saddles, spirals, centers)
- **Linearization**: Hartman-Grobman theorem for hyperbolic fixed points
- **Planar Systems**: Poincaré-Bendixson theorem foundation, limit cycles
- **Higher Dimensions**: Attractors, Lorenz system, Hamiltonian mechanics, KAM theorem

**Chaos Theory:**
- **Discrete Systems**: Logistic map, period doubling, bifurcation diagrams
- **Lyapunov Exponents**: λ > 0 ⟹ chaos, numerical computation
- **Poincaré Maps**: First return maps, periodic orbits
- **Homoclinic Chaos**: Melnikov method for chaos detection
- **Period Theory**: Sarkovskii's theorem (period 3 implies all periods)
- **Symbolic Dynamics**: Orbit encoding, admissible sequences
- **Fractals**: Box-counting dimension, strange attractors
- **Topological Chaos**: Smale horseshoe, stretch and fold mechanisms

**Applications**: Physics (pendulum, Lorenz), biology (population dynamics), engineering (nonlinear control)

### Partial Differential Equations (`maths/advanced/pde/partial_differential_equations.hpp`)
**Classical PDE theory and method of characteristics**

**PDE Classification and Fundamentals:**
- **Order and Linearity**: First/second/higher order, linear/quasi-linear/semi-linear/fully nonlinear
- **Second Order Types**: Elliptic (Δ < 0), parabolic (Δ = 0), hyperbolic (Δ > 0) via discriminant Δ = B² - AC
- **Boundary Conditions**: Dirichlet (u = g), Neumann (∂u/∂n = g), Robin (αu + β∂u/∂n = g), Cauchy
- **Superposition Principle**: Linear combinations for linear PDEs, solution space structure

**Well-Known PDEs:**
- **Heat Equation**: u_t = α u_xx (fundamental solution, diffusion, smoothing)
- **Wave Equation**: u_tt = c² u_xx (d'Alembert solution, propagation)
- **Laplace Equation**: Δu = 0 (harmonic functions, mean value property)
- **Poisson Equation**: Δu = f (with source term)
- **Transport Equation**: u_t + c·∇u = 0 (advection)

**Method of Characteristics:**
- **First Order Linear**: Constant/variable coefficients, characteristic curves dy/dx = b/a
- **Quasi-Linear Equations**: a(x,y,u) u_x + b(x,y,u) u_y = c(x,y,u), Charpit's method
- **Fully Nonlinear Equations**: F(x, y, u, u_x, u_y) = 0, complete Charpit system
- **Geometrical Interpretation**: Integral surfaces, Monge cones, characteristic directions
- **Second Order Characteristics**: A(dy)² - 2B dx dy + C(dx)² = 0, canonical forms

**Key Algorithms**: Classification via discriminant, characteristic ODE integration (Euler, RK4), Charpit solver, solution verification, boundary condition checking

**Applications**: Heat diffusion, wave propagation, fluid mechanics, electrostatics, quantum mechanics, optimal control

### PDE Solution Methods (`maths/advanced/pde/pde_solution_methods.hpp`)
**Classical solution techniques for PDEs**

**Linear Equations with Constant Coefficients:**
- **Inverse Operators**: Differential operator D = d/dx, inverse operator D⁻¹ (integration)
- **Polynomial Operators**: P(D) = aₙDⁿ + ... + a₁D + a₀, factorization methods
- **Exponential Shift**: Operator shift formula e^(ax) P(D) e^(-ax) = P(D - a)
- **Homogeneous Equations**: P(D)u = 0, complementary function from characteristic equation
- **Nonhomogeneous Equations**: P(D)u = f(x), particular solutions via inverse operators
- **Solution Structure**: General solution = complementary function + particular solution

**Orthogonal Expansions:**
- **Orthogonal Polynomials**:
  - Legendre polynomials Pₙ(x) on [-1,1]: (n+1)Pₙ₊₁ = (2n+1)xPₙ - nPₙ₋₁
  - Chebyshev polynomials Tₙ(x) = cos(n arccos x), minimal deviation property
  - Hermite polynomials Hₙ(x): Hₙ₊₁ = 2xHₙ - 2nHₙ₋₁, weight function e^(-x²)
  - Laguerre polynomials Lₙ(x): (n+1)Lₙ₊₁ = (2n+1-x)Lₙ - nLₙ₋₁, weight e^(-x)
- **Fourier Series Expansions**:
  - Trigonometric series: f(x) = a₀/2 + ∑ aₙcos(nπx/L) + bₙsin(nπx/L)
  - Half-range expansions: sine series, cosine series
  - Convergence theorems for piecewise smooth functions
- **Bessel Functions**:
  - Bessel functions of first kind Jₙ(x): series expansions, recurrence relations
  - Modified Bessel functions Iₙ(x) for imaginary arguments
  - Zeros of Bessel functions for eigenvalue problems
  - Applications to cylindrical boundary value problems

**Separation of Variables:**
- **Wave Equation Solutions** (hyperbolic): u_tt = c²u_xx
  - Series form: u(x,t) = ∑ (Aₙcos(ωₙt) + Bₙsin(ωₙt))sin(nπx/L)
  - Standing waves, normal modes, frequency spectrum
- **Heat Equation Solutions** (parabolic): u_t = αu_xx
  - Exponential decay: u(x,t) = ∑ Aₙ exp(-α(nπ/L)²t)sin(nπx/L)
  - Long-time behavior and steady states
- **Laplace Equation Solutions** (elliptic): Δu = 0
  - Rectangular domains: product solutions X(x)Y(y)
  - Dirichlet boundary value problems
- **Cylindrical Coordinate Systems**: ∇²u = u_rr + (1/r)u_r + (1/r²)u_θθ + u_zz
  - Bessel function solutions for radial dependence
  - Eigenvalues from zeros of Jₙ(x)
- **Spherical Coordinate Systems**: Laplacian in (r,θ,φ)
  - Legendre polynomial solutions for angular dependence
  - Azimuthal symmetry problems
- **Nonhomogeneous Problems**: Eigenfunction expansions, Duhamel's principle for time-dependent sources

**Key Techniques**: Orthogonality relations, eigenfunction expansions, Fourier coefficient computation, separation ansatz u(x,t) = X(x)T(t), boundary condition matching, series convergence analysis

**Applications**: Vibrating strings, heat conduction, electrostatic potential, quantum mechanics (particle in box), acoustics, diffusion processes

### PDE Transform Methods (`maths/advanced/pde/pde_transform_methods.hpp`)
**Laplace and Fourier transforms for solving PDEs**

**Laplace Transforms:**
- **Definition and Notation**: L{f(t)} = F(s) = ∫₀^∞ e^(-st) f(t) dt
- **Transform Pairs**: Exponentials, polynomials, trigonometric functions, hyperbolic functions
- **Properties**:
  - Linearity: L{af + bg} = aL{f} + bL{g}
  - First shifting theorem: L{e^(at)f(t)} = F(s-a)
  - Second shifting theorem (time delay): L{f(t-a)u(t-a)} = e^(-as)F(s)
  - Transform of derivatives: L{f'(t)} = sF(s) - f(0), L{f''(t)} = s²F(s) - sf(0) - f'(0)
  - Transform of integrals: L{∫₀ᵗ f(τ)dτ} = F(s)/s
- **Convolution Theorem**: L{f * g} = F(s)G(s)
- **Inverse Transform**: Partial fraction decomposition, residue method

**Fourier Transforms:**
- **Fourier Integral**: F{f(x)} = F(k) = ∫₋∞^∞ f(x) e^(-ikx) dx
- **Inverse Transform**: f(x) = (1/2π) ∫₋∞^∞ F(k) e^(ikx) dk
- **Transform Pairs**: Gaussian, rectangular pulse, Dirac delta, double exponential
- **Properties**:
  - Linearity, time shifting, frequency shifting, scaling
  - Differentiation: F{f'(x)} = ikF(k), F{f^(n)(x)} = (ik)^n F(k)
  - Multiplication by x: F{xf(x)} = iF'(k)
- **Parseval's Theorem**: ∫ |f(x)|² dx = (1/2π) ∫ |F(k)|² dk (energy conservation)
- **Convolution Theorem**: F{f * g} = F{f} · F{g}

**Fourier Sine and Cosine Transforms:**
- **Sine Transform**: Fs{f(x)} = ∫₀^∞ f(x) sin(kx) dx for odd extensions
- **Cosine Transform**: Fc{f(x)} = ∫₀^∞ f(x) cos(kx) dx for even extensions
- **Inverse Transforms**: f(x) = (2/π) ∫₀^∞ Fs(k) sin(kx) dk
- **Derivative Properties**: Fs{f''(x)} = -k²Fs{f(x)} - kf(0)

**Finite Fourier Transforms:**
- **Finite Sine Transform**: Fsn = ∫₀^L f(x) sin(nπx/L) dx
- **Finite Cosine Transform**: Fcn = ∫₀^L f(x) cos(nπx/L) dx
- **Applications**: Heat equation on finite intervals, boundary value problems

**Applications**: Transform methods for ODEs, heat equation, wave equation, diffusion problems, signal processing

### PDE Classification Solutions (`maths/advanced/pde/pde_classification_solutions.hpp`)
**Detailed solutions for parabolic, elliptic, and hyperbolic PDEs**

**Parabolic Equations (Heat/Diffusion):**
- **Heat Equation**: u_t = α u_xx (one-dimensional diffusion)
- **Fundamental Solution**: Heat kernel G(x,t;ξ) = 1/√(4παt) exp(-(x-ξ)²/4αt)
- **Infinite Domain Solutions**: Convolution with initial data
- **Bounded Domain Solutions**:
  - Dirichlet BC: u(x,t) = ∑ Aₙ exp(-α(nπ/L)²t) sin(nπx/L)
  - Neumann BC: u(x,t) = A₀ + ∑ Aₙ exp(-α(nπ/L)²t) cos(nπx/L)
- **Maximum Principles**: Weak and strong maximum principles
- **2D Heat Equation**: Rectangular and circular domains
- **Properties**: Infinite speed of propagation, smoothing effect, irreversibility

**Elliptic Equations (Laplace/Poisson):**
- **Laplace Equation**: Δu = 0 (harmonic functions)
- **Poisson Equation**: Δu = f (with source term)
- **Mean Value Property**: u(x₀,y₀) = (1/2πr) ∫ u on circle
- **Maximum Principles**: Maximum and minimum on boundary
- **Green's Functions**: G(x,y;ξ,η) = -(1/2π) ln(r) for 2D unbounded domain
- **Rectangular Domains**: Separation of variables with sinh/cosh solutions
- **Circular Domains**: Poisson integral formula
- **Harmonic Functions**: Solutions satisfy mean value property
- **Properties**: No time evolution, boundary value problems, smoothness

**Hyperbolic Equations (Wave):**
- **Wave Equation**: u_tt = c² u_xx (one-dimensional)
- **d'Alembert's Solution**: u(x,t) = ½[f(x+ct) + f(x-ct)] + 1/(2c) ∫ g(s) ds
- **Domain of Dependence**: Solution at (x,t) depends only on [x-ct, x+ct]
- **Standing Waves**: u(x,t) = ∑ (Aₙcos(ωₙt) + Bₙsin(ωₙt))sin(nπx/L)
- **Energy Conservation**: E = ½∫[u_t² + c²u_x²]dx is constant
- **2D Wave Equation**: Rectangular domains, eigenfrequencies ωₘₙ = c√(λₘ² + μₙ²)
- **Characteristic Cones**: Causality and light cones in spacetime
- **Finite Speed of Propagation**: Disturbances travel at speed c
- **Properties**: Reversible, energy conserving, finite propagation speed

**Key Concepts**: Well-posedness, uniqueness, regularity, stability, physical interpretation

**Applications**: Heat conduction, diffusion processes, electrostatics, membrane vibrations, acoustic waves, electromagnetic waves

### Probability & Statistics (`maths/advanced/probability/distributions.hpp`)
**Comprehensive probability distributions**

**Discrete Distributions:**
- **Bernoulli**: P(X=1) = p
- **Binomial**: C(n,k) p^k (1-p)^(n-k)
- **Poisson**: λ^k e^(-λ) / k!
- **Geometric**: Trials until first success
- **Negative Binomial**: Failures before r successes, overdispersion modeling
- **Hypergeometric**: Sampling without replacement, finite population correction

**Continuous Distributions:**
- **Uniform**: constant density on [a,b]
- **Normal (Gaussian)**: N(μ, σ²) with PDF, CDF, quantile
- **Exponential**: memoryless, rate λ
- **Gamma**: shape α, rate β
- **Beta**: on [0,1], conjugate prior
- **Chi-squared**: χ²(k) for hypothesis testing
- **F-Distribution**: Ratio of chi-squared, ANOVA, regression F-tests

**Statistical Functions:**
- PMF, PDF, CDF for all distributions
- Quantile functions (inverse CDF)
- Mean, variance, standard deviation
- Sampling with std::random
- Error function erf(x)
- Gamma function Γ(x)

**Statistical Tests:**
- One-sample t-test
- Chi-squared goodness-of-fit
- Maximum likelihood estimation

### Basic Mathematics

- **Calculus** (`maths/basic/calculus/theorems.hpp`): Numerical derivatives, integration (Simpson's rule)
- **Trigonometry** (`maths/basic/trigonometry/identities.hpp`): Computational trig identities
- **Linear Algebra** (`maths/basic/linear_algebra/`): Matrix operations, vectors, eigenvalues
- **Transforms** (`maths/basic/transforms/`): Fourier transforms, polar coordinates

### Applied Mathematics

- **Financial Mathematics** (`maths/finance/`): Options pricing, Black-Scholes, risk metrics
- **Actuarial Science** (`maths/actuarial/`): Life tables, annuities, mortality models
- **Econometrics** (`maths/econometrics/`): Time series analysis, regression models

## 🔬 Physics Modules

### Basic Physics

**Classical Mechanics:**
- **Newton's Laws** (`physics/newton_laws.hpp`): Force calculations, Newton's second law
- **Kinematics** (`physics/kinematics.hpp`): Position, velocity, acceleration equations
- **Dynamics** (`physics/dynamics.hpp`): Force systems, friction, tension
- **Energy & Momentum** (`physics/energy_momentum.hpp`): Conservation laws, work, kinetic/potential energy
- **Circular Motion** (`physics/circular_motion.hpp`): Centripetal force, angular velocity
- **Rotational Dynamics** (`physics/rotational_dynamics.hpp`): Torque, moment of inertia, angular momentum
- **Harmonic Motion** (`physics/harmonic_motion.hpp`): Simple harmonic oscillator, pendulum
- **Oscillations** (`physics/oscillations.hpp`): Damped and driven oscillations
- **Gravitation** (`physics/gravitation.hpp`): Universal gravitation, gravitational fields
- **Orbital** (`physics/orbital.hpp`): Orbital mechanics, Kepler's laws

**Electromagnetism:**
- **Electrostatics** (`physics/electrostatics.hpp`): Coulomb's law, electric fields, potential
- **Magnetism** (`physics/magnetism.hpp`): Magnetic fields, Lorentz force
- **Electric Circuits** (`physics/electric_circuits.hpp`): Ohm's law, RC/RL circuits
- **Electromagnetic Induction** (`physics/electromagnetic_induction.hpp`): Faraday's law, Lenz's law
- **Electromagnetic Waves** (`physics/electromagnetic_waves.hpp`): Wave propagation, Poynting vector
- **Maxwell Equations** (`physics/maxwell_equations.hpp`): Complete electromagnetic theory

**Waves and Optics:**
- **Wave Mechanics** (`physics/wave_mechanics.hpp`): Wave equation, interference, diffraction
- **Optics** (`physics/optics.hpp`): Reflection, refraction, lenses, mirrors

**Thermodynamics:**
- **Thermodynamics** (`physics/thermodynamics.hpp`): Laws of thermodynamics, entropy, cycles
- **Heat Transfer** (`physics/heat_transfer.hpp`): Conduction, convection, radiation
- **Thermal Expansion** (`physics/thermal_expansion.hpp`): Linear and volumetric expansion
- **Calorimetry** (`physics/calorimetry.hpp`): Specific heat, latent heat

**Fluid Mechanics:**
- **Fluid Mechanics** (`physics/fluid_mechanics.hpp`): Bernoulli's equation, continuity, viscosity
- **Surface Tension** (`physics/surface_tension.hpp`): Capillary action, contact angle

**Modern Physics:**
- **Special Relativity** (`physics/special_relativity.hpp`): Lorentz transformations, time dilation, E=mc²
- **Quantum Basics** (`physics/quantum_basics.hpp`): Planck's law, photoelectric effect, uncertainty principle

**Statics:**
- **Inclined Plane** (`physics/inclined_plane.hpp`): Forces on inclines, friction
- **Elasticity** (`physics/elasticity.hpp`): Hooke's law, Young's modulus, stress-strain

### Advanced Physics

**Advanced Classical Mechanics** (`physics/advanced/classical/`):
- **Hamiltonian Mechanics** (`hamiltonian.hpp`): Hamilton's equations, canonical transformations, generating functions
- **Phase Space** (`phase_space.hpp`): Phase space analysis, Poincaré sections, symplectic structure
- **Liouville Theorem** (`liouville.hpp`): Phase space volume conservation, statistical mechanics connection

**Cosmology** (`physics/advanced/cosmology/`):
- **Friedmann Equations** (`friedmann_equations.hpp`): FLRW metric, expansion dynamics, critical density
- **Expanding Universe** (`expanding_universe.hpp`): Hubble's law, scale factor evolution, redshift
- **Early Universe** (`early_universe.hpp`): Radiation/matter domination, recombination, nucleosynthesis
- **Energy Density** (`energy_density.hpp`): Matter, radiation, dark energy components

**Fluid Dynamics** (`physics/advanced/fluid_dynamics/`):
- **Governing Equations** (`governing_equations.hpp`): Navier-Stokes, continuity, energy equations
- **Flow Types** (`flow_types.hpp`): Laminar, turbulent, compressible, incompressible
- **Compressible Flow** (`compressible_flow.hpp`): Mach number, shock waves, supersonic flow
- **Boundary Layer** (`boundary_layer.hpp`): Boundary layer theory, separation, drag
- **Vorticity** (`vorticity.hpp`): Vorticity dynamics, circulation, Kelvin's theorem
- **Turbulence** (`turbulence.hpp`): Reynolds decomposition, energy cascade, turbulence models
- **Dimensionless Numbers** (`dimensionless_numbers.hpp`): Reynolds, Prandtl, Mach, Froude numbers

**Gauge Theory** (`physics/advanced/gauge_theory/`):
- **Gauge Invariance** (`gauge_invariance.hpp`): U(1), SU(2), SU(3) gauge symmetries
- **Higgs Mechanism** (`higgs_mechanism.hpp`): Spontaneous symmetry breaking, mass generation
- **Symmetries** (`symmetries.hpp`): Discrete and continuous symmetries, CPT theorem
- **Running Couplings** (`running_couplings.hpp`): Renormalization group, beta functions
- **Helicity** (`helicity.hpp`): Helicity conservation, polarization
- **CP Violation** (`cp_violation_kaons.hpp`): CP violation in kaon systems

**Quantum Field Theory** (`physics/advanced/qft/`):
- **Particle Physics** (`particle_physics.hpp`): Standard Model particles, interactions
- **Antiparticles** (`antiparticles.hpp`): Particle-antiparticle creation/annihilation
- **Interactions** (`interactions.hpp`): Electromagnetic, weak, strong interactions
- **Cross Sections** (`cross_sections.hpp`): Scattering amplitudes, differential cross sections
- **Decays** (`decays.hpp`): Decay rates, branching ratios, lifetime calculations
- **Spin Statistics** (`spin_statistics.hpp`): Fermi-Dirac, Bose-Einstein statistics
- **Supersymmetry** (`supersymmetry.hpp`): SUSY transformations, superpartners
- **Quark-Gluon Plasma** (`quark_gluon_plasma.hpp`): QCD matter at extreme temperatures

## 🎯 Demo Programs

All demos compile with: `g++ -std=c++17 -I./include -o demo examples/demo.cpp -lm`

### Differential Algebra Demos

1. **`differential_algebra_demo`**
   - Basic differential polynomials and fields
   - Characteristic sets and reduction
   - Differential ideals
   - Manifolds and dimension theory
   - Elimination theory
   - Low power theorem

2. **`differential_algebra_advanced_demo`**
   - Manifold intersections (Kronecker's theorem)
   - Orthonomic systems (Riquier)
   - Partial differential algebra
   - Classical PDEs: Laplace, wave, heat
   - Cauchy-Kovalevskaya existence

### Analysis Demos

3. **`fourier_analysis_demo`**
   - DFT vs FFT performance (11x speedup on N=256)
   - Circulant matrices and convolution theorem
   - Haar and Daubechies-4 wavelets
   - Fourier series and multipliers
   - STFT and spectrogram
   - Chirp signal detection

4. **`nonsmooth_algorithms_demo`**
   - Proximal operators (soft thresholding)
   - Subgradient descent
   - ISTA and FISTA comparison
   - ADMM consensus optimization
   - Sparse signal recovery (L1 minimization)

5. **`advanced_subdifferentials_demo`**
   - Clarke subdifferential computation
   - Normal and tangent cones
   - Limiting subdifferentials
   - Metric regularity testing
   - Nonsmooth optimization examples

### Stochastic Methods Demos

6. **`stochastic_methods_demo`**
   - Monte Carlo integration (basic, importance, stratified sampling)
   - Markov chains and stationary distributions
   - MCMC sampling (Metropolis-Hastings, Gibbs)
   - Brownian motion and stochastic processes
   - Boltzmann equation and DSMC simulation
   - Hamiltonian Monte Carlo
   - MCMC convergence diagnostics (ESS, R̂)

### Variational Calculus Demos

7. **`variational_calculus_demo`**
   - Contact geometry and contact structures
   - Lagrangians and Euler-Lagrange equations
   - Poincaré-Cartan forms and integral invariants
   - Legendre transformations (Lagrangian ↔ Hamiltonian)
   - Noether's theorem (symmetries → conservation laws)
   - Second variation and Jacobi fields
   - Bäcklund transformations (Sine-Gordon)
   - Field theories (Klein-Gordon, Maxwell, Yang-Mills)

## 🚀 Building and Running

### Prerequisites
- C++17 compatible compiler (g++, clang++)
- Standard library with `<random>`, `<vector>`, `<cmath>`

### Compile Individual Demo
```bash
g++ -std=c++17 -I./include -o fourier_demo examples/fourier_analysis_demo.cpp -lm
./fourier_demo
```

### Compile All Demos
```bash
for demo in examples/*_demo.cpp; do
    name=$(basename "$demo" .cpp)
    g++ -std=c++17 -I./include -o "$name" "$demo" -lm
    echo "Compiled $name"
done
```

## ✨ Features

### Design Philosophy
- **Computational**: All functions return numerical results, not strings
- **Header-Only**: Easy integration, no linking required
- **Zero Dependencies**: Only standard library
- **Well-Documented**: Detailed @param, @return, mathematical formulas
- **Tested**: Each module has comprehensive demo showing practical usage

### Performance
- FFT: 11x speedup over DFT (N=256)
- FISTA: O(1/k²) convergence vs ISTA's O(1/k)
- Reduction algorithms with termination safeguards
- Efficient circulant matrix multiplication via FFT

### Mathematical Rigor
- Based on standard texts:
  - Ritt's "Differential Algebra" (1950)
  - Clarke's "Optimization and Nonsmooth Analysis"
  - Mallat's "A Wavelet Tour of Signal Processing"
  - Arnold's "Mathematical Methods of Classical Mechanics"
  - Standard texts on stochastic processes and Monte Carlo methods

## 📊 Statistics

- **Mathematics Modules**:
  - Basic: 4 modules (calculus, linear algebra, trigonometry, transforms)
  - Advanced: 11 modules (differential algebra, Fourier analysis, subdifferentials, nonsmooth algorithms, stochastic methods, variational calculus, dynamical systems, probability, PDEs - classification, solutions, transforms)
  - Applied: 3 modules (finance, actuarial, econometrics)
- **Physics Modules**:
  - Basic: 25+ modules covering classical mechanics, E&M, thermodynamics, optics, modern physics
  - Advanced: 20+ modules in Hamiltonian mechanics, cosmology, fluid dynamics, gauge theory, QFT
- **Demos**: 7+ comprehensive demonstration programs
- **Distributions**: 13 probability distributions (Bernoulli, Binomial, Poisson, Geometric, Negative Binomial, Hypergeometric, Uniform, Normal, Exponential, Gamma, Beta, Chi-squared, F-distribution)
- **Key Algorithms**:
  - DFT, FFT (O(N log N))
  - ISTA, FISTA (O(1/k²) convergence)
  - ADMM, Ritt's algorithm
  - Monte Carlo, MCMC (Metropolis-Hastings, Gibbs, HMC)
  - Euler-Lagrange, Noether's theorem, Legendre transforms
  - RK4, Picard iteration, Floquet theory
  - Lyapunov exponents, bifurcation diagrams, fractal dimensions
  - Method of characteristics (linear, quasi-linear, nonlinear PDEs)
  - Charpit's method for fully nonlinear equations
  - PDE classification via discriminant
  - Separation of variables (wave, heat, Laplace equations)
  - Orthogonal polynomial expansions (Legendre, Chebyshev, Hermite, Laguerre)
  - Bessel function computations and zero-finding
  - Fourier series coefficient computation
  - Laplace transforms and inverse transforms
  - Fourier transforms (full, sine, cosine, finite)
  - d'Alembert's solution for wave equation
  - Heat kernel and fundamental solutions
  - Green's functions for Poisson equation
  - Poisson integral formula

## 🎓 Educational Value

Each module serves as both:
1. **Production-ready code** for numerical computations
2. **Educational reference** showing how abstract mathematics translates to algorithms
3. **Research tool** for:
   - Ordinary and partial differential equations
   - Dynamical systems, chaos theory, and bifurcation analysis
   - Nonsmooth optimization and variational calculus
   - Signal processing and Fourier analysis
   - Stochastic methods and Monte Carlo simulations
   - Classical and quantum field theory
   - Statistical mechanics and computational physics
   - PDE theory: parabolic (heat), elliptic (Laplace/Poisson), hyperbolic (wave)
   - Transform methods: Laplace and Fourier transforms for PDEs
   - Boundary value problems and initial value problems
   - Green's functions and fundamental solutions

## 📝 License

*(Add your license here)*

## 🤝 Contributing

Contributions welcome! Please ensure:
- Follow the computational pattern (no string returns)
- Include comprehensive demos
- Document with @param, @return
- Test all functionality

## 📧 Contact

*(Add contact information)*

---

**Built with C++17 | Header-Only | Zero Dependencies**
