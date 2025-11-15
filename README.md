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
│   ├── maths/                     # All mathematics modules (flattened)
│   │   ├── actuarial_life_tables.hpp
│   │   ├── advanced_subdifferentials.hpp
│   │   ├── black_scholes.hpp
│   │   ├── calculus_theorems.hpp
│   │   ├── complex_analysis.hpp   # NEW: Zeros, products, Gamma, Blaschke
│   │   ├── differential_algebra.hpp
│   │   ├── distributions.hpp
│   │   ├── econometrics_regression.hpp
│   │   ├── fourier_analysis.hpp
│   │   ├── matrices.hpp
│   │   ├── monte_carlo.hpp
│   │   ├── nonsmooth_algorithms.hpp
│   │   ├── ode_dynamical_systems.hpp
│   │   ├── partial_differential_equations.hpp
│   │   ├── pde_classification_solutions.hpp
│   │   ├── pde_numerical_methods.hpp
│   │   ├── pde_solution_methods.hpp
│   │   ├── pde_transform_methods.hpp
│   │   ├── pde_variational_methods.hpp
│   │   ├── polar_transforms.hpp
│   │   ├── stochastic_differential_equations.hpp
│   │   ├── trigonometry_identities.hpp
│   │   ├── variational_calculus.hpp
│   │   └── vectors.hpp
│   └── physics/
│       ├── (basic modules)        # Classical mechanics, waves, etc.
│       ├── advanced_quantum_mechanics.hpp      # NEW: Advanced QM topics
│       ├── quantum_chemistry.hpp               # NEW: Atomic/molecular structure
│       ├── quantum_foundations.hpp             # NEW: Historical QM development
│       ├── relativistic_quantum_mechanics.hpp  # NEW: Spin and Dirac theory
│       └── advanced/              # Advanced physics topics
│           ├── classical/         # Hamiltonian, Liouville, phase space
│           ├── cosmology/         # Friedmann equations, early universe
│           ├── fluid_dynamics/    # Turbulence, compressible flow
│           ├── gauge_theory/      # Gauge invariance, Higgs mechanism
│           ├── operator_algebras.hpp  # NEW: Von Neumann, C*-algebras
│           └── qft/              # Quantum field theory
├── examples/                      # Physics demonstration programs
└── README.md
```

## 🧮 Mathematics Modules

All mathematics modules are located in `include/maths/` with a flat structure for easy access.

### Differential Algebra (`maths/differential_algebra.hpp`)
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

### Fourier Analysis (`maths/fourier_analysis.hpp`)
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

### Advanced Subdifferentials (`maths/advanced_subdifferentials.hpp`)
**Nonsmooth analysis & variational calculus**

- **Clarke Subdifferential** ∂_C f(x)
- **Limiting (Mordukhovich) Subdifferential** ∂_L f(x)
- **Fréchet, Proximal, Graded Subdifferentials**
- **Normal & Tangent Cones**: N_C(x), T_C(x)
- **Coderivatives** D*F for set-valued mappings
- **Calculus Rules**: sum, chain, max rules
- **Metric Regularity Criterion**

### Nonsmooth Algorithms (`maths/nonsmooth_algorithms.hpp`)
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

### Stochastic Methods (`maths/monte_carlo.hpp`)
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

### Stochastic Differential Equations (`maths/stochastic_differential_equations.hpp`)
**Itô calculus and stochastic processes**

- **Itô Integrals**
  - Construction of Itô integrals: ∫ f(t, W) dW
  - Itô isometry: E[|∫ f dW|²] = E[∫ f² dt]
  - Quadratic variation and properties
  - Extensions for adapted processes

- **Itô's Lemma**
  - Change of variables formula for SDEs
  - Application to geometric Brownian motion
  - Multidimensional Itô formula
  - Derivation of Black-Scholes PDE

- **Stochastic Differential Equations**
  - Euler-Maruyama method (strong order 0.5)
  - Milstein method (strong order 1.0)
  - Ornstein-Uhlenbeck process (mean-reverting)
  - Stochastic analogs of classical ODEs

- **Filtering Problems**
  - Kalman filter (prediction and update steps)
  - Optimal state estimation
  - Linear Gaussian systems
  - Innovation and Kalman gain

- **Optimal Stopping**
  - American option pricing
  - Dynamic programming approach
  - Optimal stopping times
  - Reward maximization problems

- **Stochastic Control**
  - Hamilton-Jacobi-Bellman equation
  - Linear-Quadratic-Gaussian (LQG) control
  - Merton's portfolio optimization
  - Optimal investment strategies

- **Mathematical Finance Applications**
  - Heston stochastic volatility model
  - Cox-Ingersoll-Ross (CIR) interest rate model
  - Vasicek model for term structure
  - Multi-factor models

**Applications:** Quantitative finance, optimal control, signal processing, filtering theory, economics

### Variational Calculus (`maths/variational_calculus.hpp`)
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

### Differential Equations and Dynamical Systems (`maths/ode_dynamical_systems.hpp`)
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

### Partial Differential Equations (`maths/partial_differential_equations.hpp`)
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

### PDE Solution Methods (`maths/pde_solution_methods.hpp`)
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
- **Series of Orthogonal Functions**:
  - Parseval's identity: ∥f∥² = ∑|cₙ|² (energy conservation)
  - Bessel's inequality: ∑|cₙ|² ≤ ∥f∥² (completeness criterion)
  - Convergence rate analysis: cₙ ~ 1/nᵖ decay estimation
  - Mean square error of partial sum approximations
- **Eigenfunction Expansions**:
  - Sturm-Liouville theory: (p(x)u')' + (q(x) + λw(x))u = 0
  - Eigenvalue computation via shooting method
  - Eigenfunction orthogonality with weight w(x)
  - Function expansion: f(x) = ∑ cₙφₙ(x) with cₙ = ⟨f,φₙ⟩/⟨φₙ,φₙ⟩
  - Standard problems: Fourier sine (λₙ = n²), Bessel (cylindrical), Legendre (spherical)

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

### PDE Transform Methods (`maths/pde_transform_methods.hpp`)
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

### PDE Classification Solutions (`maths/pde_classification_solutions.hpp`)
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

**Green's Functions for All PDE Types:**
- **Parabolic (Heat)**:
  - 1D/2D heat kernels: G(x,t;ξ,τ) = 1/√(4πα(t-τ)) exp(-(x-ξ)²/(4α(t-τ)))
  - Half-space with Dirichlet/Neumann BC via method of images
  - Solution with source terms via convolution
- **Elliptic (Poisson)**:
  - 2D: G(x,y;ξ,η) = -(1/2π) ln(r), 3D: G = -1/(4πr)
  - Half-space with Dirichlet BC using method of images
  - Rectangle domain with eigenfunction expansion
  - Poisson solver: u(x,y) = ∫∫ G(x,y;ξ,η) f(ξ,η) dξdη
- **Hyperbolic (Wave)**:
  - 1D retarded: G(x,t;ξ,τ) = 1/(2c) H(t-τ-|x-ξ|/c)
  - 2D (odd dimensions): G = H(c(t-τ) - r) / (2π√(c²(t-τ)² - r²))
  - 3D with Huygens' principle: G = δ(c(t-τ) - r) / (4πr)
  - Duhamel's principle for source terms
  - Causality and light cone checking
- **Method of Images**: Image source locations for Dirichlet/Neumann BC

### PDE Variational Methods (`maths/pde_variational_methods.hpp`)
**Weak formulations and variational methods for PDEs**

**Line Integrals and Variational Notation:**
- **Line Integrals**: ∫_C F·dr along curves for variational formulations
- **Variational Derivatives**: δF/δu via Gateaux derivatives
- **Functional Derivatives**: Euler-Lagrange equations for functionals

**Multiple Integrals:**
- **Double and Triple Integrals**: Change of variables, Jacobians
- **Divergence Theorem**: ∫_Ω div(F) dV = ∫_∂Ω F·n dS
- **Green's Identities**: First, second, and third identities
  - First identity: ∫_Ω (∇u·∇v + v∇²u) dV = ∫_∂Ω v(∂u/∂n) dS
  - Second identity: ∫_Ω (v∇²u - u∇²v) dV = ∫_∂Ω (v∂u/∂n - u∂v/∂n) dS
  - Integration by parts in multiple dimensions

**Weak Variational Formulation:**
- **Test Functions**: Compact support, smoothness requirements
- **Trial Functions**: Finite-dimensional approximations
- **Weak Derivatives**: Distributional derivatives
- **Sobolev Spaces**: H¹, H₀¹ function spaces
- **Weak Solutions**: ∫_Ω ∇u·∇v dx = ∫_Ω fv dx for all test functions v

**Weighted Residual Methods (WRM):**
- **General Framework**: ∫_Ω w(x) R(x) dx = 0 where R = L[u] - f is residual
- **Test Function Selection**:
  - Hat functions (piecewise linear finite elements)
  - Polynomial basis with boundary conditions
  - Trigonometric functions (Fourier modes)
  - Bubble functions for incompressible flow
  - Completeness and linear independence criteria
- **Collocation Method**: w(x) = δ(x - xᵢ) at collocation points
  - Chebyshev nodes for spectral accuracy
  - Direct enforcement R(xᵢ) = 0
- **Subdomain Method**: w(x) = 1 on subdomain Ωᵢ, 0 elsewhere
  - Domain partitioning strategies
  - Integrated residual minimization
- **Least Squares Method**: Minimize J = ∫_Ω R² dx
  - Optimality conditions: ∂J/∂cᵢ = 0
  - Symmetric positive definite systems

**Galerkin Method:**
- **Finite Element Approximation**: Basis function expansion
- **Hat Functions**: Piecewise linear basis
- **Stiffness Matrix**: a(φᵢ, φⱼ) assembly
- **Load Vector**: L(φᵢ) computation
- **Galerkin Orthogonality**: Optimal approximation in energy norm

**Rayleigh-Ritz Method:**
- **Energy Minimization**: E[u] = ½∫(u')² dx - ∫fu dx
- **Rayleigh Quotient**: R[u] = ∫(u')² dx / ∫u² dx for eigenvalues
- **Ritz Coefficients**: Minimize energy functional
- **Upper Bounds**: Eigenvalue estimates

**Transient Problems:**
- **Semi-Discrete Methods**: Spatial discretization first
- **Time Stepping**: Backward Euler, Crank-Nicolson
- **Mass and Stiffness Matrices**: (M + Δt·K)u^(n+1) = Mu^n + Δt·F
- **Energy Stability**: ½d/dt(∫u² dx) ≤ 0

**Applications**: Finite element analysis, structural mechanics, computational fluid dynamics, elasticity

### PDE Numerical Methods (`maths/pde_numerical_methods.hpp`)
**Numerical approximation and finite difference schemes**

**Taylor Series Expansions:**
- **Forward Difference**: f'(x) ≈ (f(x+h) - f(x))/h, O(h) error
- **Backward Difference**: f'(x) ≈ (f(x) - f(x-h))/h, O(h) error
- **Central Difference**: f'(x) ≈ (f(x+h) - f(x-h))/(2h), O(h²) error
- **Second Derivative**: f''(x) ≈ (f(x+h) - 2f(x) + f(x-h))/h²
- **Higher Order Schemes**: 4th order accurate central differences
- **Truncation Error Analysis**: Leading error terms

**Successive Approximations:**
- **Picard Iteration**: u_{n+1}(t) = u₀ + ∫ f(s, u_n(s)) ds
- **Fixed Point Methods**: u_{n+1} = G(u_n), convergence criteria
- **Successive Over-Relaxation (SOR)**: ω-parameter for acceleration
- **Convergence**: Banach fixed point theorem

**Boundary Perturbations:**
- **Regular Perturbation**: u = u₀ + εu₁ + ε²u₂ + ...
- **Singular Perturbation**: Boundary layers, matched asymptotics
- **Boundary Layer Thickness**: δ ~ √ε for second order problems
- **Inner and Outer Expansions**: Composite solutions

**Perturbation Methods:**
- **Multiple Scales Analysis**: Disparate time/space scales T₀ = t, T₁ = εt
  - Solvability conditions to eliminate secular terms
  - Fast and slow scale separation
- **Matched Asymptotic Expansions**: Singular perturbation problems
  - Outer expansion (valid away from boundary)
  - Inner expansion (boundary layer with stretched coordinate ξ = x/δ(ε))
  - Van Dyke matching principle
  - Composite solutions
- **WKB Approximation**: Rapidly oscillating solutions
  - Eikonal equation: (dS/dx)² = k²(x)
  - Amplitude expansion: A = A₀ + εA₁ + ...
  - Connection formulas at turning points
- **Poincare-Lindstedt Method**: Nonlinear oscillators
  - Frequency correction: ω = ω₀ + εω₁ + ε²ω₂ + ...
  - Elimination of secular terms
  - Stretched time coordinate
- **Asymptotic Sequence Verification**: Check uₙ₊₁/uₙ → 0 as ε → 0

**Finite Difference Schemes for First Order Equations:**
- **Upwind Scheme**: Backward difference for c > 0, first order accurate
- **Lax-Friedrichs**: Central difference with averaging, stable
- **Lax-Wendroff**: Second order in space and time
- **CFL Condition**: |c|Δt/Δx ≤ 1 for stability
- **Stability Analysis**: Von Neumann stability analysis

**Finite Difference Schemes for Second Order Equations:**
- **Explicit Heat Equation**: u_i^{n+1} = u_i^n + r(u_{i+1}^n - 2u_i^n + u_{i-1}^n)
  - Stable if r = αΔt/Δx² ≤ 1/2
- **Implicit (Backward Euler)**: Unconditionally stable, first order in time
- **Crank-Nicolson**: θ = 1/2, unconditionally stable, second order in time
- **ADI (Alternating Direction Implicit)**: Efficient 2D solver
  - Step 1: Implicit in x, explicit in y
  - Step 2: Explicit in x, implicit in y
- **Stability**: Amplification factor analysis, von Neumann method
- **Tridiagonal Systems**: Thomas algorithm O(n) solution

**Key Algorithms**: Upwind, Lax-Friedrichs, Lax-Wendroff, Crank-Nicolson, ADI, SOR, Picard iteration, multiple scales, matched asymptotics, WKB

**Applications**: Computational fluid dynamics, heat transfer, wave propagation, image processing, option pricing, boundary layer problems, quantum mechanics

### Probability & Statistics (`maths/distributions.hpp`)
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
- **Student's t**: t-distribution with ν degrees of freedom, small sample inference
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

### Complex Analysis (`maths/complex_analysis.hpp`)
**Advanced complex function theory**

- **Zeros of Holomorphic Functions**
  - Argument principle: (1/2πi) ∮ f'/f dz counts zeros minus poles
  - Rouché's theorem for zero counting
  - Zero multiplicity and order computation
  - Jensen's formula relating zeros to growth

- **Infinite Products**
  - Weierstrass factorization theorem
  - Elementary factors E_n(z)
  - Canonical products and genus
  - Hadamard's theorem for entire functions

- **Ring H(D)**
  - Principal ideals in holomorphic functions
  - Common zeros and greatest common divisors
  - Maximal ideals and point evaluation
  - Identity theorem applications

- **Euler's Gamma Function**
  - Γ(z) via Weierstrass product formula
  - Reflection formula: Γ(z)Γ(1-z) = π/sin(πz)
  - Duplication and multiplication formulas
  - Stirling's approximation for large |z|
  - Beta function: B(z,w) = Γ(z)Γ(w)/Γ(z+w)
  - Pochhammer symbols (rising factorials)
  - Digamma function ψ(z) = Γ'(z)/Γ(z)

- **Divisors and Meromorphic Functions**
  - Divisor representation for meromorphic functions
  - Principal divisors and equivalence
  - Construction of meromorphic functions from divisors
  - Mittag-Leffler theorem for prescribed poles

- **Infinite Blaschke Products**
  - Blaschke condition: Σ(1 - |aₙ|) < ∞
  - Blaschke factor: b_a(z) = (a-z)/(1-āz)
  - Convergence and boundary behavior
  - Applications to Hardy spaces

- **Confluent Hypergeometric Functions**
  - Kummer's function M(a,b,z) = ₁F₁(a;b;z)
  - Kummer's U function (second solution)
  - Associated Laguerre polynomials L_n^k(x)
  - Applications to hydrogen atom radial functions

**Applications:** Analytic number theory, complex dynamics, quantum mechanics, special functions

### Basic Mathematics

- **Calculus** (`maths/calculus_theorems.hpp`): Numerical derivatives, integration (Simpson's rule)
- **Trigonometry** (`maths/trigonometry_identities.hpp`): Computational trig identities
- **Linear Algebra** (`maths/matrices.hpp`, `maths/vectors.hpp`): Matrix operations, vectors, eigenvalues
- **Transforms** (`maths/polar_transforms.hpp`): Polar coordinates and transformations

### Applied Mathematics

- **Financial Mathematics** (`maths/black_scholes.hpp`): Options pricing, Black-Scholes, risk metrics
- **Actuarial Science** (`maths/actuarial_life_tables.hpp`): Life tables, annuities, mortality models
- **Econometrics** (`maths/econometrics_regression.hpp`): Time series analysis, regression models

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

**Operator Algebras and Quantum Mechanics** (`physics/advanced/operator_algebras.hpp`):
Comprehensive functional analysis framework for quantum mechanics (~2,800 lines)

- **Hilbert Spaces**
  - Inner products ⟨ψ|φ⟩ and norms
  - Gram-Schmidt orthogonalization
  - Orthonormal bases and completeness
  - Tensor products |ψ⟩ ⊗ |φ⟩ for composite systems

- **Bounded Operators**
  - Operator norm ||A|| = sup ||Aψ||
  - Adjoints A†, self-adjoint operators
  - Unitary operators: U†U = I
  - Commutators [A,B] = AB - BA
  - Trace and partial trace operations
  - Projection operators

- **Von Neumann Algebras (Rings of Operators)**
  - Commutants A' = {B : AB = BA for all A}
  - Double commutant theorem: A'' = A (weak closure)
  - Factors (von Neumann algebras with trivial center)
  - Quantum observables as self-adjoint operators

- **Unitary Representations**
  - Group representations on Hilbert spaces
  - Irreducibility and Schur's lemma
  - Characters and orthogonality
  - Direct sums and tensor products
  - Reduction of representations

- **Murray-von Neumann Factor Classification**
  - Type I: B(H), finite/infinite dimensional
  - Type II₁: Finite factors (tracial states)
  - Type II∞: Infinite factors
  - Type III: Properly infinite factors
  - Dimension theory and comparison of projections

- **Elementary C*-Algebra Theory**
  - **Banach Algebras**: Submultiplicativity, spectral radius, Neumann series, resolvent
  - **Commutative Banach Algebras**: Characters, Gelfand transform, maximal ideals, Shilov boundary
  - **Commutative C*-Algebras**: Gelfand-Naimark theorem (isomorphism with C(X))
  - **Spectrum and Functional Calculus**: σ(A), continuous functional calculus, spectral mapping theorem
  - **Positivity**: Positive cone, order structure, square roots, polar decomposition
  - **Ideals**: Left/right/two-sided ideals, quotients, maximal ideals, simplicity
  - **States**: State functionals, vector states, pure states, faithful states, tracial states
  - **Representations and GNS Construction**: GNS theorem (states → representations), cyclic vectors, universal representation
  - **Gelfand-Naimark Theorem**: Every C*-algebra embeds in B(H), abstract vs concrete C*-algebras
  - **Complete Positivity**: CP maps, Stinespring dilation, Kraus representation, quantum channels
  - **Pure States and Irreducible Representations**: Correspondence via GNS, extremal points
  - **Compact Operators**: K(H), finite-rank operators, Calkin algebra B(H)/K(H), Fredholm operators
  - **Double Commutant Theorem**: von Neumann bicommutant A'' = Ā^SOT, Kaplansky density

**Applications:** Quantum mechanics foundations, quantum information theory, quantum field theory, statistical mechanics, mathematical physics

**Quantum Mechanics Foundations** (`physics/quantum_foundations.hpp`):
Historical development of quantum mechanics (~1,000 lines)

- **Introduction to Quantum Mechanics**
  - Failures of classical physics
  - Ultraviolet catastrophe in blackbody radiation
  - Photoelectric effect paradox
  - Atomic stability problem
  - Need for quantization

- **Planck and Quantization**
  - Planck's blackbody radiation law: u(ν,T) = (8πhν³/c³)/(e^(hν/kT) - 1)
  - Energy quantization: E = nhν
  - Planck constant h = 6.626 × 10⁻³⁴ J·s
  - Photoelectric effect: E_kinetic = hν - W
  - Einstein's photon hypothesis
  - Specific heat models (Einstein, Debye)

- **Bohr and the Hydrogen Atom**
  - Bohr model of hydrogen
  - Orbital radii: r_n = n²a₀ (Bohr radius a₀ = 0.529 Å)
  - Energy levels: E_n = -13.6 eV/n²
  - Rydberg formula: 1/λ = R_∞(1/n₁² - 1/n₂²)
  - Spectral series: Lyman, Balmer, Paschen, Brackett, Pfund
  - Angular momentum quantization: L = nℏ

- **Matrix Mechanics (Heisenberg)**
  - Heisenberg's formulation with matrices
  - Position and momentum matrices
  - Canonical commutation relation: [x,p] = iℏ
  - Ladder operators a, a† for harmonic oscillator
  - Energy eigenvalues: E_n = ℏω(n + 1/2)
  - Matrix elements and transition amplitudes

- **Uncertainty Relations**
  - Heisenberg uncertainty principle: ΔxΔp ≥ ℏ/2
  - General uncertainty: ΔAΔB ≥ ½|⟨[A,B]⟩|
  - Wave packet spreading: σ_x(t) = σ₀√(1 + (ℏt/2mσ₀²)²)
  - Energy-time uncertainty: ΔEΔt ≥ ℏ/2
  - Coherent states (minimum uncertainty)

- **Wave Mechanics (Schrödinger)**
  - Schrödinger's wave formulation
  - De Broglie relations: λ = h/p, ω = E/ℏ
  - Wave function ψ(x,t) and probability interpretation
  - Time-dependent Schrödinger equation: iℏ∂ψ/∂t = Hψ
  - Time-independent equation: Hψ = Eψ
  - Gaussian wave packets
  - Harmonic oscillator eigenstates with Hermite polynomials
  - Born rule: P(x) = |ψ(x)|²

**Applications:** Quantum mechanics education, atomic physics, quantum chemistry, historical physics

**Advanced Quantum Mechanics** (`physics/advanced_quantum_mechanics.hpp`):
Advanced topics in quantum mechanics (~1,650 lines)

- **Kummer's Confluent Hypergeometric Functions**
  - Kummer's M function: M(a,b,z) = ₁F₁(a;b;z)
  - Kummer's U function (second solution)
  - Pochhammer symbols and series expansions
  - Associated Laguerre polynomials: L_n^k(x)
  - Hydrogen radial wave functions via Laguerre polynomials

- **Hamiltonian Mechanics**
  - Hamiltonian H(q,p) = p²/2m + V(q)
  - Hamilton's equations: dq/dt = ∂H/∂p, dp/dt = -∂H/∂q
  - Poisson brackets: {f,g} = ∂f/∂q·∂g/∂p - ∂f/∂p·∂g/∂q
  - Canonical transformations and generating functions
  - Phase space trajectories
  - Liouville's theorem connection

- **Classical Harmonic Oscillator**
  - Classical solutions: x(t) = A cos(ωt + φ)
  - Phase space ellipses
  - Energy E = ½mω²A²
  - Action-angle variables
  - Correspondence with quantum oscillator

- **Mathematics of Plane Waves**
  - Plane wave solutions: ψ(x,t) = Ae^(i(kx-ωt))
  - Dispersion relations: ω(k) = ℏk²/2m
  - Fourier transforms and wave packets
  - Parseval's theorem and normalization
  - Group and phase velocities

- **Schrödinger Equation for Free Particle**
  - Free particle solutions: ψ_k(x,t) = e^(i(kx-ωt))
  - Energy-momentum relation: E = ℏ²k²/2m
  - Continuity equation: ∂ρ/∂t + ∇·j = 0
  - Probability current density
  - Time evolution operators

- **Wave Functions and Wave Packets**
  - Gaussian wave packets: ψ(x,0) = (2πσ²)^(-1/4) e^(-x²/4σ²) e^(ik₀x)
  - Wave packet spreading with time
  - Normalization and expectation values
  - Position and momentum uncertainties
  - Fourier transforms between representations

- **Quantum Tunneling**
  - Transmission coefficients for barriers
  - Rectangular barrier: T = 1/(1 + V₀²sinh²(κa)/4E(V₀-E))
  - WKB approximation: T ≈ exp(-2∫κ(x)dx)
  - Alpha decay and nuclear physics
  - Scanning tunneling microscopy (STM)
  - Tunneling time and probability

- **Perturbation Theory (Nondegenerate States)**
  - First-order energy correction: E_n^(1) = ⟨ψ_n^(0)|H'|ψ_n^(0)⟩
  - Second-order energy: E_n^(2) = Σ|⟨ψ_k|H'|ψ_n⟩|²/(E_n-E_k)
  - First-order wave function correction
  - Anharmonic oscillator perturbations
  - Convergence and validity conditions

- **Stark Effect of the Hydrogen Atom**
  - Linear Stark effect (degenerate states): ΔE ∝ E
  - Quadratic Stark effect (ground state): ΔE ∝ E²
  - Hydrogen polarizability α = (9/2)a₀³
  - Stark splitting for n=2: ΔE = 3eEa₀
  - Matrix elements and selection rules
  - Avoided crossings in Stark maps

- **Pauli's Exclusion Principle**
  - No two fermions in same quantum state
  - Quantum number uniqueness (n,l,m_l,m_s)
  - Shell filling: maximum 2n² electrons per shell
  - Subshell capacity: 2(2l+1) electrons
  - Antisymmetric wave functions for fermions
  - Slater determinants
  - Fermi energy: E_F = (ℏ²/2m)(3π²n)^(2/3)
  - Degeneracy pressure in white dwarfs

- **Electron Spin**
  - Spin quantum number s = 1/2
  - Spin angular momentum: |S| = (√3/2)ℏ
  - Spin z-component: S_z = ±ℏ/2
  - Spinors: |↑⟩ = (1,0)ᵀ, |↓⟩ = (0,1)ᵀ
  - Pauli matrices σ_x, σ_y, σ_z
  - Magnetic moment: μ = -g_e(e/2m)S
  - Zeeman effect: ΔE = g_e μ_B m_s B
  - Stern-Gerlach experiment
  - Spin-orbit coupling: H_SO ∝ L·S
  - Fine structure in hydrogen

- **Two-Electron Systems**
  - Total wave function antisymmetry
  - Singlet state (S=0): |S=0,M=0⟩ = (1/√2)(|↑↓⟩ - |↓↑⟩)
  - Triplet states (S=1): |S=1,M⟩ with M = -1,0,+1
  - Spatial symmetry requirements
  - Exchange energy: ΔE = 2K
  - Direct Coulomb integral J
  - Exchange integral K
  - Ortho and para states

- **Helium Atom**
  - Ground state energy: E₀ = -79.0 eV (experimental)
  - Independent particle approximation: -108.8 eV
  - Variational method with Z_eff
  - Optimal screening: Z_eff = Z - 5/16 ≈ 1.69
  - First ionization energy: 24.6 eV
  - Second ionization energy: 54.4 eV
  - Electron-electron repulsion corrections
  - Excited states: 1s2s configuration
  - Singlet ¹S and triplet ³S energies
  - Exchange splitting

- **Helium Atom Orbitals**
  - Hydrogenic orbital approximations
  - Radial wave functions: R_1s(r), R_2s(r), R_2p(r)
  - Effective nuclear charge Z_eff
  - Product wave functions: ψ(r₁,r₂) = ψ_a(r₁)ψ_b(r₂)
  - Symmetric spatial: ψ_+ = (1/√2)[ψ_a(1)ψ_b(2) + ψ_a(2)ψ_b(1)]
  - Antisymmetric spatial: ψ_- = (1/√2)[ψ_a(1)ψ_b(2) - ψ_a(2)ψ_b(1)]
  - Probability densities |ψ|²
  - Radial expectation values ⟨r⟩
  - Most probable radii

**Applications:** Atomic physics, quantum chemistry, spectroscopy, multi-electron systems, perturbation theory, solid-state physics

**Quantum Chemistry: Atomic and Molecular Structure** (`physics/quantum_chemistry.hpp`):
Comprehensive quantum chemistry module for atoms and molecules (~1,300 lines)

**Atomic Structure:**

- **Atomic and Molecular Wave Functions**
  - Multi-electron wave functions: ψ(1,2,...,N)
  - Product wave functions vs antisymmetrized
  - Slater determinants for fermions: ψ(1,2) = -ψ(2,1)
  - Normalization integrals ∫|ψ|² dτ = 1
  - Spin-spatial factorization: ψ(r,s) = ψ_spatial(r) × χ_spin(s)
  - Exchange symmetry verification

- **The Hartree-Fock Method**
  - Self-consistent field (SCF) theory
  - Fock operator: F = h + Σⱼ(2Jⱼ - Kⱼ)
  - Coulomb integral Jᵢⱼ: electron-electron repulsion
  - Exchange integral Kᵢⱼ: quantum exchange effects
  - Hartree-Fock energy: E_HF = Σᵢhᵢᵢ + ½ΣᵢΣⱼ(2Jᵢⱼ - Kᵢⱼ)
  - SCF iteration and convergence criteria
  - Koopmans' theorem: ionization energy ≈ -εᵢ

- **Slater Orbitals**
  - Slater-type orbitals (STOs): φₙₗₘ = N rⁿ⁻¹ e^(-ζr) Yₗₘ
  - Slater's rules for screening constants
  - Effective nuclear charge: Z_eff = Z - S
  - Slater exponents ζ = Z_eff/n*
  - Overlap integrals between STOs
  - Orbital normalization

- **Multiplet Theory**
  - Term symbols: ²ˢ⁺¹Lⱼ notation
  - L-S coupling (Russell-Saunders): L = Σlᵢ, S = Σsᵢ
  - Total angular momentum: J = L + S
  - Hund's rules for ground states:
    1. Maximize total spin S
    2. Maximize total orbital angular momentum L
    3. J = |L-S| if less than half-filled, J = L+S if more
  - Spectroscopic notation (S, P, D, F, G, ...)
  - Fine structure splitting
  - Multiplicity 2S+1

**Molecular Structure:**

- **The Born-Oppenheimer Approximation**
  - Electronic-nuclear motion separation
  - Mass ratio justification: m_e/M_n << 1
  - Wave function factorization: Ψ(r,R) ≈ ψ_el(r;R) × χ_nuc(R)
  - Electronic Hamiltonian at fixed nuclear positions
  - Adiabatic vs diabatic representations
  - Validity criterion: ω_vib << ω_el

- **Nuclear Motion of Diatomic Molecules**
  - Reduced mass: μ = m₁m₂/(m₁ + m₂)
  - Rotational energy levels: E_J = BJ(J+1)
  - Rotational constant: B = ℏ²/(2I)
  - Vibrational energy (harmonic): E_v = ℏω(v + 1/2)
  - Anharmonic corrections: -χₑℏω(v + 1/2)²
  - Rovibrational coupling: E(v,J) = E_vib + E_rot
  - Centrifugal distortion: -DJ²(J+1)²
  - Morse potential: V(R) = Dₑ[1 - e^(-a(R-Rₑ))]²
  - Selection rules: ΔJ = ±1, Δv = ±1

- **The Hydrogen Molecular Ion H₂⁺**
  - LCAO (Linear Combination of Atomic Orbitals)
  - Molecular orbitals: ψ = c₁φ_A ± c₂φ_B
  - Bonding (σ_g) and antibonding (σ_u*) orbitals
  - Bonding/antibonding energies: E_± = (H_AA ± H_AB)/(1 ± S_AB)
  - Overlap integral S_AB for 1s orbitals
  - Equilibrium bond length: R_e ≈ 2.5a₀
  - Dissociation energy: D₀ ≈ 2.8 eV
  - Energy curve E(R)

- **The Hydrogen Molecule H₂**
  - Molecular orbital configuration: (σ_g 1s)²
  - Valence bond (VB) wave function: covalent structure
  - Molecular orbital (MO) wave function
  - Heitler-London approximation: E = (Q + J)/(1 + S²)
  - Bond dissociation energy: D₀ = 4.75 eV
  - Equilibrium bond length: R_e = 0.74 Å
  - Ionic-covalent resonance: ψ = c₁ψ_covalent + c₂ψ_ionic
  - Comparison of VB and MO theories

- **The Chemical Bond**
  - Bond order: BO = (n_bonding - n_antibonding)/2
  - σ, π, and δ bonds
  - Hybridization: sp, sp², sp³, sp³d, sp³d²
  - Electronegativity and ionic character
  - Percent ionic character: 100[1 - e^(-0.25Δχ²)]
  - Bond length correlation with bond order
  - Bond energy correlation with bond order
  - Resonance structures and hybrid energies

- **Structures of Simple Polyatomic Molecules**
  - VSEPR (Valence Shell Electron Pair Repulsion) theory
  - Molecular geometries:
    - Linear (180°): 2 electron pairs
    - Trigonal planar (120°): 3 pairs, no lone pairs
    - Bent (<120°): 3 pairs with lone pairs
    - Tetrahedral (109.5°): 4 pairs, no lone pairs
    - Trigonal pyramidal (107°): 4 pairs, 1 lone pair
    - Bent (104.5°): 4 pairs, 2 lone pairs
    - Trigonal bipyramidal: 5 pairs
    - Octahedral (90°): 6 pairs
  - Walsh diagrams: orbital energy vs geometry
  - Examples: H₂O (bent), NH₃ (pyramidal), CH₄ (tetrahedral), CO₂ (linear)
  - Dipole moments: μ = Σqᵢrᵢ

- **The Hückel Molecular Orbital Method**
  - π-electron theory for conjugated systems
  - Hückel Hamiltonian matrix: H_ii = α, H_ij = β (adjacent)
  - Hückel 4n+2 aromaticity rule
  - Aromatic: benzene (6π), naphthalene (10π), cyclopentadienyl⁻ (6π)
  - Antiaromatic: cyclobutadiene (4π)
  - Total π-electron energy: E_π = Σᵢnᵢεᵢ
  - Delocalization (resonance) energy
  - Bond order: p_ij = Σₖnₖc_ikc_jk
  - Charge density: q_i = Σₖnₖ|c_ik|²
  - Aromatic stabilization energy
  - Examples: benzene resonance energy = 2β

**Applications:** Quantum chemistry, computational chemistry, molecular spectroscopy, chemical bonding theory, organic chemistry, materials science, drug design

**Relativistic Quantum Mechanics and Spin** (`physics/relativistic_quantum_mechanics.hpp`):
Comprehensive spin-1/2 theory and relativistic quantum mechanics (~1,174 lines)

**Spin and Atomic Spectra:**

- **Degenerate Position Eigenstates**
  - Degeneracy in quantum systems: g_n for various potentials
  - Hydrogen degeneracy: g_n = 2n² (including spin)
  - 3D isotropic harmonic oscillator: g_N = ½(N+1)(N+2)
  - Good quantum numbers and conserved quantities
  - Accidental degeneracy vs symmetry-based degeneracy
  - Degenerate perturbation theory framework
  - Lifting degeneracy with perturbations

- **Spin-Half Particles**
  - Pauli matrices: σ_x, σ_y, σ_z
  - Spin operators: S⃗ = (ℏ/2)σ⃗
  - Spin eigenstates: |↑⟩ = (1,0)ᵀ, |↓⟩ = (0,1)ᵀ
  - General spin states: |ψ⟩ = cos(θ/2)|↑⟩ + e^(iφ)sin(θ/2)|↓⟩
  - Spin expectation values: ⟨S_x⟩, ⟨S_y⟩, ⟨S_z⟩
  - Bloch sphere representation
  - Density matrices for mixed states: ρ = ½(I + r⃗·σ⃗)
  - Purity: Tr(ρ²), pure states vs mixed states
  - Larmor precession: ω = -γB
  - Time evolution of spin states

- **Spin Magnetic Moment (Stern-Gerlach Experiment)**
  - Magnetic moment: μ⃗ = -g_e(μ_B/ℏ)S⃗
  - Bohr magneton: μ_B = eℏ/2m_e ≈ 9.274×10⁻²⁴ J/T
  - Electron g-factor: g_e ≈ 2.00232 (QED correction)
  - Stern-Gerlach force: F_z = μ_z(∂B_z/∂z)
  - Beam deflection: Δz = (μ_B/m_e)(∂B_z/∂z)t²
  - Sequential Stern-Gerlach experiments
  - Landé g-factor: g_J = 1 + [J(J+1) - L(L+1) + S(S+1)]/[2J(J+1)]
  - Atomic magnetic moments for arbitrary J, L, S

- **Spin-Orbit Coupling**
  - Spin-orbit Hamiltonian: H_SO = (1/2m²c²r)(dV/dr)L⃗·S⃗
  - Fine structure energy: ΔE_SO ∝ Z⁴α⁴/(n³l(l+1/2)(l+1))
  - Total angular momentum: J⃗ = L⃗ + S⃗
  - Possible j values: j = l ± 1/2
  - Fine structure splitting: ΔE_fs between j levels
  - ⟨L⃗·S⃗⟩ expectation value: ½[j(j+1) - l(l+1) - s(s+1)]ℏ²
  - Thomas precession factor 1/2
  - Relativistic origin from Dirac equation

- **Zeeman Effect Revisited**
  - **Normal Zeeman** (no spin): ΔE = μ_B m_l B
  - **Anomalous Zeeman** (with spin): ΔE = g_J μ_B m_J B
  - **Paschen-Back Effect** (strong field): ΔE = μ_B(m_l + 2m_s)B
  - Transition between weak and strong field regimes
  - Selection rules: ΔJ = 0, ±1; Δm_J = 0, ±1 (Δm_J = 0 forbidden for J=0→J=0)
  - Hyperfine structure from nuclear spin
  - Hyperfine splitting: ΔE_hf = ½Ahf[F(F+1) - I(I+1) - J(J+1)]
  - Hydrogen 21cm line: F=1 → F=0 transition at 1420.405 MHz
  - Zeeman splitting of hyperfine levels

**Relativistic Quantum Mechanics:**

- **Relativistic Notation and Formalism**
  - 4-vector components: x^μ = (ct, x, y, z)
  - Minkowski metric tensor: g_μν = diag(±1, ∓1, ∓1, ∓1)
  - D'Alembertian operator: □ = ∂_μ ∂^μ = (1/c²)∂²/∂t² - ∇²
  - Lorentz invariant scalar products: x·y = x^μ y_μ
  - Natural units (ℏ = c = 1) conversions

- **The Klein-Gordon Equation** (Comprehensive Treatment)
  - Klein-Gordon equation: (□ + m²c²/ℏ²)ψ = 0
  - Dispersion relation: ω² = c²k² + (mc²/ℏ)²
  - Plane wave solutions: ψ = Ae^(i(k⃗·r⃗-ωt))
  - Energy-momentum relation: E² = (pc)² + (mc²)²
  - Positive and negative energy solutions: E = ±E_p
  - Conserved current density: j^μ = (iℏ/2m)(ψ*∂^μψ - ψ∂^μψ*)
  - Probability density (not positive definite): ρ = j^0/c
  - Klein paradox: T > 1 for V₀ > E + 2mc² (pair production)
  - Continuity equation: ∂ρ/∂t + ∇·j⃗ = 0

- **Nonrelativistic Limit (Klein-Gordon → Schrödinger)**
  - Ansatz: ψ = φ(x,t)e^(-imc²t/ℏ)
  - Schrödinger equation recovery: iℏ∂φ/∂t = -(ℏ²/2m)∇²φ
  - Relativistic corrections: E ≈ mc² + p²/2m - p⁴/8m³c² + ...
  - Velocity ratio β = v/c = pc/E
  - Validity criterion: p << mc (v << c)

- **Free Spin-0 Particles**
  - General solution: ψ(x,t) = ∫[A(k)e^(i(kx-ωt)) + B(k)e^(i(kx+ωt))]dk
  - Energy eigenvalue: E_p = √((pc)² + (mc²)²)
  - Group velocity: v_g = dω/dk = pc²/E
  - Phase velocity: v_p = ω/k = E/p (can exceed c)
  - Gaussian wave packets: ψ = exp(-x²/4σ²)exp(ik₀x)
  - Klein-Gordon inner product: (ψ₁, ψ₂) = i∫[ψ₁*∂_t ψ₂ - (∂_t ψ₁*)ψ₂]d³x

- **Energy-Momentum Tensor T^μν**
  - Energy density: T^00 = (1/2)[(∂_t ψ)² + c²(∇ψ)² + (mc²/ℏ)²ψ²]
  - Momentum density: T^0i = (∂_t ψ*)(∂_i ψ) + (∂_i ψ*)(∂_t ψ)
  - Stress tensor: T^ij = c²[(∂_i ψ*)(∂_j ψ) + c.c.] - δ^ij L
  - Conservation: ∂_μ T^μν = 0
  - Hamiltonian density: H = T^00

- **Klein-Gordon in Schrödinger Form**
  - Two-component: Ψ = (ψ, π)^T where π = ∂ψ/∂t
  - First-order evolution: iℏ∂Ψ/∂t = H_KG Ψ
  - Hamiltonian matrix: H_KG = [[0, 1], [c²∇² - (mc²)², 0]]
  - Positive-definite norm: ||Ψ||² = ∫[|π|² + c²|∇ψ|² + (mc²/ℏ)²|ψ|²]d³x

- **Charge Conjugation**
  - Charge conjugation operator: C: ψ → ψ*
  - Particle/antiparticle states (positive/negative frequency)
  - C-parity for neutral scalars: C = ±1
  - Current transformation: j^μ → -j^μ under C
  - Self-conjugate states (real scalar fields)

- **Feshbach-Villars Representation**
  - FV transformation: φ = (1/√2)(ψ + iπ/mc²), χ = (1/√2)(ψ - iπ/mc²)
  - Positive-definite density: ρ_FV = |φ|² + |χ|² ≥ 0
  - Coupled equations: iℏ∂φ/∂t = mc²φ - iℏc∇χ, iℏ∂χ/∂t = -mc²χ + iℏc∇φ
  - FV Hamiltonian: H_FV = βmc² + α⃗·(cp⃗)
  - Nonrelativistic limit: χ → 0, φ → ψ_Schrödinger

- **Klein-Gordon with Electromagnetic Field**
  - Minimal coupling: ∂_μ → D_μ = ∂_μ + (iq/ℏc)A_μ
  - Modified equation: [(∂_μ + iqA_μ)(∂^μ + iqA^μ) + (mc/ℏ)²]ψ = 0
  - Energy-momentum relation: (E - qφ)² = (p⃗ - qA⃗)²c² + (mc²)²
  - Current with field: j^μ = (iq/2m)[ψ*(D^μψ) - (D^μψ)*ψ]
  - Landau levels: E_n = √[(mc²)² + 2n|q|ℏcB]
  - Cyclotron frequency: ω_c = |q|B/(γm)

- **Gauge Invariance (U(1))**
  - Wave function transformation: ψ → ψ' = e^(iqΛ/ℏ)ψ
  - Vector potential: A_μ → A'_μ = A_μ - ∂_μΛ
  - Scalar potential: φ → φ' = φ + ∂_t Λ
  - Aharonov-Bohm phase: exp(iq/ℏ ∮A⃗·dl⃗)
  - Field strength (gauge invariant): F_μν = ∂_μA_ν - ∂_νA_μ

- **Nonrelativistic Limit with Fields**
  - Pauli equation: iℏ∂ψ/∂t = [(p⃗ - qA⃗)²/2m + qφ]ψ
  - Darwin term: H_Darwin = -(ℏ²/8m²c²)∇²V
  - Relativistic kinetic correction: -(p⃗ - qA⃗)⁴/8m³c²
  - No spin-orbit coupling (spin-0 particle)
  - Diamagnetic energy: ΔE = -(q²B²r_⊥²)/(8mc²)

- **Interpretation of One-Particle Operators**
  - Position operator: Newton-Wigner (non-local)
  - Momentum: p̂ = -iℏ∇ (well-defined)
  - Energy: Ê = iℏ∂/∂t (positive and negative eigenvalues)
  - Charge density: Not positive-definite (requires second quantization)
  - Current density: Well-defined for Klein-Gordon field
  - Angular momentum: L̂ = r⃗ × p̂ (orbital only, no spin)
  - Second quantization necessity: Negative energy states → antiparticles
  - No Zitterbewegung (unlike Dirac equation)
  - Compton wavelength: λ_C = ℏ/(mc) (localization scale)

- **The Dirac Equation**
  - 4-component Dirac spinors: ψ = (ψ₁, ψ₂, ψ₃, ψ₄)ᵀ
  - Dirac equation: iℏ∂ψ/∂t = (cα⃗·p⃗ + βmc²)ψ
  - Dirac matrices: α_i (4×4), β (4×4)
  - Gamma matrices: γ⁰, γⁱ with anticommutation {γ^μ, γ^ν} = 2g^μν
  - Free particle solutions: u(p) for positive energy, v(p) for negative
  - Positive definite probability density: ρ = ψ†ψ > 0
  - Current density: j⃗ = cψ†α⃗ψ
  - Continuity equation: ∂ρ/∂t + ∇·j⃗ = 0
  - Non-relativistic limit: Pauli equation with spin-orbit coupling
  - Antiparticle interpretation (hole theory)

- **Spin and the Dirac Particle**
  - Intrinsic spin s = 1/2 from Dirac equation
  - Spin angular momentum: |S| = (√3/2)ℏ
  - Helicity operator: h = Σ⃗·p̂ (chirality in massless limit)
  - Helicity eigenvalues: ±1 (right/left-handed)
  - Gyromagnetic ratio: g = 2 (exact prediction from Dirac)
  - QED corrections: g_e ≈ 2.00232 (Schwinger correction)
  - Anomalous magnetic moment: a_e = (g-2)/2 ≈ 0.00116
  - Zitterbewegung (trembling motion):
    - Frequency: ω = 2mc²/ℏ ≈ 10²¹ rad/s
    - Amplitude: λ_C = ℏ/(mc) ≈ 3.86×10⁻¹³ m (Compton wavelength)
  - TBMT (Thomas-Bargmann-Michel-Telegdi) equation for spin precession

- **Spin-Orbit Coupling in the Dirac Hamiltonian**
  - Automatic L⃗·S⃗ coupling from Dirac equation
  - Correct Thomas precession factor 1/2 (not 1)
  - Fine structure from relativistic corrections
  - Darwin term for s-states: ΔE_Darwin = (πℏ²/2m²c²)Z|ψ(0)|²
  - Kinetic energy correction: ΔE_kin = -p⁴/(8m³c²)
  - Total fine structure Hamiltonian
  - Non-relativistic expansion to order (v/c)²

- **The Dirac Hydrogen Atom**
  - Exact Dirac energy levels: E_nj = mc²[1 + (Zα)²/(n - j - 1/2 + √((j+1/2)² - (Zα)²))²]^(-1/2)
  - Fine structure constant: α ≈ 1/137.036
  - Quantum numbers: n (principal), j (total angular momentum), l (orbital)
  - j = l ± 1/2 for given l
  - Fine structure splitting between j states
  - n²S₁/₂, n²P₁/₂, n²P₃/₂ notation
  - Degeneracy: n²S₁/₂ and n²P₁/₂ degenerate in Dirac theory
  - Lamb shift (QED correction): 2S₁/₂ - 2P₁/₂ ≈ 1057 MHz
  - Hydrogen spectrum with fine structure and Lamb shift
  - Vacuum polarization and self-energy corrections

- **The Dirac Particle in a Magnetic Field**
  - Minimal coupling: p⃗ → p⃗ - eA⃗
  - Landau levels for Dirac particles: E_n = ±√((mc²)² + 2n|e|ℏcB)
  - Cyclotron frequency: ω_c = |e|B/(γm)
  - Magnetic length: l_B = √(ℏc/|e|B)
  - Automatic Pauli term: -μ⃗·B⃗ (no ad hoc addition needed)
  - Anomalous magnetic moment from QED
  - Quantum Hall effect foundation
  - Critical magnetic field: B_c = m²c³/(eℏ) ≈ 4.4×10¹³ Gauss
  - Pair production threshold in strong B fields
  - Synchrotron radiation power: P ∝ γ⁴B²

**Applications:** Relativistic quantum mechanics, atomic spectroscopy, spin resonance (ESR/NMR), quantum electrodynamics (QED), high-energy physics, particle physics, astrophysics (pulsars, magnetars), precision measurements (g-2 experiments), relativistic quantum chemistry

## 🚀 Usage

### Integration
Simply include the required header files in your C++ project:

```cpp
#include "maths/fourier_analysis.hpp"
#include "maths/pde_numerical_methods.hpp"
#include "physics/electrostatics.hpp"

using namespace maths;
using namespace physics;
```

### Compilation
All modules are header-only and require C++17:

```bash
g++ -std=c++17 -I./include your_program.cpp -o your_program -lm
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
  - Standard texts on complex analysis (Ahlfors, Rudin)
  - Operator algebra texts (Kadison & Ringrose, Takesaki)
  - Quantum mechanics (Griffiths, Sakurai, Cohen-Tannoudji)
  - Standard texts on stochastic processes and Monte Carlo methods

## 📊 Statistics

- **Mathematics Modules**: 24 header-only modules in flat structure
  - **Complex Analysis** (~1,650 lines): Zeros of holomorphic functions, infinite products, Gamma function, divisors, Blaschke products, Kummer's functions
  - Differential algebra, Fourier analysis, subdifferentials, nonsmooth algorithms
  - Monte Carlo & MCMC methods, stochastic differential equations (SDEs) & Itô calculus
  - Variational calculus, dynamical systems & chaos
  - Partial differential equations (6 modules: classification, solutions, transforms, variational, numerical methods, plus main PDE module)
  - Probability distributions, matrices, vectors, calculus, trigonometry
  - Financial mathematics, actuarial science, econometrics
- **Physics Modules**:
  - Basic: 25+ modules covering classical mechanics, E&M, thermodynamics, optics, modern physics
  - **Quantum Mechanics & Chemistry** (5 comprehensive modules, ~8,906 lines total):
    - **Operator Algebras** (~2,800 lines): von Neumann algebras, unitary representations, factor classification, elementary C*-algebra theory (13 classes), GNS construction
    - **Quantum Foundations** (~1,000 lines): Historical development from Planck to Schrödinger, Bohr model, matrix mechanics, uncertainty relations
    - **Advanced Quantum Mechanics** (~1,650 lines): Kummer's functions, Hamiltonian mechanics, perturbation theory, Stark effect, Pauli exclusion, electron spin, helium atom
    - **Quantum Chemistry** (~1,300 lines): Atomic structure (Hartree-Fock, Slater orbitals, multiplet theory), molecular structure (Born-Oppenheimer, diatomic molecules, H₂⁺, H₂, chemical bonding, VSEPR, Hückel MO theory)
    - **Relativistic Quantum Mechanics** (~2,156 lines): Spin-1/2 theory (Pauli matrices, Bloch sphere, Stern-Gerlach), atomic spectra (spin-orbit coupling, Zeeman effect), comprehensive Klein-Gordon equation (12 topics: notation, equation, nonrelativistic limit, free particles, energy-momentum tensor, Schrödinger form, charge conjugation, Feshbach-Villars, EM fields, gauge invariance, operators interpretation), Dirac equation, Dirac hydrogen atom, Landau levels
  - Advanced: 23+ modules in Hamiltonian mechanics, cosmology, fluid dynamics, gauge theory, QFT
- **Probability Distributions**: 14 distributions (Bernoulli, Binomial, Poisson, Geometric, Negative Binomial, Hypergeometric, Uniform, Normal, Exponential, Gamma, Beta, Chi-squared, Student's t, F-distribution)
- **Key Algorithms**:
  - DFT, FFT (O(N log N))
  - ISTA, FISTA (O(1/k²) convergence)
  - ADMM, Ritt's algorithm
  - Monte Carlo, MCMC (Metropolis-Hastings, Gibbs, HMC)
  - Itô integrals, Itô's lemma, Euler-Maruyama, Milstein methods
  - Kalman filter, optimal stopping, stochastic control
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
  - Galerkin finite element method
  - Rayleigh-Ritz energy minimization
  - Weighted Residual Methods (Galerkin, collocation, subdomain, least squares)
  - Collocation with Chebyshev nodes
  - Upwind, Lax-Friedrichs, Lax-Wendroff schemes
  - ADI (Alternating Direction Implicit)
  - Crank-Nicolson time stepping
  - SOR (Successive Over-Relaxation)
  - Picard iteration
  - Von Neumann stability analysis
  - Multiple scales analysis
  - Matched asymptotic expansions
  - WKB approximation
  - Poincare-Lindstedt method

## 🎓 Educational Value

Each module serves as both:
1. **Production-ready code** for numerical computations
2. **Educational reference** showing how abstract mathematics translates to algorithms
3. **Research tool** for:
   - **Complex Analysis**: Zeros of holomorphic functions, infinite products, special functions (Gamma, Beta), Blaschke products, Hardy spaces
   - **Operator Algebras**: Von Neumann algebras, C*-algebras, GNS construction, spectral theory, quantum observables
   - **Quantum Mechanics**: Historical development, Schrödinger equation, perturbation theory, multi-electron systems, atomic structure
   - **Quantum Chemistry**: Hartree-Fock method, molecular orbital theory, chemical bonding, VSEPR theory, Hückel aromaticity, spectroscopy
   - **Relativistic Quantum Mechanics**: Dirac equation, spin theory, Klein-Gordon equation, fine structure, Zeeman effect, Landau levels, QED corrections
   - Ordinary and partial differential equations
   - Stochastic differential equations and Itô calculus
   - Dynamical systems, chaos theory, and bifurcation analysis
   - Nonsmooth optimization and variational calculus
   - Signal processing and Fourier analysis
   - Monte Carlo methods, MCMC, and stochastic simulations
   - Filtering theory (Kalman filter), optimal stopping, stochastic control
   - Classical and quantum field theory
   - Statistical mechanics and computational physics
   - PDE theory: parabolic (heat), elliptic (Laplace/Poisson), hyperbolic (wave)
   - Transform methods: Laplace and Fourier transforms for PDEs
   - Boundary value problems and initial value problems
   - Green's functions and fundamental solutions
   - Weak formulations and variational methods
   - Finite element methods and Galerkin approximations
   - Finite difference schemes and numerical stability
   - Computational fluid dynamics and heat transfer
   - Mathematical finance and quantitative modeling

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
