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
│   │   ├── algebra/               # Differential Algebra (Ritt-Kolchin theory)
│   │   ├── analysis/              # Analysis (Fourier, subdifferentials)
│   │   ├── calculus/              # Calculus theorems and methods
│   │   ├── probability/           # Probability distributions & statistics
│   │   ├── linear_algebra/        # Matrix operations
│   │   ├── trigonometry/          # Trig functions and identities
│   │   ├── transforms/            # Integral transforms
│   │   ├── finance/               # Financial mathematics
│   │   ├── actuarial/             # Actuarial science
│   │   └── econometrics/          # Econometric models
│   └── physics/
│       └── (physics modules)
├── examples/                       # Demonstration programs
└── README.md
```

## 🧮 Mathematics Modules

### Differential Algebra (`maths/algebra/differential_algebra.hpp`)
**~1700 lines | Chapters I-IX from Ritt's "Differential Algebra"**

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

### Fourier Analysis (`maths/analysis/fourier_analysis.hpp`)
**~1200 lines | Discrete & continuous Fourier theory**

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

### Advanced Subdifferentials (`maths/analysis/advanced_subdifferentials.hpp`)
**~900 lines | Nonsmooth analysis & variational calculus**

- **Clarke Subdifferential** ∂_C f(x)
- **Limiting (Mordukhovich) Subdifferential** ∂_L f(x)
- **Fréchet, Proximal, Graded Subdifferentials**
- **Normal & Tangent Cones**: N_C(x), T_C(x)
- **Coderivatives** D*F for set-valued mappings
- **Calculus Rules**: sum, chain, max rules
- **Metric Regularity Criterion**

### Nonsmooth Algorithms (`maths/analysis/nonsmooth_algorithms.hpp`)
**~800 lines | Optimization algorithms**

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

### Probability & Statistics (`maths/probability/distributions.hpp`)
**~650 lines | Comprehensive probability distributions**

**Discrete Distributions:**
- **Bernoulli**: P(X=1) = p
- **Binomial**: C(n,k) p^k (1-p)^(n-k)
- **Poisson**: λ^k e^(-λ) / k!
- **Geometric**: trials until first success

**Continuous Distributions:**
- **Uniform**: constant density on [a,b]
- **Normal (Gaussian)**: N(μ, σ²) with PDF, CDF, quantile
- **Exponential**: memoryless, rate λ
- **Gamma**: shape α, rate β
- **Beta**: on [0,1], conjugate prior
- **Chi-squared**: χ²(k) for hypothesis testing

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

### Other Mathematics Modules

- **Calculus** (`calculus/theorems.hpp`): Numerical derivatives, integration (Simpson's rule)
- **Trigonometry** (`trigonometry/identities.hpp`): Computational trig identities
- **Linear Algebra**: Matrix operations, eigenvalues
- **Financial Mathematics**: Options pricing, risk metrics
- **Actuarial Science**: Life tables, annuities
- **Econometrics**: Time series, regression

## 🔬 Physics Modules

*(Physics modules follow the same computational pattern)*

## 🎯 Demo Programs

All demos compile with: `g++ -std=c++17 -I./include -o demo examples/demo.cpp -lm`

### Differential Algebra Demos

1. **`differential_algebra_demo`** (~400 lines)
   - Basic differential polynomials and fields
   - Characteristic sets and reduction
   - Differential ideals
   - Manifolds and dimension theory
   - Elimination theory
   - Low power theorem

2. **`differential_algebra_advanced_demo`** (~350 lines)
   - Manifold intersections (Kronecker's theorem)
   - Orthonomic systems (Riquier)
   - Partial differential algebra
   - Classical PDEs: Laplace, wave, heat
   - Cauchy-Kovalevskaya existence

### Analysis Demos

3. **`fourier_analysis_demo`** (~350 lines)
   - DFT vs FFT performance (11x speedup on N=256)
   - Circulant matrices and convolution theorem
   - Haar and Daubechies-4 wavelets
   - Fourier series and multipliers
   - STFT and spectrogram
   - Chirp signal detection

4. **`nonsmooth_algorithms_demo`** (~500 lines)
   - Proximal operators (soft thresholding)
   - Subgradient descent
   - ISTA and FISTA comparison
   - ADMM consensus optimization
   - Sparse signal recovery (L1 minimization)

5. **`advanced_subdifferentials_demo`** (~500 lines)
   - Clarke subdifferential computation
   - Normal and tangent cones
   - Limiting subdifferentials
   - Metric regularity testing
   - Nonsmooth optimization examples

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
  - Standard probability theory texts

## 📊 Statistics

- **Total Lines**: ~10,000+ lines of computational mathematics
- **Modules**: 15+ major mathematical areas
- **Demos**: 10+ comprehensive demonstration programs
- **Distributions**: 10 probability distributions with full statistics
- **Algorithms**: DFT, FFT, ISTA, FISTA, ADMM, Ritt's algorithm, and more

## 🎓 Educational Value

Each module serves as both:
1. **Production-ready code** for numerical computations
2. **Educational reference** showing how abstract mathematics translates to algorithms
3. **Research tool** for differential equations, nonsmooth optimization, signal processing

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
