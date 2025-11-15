#ifndef PHYSICS_LOOP_QUANTUM_GRAVITY_HPP
#define PHYSICS_LOOP_QUANTUM_GRAVITY_HPP

#include <vector>
#include <complex>
#include <functional>
#include <cmath>
#include <array>
#include <map>
#include <set>
#include <string>

namespace physics {
namespace loop_quantum_gravity {

using Complex = std::complex<double>;

/**
 * @brief Physical constants for quantum gravity
 */
namespace constants {
    constexpr double c = 299792458.0;           // Speed of light (m/s)
    constexpr double G = 6.67430e-11;           // Gravitational constant (m³/kg·s²)
    constexpr double hbar = 1.054571817e-34;    // Reduced Planck constant (J·s)
    constexpr double k_B = 1.380649e-23;        // Boltzmann constant (J/K)

    // Derived constants
    constexpr double l_P = 1.616255e-35;        // Planck length (m): √(ℏG/c³)
    constexpr double t_P = 5.391247e-44;        // Planck time (s): l_P/c
    constexpr double m_P = 2.176434e-8;         // Planck mass (kg): √(ℏc/G)
    constexpr double E_P = 1.956e9;             // Planck energy (J): m_P c²
    constexpr double T_P = 1.416784e32;         // Planck temperature (K): E_P/k_B

    // LQG specific
    constexpr double gamma = 0.2375;            // Barbero-Immirzi parameter (typical value)
}

/**
 * @brief Quantum Space Structure
 *
 * Overview of quantum gravity and quantum space structure
 */
class QuantumSpaceStructure {
public:
    /**
     * @brief Quantum gravity program
     *
     * Background-independent quantization of general relativity
     */
    static std::string quantum_gravity() {
        return "LQG: background-independent quantum theory of spacetime geometry";
    }

    /**
     * @brief Planck scale discreteness
     *
     * Spacetime has discrete structure at Planck scale
     */
    static std::string planck_scale() {
        return "Space quantized at l_P = √(ℏG/c³) ≈ 1.6×10⁻³⁵ m";
    }

    /**
     * @brief Planck length
     *
     * Fundamental length scale of quantum gravity
     */
    static double planck_length() {
        return std::sqrt(constants::hbar * constants::G /
                        (constants::c * constants::c * constants::c));
    }

    /**
     * @brief Planck area
     *
     * Fundamental area quantum A_P = l_P²
     */
    static double planck_area() {
        double l_p = planck_length();
        return l_p * l_p;
    }

    /**
     * @brief Planck volume
     *
     * Fundamental volume quantum V_P = l_P³
     */
    static double planck_volume() {
        double l_p = planck_length();
        return l_p * l_p * l_p;
    }

    /**
     * @brief Main features of LQG
     *
     * Background independence, diffeomorphism invariance, spin networks
     */
    static std::string main_features() {
        return "Background independence, diff invariance, discrete geometry, spin networks";
    }

    /**
     * @brief Connection dynamics
     *
     * Uses Ashtekar-Barbero connection variables
     */
    static std::string ashtekar_variables() {
        return "Variables: Ashtekar connection A^i_a, densitized triad E^a_i";
    }

    /**
     * @brief Quantum geometry
     *
     * Geometry quantized: area and volume have discrete spectra
     */
    static std::string quantum_geometry() {
        return "Area: A = 8πγl_P²√(j(j+1)), Volume: discrete spectrum";
    }

    /**
     * @brief No singularities
     *
     * Quantum geometry resolves classical singularities
     */
    static std::string singularity_resolution() {
        return "Singularities replaced by quantum bounces (Big Bounce)";
    }
};

/**
 * @brief Kinematical State Space 𝓚
 *
 * The fundamental Hilbert space before constraints
 */
class KinematicalStateSpace {
public:
    /**
     * @brief Definition of 𝓚
     *
     * Space of cylindrical functions on connections
     */
    static std::string definition() {
        return "𝓚 = space of cylindrical functions Ψ[A] on 𝒜/𝒢";
    }

    /**
     * @brief Configuration space
     *
     * 𝒜 = space of SU(2) connections on spatial manifold Σ
     */
    static std::string configuration_space() {
        return "𝒜 = {A^i_a | A: Σ → su(2) connection}";
    }

    /**
     * @brief Gauge group
     *
     * 𝒢 = group of SU(2) gauge transformations
     */
    static std::string gauge_group() {
        return "𝒢 = {g: Σ → SU(2) | gauge transformations}";
    }

    /**
     * @brief Cylindrical functions
     *
     * Functions depending on A only via holonomies along finite graph
     */
    static std::string cylindrical_functions() {
        return "Ψ_γ[A] = f(h_e₁[A], ..., h_eₙ[A]) on graph γ";
    }

    /**
     * @brief Holonomy
     *
     * Parallel transport h_e[A] ∈ SU(2) along edge e
     */
    static std::string holonomy() {
        return "h_e[A] = 𝒫 exp(∫_e A) ∈ SU(2)";
    }

    /**
     * @brief Ashtekar-Lewandowski measure
     *
     * Unique diffeomorphism-invariant measure on 𝒜/𝒢
     */
    static std::string al_measure() {
        return "dμ_AL: unique diff-invariant measure on 𝒜/𝒢";
    }

    /**
     * @brief Scalar product on 𝓚
     *
     * ⟨Ψ₁|Ψ₂⟩ = ∫ Ψ₁*[A] Ψ₂[A] dμ_AL[A]
     */
    static std::string scalar_product() {
        return "⟨Ψ₁|Ψ₂⟩ = ∫ Ψ₁*[A] Ψ₂[A] dμ_AL[A]";
    }

    /**
     * @brief Kinematical vs physical
     *
     * 𝓚 is kinematical; physical states satisfy constraints
     */
    static std::string kinematical_vs_physical() {
        return "𝓚 kinematical, 𝓗_phys ⊂ 𝓚 satisfies Gauss, diff, Hamiltonian constraints";
    }
};

/**
 * @brief Structures in 𝓚
 *
 * Mathematical structures of the kinematical Hilbert space
 */
class StructuresInK {
public:
    /**
     * @brief Graph structure
     *
     * States labeled by embedded graphs γ in Σ
     */
    static std::string graph_structure() {
        return "γ = (V, E): graph with vertices V, edges E embedded in Σ";
    }

    /**
     * @brief Orthonormal basis
     *
     * Spin network states form overcomplete basis
     */
    static std::string basis() {
        return "Spin network states |γ, j_e, i_v⟩ form orthonormal basis";
    }

    /**
     * @brief Decomposition
     *
     * 𝓚 = ⊕_γ 𝓚_γ (direct sum over graphs)
     */
    static std::string decomposition() {
        return "𝓚 = ⊕_γ 𝓚_γ where 𝓚_γ = L²(SU(2)^|E|, dμ_Haar)";
    }

    /**
     * @brief Separability
     *
     * Each 𝓚_γ is separable, but 𝓚 is non-separable
     */
    static std::string separability() {
        return "𝓚_γ separable, but 𝓚 non-separable (uncountable sum)";
    }

    /**
     * @brief Haar measure
     *
     * Induced measure on holonomies
     */
    static std::string haar_measure() {
        return "dμ_Haar: product of Haar measures on SU(2) for each edge";
    }

    /**
     * @brief Peter-Weyl theorem
     *
     * Decomposition of L²(SU(2)) into irreps
     */
    static std::string peter_weyl() {
        return "L²(SU(2)) = ⊕_{j=0,1/2,1,...} V_j ⊗ V_j* (j = spin)";
    }

    /**
     * @brief Representation labels
     *
     * j = 0, 1/2, 1, 3/2, 2, ... (SU(2) spins)
     */
    static std::string spin_labels() {
        return "j ∈ {0, 1/2, 1, 3/2, 2, ...} labels irreps of SU(2)";
    }
};

/**
 * @brief Invariances of the Scalar Product
 *
 * Symmetries preserved by the AL measure
 */
class ScalarProductInvariances {
public:
    /**
     * @brief Gauge invariance
     *
     * ⟨Ψ₁|Ψ₂⟩ invariant under SU(2) gauge transformations
     */
    static std::string gauge_invariance() {
        return "⟨U_g Ψ₁|U_g Ψ₂⟩ = ⟨Ψ₁|Ψ₂⟩ for g ∈ 𝒢";
    }

    /**
     * @brief Diffeomorphism invariance
     *
     * ⟨Ψ₁|Ψ₂⟩ invariant under spatial diffeomorphisms
     */
    static std::string diff_invariance() {
        return "⟨U_φ Ψ₁|U_φ Ψ₂⟩ = ⟨Ψ₁|Ψ₂⟩ for φ ∈ Diff(Σ)";
    }

    /**
     * @brief Uniqueness
     *
     * AL measure is unique with these invariances
     */
    static std::string uniqueness() {
        return "dμ_AL unique measure invariant under 𝒢 and Diff(Σ)";
    }

    /**
     * @brief Analytic continuation
     *
     * Extends Fock-like constructions to non-perturbative regime
     */
    static std::string non_perturbative() {
        return "AL measure: non-perturbative, no background metric needed";
    }
};

/**
 * @brief Internal Gauge Invariance and 𝓚₀
 *
 * Gauge-invariant subspace
 */
class InternalGaugeInvariance {
public:
    /**
     * @brief Gauss constraint
     *
     * Generator of SU(2) gauge transformations
     */
    static std::string gauss_constraint() {
        return "Ĝ_i[Λ]Ψ = 0 (gauge invariance at each point)";
    }

    /**
     * @brief Gauge-invariant space
     *
     * 𝓚₀ = {Ψ ∈ 𝓚 | Ĝ_i[Λ]Ψ = 0 ∀Λ}
     */
    static std::string k_zero() {
        return "𝓚₀ = space of gauge-invariant states in 𝓚";
    }

    /**
     * @brief Spin network basis
     *
     * Gauge-invariant states = spin networks with intertwiners
     */
    static std::string spin_network_basis() {
        return "Basis of 𝓚₀: spin networks |γ, j_e, i_v⟩";
    }

    /**
     * @brief Intertwiner requirement
     *
     * Vertices must carry intertwiners between edge spins
     */
    static std::string intertwiner() {
        return "i_v: Inv(⊗_{e∈v} V_{j_e}) (gauge invariant at vertex v)";
    }

    /**
     * @brief Dimension of intertwiner space
     *
     * dim Inv = 1 for 3-valent with j₁+j₂≥j₃, 0 otherwise
     */
    static int intertwiner_dimension_3valent(int twice_j1, int twice_j2, int twice_j3) {
        // Check triangle inequality
        int j1 = twice_j1, j2 = twice_j2, j3 = twice_j3;
        if (j1 + j2 >= j3 && j2 + j3 >= j1 && j3 + j1 >= j2) {
            // Check parity
            if ((j1 + j2 + j3) % 2 == 0) {
                return 1;
            }
        }
        return 0;
    }

    /**
     * @brief Recoupling theory
     *
     * For n-valent vertices: use 6j symbols, recoupling
     */
    static std::string recoupling() {
        return "n-valent: dim Inv given by recoupling theory (6j, 9j symbols)";
    }
};

/**
 * @brief Spin Network States
 *
 * Fundamental quantum states of geometry
 */
class SpinNetworkStates {
public:
    /**
     * @brief Spin network definition
     *
     * |s⟩ = |γ, j_e, i_v⟩: graph + spins + intertwiners
     */
    static std::string definition() {
        return "|s⟩ = |γ, {j_e}, {i_v}⟩: graph γ, spins j_e on edges, intertwiners i_v";
    }

    /**
     * @brief Graph γ
     *
     * Embedded graph in spatial manifold Σ
     */
    static std::string graph() {
        return "γ = (V, E) embedded in Σ (vertices V, edges E)";
    }

    /**
     * @brief Edge labels
     *
     * j_e ∈ {0, 1/2, 1, 3/2, ...} for each edge e
     */
    static std::string edge_labels() {
        return "j_e ∈ ℕ/2: SU(2) spin on edge e";
    }

    /**
     * @brief Vertex labels
     *
     * i_v: intertwiner at vertex v
     */
    static std::string vertex_labels() {
        return "i_v ∈ Inv(⊗_{e∈v} V_{j_e}): gauge invariant tensor";
    }

    /**
     * @brief Orthonormality
     *
     * ⟨s|s'⟩ = δ_{γγ'} δ_{jj'} δ_{ii'}
     */
    static std::string orthonormality() {
        return "⟨s|s'⟩ = δ_{γγ'} δ_{jj'} δ_{ii'} (discrete orthonormality)";
    }

    /**
     * @brief Completeness
     *
     * Σ_s |s⟩⟨s| = I on 𝓚₀
     */
    static std::string completeness() {
        return "Σ_{s} |s⟩⟨s| = I_𝓚₀ (overcomplete in 𝓚)";
    }

    /**
     * @brief Valence of vertex
     *
     * Number of edges meeting at vertex
     */
    static int valence(const std::vector<int>& edges_at_vertex) {
        return static_cast<int>(edges_at_vertex.size());
    }

    /**
     * @brief Physical interpretation
     *
     * Spin networks = quantum states of 3-geometry
     */
    static std::string interpretation() {
        return "Spin network = quantum excitation of spatial geometry";
    }

    /**
     * @brief Edges carry area
     *
     * Edge with spin j piercing surface S contributes to area
     */
    static std::string area_contribution() {
        return "Edge j contributes A_j = 8πγl_P²√(j(j+1)) to area";
    }

    /**
     * @brief Vertices carry volume
     *
     * Volume concentrated at vertices
     */
    static std::string volume_at_vertices() {
        return "Volume V concentrated at vertices (complex formula)";
    }
};

/**
 * @brief Details About Spin Networks
 *
 * Mathematical details and properties
 */
class SpinNetworkDetails {
public:
    /**
     * @brief Wigner D-matrices
     *
     * Matrix elements D^j_mn(g) of SU(2) representation
     */
    static std::string wigner_d() {
        return "D^j_mn(g): (2j+1)×(2j+1) matrix of irrep j";
    }

    /**
     * @brief Spin network function
     *
     * Ψ_s[A] = Tr[D^{j₁}(h_{e₁}) ⊗ ... ⊗ D^{jₙ}(h_{eₙ}) · i_v]
     */
    static std::string wave_function() {
        return "Ψ_s[A] = contraction of D^j(h_e) with intertwiners i_v";
    }

    /**
     * @brief 3j symbols
     *
     * Clebsch-Gordan coefficients for 3-valent vertices
     */
    static std::string three_j() {
        return "(j₁ j₂ j₃; m₁ m₂ m₃): CG coefficients for 3-vertex";
    }

    /**
     * @brief 6j symbols
     *
     * Recoupling coefficients for 4-valent vertices
     */
    static std::string six_j() {
        return "{j₁ j₂ j₃; j₄ j₅ j₆}: recoupling for 4-vertex";
    }

    /**
     * @brief Penrose binor calculus
     *
     * Graphical notation for spin network calculations
     */
    static std::string penrose_notation() {
        return "Penrose graphical calculus: diagrammatic spin network computation";
    }

    /**
     * @brief Fusion rules
     *
     * j₁ ⊗ j₂ = |j₁-j₂| ⊕ ... ⊕ (j₁+j₂)
     */
    static std::string fusion() {
        return "j₁ ⊗ j₂ = ⊕_{j=|j₁-j₂|}^{j₁+j₂} j (triangle inequality)";
    }

    /**
     * @brief Quantum numbers
     *
     * Magnetic quantum numbers m = -j, ..., +j
     */
    static std::string magnetic_numbers() {
        return "m ∈ {-j, -j+1, ..., j-1, j} (2j+1 values)";
    }

    /**
     * @brief Embedding dependence
     *
     * Depends on embedding of γ in Σ (knotting, linking)
     */
    static std::string embedding() {
        return "Spin network depends on embedding: knots, links matter";
    }
};

/**
 * @brief Diffeomorphism Invariance and 𝓚_Diff
 *
 * Diffeomorphism-invariant states
 */
class DiffeomorphismInvariance {
public:
    /**
     * @brief Diffeomorphism constraint
     *
     * Generator of spatial diffeomorphisms
     */
    static std::string diff_constraint() {
        return "D̂_a[N^a]Ψ = 0 (diff invariance)";
    }

    /**
     * @brief Diff-invariant space
     *
     * 𝓚_Diff = {Ψ ∈ 𝓚₀ | D̂_a[N^a]Ψ = 0 ∀N^a}
     */
    static std::string k_diff() {
        return "𝓚_Diff = space of gauge + diff invariant states";
    }

    /**
     * @brief Diffeomorphism action
     *
     * φ ∈ Diff(Σ) acts by φ*: γ → φ(γ)
     */
    static std::string diff_action() {
        return "φ*: |γ, j, i⟩ → |φ(γ), j, i⟩ (pushforward of graph)";
    }

    /**
     * @brief Quotient space
     *
     * 𝓚_Diff = 𝓚₀ / Diff(Σ)
     */
    static std::string quotient() {
        return "𝓚_Diff = 𝓚₀ / Diff(Σ) (equivalence classes)";
    }

    /**
     * @brief Knot classes
     *
     * States labeled by diffeomorphism equivalence of graphs
     */
    static std::string knot_classes() {
        return "States labeled by knot/link classes [γ]_Diff";
    }

    /**
     * @brief Abstract graphs
     *
     * Only combinatorial structure matters, not embedding
     */
    static std::string abstract() {
        return "Up to diff: only abstract graph + coloring matters";
    }

    /**
     * @brief Separability
     *
     * 𝓚_Diff is separable (countable basis)
     */
    static std::string separability() {
        return "𝓚_Diff separable: countable diff classes [γ]";
    }
};

/**
 * @brief Knots and s-Knot States
 *
 * Knot theory in LQG
 */
class KnotStates {
public:
    /**
     * @brief s-knot definition
     *
     * Spin network up to diffeomorphisms
     */
    static std::string s_knot() {
        return "s-knot = [γ, j, i]_Diff (diff equivalence class)";
    }

    /**
     * @brief Knot invariants
     *
     * Colored knot/link invariants
     */
    static std::string knot_invariants() {
        return "s-knots → colored Jones polynomials, Kauffman brackets";
    }

    /**
     * @brief Ambient isotopy
     *
     * Graphs related by continuous deformation
     */
    static std::string isotopy() {
        return "Ambient isotopy: continuous deformation in Σ";
    }

    /**
     * @brief Abstract spin networks
     *
     * Combinatorial data: graph + spins + intertwiners
     */
    static std::string abstract_sn() {
        return "Abstract SN: combinatorial graph, no embedding";
    }

    /**
     * @brief Linking and knotting
     *
     * Different embeddings → different physical states
     */
    static std::string linking() {
        return "Linking/knotting: physical (different quantum states)";
    }

    /**
     * @brief Trivalent decomposition
     *
     * Any graph can be decomposed into trivalent vertices
     */
    static std::string trivalent() {
        return "Any n-valent → decompose to 3-valent (standard form)";
    }

    /**
     * @brief Turaev-Viro model
     *
     * Connection to 3D TQFT
     */
    static std::string turaev_viro() {
        return "Related to Turaev-Viro 3D TQFT with q = root of unity";
    }
};

/**
 * @brief The Hilbert Space 𝓚_Diff is Separable
 *
 * Proof of separability
 */
class KDiffSeparability {
public:
    /**
     * @brief Countability argument
     *
     * Countably many abstract graphs
     */
    static std::string countability() {
        return "Abstract graphs: countable (finite V, E, combinatorial)";
    }

    /**
     * @brief Coloring count
     *
     * Finitely many colorings for bounded total spin
     */
    static std::string colorings() {
        return "Colorings: finite for Σ_e j_e ≤ J (bounded)";
    }

    /**
     * @brief Basis enumeration
     *
     * Enumerate by total spin J
     */
    static std::string enumeration() {
        return "Basis: ⋃_{J=0}^∞ {s-knots with Σj ≤ J} (countable)";
    }

    /**
     * @brief Physical significance
     *
     * Separable → standard quantum mechanics applies
     */
    static std::string physical() {
        return "Separable → well-defined quantum theory";
    }

    /**
     * @brief Contrast with 𝓚
     *
     * 𝓚 non-separable, but 𝓚_Diff separable
     */
    static std::string contrast() {
        return "𝓚 non-separable, but 𝓚_Diff separable (diff constraint!)";
    }
};

/**
 * @brief The Connection Operator Â
 *
 * Connection as quantum operator
 */
class ConnectionOperator {
public:
    /**
     * @brief Configuration variable
     *
     * A^i_a: Ashtekar-Barbero connection
     */
    static std::string variable() {
        return "A^i_a: SU(2) connection (configuration variable)";
    }

    /**
     * @brief Operator action
     *
     * Multiplication operator on wave functionals
     */
    static std::string action() {
        return "Â^i_a(x) Ψ[A] = A^i_a(x) Ψ[A] (multiplication)";
    }

    /**
     * @brief Non-well-defined
     *
     * A(x) not well-defined on 𝓚 (distributional)
     */
    static std::string distributional() {
        return "A(x) ill-defined on 𝓚 (distribution, not function)";
    }

    /**
     * @brief Smeared operator
     *
     * Well-defined: A(S) = ∫_S A (surface integral)
     */
    static std::string smeared() {
        return "Â(S) = ∫_S A^i_a ε^a dΣ well-defined on spin networks";
    }

    /**
     * @brief Holonomy as fundamental
     *
     * h_e[A] = 𝒫 exp(∫_e A) is well-defined
     */
    static std::string holonomy() {
        return "h_e[A] ∈ SU(2) well-defined (holonomy = fundamental)";
    }

    /**
     * @brief Polymer representation
     *
     * States concentrated on graphs (polymer-like)
     */
    static std::string polymer() {
        return "Polymer rep: states on graphs (not Fock representation)";
    }
};

/**
 * @brief The Conjugate Momentum Operator Ê
 *
 * Densitized triad as momentum
 */
class ConjugateMomentumOperator {
public:
    /**
     * @brief Momentum variable
     *
     * E^a_i: densitized triad (momentum conjugate to A)
     */
    static std::string variable() {
        return "E^a_i: densitized triad E = det(e) e^a_i (momentum)";
    }

    /**
     * @brief Poisson bracket
     *
     * {A^i_a(x), E^b_j(y)} = δ^i_j δ^b_a δ³(x,y)
     */
    static std::string poisson_bracket() {
        return "{A^i_a(x), E^b_j(y)} = δ^i_j δ^b_a δ³(x-y)";
    }

    /**
     * @brief Quantum operator
     *
     * Ê = -iℏ δ/δA (functional derivative)
     */
    static std::string operator_form() {
        return "Ê^a_i(x) = -iℏ δ/δA^i_a(x) on Ψ[A]";
    }

    /**
     * @brief Well-defined on graphs
     *
     * E acts by derivative on holonomies
     */
    static std::string graph_action() {
        return "Ê acts by derivative: d/dh on holonomies h_e";
    }

    /**
     * @brief Flux operator
     *
     * E(S) = ∫_S E^a_i n_a well-defined
     */
    static std::string flux() {
        return "Flux Ê(S,f) = ∫_S E^a_i f^i n_a through surface S";
    }

    /**
     * @brief Commutation relation
     *
     * [Â, Ê] = iℏ (canonical quantization)
     */
    static std::string commutator() {
        return "[Â(S), Ê(S')] ∼ iℏ (depends on S ∩ S')";
    }

    /**
     * @brief Geometric meaning
     *
     * E determines metric: q_ab = E^i_a E^i_b / det(E)
     */
    static std::string metric() {
        return "Metric: q_ab ∼ E^i_a E_i^b (E encodes 3-geometry)";
    }
};

/**
 * @brief The Operator Â(S)
 *
 * Connection operator on surfaces
 */
class ConnectionOperatorOnSurface {
public:
    /**
     * @brief Definition
     *
     * Â(S) = ∫_S A ∧ dΣ (surface integral)
     */
    static std::string definition() {
        return "Â(S,f) = ∫_S A^i_a f^i ε^a dΣ (smeared connection)";
    }

    /**
     * @brief Action on spin networks
     *
     * Inserts punctures where S intersects γ
     */
    static std::string action_on_sn() {
        return "Â(S)|γ,j,i⟩ ∼ Σ_{p∈S∩γ} τ^i |γ,j,i⟩ (Pauli matrices at punctures)";
    }

    /**
     * @brief Creates edges
     *
     * Can create new edges piercing S
     */
    static std::string creates_edges() {
        return "Â(S) can create edges: |γ⟩ → |γ ∪ e_new⟩";
    }

    /**
     * @brief Gauge transformation
     *
     * A(S) generates gauge transformations on S
     */
    static std::string gauge_generator() {
        return "Generates SU(2) rotation of spins at S ∩ γ";
    }
};

/**
 * @brief Quanta of Area
 *
 * Area operator eigenvalues
 */
class QuantaOfArea {
public:
    /**
     * @brief Area operator
     *
     * Â(S) = Σ_{p∈S∩γ} √(Ê(S,p)² + ...)
     */
    static std::string operator_def() {
        return "Â(S) = Σ_{p∈S∩γ} √(E^a_i(p) E^b_j(p) n_a n_b) (sum over punctures)";
    }

    /**
     * @brief Eigenvalue formula
     *
     * A = 8πγl_P² Σ_p √(j_p(j_p+1))
     */
    static double area_eigenvalue(const std::vector<int>& twice_j_punctures, double gamma_IP) {
        double l_P = constants::l_P;
        double area = 0.0;
        for (int twice_j : twice_j_punctures) {
            double j = twice_j / 2.0;
            area += std::sqrt(j * (j + 1.0));
        }
        return 8.0 * M_PI * gamma_IP * l_P * l_P * area;
    }

    /**
     * @brief Discrete spectrum
     *
     * Area has discrete eigenvalues (quantum geometry!)
     */
    static std::string discrete_spectrum() {
        return "Area spectrum discrete: A_j = 8πγl_P²√(j(j+1))";
    }

    /**
     * @brief Minimal area
     *
     * Smallest non-zero area (j = 1/2)
     */
    static double minimal_area(double gamma_IP) {
        double l_P = constants::l_P;
        double j = 0.5;
        return 8.0 * M_PI * gamma_IP * l_P * l_P * std::sqrt(j * (j + 1.0));
    }

    /**
     * @brief Area gap
     *
     * ΔA ∼ l_P² (Planck area)
     */
    static std::string area_gap() {
        return "Area gap ΔA ∼ l_P² (quantized in Planck units)";
    }

    /**
     * @brief Black hole entropy
     *
     * S_BH = A/(4l_P²) ~ number of punctures
     */
    static std::string black_hole_entropy() {
        return "S_BH = (A_horizon)/(4γl_P²) ~ N_punctures (Bekenstein-Hawking)";
    }

    /**
     * @brief Barbero-Immirzi parameter
     *
     * γ fixed by black hole entropy
     */
    static std::string immirzi() {
        return "γ ≈ 0.2375 fixed by matching S_BH = A/(4l_P²)";
    }
};

/**
 * @brief n-hand Operators and Recoupling Theory
 *
 * Multi-valent vertices and recoupling
 */
class RecouplingTheory {
public:
    /**
     * @brief n-valent vertices
     *
     * Vertices with n edges (n ≥ 3)
     */
    static std::string n_valent() {
        return "n-valent vertex: n edges meeting at v";
    }

    /**
     * @brief Intertwiner space dimension
     *
     * dim Inv(j₁ ⊗ ... ⊗ jₙ)
     */
    static std::string intertwiner_space() {
        return "dim Inv(V_{j₁} ⊗ ... ⊗ V_{jₙ}) via recoupling";
    }

    /**
     * @brief 6j symbols
     *
     * Wigner 6j: {j₁ j₂ j₃; j₄ j₅ j₆}
     */
    static std::string six_j() {
        return "{j₁ j₂ j₃; j₄ j₅ j₆}: recoupling coefficient for 4-valent";
    }

    /**
     * @brief 9j symbols
     *
     * Used for higher-valent vertices
     */
    static std::string nine_j() {
        return "{j₁ j₂ j₃; j₄ j₅ j₆; j₇ j₈ j₉}: 5-valent and beyond";
    }

    /**
     * @brief Recoupling basis
     *
     * Different ways to couple spins
     */
    static std::string basis_choice() {
        return "Intertwiner basis: choice of recoupling tree";
    }

    /**
     * @brief Racah formula
     *
     * Explicit expressions for 6j symbols
     */
    static std::string racah() {
        return "6j computed via Racah formula (sums over factorials)";
    }

    /**
     * @brief Tetrahedral symmetry
     *
     * 6j symbols have tetrahedral symmetry
     */
    static std::string symmetry() {
        return "6j: invariant under 24 tetrahedral permutations";
    }
};

/**
 * @brief Degenerate Sector
 *
 * Degenerate eigenvalues and gauge fixing
 */
class DegenerateSector {
public:
    /**
     * @brief Area degeneracy
     *
     * Different {j_p} giving same area
     */
    static std::string area_degeneracy() {
        return "Many spin configurations → same area eigenvalue";
    }

    /**
     * @brief Volume degeneracy
     *
     * Huge degeneracy in volume spectrum
     */
    static std::string volume_degeneracy() {
        return "Volume: exponentially large degeneracy";
    }

    /**
     * @brief Gauge degrees of freedom
     *
     * Intertwiners = gauge-invariant DOF
     */
    static std::string gauge_dof() {
        return "Intertwiner quantum numbers resolve degeneracy";
    }

    /**
     * @brief Entropy from degeneracy
     *
     * S = ln(degeneracy) (statistical entropy)
     */
    static std::string statistical_entropy() {
        return "S = k_B ln Ω (Ω = degeneracy of macrostate)";
    }

    /**
     * @brief Black hole microstates
     *
     * Horizon punctures → entropy
     */
    static std::string bh_microstates() {
        return "BH microstates: spin network punctures on horizon";
    }
};

/**
 * @brief Quanta of Volume
 *
 * Volume operator eigenvalues
 */
class QuantaOfVolume {
public:
    /**
     * @brief Volume operator
     *
     * V̂(R) for region R (complicated formula)
     */
    static std::string operator_def() {
        return "V̂(R) = Σ_{v∈R} V̂_v (sum over vertices in R)";
    }

    /**
     * @brief Volume at vertex
     *
     * V_v depends on spins {j_e} meeting at v
     */
    static std::string vertex_volume() {
        return "V_v = f(j₁,...,j_n,i_v) (complex function of spins/intertwiners)";
    }

    /**
     * @brief Discrete spectrum
     *
     * Volume has discrete eigenvalues V_n ∼ l_P³
     */
    static std::string discrete_spectrum() {
        return "Volume spectrum discrete: V ∼ n l_P³ (n depends on graph)";
    }

    /**
     * @brief Minimal volume
     *
     * Smallest quantum: V_min ∼ l_P³
     */
    static double minimal_volume() {
        return constants::l_P * constants::l_P * constants::l_P;
    }

    /**
     * @brief Volume gap
     *
     * ΔV ∼ l_P³ (Planck volume)
     */
    static std::string volume_gap() {
        return "Volume gap ΔV ∼ l_P³ (no arbitrary small volumes)";
    }

    /**
     * @brief Singularity resolution
     *
     * V > 0 always (no V = 0 classical singularity)
     */
    static std::string no_singularity() {
        return "V bounded below: no classical singularities (Big Bounce)";
    }

    /**
     * @brief Rovelli-Smolin formula
     *
     * Original volume formula (simplified)
     */
    static std::string rovelli_smolin() {
        return "V_v ∼ (l_P³/96√2)|ε^{ijk} Σ_e j_e^(i) j_e^(j) j_e^(k)|";
    }
};

/**
 * @brief Quantum Geometry
 *
 * The geometric interpretation of spin networks
 */
class QuantumGeometry {
public:
    /**
     * @brief Discrete geometry
     *
     * Geometry built from quanta
     */
    static std::string discrete() {
        return "3-geometry = discrete chunks (area/volume quanta)";
    }

    /**
     * @brief Graph as geometry
     *
     * Spin network = quantum state of 3-geometry
     */
    static std::string graph_geometry() {
        return "Graph γ = skeleton of quantum geometry";
    }

    /**
     * @brief Edges = area quanta
     *
     * Edges carry quantized area
     */
    static std::string edges() {
        return "Edges: area quanta A_j = 8πγl_P²√(j(j+1))";
    }

    /**
     * @brief Vertices = volume quanta
     *
     * Vertices carry quantized volume
     */
    static std::string vertices() {
        return "Vertices: volume quanta V_v ∼ l_P³ f({j_e})";
    }

    /**
     * @brief Continuum limit
     *
     * Smooth geometry from fine-grained networks
     */
    static std::string continuum() {
        return "Continuum limit: dense spin network → smooth metric";
    }

    /**
     * @brief Weave states
     *
     * Semiclassical states approximating smooth geometry
     */
    static std::string weaves() {
        return "Weave: fine-grained network approximating classical q_ab";
    }

    /**
     * @brief Polymer-like
     *
     * Space has polymer structure at Planck scale
     */
    static std::string polymer() {
        return "Polymer structure: space = network of Planck-scale chunks";
    }

    /**
     * @brief No background
     *
     * No pre-existing space; geometry IS quantum state
     */
    static std::string background_independence() {
        return "Background independent: no a priori space";
    }
};

/**
 * @brief The Texture of Space: Weaves
 *
 * Weave states as semiclassical approximations
 */
class Weaves {
public:
    /**
     * @brief Weave definition
     *
     * Fine-grained spin network approximating smooth geometry
     */
    static std::string definition() {
        return "Weave: spin network with fine mesh ~ l_P, many edges";
    }

    /**
     * @brief Classical limit
     *
     * Weave state → classical 3-metric in semiclassical limit
     */
    static std::string classical_limit() {
        return "Weave → smooth q_ab as l_P/L → 0 (coarse graining)";
    }

    /**
     * @brief Coarse graining
     *
     * Average over regions ΔV >> l_P³
     */
    static std::string coarse_graining() {
        return "Coarse grain: ⟨q_ab⟩_ΔV ≈ classical metric (for ΔV >> l_P³)";
    }

    /**
     * @brief Mesh size
     *
     * Characteristic spacing ε between edges
     */
    static std::string mesh_size() {
        return "Mesh ε: typical edge separation (ε ~ l_P for Planck-scale weave)";
    }

    /**
     * @brief Coherent states
     *
     * Weaves as coherent states (peaked on classical geometry)
     */
    static std::string coherent_states() {
        return "Weaves ~ coherent states for quantum geometry";
    }

    /**
     * @brief Fluctuations
     *
     * Quantum fluctuations δq ~ l_P²/ε²
     */
    static std::string fluctuations() {
        return "Quantum fluctuations: δq_ab ~ (l_P/ε)² (suppressed for ε >> l_P)";
    }

    /**
     * @brief Polymer Planck lattice
     *
     * Dense regular weave ≈ Planck lattice
     */
    static std::string planck_lattice() {
        return "Regular weave ≈ cubic lattice at scale l_P (polymer)";
    }

    /**
     * @brief Embedding in continuum
     *
     * Weave embedded in smooth manifold Σ
     */
    static std::string embedding() {
        return "Weave γ embedded in topological Σ (diff invariance later)";
    }

    /**
     * @brief Effective continuum
     *
     * For ε << L: effective continuum description
     */
    static std::string effective_continuum() {
        return "ε << L: effective GR with quantum corrections ~ (l_P/L)²";
    }
};

/**
 * @brief Loop Quantum Cosmology (LQC)
 *
 * Application of LQG to cosmology
 */
class LoopQuantumCosmology {
public:
    /**
     * @brief Big Bounce
     *
     * Quantum geometry replaces Big Bang singularity with bounce
     */
    static std::string big_bounce() {
        return "Big Bang singularity → Big Bounce (ρ ≤ ρ_max ~ ρ_Planck)";
    }

    /**
     * @brief Maximum density
     *
     * ρ_max ~ 0.41 ρ_Planck (quantum geometry bound)
     */
    static double max_density_ratio() {
        return 0.41;  // ρ_max/ρ_Planck
    }

    /**
     * @brief Friedmann equation modification
     *
     * H² = (8πG/3)ρ(1 - ρ/ρ_crit) (LQC correction)
     */
    static std::string modified_friedmann() {
        return "H² = (8πG/3)ρ(1 - ρ/ρ_crit) (bounce when ρ = ρ_crit)";
    }

    /**
     * @brief Volume quantization
     *
     * Universe volume quantized: V_n = n V_0
     */
    static std::string volume_quantization() {
        return "V_universe = n × V_Planck (discrete quantum geometry)";
    }

    /**
     * @brief Effective dynamics
     *
     * Quantum corrections to classical cosmology
     */
    static std::string effective_dynamics() {
        return "Effective equation: quantum corrections ∝ ρ/ρ_Planck";
    }

    /**
     * @brief Pre-big-bang
     *
     * Universe existed before bounce (contracting phase)
     */
    static std::string pre_big_bang() {
        return "Contracting universe → bounce → expanding universe";
    }

    /**
     * @brief Observational signatures
     *
     * CMB anomalies, tensor-to-scalar ratio
     */
    static std::string observations() {
        return "CMB: suppressed power at large scales, r < 0.01";
    }

    /**
     * @brief Compute Planck density
     *
     * ρ_Planck = c⁵/(ℏG²) (maximum density scale)
     */
    static double planck_density() {
        return (constants::c * constants::c * constants::c * constants::c * constants::c) /
               (constants::hbar * constants::G * constants::G);
    }

    /**
     * @brief Compute critical density for bounce
     *
     * ρ_crit ≈ 0.41 ρ_Planck
     */
    static double critical_density() {
        return max_density_ratio() * planck_density();
    }

    /**
     * @brief Compute Hubble parameter with LQC corrections
     *
     * H² = (8πG/3)ρ(1 - ρ/ρ_crit)
     */
    static double hubble_parameter_lqc(double density) {
        double rho_crit = critical_density();
        double factor = (8.0 * M_PI * constants::G / 3.0) * density * (1.0 - density / rho_crit);
        return std::sqrt(std::max(0.0, factor));  // Ensure non-negative
    }

    /**
     * @brief Compute scale factor evolution near bounce
     *
     * a(t) for given density (assuming flat FRW)
     */
    static double scale_factor_near_bounce(double density, double density_today = 1.0) {
        // a ∝ ρ^(-1/2) in radiation era near bounce
        return std::sqrt(density_today / std::max(density, 1e-100));
    }
};

/**
 * @brief Inflation in LQC
 *
 * Inflationary cosmology in loop quantum framework
 */
class InflationLQC {
public:
    /**
     * @brief Quantum bounce inflation
     *
     * Bounce can seed inflation naturally
     */
    static std::string bounce_inflation() {
        return "Bounce → high energy density → slow-roll inflation";
    }

    /**
     * @brief Slow-roll conditions
     *
     * Modified in LQC: ε, η slow-roll parameters
     */
    static std::string slow_roll() {
        return "Slow-roll: ε = (1/2)(V'/V)² << 1, η = V''/V << 1";
    }

    /**
     * @brief Power spectrum
     *
     * P(k) = (H²/2π)² (1/2ε) (quantum corrections)
     */
    static std::string power_spectrum() {
        return "P(k) modified: LQC corrections at trans-Planckian scales";
    }

    /**
     * @brief Tensor-to-scalar ratio
     *
     * r = 16ε (observable prediction)
     */
    static double tensor_to_scalar(double epsilon) {
        return 16.0 * epsilon;
    }

    /**
     * @brief Pre-inflationary dynamics
     *
     * Quantum geometry before inflation
     */
    static std::string pre_inflation() {
        return "Pre-inflation: quantum bounce sets initial conditions";
    }

    /**
     * @brief Trans-Planckian problem
     *
     * LQC provides UV completion
     */
    static std::string trans_planckian() {
        return "Trans-Planckian modes: LQC discrete geometry = natural cutoff";
    }

    /**
     * @brief Graceful exit
     *
     * Transition from inflation to radiation era
     */
    static std::string graceful_exit() {
        return "Reheating after inflation (quantum corrections small)";
    }

    /**
     * @brief Compute slow-roll parameter ε
     *
     * ε = (1/2)(V'/V)² where V is potential, V' is derivative
     */
    static double epsilon_slow_roll(double V, double V_prime) {
        if (std::abs(V) < 1e-100) return 1.0;  // Avoid division by zero
        double ratio = V_prime / V;
        return 0.5 * ratio * ratio;
    }

    /**
     * @brief Compute slow-roll parameter η
     *
     * η = V''/V where V'' is second derivative
     */
    static double eta_slow_roll(double V, double V_double_prime) {
        if (std::abs(V) < 1e-100) return 1.0;  // Avoid division by zero
        return V_double_prime / V;
    }

    /**
     * @brief Compute scalar power spectrum amplitude
     *
     * P_s(k) = (H²/2π)² (1/2ε) at horizon crossing
     */
    static double scalar_power_spectrum(double hubble, double epsilon) {
        if (epsilon < 1e-100) return 0.0;  // Avoid division by zero
        double H_over_2pi = hubble / (2.0 * M_PI);
        return (H_over_2pi * H_over_2pi) / (2.0 * epsilon);
    }

    /**
     * @brief Compute tensor power spectrum amplitude
     *
     * P_t(k) = 2(H/2π)²
     */
    static double tensor_power_spectrum(double hubble) {
        double H_over_2pi = hubble / (2.0 * M_PI);
        return 2.0 * H_over_2pi * H_over_2pi;
    }

    /**
     * @brief Compute spectral index
     *
     * n_s - 1 ≈ -6ε + 2η
     */
    static double spectral_index(double epsilon, double eta) {
        return 1.0 - 6.0 * epsilon + 2.0 * eta;
    }

    /**
     * @brief Compute number of e-folds
     *
     * N = ∫ H dt ≈ ∫ (V/V') dφ (slow-roll approximation)
     */
    static double efolds_estimate(double V_initial, double V_final, double V_avg, double V_prime_avg) {
        if (std::abs(V_prime_avg) < 1e-100) return 0.0;
        // Simple estimate: N ≈ (V/V') Δφ
        return (V_avg / V_prime_avg) * std::log(V_initial / V_final);
    }
};

/**
 * @brief Black Hole Thermodynamics - Statistical Ensemble
 *
 * Microcanonical ensemble for black holes
 */
class BHStatisticalEnsemble {
public:
    /**
     * @brief Horizon as boundary
     *
     * Isolated horizon: inner boundary condition
     */
    static std::string isolated_horizon() {
        return "Isolated horizon Δ: (null, non-expanding, weakly isolated)";
    }

    /**
     * @brief Horizon area
     *
     * Area A = 4πr²_s (Schwarzschild radius)
     */
    static double schwarzschild_area(double mass) {
        double r_s = 2.0 * constants::G * mass / (constants::c * constants::c);
        return 4.0 * M_PI * r_s * r_s;
    }

    /**
     * @brief Chern-Simons theory
     *
     * Horizon degrees of freedom from SU(2) CS theory
     */
    static std::string chern_simons() {
        return "Horizon: U(1) CS theory → punctures with spins";
    }

    /**
     * @brief Microstate counting
     *
     * Ω(A) = number of spin network punctures on horizon
     */
    static std::string microstates() {
        return "Microstates: spin network punctures piercing horizon";
    }

    /**
     * @brief Entropy formula
     *
     * S = k_B ln Ω(A)
     */
    static std::string entropy() {
        return "S = k_B ln Ω (Boltzmann entropy from microstate counting)";
    }

    /**
     * @brief Quantum geometry on horizon
     *
     * Horizon = 2D quantum surface
     */
    static std::string quantum_horizon() {
        return "Horizon: 2D spin network with quantized area";
    }
};

/**
 * @brief Derivation of Bekenstein-Hawking Entropy
 *
 * S_BH = A/(4l_P²) from LQG
 */
class BekensteinHawkingEntropy {
public:
    /**
     * @brief Area constraint
     *
     * Horizon area A = Σ_p a_p where a_p = 8πγl_P²√(j_p(j_p+1))
     */
    static std::string area_constraint() {
        return "A_horizon = Σ_p 8πγl_P²√(j_p(j_p+1)) (sum over punctures)";
    }

    /**
     * @brief Number of punctures
     *
     * N ~ A/(8πγl_P²√(j(j+1))) for dominant spin j
     */
    static double number_of_punctures(double area, double gamma_IP, int twice_j) {
        double j = twice_j / 2.0;
        double a_j = 8.0 * M_PI * gamma_IP * constants::l_P * constants::l_P *
                     std::sqrt(j * (j + 1.0));
        return area / a_j;
    }

    /**
     * @brief Counting formula
     *
     * Ω(A,N) = combinatorial count of configurations
     */
    static std::string counting() {
        return "Ω(A,N): count spin assignments {j_p} with Σ a_p = A";
    }

    /**
     * @brief Dominant contribution
     *
     * j = 1/2 punctures dominate for large A
     */
    static std::string dominant_spin() {
        return "Dominant: j = 1/2 (minimal area quanta)";
    }

    /**
     * @brief Entropy calculation
     *
     * S = k_B ln Ω ≈ γ₀ A/(4l_P²) where γ₀ ≈ 0.2375
     */
    static double entropy(double area, double gamma_IP) {
        return constants::k_B * area / (4.0 * gamma_IP * constants::l_P * constants::l_P);
    }

    /**
     * @brief Immirzi parameter fixing
     *
     * γ fixed by matching S = A/(4l_P²)
     */
    static std::string immirzi_fixing() {
        return "γ ≈ ln(2)/(π√3) ≈ 0.2375 from S_BH match";
    }

    /**
     * @brief Bekenstein-Hawking formula
     *
     * S_BH = k_B A/(4l_P²) = k_B c³A/(4ℏG)
     */
    static std::string bh_formula() {
        return "S_BH = k_B c³ A/(4ℏG) (exact agreement with Hawking!)";
    }

    /**
     * @brief Quantum corrections
     *
     * Leading correction: logarithmic in A
     */
    static std::string corrections() {
        return "S = A/(4γl_P²) - (3/2)ln(A/l_P²) + O(1) (quantum corrections)";
    }
};

/**
 * @brief Ringing Modes Frequencies
 *
 * Quasi-normal modes of black holes from LQG
 */
class RingingModes {
public:
    /**
     * @brief Quasi-normal modes
     *
     * Damped oscillations: ω = ω_R - iω_I
     */
    static std::string qnm() {
        return "QNM: h(t) ~ e^(-ω_I t) e^(iω_R t) (ringdown)";
    }

    /**
     * @brief Bohr correspondence
     *
     * ℏω ~ ΔE between horizon quantum states
     */
    static std::string bohr() {
        return "ℏω_R ~ ΔA/A (area quantum transitions)";
    }

    /**
     * @brief Area spectrum
     *
     * A_n = 8πγl_P² Σ_{i=1}^n √(j_i(j_i+1))
     */
    static std::string area_spectrum() {
        return "Discrete area: ΔA_min = 8πγl_P²√(j(j+1))";
    }

    /**
     * @brief Frequency formula
     *
     * ω ~ c/(r_s) × (area quantum)
     */
    static double frequency_estimate(double mass) {
        double r_s = 2.0 * constants::G * mass / (constants::c * constants::c);
        return constants::c / r_s;
    }

    /**
     * @brief Damping time
     *
     * τ = 1/ω_I ~ r_s/c (light crossing time)
     */
    static std::string damping() {
        return "τ_damp ~ r_s/c (horizon crossing time)";
    }

    /**
     * @brief Overtones
     *
     * Multiple frequencies from different j transitions
     */
    static std::string overtones() {
        return "Overtone spectrum: ω_n from different Δj transitions";
    }

    /**
     * @brief Observable via gravitational waves
     *
     * LIGO/Virgo detection of ringdown
     */
    static std::string observability() {
        return "GW ringdown: test LQG via QNM spectrum";
    }

    /**
     * @brief Compute damping time
     *
     * τ = r_s/c (light crossing time of horizon)
     */
    static double damping_time(double mass) {
        double r_s = 2.0 * constants::G * mass / (constants::c * constants::c);
        return r_s / constants::c;
    }

    /**
     * @brief Compute QNM quality factor
     *
     * Q = ω_R/(2ω_I) = π(ω_R τ)
     */
    static double quality_factor(double mass) {
        double omega_R = frequency_estimate(mass);
        double tau = damping_time(mass);
        return M_PI * omega_R * tau;
    }
};

/**
 * @brief Bekenstein-Mukhanov Effect
 *
 * Discrete area spectrum → discrete entropy
 */
class BekensteinMukhanovEffect {
public:
    /**
     * @brief Area quantization
     *
     * Area comes in discrete quanta
     */
    static std::string quantization() {
        return "ΔA = 8πγl_P²√(j(j+1)) (area quantum)";
    }

    /**
     * @brief Entropy spacing
     *
     * ΔS ~ k_B (entropy increases by discrete steps)
     */
    static std::string entropy_spacing() {
        return "ΔS ~ k_B ΔA/(4γl_P²) ~ k_B (discrete entropy!)";
    }

    /**
     * @brief Evaporation discrete
     *
     * Hawking radiation emitted in quanta
     */
    static std::string discrete_evaporation() {
        return "BH evaporation: discrete jumps ΔA (not continuous)";
    }

    /**
     * @brief Mukhanov proposal
     *
     * Entropy eigenvalue spacing constant
     */
    static std::string mukhanov() {
        return "Mukhanov: ΔS = k_B ln(n) for some integer n";
    }

    /**
     * @brief LQG realization
     *
     * Adding one j=1/2 puncture
     */
    static std::string lqg_realization() {
        return "LQG: ΔS from adding j=1/2 puncture to horizon";
    }

    /**
     * @brief Observational prospects
     *
     * Final stages of BH evaporation
     */
    static std::string observations() {
        return "Observable: Planck-mass BH evaporation (discrete bursts)";
    }

    /**
     * @brief Compute area quantum for given spin
     *
     * ΔA_j = 8πγl_P²√(j(j+1))
     */
    static double area_quantum(int twice_j, double gamma_IP = constants::gamma) {
        double j = twice_j / 2.0;
        return 8.0 * M_PI * gamma_IP * constants::l_P * constants::l_P *
               std::sqrt(j * (j + 1.0));
    }

    /**
     * @brief Compute entropy quantum
     *
     * ΔS = k_B ΔA/(4γl_P²)
     */
    static double entropy_quantum(int twice_j, double gamma_IP = constants::gamma) {
        double dA = area_quantum(twice_j, gamma_IP);
        return constants::k_B * dA / (4.0 * gamma_IP * constants::l_P * constants::l_P);
    }

    /**
     * @brief Compute minimal area change
     *
     * Minimal j=1/2 puncture
     */
    static double minimal_area_change(double gamma_IP = constants::gamma) {
        return area_quantum(1, gamma_IP);  // j=1/2, so twice_j=1
    }

    /**
     * @brief Compute minimal entropy change
     *
     * ΔS_min for j=1/2
     */
    static double minimal_entropy_change(double gamma_IP = constants::gamma) {
        return entropy_quantum(1, gamma_IP);  // j=1/2
    }
};

/**
 * @brief Observable Effects
 *
 * Testable predictions of LQG
 */
class ObservableEffects {
public:
    /**
     * @brief Modified dispersion relations
     *
     * E² = p²c² + m²c⁴ + ξ(l_P/λ)E³ (Planck-scale corrections)
     */
    static std::string dispersion() {
        return "Modified: E² ≈ p²c² + α(l_P/λ)E³ (Lorentz violation at l_P)";
    }

    /**
     * @brief Time-of-flight delays
     *
     * Different energies arrive at different times
     */
    static std::string time_delay() {
        return "Δt ~ ΔE × l_P/c × D (D = distance to source)";
    }

    /**
     * @brief Gamma-ray bursts
     *
     * High-energy photons from GRBs test quantum gravity
     */
    static std::string grb() {
        return "GRB: E ~ 10 GeV, D ~ Gpc → Δt ~ μs (testable!)";
    }

    /**
     * @brief CMB anomalies
     *
     * LQC predicts large-scale power suppression
     */
    static std::string cmb() {
        return "CMB: suppressed power at l < 30 (LQC bounce signature)";
    }

    /**
     * @brief Black hole shadows
     *
     * Quantum corrections to photon sphere
     */
    static std::string bh_shadow() {
        return "BH shadow: quantum corrections Δr/r ~ (l_P/r_s)²";
    }

    /**
     * @brief Gravitational wave echoes
     *
     * Quantum structure near horizon
     */
    static std::string gw_echoes() {
        return "GW echoes: reflections from quantum horizon structure";
    }

    /**
     * @brief Primordial gravitational waves
     *
     * r < 0.01 (tensor-to-scalar ratio from LQC)
     */
    static std::string primordial_gw() {
        return "Primordial GW: r < 0.01 (LQC bounce → low tensor modes)";
    }

    /**
     * @brief Current constraints
     *
     * No Lorentz violation detected to 10^-17 eV
     */
    static std::string constraints() {
        return "Current limits: Lorentz violation ξ < 10⁻² (Fermi-LAT)";
    }

    /**
     * @brief Compute time-of-flight delay
     *
     * Δt ≈ α × (ΔE/E_Planck) × (l_P/c) × D
     * where α is suppression parameter, ΔE is energy difference, D is distance
     */
    static double time_delay(double energy_diff_eV, double distance_m, double alpha = 1.0) {
        double E_Planck_eV = constants::E_P / 1.602176634e-19;  // Convert J to eV
        double factor = alpha * (energy_diff_eV / E_Planck_eV) * (constants::l_P / constants::c);
        return factor * distance_m;
    }

    /**
     * @brief Compute modified energy-momentum relation
     *
     * E² ≈ p²c² + m²c⁴ + α(l_P E)²E (first-order LQG correction)
     */
    static double modified_energy(double momentum, double mass, double alpha = 1.0) {
        double E_classical_sq = momentum * momentum * constants::c * constants::c +
                               mass * mass * constants::c * constants::c * constants::c * constants::c;
        double E_classical = std::sqrt(E_classical_sq);
        // Perturbative correction
        double correction = alpha * constants::l_P * E_classical * E_classical;
        return E_classical + correction;
    }

    /**
     * @brief Compute velocity modification
     *
     * v/c = ∂E/∂p ≈ pc²/E + (LQG corrections)
     */
    static double modified_velocity(double momentum, double energy) {
        if (energy < 1e-100) return 0.0;
        return (momentum * constants::c * constants::c) / energy;
    }

    /**
     * @brief Compute photon sphere correction
     *
     * Δr/r_s ~ (l_P/r_s)²
     */
    static double photon_sphere_correction(double mass) {
        double r_s = 2.0 * constants::G * mass / (constants::c * constants::c);
        double ratio = constants::l_P / r_s;
        return ratio * ratio;
    }

    /**
     * @brief Compute GW echo time delay
     *
     * Δt_echo ~ r_s/c × ln(r_s/l_P) (quantum horizon effect)
     */
    static double gw_echo_delay(double mass) {
        double r_s = 2.0 * constants::G * mass / (constants::c * constants::c);
        double ln_factor = std::log(r_s / constants::l_P);
        return (r_s / constants::c) * ln_factor;
    }

    /**
     * @brief Estimate GRB time delay for typical parameters
     *
     * E_high ~ 10 GeV, E_low ~ 1 GeV, D ~ 1 Gpc
     */
    static double grb_delay_estimate() {
        double E_high = 10.0e9 * 1.602176634e-19;  // 10 GeV in Joules
        double E_low = 1.0e9 * 1.602176634e-19;    // 1 GeV in Joules
        double distance = 1.0e9 * 3.086e22;         // 1 Gpc in meters
        double dE = E_high - E_low;
        return time_delay(dE / 1.602176634e-19, distance, 1.0);  // Convert back to eV
    }
};

/**
 * @brief From Loops to Spinfoams
 *
 * Transition from canonical to covariant formulation
 */
class FromLoopsToSpinfoams {
public:
    /**
     * @brief Canonical vs covariant
     *
     * LQG (3+1) vs spinfoams (4D covariant)
     */
    static std::string canonical_vs_covariant() {
        return "Canonical LQG: evolution in time; Spinfoams: spacetime histories";
    }

    /**
     * @brief Path integral
     *
     * Z = Σ_geometries e^(iS/ℏ) (sum over 4-geometries)
     */
    static std::string path_integral() {
        return "Spinfoam: Z = Σ_σ A(σ) (sum over 2-complexes σ)";
    }

    /**
     * @brief Spacetime foam
     *
     * 4D analogue of Feynman path integral
     */
    static std::string foam() {
        return "Spinfoam = quantum 4-geometry (spacetime histories)";
    }

    /**
     * @brief Spin networks as boundaries
     *
     * Spin networks = boundary states of spinfoams
     */
    static std::string boundaries() {
        return "∂(spinfoam) = spin network (3D boundary of 4D geometry)";
    }

    /**
     * @brief Amplitude
     *
     * A(σ) = transition amplitude between spin network states
     */
    static std::string amplitude() {
        return "A(σ): ⟨s_f|e^(-iĤt)|s_i⟩ (boundary spin networks s_i, s_f)";
    }

    /**
     * @brief Wheeler-DeWitt
     *
     * Ĥ|Ψ⟩ = 0 → spinfoam sum
     */
    static std::string wheeler_dewitt() {
        return "Hamiltonian constraint → sum over spacetime histories";
    }

    /**
     * @brief Discrete vs continuum
     *
     * Spinfoam = discrete approximation to path integral
     */
    static std::string discrete_continuum() {
        return "Spinfoam: discretized path integral (Regge-like)";
    }
};

/**
 * @brief Spinfoam Formalism
 *
 * Mathematical structure of spinfoams
 */
class SpinfoamFormalism {
public:
    /**
     * @brief 2-complex
     *
     * Σ = (faces F, edges E, vertices V) dual to triangulation
     */
    static std::string two_complex() {
        return "2-complex σ: vertices V, edges E, faces F (dual to triangulation)";
    }

    /**
     * @brief Labeling
     *
     * Faces → spins j_f, edges → intertwiners i_e
     */
    static std::string labeling() {
        return "Faces: spins j_f, Edges: intertwiners i_e, Vertices: amplitudes";
    }

    /**
     * @brief Amplitude formula
     *
     * A(σ) = Σ_{j,i} Π_f (2j_f+1) Π_v A_v({j,i})
     */
    static std::string amplitude() {
        return "A(σ) = Σ_{j,i} Π_f d_j Π_v A_v (vertex amplitude A_v)";
    }

    /**
     * @brief Vertex amplitude
     *
     * A_v from 15j symbol (10 faces, 15 edges around vertex)
     */
    static std::string vertex_amplitude() {
        return "A_v = {15j symbol} (Wigner symbol for 4-simplex)";
    }

    /**
     * @brief Face amplitude
     *
     * (2j+1) = dimension of SU(2) irrep
     */
    static std::string face_amplitude() {
        return "Face: d_j = 2j+1 (dimension weight)";
    }

    /**
     * @brief Partition function
     *
     * Z = Σ_σ A(σ) (sum over 2-complexes)
     */
    static std::string partition_function() {
        return "Z = Σ_σ A(σ) (sum over all labelings)";
    }

    /**
     * @brief Transition amplitude
     *
     * ⟨s_f|s_i⟩ = Σ_σ A(σ) with ∂σ = s_i ∪ s_f
     */
    static std::string transition() {
        return "⟨s_f|s_i⟩ = Σ_{σ:∂σ=s_i∪s_f} A(σ)";
    }
};

/**
 * @brief Spinfoam Boundaries
 *
 * Boundary conditions for spinfoams
 */
class SpinfoamBoundaries {
public:
    /**
     * @brief Boundary spin networks
     *
     * ∂σ consists of spin networks on initial/final slices
     */
    static std::string boundary_sn() {
        return "∂σ = s_initial ∪ s_final (3D spin networks)";
    }

    /**
     * @brief Gluing
     *
     * Spinfoams glued along common boundaries
     */
    static std::string gluing() {
        return "Gluing: σ₁ ∪_{s} σ₂ (compose along shared boundary s)";
    }

    /**
     * @brief Cylindrical consistency
     *
     * ⟨s|s⟩ = 1 (normalization)
     */
    static std::string cylindrical() {
        return "Cylindrical: ⟨s|s⟩ = 1 (probability conservation)";
    }

    /**
     * @brief No boundary
     *
     * Closed spinfoam: ∂σ = ∅ (universe with no boundary)
     */
    static std::string no_boundary() {
        return "∂σ = ∅: closed universe (Hartle-Hawking state)";
    }

    /**
     * @brief Time evolution
     *
     * Hamiltonian from boundary deformation
     */
    static std::string evolution() {
        return "Ĥ generates boundary deformation (discrete time step)";
    }
};

/**
 * @brief 3D Quantum Gravity
 *
 * Toy model for spinfoams
 */
class ThreeDQuantumGravity {
public:
    /**
     * @brief Topological theory
     *
     * 3D gravity has no local degrees of freedom
     */
    static std::string topological() {
        return "3D gravity: topological (no local DOF, only global)";
    }

    /**
     * @brief Ponzano-Regge model
     *
     * 3D spinfoam: 6j symbols at tetrahedra
     */
    static std::string ponzano_regge() {
        return "Z = Σ_{j} Π_tetrahedra {6j symbols} (Ponzano-Regge)";
    }

    /**
     * @brief Turaev-Viro model
     *
     * Quantum group version (q = root of unity)
     */
    static std::string turaev_viro() {
        return "Turaev-Viro: quantum 6j at q^k = 1 (finite sum)";
    }

    /**
     * @brief Exactly solvable
     *
     * 3D gravity exactly solvable (pedagogical model)
     */
    static std::string solvable() {
        return "Exactly solvable: partition function computable";
    }

    /**
     * @brief BTZ black hole
     *
     * 3D black hole solution
     */
    static std::string btz() {
        return "BTZ: 3D rotating BH, horizon entropy from Chern-Simons";
    }
};

/**
 * @brief BF Theory
 *
 * Topological BF theory as starting point
 */
class BFTheory {
public:
    /**
     * @brief BF action
     *
     * S = ∫ Tr(B ∧ F) where F = dA + A ∧ A
     */
    static std::string action() {
        return "S_BF = ∫ Tr(B ∧ F) (topological, gauge invariant)";
    }

    /**
     * @brief Plebanski formulation
     *
     * GR = constrained BF theory
     */
    static std::string plebanski() {
        return "GR: BF + simplicity constraints (B ~ e ∧ e)";
    }

    /**
     * @brief Simplicity constraints
     *
     * B^IJ = *(e^I ∧ e^J) (Plebanski constraint)
     */
    static std::string simplicity() {
        return "Simplicity: B^IJ ~ ε^IJKL e_K ∧ e_L (relates B to metric)";
    }

    /**
     * @brief Quantization
     *
     * BF theory quantizes to TQFT
     */
    static std::string quantization() {
        return "Quantum BF: exactly solvable TQFT (state sum)";
    }

    /**
     * @brief Spinfoam from BF
     *
     * Impose simplicity → gravity spinfoam
     */
    static std::string to_spinfoam() {
        return "BF + simplicity (quantum) → gravity spinfoam models";
    }
};

/**
 * @brief Spinfoam/GFT Duality
 *
 * Group field theory formulation
 */
class SpinfoamGFTDuality {
public:
    /**
     * @brief Group field theory
     *
     * Field φ on G^×4 (four copies of SU(2))
     */
    static std::string gft() {
        return "GFT: field φ(g₁,g₂,g₃,g₄) on SU(2)^×4";
    }

    /**
     * @brief Feynman diagrams = spinfoams
     *
     * GFT Feynman graphs are spinfoams
     */
    static std::string duality() {
        return "GFT Feynman diagrams ↔ spinfoams (dual description)";
    }

    /**
     * @brief Action
     *
     * S_GFT = kinetic + interaction terms
     */
    static std::string action() {
        return "S = ∫ φ̄ K φ + λ ∫ φ⁵ + ... (field theory action)";
    }

    /**
     * @brief Vertices = simplices
     *
     * GFT interaction vertex = 4-simplex
     */
    static std::string vertices() {
        return "5-valent vertex in GFT = 4-simplex in spinfoam";
    }

    /**
     * @brief Condensates
     *
     * GFT condensate → macroscopic geometry
     */
    static std::string condensate() {
        return "⟨φ⟩ ≠ 0: condensate → emergent continuum spacetime";
    }

    /**
     * @brief Cosmology from GFT
     *
     * GFT condensate cosmology
     */
    static std::string cosmology() {
        return "GFT condensate → FRW cosmology (emergent)";
    }
};

/**
 * @brief BC (Barrett-Crane) Models
 *
 * Early spinfoam models
 */
class BCModels {
public:
    /**
     * @brief Barrett-Crane model
     *
     * Euclidean 4D spinfoam with SO(4) gauge group
     */
    static std::string barrett_crane() {
        return "BC: vertex = 10j symbol from SO(4) recoupling";
    }

    /**
     * @brief Euclidean signature
     *
     * SO(4) = SU(2) × SU(2)
     */
    static std::string euclidean() {
        return "SO(4): j_f = (j_+, j_-) (self-dual/anti-self-dual)";
    }

    /**
     * @brief Simplicity constraints
     *
     * j_+ = j_- (simple representation)
     */
    static std::string simplicity() {
        return "Simplicity: j_+ = j_- (Barrett-Crane constraint)";
    }

    /**
     * @brief Problems
     *
     * BC model has issues (wrong classical limit)
     */
    static std::string problems() {
        return "Problems: no propagating DOF, wrong n-point functions";
    }

    /**
     * @brief Superseded
     *
     * Replaced by EPRL/FK models
     */
    static std::string superseded() {
        return "BC superseded by EPRL/FK models (correct semiclassical limit)";
    }
};

/**
 * @brief Group Field Theory
 *
 * Detailed GFT formalism
 */
class GroupFieldTheory {
public:
    /**
     * @brief Field definition
     *
     * φ: G^×n → ℂ (field on group manifold)
     */
    static std::string field() {
        return "φ(g₁,...,g_n): field on SU(2)^×n (n = 3 or 4)";
    }

    /**
     * @brief Gauge invariance
     *
     * φ(g₁h,...,g_nh) = φ(g₁,...,g_n) for h ∈ G
     */
    static std::string gauge_invariance() {
        return "Right invariance: φ(g_i h) = φ(g_i)";
    }

    /**
     * @brief Kinetic term
     *
     * ∫ φ̄(g) K(g,g') φ(g') (Laplacian on G)
     */
    static std::string kinetic() {
        return "Kinetic: ∫ φ̄ (Δ_G + m²) φ (Laplacian on group)";
    }

    /**
     * @brief Interaction
     *
     * Vertex with n strands (n = d+1 in d dimensions)
     */
    static std::string interaction() {
        return "Interaction: ∫ φ^{d+1} (d = spacetime dimension)";
    }

    /**
     * @brief Propagator
     *
     * ⟨φ(g)φ̄(g')⟩ (Green's function on group)
     */
    static std::string propagator() {
        return "Propagator: ⟨φφ̄⟩ = Σ_j d_j χ_j(gg'^{-1}) (Peter-Weyl)";
    }

    /**
     * @brief Renormalization
     *
     * GFT renormalization program (ongoing)
     */
    static std::string renormalization() {
        return "Renormalization: scale-dependent couplings (active research)";
    }

    /**
     * @brief Mean field
     *
     * Saddle point approximation
     */
    static std::string mean_field() {
        return "Mean field: δS/δφ = 0 → classical field equation";
    }
};

/**
 * @brief Lorentzian Models
 *
 * Lorentzian signature spinfoams (physical spacetime)
 */
class LorentzianModels {
public:
    /**
     * @brief EPRL model
     *
     * Engle-Pereira-Rovelli-Livine model (Lorentzian)
     */
    static std::string eprl() {
        return "EPRL: SL(2,C) spinfoam with Immirzi parameter γ";
    }

    /**
     * @brief FK model
     *
     * Freidel-Krasnov model (alternative Lorentzian)
     */
    static std::string fk() {
        return "FK: closely related to EPRL (same semiclassical limit)";
    }

    /**
     * @brief Gauge group
     *
     * SL(2,C) for Lorentzian signature
     */
    static std::string gauge_group() {
        return "SL(2,C): double cover of SO↑(1,3) (Lorentz group)";
    }

    /**
     * @brief Representations
     *
     * (ρ, k) labels for SL(2,C) irreps
     */
    static std::string representations() {
        return "SL(2,C) reps: (ρ,k) where ρ ∈ ℝ⁺, k ∈ ℤ/2";
    }

    /**
     * @brief Vertex amplitude
     *
     * A_v from SL(2,C) 15j symbol (complex)
     */
    static std::string vertex() {
        return "A_v: SL(2,C) {15j} symbol (Lorentzian geometry)";
    }

    /**
     * @brief Semiclassical limit
     *
     * Large spins → Regge geometry
     */
    static std::string semiclassical() {
        return "j → ∞: EPRL/FK → Regge action (correct limit!)";
    }

    /**
     * @brief Asymptotics
     *
     * Stationary phase approximation
     */
    static std::string asymptotics() {
        return "A_v ~ e^(iS_Regge/ℏ) for large j (WKB)";
    }

    /**
     * @brief n-point functions
     *
     * Graviton propagator from boundary correlators
     */
    static std::string n_point() {
        return "⟨h_μν(x) h_ρσ(y)⟩: graviton 2-point function";
    }
};

/**
 * @brief Physics from Spinfoams
 *
 * Physical observables and predictions
 */
class PhysicsFromSpinfoams {
public:
    /**
     * @brief Graviton propagator
     *
     * 2-point function of metric perturbations
     */
    static std::string graviton_propagator() {
        return "⟨h(x)h(y)⟩ ~ 1/|x-y|² (from boundary correlators)";
    }

    /**
     * @brief Particle scattering
     *
     * Matter coupled to quantum geometry
     */
    static std::string scattering() {
        return "S-matrix: ⟨out|in⟩ from spinfoam with matter insertions";
    }

    /**
     * @brief Minkowski vacuum
     *
     * Flat space as spinfoam state
     */
    static std::string minkowski() {
        return "η_μν: sum over flat spinfoams (coherent state)";
    }

    /**
     * @brief Coherent states
     *
     * Semiclassical geometries as coherent superpositions
     */
    static std::string coherent() {
        return "|g_μν⟩ ~ Σ_σ e^(-||σ-g||²) |σ⟩ (peaked on classical g)";
    }

    /**
     * @brief Quantum corrections
     *
     * Deviations from GR at Planck scale
     */
    static std::string corrections() {
        return "⟨O⟩ = ⟨O⟩_GR + ℏ ⟨O⟩_(1) + ℏ² ⟨O⟩_(2) + ...";
    }

    /**
     * @brief Cosmological constant
     *
     * Λ from spinfoam asymptotics?
     */
    static std::string cosmological_constant() {
        return "Λ_eff from large-scale spinfoam structure (speculative)";
    }

    /**
     * @brief Emergence of locality
     *
     * How does local spacetime emerge?
     */
    static std::string locality() {
        return "Locality emerges from fine-grained spinfoam (coarse graining)";
    }

    /**
     * @brief Continuum limit
     *
     * Refinement limit of triangulation
     */
    static std::string continuum() {
        return "Continuum: ε → 0 limit of spinfoam (triangulation refined)";
    }
};

} // namespace loop_quantum_gravity
} // namespace physics

#endif // PHYSICS_LOOP_QUANTUM_GRAVITY_HPP

