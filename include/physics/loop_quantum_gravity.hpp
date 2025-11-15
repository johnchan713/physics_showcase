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

} // namespace loop_quantum_gravity
} // namespace physics

#endif // PHYSICS_LOOP_QUANTUM_GRAVITY_HPP
