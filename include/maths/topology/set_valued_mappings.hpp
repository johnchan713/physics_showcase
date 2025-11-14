#ifndef MATHS_TOPOLOGY_SET_VALUED_MAPPINGS_HPP
#define MATHS_TOPOLOGY_SET_VALUED_MAPPINGS_HPP

#include <vector>
#include <set>
#include <functional>
#include <algorithm>
#include <string>

/**
 * @file set_valued_mappings.hpp
 * @brief Set-valued mappings (multimaps/correspondences) and their continuity
 *
 * Implements:
 * - Set-valued mappings F: X ⇉ Y
 * - Upper and lower hemicontinuity
 * - Continuity of correspondences
 * - Kakutani Fixed-Point Theorem
 * - Michael Selection Theorem
 */

namespace maths::topology {

/**
 * @class SetValuedMappings
 * @brief Theory of multimaps/correspondences
 */
class SetValuedMappings {
public:
    /**
     * @brief Basic definitions
     */
    static std::string definitions() {
        return "Set-Valued Mappings (Multimaps/Correspondences):\n"
               "\n"
               "A set-valued mapping F: X ⇉ Y assigns to each x ∈ X\n"
               "a subset F(x) ⊆ Y (possibly empty).\n"
               "\n"
               "Domain: dom(F) = {x ∈ X : F(x) ≠ ∅}\n"
               "Graph: graph(F) = {(x, y) : y ∈ F(x)} ⊆ X × Y\n"
               "Image: F(A) = ⋃_{x ∈ A} F(x)\n"
               "\n"
               "Inverse/Preimage:\n"
               "- F⁻¹(B) = {x : F(x) ∩ B ≠ ∅} (small preimage)\n"
               "- F⁺(B) = {x : F(x) ⊆ B} (large preimage)\n"
               "\n"
               "Note: F⁺(B) ⊆ F⁻¹(B) always.\n"
               "\n"
               "Examples:\n"
               "- Subdifferential ∂f: ℝⁿ ⇉ ℝⁿ\n"
               "- Budget correspondence in economics\n"
               "- Solution maps for equations";
    }

    /**
     * @brief Continuity concepts
     */
    static std::string continuityConcepts() {
        return "Continuity of Set-Valued Mappings:\n"
               "\n"
               "F: X ⇉ Y is:\n"
               "\n"
               "1. Upper hemicontinuous (uhc) at x₀ if:\n"
               "   For every open V ⊇ F(x₀), ∃ neighborhood U of x₀:\n"
               "   F(x) ⊆ V for all x ∈ U\n"
               "   \n"
               "   Equivalently: F⁺(V) is open for every open V\n"
               "\n"
               "2. Lower hemicontinuous (lhc) at x₀ if:\n"
               "   For every open V with V ∩ F(x₀) ≠ ∅,\n"
               "   ∃ neighborhood U of x₀: V ∩ F(x) ≠ ∅ for all x ∈ U\n"
               "   \n"
               "   Equivalently: F⁻¹(V) is open for every open V\n"
               "\n"
               "3. Continuous if both uhc and lhc\n"
               "\n"
               "Note: For single-valued maps (F(x) = {f(x)}),\n"
               "      uhc + lhc = continuity of f";
    }

    /**
     * @brief Sequential characterization
     */
    static std::string sequentialCharacterization() {
        return "Sequential Characterization:\n"
               "\n"
               "For metric spaces:\n"
               "\n"
               "F is uhc at x₀ ⟺\n"
               "  ∀xₙ → x₀, ∀yₙ ∈ F(xₙ) convergent subsequence,\n"
               "  lim yₙₖ ∈ F(x₀)\n"
               "\n"
               "  (if F has compact values: every sequence yₙ ∈ F(xₙ)\n"
               "   has subsequence converging to some y ∈ F(x₀))\n"
               "\n"
               "F is lhc at x₀ ⟺\n"
               "  ∀xₙ → x₀, ∀y ∈ F(x₀), ∃yₙ ∈ F(xₙ): yₙ → y\n"
               "\n"
               "  (can approximate limits from F(x₀) by elements\n"
               "   from nearby F(xₙ))\n"
               "\n"
               "Intuition:\n"
               "- uhc: graph doesn't \"blow up\" at x₀\n"
               "- lhc: no \"sudden disappearance\" of values";
    }

    /**
     * @brief Closed graph theorem for multimaps
     */
    static std::string closedGraph() {
        return "Closed Graph Property:\n"
               "\n"
               "F has closed graph if:\n"
               "  xₙ → x, yₙ → y, yₙ ∈ F(xₙ) ⇒ y ∈ F(x)\n"
               "\n"
               "Theorem: If Y is compact and F has compact values,\n"
               "then F is uhc ⟺ F has closed graph\n"
               "\n"
               "Warning: Closed graph ⇏ uhc in general!\n"
               "         Need compactness assumption.\n"
               "\n"
               "Maximum Theorem: If F: X ⇉ Y is uhc with compact values,\n"
               "f: X × Y → ℝ continuous, then:\n"
               "  M(x) = max_{y ∈ F(x)} f(x, y)\n"
               "is continuous and argmax correspondence\n"
               "  G(x) = {y ∈ F(x) : f(x, y) = M(x)}\n"
               "is uhc with compact values.";
    }
};

/**
 * @class FixedPointTheorems
 * @brief Fixed-point theorems for set-valued maps
 */
class FixedPointTheorems {
public:
    /**
     * @brief Kakutani Fixed-Point Theorem
     */
    static std::string kakutani() {
        return "Kakutani Fixed-Point Theorem:\n"
               "\n"
               "Let K ⊆ ℝⁿ be non-empty, compact, convex.\n"
               "Let F: K ⇉ K satisfy:\n"
               "  1. F(x) is non-empty, convex for all x\n"
               "  2. F is upper hemicontinuous\n"
               "  3. F has closed graph\n"
               "\n"
               "Then F has a fixed point: ∃x* ∈ K with x* ∈ F(x*).\n"
               "\n"
               "Generalizes Brouwer Fixed-Point Theorem\n"
               "(Brouwer = Kakutani for single-valued maps)\n"
               "\n"
               "Applications:\n"
               "- Nash equilibrium existence (game theory)\n"
               "- Existence of Walrasian equilibrium (economics)\n"
               "- Variational inequalities\n"
               "- Differential inclusions\n"
               "\n"
               "Proof sketch: Approximate by continuous single-valued\n"
               "selections, apply Brouwer, take limit.";
    }

    /**
     * @brief Michael Selection Theorem
     */
    static std::string michaelSelection() {
        return "Michael Selection Theorem:\n"
               "\n"
               "Let X be paracompact, Y Banach space.\n"
               "Let F: X ⇉ Y be lhc with closed convex values.\n"
               "\n"
               "Then F admits a continuous selection:\n"
               "  ∃f: X → Y continuous with f(x) ∈ F(x) for all x\n"
               "\n"
               "Key points:\n"
               "- Needs lower hemicontinuity (not upper!)\n"
               "- Requires convexity of values\n"
               "- Paracompactness technical (holds for metric spaces)\n"
               "\n"
               "Applications:\n"
               "- Constructing continuous solutions to inclusions\n"
               "- Approximation theory\n"
               "- Continuous choice functions\n"
               "\n"
               "Non-convex case: Generally no continuous selection\n"
               "(e.g., F(x) = {-1, 1} on ℝ is lhc but has no\n"
               "continuous selection)";
    }

    /**
     * @brief Other fixed-point results
     */
    static std::string otherResults() {
        return "Other Fixed-Point Theorems:\n"
               "\n"
               "Banach Fixed-Point (Contraction Mapping):\n"
               "- F: X → X contraction on complete metric space\n"
               "- Has unique fixed point\n"
               "- Iterates converge to fixed point\n"
               "\n"
               "Brouwer Fixed-Point:\n"
               "- f: K → K continuous, K ⊆ ℝⁿ compact convex\n"
               "- Has fixed point\n"
               "\n"
               "Schauder Fixed-Point:\n"
               "- f: K → K continuous, K ⊆ Banach space\n"
               "- K compact convex\n"
               "- Has fixed point\n"
               "\n"
               "Fan-Glicksberg:\n"
               "- F: K ⇉ K uhc with compact convex values\n"
               "- K locally convex, compact\n"
               "- Has fixed point\n"
               "\n"
               "All proved using tools from algebraic topology\n"
               "(degree theory, simplicial approximation, etc.)";
    }
};

/**
 * @class Measurability
 * @brief Measurable multifunctions
 */
class Measurability {
public:
    /**
     * @brief Measurable selections
     */
    static std::string measurableSelections() {
        return "Measurable Set-Valued Mappings:\n"
               "\n"
               "F: (X, Σ) ⇉ Y is measurable if:\n"
               "  F⁻¹(U) = {x : F(x) ∩ U ≠ ∅} ∈ Σ\n"
               "for every open U ⊆ Y.\n"
               "\n"
               "Selection: measurable f: X → Y with f(x) ∈ F(x)\n"
               "\n"
               "Kuratowski-Ryll-Nardzewski Theorem:\n"
               "If (Y, d) is separable metric and F: X ⇉ Y measurable\n"
               "with closed values, then F admits measurable selection.\n"
               "\n"
               "Aumann Measurable Selection:\n"
               "- Weaker than KRN\n"
               "- For integrable selections\n"
               "\n"
               "Applications:\n"
               "- Stochastic processes\n"
               "- Random sets\n"
               "- Optimal control with uncertainty";
    }
};

} // namespace maths::topology

#endif // MATHS_TOPOLOGY_SET_VALUED_MAPPINGS_HPP
