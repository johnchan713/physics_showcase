#ifndef MATHS_LINEAR_ALGEBRA_VECTORS_HPP
#define MATHS_LINEAR_ALGEBRA_VECTORS_HPP

#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <sstream>

/**
 * @file vectors.hpp
 * @brief Vector operations and properties
 *
 * Implements:
 * - Vector arithmetic (addition, scalar multiplication)
 * - Dot product and cross product
 * - Vector norms and normalization
 * - Orthogonality and projections
 * - Linear independence
 */

namespace maths::linear_algebra {

/**
 * @class Vector
 * @brief N-dimensional vector with operations
 */
class Vector {
private:
    std::vector<double> components;

public:
    /**
     * @brief Construct vector from components
     */
    explicit Vector(const std::vector<double>& comps) : components(comps) {}

    /**
     * @brief Construct zero vector of dimension n
     */
    explicit Vector(size_t n) : components(n, 0.0) {}

    /**
     * @brief Get dimension
     */
    size_t dimension() const {
        return components.size();
    }

    /**
     * @brief Access component (mutable)
     */
    double& operator[](size_t i) {
        if (i >= components.size()) {
            throw std::out_of_range("Vector index out of range");
        }
        return components[i];
    }

    /**
     * @brief Access component (const)
     */
    const double& operator[](size_t i) const {
        if (i >= components.size()) {
            throw std::out_of_range("Vector index out of range");
        }
        return components[i];
    }

    /**
     * @brief Get all components
     */
    const std::vector<double>& getComponents() const {
        return components;
    }

    /**
     * @brief Vector addition
     */
    Vector operator+(const Vector& other) const {
        if (dimension() != other.dimension()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }

        std::vector<double> result(dimension());
        for (size_t i = 0; i < dimension(); ++i) {
            result[i] = components[i] + other[i];
        }
        return Vector(result);
    }

    /**
     * @brief Vector subtraction
     */
    Vector operator-(const Vector& other) const {
        if (dimension() != other.dimension()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }

        std::vector<double> result(dimension());
        for (size_t i = 0; i < dimension(); ++i) {
            result[i] = components[i] - other[i];
        }
        return Vector(result);
    }

    /**
     * @brief Scalar multiplication
     */
    Vector operator*(double scalar) const {
        std::vector<double> result(dimension());
        for (size_t i = 0; i < dimension(); ++i) {
            result[i] = components[i] * scalar;
        }
        return Vector(result);
    }

    /**
     * @brief Negation
     */
    Vector operator-() const {
        return (*this) * (-1.0);
    }

    /**
     * @brief Dot product (inner product)
     *
     * v · w = Σ vᵢwᵢ
     */
    double dot(const Vector& other) const {
        if (dimension() != other.dimension()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }

        double result = 0.0;
        for (size_t i = 0; i < dimension(); ++i) {
            result += components[i] * other[i];
        }
        return result;
    }

    /**
     * @brief Euclidean norm (magnitude, length)
     *
     * ‖v‖ = √(v₁² + v₂² + ... + vₙ²)
     */
    double norm() const {
        return std::sqrt(dot(*this));
    }

    /**
     * @brief Normalize vector (unit vector in same direction)
     *
     * v̂ = v / ‖v‖
     */
    Vector normalize() const {
        double n = norm();
        if (n == 0.0) {
            throw std::runtime_error("Cannot normalize zero vector");
        }
        return (*this) * (1.0 / n);
    }

    /**
     * @brief Distance between two vectors
     *
     * d(v, w) = ‖v - w‖
     */
    double distance(const Vector& other) const {
        return (*this - other).norm();
    }

    /**
     * @brief Check if zero vector
     */
    bool isZero(double tolerance = 1e-10) const {
        return norm() < tolerance;
    }

    /**
     * @brief String representation
     */
    std::string toString() const {
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < dimension(); ++i) {
            oss << components[i];
            if (i < dimension() - 1) oss << ", ";
        }
        oss << "]";
        return oss.str();
    }
};

/**
 * @brief Scalar multiplication (left side)
 */
inline Vector operator*(double scalar, const Vector& v) {
    return v * scalar;
}

/**
 * @class CrossProduct
 * @brief Cross product for 3D vectors
 */
class CrossProduct {
public:
    /**
     * @brief Cross product v × w (3D only)
     *
     * v × w = (v₂w₃ - v₃w₂, v₃w₁ - v₁w₃, v₁w₂ - v₂w₁)
     *
     * Result is perpendicular to both v and w
     */
    static Vector compute(const Vector& v, const Vector& w) {
        if (v.dimension() != 3 || w.dimension() != 3) {
            throw std::invalid_argument("Cross product requires 3D vectors");
        }

        std::vector<double> result(3);
        result[0] = v[1] * w[2] - v[2] * w[1];
        result[1] = v[2] * w[0] - v[0] * w[2];
        result[2] = v[0] * w[1] - v[1] * w[0];

        return Vector(result);
    }

    /**
     * @brief Properties of cross product
     */
    static std::string properties() {
        return "Cross Product Properties (v × w):\n"
               "\n"
               "1. Anticommutative: v × w = -(w × v)\n"
               "2. Distributive: v × (w + u) = v × w + v × u\n"
               "3. Scalar mult: (cv) × w = c(v × w) = v × (cw)\n"
               "4. Perpendicular: (v × w) · v = 0, (v × w) · w = 0\n"
               "5. Magnitude: ‖v × w‖ = ‖v‖‖w‖ sin(θ)\n"
               "6. Parallel: v × w = 0 ⟺ v ∥ w\n"
               "7. Right-hand rule determines direction";
    }

    /**
     * @brief Geometric interpretation
     */
    static std::string geometricMeaning() {
        return "Geometric interpretation:\n"
               "\n"
               "‖v × w‖ = area of parallelogram spanned by v and w\n"
               "\n"
               "Direction: perpendicular to plane containing v and w\n"
               "           (right-hand rule)\n"
               "\n"
               "If v and w are parallel: v × w = 0";
    }
};

/**
 * @class VectorProjection
 * @brief Projections and orthogonal components
 */
class VectorProjection {
public:
    /**
     * @brief Scalar projection of v onto w
     *
     * comp_w(v) = (v · w) / ‖w‖
     */
    static double scalarProjection(const Vector& v, const Vector& w) {
        double w_norm = w.norm();
        if (w_norm == 0.0) {
            throw std::runtime_error("Cannot project onto zero vector");
        }
        return v.dot(w) / w_norm;
    }

    /**
     * @brief Vector projection of v onto w
     *
     * proj_w(v) = [(v · w) / (w · w)] w
     */
    static Vector vectorProjection(const Vector& v, const Vector& w) {
        double w_dot_w = w.dot(w);
        if (w_dot_w == 0.0) {
            throw std::runtime_error("Cannot project onto zero vector");
        }
        double scalar = v.dot(w) / w_dot_w;
        return w * scalar;
    }

    /**
     * @brief Orthogonal component (rejection)
     *
     * perp_w(v) = v - proj_w(v)
     */
    static Vector orthogonalComponent(const Vector& v, const Vector& w) {
        return v - vectorProjection(v, w);
    }

    /**
     * @brief Verify v = proj_w(v) + perp_w(v)
     */
    static bool verifyDecomposition(const Vector& v, const Vector& w,
                                    double tolerance = 1e-10) {
        Vector proj = vectorProjection(v, w);
        Vector perp = orthogonalComponent(v, w);
        Vector sum = proj + perp;

        return (v - sum).norm() < tolerance;
    }

    /**
     * @brief Properties
     */
    static std::string properties() {
        return "Projection Properties:\n"
               "\n"
               "1. v = proj_w(v) + perp_w(v) (orthogonal decomposition)\n"
               "2. proj_w(v) · perp_w(v) = 0 (orthogonal components)\n"
               "3. proj_w(v) is parallel to w\n"
               "4. perp_w(v) is perpendicular to w\n"
               "5. ‖v‖² = ‖proj_w(v)‖² + ‖perp_w(v)‖² (Pythagoras)";
    }
};

/**
 * @class Orthogonality
 * @brief Orthogonal and orthonormal vectors
 */
class Orthogonality {
public:
    /**
     * @brief Check if two vectors are orthogonal
     *
     * v ⊥ w ⟺ v · w = 0
     */
    static bool areOrthogonal(const Vector& v, const Vector& w,
                             double tolerance = 1e-10) {
        return std::abs(v.dot(w)) < tolerance;
    }

    /**
     * @brief Check if vector is unit vector
     *
     * ‖v‖ = 1
     */
    static bool isUnitVector(const Vector& v, double tolerance = 1e-10) {
        return std::abs(v.norm() - 1.0) < tolerance;
    }

    /**
     * @brief Check if vectors are orthonormal
     *
     * Orthonormal: pairwise orthogonal and all unit vectors
     */
    static bool areOrthonormal(const std::vector<Vector>& vectors,
                              double tolerance = 1e-10) {
        for (size_t i = 0; i < vectors.size(); ++i) {
            // Check unit length
            if (!isUnitVector(vectors[i], tolerance)) {
                return false;
            }

            // Check orthogonality with all other vectors
            for (size_t j = i + 1; j < vectors.size(); ++j) {
                if (!areOrthogonal(vectors[i], vectors[j], tolerance)) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Gram-Schmidt orthogonalization
     *
     * Convert linearly independent vectors to orthogonal vectors
     */
    static std::vector<Vector> gramSchmidt(const std::vector<Vector>& vectors) {
        if (vectors.empty()) return {};

        std::vector<Vector> orthogonal;
        orthogonal.reserve(vectors.size());

        // First vector unchanged
        orthogonal.push_back(vectors[0]);

        // Orthogonalize subsequent vectors
        for (size_t i = 1; i < vectors.size(); ++i) {
            Vector v = vectors[i];

            // Subtract projections onto all previous orthogonal vectors
            for (size_t j = 0; j < orthogonal.size(); ++j) {
                v = v - VectorProjection::vectorProjection(v, orthogonal[j]);
            }

            // Check if resulting vector is nonzero
            if (!v.isZero()) {
                orthogonal.push_back(v);
            }
        }

        return orthogonal;
    }

    /**
     * @brief Gram-Schmidt with normalization (orthonormal basis)
     */
    static std::vector<Vector> gramSchmidtOrthonormal(const std::vector<Vector>& vectors) {
        auto orthogonal = gramSchmidt(vectors);

        std::vector<Vector> orthonormal;
        orthonormal.reserve(orthogonal.size());

        for (const auto& v : orthogonal) {
            orthonormal.push_back(v.normalize());
        }

        return orthonormal;
    }

    /**
     * @brief Standard basis vectors
     *
     * e₁ = (1, 0, 0, ..., 0)
     * e₂ = (0, 1, 0, ..., 0)
     * etc.
     */
    static std::vector<Vector> standardBasis(size_t dimension) {
        std::vector<Vector> basis;
        basis.reserve(dimension);

        for (size_t i = 0; i < dimension; ++i) {
            std::vector<double> components(dimension, 0.0);
            components[i] = 1.0;
            basis.push_back(Vector(components));
        }

        return basis;
    }
};

/**
 * @class LinearIndependence
 * @brief Check linear independence of vectors
 */
class LinearIndependence {
public:
    /**
     * @brief Check if vectors are linearly independent
     *
     * Vectors v₁, ..., vₖ are linearly independent if
     * c₁v₁ + ... + cₖvₖ = 0 implies c₁ = ... = cₖ = 0
     */
    static bool areIndependent(const std::vector<Vector>& vectors,
                               double tolerance = 1e-10) {
        if (vectors.empty()) return true;

        // Use Gram-Schmidt: if any vector becomes zero, they're dependent
        auto orthogonal = Orthogonality::gramSchmidt(vectors);

        // If we lost vectors, they were linearly dependent
        return orthogonal.size() == vectors.size();
    }

    /**
     * @brief Find dimension of span (rank)
     *
     * Dimension of space spanned by vectors
     */
    static size_t spanDimension(const std::vector<Vector>& vectors) {
        auto orthogonal = Orthogonality::gramSchmidt(vectors);
        return orthogonal.size();
    }

    /**
     * @brief Extract linearly independent subset
     */
    static std::vector<Vector> extractIndependentSet(const std::vector<Vector>& vectors,
                                                     double tolerance = 1e-10) {
        return Orthogonality::gramSchmidt(vectors);
    }

    /**
     * @brief Properties
     */
    static std::string properties() {
        return "Linear Independence:\n"
               "\n"
               "Vectors v₁, ..., vₖ are linearly independent if:\n"
               "c₁v₁ + ... + cₖvₖ = 0 ⟹ c₁ = ... = cₖ = 0\n"
               "\n"
               "Equivalently:\n"
               "- No vector is a linear combination of others\n"
               "- Gram-Schmidt produces k nonzero vectors\n"
               "- Matrix with columns v₁, ..., vₖ has rank k\n"
               "\n"
               "Maximum number of independent vectors in ℝⁿ: n";
    }
};

/**
 * @class VectorSpaces
 * @brief Vector space properties and operations
 */
class VectorSpaces {
public:
    /**
     * @brief Axioms of vector space
     */
    static std::string axioms() {
        return "Vector Space Axioms:\n"
               "\n"
               "Addition properties:\n"
               "1. Closure: u + v ∈ V\n"
               "2. Commutativity: u + v = v + u\n"
               "3. Associativity: (u + v) + w = u + (v + w)\n"
               "4. Zero vector: ∃0 such that v + 0 = v\n"
               "5. Additive inverse: ∃(-v) such that v + (-v) = 0\n"
               "\n"
               "Scalar multiplication properties:\n"
               "6. Closure: cv ∈ V\n"
               "7. Associativity: c(dv) = (cd)v\n"
               "8. Identity: 1v = v\n"
               "9. Distributivity: c(u + v) = cu + cv\n"
               "10. Distributivity: (c + d)v = cv + dv";
    }

    /**
     * @brief Subspace definition
     */
    static std::string subspaceDefinition() {
        return "Subspace:\n"
               "\n"
               "W ⊆ V is a subspace if:\n"
               "1. 0 ∈ W (contains zero vector)\n"
               "2. u, v ∈ W ⟹ u + v ∈ W (closed under addition)\n"
               "3. v ∈ W, c ∈ ℝ ⟹ cv ∈ W (closed under scalar mult)\n"
               "\n"
               "Examples:\n"
               "- Span of vectors\n"
               "- Kernel (null space) of linear transformation\n"
               "- Range (column space) of matrix";
    }

    /**
     * @brief Basis definition
     */
    static std::string basisDefinition() {
        return "Basis:\n"
               "\n"
               "Set {v₁, ..., vₙ} is a basis for V if:\n"
               "1. Linearly independent\n"
               "2. Span(v₁, ..., vₙ) = V\n"
               "\n"
               "Properties:\n"
               "- Every v ∈ V has unique representation v = c₁v₁ + ... + cₙvₙ\n"
               "- All bases of V have same number of vectors (dimension)\n"
               "- Standard basis for ℝⁿ: {e₁, ..., eₙ}";
    }

    /**
     * @brief Dimension theorem
     */
    static std::string dimensionTheorem() {
        return "Dimension Theorems:\n"
               "\n"
               "1. All bases of V have same size = dim(V)\n"
               "\n"
               "2. If dim(V) = n:\n"
               "   - n independent vectors form a basis\n"
               "   - n spanning vectors form a basis\n"
               "   - > n vectors must be dependent\n"
               "   - < n vectors cannot span V\n"
               "\n"
               "3. Rank-Nullity Theorem:\n"
               "   dim(domain) = dim(range) + dim(kernel)";
    }
};

} // namespace maths::linear_algebra

#endif // MATHS_LINEAR_ALGEBRA_VECTORS_HPP
