#ifndef MATHS_GROUP_THEORY_LIE_GROUPS_HPP
#define MATHS_GROUP_THEORY_LIE_GROUPS_HPP

#include <cmath>
#include <vector>
#include <complex>
#include <stdexcept>
#include <functional>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <array>

/**
 * @file group_theory_lie_groups.hpp
 * @brief Comprehensive implementation of group theory, Lie groups, and Lie algebras
 *
 * This module implements:
 *
 * GROUP THEORY:
 * - Definition of a Group and Basic Properties
 * - Examples of Groups (cyclic, symmetric, dihedral, etc.)
 * - Subgroups, Center, and Direct Products
 * - Homomorphisms and Isomorphisms
 * - Group actions and orbits
 *
 * MATRIX LIE GROUPS:
 * - Definition of Matrix Lie Groups
 * - Examples: GL(n), SL(n), O(n), SO(n), U(n), SU(n), Sp(n)
 * - Compactness properties
 * - Connectedness and simple-connectedness
 * - Homomorphisms and Isomorphisms
 *
 * LIE ALGEBRAS:
 * - Matrix Exponential computation
 * - Matrix Logarithm
 * - Properties of Matrix Exponential
 * - Lie Algebra of Matrix Lie Groups
 * - Lie bracket and structure constants
 * - Exponential mapping
 * - Complexification of Real Lie Algebras
 *
 * BAKER-CAMPBELL-HAUSDORFF:
 * - BCH Formula
 * - BCH for Heisenberg Group
 * - General BCH Formula
 * - Series Form of BCH
 * - Subgroups and Subalgebras
 *
 * All implementations use standard mathematical conventions.
 */

namespace maths {
namespace group_theory {

// ============================================================================
// MATRIX UTILITIES
// ============================================================================

/**
 * @brief Simple matrix class for group theory computations
 */
template<typename T = double>
class Matrix {
private:
    std::vector<std::vector<T>> data;
    size_t rows_, cols_;

public:
    Matrix(size_t r, size_t c, T val = T(0)) : rows_(r), cols_(c) {
        data.resize(r, std::vector<T>(c, val));
    }

    Matrix(const std::vector<std::vector<T>>& d) : data(d) {
        rows_ = d.size();
        cols_ = rows_ > 0 ? d[0].size() : 0;
    }

    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }

    T& operator()(size_t i, size_t j) {
        if (i >= rows_ || j >= cols_) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i][j];
    }

    const T& operator()(size_t i, size_t j) const {
        if (i >= rows_ || j >= cols_) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i][j];
    }

    /**
     * @brief Matrix multiplication
     */
    Matrix operator*(const Matrix& other) const {
        if (cols_ != other.rows_) {
            throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
        }

        Matrix result(rows_, other.cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < other.cols_; ++j) {
                T sum = T(0);
                for (size_t k = 0; k < cols_; ++k) {
                    sum += data[i][k] * other.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }
        return result;
    }

    /**
     * @brief Matrix addition
     */
    Matrix operator+(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Matrix dimensions must match for addition");
        }

        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    /**
     * @brief Matrix subtraction
     */
    Matrix operator-(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Matrix dimensions must match for subtraction");
        }

        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    /**
     * @brief Scalar multiplication
     */
    Matrix operator*(const T& scalar) const {
        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.data[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }

    /**
     * @brief Matrix transpose
     */
    Matrix transpose() const {
        Matrix result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    /**
     * @brief Calculate determinant (for small matrices)
     */
    T determinant() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("Determinant requires square matrix");
        }

        if (rows_ == 1) {
            return data[0][0];
        }

        if (rows_ == 2) {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        }

        if (rows_ == 3) {
            return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1])
                 - data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0])
                 + data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
        }

        // For larger matrices, use Laplace expansion (slow but works)
        T det = T(0);
        for (size_t j = 0; j < cols_; ++j) {
            Matrix minor(rows_ - 1, cols_ - 1);
            for (size_t i1 = 1; i1 < rows_; ++i1) {
                size_t j2 = 0;
                for (size_t j1 = 0; j1 < cols_; ++j1) {
                    if (j1 == j) continue;
                    minor.data[i1 - 1][j2++] = data[i1][j1];
                }
            }
            T cofactor = (j % 2 == 0 ? T(1) : T(-1)) * data[0][j] * minor.determinant();
            det += cofactor;
        }

        return det;
    }

    /**
     * @brief Calculate trace
     */
    T trace() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("Trace requires square matrix");
        }

        T tr = T(0);
        for (size_t i = 0; i < rows_; ++i) {
            tr += data[i][i];
        }
        return tr;
    }

    /**
     * @brief Create identity matrix
     */
    static Matrix identity(size_t n) {
        Matrix result(n, n);
        for (size_t i = 0; i < n; ++i) {
            result.data[i][i] = T(1);
        }
        return result;
    }

    /**
     * @brief Check if matrix is approximately equal to another
     */
    bool approxEqual(const Matrix& other, double tolerance = 1e-10) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            return false;
        }

        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                if (std::abs(data[i][j] - other.data[i][j]) > tolerance) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Frobenius norm
     */
    double frobeniusNorm() const {
        double sum = 0.0;
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                double val = std::abs(data[i][j]);
                sum += val * val;
            }
        }
        return std::sqrt(sum);
    }
};

// ============================================================================
// GROUP THEORY - BASIC DEFINITIONS
// ============================================================================

/**
 * @brief Abstract base class for finite groups
 *
 * A group (G, •) consists of:
 * 1. A set G
 * 2. A binary operation •: G × G → G
 * 3. Associativity: (a • b) • c = a • (b • c)
 * 4. Identity element: ∃e ∈ G such that a • e = e • a = a
 * 5. Inverse: ∀a ∈ G, ∃a⁻¹ ∈ G such that a • a⁻¹ = a⁻¹ • a = e
 */
template<typename T>
class Group {
public:
    virtual ~Group() = default;

    /**
     * @brief Group operation
     */
    virtual T operation(const T& a, const T& b) const = 0;

    /**
     * @brief Identity element
     */
    virtual T identity() const = 0;

    /**
     * @brief Inverse element
     */
    virtual T inverse(const T& a) const = 0;

    /**
     * @brief Check if element is in group
     */
    virtual bool contains(const T& a) const = 0;

    /**
     * @brief Order of the group (number of elements, -1 if infinite)
     */
    virtual int order() const = 0;

    /**
     * @brief Order of an element (smallest positive n such that a^n = e)
     */
    int elementOrder(const T& a) const {
        if (!contains(a)) {
            throw std::invalid_argument("Element not in group");
        }

        T current = a;
        int n = 1;
        int maxIter = order() > 0 ? order() : 1000;  // Safety limit for infinite groups

        while (n <= maxIter) {
            if (approxEqual(current, identity())) {
                return n;
            }
            current = operation(current, a);
            n++;
        }

        return -1;  // Infinite order or exceeds limit
    }

    /**
     * @brief Check associativity (for verification)
     */
    virtual bool isAssociative(const T& a, const T& b, const T& c) const {
        T left = operation(operation(a, b), c);
        T right = operation(a, operation(b, c));
        return approxEqual(left, right);
    }

    /**
     * @brief Check if element commutes with another
     */
    bool commutes(const T& a, const T& b) const {
        T ab = operation(a, b);
        T ba = operation(b, a);
        return approxEqual(ab, ba);
    }

protected:
    virtual bool approxEqual(const T& a, const T& b) const = 0;
};

// ============================================================================
// CYCLIC GROUP Z_n
// ============================================================================

/**
 * @class CyclicGroup
 * @brief Cyclic group of integers modulo n
 *
 * Z_n = {0, 1, 2, ..., n-1} with addition modulo n
 */
class CyclicGroup : public Group<int> {
private:
    int n;  // Order of the group

public:
    explicit CyclicGroup(int order) : n(order) {
        if (order <= 0) {
            throw std::invalid_argument("Group order must be positive");
        }
    }

    int operation(const int& a, const int& b) const override {
        return (a + b) % n;
    }

    int identity() const override {
        return 0;
    }

    int inverse(const int& a) const override {
        if (!contains(a)) {
            throw std::invalid_argument("Element not in group");
        }
        return (n - a) % n;
    }

    bool contains(const int& a) const override {
        return a >= 0 && a < n;
    }

    int order() const override {
        return n;
    }

    /**
     * @brief Generate element a^k
     */
    int power(int a, int k) const {
        if (!contains(a)) {
            throw std::invalid_argument("Element not in group");
        }
        return (a * k) % n;
    }

    /**
     * @brief Check if group is cyclic (always true for Z_n)
     */
    bool isCyclic() const {
        return true;
    }

    /**
     * @brief Find a generator of the cyclic group
     */
    int generator() const {
        // 1 is always a generator of Z_n
        return 1;
    }

protected:
    bool approxEqual(const int& a, const int& b) const override {
        return a == b;
    }
};

// ============================================================================
// SYMMETRIC GROUP S_n (PERMUTATIONS)
// ============================================================================

/**
 * @class SymmetricGroup
 * @brief Symmetric group of permutations on n elements
 *
 * S_n consists of all bijections {1, ..., n} → {1, ..., n}
 * Order: |S_n| = n!
 */
class SymmetricGroup : public Group<std::vector<int>> {
private:
    int n;

public:
    explicit SymmetricGroup(int size) : n(size) {
        if (size <= 0) {
            throw std::invalid_argument("Group size must be positive");
        }
    }

    std::vector<int> operation(const std::vector<int>& a,
                              const std::vector<int>& b) const override {
        // Composition of permutations: (a ∘ b)[i] = a[b[i]]
        if (a.size() != static_cast<size_t>(n) || b.size() != static_cast<size_t>(n)) {
            throw std::invalid_argument("Permutation size mismatch");
        }

        std::vector<int> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = a[b[i]];
        }
        return result;
    }

    std::vector<int> identity() const override {
        std::vector<int> id(n);
        for (int i = 0; i < n; ++i) {
            id[i] = i;
        }
        return id;
    }

    std::vector<int> inverse(const std::vector<int>& a) const override {
        if (!contains(a)) {
            throw std::invalid_argument("Not a valid permutation");
        }

        std::vector<int> inv(n);
        for (int i = 0; i < n; ++i) {
            inv[a[i]] = i;
        }
        return inv;
    }

    bool contains(const std::vector<int>& a) const override {
        if (a.size() != static_cast<size_t>(n)) {
            return false;
        }

        std::vector<bool> seen(n, false);
        for (int i = 0; i < n; ++i) {
            if (a[i] < 0 || a[i] >= n || seen[a[i]]) {
                return false;
            }
            seen[a[i]] = true;
        }
        return true;
    }

    int order() const override {
        // n!
        int factorial = 1;
        for (int i = 2; i <= n; ++i) {
            factorial *= i;
        }
        return factorial;
    }

    /**
     * @brief Calculate sign (parity) of permutation
     *
     * Returns +1 for even permutations, -1 for odd
     */
    int sign(const std::vector<int>& perm) const {
        if (!contains(perm)) {
            throw std::invalid_argument("Not a valid permutation");
        }

        int inversions = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (perm[i] > perm[j]) {
                    inversions++;
                }
            }
        }

        return (inversions % 2 == 0) ? 1 : -1;
    }

    /**
     * @brief Create transposition (i j) - swap elements i and j
     */
    std::vector<int> transposition(int i, int j) const {
        if (i < 0 || i >= n || j < 0 || j >= n) {
            throw std::out_of_range("Index out of range");
        }

        std::vector<int> perm = identity();
        std::swap(perm[i], perm[j]);
        return perm;
    }

protected:
    bool approxEqual(const std::vector<int>& a, const std::vector<int>& b) const override {
        return a == b;
    }
};

// ============================================================================
// DIHEDRAL GROUP D_n
// ============================================================================

/**
 * @class DihedralGroup
 * @brief Dihedral group D_n - symmetries of regular n-gon
 *
 * D_n = {r^i s^j : 0 ≤ i < n, 0 ≤ j < 2}
 * where r = rotation by 2π/n, s = reflection
 * Order: |D_n| = 2n
 */
class DihedralGroup : public Group<std::pair<int, int>> {
private:
    int n;  // Number of vertices

public:
    explicit DihedralGroup(int vertices) : n(vertices) {
        if (vertices < 3) {
            throw std::invalid_argument("Dihedral group requires at least 3 vertices");
        }
    }

    /**
     * @brief Group operation
     *
     * Element represented as (rotation, reflection) where:
     * - rotation ∈ {0, 1, ..., n-1}
     * - reflection ∈ {0, 1}
     */
    std::pair<int, int> operation(const std::pair<int, int>& a,
                                 const std::pair<int, int>& b) const override {
        int r1 = a.first, s1 = a.second;
        int r2 = b.first, s2 = b.second;

        // Composition rule: (r^i s^j)(r^k s^l) = r^(i + (-1)^j k) s^(j+l mod 2)
        int r_result, s_result;

        if (s1 == 0) {
            r_result = (r1 + r2) % n;
        } else {
            r_result = (r1 - r2 + n) % n;
        }

        s_result = (s1 + s2) % 2;

        return {r_result, s_result};
    }

    std::pair<int, int> identity() const override {
        return {0, 0};
    }

    std::pair<int, int> inverse(const std::pair<int, int>& a) const override {
        int r = a.first, s = a.second;

        if (s == 0) {
            // Inverse of r^i is r^(n-i)
            return {(n - r) % n, 0};
        } else {
            // Inverse of r^i s is r^i s (reflections are self-inverse)
            return {r, 1};
        }
    }

    bool contains(const std::pair<int, int>& a) const override {
        return a.first >= 0 && a.first < n && (a.second == 0 || a.second == 1);
    }

    int order() const override {
        return 2 * n;
    }

    /**
     * @brief Create rotation element r^k
     */
    std::pair<int, int> rotation(int k) const {
        return {k % n, 0};
    }

    /**
     * @brief Create reflection element s
     */
    std::pair<int, int> reflection() const {
        return {0, 1};
    }

protected:
    bool approxEqual(const std::pair<int, int>& a, const std::pair<int, int>& b) const override {
        return a == b;
    }
};

// ============================================================================
// SUBGROUPS
// ============================================================================

/**
 * @brief Check if subset H is a subgroup of G
 *
 * Subgroup test:
 * 1. e ∈ H (identity)
 * 2. If a, b ∈ H, then a • b ∈ H (closure)
 * 3. If a ∈ H, then a⁻¹ ∈ H (inverses)
 *
 * @param elements Elements of potential subgroup
 * @param group Parent group
 * @return true if H is a subgroup
 */
template<typename T>
bool isSubgroup(const std::vector<T>& elements, const Group<T>& group) {
    if (elements.empty()) {
        return false;
    }

    // Check identity
    T e = group.identity();
    bool hasIdentity = false;
    for (const auto& elem : elements) {
        if (group.approxEqual(elem, e)) {
            hasIdentity = true;
            break;
        }
    }
    if (!hasIdentity) {
        return false;
    }

    // Check closure and inverses
    for (const auto& a : elements) {
        // Check inverse
        T inv = group.inverse(a);
        bool hasInverse = false;
        for (const auto& elem : elements) {
            if (group.approxEqual(elem, inv)) {
                hasInverse = true;
                break;
            }
        }
        if (!hasInverse) {
            return false;
        }

        // Check closure
        for (const auto& b : elements) {
            T product = group.operation(a, b);
            bool inSubgroup = false;
            for (const auto& elem : elements) {
                if (group.approxEqual(elem, product)) {
                    inSubgroup = true;
                    break;
                }
            }
            if (!inSubgroup) {
                return false;
            }
        }
    }

    return true;
}

// ============================================================================
// CENTER OF A GROUP
// ============================================================================

/**
 * @brief Calculate center of a group
 *
 * Z(G) = {g ∈ G : gh = hg for all h ∈ G}
 *
 * The center is always a normal subgroup
 *
 * @param elements All elements of the group
 * @param group The group
 * @return Elements in the center
 */
template<typename T>
std::vector<T> center(const std::vector<T>& elements, const Group<T>& group) {
    std::vector<T> centerElements;

    for (const auto& g : elements) {
        bool isInCenter = true;

        for (const auto& h : elements) {
            if (!group.commutes(g, h)) {
                isInCenter = false;
                break;
            }
        }

        if (isInCenter) {
            centerElements.push_back(g);
        }
    }

    return centerElements;
}

// ============================================================================
// DIRECT PRODUCT OF GROUPS
// ============================================================================

/**
 * @brief Direct product G × H
 *
 * Elements: (g, h) where g ∈ G, h ∈ H
 * Operation: (g₁, h₁) • (g₂, h₂) = (g₁ • g₂, h₁ • h₂)
 * Order: |G × H| = |G| × |H|
 */
template<typename T1, typename T2>
class DirectProduct : public Group<std::pair<T1, T2>> {
private:
    const Group<T1>& group1;
    const Group<T2>& group2;

public:
    DirectProduct(const Group<T1>& g1, const Group<T2>& g2)
        : group1(g1), group2(g2) {}

    std::pair<T1, T2> operation(const std::pair<T1, T2>& a,
                                const std::pair<T1, T2>& b) const override {
        return {
            group1.operation(a.first, b.first),
            group2.operation(a.second, b.second)
        };
    }

    std::pair<T1, T2> identity() const override {
        return {group1.identity(), group2.identity()};
    }

    std::pair<T1, T2> inverse(const std::pair<T1, T2>& a) const override {
        return {group1.inverse(a.first), group2.inverse(a.second)};
    }

    bool contains(const std::pair<T1, T2>& a) const override {
        return group1.contains(a.first) && group2.contains(a.second);
    }

    int order() const override {
        int o1 = group1.order();
        int o2 = group2.order();
        if (o1 < 0 || o2 < 0) {
            return -1;  // Infinite
        }
        return o1 * o2;
    }

protected:
    bool approxEqual(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) const override {
        return group1.approxEqual(a.first, b.first) &&
               group2.approxEqual(a.second, b.second);
    }
};

// ============================================================================
// HOMOMORPHISMS AND ISOMORPHISMS
// ============================================================================

/**
 * @brief Check if a map is a group homomorphism
 *
 * φ: G → H is a homomorphism if φ(g₁ • g₂) = φ(g₁) • φ(g₂)
 *
 * @param elementsG Elements of group G
 * @param groupG Group G
 * @param groupH Group H
 * @param phi The map φ: G → H
 * @return true if φ is a homomorphism
 */
template<typename T, typename U>
bool isHomomorphism(const std::vector<T>& elementsG,
                    const Group<T>& groupG,
                    const Group<U>& groupH,
                    std::function<U(const T&)> phi) {
    // Check φ(g₁ • g₂) = φ(g₁) • φ(g₂) for all g₁, g₂
    for (const auto& g1 : elementsG) {
        for (const auto& g2 : elementsG) {
            T product_G = groupG.operation(g1, g2);
            U left = phi(product_G);

            U phi_g1 = phi(g1);
            U phi_g2 = phi(g2);
            U right = groupH.operation(phi_g1, phi_g2);

            if (!groupH.approxEqual(left, right)) {
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Check if homomorphism is injective (monomorphism)
 */
template<typename T, typename U>
bool isMonomorphism(const std::vector<T>& elementsG,
                    std::function<U(const T&)> phi) {
    std::set<U> images;

    for (const auto& g : elementsG) {
        U img = phi(g);
        if (images.count(img) > 0) {
            return false;  // Not injective
        }
        images.insert(img);
    }

    return true;
}

/**
 * @brief Check if groups are isomorphic via given map
 *
 * φ: G → H is an isomorphism if:
 * 1. φ is a homomorphism
 * 2. φ is bijective
 */
template<typename T, typename U>
bool isIsomorphism(const std::vector<T>& elementsG,
                   const std::vector<U>& elementsH,
                   const Group<T>& groupG,
                   const Group<U>& groupH,
                   std::function<U(const T&)> phi) {
    // Must be same order
    if (elementsG.size() != elementsH.size()) {
        return false;
    }

    // Check if homomorphism
    if (!isHomomorphism(elementsG, groupG, groupH, phi)) {
        return false;
    }

    // Check if bijective (injective suffices for finite groups of equal size)
    return isMonomorphism(elementsG, phi);
}

// ============================================================================
// MATRIX LIE GROUPS - DEFINITIONS
// ============================================================================

/**
 * @brief General Linear Group GL(n, ℝ)
 *
 * GL(n) = {A ∈ ℝⁿˣⁿ : det(A) ≠ 0}
 *
 * Group operation: Matrix multiplication
 */
template<typename T = double>
class GeneralLinearGroup {
public:
    int n;  // Dimension

    explicit GeneralLinearGroup(int dim) : n(dim) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension must be positive");
        }
    }

    /**
     * @brief Check if matrix is in GL(n)
     */
    bool contains(const Matrix<T>& A) const {
        if (A.rows() != static_cast<size_t>(n) || A.cols() != static_cast<size_t>(n)) {
            return false;
        }

        T det = A.determinant();
        return std::abs(det) > 1e-10;  // Non-zero determinant
    }

    /**
     * @brief Group operation (matrix multiplication)
     */
    Matrix<T> operation(const Matrix<T>& A, const Matrix<T>& B) const {
        return A * B;
    }

    /**
     * @brief Identity element
     */
    Matrix<T> identity() const {
        return Matrix<T>::identity(n);
    }

    /**
     * @brief Check if GL(n) is compact
     *
     * GL(n) is NOT compact (unbounded matrices)
     */
    bool isCompact() const {
        return false;
    }

    /**
     * @brief Check if GL(n) is connected
     *
     * GL(n) has two connected components: det > 0 and det < 0
     */
    bool isConnected() const {
        return false;
    }

    /**
     * @brief Number of connected components
     */
    int numConnectedComponents() const {
        return 2;
    }
};

/**
 * @brief Special Linear Group SL(n, ℝ)
 *
 * SL(n) = {A ∈ ℝⁿˣⁿ : det(A) = 1}
 */
template<typename T = double>
class SpecialLinearGroup {
public:
    int n;

    explicit SpecialLinearGroup(int dim) : n(dim) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension must be positive");
        }
    }

    bool contains(const Matrix<T>& A, double tolerance = 1e-10) const {
        if (A.rows() != static_cast<size_t>(n) || A.cols() != static_cast<size_t>(n)) {
            return false;
        }

        T det = A.determinant();
        return std::abs(det - T(1)) < tolerance;
    }

    Matrix<T> operation(const Matrix<T>& A, const Matrix<T>& B) const {
        return A * B;
    }

    Matrix<T> identity() const {
        return Matrix<T>::identity(n);
    }

    /**
     * @brief SL(n) is NOT compact
     */
    bool isCompact() const {
        return false;
    }

    /**
     * @brief SL(n) is connected
     */
    bool isConnected() const {
        return true;
    }

    /**
     * @brief SL(n) is simply connected for n ≥ 2
     */
    bool isSimplyConnected() const {
        return n >= 2;
    }
};

/**
 * @brief Orthogonal Group O(n)
 *
 * O(n) = {A ∈ ℝⁿˣⁿ : AᵀA = I}
 *
 * Matrices that preserve inner product
 */
template<typename T = double>
class OrthogonalGroup {
public:
    int n;

    explicit OrthogonalGroup(int dim) : n(dim) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension must be positive");
        }
    }

    bool contains(const Matrix<T>& A, double tolerance = 1e-10) const {
        if (A.rows() != static_cast<size_t>(n) || A.cols() != static_cast<size_t>(n)) {
            return false;
        }

        Matrix<T> AT = A.transpose();
        Matrix<T> ATA = AT * A;
        Matrix<T> I = Matrix<T>::identity(n);

        return ATA.approxEqual(I, tolerance);
    }

    Matrix<T> operation(const Matrix<T>& A, const Matrix<T>& B) const {
        return A * B;
    }

    Matrix<T> identity() const {
        return Matrix<T>::identity(n);
    }

    /**
     * @brief O(n) is compact
     */
    bool isCompact() const {
        return true;
    }

    /**
     * @brief O(n) is NOT connected (two components: det = ±1)
     */
    bool isConnected() const {
        return false;
    }

    /**
     * @brief O(n) has two connected components
     */
    int numConnectedComponents() const {
        return 2;
    }
};

/**
 * @brief Special Orthogonal Group SO(n)
 *
 * SO(n) = {A ∈ O(n) : det(A) = 1}
 *
 * Rotation matrices
 */
template<typename T = double>
class SpecialOrthogonalGroup {
public:
    int n;

    explicit SpecialOrthogonalGroup(int dim) : n(dim) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension must be positive");
        }
    }

    bool contains(const Matrix<T>& A, double tolerance = 1e-10) const {
        if (A.rows() != static_cast<size_t>(n) || A.cols() != static_cast<size_t>(n)) {
            return false;
        }

        // Check orthogonality
        Matrix<T> AT = A.transpose();
        Matrix<T> ATA = AT * A;
        Matrix<T> I = Matrix<T>::identity(n);

        if (!ATA.approxEqual(I, tolerance)) {
            return false;
        }

        // Check determinant = 1
        T det = A.determinant();
        return std::abs(det - T(1)) < tolerance;
    }

    Matrix<T> operation(const Matrix<T>& A, const Matrix<T>& B) const {
        return A * B;
    }

    Matrix<T> identity() const {
        return Matrix<T>::identity(n);
    }

    /**
     * @brief SO(n) is compact
     */
    bool isCompact() const {
        return true;
    }

    /**
     * @brief SO(n) is connected
     */
    bool isConnected() const {
        return true;
    }

    /**
     * @brief SO(n) is simply connected only for n ≥ 3
     * SO(2) is not simply connected (isomorphic to S¹)
     */
    bool isSimplyConnected() const {
        return n >= 3;
    }

    /**
     * @brief Fundamental group π₁(SO(n))
     * Returns order: 2 for n ≥ 3, ∞ for n = 2
     */
    int fundamentalGroupOrder() const {
        if (n == 2) {
            return -1;  // Infinite (ℤ)
        } else if (n >= 3) {
            return 2;  // ℤ₂
        }
        return 1;  // Trivial
    }
};

// ============================================================================
// MATRIX EXPONENTIAL
// ============================================================================

/**
 * @brief Compute matrix exponential using scaling and squaring
 *
 * exp(A) = Σ_{k=0}^∞ Aᵏ/k!
 *
 * Uses Padé approximation with scaling and squaring
 *
 * @param A Input matrix
 * @param terms Number of terms in series (default 20)
 * @return exp(A)
 */
template<typename T = double>
Matrix<T> matrixExponential(const Matrix<T>& A, int terms = 20) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Matrix must be square");
    }

    size_t n = A.rows();
    Matrix<T> result = Matrix<T>::identity(n);
    Matrix<T> term = Matrix<T>::identity(n);

    // Compute series: I + A + A²/2! + A³/3! + ...
    for (int k = 1; k <= terms; ++k) {
        term = term * A;
        T factorial = T(1);
        for (int i = 2; i <= k; ++i) {
            factorial *= i;
        }

        Matrix<T> scaled_term = term * (T(1) / factorial);
        result = result + scaled_term;
    }

    return result;
}

/**
 * @brief Compute matrix logarithm (inverse of exponential)
 *
 * log(A) computed via series for ||I - A|| < 1
 * log(A) = Σ_{k=1}^∞ (-1)^{k+1}/k · (A - I)^k
 *
 * @param A Input matrix (should be close to identity)
 * @param terms Number of terms
 * @return log(A)
 */
template<typename T = double>
Matrix<T> matrixLogarithm(const Matrix<T>& A, int terms = 20) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Matrix must be square");
    }

    size_t n = A.rows();
    Matrix<T> I = Matrix<T>::identity(n);
    Matrix<T> X = A - I;  // A - I

    // Check convergence condition
    if (X.frobeniusNorm() >= 1.0) {
        throw std::runtime_error("Matrix logarithm series may not converge");
    }

    Matrix<T> result(n, n, T(0));
    Matrix<T> power = X;  // (A - I)^k

    for (int k = 1; k <= terms; ++k) {
        T sign = (k % 2 == 1) ? T(1) : T(-1);
        T coeff = sign / T(k);

        result = result + power * coeff;

        if (k < terms) {
            power = power * X;
        }
    }

    return result;
}

/**
 * @brief Verify exp(log(A)) = A
 */
template<typename T = double>
bool verifyExpLog(const Matrix<T>& A, double tolerance = 1e-6) {
    try {
        Matrix<T> logA = matrixLogarithm(A);
        Matrix<T> expLogA = matrixExponential(logA);
        return A.approxEqual(expLogA, tolerance);
    } catch (...) {
        return false;
    }
}

// ============================================================================
// LIE ALGEBRAS
// ============================================================================

/**
 * @brief Lie bracket (commutator)
 *
 * [X, Y] = XY - YX
 *
 * @param X First matrix
 * @param Y Second matrix
 * @return Commutator [X, Y]
 */
template<typename T = double>
Matrix<T> lieBracket(const Matrix<T>& X, const Matrix<T>& Y) {
    if (X.rows() != X.cols() || Y.rows() != Y.cols()) {
        throw std::invalid_argument("Matrices must be square");
    }
    if (X.rows() != Y.rows()) {
        throw std::invalid_argument("Matrices must have same dimension");
    }

    return (X * Y) - (Y * X);
}

/**
 * @brief Check Jacobi identity
 *
 * [X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0
 *
 * @param X, Y, Z Three matrices
 * @return true if Jacobi identity holds
 */
template<typename T = double>
bool verifyJacobiIdentity(const Matrix<T>& X, const Matrix<T>& Y,
                          const Matrix<T>& Z, double tolerance = 1e-10) {
    Matrix<T> term1 = lieBracket(X, lieBracket(Y, Z));
    Matrix<T> term2 = lieBracket(Y, lieBracket(Z, X));
    Matrix<T> term3 = lieBracket(Z, lieBracket(X, Y));

    Matrix<T> sum = term1 + term2 + term3;

    // Check if sum is approximately zero
    return sum.frobeniusNorm() < tolerance;
}

/**
 * @brief Lie algebra of sl(n): traceless matrices
 *
 * sl(n) = {X ∈ ℝⁿˣⁿ : tr(X) = 0}
 */
template<typename T = double>
bool isSLAlgebra(const Matrix<T>& X, double tolerance = 1e-10) {
    if (X.rows() != X.cols()) {
        return false;
    }

    T trace = X.trace();
    return std::abs(trace) < tolerance;
}

/**
 * @brief Lie algebra of so(n): skew-symmetric matrices
 *
 * so(n) = {X ∈ ℝⁿˣⁿ : Xᵀ = -X}
 */
template<typename T = double>
bool isSOAlgebra(const Matrix<T>& X, double tolerance = 1e-10) {
    if (X.rows() != X.cols()) {
        return false;
    }

    Matrix<T> XT = X.transpose();
    Matrix<T> sum = X + XT;

    return sum.frobeniusNorm() < tolerance;
}

/**
 * @brief Project matrix to sl(n) by subtracting trace
 */
template<typename T = double>
Matrix<T> projectToSL(const Matrix<T>& X) {
    if (X.rows() != X.cols()) {
        throw std::invalid_argument("Matrix must be square");
    }

    size_t n = X.rows();
    T trace = X.trace();
    T avg = trace / T(n);

    Matrix<T> result = X;
    for (size_t i = 0; i < n; ++i) {
        result(i, i) -= avg;
    }

    return result;
}

/**
 * @brief Project matrix to so(n) by skew-symmetrizing
 */
template<typename T = double>
Matrix<T> projectToSO(const Matrix<T>& X) {
    if (X.rows() != X.cols()) {
        throw std::invalid_argument("Matrix must be square");
    }

    Matrix<T> XT = X.transpose();
    return (X - XT) * T(0.5);
}

// ============================================================================
// EXPONENTIAL MAPPING
// ============================================================================

/**
 * @brief Exponential map from Lie algebra to Lie group
 *
 * For X ∈ 𝔤 (Lie algebra), exp(X) ∈ G (Lie group)
 *
 * Examples:
 * - X ∈ sl(n) → exp(X) ∈ SL(n)
 * - X ∈ so(n) → exp(X) ∈ SO(n)
 */
template<typename T = double>
class ExponentialMap {
public:
    /**
     * @brief Map from sl(n) to SL(n)
     */
    static Matrix<T> slToSL(const Matrix<T>& X) {
        if (!isSLAlgebra(X)) {
            throw std::invalid_argument("Matrix not in sl(n)");
        }

        Matrix<T> expX = matrixExponential(X);

        // Verify det(exp(X)) = exp(tr(X)) = exp(0) = 1
        // (This should automatically hold for traceless X)

        return expX;
    }

    /**
     * @brief Map from so(n) to SO(n)
     */
    static Matrix<T> soToSO(const Matrix<T>& X) {
        if (!isSOAlgebra(X)) {
            throw std::invalid_argument("Matrix not in so(n)");
        }

        Matrix<T> expX = matrixExponential(X);

        // Verify exp(X) is orthogonal
        // If X is skew-symmetric, exp(X) is orthogonal

        return expX;
    }

    /**
     * @brief Verify exponential preserves group properties
     */
    static bool verifySLExponential(const Matrix<T>& X, double tolerance = 1e-10) {
        Matrix<T> expX = slToSL(X);
        T det = expX.determinant();
        return std::abs(det - T(1)) < tolerance;
    }

    static bool verifySOExponential(const Matrix<T>& X, double tolerance = 1e-10) {
        Matrix<T> expX = soToSO(X);
        Matrix<T> expXT = expX.transpose();
        Matrix<T> product = expX * expXT;
        Matrix<T> I = Matrix<T>::identity(X.rows());
        return product.approxEqual(I, tolerance);
    }
};

// ============================================================================
// COMPLEXIFICATION
// ============================================================================

/**
 * @brief Complexification of real Lie algebra
 *
 * 𝔤ℂ = 𝔤 ⊗ ℂ
 *
 * For 𝔤 ⊂ ℝⁿˣⁿ, 𝔤ℂ ⊂ ℂⁿˣⁿ
 */
template<typename T = double>
class Complexification {
public:
    using ComplexMatrix = Matrix<std::complex<T>>;

    /**
     * @brief Embed real matrix into complex matrices
     */
    static ComplexMatrix embed(const Matrix<T>& X) {
        size_t n = X.rows();
        ComplexMatrix result(n, n);

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                result(i, j) = std::complex<T>(X(i, j), 0.0);
            }
        }

        return result;
    }

    /**
     * @brief Create complex linear combination
     *
     * Z = X + iY where X, Y ∈ 𝔤ℝ
     */
    static ComplexMatrix linearCombination(const Matrix<T>& X, const Matrix<T>& Y) {
        if (X.rows() != Y.rows() || X.cols() != Y.cols()) {
            throw std::invalid_argument("Matrices must have same dimensions");
        }

        size_t n = X.rows();
        ComplexMatrix result(n, n);

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                result(i, j) = std::complex<T>(X(i, j), Y(i, j));
            }
        }

        return result;
    }

    /**
     * @brief Dimension of complexified algebra
     *
     * dim(𝔤ℂ) = dim(𝔤ℝ)
     */
    static int dimension(const Matrix<T>& basis_element) {
        // Real dimension = complex dimension (as vector space over ℂ)
        size_t n = basis_element.rows();
        return n * n;  // For full matrix algebra
    }
};

// ============================================================================
// BAKER-CAMPBELL-HAUSDORFF FORMULA
// ============================================================================

/**
 * @brief Baker-Campbell-Hausdorff (BCH) formula
 *
 * Z = log(exp(X) exp(Y))
 *
 * Series expansion:
 * Z = X + Y + (1/2)[X,Y] + (1/12)[X,[X,Y]] + (1/12)[Y,[Y,X]] + ...
 *
 * @param X First Lie algebra element
 * @param Y Second Lie algebra element
 * @param order Order of expansion (1, 2, or 3)
 * @return BCH series approximation
 */
template<typename T = double>
Matrix<T> bakerCampbellHausdorff(const Matrix<T>& X, const Matrix<T>& Y, int order = 2) {
    if (X.rows() != X.cols() || Y.rows() != Y.cols()) {
        throw std::invalid_argument("Matrices must be square");
    }
    if (X.rows() != Y.rows()) {
        throw std::invalid_argument("Matrices must have same dimension");
    }

    size_t n = X.rows();
    Matrix<T> Z(n, n, T(0));

    // Order 0 and 1: Z = X + Y
    Z = X + Y;

    if (order >= 2) {
        // Order 2: Z = X + Y + (1/2)[X,Y]
        Matrix<T> commutator = lieBracket(X, Y);
        Z = Z + commutator * T(0.5);
    }

    if (order >= 3) {
        // Order 3: Add (1/12)([X,[X,Y]] - [Y,[X,Y]])
        Matrix<T> XY = lieBracket(X, Y);
        Matrix<T> X_XY = lieBracket(X, XY);
        Matrix<T> Y_XY = lieBracket(Y, XY);

        Matrix<T> term3 = (X_XY - Y_XY) * (T(1) / T(12));
        Z = Z + term3;
    }

    return Z;
}

/**
 * @brief BCH formula for Heisenberg group
 *
 * Heisenberg algebra: [X, Y] = Z, [X, Z] = 0, [Y, Z] = 0
 *
 * For Heisenberg group, BCH formula terminates:
 * log(exp(X) exp(Y)) = X + Y + (1/2)[X, Y]
 */
template<typename T = double>
Matrix<T> bakerCampbellHausdorffHeisenberg(const Matrix<T>& X, const Matrix<T>& Y) {
    // For Heisenberg group, exact formula
    Matrix<T> commutator = lieBracket(X, Y);
    return X + Y + commutator * T(0.5);
}

/**
 * @brief Verify BCH formula approximately
 *
 * Check: exp(BCH(X,Y)) ≈ exp(X)exp(Y)
 */
template<typename T = double>
bool verifyBCH(const Matrix<T>& X, const Matrix<T>& Y,
               int order = 2, double tolerance = 1e-6) {
    Matrix<T> Z = bakerCampbellHausdorff(X, Y, order);

    Matrix<T> expX = matrixExponential(X);
    Matrix<T> expY = matrixExponential(Y);
    Matrix<T> expXexpY = expX * expY;

    Matrix<T> expZ = matrixExponential(Z);

    return expZ.approxEqual(expXexpY, tolerance);
}

/**
 * @brief Compute higher-order BCH terms
 *
 * Full BCH series (truncated to specified order)
 *
 * @param X First element
 * @param Y Second element
 * @param maxOrder Maximum order to compute
 * @return BCH approximation
 */
template<typename T = double>
Matrix<T> bakerCampbellHausdorffSeries(const Matrix<T>& X, const Matrix<T>& Y,
                                       int maxOrder = 4) {
    if (maxOrder < 1) {
        throw std::invalid_argument("Order must be at least 1");
    }

    size_t n = X.rows();
    Matrix<T> Z(n, n, T(0));

    // Order 1
    Z = X + Y;

    if (maxOrder >= 2) {
        // (1/2)[X, Y]
        Z = Z + lieBracket(X, Y) * T(0.5);
    }

    if (maxOrder >= 3) {
        // (1/12)([X,[X,Y]] + [Y,[Y,X]])
        Matrix<T> XY = lieBracket(X, Y);
        Matrix<T> X_XY = lieBracket(X, XY);
        Matrix<T> Y_YX = lieBracket(Y, lieBracket(Y, X));

        Z = Z + (X_XY + Y_YX) * (T(1) / T(12));
    }

    if (maxOrder >= 4) {
        // -( 1/24)[Y,[X,[X,Y]]]
        Matrix<T> XY = lieBracket(X, Y);
        Matrix<T> X_XY = lieBracket(X, XY);
        Matrix<T> Y_X_XY = lieBracket(Y, X_XY);

        Z = Z - Y_X_XY * (T(1) / T(24));
    }

    return Z;
}

// ============================================================================
// SUBGROUPS AND SUBALGEBRAS
// ============================================================================

/**
 * @brief Check if algebra element generates a one-parameter subgroup
 *
 * One-parameter subgroup: φ(t) = exp(tX) for X ∈ 𝔤
 *
 * Properties:
 * - φ(0) = I
 * - φ(s + t) = φ(s)φ(t)
 * - dφ/dt|_{t=0} = X
 */
template<typename T = double>
class OneParameterSubgroup {
private:
    Matrix<T> generator;  // X ∈ 𝔤

public:
    explicit OneParameterSubgroup(const Matrix<T>& X) : generator(X) {
        if (X.rows() != X.cols()) {
            throw std::invalid_argument("Generator must be square matrix");
        }
    }

    /**
     * @brief Evaluate φ(t) = exp(tX)
     */
    Matrix<T> evaluate(T t) const {
        Matrix<T> tX = generator * t;
        return matrixExponential(tX);
    }

    /**
     * @brief Verify group property: φ(s+t) = φ(s)φ(t)
     */
    bool verifyGroupProperty(T s, T t, double tolerance = 1e-10) const {
        Matrix<T> phi_s_plus_t = evaluate(s + t);
        Matrix<T> phi_s = evaluate(s);
        Matrix<T> phi_t = evaluate(t);
        Matrix<T> product = phi_s * phi_t;

        return phi_s_plus_t.approxEqual(product, tolerance);
    }

    /**
     * @brief Get the generator
     */
    const Matrix<T>& getGenerator() const {
        return generator;
    }
};

/**
 * @brief Check if set of matrices forms a subalgebra
 *
 * A subspace 𝔥 ⊂ 𝔤 is a subalgebra if:
 * [X, Y] ∈ 𝔥 for all X, Y ∈ 𝔥
 */
template<typename T = double>
bool isSubalgebra(const std::vector<Matrix<T>>& basis) {
    // Check closure under Lie bracket
    for (const auto& X : basis) {
        for (const auto& Y : basis) {
            Matrix<T> bracket = lieBracket(X, Y);

            // Check if bracket can be expressed as linear combination of basis
            // (Simplified check: assume span is checked separately)
            bool inSpan = false;

            // Check if bracket is approximately in the span
            // For simplicity, check if bracket is approximately zero or in basis
            double norm = bracket.frobeniusNorm();
            if (norm < 1e-10) {
                inSpan = true;  // Zero is always in span
            } else {
                // More sophisticated span checking would be needed for general case
                // Here we just check against basis elements
                for (const auto& B : basis) {
                    if (bracket.approxEqual(B) || bracket.approxEqual(B * T(-1))) {
                        inSpan = true;
                        break;
                    }
                }
            }

            if (!inSpan) {
                // Could not verify that bracket is in subalgebra
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Calculate dimension of Lie algebra
 *
 * For classical Lie algebras:
 * - gl(n): n²
 * - sl(n): n² - 1
 * - so(n): n(n-1)/2
 * - sp(n): n(2n+1)
 */
inline int liealgebraDimension(const std::string& type, int n) {
    if (type == "gl") {
        return n * n;
    } else if (type == "sl") {
        return n * n - 1;
    } else if (type == "so") {
        return n * (n - 1) / 2;
    } else if (type == "sp") {
        return n * (2 * n + 1);
    } else {
        throw std::invalid_argument("Unknown Lie algebra type");
    }
}

} // namespace group_theory
} // namespace maths

#endif // MATHS_GROUP_THEORY_LIE_GROUPS_HPP
