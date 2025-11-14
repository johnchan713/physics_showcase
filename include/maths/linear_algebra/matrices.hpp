#ifndef MATHS_LINEAR_ALGEBRA_MATRICES_HPP
#define MATHS_LINEAR_ALGEBRA_MATRICES_HPP

#include "vectors.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <sstream>

/**
 * @file matrices.hpp
 * @brief Matrix operations and properties
 *
 * Implements:
 * - Matrix arithmetic
 * - Matrix multiplication
 * - Determinant and trace
 * - Inverse and transpose
 * - Eigenvalues and eigenvectors (basic)
 */

namespace maths::linear_algebra {

/**
 * @class Matrix
 * @brief m × n matrix with operations
 */
class Matrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows_;
    size_t cols_;

public:
    /**
     * @brief Construct matrix from 2D vector
     */
    explicit Matrix(const std::vector<std::vector<double>>& d)
        : data(d), rows_(d.size()), cols_(d.empty() ? 0 : d[0].size()) {
        // Verify all rows have same length
        for (const auto& row : data) {
            if (row.size() != cols_) {
                throw std::invalid_argument("All rows must have same length");
            }
        }
    }

    /**
     * @brief Construct zero matrix of size m × n
     */
    Matrix(size_t rows, size_t cols)
        : data(rows, std::vector<double>(cols, 0.0)),
          rows_(rows), cols_(cols) {}

    /**
     * @brief Get dimensions
     */
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }

    /**
     * @brief Access element (mutable)
     */
    double& operator()(size_t i, size_t j) {
        if (i >= rows_ || j >= cols_) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i][j];
    }

    /**
     * @brief Access element (const)
     */
    const double& operator()(size_t i, size_t j) const {
        if (i >= rows_ || j >= cols_) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i][j];
    }

    /**
     * @brief Matrix addition
     */
    Matrix operator+(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Matrices must have same dimensions");
        }

        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = data[i][j] + other(i, j);
            }
        }
        return result;
    }

    /**
     * @brief Matrix subtraction
     */
    Matrix operator-(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Matrices must have same dimensions");
        }

        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = data[i][j] - other(i, j);
            }
        }
        return result;
    }

    /**
     * @brief Scalar multiplication
     */
    Matrix operator*(double scalar) const {
        Matrix result(rows_, cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = data[i][j] * scalar;
            }
        }
        return result;
    }

    /**
     * @brief Matrix multiplication
     *
     * (m × n) × (n × p) = (m × p)
     */
    Matrix operator*(const Matrix& other) const {
        if (cols_ != other.rows_) {
            throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
        }

        Matrix result(rows_, other.cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < other.cols_; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < cols_; ++k) {
                    sum += data[i][k] * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    /**
     * @brief Matrix-vector multiplication
     */
    Vector operator*(const Vector& v) const {
        if (cols_ != v.dimension()) {
            throw std::invalid_argument("Matrix columns must match vector dimension");
        }

        std::vector<double> result(rows_);
        for (size_t i = 0; i < rows_; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < cols_; ++j) {
                sum += data[i][j] * v[j];
            }
            result[i] = sum;
        }
        return Vector(result);
    }

    /**
     * @brief Transpose
     */
    Matrix transpose() const {
        Matrix result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(j, i) = data[i][j];
            }
        }
        return result;
    }

    /**
     * @brief Trace (sum of diagonal elements, square matrices only)
     */
    double trace() const {
        if (rows_ != cols_) {
            throw std::runtime_error("Trace requires square matrix");
        }

        double sum = 0.0;
        for (size_t i = 0; i < rows_; ++i) {
            sum += data[i][i];
        }
        return sum;
    }

    /**
     * @brief Determinant (square matrices only, using Laplace expansion)
     */
    double determinant() const {
        if (rows_ != cols_) {
            throw std::runtime_error("Determinant requires square matrix");
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

        // Laplace expansion for larger matrices (inefficient but general)
        double det = 0.0;
        for (size_t j = 0; j < cols_; ++j) {
            double cofactor = ((j % 2 == 0) ? 1.0 : -1.0) * data[0][j];
            det += cofactor * minor(0, j).determinant();
        }
        return det;
    }

    /**
     * @brief Get minor (submatrix by removing row i and column j)
     */
    Matrix minor(size_t row, size_t col) const {
        Matrix result(rows_ - 1, cols_ - 1);
        size_t r = 0;
        for (size_t i = 0; i < rows_; ++i) {
            if (i == row) continue;
            size_t c = 0;
            for (size_t j = 0; j < cols_; ++j) {
                if (j == col) continue;
                result(r, c) = data[i][j];
                ++c;
            }
            ++r;
        }
        return result;
    }

    /**
     * @brief Identity matrix
     */
    static Matrix identity(size_t n) {
        Matrix I(n, n);
        for (size_t i = 0; i < n; ++i) {
            I(i, i) = 1.0;
        }
        return I;
    }

    /**
     * @brief Check if matrix is square
     */
    bool isSquare() const {
        return rows_ == cols_;
    }

    /**
     * @brief Check if matrix is symmetric (A = Aᵀ)
     */
    bool isSymmetric(double tolerance = 1e-10) const {
        if (!isSquare()) return false;

        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = i + 1; j < cols_; ++j) {
                if (std::abs(data[i][j] - data[j][i]) > tolerance) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Check if matrix is diagonal
     */
    bool isDiagonal(double tolerance = 1e-10) const {
        if (!isSquare()) return false;

        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                if (i != j && std::abs(data[i][j]) > tolerance) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Check if matrix is identity
     */
    bool isIdentity(double tolerance = 1e-10) const {
        if (!isSquare()) return false;

        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                double expected = (i == j) ? 1.0 : 0.0;
                if (std::abs(data[i][j] - expected) > tolerance) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief String representation
     */
    std::string toString() const {
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < rows_; ++i) {
            if (i > 0) oss << " ";
            oss << "[";
            for (size_t j = 0; j < cols_; ++j) {
                oss << data[i][j];
                if (j < cols_ - 1) oss << ", ";
            }
            oss << "]";
            if (i < rows_ - 1) oss << "\n";
        }
        oss << "]";
        return oss.str();
    }
};

/**
 * @brief Scalar multiplication (left side)
 */
inline Matrix operator*(double scalar, const Matrix& M) {
    return M * scalar;
}

/**
 * @class MatrixProperties
 * @brief Properties and theorems about matrices
 */
class MatrixProperties {
public:
    /**
     * @brief Determinant properties
     */
    static std::string determinantProperties() {
        return "Determinant Properties:\n"
               "\n"
               "1. det(AB) = det(A)det(B)\n"
               "2. det(Aᵀ) = det(A)\n"
               "3. det(A⁻¹) = 1/det(A) (if A invertible)\n"
               "4. det(cA) = cⁿ det(A) (n×n matrix)\n"
               "5. det(I) = 1\n"
               "6. Row swap: multiplies det by -1\n"
               "7. Row scaling: multiplies det by scale\n"
               "8. det(A) = 0 ⟺ A singular (not invertible)";
    }

    /**
     * @brief Trace properties
     */
    static std::string traceProperties() {
        return "Trace Properties:\n"
               "\n"
               "1. tr(A + B) = tr(A) + tr(B)\n"
               "2. tr(cA) = c·tr(A)\n"
               "3. tr(Aᵀ) = tr(A)\n"
               "4. tr(AB) = tr(BA) (cyclic property)\n"
               "5. tr(ABC) = tr(CAB) = tr(BCA)\n"
               "6. tr(I) = n (for n×n identity)\n"
               "7. tr(A) = sum of eigenvalues";
    }

    /**
     * @brief Transpose properties
     */
    static std::string transposeProperties() {
        return "Transpose Properties:\n"
               "\n"
               "1. (Aᵀ)ᵀ = A\n"
               "2. (A + B)ᵀ = Aᵀ + Bᵀ\n"
               "3. (cA)ᵀ = cAᵀ\n"
               "4. (AB)ᵀ = BᵀAᵀ (reverse order!)\n"
               "5. (ABC)ᵀ = CᵀBᵀAᵀ\n"
               "6. (A⁻¹)ᵀ = (Aᵀ)⁻¹";
    }

    /**
     * @brief Matrix multiplication properties
     */
    static std::string multiplicationProperties() {
        return "Matrix Multiplication:\n"
               "\n"
               "1. Associative: (AB)C = A(BC)\n"
               "2. Distributive: A(B + C) = AB + AC\n"
               "3. NOT commutative: AB ≠ BA (in general)\n"
               "4. Identity: AI = IA = A\n"
               "5. Zero: A·0 = 0·A = 0\n"
               "\n"
               "Special cases where AB = BA:\n"
               "- A = cI (scalar matrix)\n"
               "- B = A⁻¹\n"
               "- Both diagonal\n"
               "- Certain symmetric matrices";
    }

    /**
     * @brief Invertibility conditions
     */
    static std::string invertibilityConditions() {
        return "Matrix Invertibility (n×n matrix A):\n"
               "\n"
               "The following are equivalent:\n"
               "1. A is invertible (∃A⁻¹ such that AA⁻¹ = I)\n"
               "2. det(A) ≠ 0\n"
               "3. rank(A) = n (full rank)\n"
               "4. Columns of A are linearly independent\n"
               "5. Rows of A are linearly independent\n"
               "6. Ax = 0 has only trivial solution x = 0\n"
               "7. Ax = b has unique solution for all b\n"
               "8. 0 is not an eigenvalue of A";
    }
};

/**
 * @class SpecialMatrices
 * @brief Special types of matrices
 */
class SpecialMatrices {
public:
    /**
     * @brief Create diagonal matrix from vector
     */
    static Matrix diagonal(const Vector& diag) {
        size_t n = diag.dimension();
        Matrix D(n, n);
        for (size_t i = 0; i < n; ++i) {
            D(i, i) = diag[i];
        }
        return D;
    }

    /**
     * @brief Symmetric matrix definition
     */
    static std::string symmetricMatrices() {
        return "Symmetric Matrices (A = Aᵀ):\n"
               "\n"
               "Properties:\n"
               "1. aᵢⱼ = aⱼᵢ for all i, j\n"
               "2. All eigenvalues are real\n"
               "3. Eigenvectors from different eigenvalues are orthogonal\n"
               "4. Can be diagonalized by orthogonal matrix:\n"
               "   A = QΛQᵀ where Q orthogonal, Λ diagonal\n"
               "\n"
               "Examples: covariance matrices, Hessian matrices";
    }

    /**
     * @brief Orthogonal matrix definition
     */
    static std::string orthogonalMatrices() {
        return "Orthogonal Matrices (QᵀQ = I):\n"
               "\n"
               "Properties:\n"
               "1. Q⁻¹ = Qᵀ (inverse = transpose)\n"
               "2. Columns form orthonormal basis\n"
               "3. Rows form orthonormal basis\n"
               "4. Preserves lengths: ‖Qx‖ = ‖x‖\n"
               "5. Preserves angles: (Qx)·(Qy) = x·y\n"
               "6. det(Q) = ±1\n"
               "\n"
               "Examples: rotation matrices, reflection matrices";
    }

    /**
     * @brief Positive definite matrices
     */
    static std::string positiveDefiniteMatrices() {
        return "Positive Definite Matrices:\n"
               "\n"
               "Symmetric matrix A is positive definite if:\n"
               "xᵀAx > 0 for all x ≠ 0\n"
               "\n"
               "Equivalent conditions:\n"
               "1. All eigenvalues > 0\n"
               "2. All leading principal minors > 0\n"
               "3. A = BᵀB for some invertible B\n"
               "4. Cholesky decomposition exists: A = LLᵀ\n"
               "\n"
               "Positive semidefinite: xᵀAx ≥ 0 (allows zero)";
    }

    /**
     * @brief Upper/lower triangular matrices
     */
    static std::string triangularMatrices() {
        return "Triangular Matrices:\n"
               "\n"
               "Upper triangular: aᵢⱼ = 0 for i > j\n"
               "[a₁₁ a₁₂ a₁₃]\n"
               "[ 0  a₂₂ a₂₃]\n"
               "[ 0   0  a₃₃]\n"
               "\n"
               "Lower triangular: aᵢⱼ = 0 for i < j\n"
               "[a₁₁  0   0 ]\n"
               "[a₂₁ a₂₂  0 ]\n"
               "[a₃₁ a₃₂ a₃₃]\n"
               "\n"
               "Properties:\n"
               "- det(A) = product of diagonal elements\n"
               "- Product of triangular matrices is triangular\n"
               "- Eigenvalues = diagonal elements";
    }
};

/**
 * @class MatrixDecompositions
 * @brief Matrix factorizations
 */
class MatrixDecompositions {
public:
    /**
     * @brief LU decomposition
     */
    static std::string luDecomposition() {
        return "LU Decomposition:\n"
               "\n"
               "A = LU\n"
               "\n"
               "where L is lower triangular, U is upper triangular\n"
               "\n"
               "Uses:\n"
               "- Solve Ax = b: Ly = b (forward sub), Ux = y (backward sub)\n"
               "- Compute determinant: det(A) = det(L)det(U)\n"
               "- Find inverse\n"
               "\n"
               "With pivoting: PA = LU (P permutation matrix)";
    }

    /**
     * @brief QR decomposition
     */
    static std::string qrDecomposition() {
        return "QR Decomposition:\n"
               "\n"
               "A = QR\n"
               "\n"
               "where Q is orthogonal (QᵀQ = I), R is upper triangular\n"
               "\n"
               "Methods:\n"
               "- Gram-Schmidt orthogonalization\n"
               "- Householder reflections\n"
               "- Givens rotations\n"
               "\n"
               "Uses:\n"
               "- Least squares: Ax ≈ b\n"
               "- Eigenvalue algorithms (QR iteration)\n"
               "- Stable computation";
    }

    /**
     * @brief Singular Value Decomposition (SVD)
     */
    static std::string svd() {
        return "Singular Value Decomposition:\n"
               "\n"
               "A = UΣVᵀ\n"
               "\n"
               "where U, V orthogonal, Σ diagonal (singular values)\n"
               "\n"
               "Properties:\n"
               "- Exists for any m×n matrix\n"
               "- Singular values σᵢ ≥ 0\n"
               "- rank(A) = number of nonzero singular values\n"
               "- ‖A‖₂ = σ_max (largest singular value)\n"
               "\n"
               "Uses:\n"
               "- Pseudoinverse: A⁺ = VΣ⁺Uᵀ\n"
               "- Low-rank approximation\n"
               "- Principal component analysis (PCA)";
    }

    /**
     * @brief Eigenvalue decomposition
     */
    static std::string eigenvalueDecomposition() {
        return "Eigenvalue Decomposition:\n"
               "\n"
               "A = PDP⁻¹\n"
               "\n"
               "where D diagonal (eigenvalues), P (eigenvectors)\n"
               "\n"
               "Requirements:\n"
               "- A must be diagonalizable\n"
               "- n linearly independent eigenvectors\n"
               "\n"
               "For symmetric A:\n"
               "A = QΛQᵀ (Q orthogonal, Λ diagonal)\n"
               "Spectral theorem: symmetric matrices always diagonalizable\n"
               "\n"
               "Uses:\n"
               "- Matrix powers: Aⁿ = PDⁿP⁻¹\n"
               "- Matrix exponential: e^A\n"
               "- Differential equations";
    }

    /**
     * @brief Cholesky decomposition
     */
    static std::string choleskyDecomposition() {
        return "Cholesky Decomposition:\n"
               "\n"
               "A = LLᵀ\n"
               "\n"
               "where L is lower triangular\n"
               "\n"
               "Requirements:\n"
               "- A must be symmetric positive definite\n"
               "\n"
               "Properties:\n"
               "- Unique (if diagonal elements of L positive)\n"
               "- More efficient than LU (half the operations)\n"
               "- Numerically stable\n"
               "\n"
               "Uses:\n"
               "- Solve symmetric positive definite systems\n"
               "- Monte Carlo simulation (correlated random variables)\n"
               "- Kalman filtering";
    }
};

/**
 * @class LinearTransformations
 * @brief Matrices as linear transformations
 */
class LinearTransformations {
public:
    /**
     * @brief Definition of linear transformation
     */
    static std::string definition() {
        return "Linear Transformation T: ℝⁿ → ℝᵐ:\n"
               "\n"
               "Properties:\n"
               "1. T(u + v) = T(u) + T(v) (additive)\n"
               "2. T(cu) = cT(u) (homogeneous)\n"
               "\n"
               "Equivalently: T(cu + dv) = cT(u) + dT(v)\n"
               "\n"
               "Matrix representation:\n"
               "T(x) = Ax\n"
               "\n"
               "Every linear transformation can be represented by matrix!";
    }

    /**
     * @brief Kernel (null space)
     */
    static std::string kernel() {
        return "Kernel (Null Space):\n"
               "\n"
               "ker(A) = {x : Ax = 0}\n"
               "\n"
               "Properties:\n"
               "- Subspace of domain\n"
               "- dim(ker(A)) = nullity(A)\n"
               "- A injective ⟺ ker(A) = {0}\n"
               "\n"
               "Finding kernel: row reduce [A|0] to RREF";
    }

    /**
     * @brief Range (column space)
     */
    static std::string range() {
        return "Range (Column Space):\n"
               "\n"
               "range(A) = {Ax : x ∈ ℝⁿ}\n"
               "         = span of columns of A\n"
               "\n"
               "Properties:\n"
               "- Subspace of codomain\n"
               "- dim(range(A)) = rank(A)\n"
               "- A surjective ⟺ range(A) = ℝᵐ\n"
               "\n"
               "Finding range: identify pivot columns of RREF";
    }

    /**
     * @brief Rank-Nullity Theorem
     */
    static std::string rankNullityTheorem() {
        return "Rank-Nullity Theorem:\n"
               "\n"
               "For A: ℝⁿ → ℝᵐ:\n"
               "\n"
               "rank(A) + nullity(A) = n\n"
               "\n"
               "i.e., dim(range(A)) + dim(ker(A)) = dim(domain)\n"
               "\n"
               "Interpretation:\n"
               "Domain splits into kernel (maps to 0)\n"
               "and complement (isomorphic to range)";
    }
};

} // namespace maths::linear_algebra

#endif // MATHS_LINEAR_ALGEBRA_MATRICES_HPP
