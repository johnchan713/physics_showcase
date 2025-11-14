#include "../include/maths/linear_algebra/matrices.hpp"
#include "../include/maths/linear_algebra/vectors.hpp"
#include <iostream>
#include <iomanip>

using namespace maths::linear_algebra;

void printMatrix(const std::string& name, const Matrix& M) {
    std::cout << "\n" << name << ":\n" << M.toString() << "\n";
}

void printValue(const std::string& name, double value) {
    std::cout << name << ": " << value << "\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(4);

    // ========================================
    // 1. Creating Matrices
    // ========================================
    std::cout << "=== MATRIX CREATION ===\n";

    Matrix A({{1, 2, 3},
              {4, 5, 6},
              {7, 8, 9}});
    printMatrix("Matrix A (3x3)", A);

    Matrix B({{9, 8, 7},
              {6, 5, 4},
              {3, 2, 1}});
    printMatrix("Matrix B (3x3)", B);

    Matrix I = Matrix::identity(3);
    printMatrix("Identity Matrix I (3x3)", I);

    // ========================================
    // 2. Basic Arithmetic
    // ========================================
    std::cout << "\n=== BASIC ARITHMETIC ===\n";

    Matrix sum = A + B;
    printMatrix("A + B", sum);

    Matrix diff = A - B;
    printMatrix("A - B", diff);

    Matrix scaled = A * 2.5;
    printMatrix("A * 2.5", scaled);

    // ========================================
    // 3. Matrix Multiplication
    // ========================================
    std::cout << "\n=== MATRIX MULTIPLICATION ===\n";

    Matrix C({{1, 2},
              {3, 4},
              {5, 6}});

    Matrix D({{7, 8, 9},
              {10, 11, 12}});

    printMatrix("Matrix C (3x2)", C);
    printMatrix("Matrix D (2x3)", D);

    Matrix CD = C * D;
    printMatrix("C * D (3x3)", CD);

    Matrix DC = D * C;
    printMatrix("D * C (2x2)", DC);

    // Matrix-vector multiplication
    Vector v({1, 2, 3});
    std::cout << "\nVector v: " << v.toString() << "\n";
    Vector Av = A * v;
    std::cout << "A * v: " << Av.toString() << "\n";

    // ========================================
    // 4. Transpose
    // ========================================
    std::cout << "\n=== TRANSPOSE ===\n";

    printMatrix("Original A", A);
    Matrix At = A.transpose();
    printMatrix("A^T (transpose)", At);

    // Verify (A^T)^T = A
    Matrix Att = At.transpose();
    printMatrix("(A^T)^T", Att);

    // ========================================
    // 5. Trace and Determinant
    // ========================================
    std::cout << "\n=== TRACE AND DETERMINANT ===\n";

    printValue("trace(A)", A.trace());
    printValue("det(A)", A.determinant());
    printValue("det(I)", I.determinant());

    // 2x2 determinant example
    Matrix M2({{3, 8},
               {4, 6}});
    printMatrix("Matrix M2 (2x2)", M2);
    printValue("det(M2) = 3*6 - 8*4", M2.determinant());

    // 3x3 determinant example
    Matrix M3({{6, 1, 1},
               {4, -2, 5},
               {2, 8, 7}});
    printMatrix("Matrix M3 (3x3)", M3);
    printValue("det(M3)", M3.determinant());

    // ========================================
    // 6. Matrix Properties
    // ========================================
    std::cout << "\n=== MATRIX PROPERTIES ===\n";

    Matrix symmetric({{1, 2, 3},
                      {2, 4, 5},
                      {3, 5, 6}});
    printMatrix("Symmetric Matrix", symmetric);
    std::cout << "Is symmetric? " << (symmetric.isSymmetric() ? "Yes" : "No") << "\n";
    std::cout << "Is diagonal? " << (symmetric.isDiagonal() ? "Yes" : "No") << "\n";

    Matrix diagonal({{5, 0, 0},
                     {0, 3, 0},
                     {0, 0, -2}});
    printMatrix("Diagonal Matrix", diagonal);
    std::cout << "Is symmetric? " << (diagonal.isSymmetric() ? "Yes" : "No") << "\n";
    std::cout << "Is diagonal? " << (diagonal.isDiagonal() ? "Yes" : "No") << "\n";
    printValue("trace(diagonal)", diagonal.trace());
    printValue("det(diagonal) = 5*3*(-2)", diagonal.determinant());

    std::cout << "\nIs I identity? " << (I.isIdentity() ? "Yes" : "No") << "\n";

    // ========================================
    // 7. Special Matrices
    // ========================================
    std::cout << "\n=== SPECIAL MATRICES ===\n";

    Vector diagVec({2, 3, 5, 7});
    Matrix diagMat = SpecialMatrices::diagonal(diagVec);
    printMatrix("Diagonal from vector [2,3,5,7]", diagMat);

    // ========================================
    // 8. Rotation Matrix Example
    // ========================================
    std::cout << "\n=== ROTATION MATRIX (2D) ===\n";

    double angle = M_PI / 4;  // 45 degrees
    Matrix rotation({{std::cos(angle), -std::sin(angle)},
                     {std::sin(angle), std::cos(angle)}});
    printMatrix("Rotation matrix (45°)", rotation);
    printValue("det(rotation)", rotation.determinant());  // Should be 1

    Vector point({1, 0});
    std::cout << "\nOriginal point: " << point.toString() << "\n";
    Vector rotated = rotation * point;
    std::cout << "After 45° rotation: " << rotated.toString() << "\n";

    // ========================================
    // 9. Matrix Powers (via repeated multiplication)
    // ========================================
    std::cout << "\n=== MATRIX POWERS ===\n";

    Matrix M({{1, 1},
              {1, 0}});  // Fibonacci matrix
    printMatrix("Fibonacci Matrix M", M);

    Matrix Msquared = M * M;
    printMatrix("M^2", Msquared);

    Matrix Mcubed = Msquared * M;
    printMatrix("M^3", Mcubed);

    Matrix Mfifth = Mcubed * Msquared;
    printMatrix("M^5", Mfifth);
    std::cout << "Note: M^n generates Fibonacci numbers!\n";

    // ========================================
    // 10. Properties Documentation
    // ========================================
    std::cout << "\n=== DETERMINANT PROPERTIES ===\n";
    std::cout << MatrixProperties::determinantProperties() << "\n";

    std::cout << "\n=== TRACE PROPERTIES ===\n";
    std::cout << MatrixProperties::traceProperties() << "\n";

    std::cout << "\n=== ORTHOGONAL MATRICES ===\n";
    std::cout << SpecialMatrices::orthogonalMatrices() << "\n";

    std::cout << "\n=== DECOMPOSITIONS (Documentation) ===\n";
    std::cout << MatrixDecompositions::luDecomposition() << "\n\n";
    std::cout << MatrixDecompositions::qrDecomposition() << "\n\n";
    std::cout << MatrixDecompositions::svd() << "\n";

    return 0;
}
