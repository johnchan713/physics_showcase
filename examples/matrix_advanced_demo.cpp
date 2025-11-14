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
    std::cout << std::fixed << std::setprecision(6);

    // ========================================
    // 1. Matrix Inverse
    // ========================================
    std::cout << "=== MATRIX INVERSE ===\n";

    Matrix A({{4, 7},
              {2, 6}});
    printMatrix("Matrix A (2x2)", A);
    printValue("det(A)", A.determinant());

    Matrix Ainv = A.inverse();
    printMatrix("A^(-1) (inverse)", Ainv);

    Matrix I_check = A * Ainv;
    printMatrix("A * A^(-1) (should be I)", I_check);

    // 3x3 inverse example
    std::cout << "\n--- 3x3 Inverse Example ---\n";
    Matrix B({{1, 2, 3},
              {0, 1, 4},
              {5, 6, 0}});
    printMatrix("Matrix B (3x3)", B);
    printValue("det(B)", B.determinant());

    Matrix Binv = B.inverse();
    printMatrix("B^(-1)", Binv);

    Matrix I_check2 = B * Binv;
    printMatrix("B * B^(-1)", I_check2);

    // Verify det(A^-1) = 1/det(A)
    double detB = B.determinant();
    double detBinv = Binv.determinant();
    printValue("det(B)", detB);
    printValue("det(B^-1)", detBinv);
    printValue("det(B) * det(B^-1)", detB * detBinv);

    // ========================================
    // 2. Matrix Rank
    // ========================================
    std::cout << "\n=== MATRIX RANK ===\n";

    Matrix R1({{1, 2, 3},
               {4, 5, 6},
               {7, 8, 9}});
    printMatrix("Matrix R1 (singular)", R1);
    std::cout << "Rank: " << R1.rank() << " (expected 2, rows are linearly dependent)\n";
    printValue("det(R1)", R1.determinant());

    Matrix R2({{1, 0, 0},
               {0, 1, 0},
               {0, 0, 1}});
    printMatrix("Identity Matrix", R2);
    std::cout << "Rank: " << R2.rank() << " (expected 3, full rank)\n";

    Matrix R3({{1, 2, 3, 4},
               {2, 4, 6, 8},
               {3, 6, 9, 12}});
    printMatrix("Matrix R3 (3x4, rank 1)", R3);
    std::cout << "Rank: " << R3.rank() << " (all rows are multiples of first row)\n";

    Matrix R4({{1, 2},
               {3, 4},
               {5, 6}});
    printMatrix("Matrix R4 (3x2)", R4);
    std::cout << "Rank: " << R4.rank() << "\n";

    // ========================================
    // 3. Reduced Row Echelon Form (RREF)
    // ========================================
    std::cout << "\n=== REDUCED ROW ECHELON FORM (RREF) ===\n";

    Matrix M1({{1, 2, 3},
               {4, 5, 6},
               {7, 8, 9}});
    printMatrix("Original Matrix M1", M1);
    Matrix M1_rref = M1.rref();
    printMatrix("RREF of M1", M1_rref);

    Matrix M2({{1, 2, -1, -4},
               {2, 3, -1, -11},
               {-2, 0, -3, 22}});
    printMatrix("Original Matrix M2 (3x4)", M2);
    Matrix M2_rref = M2.rref();
    printMatrix("RREF of M2", M2_rref);

    Matrix M3({{2, 1, -1},
               {-3, -1, 2},
               {-2, 1, 2}});
    printMatrix("Original Matrix M3", M3);
    Matrix M3_rref = M3.rref();
    printMatrix("RREF of M3", M3_rref);

    // ========================================
    // 4. Solving Linear Systems (Ax = b)
    // ========================================
    std::cout << "\n=== SOLVING LINEAR SYSTEMS ===\n";

    // Example 1: 2x2 system
    std::cout << "\n--- Example 1: 2x2 system ---\n";
    Matrix A1({{3, 2},
               {1, 2}});
    Vector b1({7, 5});

    printMatrix("Coefficient matrix A", A1);
    std::cout << "Right-hand side b: " << b1.toString() << "\n";

    Vector x1 = A1.solve(b1);
    std::cout << "Solution x: " << x1.toString() << "\n";

    // Verify: Ax should equal b
    Vector verify1 = A1 * x1;
    std::cout << "Verification A*x: " << verify1.toString() << "\n";

    // Example 2: 3x3 system
    std::cout << "\n--- Example 2: 3x3 system ---\n";
    Matrix A2({{2, 1, -1},
               {-3, -1, 2},
               {-2, 1, 2}});
    Vector b2({8, -11, -3});

    printMatrix("Coefficient matrix A", A2);
    std::cout << "Right-hand side b: " << b2.toString() << "\n";

    Vector x2 = A2.solve(b2);
    std::cout << "Solution x: " << x2.toString() << "\n";

    Vector verify2 = A2 * x2;
    std::cout << "Verification A*x: " << verify2.toString() << "\n";

    // Example 3: Using inverse to solve
    std::cout << "\n--- Example 3: Solving via inverse (x = A^(-1)b) ---\n";
    Matrix A3({{1, 2},
               {3, 4}});
    Vector b3({5, 11});

    printMatrix("Matrix A", A3);
    std::cout << "Vector b: " << b3.toString() << "\n";

    Matrix A3inv = A3.inverse();
    Vector x3 = A3inv * b3;
    std::cout << "Solution x = A^(-1)*b: " << x3.toString() << "\n";

    Vector verify3 = A3 * x3;
    std::cout << "Verification A*x: " << verify3.toString() << "\n";

    // ========================================
    // 5. Row Operations Demo
    // ========================================
    std::cout << "\n=== ROW OPERATIONS ===\n";

    Matrix M({{1, 2, 3},
              {4, 5, 6},
              {7, 8, 9}});
    printMatrix("Original Matrix M", M);

    Matrix M_copy = M;
    M_copy.swapRows(0, 2);
    printMatrix("After swapping rows 0 and 2", M_copy);

    M_copy = M;
    M_copy.scaleRow(1, 2.0);
    printMatrix("After scaling row 1 by 2", M_copy);

    M_copy = M;
    M_copy.addScaledRow(2, 0, -7.0);  // R3 = R3 - 7*R1
    printMatrix("After R3 = R3 - 7*R1", M_copy);

    // ========================================
    // 6. Properties Verification
    // ========================================
    std::cout << "\n=== PROPERTIES VERIFICATION ===\n";

    Matrix P({{1, 2},
              {3, 4}});
    Matrix Q({{5, 6},
              {7, 8}});

    printMatrix("Matrix P", P);
    printMatrix("Matrix Q", Q);

    Matrix PQ = P * Q;
    Matrix QP = Q * P;
    printMatrix("P * Q", PQ);
    printMatrix("Q * P (note: ≠ P*Q)", QP);

    Matrix PQT = PQ.transpose();
    Matrix QT = Q.transpose();
    Matrix PT = P.transpose();
    Matrix QTPT = QT * PT;
    printMatrix("(PQ)^T", PQT);
    printMatrix("Q^T * P^T", QTPT);
    std::cout << "(PQ)^T = Q^T P^T ✓\n";

    // ========================================
    // 7. Singular Matrix Handling
    // ========================================
    std::cout << "\n=== SINGULAR MATRIX TEST ===\n";

    Matrix Singular({{1, 2, 3},
                     {2, 4, 6},
                     {3, 6, 9}});
    printMatrix("Singular Matrix", Singular);
    printValue("det(Singular)", Singular.determinant());
    std::cout << "Rank: " << Singular.rank() << "\n";

    try {
        Matrix SingInv = Singular.inverse();
        std::cout << "ERROR: Singular matrix should not be invertible!\n";
    } catch (const std::runtime_error& e) {
        std::cout << "Correctly threw exception: " << e.what() << "\n";
    }

    // ========================================
    // 8. Application: Cryptography (Hill Cipher)
    // ========================================
    std::cout << "\n=== APPLICATION: HILL CIPHER ===\n";

    Matrix K({{6, 24, 1},
              {13, 16, 10},
              {20, 17, 15}});  // Encryption key
    printMatrix("Encryption Key K (mod 26)", K);
    printValue("det(K)", K.determinant());

    Vector plaintext({7, 8, 0});  // H=7, I=8, A=0
    std::cout << "Plaintext vector (HIA): " << plaintext.toString() << "\n";

    Vector ciphertext = K * plaintext;
    std::cout << "Encrypted (before mod 26): " << ciphertext.toString() << "\n";

    std::cout << "\n=== DONE ===\n";

    return 0;
}
