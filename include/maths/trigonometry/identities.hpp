#ifndef MATHS_TRIGONOMETRY_IDENTITIES_HPP
#define MATHS_TRIGONOMETRY_IDENTITIES_HPP

#include <cmath>
#include <string>
#include <vector>

/**
 * @file identities.hpp
 * @brief Trigonometric functions and identities
 *
 * Implements:
 * - Fundamental trig functions
 * - Pythagorean identities
 * - Angle sum/difference formulas
 * - Double and half-angle formulas
 * - Product-to-sum and sum-to-product
 * - Inverse trig functions
 */

namespace maths::trigonometry {

/**
 * @class TrigonometricFunctions
 * @brief Basic trigonometric functions and relationships
 */
class TrigonometricFunctions {
public:
    /**
     * @brief Convert degrees to radians
     */
    static double degreesToRadians(double degrees) {
        return degrees * M_PI / 180.0;
    }

    /**
     * @brief Convert radians to degrees
     */
    static double radiansToDegrees(double radians) {
        return radians * 180.0 / M_PI;
    }

    /**
     * @brief Six trigonometric functions
     */
    static double sine(double angle_rad) {
        return std::sin(angle_rad);
    }

    static double cosine(double angle_rad) {
        return std::cos(angle_rad);
    }

    static double tangent(double angle_rad) {
        return std::tan(angle_rad);
    }

    static double cotangent(double angle_rad) {
        return 1.0 / std::tan(angle_rad);
    }

    static double secant(double angle_rad) {
        return 1.0 / std::cos(angle_rad);
    }

    static double cosecant(double angle_rad) {
        return 1.0 / std::sin(angle_rad);
    }

    /**
     * @brief Definitions in terms of unit circle
     */
    static std::string unitCircleDefinitions() {
        return "Unit circle definitions (radius r = 1):\n"
               "\n"
               "sin(θ) = y (vertical coordinate)\n"
               "cos(θ) = x (horizontal coordinate)\n"
               "tan(θ) = y/x = sin(θ)/cos(θ)\n"
               "\n"
               "cot(θ) = x/y = cos(θ)/sin(θ)\n"
               "sec(θ) = 1/x = 1/cos(θ)\n"
               "csc(θ) = 1/y = 1/sin(θ)";
    }

    /**
     * @brief Right triangle definitions
     */
    static std::string rightTriangleDefinitions() {
        return "Right triangle definitions (angle θ):\n"
               "\n"
               "sin(θ) = opposite / hypotenuse\n"
               "cos(θ) = adjacent / hypotenuse\n"
               "tan(θ) = opposite / adjacent\n"
               "\n"
               "cot(θ) = adjacent / opposite\n"
               "sec(θ) = hypotenuse / adjacent\n"
               "csc(θ) = hypotenuse / opposite\n"
               "\n"
               "Mnemonic: SOH-CAH-TOA";
    }

    /**
     * @brief Special angles (exact values)
     */
    static std::string specialAngles() {
        return "Special angles (degrees → radians):\n"
               "\n"
               "     θ°   θ(rad)   sin     cos     tan\n"
               "    0°     0       0       1       0\n"
               "   30°    π/6      1/2     √3/2    √3/3\n"
               "   45°    π/4      √2/2    √2/2    1\n"
               "   60°    π/3      √3/2    1/2     √3\n"
               "   90°    π/2      1       0       ∞\n"
               "  180°    π        0      -1       0\n"
               "  270°    3π/2    -1       0       ∞\n"
               "  360°    2π       0       1       0";
    }
};

/**
 * @class PythagoreanIdentities
 * @brief Fundamental Pythagorean identities
 */
class PythagoreanIdentities {
public:
    /**
     * @brief Fundamental identity: sin²θ + cos²θ = 1
     */
    static std::string fundamentalIdentity() {
        return "Fundamental Pythagorean Identity:\n"
               "\n"
               "sin²(θ) + cos²(θ) = 1\n"
               "\n"
               "Proof: From unit circle x² + y² = 1\n"
               "where x = cos(θ), y = sin(θ)";
    }

    /**
     * @brief Verify sin²θ + cos²θ = 1
     */
    static bool verifyFundamental(double theta, double tolerance = 1e-10) {
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sum = sin_theta * sin_theta + cos_theta * cos_theta;
        return std::abs(sum - 1.0) < tolerance;
    }

    /**
     * @brief Tangent identity: 1 + tan²θ = sec²θ
     */
    static std::string tangentIdentity() {
        return "Tangent Pythagorean Identity:\n"
               "\n"
               "1 + tan²(θ) = sec²(θ)\n"
               "\n"
               "Proof: Divide sin²θ + cos²θ = 1 by cos²θ:\n"
               "(sin²θ/cos²θ) + 1 = 1/cos²θ\n"
               "tan²θ + 1 = sec²θ";
    }

    /**
     * @brief Verify 1 + tan²θ = sec²θ
     */
    static bool verifyTangent(double theta, double tolerance = 1e-10) {
        if (std::abs(std::cos(theta)) < 1e-10) {
            return false;  // Undefined
        }
        double tan_theta = std::tan(theta);
        double sec_theta = 1.0 / std::cos(theta);
        double lhs = 1.0 + tan_theta * tan_theta;
        double rhs = sec_theta * sec_theta;
        return std::abs(lhs - rhs) < tolerance;
    }

    /**
     * @brief Cotangent identity: 1 + cot²θ = csc²θ
     */
    static std::string cotangentIdentity() {
        return "Cotangent Pythagorean Identity:\n"
               "\n"
               "1 + cot²(θ) = csc²(θ)\n"
               "\n"
               "Proof: Divide sin²θ + cos²θ = 1 by sin²θ:\n"
               "1 + (cos²θ/sin²θ) = 1/sin²θ\n"
               "1 + cot²θ = csc²θ";
    }

    /**
     * @brief All three Pythagorean identities
     */
    static std::string allIdentities() {
        return "Three Pythagorean Identities:\n"
               "\n"
               "1. sin²(θ) + cos²(θ) = 1\n"
               "2. 1 + tan²(θ) = sec²(θ)\n"
               "3. 1 + cot²(θ) = csc²(θ)\n"
               "\n"
               "All derived from x² + y² = 1 on unit circle";
    }
};

/**
 * @class AngleFormulas
 * @brief Sum, difference, double, and half-angle formulas
 */
class AngleFormulas {
public:
    /**
     * @brief Sine of sum: sin(α + β)
     */
    static double sineSum(double alpha, double beta) {
        return std::sin(alpha) * std::cos(beta) + std::cos(alpha) * std::sin(beta);
    }

    /**
     * @brief Sine of difference: sin(α - β)
     */
    static double sineDifference(double alpha, double beta) {
        return std::sin(alpha) * std::cos(beta) - std::cos(alpha) * std::sin(beta);
    }

    /**
     * @brief Cosine of sum: cos(α + β)
     */
    static double cosineSum(double alpha, double beta) {
        return std::cos(alpha) * std::cos(beta) - std::sin(alpha) * std::sin(beta);
    }

    /**
     * @brief Cosine of difference: cos(α - β)
     */
    static double cosineDifference(double alpha, double beta) {
        return std::cos(alpha) * std::cos(beta) + std::sin(alpha) * std::sin(beta);
    }

    /**
     * @brief Tangent of sum: tan(α + β)
     */
    static double tangentSum(double alpha, double beta) {
        return (std::tan(alpha) + std::tan(beta)) /
               (1.0 - std::tan(alpha) * std::tan(beta));
    }

    /**
     * @brief Tangent of difference: tan(α - β)
     */
    static double tangentDifference(double alpha, double beta) {
        return (std::tan(alpha) - std::tan(beta)) /
               (1.0 + std::tan(alpha) * std::tan(beta));
    }

    /**
     * @brief Angle sum/difference formulas
     */
    static std::string sumDifferenceFormulas() {
        return "Angle Sum and Difference Formulas:\n"
               "\n"
               "sin(α ± β) = sin(α)cos(β) ± cos(α)sin(β)\n"
               "cos(α ± β) = cos(α)cos(β) ∓ sin(α)sin(β)\n"
               "tan(α ± β) = [tan(α) ± tan(β)] / [1 ∓ tan(α)tan(β)]\n"
               "\n"
               "Note: Signs are opposite for cosine!";
    }

    /**
     * @brief Double angle formulas
     */
    static std::string doubleAngleFormulas() {
        return "Double Angle Formulas (θ = 2α):\n"
               "\n"
               "sin(2α) = 2sin(α)cos(α)\n"
               "\n"
               "cos(2α) = cos²(α) - sin²(α)\n"
               "        = 2cos²(α) - 1\n"
               "        = 1 - 2sin²(α)\n"
               "\n"
               "tan(2α) = 2tan(α) / [1 - tan²(α)]";
    }

    /**
     * @brief Half angle formulas
     */
    static std::string halfAngleFormulas() {
        return "Half Angle Formulas:\n"
               "\n"
               "sin(α/2) = ±√[(1 - cos(α))/2]\n"
               "cos(α/2) = ±√[(1 + cos(α))/2]\n"
               "tan(α/2) = ±√[(1 - cos(α))/(1 + cos(α))]\n"
               "         = sin(α) / [1 + cos(α)]\n"
               "         = [1 - cos(α)] / sin(α)\n"
               "\n"
               "Sign depends on quadrant of α/2";
    }

    /**
     * @brief Triple angle formulas
     */
    static std::string tripleAngleFormulas() {
        return "Triple Angle Formulas:\n"
               "\n"
               "sin(3α) = 3sin(α) - 4sin³(α)\n"
               "cos(3α) = 4cos³(α) - 3cos(α)\n"
               "tan(3α) = [3tan(α) - tan³(α)] / [1 - 3tan²(α)]";
    }
};

/**
 * @class ProductSumFormulas
 * @brief Product-to-sum and sum-to-product identities
 */
class ProductSumFormulas {
public:
    /**
     * @brief Product-to-sum formulas
     */
    static std::string productToSum() {
        return "Product-to-Sum Formulas:\n"
               "\n"
               "sin(α)sin(β) = [cos(α-β) - cos(α+β)] / 2\n"
               "cos(α)cos(β) = [cos(α-β) + cos(α+β)] / 2\n"
               "sin(α)cos(β) = [sin(α+β) + sin(α-β)] / 2\n"
               "cos(α)sin(β) = [sin(α+β) - sin(α-β)] / 2\n"
               "\n"
               "Useful for integration and signal processing";
    }

    /**
     * @brief Sum-to-product formulas
     */
    static std::string sumToProduct() {
        return "Sum-to-Product Formulas:\n"
               "\n"
               "sin(α) + sin(β) = 2sin[(α+β)/2]cos[(α-β)/2]\n"
               "sin(α) - sin(β) = 2cos[(α+β)/2]sin[(α-β)/2]\n"
               "cos(α) + cos(β) = 2cos[(α+β)/2]cos[(α-β)/2]\n"
               "cos(α) - cos(β) = -2sin[(α+β)/2]sin[(α-β)/2]\n"
               "\n"
               "Useful for solving trig equations";
    }

    /**
     * @brief Power reduction formulas
     */
    static std::string powerReduction() {
        return "Power Reduction Formulas:\n"
               "\n"
               "sin²(α) = [1 - cos(2α)] / 2\n"
               "cos²(α) = [1 + cos(2α)] / 2\n"
               "tan²(α) = [1 - cos(2α)] / [1 + cos(2α)]\n"
               "\n"
               "Derived from double angle formulas\n"
               "Useful for integration of trig powers";
    }
};

/**
 * @class InverseTrigFunctions
 * @brief Inverse trigonometric functions
 */
class InverseTrigFunctions {
public:
    /**
     * @brief Arcsine (inverse sine)
     *
     * @param x Value in [-1, 1]
     * @return Angle in [-π/2, π/2]
     */
    static double arcsin(double x) {
        if (x < -1.0 || x > 1.0) {
            throw std::invalid_argument("arcsin domain: [-1, 1]");
        }
        return std::asin(x);
    }

    /**
     * @brief Arccosine (inverse cosine)
     *
     * @param x Value in [-1, 1]
     * @return Angle in [0, π]
     */
    static double arccos(double x) {
        if (x < -1.0 || x > 1.0) {
            throw std::invalid_argument("arccos domain: [-1, 1]");
        }
        return std::acos(x);
    }

    /**
     * @brief Arctangent (inverse tangent)
     *
     * @param x Any real value
     * @return Angle in (-π/2, π/2)
     */
    static double arctan(double x) {
        return std::atan(x);
    }

    /**
     * @brief Two-argument arctangent (atan2)
     *
     * @param y Numerator
     * @param x Denominator
     * @return Angle in (-π, π]
     */
    static double arctan2(double y, double x) {
        return std::atan2(y, x);
    }

    /**
     * @brief Domains and ranges
     */
    static std::string domainsAndRanges() {
        return "Inverse Trig Functions:\n"
               "\n"
               "Function    Domain      Range\n"
               "arcsin(x)   [-1, 1]     [-π/2, π/2]\n"
               "arccos(x)   [-1, 1]     [0, π]\n"
               "arctan(x)   (-∞, ∞)     (-π/2, π/2)\n"
               "arccot(x)   (-∞, ∞)     (0, π)\n"
               "arcsec(x)   (-∞,-1]∪[1,∞) [0,π], x≠π/2\n"
               "arccsc(x)   (-∞,-1]∪[1,∞) [-π/2,π/2], x≠0";
    }

    /**
     * @brief Inverse trig identities
     */
    static std::string inverseIdentities() {
        return "Inverse Trig Identities:\n"
               "\n"
               "arcsin(x) + arccos(x) = π/2\n"
               "arctan(x) + arccot(x) = π/2\n"
               "\n"
               "arcsin(-x) = -arcsin(x)\n"
               "arccos(-x) = π - arccos(x)\n"
               "arctan(-x) = -arctan(x)\n"
               "\n"
               "arctan(1/x) = π/2 - arctan(x)  (x > 0)\n"
               "arctan(1/x) = -π/2 - arctan(x) (x < 0)";
    }

    /**
     * @brief Derivatives of inverse trig functions
     */
    static std::string derivatives() {
        return "Derivatives of Inverse Trig Functions:\n"
               "\n"
               "d/dx [arcsin(x)] = 1/√(1-x²)\n"
               "d/dx [arccos(x)] = -1/√(1-x²)\n"
               "d/dx [arctan(x)] = 1/(1+x²)\n"
               "d/dx [arccot(x)] = -1/(1+x²)\n"
               "d/dx [arcsec(x)] = 1/(|x|√(x²-1))\n"
               "d/dx [arccsc(x)] = -1/(|x|√(x²-1))";
    }
};

/**
 * @class LawsOfTriangles
 * @brief Law of sines and law of cosines
 */
class LawsOfTriangles {
public:
    /**
     * @brief Law of Sines
     *
     * a/sin(A) = b/sin(B) = c/sin(C) = 2R
     *
     * where R is circumradius
     */
    static std::string lawOfSines() {
        return "Law of Sines:\n"
               "\n"
               "a/sin(A) = b/sin(B) = c/sin(C) = 2R\n"
               "\n"
               "where a, b, c are sides opposite angles A, B, C\n"
               "      R is circumradius of triangle\n"
               "\n"
               "Use when: Given 2 angles + 1 side, or 2 sides + non-included angle";
    }

    /**
     * @brief Solve triangle using law of sines
     *
     * Given angle A, angle B, and side a, find side b
     */
    static double solveWithSines(double angleA_rad, double angleB_rad, double sideA) {
        // b = a * sin(B) / sin(A)
        return sideA * std::sin(angleB_rad) / std::sin(angleA_rad);
    }

    /**
     * @brief Law of Cosines
     *
     * c² = a² + b² - 2ab·cos(C)
     *
     * Generalizes Pythagorean theorem
     */
    static std::string lawOfCosines() {
        return "Law of Cosines:\n"
               "\n"
               "c² = a² + b² - 2ab·cos(C)\n"
               "b² = a² + c² - 2ac·cos(B)\n"
               "a² = b² + c² - 2bc·cos(A)\n"
               "\n"
               "Generalizes Pythagorean theorem (when C = 90°, cos(C) = 0)\n"
               "\n"
               "Use when: Given 3 sides, or 2 sides + included angle";
    }

    /**
     * @brief Solve for side using law of cosines
     *
     * Given sides a, b and included angle C, find side c
     */
    static double solveWithCosines(double sideA, double sideB, double angleC_rad) {
        // c² = a² + b² - 2ab·cos(C)
        double c_squared = sideA * sideA + sideB * sideB -
                          2.0 * sideA * sideB * std::cos(angleC_rad);
        return std::sqrt(c_squared);
    }

    /**
     * @brief Find angle using law of cosines
     *
     * Given all three sides a, b, c, find angle C
     */
    static double findAngleWithCosines(double sideA, double sideB, double sideC) {
        // cos(C) = (a² + b² - c²) / (2ab)
        double cos_C = (sideA * sideA + sideB * sideB - sideC * sideC) /
                       (2.0 * sideA * sideB);
        return std::acos(cos_C);
    }

    /**
     * @brief Area formulas
     */
    static std::string areaFormulas() {
        return "Triangle Area Formulas:\n"
               "\n"
               "Given base b and height h:\n"
               "  A = (1/2)bh\n"
               "\n"
               "Given two sides and included angle:\n"
               "  A = (1/2)ab·sin(C)\n"
               "\n"
               "Heron's formula (three sides a, b, c):\n"
               "  s = (a+b+c)/2 (semiperimeter)\n"
               "  A = √[s(s-a)(s-b)(s-c)]\n"
               "\n"
               "Using circumradius R:\n"
               "  A = abc/(4R)";
    }
};

} // namespace maths::trigonometry

#endif // MATHS_TRIGONOMETRY_IDENTITIES_HPP
