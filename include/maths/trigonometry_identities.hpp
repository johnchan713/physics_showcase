#ifndef MATHS_TRIGONOMETRY_IDENTITIES_HPP
#define MATHS_TRIGONOMETRY_IDENTITIES_HPP

#include <cmath>
#include <stdexcept>

/**
 * @file identities.hpp
 * @brief Trigonometric functions and computational identities
 *
 * Implements computational trigonometric functions following physics module pattern.
 */

namespace maths::trigonometry {

/**
 * @class TrigonometricFunctions
 * @brief Basic trigonometric functions and conversions
 */
class TrigonometricFunctions {
public:
    /**
     * @brief Convert degrees to radians
     * @param degrees Angle in degrees
     * @return Angle in radians
     */
    static double degreesToRadians(double degrees) {
        return degrees * M_PI / 180.0;
    }

    /**
     * @brief Convert radians to degrees
     * @param radians Angle in radians
     * @return Angle in degrees
     */
    static double radiansToDegrees(double radians) {
        return radians * 180.0 / M_PI;
    }

    /**
     * @brief Sine function
     * @param angle_rad Angle in radians
     * @return sin(angle)
     */
    static double sine(double angle_rad) {
        return std::sin(angle_rad);
    }

    /**
     * @brief Cosine function
     * @param angle_rad Angle in radians
     * @return cos(angle)
     */
    static double cosine(double angle_rad) {
        return std::cos(angle_rad);
    }

    /**
     * @brief Tangent function
     * @param angle_rad Angle in radians
     * @return tan(angle)
     */
    static double tangent(double angle_rad) {
        return std::tan(angle_rad);
    }

    /**
     * @brief Cotangent function
     * @param angle_rad Angle in radians
     * @return cot(angle) = 1/tan(angle)
     */
    static double cotangent(double angle_rad) {
        return 1.0 / std::tan(angle_rad);
    }

    /**
     * @brief Secant function
     * @param angle_rad Angle in radians
     * @return sec(angle) = 1/cos(angle)
     */
    static double secant(double angle_rad) {
        return 1.0 / std::cos(angle_rad);
    }

    /**
     * @brief Cosecant function
     * @param angle_rad Angle in radians
     * @return csc(angle) = 1/sin(angle)
     */
    static double cosecant(double angle_rad) {
        return 1.0 / std::sin(angle_rad);
    }
};

/**
 * @class PythagoreanIdentities
 * @brief Computational verification of Pythagorean identities
 */
class PythagoreanIdentities {
public:
    /**
     * @brief Verify sin²θ + cos²θ = 1
     * @param theta Angle in radians
     * @param tolerance Acceptable error (default: 1e-10)
     * @return true if identity holds within tolerance
     */
    static bool verifyFundamental(double theta, double tolerance = 1e-10) {
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sum = sin_theta * sin_theta + cos_theta * cos_theta;
        return std::abs(sum - 1.0) < tolerance;
    }

    /**
     * @brief Verify 1 + tan²θ = sec²θ
     * @param theta Angle in radians
     * @param tolerance Acceptable error (default: 1e-10)
     * @return true if identity holds within tolerance
     */
    static bool verifyTangent(double theta, double tolerance = 1e-10) {
        if (std::abs(std::cos(theta)) < 1e-10) {
            return false;  // Undefined at cos(θ) = 0
        }
        double tan_theta = std::tan(theta);
        double sec_theta = 1.0 / std::cos(theta);
        double lhs = 1.0 + tan_theta * tan_theta;
        double rhs = sec_theta * sec_theta;
        return std::abs(lhs - rhs) < tolerance;
    }

    /**
     * @brief Verify 1 + cot²θ = csc²θ
     * @param theta Angle in radians
     * @param tolerance Acceptable error (default: 1e-10)
     * @return true if identity holds within tolerance
     */
    static bool verifyCotangent(double theta, double tolerance = 1e-10) {
        if (std::abs(std::sin(theta)) < 1e-10) {
            return false;  // Undefined at sin(θ) = 0
        }
        double cot_theta = 1.0 / std::tan(theta);
        double csc_theta = 1.0 / std::sin(theta);
        double lhs = 1.0 + cot_theta * cot_theta;
        double rhs = csc_theta * csc_theta;
        return std::abs(lhs - rhs) < tolerance;
    }
};

/**
 * @class AngleFormulas
 * @brief Sum, difference, double, and half-angle computations
 */
class AngleFormulas {
public:
    /**
     * @brief Compute sin(α + β)
     * @param alpha First angle in radians
     * @param beta Second angle in radians
     * @return sin(α + β)
     */
    static double sineSum(double alpha, double beta) {
        return std::sin(alpha) * std::cos(beta) + std::cos(alpha) * std::sin(beta);
    }

    /**
     * @brief Compute sin(α - β)
     * @param alpha First angle in radians
     * @param beta Second angle in radians
     * @return sin(α - β)
     */
    static double sineDifference(double alpha, double beta) {
        return std::sin(alpha) * std::cos(beta) - std::cos(alpha) * std::sin(beta);
    }

    /**
     * @brief Compute cos(α + β)
     * @param alpha First angle in radians
     * @param beta Second angle in radians
     * @return cos(α + β)
     */
    static double cosineSum(double alpha, double beta) {
        return std::cos(alpha) * std::cos(beta) - std::sin(alpha) * std::sin(beta);
    }

    /**
     * @brief Compute cos(α - β)
     * @param alpha First angle in radians
     * @param beta Second angle in radians
     * @return cos(α - β)
     */
    static double cosineDifference(double alpha, double beta) {
        return std::cos(alpha) * std::cos(beta) + std::sin(alpha) * std::sin(beta);
    }

    /**
     * @brief Compute tan(α + β)
     * @param alpha First angle in radians
     * @param beta Second angle in radians
     * @return tan(α + β)
     */
    static double tangentSum(double alpha, double beta) {
        return (std::tan(alpha) + std::tan(beta)) /
               (1.0 - std::tan(alpha) * std::tan(beta));
    }

    /**
     * @brief Compute tan(α - β)
     * @param alpha First angle in radians
     * @param beta Second angle in radians
     * @return tan(α - β)
     */
    static double tangentDifference(double alpha, double beta) {
        return (std::tan(alpha) - std::tan(beta)) /
               (1.0 + std::tan(alpha) * std::tan(beta));
    }

    /**
     * @brief Compute sin(2α)
     * @param alpha Angle in radians
     * @return sin(2α) = 2sin(α)cos(α)
     */
    static double sineDouble(double alpha) {
        return 2.0 * std::sin(alpha) * std::cos(alpha);
    }

    /**
     * @brief Compute cos(2α) using cos²(α) - sin²(α)
     * @param alpha Angle in radians
     * @return cos(2α)
     */
    static double cosineDouble(double alpha) {
        double cos_alpha = std::cos(alpha);
        double sin_alpha = std::sin(alpha);
        return cos_alpha * cos_alpha - sin_alpha * sin_alpha;
    }

    /**
     * @brief Compute tan(2α)
     * @param alpha Angle in radians
     * @return tan(2α) = 2tan(α)/(1-tan²(α))
     */
    static double tangentDouble(double alpha) {
        double tan_alpha = std::tan(alpha);
        return (2.0 * tan_alpha) / (1.0 - tan_alpha * tan_alpha);
    }

    /**
     * @brief Compute sin(α/2) (positive root)
     * @param alpha Angle in radians
     * @return sin(α/2) = √[(1-cos(α))/2]
     */
    static double sineHalf(double alpha) {
        return std::sqrt((1.0 - std::cos(alpha)) / 2.0);
    }

    /**
     * @brief Compute cos(α/2) (positive root)
     * @param alpha Angle in radians
     * @return cos(α/2) = √[(1+cos(α))/2]
     */
    static double cosineHalf(double alpha) {
        return std::sqrt((1.0 + std::cos(alpha)) / 2.0);
    }

    /**
     * @brief Compute tan(α/2) using sin/(1+cos) formula
     * @param alpha Angle in radians
     * @return tan(α/2) = sin(α)/(1+cos(α))
     */
    static double tangentHalf(double alpha) {
        return std::sin(alpha) / (1.0 + std::cos(alpha));
    }
};

/**
 * @class ProductSumFormulas
 * @brief Product-to-sum computations
 */
class ProductSumFormulas {
public:
    /**
     * @brief Compute sin(α)sin(β) as sum of cosines
     * @return [cos(α-β) - cos(α+β)]/2
     */
    static double sineSineProduct(double alpha, double beta) {
        return (std::cos(alpha - beta) - std::cos(alpha + beta)) / 2.0;
    }

    /**
     * @brief Compute cos(α)cos(β) as sum of cosines
     * @return [cos(α-β) + cos(α+β)]/2
     */
    static double cosineCosineProduct(double alpha, double beta) {
        return (std::cos(alpha - beta) + std::cos(alpha + beta)) / 2.0;
    }

    /**
     * @brief Compute sin(α)cos(β) as sum of sines
     * @return [sin(α+β) + sin(α-β)]/2
     */
    static double sineCosineProduct(double alpha, double beta) {
        return (std::sin(alpha + beta) + std::sin(alpha - beta)) / 2.0;
    }
};

/**
 * @class InverseTrigFunctions
 * @brief Inverse trigonometric functions with domain checking
 */
class InverseTrigFunctions {
public:
    /**
     * @brief Arcsine (inverse sine)
     * @param x Value in [-1, 1]
     * @return Angle in [-π/2, π/2]
     * @throws std::invalid_argument if x not in [-1, 1]
     */
    static double arcsin(double x) {
        if (x < -1.0 || x > 1.0) {
            throw std::invalid_argument("arcsin domain: [-1, 1]");
        }
        return std::asin(x);
    }

    /**
     * @brief Arccosine (inverse cosine)
     * @param x Value in [-1, 1]
     * @return Angle in [0, π]
     * @throws std::invalid_argument if x not in [-1, 1]
     */
    static double arccos(double x) {
        if (x < -1.0 || x > 1.0) {
            throw std::invalid_argument("arccos domain: [-1, 1]");
        }
        return std::acos(x);
    }

    /**
     * @brief Arctangent (inverse tangent)
     * @param x Any real value
     * @return Angle in (-π/2, π/2)
     */
    static double arctan(double x) {
        return std::atan(x);
    }

    /**
     * @brief Two-argument arctangent
     * @param y Numerator
     * @param x Denominator
     * @return Angle in (-π, π]
     */
    static double arctan2(double y, double x) {
        return std::atan2(y, x);
    }

    /**
     * @brief Verify arcsin(x) + arccos(x) = π/2
     * @param x Value in [-1, 1]
     * @param tolerance Acceptable error (default: 1e-10)
     * @return true if identity holds
     */
    static bool verifyArcsinArccos(double x, double tolerance = 1e-10) {
        if (x < -1.0 || x > 1.0) return false;
        double sum = std::asin(x) + std::acos(x);
        return std::abs(sum - M_PI / 2.0) < tolerance;
    }
};

/**
 * @class LawsOfTriangles
 * @brief Computational triangle solving
 */
class LawsOfTriangles {
public:
    /**
     * @brief Solve for side using law of sines
     *
     * Given angle A, angle B, and side a, compute side b
     * Formula: b = a·sin(B)/sin(A)
     *
     * @param angleA_rad Angle A in radians
     * @param angleB_rad Angle B in radians
     * @param sideA Length of side a (opposite angle A)
     * @return Length of side b (opposite angle B)
     */
    static double solveWithSines(double angleA_rad, double angleB_rad, double sideA) {
        return sideA * std::sin(angleB_rad) / std::sin(angleA_rad);
    }

    /**
     * @brief Solve for side using law of cosines
     *
     * Given sides a, b and included angle C, compute side c
     * Formula: c² = a² + b² - 2ab·cos(C)
     *
     * @param sideA Length of side a
     * @param sideB Length of side b
     * @param angleC_rad Included angle C in radians
     * @return Length of side c
     */
    static double solveWithCosines(double sideA, double sideB, double angleC_rad) {
        double c_squared = sideA * sideA + sideB * sideB -
                          2.0 * sideA * sideB * std::cos(angleC_rad);
        return std::sqrt(c_squared);
    }

    /**
     * @brief Find angle using law of cosines
     *
     * Given all three sides a, b, c, compute angle C
     * Formula: cos(C) = (a² + b² - c²)/(2ab)
     *
     * @param sideA Length of side a
     * @param sideB Length of side b
     * @param sideC Length of side c (opposite angle C)
     * @return Angle C in radians
     */
    static double findAngleWithCosines(double sideA, double sideB, double sideC) {
        double cos_C = (sideA * sideA + sideB * sideB - sideC * sideC) /
                       (2.0 * sideA * sideB);
        return std::acos(cos_C);
    }

    /**
     * @brief Compute triangle area using two sides and included angle
     *
     * Formula: A = (1/2)ab·sin(C)
     *
     * @param sideA Length of side a
     * @param sideB Length of side b
     * @param angleC_rad Included angle C in radians
     * @return Area of triangle
     */
    static double areaTwoSidesAngle(double sideA, double sideB, double angleC_rad) {
        return 0.5 * sideA * sideB * std::sin(angleC_rad);
    }

    /**
     * @brief Compute triangle area using Heron's formula
     *
     * Given three sides a, b, c:
     * s = (a+b+c)/2
     * A = √[s(s-a)(s-b)(s-c)]
     *
     * @param sideA Length of side a
     * @param sideB Length of side b
     * @param sideC Length of side c
     * @return Area of triangle
     */
    static double areaHeron(double sideA, double sideB, double sideC) {
        double s = (sideA + sideB + sideC) / 2.0;  // Semiperimeter
        return std::sqrt(s * (s - sideA) * (s - sideB) * (s - sideC));
    }
};

} // namespace maths::trigonometry

#endif // MATHS_TRIGONOMETRY_IDENTITIES_HPP
