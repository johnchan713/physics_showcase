#ifndef MATHS_TRANSFORMS_POLAR_HPP
#define MATHS_TRANSFORMS_POLAR_HPP

#include <cmath>
#include <complex>
#include <string>
#include <stdexcept>

/**
 * @file polar.hpp
 * @brief Polar coordinates and transformations
 *
 * Implements:
 * - Cartesian ↔ Polar coordinate transformations
 * - Polar form of complex numbers
 * - Jacobian for polar integrals
 * - Applications in calculus and physics
 */

namespace maths::transforms {

/**
 * @class PolarCoordinates
 * @brief 2D Polar coordinate system (r, θ)
 */
class PolarCoordinates {
public:
    /**
     * @brief Convert Cartesian (x, y) to Polar (r, θ)
     *
     * r = √(x² + y²)
     * θ = atan2(y, x)
     *
     * @param x Cartesian x-coordinate
     * @param y Cartesian y-coordinate
     * @return {r, theta} where theta in [-π, π]
     */
    static std::pair<double, double> fromCartesian(double x, double y) {
        double r = std::sqrt(x * x + y * y);
        double theta = std::atan2(y, x);  // Returns angle in [-π, π]
        return {r, theta};
    }

    /**
     * @brief Convert Polar (r, θ) to Cartesian (x, y)
     *
     * x = r cos(θ)
     * y = r sin(θ)
     *
     * @param r Radial distance (r ≥ 0)
     * @param theta Angle in radians
     * @return {x, y}
     */
    static std::pair<double, double> toCartesian(double r, double theta) {
        if (r < 0.0) {
            throw std::invalid_argument("Radius must be non-negative");
        }

        double x = r * std::cos(theta);
        double y = r * std::sin(theta);
        return {x, y};
    }

    /**
     * @brief Definition and conventions
     */
    static std::string definition() {
        return "Polar Coordinates (r, θ):\n"
               "\n"
               "r: radial distance from origin (r ≥ 0)\n"
               "θ: angle from positive x-axis (counterclockwise)\n"
               "\n"
               "Transformation:\n"
               "Cartesian → Polar:\n"
               "  r = √(x² + y²)\n"
               "  θ = atan2(y, x)\n"
               "\n"
               "Polar → Cartesian:\n"
               "  x = r cos(θ)\n"
               "  y = r sin(θ)\n"
               "\n"
               "Angle range: typically θ ∈ [-π, π] or [0, 2π)";
    }

    /**
     * @brief Distance between two points in polar coordinates
     *
     * d² = r₁² + r₂² - 2r₁r₂ cos(θ₂ - θ₁)
     *
     * Law of cosines!
     */
    static double distance(double r1, double theta1, double r2, double theta2) {
        if (r1 < 0.0 || r2 < 0.0) {
            throw std::invalid_argument("Radii must be non-negative");
        }

        double d_squared = r1*r1 + r2*r2 - 2.0*r1*r2*std::cos(theta2 - theta1);
        return std::sqrt(d_squared);
    }

    /**
     * @brief Area element in polar coordinates
     *
     * dA = r dr dθ
     *
     * Jacobian factor: r
     */
    static std::string areaElement() {
        return "Area Element in Polar Coordinates:\n"
               "\n"
               "dA = r dr dθ\n"
               "\n"
               "Jacobian:\n"
               "J = |∂(x,y)/∂(r,θ)| = r\n"
               "\n"
               "Double integral:\n"
               "∬_R f(x,y) dx dy = ∬_R f(r cos θ, r sin θ) r dr dθ\n"
               "\n"
               "Example: Circle of radius R\n"
               "Area = ∫₀^(2π) ∫₀^R r dr dθ = πR²";
    }

    /**
     * @brief Jacobian determinant
     *
     * J = |∂x/∂r  ∂x/∂θ|   |cos θ  -r sin θ|
     *     |∂y/∂r  ∂y/∂θ| = |sin θ   r cos θ| = r
     */
    static double jacobian(double r) {
        return r;
    }
};

/**
 * @class CylindricalCoordinates
 * @brief 3D Cylindrical coordinate system (r, θ, z)
 */
class CylindricalCoordinates {
public:
    /**
     * @brief Convert Cartesian (x, y, z) to Cylindrical (r, θ, z)
     */
    static std::tuple<double, double, double> fromCartesian(double x, double y, double z) {
        double r = std::sqrt(x * x + y * y);
        double theta = std::atan2(y, x);
        return {r, theta, z};  // z unchanged
    }

    /**
     * @brief Convert Cylindrical (r, θ, z) to Cartesian (x, y, z)
     */
    static std::tuple<double, double, double> toCartesian(double r, double theta, double z) {
        if (r < 0.0) {
            throw std::invalid_argument("Radius must be non-negative");
        }

        double x = r * std::cos(theta);
        double y = r * std::sin(theta);
        return {x, y, z};
    }

    /**
     * @brief Volume element in cylindrical coordinates
     *
     * dV = r dr dθ dz
     */
    static std::string volumeElement() {
        return "Volume Element in Cylindrical Coordinates:\n"
               "\n"
               "dV = r dr dθ dz\n"
               "\n"
               "Jacobian:\n"
               "J = |∂(x,y,z)/∂(r,θ,z)| = r\n"
               "\n"
               "Triple integral:\n"
               "∭_V f(x,y,z) dx dy dz = ∭_V f(r,θ,z) r dr dθ dz\n"
               "\n"
               "Example: Cylinder of radius R, height h\n"
               "Volume = ∫₀^(2π) ∫₀^R ∫₀^h r dr dθ dz = πR²h";
    }

    /**
     * @brief Jacobian determinant
     */
    static double jacobian(double r) {
        return r;
    }
};

/**
 * @class SphericalCoordinates
 * @brief 3D Spherical coordinate system (ρ, θ, φ)
 */
class SphericalCoordinates {
public:
    /**
     * @brief Convert Cartesian (x, y, z) to Spherical (ρ, θ, φ)
     *
     * ρ: radial distance
     * θ: azimuthal angle (from +x axis in xy-plane)
     * φ: polar angle (from +z axis)
     */
    static std::tuple<double, double, double> fromCartesian(double x, double y, double z) {
        double rho = std::sqrt(x*x + y*y + z*z);
        double theta = std::atan2(y, x);
        double phi = (rho > 0.0) ? std::acos(z / rho) : 0.0;
        return {rho, theta, phi};
    }

    /**
     * @brief Convert Spherical (ρ, θ, φ) to Cartesian (x, y, z)
     */
    static std::tuple<double, double, double> toCartesian(double rho, double theta, double phi) {
        if (rho < 0.0) {
            throw std::invalid_argument("Radial distance must be non-negative");
        }

        double x = rho * std::sin(phi) * std::cos(theta);
        double y = rho * std::sin(phi) * std::sin(theta);
        double z = rho * std::cos(phi);
        return {x, y, z};
    }

    /**
     * @brief Volume element in spherical coordinates
     *
     * dV = ρ² sin(φ) dρ dθ dφ
     */
    static std::string volumeElement() {
        return "Volume Element in Spherical Coordinates:\n"
               "\n"
               "dV = ρ² sin(φ) dρ dθ dφ\n"
               "\n"
               "Jacobian:\n"
               "J = |∂(x,y,z)/∂(ρ,θ,φ)| = ρ² sin(φ)\n"
               "\n"
               "Triple integral:\n"
               "∭_V f(x,y,z) dx dy dz = ∭_V f(ρ,θ,φ) ρ² sin(φ) dρ dθ dφ\n"
               "\n"
               "Example: Sphere of radius R\n"
               "Volume = ∫₀^(2π) ∫₀^π ∫₀^R ρ² sin φ dρ dφ dθ = (4/3)πR³";
    }

    /**
     * @brief Jacobian determinant
     */
    static double jacobian(double rho, double phi) {
        return rho * rho * std::sin(phi);
    }

    /**
     * @brief Convention note
     */
    static std::string convention() {
        return "Spherical Coordinate Conventions:\n"
               "\n"
               "Physics convention (ISO):\n"
               "  ρ (or r): radial distance [0, ∞)\n"
               "  θ: azimuthal angle [0, 2π)\n"
               "  φ: polar angle from +z axis [0, π]\n"
               "\n"
               "Mathematics convention:\n"
               "  Sometimes θ and φ are swapped!\n"
               "  Always check convention being used.";
    }
};

/**
 * @class ComplexPolarForm
 * @brief Polar form of complex numbers
 */
class ComplexPolarForm {
public:
    /**
     * @brief Convert complex to polar form
     *
     * z = x + iy → z = r e^(iθ) = r(cos θ + i sin θ)
     *
     * @return {r, theta} where r = |z|, theta = arg(z)
     */
    static std::pair<double, double> toPolar(const std::complex<double>& z) {
        double r = std::abs(z);
        double theta = std::arg(z);
        return {r, theta};
    }

    /**
     * @brief Convert polar to complex
     *
     * r e^(iθ) → r cos θ + i r sin θ
     */
    static std::complex<double> fromPolar(double r, double theta) {
        if (r < 0.0) {
            throw std::invalid_argument("Modulus must be non-negative");
        }
        return std::polar(r, theta);
    }

    /**
     * @brief Euler's formula
     */
    static std::string eulersFormula() {
        return "Euler's Formula:\n"
               "\n"
               "e^(iθ) = cos(θ) + i sin(θ)\n"
               "\n"
               "Special cases:\n"
               "e^(iπ) = -1  (Euler's identity)\n"
               "e^(iπ/2) = i\n"
               "e^(-iθ) = cos(θ) - i sin(θ) = [e^(iθ)]*\n"
               "\n"
               "Polar form of complex number:\n"
               "z = r e^(iθ)\n"
               "\n"
               "where r = |z| (modulus), θ = arg(z) (argument)";
    }

    /**
     * @brief Multiplication in polar form
     *
     * z₁ z₂ = r₁ e^(iθ₁) · r₂ e^(iθ₂) = (r₁r₂) e^(i(θ₁+θ₂))
     *
     * Multiply moduli, add arguments!
     */
    static std::string multiplicationRule() {
        return "Multiplication in Polar Form:\n"
               "\n"
               "z₁ = r₁ e^(iθ₁), z₂ = r₂ e^(iθ₂)\n"
               "\n"
               "z₁ · z₂ = (r₁r₂) e^(i(θ₁+θ₂))\n"
               "\n"
               "Rule: Multiply moduli, add arguments\n"
               "\n"
               "Geometric interpretation:\n"
               "- Scale by r₂\n"
               "- Rotate by θ₂";
    }

    /**
     * @brief De Moivre's theorem
     *
     * (cos θ + i sin θ)^n = cos(nθ) + i sin(nθ)
     */
    static std::string deMoivresTheorem() {
        return "De Moivre's Theorem:\n"
               "\n"
               "(cos θ + i sin θ)^n = cos(nθ) + i sin(nθ)\n"
               "\n"
               "Or in exponential form:\n"
               "(e^(iθ))^n = e^(inθ)\n"
               "\n"
               "Applications:\n"
               "1. Compute powers of complex numbers\n"
               "2. Find nth roots\n"
               "3. Derive trig identities\n"
               "\n"
               "Example: (1 + i)^10 in polar form";
    }

    /**
     * @brief nth roots of complex number
     *
     * If z = r e^(iθ), then nth roots are:
     * z^(1/n) = r^(1/n) e^(i(θ + 2πk)/n) for k = 0, 1, ..., n-1
     */
    static std::vector<std::complex<double>> nthRoots(const std::complex<double>& z, int n) {
        if (n <= 0) {
            throw std::invalid_argument("n must be positive");
        }

        auto [r, theta] = toPolar(z);
        double r_root = std::pow(r, 1.0 / n);

        std::vector<std::complex<double>> roots;
        for (int k = 0; k < n; ++k) {
            double angle = (theta + 2.0 * M_PI * k) / n;
            roots.push_back(fromPolar(r_root, angle));
        }

        return roots;
    }

    /**
     * @brief Roots of unity
     *
     * nth roots of 1: e^(2πik/n) for k = 0, 1, ..., n-1
     */
    static std::vector<std::complex<double>> rootsOfUnity(int n) {
        return nthRoots(std::complex<double>(1.0, 0.0), n);
    }

    static std::string rootsOfUnityDescription() {
        return "Roots of Unity:\n"
               "\n"
               "nth roots of 1:\n"
               "ωₖ = e^(2πik/n) for k = 0, 1, ..., n-1\n"
               "\n"
               "Properties:\n"
               "1. Evenly spaced around unit circle\n"
               "2. ω₀ = 1 (always one root at 1)\n"
               "3. Sum of all roots = 0 (for n > 1)\n"
               "4. Form cyclic group under multiplication\n"
               "\n"
               "Example (n=3): 1, e^(2πi/3), e^(4πi/3)\n"
               "           = 1, -1/2 + i√3/2, -1/2 - i√3/2";
    }
};

/**
 * @class PolarCurves
 * @brief Common curves in polar coordinates
 */
class PolarCurves {
public:
    /**
     * @brief Circle: r = R (constant)
     */
    static std::string circle() {
        return "Circle:\n"
               "r = R (constant)\n"
               "\n"
               "Circle of radius R centered at origin";
    }

    /**
     * @brief Cardioid: r = a(1 + cos θ)
     */
    static std::string cardioid() {
        return "Cardioid:\n"
               "r = a(1 + cos θ)\n"
               "\n"
               "Heart-shaped curve\n"
               "Special case of limaçon";
    }

    /**
     * @brief Rose curve: r = a cos(nθ)
     */
    static std::string rose() {
        return "Rose Curve:\n"
               "r = a cos(nθ) or r = a sin(nθ)\n"
               "\n"
               "n odd: n petals\n"
               "n even: 2n petals\n"
               "\n"
               "Example: r = cos(3θ) has 3 petals";
    }

    /**
     * @brief Spiral of Archimedes: r = aθ
     */
    static std::string spiral() {
        return "Archimedean Spiral:\n"
               "r = aθ\n"
               "\n"
               "Uniform spacing between turns\n"
               "\n"
               "Logarithmic spiral: r = ae^(bθ)\n"
               "(appears in nature: nautilus shell)";
    }

    /**
     * @brief Lemniscate: r² = a² cos(2θ)
     */
    static std::string lemniscate() {
        return "Lemniscate of Bernoulli:\n"
               "r² = a² cos(2θ)\n"
               "\n"
               "Figure-eight shape (∞)\n"
               "r = ±a√(cos 2θ) for |θ| ≤ π/4";
    }

    /**
     * @brief Area enclosed by polar curve
     *
     * A = (1/2) ∫_α^β r² dθ
     */
    static std::string areaFormula() {
        return "Area Enclosed by Polar Curve:\n"
               "\n"
               "A = (1/2) ∫_α^β r(θ)² dθ\n"
               "\n"
               "Example: Cardioid r = a(1 + cos θ)\n"
               "A = (1/2) ∫₀^(2π) a²(1 + cos θ)² dθ = (3/2)πa²";
    }

    /**
     * @brief Arc length in polar coordinates
     *
     * s = ∫_α^β √(r² + (dr/dθ)²) dθ
     */
    static std::string arcLength() {
        return "Arc Length in Polar Coordinates:\n"
               "\n"
               "s = ∫_α^β √(r² + (dr/dθ)²) dθ\n"
               "\n"
               "Example: Circle r = R\n"
               "dr/dθ = 0\n"
               "s = ∫₀^(2π) R dθ = 2πR";
    }
};

} // namespace maths::transforms

#endif // MATHS_TRANSFORMS_POLAR_HPP
