#ifndef PHYSICS_FLUID_MECHANICS_HPP
#define PHYSICS_FLUID_MECHANICS_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace fluid_mechanics {

/**
 * @brief Fluid Mechanics
 *
 * Functions for analyzing fluid statics and dynamics including:
 * - Pressure in fluids
 * - Barometer measurements
 * - Fluid flow and continuity
 * - Bernoulli's equation
 * - Efflux and jets
 * - Pipe flow and friction
 *
 * Key principles:
 * - Continuity equation: A₁v₁ = A₂v₂
 * - Bernoulli's equation: P + ½ρv² + ρgh = constant
 * - Torricelli's theorem: v = √(2gh)
 */

// ============================================================================
// Physical Constants
// ============================================================================

namespace constants {
    constexpr double WATER_DENSITY = 1000.0;     // kg/m³
    constexpr double MERCURY_DENSITY = 13600.0;  // kg/m³
    constexpr double AIR_DENSITY_STP = 1.225;    // kg/m³ at STP
    constexpr double ATM_PRESSURE = 101325.0;    // Pa (1 atm)
    constexpr double GRAVITY = 9.81;             // m/s²
}

// ============================================================================
// Pressure in Fluids
// ============================================================================

/**
 * @brief Calculate pressure at depth in a fluid
 *
 * P = P₀ + ρgh
 * where P₀ is atmospheric pressure, ρ is fluid density, h is depth
 *
 * @param depth Depth below surface (in meters, must be >= 0)
 * @param fluidDensity Fluid density (in kg/m³, must be > 0)
 * @param atmosphericPressure Surface pressure (in Pa, default: 101325 Pa)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Pressure at depth (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double pressureAtDepth(double depth, double fluidDensity,
                              double atmosphericPressure = constants::ATM_PRESSURE,
                              double gravity = constants::GRAVITY) {
    if (depth < 0) {
        throw std::invalid_argument("Depth must be non-negative");
    }
    if (fluidDensity <= 0) {
        throw std::invalid_argument("Fluid density must be positive");
    }
    return atmosphericPressure + fluidDensity * gravity * depth;
}

/**
 * @brief Calculate gauge pressure (pressure above atmospheric)
 *
 * P_gauge = ρgh
 *
 * @param depth Depth below surface (in meters, must be >= 0)
 * @param fluidDensity Fluid density (in kg/m³, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Gauge pressure (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double gaugePressure(double depth, double fluidDensity, double gravity = constants::GRAVITY) {
    if (depth < 0) {
        throw std::invalid_argument("Depth must be non-negative");
    }
    if (fluidDensity <= 0) {
        throw std::invalid_argument("Fluid density must be positive");
    }
    return fluidDensity * gravity * depth;
}

// ============================================================================
// Measurement of Heights by Barometer
// ============================================================================

/**
 * @brief Calculate atmospheric pressure from barometer height
 *
 * P_atm = ρ_mercury × g × h
 *
 * Standard atmospheric pressure supports 760 mm of mercury
 *
 * @param mercuryHeight Height of mercury column (in meters)
 * @param mercuryDensity Density of mercury (default: 13600 kg/m³)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Atmospheric pressure (in Pa)
 * @throws std::invalid_argument if mercuryHeight < 0 or mercuryDensity <= 0
 */
inline double barometerPressure(double mercuryHeight,
                                double mercuryDensity = constants::MERCURY_DENSITY,
                                double gravity = constants::GRAVITY) {
    if (mercuryHeight < 0) {
        throw std::invalid_argument("Mercury height must be non-negative");
    }
    if (mercuryDensity <= 0) {
        throw std::invalid_argument("Mercury density must be positive");
    }
    return mercuryDensity * gravity * mercuryHeight;
}

/**
 * @brief Calculate altitude from pressure difference
 *
 * Using barometric formula (for small altitude changes):
 * Δh ≈ (P₁ - P₂)/(ρ_air × g)
 *
 * @param pressure1 Pressure at lower altitude (in Pa, must be > 0)
 * @param pressure2 Pressure at higher altitude (in Pa, must be > 0)
 * @param airDensity Density of air (default: 1.225 kg/m³)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Altitude difference (in meters)
 * @throws std::invalid_argument if parameters out of range
 */
inline double altitudeFromPressure(double pressure1, double pressure2,
                                   double airDensity = constants::AIR_DENSITY_STP,
                                   double gravity = constants::GRAVITY) {
    if (pressure1 <= 0 || pressure2 <= 0) {
        throw std::invalid_argument("Pressures must be positive");
    }
    if (airDensity <= 0) {
        throw std::invalid_argument("Air density must be positive");
    }
    return (pressure1 - pressure2) / (airDensity * gravity);
}

// ============================================================================
// Fluid Steady Flow and Continuity
// ============================================================================

/**
 * @brief Continuity equation: Calculate velocity at second point
 *
 * For incompressible fluid: A₁v₁ = A₂v₂
 * Mass flow rate is constant: ρA₁v₁ = ρA₂v₂
 *
 * @param area1 Cross-sectional area at point 1 (in m², must be > 0)
 * @param velocity1 Velocity at point 1 (in m/s, must be >= 0)
 * @param area2 Cross-sectional area at point 2 (in m², must be > 0)
 * @return Velocity at point 2 (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double continuityEquation(double area1, double velocity1, double area2) {
    if (area1 <= 0 || area2 <= 0) {
        throw std::invalid_argument("Areas must be positive");
    }
    if (velocity1 < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    return (area1 * velocity1) / area2;
}

/**
 * @brief Calculate volume flow rate (discharge)
 *
 * Q = A × v (volume per unit time)
 *
 * @param area Cross-sectional area (in m², must be > 0)
 * @param velocity Flow velocity (in m/s, must be >= 0)
 * @return Volume flow rate (in m³/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double volumeFlowRate(double area, double velocity) {
    if (area <= 0) {
        throw std::invalid_argument("Area must be positive");
    }
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    return area * velocity;
}

/**
 * @brief Calculate mass flow rate
 *
 * ṁ = ρAv (mass per unit time)
 *
 * @param density Fluid density (in kg/m³, must be > 0)
 * @param area Cross-sectional area (in m², must be > 0)
 * @param velocity Flow velocity (in m/s, must be >= 0)
 * @return Mass flow rate (in kg/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double massFlowRate(double density, double area, double velocity) {
    if (density <= 0) {
        throw std::invalid_argument("Density must be positive");
    }
    if (area <= 0) {
        throw std::invalid_argument("Area must be positive");
    }
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    return density * area * velocity;
}

// ============================================================================
// Bernoulli's Equation
// ============================================================================

/**
 * @brief Energy due to pressure (pressure energy per unit volume)
 *
 * Energy density = P (pressure itself represents energy per unit volume)
 *
 * @param pressure Pressure (in Pa, must be >= 0)
 * @return Pressure energy per unit volume (in J/m³ = Pa)
 */
inline double pressureEnergy(double pressure) {
    if (pressure < 0) {
        throw std::invalid_argument("Pressure must be non-negative");
    }
    return pressure;
}

/**
 * @brief Bernoulli's equation: Calculate pressure at second point
 *
 * P₁ + ½ρv₁² + ρgh₁ = P₂ + ½ρv₂² + ρgh₂
 * Therefore: P₂ = P₁ + ½ρ(v₁² - v₂²) + ρg(h₁ - h₂)
 *
 * @param pressure1 Pressure at point 1 (in Pa, must be >= 0)
 * @param velocity1 Velocity at point 1 (in m/s, must be >= 0)
 * @param height1 Height at point 1 (in meters)
 * @param velocity2 Velocity at point 2 (in m/s, must be >= 0)
 * @param height2 Height at point 2 (in meters)
 * @param density Fluid density (in kg/m³, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Pressure at point 2 (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double bernoulliPressure(double pressure1, double velocity1, double height1,
                                double velocity2, double height2, double density,
                                double gravity = constants::GRAVITY) {
    if (pressure1 < 0) {
        throw std::invalid_argument("Pressure must be non-negative");
    }
    if (velocity1 < 0 || velocity2 < 0) {
        throw std::invalid_argument("Velocities must be non-negative");
    }
    if (density <= 0) {
        throw std::invalid_argument("Density must be positive");
    }

    double kineticTerm = 0.5 * density * (velocity1 * velocity1 - velocity2 * velocity2);
    double potentialTerm = density * gravity * (height1 - height2);

    return pressure1 + kineticTerm + potentialTerm;
}

/**
 * @brief Calculate total mechanical energy per unit volume (Bernoulli constant)
 *
 * E/V = P + ½ρv² + ρgh
 *
 * @param pressure Pressure (in Pa, must be >= 0)
 * @param velocity Flow velocity (in m/s, must be >= 0)
 * @param height Height (in meters)
 * @param density Fluid density (in kg/m³, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Total energy per unit volume (in J/m³)
 * @throws std::invalid_argument if parameters out of range
 */
inline double bernoulliConstant(double pressure, double velocity, double height, double density,
                                double gravity = constants::GRAVITY) {
    if (pressure < 0) {
        throw std::invalid_argument("Pressure must be non-negative");
    }
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    if (density <= 0) {
        throw std::invalid_argument("Density must be positive");
    }

    return pressure + 0.5 * density * velocity * velocity + density * gravity * height;
}

// ============================================================================
// Velocity of a Jet / Efflux (Torricelli's Theorem)
// ============================================================================

/**
 * @brief Torricelli's theorem: Velocity of efflux from tank
 *
 * v = √(2gh)
 * where h is the depth below the free surface
 *
 * This is the speed of fluid flowing out of an orifice
 *
 * @param depth Depth below free surface (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Efflux velocity (in m/s)
 * @throws std::invalid_argument if depth <= 0
 */
inline double effluxVelocity(double depth, double gravity = constants::GRAVITY) {
    if (depth <= 0) {
        throw std::invalid_argument("Depth must be positive");
    }
    return std::sqrt(2.0 * gravity * depth);
}

/**
 * @brief Calculate jet velocity considering pressure difference
 *
 * For pressurized tank: v = √(2(P₁ - P₂)/ρ + 2gh)
 * where P₁ is tank pressure, P₂ is outside pressure
 *
 * @param pressureDiff Pressure difference P₁ - P₂ (in Pa)
 * @param depth Depth below surface (in meters, must be >= 0)
 * @param density Fluid density (in kg/m³, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Jet velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double jetVelocityPressurized(double pressureDiff, double depth, double density,
                                     double gravity = constants::GRAVITY) {
    if (depth < 0) {
        throw std::invalid_argument("Depth must be non-negative");
    }
    if (density <= 0) {
        throw std::invalid_argument("Density must be positive");
    }

    return std::sqrt(2.0 * pressureDiff / density + 2.0 * gravity * depth);
}

/**
 * @brief Calculate horizontal range of jet
 *
 * For jet exiting at depth h below surface, height H above ground:
 * Range = 2√(h(H - h))
 *
 * Maximum range occurs when h = H/2
 *
 * @param depthBelowSurface Depth of orifice below free surface (in meters, must be > 0)
 * @param heightAboveGround Height of orifice above ground (in meters, must be > 0)
 * @return Horizontal range (in meters)
 * @throws std::invalid_argument if parameters out of range
 */
inline double jetRange(double depthBelowSurface, double heightAboveGround) {
    if (depthBelowSurface <= 0) {
        throw std::invalid_argument("Depth below surface must be positive");
    }
    if (heightAboveGround <= 0) {
        throw std::invalid_argument("Height above ground must be positive");
    }

    return 2.0 * std::sqrt(depthBelowSurface * heightAboveGround);
}

// ============================================================================
// Efflux of Gases
// ============================================================================

/**
 * @brief Calculate efflux velocity of gas (isothermal)
 *
 * For isothermal expansion: v = √(2P/ρ)
 * where P is pressure difference, ρ is gas density
 *
 * @param pressureDiff Pressure difference (in Pa, must be > 0)
 * @param gasDensity Gas density (in kg/m³, must be > 0)
 * @return Efflux velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double gasEffluxIsothermal(double pressureDiff, double gasDensity) {
    if (pressureDiff <= 0) {
        throw std::invalid_argument("Pressure difference must be positive");
    }
    if (gasDensity <= 0) {
        throw std::invalid_argument("Gas density must be positive");
    }
    return std::sqrt(2.0 * pressureDiff / gasDensity);
}

/**
 * @brief Calculate efflux velocity of gas (adiabatic)
 *
 * For adiabatic expansion: v = √(2γP/(ρ(γ-1)))
 * where γ is heat capacity ratio
 *
 * @param pressureDiff Pressure difference (in Pa, must be > 0)
 * @param gasDensity Gas density (in kg/m³, must be > 0)
 * @param gamma Heat capacity ratio γ = Cp/Cv (must be > 1)
 * @return Efflux velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double gasEffluxAdiabatic(double pressureDiff, double gasDensity, double gamma) {
    if (pressureDiff <= 0) {
        throw std::invalid_argument("Pressure difference must be positive");
    }
    if (gasDensity <= 0) {
        throw std::invalid_argument("Gas density must be positive");
    }
    if (gamma <= 1.0) {
        throw std::invalid_argument("Gamma must be greater than 1");
    }
    return std::sqrt(2.0 * gamma * pressureDiff / (gasDensity * (gamma - 1.0)));
}

// ============================================================================
// Friction in Pipes
// ============================================================================

/**
 * @brief Calculate pressure drop due to friction in pipe (Darcy-Weisbach)
 *
 * ΔP = f × (L/D) × (ρv²/2)
 * where f is friction factor, L is length, D is diameter
 *
 * @param frictionFactor Darcy friction factor (dimensionless, must be > 0)
 * @param length Pipe length (in meters, must be > 0)
 * @param diameter Pipe diameter (in meters, must be > 0)
 * @param density Fluid density (in kg/m³, must be > 0)
 * @param velocity Flow velocity (in m/s, must be >= 0)
 * @return Pressure drop (in Pa)
 * @throws std::invalid_argument if parameters out of range
 */
inline double pipeFrictionPressureDrop(double frictionFactor, double length, double diameter,
                                       double density, double velocity) {
    if (frictionFactor <= 0) {
        throw std::invalid_argument("Friction factor must be positive");
    }
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (diameter <= 0) {
        throw std::invalid_argument("Diameter must be positive");
    }
    if (density <= 0) {
        throw std::invalid_argument("Density must be positive");
    }
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }

    return frictionFactor * (length / diameter) * (0.5 * density * velocity * velocity);
}

/**
 * @brief Calculate head loss due to friction
 *
 * h_f = f × (L/D) × (v²/(2g))
 *
 * @param frictionFactor Darcy friction factor (dimensionless, must be > 0)
 * @param length Pipe length (in meters, must be > 0)
 * @param diameter Pipe diameter (in meters, must be > 0)
 * @param velocity Flow velocity (in m/s, must be >= 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Head loss (in meters of fluid column)
 * @throws std::invalid_argument if parameters out of range
 */
inline double headLossFriction(double frictionFactor, double length, double diameter,
                               double velocity, double gravity = constants::GRAVITY) {
    if (frictionFactor <= 0) {
        throw std::invalid_argument("Friction factor must be positive");
    }
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (diameter <= 0) {
        throw std::invalid_argument("Diameter must be positive");
    }
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }

    return frictionFactor * (length / diameter) * (velocity * velocity) / (2.0 * gravity);
}

/**
 * @brief Calculate Reynolds number for flow characterization
 *
 * Re = ρvD/μ = vD/ν
 * where μ is dynamic viscosity, ν is kinematic viscosity
 *
 * Re < 2300: Laminar flow
 * 2300 < Re < 4000: Transition
 * Re > 4000: Turbulent flow
 *
 * @param velocity Flow velocity (in m/s, must be >= 0)
 * @param diameter Characteristic length (pipe diameter) (in meters, must be > 0)
 * @param kinematicViscosity Kinematic viscosity ν (in m²/s, must be > 0)
 * @return Reynolds number (dimensionless)
 * @throws std::invalid_argument if parameters out of range
 */
inline double reynoldsNumber(double velocity, double diameter, double kinematicViscosity) {
    if (velocity < 0) {
        throw std::invalid_argument("Velocity must be non-negative");
    }
    if (diameter <= 0) {
        throw std::invalid_argument("Diameter must be positive");
    }
    if (kinematicViscosity <= 0) {
        throw std::invalid_argument("Kinematic viscosity must be positive");
    }

    return (velocity * diameter) / kinematicViscosity;
}

/**
 * @brief Power loss due to friction in pipe
 *
 * Power = ΔP × Q = ΔP × A × v
 * where Q is volume flow rate
 *
 * @param pressureDrop Pressure drop due to friction (in Pa, must be >= 0)
 * @param volumeFlowRate Volume flow rate (in m³/s, must be >= 0)
 * @return Power loss (in Watts)
 * @throws std::invalid_argument if parameters out of range
 */
inline double powerLossFriction(double pressureDrop, double volumeFlowRate) {
    if (pressureDrop < 0) {
        throw std::invalid_argument("Pressure drop must be non-negative");
    }
    if (volumeFlowRate < 0) {
        throw std::invalid_argument("Volume flow rate must be non-negative");
    }

    return pressureDrop * volumeFlowRate;
}

} // namespace fluid_mechanics
} // namespace physics

#endif // PHYSICS_FLUID_MECHANICS_HPP
