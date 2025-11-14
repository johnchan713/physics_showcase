#ifndef PHYSICS_UNITS_HPP
#define PHYSICS_UNITS_HPP

#include <cmath>

namespace physics {
namespace units {

/**
 * @brief Unit Conversions and Relationships
 *
 * This namespace provides functions for converting between different unit systems
 * and understanding relationships between fundamental quantities.
 *
 * Common unit systems:
 * - SI (International System): meter, kilogram, second, Newton
 * - CGS (Centimeter-Gram-Second): centimeter, gram, second, Dyne
 * - Imperial: foot, pound, second
 */

// ============================================================================
// CGS and SI Unit Conversions
// ============================================================================

/**
 * @brief Convert force from Newtons to Dynes
 *
 * Relationship: 1 Newton = 10^5 Dynes
 *
 * The Dyne is the CGS unit of force, defined as the force required to
 * accelerate a mass of one gram at a rate of one centimeter per second squared.
 *
 * @param newtons Force in Newtons (SI unit)
 * @return Force in Dynes (CGS unit)
 */
inline double newtonsToДynes(double newtons) {
    return newtons * 1.0e5;
}

/**
 * @brief Convert force from Dynes to Newtons
 *
 * Relationship: 1 Dyne = 10^-5 Newtons
 *
 * @param dynes Force in Dynes (CGS unit)
 * @return Force in Newtons (SI unit)
 */
inline double dynesToNewtons(double dynes) {
    return dynes * 1.0e-5;
}

/**
 * @brief Relationship between Dyne and Gram
 *
 * Calculate force in Dynes from mass in grams and acceleration in cm/s²
 * F(dyne) = m(g) × a(cm/s²)
 *
 * This demonstrates the definition: 1 Dyne = 1 gram × 1 cm/s²
 *
 * @param massGrams Mass in grams
 * @param accelerationCmPerSec2 Acceleration in cm/s²
 * @return Force in Dynes
 */
inline double calculateForceDynes(double massGrams, double accelerationCmPerSec2) {
    return massGrams * accelerationCmPerSec2;
}

/**
 * @brief Convert mass from kilograms to grams
 *
 * Relationship: 1 kg = 1000 g
 *
 * @param kilograms Mass in kilograms
 * @return Mass in grams
 */
inline double kilogramsToGrams(double kilograms) {
    return kilograms * 1000.0;
}

/**
 * @brief Convert mass from grams to kilograms
 *
 * Relationship: 1 g = 0.001 kg
 *
 * @param grams Mass in grams
 * @return Mass in kilograms
 */
inline double gramsToKilograms(double grams) {
    return grams / 1000.0;
}

/**
 * @brief Demonstrate equivalence: 1N = 1kg⋅m/s² = 10^5 g⋅cm/s² = 10^5 Dynes
 *
 * Given a force in SI units (mass in kg, acceleration in m/s²),
 * calculate the equivalent force in CGS units (Dynes)
 *
 * @param massKg Mass in kilograms
 * @param accelMPerSec2 Acceleration in m/s²
 * @return Force in Dynes
 */
inline double forceInDynesFromSI(double massKg, double accelMPerSec2) {
    double forceNewtons = massKg * accelMPerSec2;
    return newtonsToДynes(forceNewtons);
}

// ============================================================================
// Length Conversions
// ============================================================================

/**
 * @brief Convert meters to centimeters
 *
 * Relationship: 1 m = 100 cm
 *
 * @param meters Length in meters
 * @return Length in centimeters
 */
inline double metersToCentimeters(double meters) {
    return meters * 100.0;
}

/**
 * @brief Convert centimeters to meters
 *
 * Relationship: 1 cm = 0.01 m
 *
 * @param centimeters Length in centimeters
 * @return Length in meters
 */
inline double centimetersToMeters(double centimeters) {
    return centimeters / 100.0;
}

/**
 * @brief Convert meters to feet
 *
 * Relationship: 1 m ≈ 3.28084 ft
 *
 * @param meters Length in meters
 * @return Length in feet
 */
inline double metersToFeet(double meters) {
    return meters * 3.28084;
}

/**
 * @brief Convert feet to meters
 *
 * Relationship: 1 ft ≈ 0.3048 m
 *
 * @param feet Length in feet
 * @return Length in meters
 */
inline double feetToMeters(double feet) {
    return feet * 0.3048;
}

// ============================================================================
// Velocity Conversions
// ============================================================================

/**
 * @brief Convert velocity from m/s to cm/s
 *
 * Relationship: 1 m/s = 100 cm/s
 *
 * @param metersPerSecond Velocity in m/s
 * @return Velocity in cm/s
 */
inline double mpsToЛmps(double metersPerSecond) {
    return metersPerSecond * 100.0;
}

/**
 * @brief Convert velocity from cm/s to m/s
 *
 * Relationship: 1 cm/s = 0.01 m/s
 *
 * @param centimetersPerSecond Velocity in cm/s
 * @return Velocity in m/s
 */
inline double cmpsToMps(double centimetersPerSecond) {
    return centimetersPerSecond / 100.0;
}

/**
 * @brief Convert velocity from m/s to km/h
 *
 * Relationship: 1 m/s = 3.6 km/h
 *
 * @param metersPerSecond Velocity in m/s
 * @return Velocity in km/h
 */
inline double mpsToKmph(double metersPerSecond) {
    return metersPerSecond * 3.6;
}

/**
 * @brief Convert velocity from km/h to m/s
 *
 * Relationship: 1 km/h ≈ 0.27778 m/s
 *
 * @param kilometersPerHour Velocity in km/h
 * @return Velocity in m/s
 */
inline double kmphToMps(double kilometersPerHour) {
    return kilometersPerHour / 3.6;
}

/**
 * @brief Convert velocity from m/s to mph (miles per hour)
 *
 * Relationship: 1 m/s ≈ 2.23694 mph
 *
 * @param metersPerSecond Velocity in m/s
 * @return Velocity in mph
 */
inline double mpsToMph(double metersPerSecond) {
    return metersPerSecond * 2.23694;
}

/**
 * @brief Convert velocity from mph to m/s
 *
 * Relationship: 1 mph ≈ 0.44704 m/s
 *
 * @param milesPerHour Velocity in mph
 * @return Velocity in m/s
 */
inline double mphToMps(double milesPerHour) {
    return milesPerHour * 0.44704;
}

// ============================================================================
// Energy Conversions
// ============================================================================

/**
 * @brief Convert energy from Joules to ergs
 *
 * Relationship: 1 Joule = 10^7 ergs
 * The erg is the CGS unit of energy
 *
 * @param joules Energy in Joules (SI)
 * @return Energy in ergs (CGS)
 */
inline double joulesToErgs(double joules) {
    return joules * 1.0e7;
}

/**
 * @brief Convert energy from ergs to Joules
 *
 * Relationship: 1 erg = 10^-7 Joules
 *
 * @param ergs Energy in ergs (CGS)
 * @return Energy in Joules (SI)
 */
inline double ergsToJoules(double ergs) {
    return ergs * 1.0e-7;
}

/**
 * @brief Convert energy from Joules to calories
 *
 * Relationship: 1 Joule ≈ 0.239006 calories
 *
 * @param joules Energy in Joules
 * @return Energy in calories
 */
inline double joulesToCalories(double joules) {
    return joules * 0.239006;
}

/**
 * @brief Convert energy from calories to Joules
 *
 * Relationship: 1 calorie ≈ 4.184 Joules
 *
 * @param calories Energy in calories
 * @return Energy in Joules
 */
inline double caloriesToJoules(double calories) {
    return calories * 4.184;
}

// ============================================================================
// Acceleration Conversions
// ============================================================================

/**
 * @brief Convert acceleration from m/s² to cm/s²
 *
 * Relationship: 1 m/s² = 100 cm/s²
 *
 * @param metersPerSecondSquared Acceleration in m/s²
 * @return Acceleration in cm/s²
 */
inline double accelMpsToĹmps(double metersPerSecondSquared) {
    return metersPerSecondSquared * 100.0;
}

/**
 * @brief Convert acceleration from cm/s² to m/s²
 *
 * Relationship: 1 cm/s² = 0.01 m/s²
 *
 * @param centimetersPerSecondSquared Acceleration in cm/s²
 * @return Acceleration in m/s²
 */
inline double accelCmpsToMps(double centimetersPerSecondSquared) {
    return centimetersPerSecondSquared / 100.0;
}

/**
 * @brief Express acceleration as multiples of g (Earth's gravity)
 *
 * Useful for expressing high accelerations (e.g., in aircraft, rockets)
 *
 * @param accelerationMps2 Acceleration in m/s²
 * @param g Acceleration due to gravity (default: 9.81 m/s²)
 * @return Acceleration in multiples of g
 */
inline double accelerationInGs(double accelerationMps2, double g = 9.81) {
    return accelerationMps2 / g;
}

// ============================================================================
// Angle Conversions
// ============================================================================

/**
 * @brief Convert degrees to radians
 *
 * Relationship: 1 degree = π/180 radians
 *
 * @param degrees Angle in degrees
 * @return Angle in radians
 */
inline double degreesToRadians(double degrees) {
    return degrees * M_PI / 180.0;
}

/**
 * @brief Convert radians to degrees
 *
 * Relationship: 1 radian = 180/π degrees
 *
 * @param radians Angle in radians
 * @return Angle in degrees
 */
inline double radiansToDegrees(double radians) {
    return radians * 180.0 / M_PI;
}

} // namespace units
} // namespace physics

#endif // PHYSICS_UNITS_HPP
