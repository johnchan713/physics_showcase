#include <iostream>
#include <iomanip>
#include "../include/physics/units.hpp"
#include "../include/physics/inclined_plane.hpp"
#include "../include/physics/energy_momentum.hpp"
#include "../include/physics/projectile.hpp"
#include "../include/physics/circular_motion.hpp"
#include "../include/physics/orbital.hpp"

// Helper function to print section headers
void printSection(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(70, '=') << "\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(2);

    std::cout << "PHYSICS SHOWCASE - ADVANCED TOPICS\n";
    std::cout << "Unit Conversions, Energy, Momentum, Projectiles, Orbits, and More\n";

    // ========================================================================
    // UNIT CONVERSIONS
    // ========================================================================

    printSection("UNIT CONVERSIONS");

    std::cout << "\n--- Relation between Dyne and Gram ---\n";
    std::cout << "In CGS system: 1 Dyne = 1 gram × 1 cm/s²\n";

    double massGrams = 100.0;
    double accelCmPerS2 = 50.0;
    double forceDynes = physics::units::calculateForceDynes(massGrams, accelCmPerS2);

    std::cout << "Mass: " << massGrams << " g\n";
    std::cout << "Acceleration: " << accelCmPerS2 << " cm/s²\n";
    std::cout << "Force: " << forceDynes << " Dynes\n";

    std::cout << "\nConversion to SI units:\n";
    double massKg = physics::units::gramsToKilograms(massGrams);
    double accelMPerS2 = physics::units::accelCmpsToMps(accelCmPerS2);
    double forceNewtons = physics::units::dynesToNewtons(forceDynes);

    std::cout << "Mass: " << massKg << " kg\n";
    std::cout << "Acceleration: " << accelMPerS2 << " m/s²\n";
    std::cout << "Force: " << forceNewtons << " N\n";

    std::cout << "\nVerification: 1 N = 10^5 Dynes\n";
    std::cout << "100 N = " << physics::units::newtonsToДynes(100.0) << " Dynes\n";
    std::cout << "1,000,000 Dynes = " << physics::units::dynesToNewtons(1000000.0) << " N\n";

    std::cout << "\n--- Other Unit Conversions ---\n";
    double speed = 100.0; // km/h
    std::cout << speed << " km/h = " << physics::units::kmphToMps(speed) << " m/s\n";

    double angle = 45.0; // degrees
    std::cout << angle << " degrees = " << physics::units::degreesToRadians(angle) << " radians\n";

    // ========================================================================
    // INCLINED PLANE - VELOCITY AT FOOT
    // ========================================================================

    printSection("INCLINED PLANE - VELOCITY AT FOOT");

    std::cout << "\n--- Object Sliding Down Frictionless Plane ---\n";
    double planeLength = 10.0;      // meters
    double angleRad = physics::units::degreesToRadians(30.0);

    std::cout << "Incline length: " << planeLength << " m\n";
    std::cout << "Angle: 30 degrees\n";

    double velocityAtFoot = physics::inclined_plane::calculateVelocityAtFootFrictionless(
        planeLength, angleRad);
    double timeToSlide = physics::inclined_plane::calculateTimeToSlide(planeLength, angleRad);

    std::cout << "Velocity at foot: " << velocityAtFoot << " m/s\n";
    std::cout << "Time to slide down: " << timeToSlide << " s\n";

    double height = planeLength * std::sin(angleRad);
    std::cout << "\nVerification using energy: height = " << height << " m\n";
    double velocityFromEnergy = physics::inclined_plane::calculateVelocityAtFootFromHeight(height);
    std::cout << "Velocity from energy conservation: " << velocityFromEnergy << " m/s\n";

    std::cout << "\n--- With Friction ---\n";
    double frictionCoeff = 0.1;
    std::cout << "Coefficient of friction: " << frictionCoeff << "\n";

    double velocityWithFriction = physics::inclined_plane::calculateVelocityAtFootWithFriction(
        planeLength, angleRad, frictionCoeff);
    double timeWithFriction = physics::inclined_plane::calculateTimeToSlideWithFriction(
        planeLength, angleRad, frictionCoeff);

    std::cout << "Velocity at foot: " << velocityWithFriction << " m/s\n";
    std::cout << "Time to slide down: " << timeWithFriction << " s\n";

    // ========================================================================
    // KINETIC ENERGY AND MOMENTUM COMPARED
    // ========================================================================

    printSection("KINETIC ENERGY AND MOMENTUM COMPARED");

    std::cout << "\n--- Two Objects with Same Momentum ---\n";
    double momentum = 100.0; // kg⋅m/s
    double mass1 = 10.0;     // kg
    double mass2 = 20.0;     // kg

    std::cout << "Both objects have momentum: " << momentum << " kg⋅m/s\n";
    std::cout << "Object 1 mass: " << mass1 << " kg\n";
    std::cout << "Object 2 mass: " << mass2 << " kg\n";

    double v1 = physics::energy_momentum::calculateVelocityFromMomentum(momentum, mass1);
    double v2 = physics::energy_momentum::calculateVelocityFromMomentum(momentum, mass2);
    double ke1 = physics::energy_momentum::calculateKineticEnergy(mass1, v1);
    double ke2 = physics::energy_momentum::calculateKineticEnergy(mass2, v2);

    std::cout << "\nObject 1: velocity = " << v1 << " m/s, KE = " << ke1 << " J\n";
    std::cout << "Object 2: velocity = " << v2 << " m/s, KE = " << ke2 << " J\n";
    std::cout << "Lighter object has MORE kinetic energy for same momentum!\n";
    std::cout << "KE ratio (1/2): " << (ke1 / ke2) << " = mass ratio (2/1): " << (mass2 / mass1) << "\n";

    std::cout << "\n--- Two Objects with Same Kinetic Energy ---\n";
    double kineticEnergy = 500.0; // Joules

    std::cout << "Both objects have KE: " << kineticEnergy << " J\n";
    std::cout << "Object 1 mass: " << mass1 << " kg\n";
    std::cout << "Object 2 mass: " << mass2 << " kg\n";

    v1 = physics::energy_momentum::calculateVelocityFromKE(kineticEnergy, mass1);
    v2 = physics::energy_momentum::calculateVelocityFromKE(kineticEnergy, mass2);
    double p1 = physics::energy_momentum::calculateMomentum(mass1, v1);
    double p2 = physics::energy_momentum::calculateMomentum(mass2, v2);

    std::cout << "\nObject 1: velocity = " << v1 << " m/s, momentum = " << p1 << " kg⋅m/s\n";
    std::cout << "Object 2: velocity = " << v2 << " m/s, momentum = " << p2 << " kg⋅m/s\n";
    std::cout << "Heavier object has MORE momentum for same kinetic energy!\n";
    std::cout << "Momentum ratio (2/1): " << (p2 / p1) << " = √(mass ratio) = " << std::sqrt(mass2 / mass1) << "\n";

    std::cout << "\n--- Relationship: KE = p²/(2m) ---\n";
    double testMomentum = 50.0;
    double testMass = 5.0;
    double keFromMomentum = physics::energy_momentum::calculateKEFromMomentum(testMomentum, testMass);
    double velocityFromP = testMomentum / testMass;
    double keDirect = 0.5 * testMass * velocityFromP * velocityFromP;

    std::cout << "Momentum: " << testMomentum << " kg⋅m/s, Mass: " << testMass << " kg\n";
    std::cout << "KE from p²/(2m): " << keFromMomentum << " J\n";
    std::cout << "KE from (1/2)mv²: " << keDirect << " J\n";

    // ========================================================================
    // PROJECTILE MOTION
    // ========================================================================

    printSection("MOTION OF A PROJECTILE");

    std::cout << "\n--- Projectile Launched at 45° ---\n";
    double initialVel = 50.0;       // m/s
    double launchAngle = physics::units::degreesToRadians(45.0);

    std::cout << "Initial velocity: " << initialVel << " m/s\n";
    std::cout << "Launch angle: 45 degrees\n";

    double range = physics::projectile::calculateRange(initialVel, launchAngle);
    double maxHeight = physics::projectile::calculateMaxHeight(initialVel, launchAngle);
    double flightTime = physics::projectile::calculateTimeOfFlight(initialVel, launchAngle);

    std::cout << "Range: " << range << " m\n";
    std::cout << "Maximum height: " << maxHeight << " m\n";
    std::cout << "Flight time: " << flightTime << " s\n";

    std::cout << "\n--- Position and Velocity at t = 2s ---\n";
    double t = 2.0;
    double xPos = physics::projectile::calculateHorizontalPosition(initialVel, launchAngle, t);
    double yPos = physics::projectile::calculateVerticalPosition(initialVel, launchAngle, t);
    double vx = physics::projectile::getHorizontalVelocityAtTime(initialVel, launchAngle);
    double vy = physics::projectile::getVerticalVelocityAtTime(initialVel, launchAngle, t);
    double speed = physics::projectile::calculateSpeedAtTime(initialVel, launchAngle, t);

    std::cout << "Position: (" << xPos << ", " << yPos << ") m\n";
    std::cout << "Velocity: (" << vx << ", " << vy << ") m/s\n";
    std::cout << "Speed: " << speed << " m/s\n";

    std::cout << "\n--- Maximum Range ---\n";
    double maxRange = physics::projectile::calculateMaximumRange(initialVel);
    std::cout << "Maximum possible range: " << maxRange << " m (at 45°)\n";
    std::cout << "Optimal angle: " << physics::units::radiansToDegrees(
        physics::projectile::getAngleForMaxRange()) << " degrees\n";

    // ========================================================================
    // CIRCULAR MOTION
    // ========================================================================

    printSection("MOTION IN ANY CIRCLE WITH CONSTANT SPEED");

    std::cout << "\n--- Object Moving in Horizontal Circle ---\n";
    double radius = 5.0;            // meters
    double tangentialVel = 10.0;    // m/s

    std::cout << "Radius: " << radius << " m\n";
    std::cout << "Tangential velocity: " << tangentialVel << " m/s\n";

    double centripetalAccel = physics::circular_motion::calculateCentripetalAcceleration(
        tangentialVel, radius);
    double angularVel = physics::circular_motion::calculateAngularVelocity(tangentialVel, radius);
    double period = physics::circular_motion::calculatePeriod(tangentialVel, radius);
    double frequency = physics::circular_motion::calculateFrequency(tangentialVel, radius);

    std::cout << "Centripetal acceleration: " << centripetalAccel << " m/s²\n";
    std::cout << "Angular velocity: " << angularVel << " rad/s\n";
    std::cout << "Period: " << period << " s\n";
    std::cout << "Frequency: " << frequency << " Hz\n";

    std::cout << "\n--- Centripetal Force Required ---\n";
    double objectMass = 2.0; // kg
    std::cout << "Mass: " << objectMass << " kg\n";

    double centripetalForce = physics::circular_motion::calculateCentripetalForce(
        objectMass, tangentialVel, radius);
    std::cout << "Centripetal force required: " << centripetalForce << " N\n";

    std::cout << "\n--- Banked Curve (Car Turning) ---\n";
    double carSpeed = 20.0;         // m/s
    double curveRadius = 50.0;      // m

    std::cout << "Car speed: " << carSpeed << " m/s (" <<
        physics::units::mpsToKmph(carSpeed) << " km/h)\n";
    std::cout << "Curve radius: " << curveRadius << " m\n";

    double bankAngle = physics::circular_motion::calculateBankingAngle(carSpeed, curveRadius);
    std::cout << "Required banking angle: " <<
        physics::units::radiansToDegrees(bankAngle) << " degrees\n";

    std::cout << "\n--- Vertical Circle (Loop-the-Loop) ---\n";
    double loopRadius = 10.0; // m

    std::cout << "Loop radius: " << loopRadius << " m\n";
    double minSpeedTop = physics::circular_motion::calculateMinVelocityTopOfLoop(loopRadius);
    std::cout << "Minimum speed at top to maintain contact: " << minSpeedTop << " m/s\n";

    double speedBottom = 20.0; // m/s
    double objectMass2 = 50.0; // kg
    std::cout << "\nObject with mass " << objectMass2 << " kg at speed " << speedBottom << " m/s\n";

    double tensionBottom = physics::circular_motion::calculateTensionAtBottom(
        objectMass2, speedBottom, loopRadius);
    std::cout << "Tension at bottom: " << tensionBottom << " N\n";

    // ========================================================================
    // ORBITAL MECHANICS - MOTION AROUND THE EARTH
    // ========================================================================

    printSection("MOTION AROUND THE EARTH - ORBITAL MECHANICS");

    std::cout << "\n--- Low Earth Orbit (LEO) ---\n";
    double leoAltitude = 400000.0; // 400 km in meters

    std::cout << "Altitude: " << (leoAltitude / 1000.0) << " km\n";

    double orbitalVel = physics::orbital::calculateOrbitalVelocityFromAltitude(leoAltitude);
    double orbitalPeriod = physics::orbital::calculateOrbitalPeriodFromAltitude(leoAltitude);
    double orbitalAccel = physics::orbital::calculateOrbitalAcceleration(
        physics::orbital::constants::EARTH_RADIUS + leoAltitude);

    std::cout << "Orbital velocity: " << orbitalVel << " m/s (" <<
        physics::units::mpsToKmph(orbitalVel) << " km/h)\n";
    std::cout << "Orbital period: " << (orbitalPeriod / 60.0) << " minutes\n";
    std::cout << "Orbital acceleration: " << orbitalAccel << " m/s²\n";
    std::cout << "  (Compare to surface g = 9.81 m/s²)\n";

    std::cout << "\n--- Geostationary Orbit ---\n";
    double geoRadius = physics::orbital::calculateGeostationaryOrbitRadius();
    double geoAltitude = geoRadius - physics::orbital::constants::EARTH_RADIUS;

    std::cout << "Orbital radius: " << (geoRadius / 1000.0) << " km from Earth's center\n";
    std::cout << "Altitude: " << (geoAltitude / 1000.0) << " km above surface\n";

    double geoVel = physics::orbital::calculateOrbitalVelocity(geoRadius);
    std::cout << "Orbital velocity: " << geoVel << " m/s\n";
    std::cout << "Period: 24 hours (matches Earth's rotation)\n";

    std::cout << "\n--- Escape Velocity ---\n";
    double escapeVel = physics::orbital::calculateEscapeVelocityFromSurface();
    double surfaceOrbitalVel = physics::orbital::calculateOrbitalVelocity(
        physics::orbital::constants::EARTH_RADIUS);

    std::cout << "Escape velocity from Earth's surface: " << escapeVel << " m/s\n";
    std::cout << "  = " << physics::units::mpsToKmph(escapeVel) << " km/h\n";
    std::cout << "\nSurface orbital velocity (hypothetical): " << surfaceOrbitalVel << " m/s\n";
    std::cout << "Ratio (escape/orbital): " <<
        physics::orbital::escapeToOrbitalVelocityRatio(physics::orbital::constants::EARTH_RADIUS)
        << " = √2\n";

    std::cout << "\n--- Orbital Energy ---\n";
    double satelliteMass = 1000.0; // kg
    double satAltitude = 500000.0; // 500 km

    std::cout << "Satellite mass: " << satelliteMass << " kg\n";
    std::cout << "Altitude: " << (satAltitude / 1000.0) << " km\n";

    double orbRadius = physics::orbital::constants::EARTH_RADIUS + satAltitude;
    double ke = physics::orbital::calculateOrbitalKineticEnergy(satelliteMass, orbRadius);
    double pe = physics::orbital::calculateOrbitalPotentialEnergy(satelliteMass, orbRadius);
    double totalEnergy = physics::orbital::calculateTotalOrbitalEnergy(satelliteMass, orbRadius);

    std::cout << "Kinetic energy: " << (ke / 1e9) << " GJ\n";
    std::cout << "Potential energy: " << (pe / 1e9) << " GJ\n";
    std::cout << "Total energy: " << (totalEnergy / 1e9) << " GJ\n";
    std::cout << "Note: Total energy = -KE (virial theorem)\n";

    std::cout << "\n--- Moon's Orbit ---\n";
    double moonDistance = physics::orbital::constants::MOON_ORBITAL_RADIUS;
    double moonOrbitalVel = physics::orbital::calculateOrbitalVelocity(moonDistance);
    double moonPeriod = physics::orbital::calculateOrbitalPeriod(moonDistance);

    std::cout << "Moon's orbital radius: " << (moonDistance / 1e6) << " × 10^6 m\n";
    std::cout << "Moon's orbital velocity: " << moonOrbitalVel << " m/s\n";
    std::cout << "Moon's orbital period: " << (moonPeriod / 86400.0) << " days\n";
    std::cout << "  (Actual period ≈ 27.3 days)\n";

    printSection("END OF ADVANCED PHYSICS SHOWCASE");
    std::cout << "\nAll advanced calculations completed successfully!\n\n";

    return 0;
}
