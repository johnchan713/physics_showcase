#include <iostream>
#include <iomanip>
#include <vector>
#include "../include/physics/newton_laws.hpp"
#include "../include/physics/kinematics.hpp"
#include "../include/physics/dynamics.hpp"

// Helper function to print section headers
void printSection(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(70, '=') << "\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(2);

    std::cout << "PHYSICS SHOWCASE - C++ Implementation\n";
    std::cout << "Demonstrating Newton's Laws, Kinematics, and Dynamics\n";

    // ========================================================================
    // NEWTON'S LAWS OF MOTION
    // ========================================================================

    printSection("NEWTON'S LAWS OF MOTION");

    // First Law - Law of Inertia
    std::cout << "\n--- First Law: Law of Inertia ---\n";
    std::vector<double> forces1 = {10.0, -5.0, -5.0}; // Balanced forces
    std::vector<double> forces2 = {10.0, -5.0, -3.0}; // Unbalanced forces

    std::cout << "Forces on object 1: 10N, -5N, -5N\n";
    std::cout << "Net force: " << physics::newton::calculateNetForce(forces1) << " N\n";
    std::cout << "In equilibrium? " << (physics::newton::isInEquilibrium(forces1) ? "Yes" : "No") << "\n";

    std::cout << "\nForces on object 2: 10N, -5N, -3N\n";
    std::cout << "Net force: " << physics::newton::calculateNetForce(forces2) << " N\n";
    std::cout << "In equilibrium? " << (physics::newton::isInEquilibrium(forces2) ? "Yes" : "No") << "\n";

    // Second Law - F = ma
    std::cout << "\n--- Second Law: F = ma ---\n";
    double mass = 10.0;     // kg
    double accel = 5.0;     // m/s²
    double force = physics::newton::calculateForce(mass, accel);

    std::cout << "Mass: " << mass << " kg\n";
    std::cout << "Acceleration: " << accel << " m/s²\n";
    std::cout << "Force required: " << force << " N\n";

    std::cout << "\nGiven force of 50N and mass of 10kg:\n";
    std::cout << "Acceleration: " << physics::newton::calculateAcceleration(50.0, 10.0) << " m/s²\n";

    // Third Law - Action-Reaction
    std::cout << "\n--- Third Law: Action-Reaction ---\n";
    double actionForce = 100.0; // N
    double reactionForce = physics::newton::calculateReactionForce(actionForce);

    std::cout << "Action force: " << actionForce << " N\n";
    std::cout << "Reaction force: " << reactionForce << " N\n";
    std::cout << "Valid action-reaction pair? "
              << (physics::newton::verifyActionReactionPair(actionForce, reactionForce) ? "Yes" : "No") << "\n";

    // Weight calculation
    std::cout << "\n--- Weight Calculation (Application of F = ma) ---\n";
    double objectMass = 75.0; // kg
    double weight = physics::newton::calculateWeight(objectMass);
    std::cout << "Mass: " << objectMass << " kg\n";
    std::cout << "Weight on Earth: " << weight << " N\n";
    std::cout << "Weight on Moon (g=1.62 m/s²): "
              << physics::newton::calculateWeight(objectMass, 1.62) << " N\n";

    // ========================================================================
    // KINEMATICS - MOTION WITH CONSTANT ACCELERATION
    // ========================================================================

    printSection("KINEMATICS - MOTION WITH CONSTANT ACCELERATION");

    // Example 1: Car accelerating
    std::cout << "\n--- Example 1: Car Acceleration ---\n";
    double v0 = 0.0;        // m/s (starting from rest)
    double a = 3.0;         // m/s²
    double t = 10.0;        // s

    std::cout << "Initial velocity: " << v0 << " m/s\n";
    std::cout << "Acceleration: " << a << " m/s²\n";
    std::cout << "Time: " << t << " s\n";

    double vf = physics::kinematics::calculateFinalVelocity(v0, a, t);
    double displacement = physics::kinematics::calculateDisplacement(v0, a, t);

    std::cout << "Final velocity: " << vf << " m/s\n";
    std::cout << "Displacement: " << displacement << " m\n";

    // Example 2: Free fall
    std::cout << "\n--- Example 2: Free Fall ---\n";
    double initialVel = 0.0;    // m/s (dropped, not thrown)
    double gravity = 9.81;      // m/s²
    double fallTime = 3.0;      // s

    std::cout << "Object dropped from rest\n";
    std::cout << "Fall time: " << fallTime << " s\n";

    double fallVelocity = physics::kinematics::calculateFinalVelocity(initialVel, gravity, fallTime);
    double fallDistance = physics::kinematics::calculateDisplacement(initialVel, gravity, fallTime);

    std::cout << "Final velocity: " << fallVelocity << " m/s\n";
    std::cout << "Distance fallen: " << fallDistance << " m\n";

    // Example 3: Braking car
    std::cout << "\n--- Example 3: Braking Car ---\n";
    double carSpeed = 30.0;         // m/s (108 km/h)
    double deceleration = 6.0;      // m/s²

    std::cout << "Initial speed: " << carSpeed << " m/s (108 km/h)\n";
    std::cout << "Deceleration: " << deceleration << " m/s²\n";

    double stopDist = physics::kinematics::calculateStoppingDistance(carSpeed, deceleration);
    double stopTime = physics::kinematics::calculateStoppingTime(carSpeed, deceleration);

    std::cout << "Stopping distance: " << stopDist << " m\n";
    std::cout << "Stopping time: " << stopTime << " s\n";

    // Example 4: Using v² = v₀² + 2as
    std::cout << "\n--- Example 4: Velocity from Displacement ---\n";
    double initialSpeed = 10.0;     // m/s
    double acceleration = 2.0;      // m/s²
    double distance = 50.0;         // m

    std::cout << "Initial speed: " << initialSpeed << " m/s\n";
    std::cout << "Acceleration: " << acceleration << " m/s²\n";
    std::cout << "Distance traveled: " << distance << " m\n";

    double finalSpeed = physics::kinematics::calculateFinalVelocityFromDisplacement(
        initialSpeed, acceleration, distance);

    std::cout << "Final speed: " << finalSpeed << " m/s\n";

    // ========================================================================
    // DYNAMICS - FORCE CAUSING MOTION
    // ========================================================================

    printSection("DYNAMICS - FORCE CAUSING RECTILINEAR MOTION");

    // Example 1: Force causing acceleration
    std::cout << "\n--- Example 1: Force Applied to Object ---\n";
    double objMass = 50.0;          // kg
    double appliedForce = 200.0;    // N
    double timeApplied = 5.0;       // s
    double initialV = 0.0;          // m/s

    std::cout << "Object mass: " << objMass << " kg\n";
    std::cout << "Applied force: " << appliedForce << " N\n";
    std::cout << "Time force applied: " << timeApplied << " s\n";
    std::cout << "Initial velocity: " << initialV << " m/s\n";

    double objAccel = physics::dynamics::calculateAccelerationFromForce(appliedForce, objMass);
    double objFinalV = physics::dynamics::calculateFinalVelocityFromForce(
        objMass, appliedForce, initialV, timeApplied);
    double objDisplacement = physics::dynamics::calculateDisplacementFromForce(
        objMass, appliedForce, initialV, timeApplied);

    std::cout << "Resulting acceleration: " << objAccel << " m/s²\n";
    std::cout << "Final velocity: " << objFinalV << " m/s\n";
    std::cout << "Displacement: " << objDisplacement << " m\n";

    // Example 2: Multiple forces
    std::cout << "\n--- Example 2: Multiple Forces Acting on Object ---\n";
    std::vector<double> multipleForces = {100.0, -30.0, 50.0, -20.0}; // N
    double objectMass = 25.0; // kg

    std::cout << "Forces acting: 100N, -30N, 50N, -20N\n";
    std::cout << "Object mass: " << objectMass << " kg\n";

    double netF = physics::dynamics::calculateNetForce(multipleForces);
    double resultAccel = physics::dynamics::calculateAccelerationFromForce(netF, objectMass);

    std::cout << "Net force: " << netF << " N\n";
    std::cout << "Acceleration: " << resultAccel << " m/s²\n";

    // Example 3: Braking force
    std::cout << "\n--- Example 3: Braking Force Required ---\n";
    double vehicleMass = 1500.0;    // kg
    double vehicleSpeed = 25.0;     // m/s
    double brakingDist = 40.0;      // m

    std::cout << "Vehicle mass: " << vehicleMass << " kg\n";
    std::cout << "Initial speed: " << vehicleSpeed << " m/s (90 km/h)\n";
    std::cout << "Available braking distance: " << brakingDist << " m\n";

    double brakingF = physics::dynamics::calculateBrakingForce(vehicleMass, vehicleSpeed, brakingDist);
    double brakingT = physics::dynamics::calculateStoppingTimeFromForce(
        vehicleMass, vehicleSpeed, std::abs(brakingF));

    std::cout << "Required braking force: " << brakingF << " N\n";
    std::cout << "Braking time: " << brakingT << " s\n";

    // Example 4: Friction
    std::cout << "\n--- Example 4: Motion with Friction ---\n";
    double boxMass = 20.0;              // kg
    double pushForce = 100.0;           // N
    double frictionCoeff = 0.25;        // dimensionless

    std::cout << "Box mass: " << boxMass << " kg\n";
    std::cout << "Applied force: " << pushForce << " N\n";
    std::cout << "Coefficient of friction: " << frictionCoeff << "\n";

    double normalF = boxMass * 9.81;    // On horizontal surface
    double frictionF = physics::dynamics::calculateFrictionForce(frictionCoeff, normalF);
    double accelWithFriction = physics::dynamics::calculateAccelerationWithFriction(
        boxMass, pushForce, frictionCoeff);

    std::cout << "Normal force: " << normalF << " N\n";
    std::cout << "Friction force: " << frictionF << " N\n";
    std::cout << "Net acceleration: " << accelWithFriction << " m/s²\n";

    // Check minimum force to overcome static friction
    double staticCoeff = 0.4;
    double minForce = physics::dynamics::calculateMinimumForceToOvercomeFriction(
        boxMass, staticCoeff);
    std::cout << "\nWith static friction coefficient: " << staticCoeff << "\n";
    std::cout << "Minimum force to start motion: " << minForce << " N\n";

    // Example 5: Work and Power
    std::cout << "\n--- Example 5: Work and Power ---\n";
    double constantForce = 50.0;    // N
    double moveDistance = 10.0;     // m
    double velocity = 5.0;          // m/s

    std::cout << "Force applied: " << constantForce << " N\n";
    std::cout << "Displacement: " << moveDistance << " m\n";
    std::cout << "Velocity: " << velocity << " m/s\n";

    double work = physics::dynamics::calculateWork(constantForce, moveDistance);
    double power = physics::dynamics::calculatePower(constantForce, velocity);

    std::cout << "Work done: " << work << " J (Joules)\n";
    std::cout << "Power: " << power << " W (Watts)\n";

    // ========================================================================
    // INTEGRATED EXAMPLE
    // ========================================================================

    printSection("INTEGRATED EXAMPLE: Rocket Launch Simulation");

    std::cout << "\nA rocket with mass 1000 kg is launched vertically upward.\n";
    std::cout << "Engine provides thrust of 15000 N for 10 seconds.\n";

    double rocketMass = 1000.0;         // kg
    double thrust = 15000.0;            // N
    double g = 9.81;                    // m/s²
    double burnTime = 10.0;             // s

    // During powered ascent
    double weight = physics::newton::calculateWeight(rocketMass, g);
    double netForceUp = thrust - weight;

    std::cout << "\n--- Powered Ascent (0-10s) ---\n";
    std::cout << "Thrust force: " << thrust << " N (upward)\n";
    std::cout << "Weight: " << weight << " N (downward)\n";
    std::cout << "Net force: " << netForceUp << " N (upward)\n";

    double ascentAccel = physics::dynamics::calculateAccelerationFromForce(netForceUp, rocketMass);
    double velocityAtBurnout = physics::kinematics::calculateFinalVelocity(0.0, ascentAccel, burnTime);
    double altitudeAtBurnout = physics::kinematics::calculateDisplacement(0.0, ascentAccel, burnTime);

    std::cout << "Acceleration: " << ascentAccel << " m/s²\n";
    std::cout << "Velocity at engine cutoff: " << velocityAtBurnout << " m/s\n";
    std::cout << "Altitude at engine cutoff: " << altitudeAtBurnout << " m\n";

    // After engine cutoff (free fall with initial upward velocity)
    std::cout << "\n--- Coasting Phase (after 10s) ---\n";
    std::cout << "Engine off, only gravity acting (deceleration: " << g << " m/s²)\n";

    double additionalHeight = physics::kinematics::calculateDisplacementFromVelocities(
        velocityAtBurnout, 0.0, -g);
    double timeToApex = physics::kinematics::calculateTimeFromVelocities(
        velocityAtBurnout, 0.0, -g);
    double maxAltitude = altitudeAtBurnout + additionalHeight;

    std::cout << "Time to reach maximum altitude: " << timeToApex << " s (after burnout)\n";
    std::cout << "Additional height gained: " << additionalHeight << " m\n";
    std::cout << "Maximum altitude: " << maxAltitude << " m\n";
    std::cout << "Total time to max altitude: " << (burnTime + timeToApex) << " s\n";

    printSection("END OF PHYSICS SHOWCASE");
    std::cout << "\nAll calculations completed successfully!\n\n";

    return 0;
}
