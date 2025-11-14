#include <iostream>
#include <iomanip>
#include "../include/physics/thermodynamics.hpp"
#include "../include/physics/fluid_mechanics.hpp"
#include "../include/physics/elasticity.hpp"
#include "../include/physics/surface_tension.hpp"
#include "../include/physics/wave_mechanics.hpp"

// Helper function to print section headers
void printSection(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(70, '=') << "\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(2);

    std::cout << "PHYSICS SHOWCASE - SCIENTIFIC APPLICATIONS\n";
    std::cout << "Thermodynamics, Fluids, Elasticity, Surface Tension, Waves\n";

    // ========================================================================
    // THERMODYNAMICS
    // ========================================================================

    printSection("THERMODYNAMICS - GAS LAWS");

    std::cout << "\n--- Boyle's Law (Isothermal Process) ---\n";
    std::cout << "Compressing a gas at constant temperature\n";

    double initialP = 101325.0;    // Pa (1 atm)
    double initialV = 2.0;         // m³
    double finalV = 0.5;           // m³

    std::cout << "Initial: P = " << (initialP / 1000.0) << " kPa, V = " << initialV << " m³\n";
    std::cout << "Final volume: " << finalV << " m³\n";

    double finalP = physics::thermodynamics::boylesLaw(initialP, initialV, finalV);
    std::cout << "Final pressure: " << (finalP / 1000.0) << " kPa\n";
    std::cout << "Pressure increased by factor: " << (finalP / initialP) << "\n";

    std::cout << "\n--- Ideal Gas Law ---\n";
    double moles = 2.0;            // mol
    double temp = 300.0;           // K (≈27°C)
    double volume = 1.0;           // m³

    std::cout << "Amount: " << moles << " mol\n";
    std::cout << "Temperature: " << temp << " K (" <<
        physics::thermodynamics::kelvinToCelsius(temp) << " °C)\n";
    std::cout << "Volume: " << volume << " m³\n";

    double pressure = physics::thermodynamics::idealGasLawPressure(moles, volume, temp);
    std::cout << "Pressure: " << (pressure / 1000.0) << " kPa\n";

    std::cout << "\n--- Barometer Height Measurement ---\n";
    double mercuryHeight = 0.76; // m (76 cm)

    double atmPressure = physics::thermodynamics::barometerPressure(mercuryHeight);
    std::cout << "Mercury column height: " << mercuryHeight << " m\n";
    std::cout << "Atmospheric pressure: " << (atmPressure / 1000.0) << " kPa\n";
    std::cout << "  = " << (atmPressure / 101325.0) << " atm\n";

    std::cout << "\n--- RMS Velocity of Gas Molecules ---\n";
    double gasTemp = 300.0;        // K
    double nitrogenMass = 4.65e-26; // kg (N₂ molecule)

    double rmsVel = physics::thermodynamics::rmsVelocity(gasTemp, nitrogenMass);
    std::cout << "Temperature: " << gasTemp << " K\n";
    std::cout << "Nitrogen molecule (N₂)\n";
    std::cout << "RMS velocity: " << rmsVel << " m/s\n";

    // ========================================================================
    // FLUID MECHANICS
    // ========================================================================

    printSection("FLUID MECHANICS");

    std::cout << "\n--- Pressure at Depth ---\n";
    double depth = 100.0; // meters

    std::cout << "Water depth: " << depth << " m\n";
    double waterPressure = physics::fluid_mechanics::pressureAtDepth(depth);
    std::cout << "Total pressure: " << (waterPressure / 1000.0) << " kPa\n";
    std::cout << "  = " << (waterPressure / 101325.0) << " atm\n";

    double gaugePressure = physics::fluid_mechanics::gaugePressure(depth);
    std::cout << "Gauge pressure: " << (gaugePressure / 1000.0) << " kPa\n";

    std::cout << "\n--- Continuity Equation (Flow in Pipes) ---\n";
    double area1 = 0.01;   // m² (100 cm²)
    double vel1 = 2.0;     // m/s
    double area2 = 0.005;  // m² (50 cm²)

    std::cout << "Pipe 1: Area = " << (area1 * 10000) << " cm², Velocity = " << vel1 << " m/s\n";
    std::cout << "Pipe 2: Area = " << (area2 * 10000) << " cm²\n";

    double vel2 = physics::fluid_mechanics::continuityEquation(area1, vel1, area2);
    std::cout << "Velocity in narrower pipe: " << vel2 << " m/s\n";

    double flowRate = physics::fluid_mechanics::volumeFlowRate(area1, vel1);
    std::cout << "Volume flow rate: " << flowRate << " m³/s = " <<
        (flowRate * 1000.0) << " L/s\n";

    std::cout << "\n--- Torricelli's Theorem (Velocity of Efflux) ---\n";
    double tankDepth = 5.0; // m

    std::cout << "Water tank depth: " << tankDepth << " m\n";
    double effluxVel = physics::fluid_mechanics::effluxVelocity(tankDepth);
    std::cout << "Efflux velocity: " << effluxVel << " m/s\n";

    double holeArea = 0.001; // m² (10 cm²)
    double effluxRate = physics::fluid_mechanics::volumetricEffluxRate(tankDepth, holeArea);
    std::cout << "Flow rate (hole area " << (holeArea * 10000) << " cm²): " <<
        (effluxRate * 1000.0) << " L/s\n";

    std::cout << "\n--- Bernoulli's Equation (Jet Range) ---\n";
    double jetDepth = 10.0; // m
    double jetHeight = 2.0; // m

    std::cout << "Tank depth: " << jetDepth << " m\n";
    std::cout << "Orifice height above ground: " << jetHeight << " m\n";

    double jetRange = physics::fluid_mechanics::horizontalJetRange(jetDepth, jetHeight);
    std::cout << "Horizontal jet range: " << jetRange << " m\n";

    std::cout << "\n--- Pipe Friction (Darcy-Weisbach) ---\n";
    double pipeLength = 100.0;     // m
    double pipeDiameter = 0.1;     // m
    double pipeVelocity = 2.0;     // m/s
    double frictionFactor = 0.02;  // typical for turbulent flow

    std::cout << "Pipe: Length = " << pipeLength << " m, Diameter = " <<
        (pipeDiameter * 100) << " cm\n";
    std::cout << "Flow velocity: " << pipeVelocity << " m/s\n";
    std::cout << "Friction factor: " << frictionFactor << "\n";

    double pressureDrop = physics::fluid_mechanics::pipeFrictionPressureDrop(
        frictionFactor, pipeLength, pipeDiameter,
        physics::fluid_mechanics::constants::WATER_DENSITY, pipeVelocity);
    std::cout << "Pressure drop due to friction: " << (pressureDrop / 1000.0) << " kPa\n";

    double headLoss = physics::fluid_mechanics::headLoss(
        frictionFactor, pipeLength, pipeDiameter, pipeVelocity);
    std::cout << "Head loss: " << headLoss << " m\n";

    // ========================================================================
    // ELASTICITY
    // ========================================================================

    printSection("ELASTICITY AND MATERIAL DEFORMATION");

    std::cout << "\n--- Young's Modulus (Steel Wire Stretching) ---\n";
    double wireLength = 2.0;        // m
    double wireDiameter = 0.002;    // m (2 mm)
    double appliedForce = 1000.0;   // N

    std::cout << "Steel wire: Length = " << wireLength << " m, Diameter = " <<
        (wireDiameter * 1000) << " mm\n";
    std::cout << "Applied force: " << appliedForce << " N\n";

    double wireArea = M_PI * wireDiameter * wireDiameter / 4.0;
    double elongation = physics::elasticity::calculateElongation(
        appliedForce, wireLength, wireArea,
        physics::elasticity::constants::STEEL_YOUNGS_MODULUS);

    std::cout << "Elongation: " << (elongation * 1000.0) << " mm\n";

    double stress = physics::elasticity::calculateStress(appliedForce, wireArea);
    double strain = physics::elasticity::calculateStrain(elongation, wireLength);
    std::cout << "Stress: " << (stress / 1e6) << " MPa\n";
    std::cout << "Strain: " << (strain * 100.0) << " %\n";

    std::cout << "\n--- Beam Deflection (Simply Supported) ---\n";
    double beamLength = 4.0;        // m
    double beamWidth = 0.1;         // m
    double beamHeight = 0.2;        // m
    double centerLoad = 5000.0;     // N

    std::cout << "Beam: Length = " << beamLength << " m, " <<
        (beamWidth * 100) << " cm × " << (beamHeight * 100) << " cm\n";
    std::cout << "Center load: " << centerLoad << " N\n";

    double momentOfInertia = physics::elasticity::calculateRectangularMomentOfInertia(
        beamWidth, beamHeight);
    double deflection = physics::elasticity::calculateBeamDeflectionCenterLoad(
        centerLoad, beamLength, physics::elasticity::constants::STEEL_YOUNGS_MODULUS,
        momentOfInertia);

    std::cout << "Maximum deflection: " << (deflection * 1000.0) << " mm\n";

    std::cout << "\n--- Bulk Modulus (Water Compression) ---\n";
    double waterVolume = 1.0;       // m³
    double pressureIncrease = 1e6;  // Pa (10 atm)

    std::cout << "Initial volume: " << waterVolume << " m³\n";
    std::cout << "Pressure increase: " << (pressureIncrease / 1e6) << " MPa\n";

    double volChange = physics::elasticity::calculateVolumeChange(
        waterVolume, pressureIncrease,
        physics::elasticity::constants::WATER_BULK_MODULUS);

    std::cout << "Volume change: " << (volChange * 1e6) << " cm³\n";
    std::cout << "Relative change: " << (std::abs(volChange) / waterVolume * 100.0) << " %\n";

    std::cout << "\n--- Cantilever Beam Deflection ---\n";
    double cantileverLength = 1.0;  // m
    double endLoad = 100.0;         // N
    double cantileverDiam = 0.01;   // m (1 cm)

    std::cout << "Cantilever: Length = " << cantileverLength << " m, Diameter = " <<
        (cantileverDiam * 100) << " cm\n";
    std::cout << "End load: " << endLoad << " N\n";

    double circularI = physics::elasticity::calculateCircularMomentOfInertia(cantileverDiam);
    double cantileverDeflection = physics::elasticity::calculateCantileverDeflection(
        endLoad, cantileverLength,
        physics::elasticity::constants::STEEL_YOUNGS_MODULUS, circularI);

    std::cout << "Deflection at free end: " << (cantileverDeflection * 1000.0) << " mm\n";

    // ========================================================================
    // SURFACE TENSION
    // ========================================================================

    printSection("SURFACE TENSION");

    std::cout << "\n--- Capillary Rise in Glass Tube ---\n";
    double tubeRadius = 0.0005;     // m (0.5 mm)

    std::cout << "Glass tube radius: " << (tubeRadius * 1000.0) << " mm\n";

    double capillaryHeight = physics::surface_tension::calculateWaterCapillaryRise(tubeRadius);
    std::cout << "Water rises to: " << (capillaryHeight * 100.0) << " cm\n";

    std::cout << "\n--- Capillary Depression (Mercury) ---\n";
    double mercuryTubeRadius = 0.001; // m (1 mm)

    std::cout << "Glass tube radius: " << (mercuryTubeRadius * 1000.0) << " mm\n";

    double mercuryDepression = physics::surface_tension::calculateMercuryCapillaryDepression(
        mercuryTubeRadius);
    std::cout << "Mercury depresses to: " << (mercuryDepression * 100.0) << " cm\n";
    std::cout << "  (negative value indicates depression)\n";

    std::cout << "\n--- Pressure Inside Water Droplet ---\n";
    double dropletRadius = 0.001; // m (1 mm)

    std::cout << "Droplet radius: " << (dropletRadius * 1000.0) << " mm\n";

    double excessPressure = physics::surface_tension::calculateDropletPressure(
        physics::surface_tension::constants::WATER_SURFACE_TENSION, dropletRadius);
    std::cout << "Excess pressure inside: " << excessPressure << " Pa\n";

    std::cout << "\n--- Soap Bubble Pressure ---\n";
    double bubbleRadius = 0.05; // m (5 cm)

    std::cout << "Bubble radius: " << (bubbleRadius * 100.0) << " cm\n";

    double bubblePressure = physics::surface_tension::calculateBubblePressure(
        physics::surface_tension::constants::SOAP_SOLUTION_SURFACE_TENSION, bubbleRadius);
    std::cout << "Excess pressure: " << bubblePressure << " Pa\n";
    std::cout << "  (Factor of 2 higher than droplet due to two surfaces)\n";

    std::cout << "\n--- Two Connected Bubbles ---\n";
    double bubble1Radius = 0.02; // m (2 cm)
    double bubble2Radius = 0.05; // m (5 cm)

    std::cout << "Bubble 1 radius: " << (bubble1Radius * 100.0) << " cm\n";
    std::cout << "Bubble 2 radius: " << (bubble2Radius * 100.0) << " cm\n";

    double pressureDiff = physics::surface_tension::calculateBubblePressureDifference(
        physics::surface_tension::constants::SOAP_SOLUTION_SURFACE_TENSION,
        bubble1Radius, bubble2Radius);
    std::cout << "Pressure difference: " << pressureDiff << " Pa\n";
    std::cout << "Air flows from smaller to larger bubble\n";

    std::cout << "\n--- Surface Energy and Droplet Splitting ---\n";
    double largeDroplet = 0.01; // m (1 cm)
    int numDroplets = 8;

    std::cout << "Large droplet radius: " << (largeDroplet * 100.0) << " cm\n";
    std::cout << "Splits into: " << numDroplets << " equal droplets\n";

    double smallRadius = physics::surface_tension::calculateSplitDropletRadius(
        largeDroplet, numDroplets);
    double energyRequired = physics::surface_tension::calculateDropletSplittingEnergy(
        physics::surface_tension::constants::WATER_SURFACE_TENSION,
        largeDroplet, smallRadius, numDroplets);

    std::cout << "Each small droplet radius: " << (smallRadius * 100.0) << " cm\n";
    std::cout << "Energy required for splitting: " << (energyRequired * 1e6) << " μJ\n";

    // ========================================================================
    // WAVE MECHANICS
    // ========================================================================

    printSection("WAVE MECHANICS AND ACOUSTICS");

    std::cout << "\n--- Sound Wave Properties ---\n";
    double soundFreq = 440.0; // Hz (A4 note)
    double soundSpeed = physics::wave_mechanics::soundVelocityInAir(20.0); // 20°C

    std::cout << "Frequency: " << soundFreq << " Hz (A4 musical note)\n";
    std::cout << "Temperature: 20°C\n";
    std::cout << "Speed of sound: " << soundSpeed << " m/s\n";

    double wavelength = physics::wave_mechanics::calculateWavelength(soundSpeed, soundFreq);
    std::cout << "Wavelength: " << wavelength << " m\n";

    std::cout << "\n--- Sound Intensity and Decibels ---\n";
    double soundPower = 0.01; // W (10 mW)
    double distance = 10.0;   // m

    std::cout << "Sound source power: " << (soundPower * 1000.0) << " mW\n";
    std::cout << "Distance from source: " << distance << " m\n";

    double intensity = physics::wave_mechanics::sphericalWaveIntensity(soundPower, distance);
    double decibelLevel = physics::wave_mechanics::calculateSoundLevelDecibels(intensity);

    std::cout << "Sound intensity: " << intensity << " W/m²\n";
    std::cout << "Sound level: " << decibelLevel << " dB\n";

    std::cout << "\n--- Doppler Effect (Ambulance Siren) ---\n";
    double sirenFreq = 1000.0;      // Hz
    double ambulanceSpeed = 30.0;   // m/s

    std::cout << "Siren frequency: " << sirenFreq << " Hz\n";
    std::cout << "Ambulance speed: " << ambulanceSpeed << " m/s\n";
    std::cout << "Sound speed: " << soundSpeed << " m/s\n";

    double freqApproaching = physics::wave_mechanics::dopplerSourceApproaching(
        sirenFreq, soundSpeed, ambulanceSpeed);
    double freqReceding = physics::wave_mechanics::dopplerSourceReceding(
        sirenFreq, soundSpeed, ambulanceSpeed);

    std::cout << "Frequency heard (approaching): " << freqApproaching << " Hz\n";
    std::cout << "Frequency heard (receding): " << freqReceding << " Hz\n";
    std::cout << "Frequency shift: " << (freqApproaching - freqReceding) << " Hz\n";

    std::cout << "\n--- Vibrating String (Guitar String) ---\n";
    double stringLength = 0.65;     // m
    double stringTension = 100.0;   // N
    double stringMass = 0.005;      // kg

    std::cout << "String: Length = " << stringLength << " m, Mass = " <<
        (stringMass * 1000.0) << " g\n";
    std::cout << "Tension: " << stringTension << " N\n";

    double linearDensity = physics::wave_mechanics::calculateLinearDensity(
        stringMass, stringLength);
    double fundamentalFreq = physics::wave_mechanics::stringFundamentalFrequency(
        stringLength, stringTension, linearDensity);

    std::cout << "Linear density: " << (linearDensity * 1000.0) << " g/m\n";
    std::cout << "Fundamental frequency: " << fundamentalFreq << " Hz\n";

    std::cout << "\nFirst 5 harmonics:\n";
    for (int n = 1; n <= 5; n++) {
        double harmonic = physics::wave_mechanics::stringHarmonicFrequency(
            n, stringLength, stringTension, linearDensity);
        double wavelength = physics::wave_mechanics::stringHarmonicWavelength(
            stringLength, n);
        std::cout << "  n=" << n << ": " << harmonic << " Hz, λ = " <<
            wavelength << " m\n";
    }

    std::cout << "\n--- Velocity of Sound in Different Media ---\n";
    std::cout << "Speed of sound:\n";
    std::cout << "  Air (20°C): " <<
        physics::wave_mechanics::soundVelocityInAir(20.0) << " m/s\n";
    std::cout << "  Water: " <<
        physics::wave_mechanics::constants::SOUND_SPEED_WATER << " m/s\n";
    std::cout << "  Steel: " <<
        physics::wave_mechanics::constants::SOUND_SPEED_STEEL << " m/s\n";
    std::cout << "  Aluminum: " <<
        physics::wave_mechanics::constants::SOUND_SPEED_ALUMINUM << " m/s\n";

    std::cout << "\n--- Newton vs Laplace Formula for Sound in Air ---\n";
    double airDensity = 1.225;      // kg/m³ at STP
    double airPressure = 101325.0;  // Pa

    std::cout << "Air at STP: P = " << (airPressure / 1000.0) << " kPa, ρ = " <<
        airDensity << " kg/m³\n";

    // Newton's formula (incorrect - assumes isothermal)
    double newtonSpeed = physics::wave_mechanics::newtonsFormulaSound(airPressure, airDensity);
    // Laplace's correction (correct - assumes adiabatic)
    double laplaceSpeed = physics::wave_mechanics::laplaceFormulaSound(
        physics::wave_mechanics::constants::GAMMA_AIR, airPressure, airDensity);

    std::cout << "Newton's formula: " << newtonSpeed << " m/s (incorrect)\n";
    std::cout << "Laplace's correction: " << laplaceSpeed << " m/s (correct)\n";
    std::cout << "Ratio (Laplace/Newton): " << (laplaceSpeed / newtonSpeed) <<
        " = √γ = " << std::sqrt(physics::wave_mechanics::constants::GAMMA_AIR) << "\n";

    std::cout << "\n--- Resonance in Open Tube (Organ Pipe) ---\n";
    double pipeLength = 1.0; // m

    std::cout << "Open pipe length: " << pipeLength << " m\n";
    std::cout << "First 3 harmonics:\n";

    for (int n = 1; n <= 3; n++) {
        double resonantFreq = physics::wave_mechanics::openTubeResonance(
            n, pipeLength, soundSpeed);
        std::cout << "  n=" << n << ": " << resonantFreq << " Hz\n";
    }

    printSection("END OF SCIENTIFIC PHYSICS SHOWCASE");
    std::cout << "\nAll scientific calculations completed successfully!\n\n";

    return 0;
}
