# Physics Showcase - C++ Implementation

A comprehensive C++ library implementing fundamental physics concepts with well-documented standalone functions. This showcase demonstrates classical mechanics principles through clean, reusable code organized in header-only libraries.

## Overview

This project implements core physics functionality as standalone C++ functions, categorized into nine comprehensive modules:

### Basic Mechanics
1. **Newton's Laws of Motion** (`newton_laws.hpp`)
2. **Kinematics** - Motion with constant acceleration (`kinematics.hpp`)
3. **Dynamics** - Force causing rectilinear motion (`dynamics.hpp`)

### Advanced Topics
4. **Unit Conversions** (`units.hpp`) - CGS, SI, Imperial unit conversions
5. **Inclined Plane** (`inclined_plane.hpp`) - Motion on slopes with/without friction
6. **Energy & Momentum** (`energy_momentum.hpp`) - KE, PE, momentum relationships
7. **Projectile Motion** (`projectile.hpp`) - 2D motion under gravity
8. **Circular Motion** (`circular_motion.hpp`) - Uniform circular motion
9. **Orbital Mechanics** (`orbital.hpp`) - Satellite orbits and escape velocity

All functions are thoroughly documented with parameter descriptions, units, return values, and exception handling.

## Project Structure

```
physics_showcase/
├── include/
│   └── physics/
│       ├── newton_laws.hpp      # Newton's three laws of motion
│       ├── kinematics.hpp       # Motion with constant acceleration
│       ├── dynamics.hpp         # Force and motion integration
│       ├── units.hpp            # Unit conversions (CGS, SI, Imperial)
│       ├── inclined_plane.hpp   # Velocity at foot of inclined plane
│       ├── energy_momentum.hpp  # Energy and momentum comparisons
│       ├── projectile.hpp       # Projectile motion (2D)
│       ├── circular_motion.hpp  # Circular motion with constant speed
│       └── orbital.hpp          # Motion around Earth (orbital mechanics)
├── examples/
│   ├── main.cpp                 # Basic mechanics demonstrations
│   └── advanced_demo.cpp        # Advanced topics demonstrations
├── Makefile                     # Build system
├── CMakeLists.txt              # CMake configuration
├── LICENSE
└── README.md
```

## Features

### Newton's Laws of Motion (`newton_laws.hpp`)

#### First Law - Law of Inertia
- `isInEquilibrium()` - Check if object is in equilibrium
- `calculateNetForce()` - Calculate net force from multiple forces

#### Second Law - F = ma
- `calculateForce()` - Calculate force from mass and acceleration
- `calculateAcceleration()` - Calculate acceleration from force and mass
- `calculateMass()` - Calculate mass from force and acceleration
- `calculateWeight()` - Calculate gravitational force

#### Third Law - Action-Reaction
- `calculateReactionForce()` - Get reaction force for given action
- `verifyActionReactionPair()` - Verify if two forces form valid pair

### Kinematics (`kinematics.hpp`)

Motion in a straight line with constant acceleration using all kinematic equations:

#### First Equation: v = v₀ + at
- `calculateFinalVelocity()` - Calculate final velocity
- `calculateAccelerationFromVelocities()` - Calculate acceleration
- `calculateTimeFromVelocities()` - Calculate time required

#### Second Equation: s = v₀t + ½at²
- `calculateDisplacement()` - Calculate displacement
- `calculateAccelerationFromDisplacement()` - Calculate acceleration

#### Third Equation: v² = v₀² + 2as
- `calculateFinalVelocityFromDisplacement()` - Calculate final velocity
- `calculateAccelerationFromVelocitySquared()` - Calculate acceleration
- `calculateDisplacementFromVelocities()` - Calculate displacement

#### Fourth Equation: s = ((v + v₀) / 2) × t
- `calculateDisplacementFromAverageVelocity()` - Using average velocity
- `calculateAverageVelocity()` - Calculate average velocity
- `calculateTimeFromAverageVelocity()` - Calculate time

#### Utility Functions
- `calculateDistance()` - Calculate distance traveled (always positive)
- `calculateStoppingDistance()` - Distance required to stop
- `calculateStoppingTime()` - Time required to stop

### Dynamics (`dynamics.hpp`)

Integration of forces with motion - how forces cause acceleration and resulting motion:

#### Force-Acceleration Relationships
- `calculateNetForce()` - Sum multiple forces
- `calculateAccelerationFromForce()` - F = ma → a = F/m
- `calculateRequiredForce()` - Force needed for desired acceleration

#### Force-Motion Integration: Velocity
- `calculateFinalVelocityFromForce()` - Final velocity after force acts
- `calculateVelocityChange()` - Change in velocity (Δv)
- `calculateTimeForVelocityChange()` - Time for velocity change

#### Force-Motion Integration: Displacement
- `calculateDisplacementFromForce()` - Displacement when force acts
- `calculateFinalVelocityFromForceAndDisplacement()` - Final velocity from force and displacement

#### Stopping Problems
- `calculateBrakingForce()` - Force needed to stop in given distance
- `calculateStoppingDistanceFromForce()` - Distance to stop with given force
- `calculateStoppingTimeFromForce()` - Time to stop with given force

#### Friction
- `calculateFrictionForce()` - Frictional force (f = μN)
- `calculateAccelerationWithFriction()` - Net acceleration with friction
- `calculateMinimumForceToOvercomeFriction()` - Minimum force to start motion

#### Work and Power
- `calculateWork()` - Work done by force (W = F·s)
- `calculatePower()` - Power (rate of work, P = F·v)

### Unit Conversions (`units.hpp`)

Convert between different unit systems (CGS, SI, Imperial):

#### Force and Mass
- `newtonsToДynes()` / `dynesToNewtons()` - Force conversions (1 N = 10⁵ Dynes)
- `kilogramsToGrams()` / `gramsToKilograms()` - Mass conversions
- `calculateForceDynes()` - Calculate force in Dynes from grams and cm/s²

#### Length and Velocity
- `metersToCentimeters()` / `centimetersToMeters()` - Length conversions
- `metersToFeet()` / `feetToMeters()` - Metric to Imperial
- `mpsToKmph()` / `kmphToMps()` - Velocity conversions
- `mpsToMph()` / `mphToMps()` - m/s to miles per hour

#### Energy and Angles
- `joulesToErgs()` / `ergsToJoules()` - Energy conversions (1 J = 10⁷ ergs)
- `degreesToRadians()` / `radiansToDegrees()` - Angle conversions

### Inclined Plane (`inclined_plane.hpp`)

Analyze motion on inclined surfaces:

#### Forces on Incline
- `calculateParallelForce()` - Force component down the slope (mg sin θ)
- `calculateNormalForce()` - Force perpendicular to slope (mg cos θ)

#### Acceleration and Velocity
- `calculateAccelerationFrictionless()` - a = g sin θ
- `calculateAccelerationWithFriction()` - a = g(sin θ - μ cos θ)
- `calculateVelocityAtFootFrictionless()` - Speed at bottom (frictionless)
- `calculateVelocityAtFootWithFriction()` - Speed at bottom (with friction)
- `calculateVelocityAtFootFromHeight()` - Using energy: v = √(2gh)

#### Additional Analysis
- `calculateMinimumAngleToSlide()` - Minimum angle to overcome static friction
- `calculateTimeToSlide()` - Time to slide down
- `calculateWorkAgainstFriction()` - Energy lost to friction

### Energy and Momentum (`energy_momentum.hpp`)

Compare and analyze kinetic energy and momentum:

#### Kinetic Energy
- `calculateKineticEnergy()` - KE = (1/2)mv²
- `calculateVelocityFromKE()` - v = √(2KE/m)

#### Momentum
- `calculateMomentum()` - p = mv
- `calculateVelocityFromMomentum()` - v = p/m

#### Relationships
- `calculateKEFromMomentum()` - KE = p²/(2m)
- `calculateMomentumFromKE()` - p = √(2m·KE)
- `kineticEnergyRatioConstantMomentum()` - How mass affects KE for same p
- `momentumRatioConstantKE()` - How mass affects p for same KE

#### Potential Energy
- `calculatePotentialEnergy()` - PE = mgh
- `calculateSpringPotentialEnergy()` - PE = (1/2)kx²
- `velocityFromFall()` - Final velocity from height (energy conservation)

#### Impulse
- `calculateImpulse()` - Change in momentum (Δp = mΔv)
- `calculateChangeInKE()` - Change in kinetic energy
- `calculateAverageForceFromImpulse()` - F = Impulse/Δt

### Projectile Motion (`projectile.hpp`)

Analyze 2D motion under gravity:

#### Initial Components
- `calculateHorizontalVelocity()` - v₀ₓ = v₀ cos θ
- `calculateVerticalVelocity()` - v₀ᵧ = v₀ sin θ

#### Flight Parameters
- `calculateRange()` - Horizontal range: R = v₀² sin(2θ) / g
- `calculateMaxHeight()` - Maximum height: H = (v₀ sin θ)² / (2g)
- `calculateTimeOfFlight()` - Total flight time
- `calculateTimeToMaxHeight()` - Time to reach peak

#### Position and Velocity
- `calculateHorizontalPosition()` - x(t) = v₀ₓ × t
- `calculateVerticalPosition()` - y(t) = v₀ᵧt - (1/2)gt²
- `getVerticalVelocityAtTime()` - vᵧ(t) = v₀ᵧ - gt
- `calculateSpeedAtTime()` - |v(t)| = √(vₓ² + vᵧ²)

#### Trajectory Analysis
- `calculateHeightAtDistance()` - y(x) parabolic trajectory equation
- `getAngleForMaxRange()` - Returns 45° for maximum range
- `calculateMaximumRange()` - Maximum range at 45°
- `calculateLaunchAngleForRange()` - Find angle for desired range

### Circular Motion (`circular_motion.hpp`)

Motion in a circle with constant speed:

#### Centripetal Quantities
- `calculateCentripetalAcceleration()` - a_c = v²/r
- `calculateCentripetalForce()` - F_c = mv²/r
- `calculateVelocityFromAcceleration()` - v = √(a_c × r)

#### Angular Motion
- `calculateAngularVelocity()` - ω = v/r
- `calculateTangentialVelocity()` - v = ωr
- `calculateCentripetalAccelFromAngular()` - a_c = ω²r

#### Period and Frequency
- `calculatePeriod()` - T = 2πr/v
- `calculateFrequency()` - f = v/(2πr)
- `calculateAngularVelocityFromPeriod()` - ω = 2π/T
- `rpmToRadPerSec()` / `radPerSecToRpm()` - RPM conversions

#### Applications
- `calculateBankingAngle()` - Banked curve angle for frictionless turn
- `calculateMinVelocityTopOfLoop()` - Minimum speed at top of vertical loop
- `calculateTensionAtBottom()` / `calculateTensionAtTop()` - String tension in vertical circle

### Orbital Mechanics (`orbital.hpp`)

Motion around Earth and celestial bodies:

#### Gravitational Force
- `calculateGravitationalForce()` - F = GMm/r²
- `calculateGravitationalAcceleration()` - g = GM/r²

#### Orbital Velocity and Period
- `calculateOrbitalVelocity()` - v = √(GM/r)
- `calculateOrbitalVelocityFromAltitude()` - Velocity at height above surface
- `calculateOrbitalPeriod()` - T = 2π√(r³/GM) (Kepler's 3rd Law)
- `calculateGeostationaryOrbitRadius()` - Radius for 24-hour period

#### Escape Velocity
- `calculateEscapeVelocity()` - v_esc = √(2GM/r)
- `calculateEscapeVelocityFromSurface()` - Earth surface escape velocity
- `escapeToOrbitalVelocityRatio()` - Verify v_esc = √2 × v_orbital

#### Orbital Energy
- `calculateOrbitalPotentialEnergy()` - PE = -GMm/r
- `calculateOrbitalKineticEnergy()` - KE = GMm/(2r)
- `calculateTotalOrbitalEnergy()` - E = -GMm/(2r)

#### Special Orbits
- `calculateLEOVelocity()` - Low Earth Orbit parameters
- `calculateWeightInOrbit()` - Gravitational force in orbit

## Units

All functions use SI units:
- **Mass**: kilograms (kg)
- **Distance/Displacement**: meters (m)
- **Time**: seconds (s)
- **Velocity**: meters per second (m/s)
- **Acceleration**: meters per second squared (m/s²)
- **Force**: Newtons (N)
- **Work/Energy**: Joules (J)
- **Power**: Watts (W)

## Building and Running

### Prerequisites
- C++ compiler with C++11 support or later (g++, clang++, MSVC)
- Make (optional, for using Makefile)
- CMake 3.10+ (optional, for using CMake)

### Option 1: Direct Compilation

```bash
# Build basic demo
g++ -std=c++11 -I./include examples/main.cpp -o physics_demo
./physics_demo

# Build advanced demo
g++ -std=c++11 -I./include examples/advanced_demo.cpp -o advanced_demo
./advanced_demo
```

### Option 2: Using Make

```bash
# Build both demos
make

# Run basic demo (Newton's Laws, Kinematics, Dynamics)
make run

# Run advanced demo (Units, Inclined Plane, Energy, Projectiles, Orbits)
make run-advanced

# Run both demos sequentially
make run-all
```

### Option 3: Using CMake

```bash
mkdir build
cd build
cmake ..
make

# Run basic demo
./physics_demo

# Run advanced demo
./advanced_demo
```

## Usage Examples

### Example 1: Newton's Second Law

```cpp
#include "physics/newton_laws.hpp"

// Calculate force required to accelerate 10kg object at 5 m/s²
double mass = 10.0;           // kg
double acceleration = 5.0;    // m/s²
double force = physics::newton::calculateForce(mass, acceleration);
// Result: 50.0 N
```

### Example 2: Kinematics - Free Fall

```cpp
#include "physics/kinematics.hpp"

// Object dropped from rest, falling for 3 seconds
double initialVelocity = 0.0;     // m/s
double gravity = 9.81;             // m/s²
double time = 3.0;                 // s

double finalVelocity = physics::kinematics::calculateFinalVelocity(
    initialVelocity, gravity, time);
// Result: 29.43 m/s

double distance = physics::kinematics::calculateDisplacement(
    initialVelocity, gravity, time);
// Result: 44.145 m
```

### Example 3: Dynamics - Braking Force

```cpp
#include "physics/dynamics.hpp"

// Calculate braking force to stop 1500kg car in 40m
double mass = 1500.0;           // kg
double speed = 25.0;            // m/s (90 km/h)
double distance = 40.0;         // m

double brakingForce = physics::dynamics::calculateBrakingForce(
    mass, speed, distance);
// Result: -11718.75 N (negative = opposite to motion)
```

### Example 4: Motion with Friction

```cpp
#include "physics/dynamics.hpp"

// Box pushed across floor with friction
double mass = 20.0;                  // kg
double appliedForce = 100.0;         // N
double frictionCoefficient = 0.25;   // dimensionless

double acceleration = physics::dynamics::calculateAccelerationWithFriction(
    mass, appliedForce, frictionCoefficient);
// Result: 2.5475 m/s²
```

## Error Handling

All functions validate input parameters and throw `std::invalid_argument` exceptions for invalid inputs:
- Mass must be greater than zero
- Time must be non-negative
- Certain denominators must be non-zero
- Physical constraints (e.g., v² cannot be negative in real solutions)

Example:
```cpp
try {
    double force = physics::newton::calculateForce(-10.0, 5.0);
} catch (const std::invalid_argument& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    // Output: "Error: Mass must be greater than zero"
}
```

## Function Documentation

All functions include detailed Doxygen-style comments:
- **@brief**: Brief description of what the function does
- **@param**: Description of each parameter with units
- **@return**: Description of return value with units
- **@throws**: Exceptions that may be thrown

Example:
```cpp
/**
 * @brief Calculate force using Newton's Second Law: F = ma
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param acceleration Acceleration of the object (in m/s²)
 * @return Force acting on the object (in Newtons)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateForce(double mass, double acceleration);
```

## Design Principles

1. **Header-Only Library**: All functions are inline in headers for easy integration
2. **Namespace Organization**: Functions organized in logical namespaces:
   - `physics::newton` - Newton's laws
   - `physics::kinematics` - Motion equations
   - `physics::dynamics` - Force and motion
   - `physics::units` - Unit conversions
   - `physics::inclined_plane` - Inclined plane mechanics
   - `physics::energy_momentum` - Energy and momentum
   - `physics::projectile` - Projectile motion
   - `physics::circular_motion` - Circular motion
   - `physics::orbital` - Orbital mechanics
3. **Comprehensive Documentation**: Every function fully documented with units and constraints
4. **Error Handling**: Input validation with meaningful exception messages
5. **SI Units**: Consistent use of SI units throughout (with conversion utilities)
6. **Pure Functions**: No side effects, easy to test and reason about

## Current Features

This showcase currently implements:
- ✅ Newton's Three Laws of Motion
- ✅ Kinematics (constant acceleration)
- ✅ Dynamics (force causing motion)
- ✅ Unit Conversions (CGS, SI, Imperial)
- ✅ Inclined Plane Mechanics
- ✅ Energy and Momentum
- ✅ 2D Projectile Motion
- ✅ Uniform Circular Motion
- ✅ Orbital Mechanics (Earth satellites)

## Future Extensions

This showcase can be further extended with:
- 3D motion with vectors
- Non-uniform circular motion (angular acceleration)
- Rigid body dynamics and rotation
- Simple harmonic motion and oscillations
- Wave mechanics
- Fluid mechanics
- Thermodynamics and heat transfer
- Relativistic mechanics
- Quantum mechanics basics

## License

This project is licensed under the terms specified in the LICENSE file.

## Contributing

When adding new physics functions:
1. Follow the existing documentation style
2. Include parameter descriptions with units
3. Add input validation with appropriate exceptions
4. Use SI units consistently
5. Add usage examples to the main.cpp
6. Update this README with new functionality

## References

The physics equations implemented are based on standard classical mechanics:
- Newton's Laws of Motion
- Kinematic equations for constant acceleration
- Dynamics combining forces and motion
- Inclined plane mechanics
- Energy and momentum conservation
- Projectile motion (parabolic trajectories)
- Uniform circular motion and centripetal force
- Kepler's Laws of Planetary Motion
- Universal Law of Gravitation
- SI unit system and CGS unit system
