# Physics Showcase - C++ Implementation

A comprehensive C++ library implementing fundamental physics concepts with well-documented standalone functions. This showcase demonstrates classical mechanics principles through clean, reusable code organized in header-only libraries.

## Overview

This project implements core physics functionality as standalone C++ functions, categorized into seventeen comprehensive modules:

### Basic Mechanics
1. **Newton's Laws of Motion** (`newton_laws.hpp`)
2. **Kinematics** - Motion with constant acceleration (`kinematics.hpp`)
3. **Dynamics** - Force causing rectilinear motion (`dynamics.hpp`)

### Motion and Energy
4. **Unit Conversions** (`units.hpp`) - CGS, SI, Imperial conversions (Dyne-Gram relation)
5. **Inclined Plane** (`inclined_plane.hpp`) - Velocity at foot, friction effects
6. **Energy & Momentum** (`energy_momentum.hpp`) - KE/momentum compared
7. **Projectile Motion** (`projectile.hpp`) - 2D parabolic trajectories

### Circular and Oscillatory Motion
8. **Circular Motion** (`circular_motion.hpp`) - Centripetal force expressions
9. **Simple Harmonic Motion** (`harmonic_motion.hpp`) - SHM, pendulums, vibrating masses

### Rotational Mechanics
10. **Rotational Dynamics** (`rotational_dynamics.hpp`) - Torque, angular momentum, moments of inertia

### Gravitation and Orbits
11. **Orbital Mechanics** (`orbital.hpp`) - Satellite motion around Earth
12. **Universal Gravitation** (`gravitation.hpp`) - Apple and Moon connection, Kepler's Laws

### Thermodynamics and Gases
13. **Thermodynamics** (`thermodynamics.hpp`) - Boyle's Law, gas laws, pressure, RMS velocity

### Fluids
14. **Fluid Mechanics** (`fluid_mechanics.hpp`) - Steady flow, continuity, Bernoulli, pipe friction

### Materials and Elasticity
15. **Elasticity** (`elasticity.hpp`) - Young's modulus, beams, bulk modulus, stress-strain

### Surface Phenomena
16. **Surface Tension** (`surface_tension.hpp`) - Capillary rise, droplet pressure

### Waves and Acoustics
17. **Wave Mechanics** (`wave_mechanics.hpp`) - Sound, Doppler effect, string vibrations

All functions are thoroughly documented with parameter descriptions, units, return values, and exception handling.

## Project Structure

```
physics_showcase/
├── include/
│   └── physics/
│       ├── newton_laws.hpp         # Newton's three laws of motion
│       ├── kinematics.hpp          # Motion with constant acceleration
│       ├── dynamics.hpp            # Force and motion integration
│       ├── units.hpp               # Unit conversions (Dyne-Gram, CGS, SI)
│       ├── inclined_plane.hpp      # Velocity at foot of inclined plane
│       ├── energy_momentum.hpp     # Energy and momentum comparisons
│       ├── projectile.hpp          # Projectile motion (2D)
│       ├── circular_motion.hpp     # Circular motion + centripetal force expressions
│       ├── harmonic_motion.hpp     # Simple harmonic motion and pendulums
│       ├── rotational_dynamics.hpp # Torque, angular momentum, moment of inertia
│       ├── orbital.hpp             # Motion around Earth (satellites)
│       ├── gravitation.hpp         # Universal gravitation and Kepler's Laws
│       ├── thermodynamics.hpp      # Boyle's Law, gas laws, RMS velocity
│       ├── fluid_mechanics.hpp     # Flow, continuity, Bernoulli, pipe friction
│       ├── elasticity.hpp          # Young's modulus, beams, bulk modulus
│       ├── surface_tension.hpp     # Capillary rise, droplet/bubble pressure
│       └── wave_mechanics.hpp      # Sound, Doppler, string vibrations
├── examples/
│   ├── main.cpp                 # Basic mechanics demonstrations
│   ├── advanced_demo.cpp        # Advanced topics demonstrations
│   └── scientific_demo.cpp      # Scientific applications demonstrations
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

#### Other Expressions for Centripetal Force
- `calculateCentripetalForceAngular()` - F_c = mω²r
- `calculateCentripetalForceFromPeriod()` - F_c = 4π²mr/T²
- `calculateCentripetalForceFromFrequency()` - F_c = 4π²mf²r

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

### Simple Harmonic Motion (`harmonic_motion.hpp`)

Oscillatory motion with restoring force proportional to displacement:

#### Acceleration and Force
- `calculateAcceleration()` - a = -ω²x
- `calculateMaxAcceleration()` - a_max = ω²A
- `calculateRestoringForce()` - F = -kx (Hooke's Law)
- `calculateMaxForce()` - F_max = kA

#### Angular Frequency and Period
- `calculateAngularFrequency()` - ω = √(k/m)
- `calculatePeriod()` - T = 2π√(m/k)
- `calculateFrequency()` - f = (1/2π)√(k/m)

#### Energy of Vibrating Mass
- `calculateKineticEnergy()` - KE = (1/2)k(A² - x²)
- `calculatePotentialEnergy()` - PE = (1/2)kx²
- `calculateTotalEnergy()` - E_total = (1/2)kA²
- `verifyEnergyConservation()` - Check KE + PE = constant

#### Displacement and Velocity
- `calculateDisplacement()` - x(t) = A cos(ωt + φ)
- `calculateVelocity()` - v(t) = -Aω sin(ωt + φ)
- `calculateMaxVelocity()` - v_max = Aω
- `calculateVelocityAtDisplacement()` - v = ω√(A² - x²)

#### Simple Pendulum
- `calculatePendulumPeriod()` - T = 2π√(L/g)
- `calculatePendulumFrequency()` - f = (1/2π)√(g/L)
- `calculatePendulumLength()` - L for desired period
- `calculatePendulumMaxVelocity()` - v_max = θ₀√(gL)
- `calculatePendulumEnergy()` - E = (1/2)mgLθ₀²

### Rotational Dynamics (`rotational_dynamics.hpp`)

Rotational analogs of linear motion:

#### Angular Acceleration
- `calculateAngularAcceleration()` - α = Δω/Δt
- `calculateFinalAngularVelocity()` - ω_f = ω_i + αt
- `calculateAngularDisplacement()` - θ = ω₀t + (1/2)αt²
- `calculateTangentialAcceleration()` - a_t = rα

#### Torque
- `calculateTorque()` - τ = rF (perpendicular force)
- `calculateTorqueWithAngle()` - τ = rF sin(θ)

#### Angular Acceleration Caused by Torque
- `calculateAngularAccelFromTorque()` - α = τ/I
- `calculateRequiredTorque()` - τ = Iα

#### Angular Momentum
- `calculateAngularMomentum()` - L = Iω
- `calculatePointMassAngularMomentum()` - L = mvr
- `calculateAngularMomentumChange()` - ΔL = τ·Δt

#### Kinetic Energy of Rotating Body
- `calculateRotationalKE()` - KE = (1/2)Iω²
- `calculateAngularVelocityFromKE()` - ω = √(2KE/I)

#### Formulas for Moment of Inertia
- `momentOfInertiaPointMass()` - I = mr²
- `momentOfInertiaRodCenter()` - I = (1/12)ML² (rod about center)
- `momentOfInertiaRodEnd()` - I = (1/3)ML² (rod about end)
- `momentOfInertiaCylinder()` - I = (1/2)MR² (solid cylinder)
- `momentOfInertiaHollowCylinder()` - I = MR²
- `momentOfInertiaSolidSphere()` - I = (2/5)MR²
- `momentOfInertiaHollowSphere()` - I = (2/3)MR²
- `momentOfInertiaDisk()` - I = (1/2)MR²
- `momentOfInertiaHoop()` - I = MR²

#### Moment of Inertia about Parallel Axis
- `parallelAxisTheorem()` - I = I_cm + Md²

#### Compound Pendulum
- `calculateCompoundPendulumPeriod()` - T = 2π√(I/(mgd))
- `calculateEquivalentLength()` - L_eq = I/(md)

#### Center of Percussion
- `calculateCenterOfPercussion()` - Sweet spot location
- `calculateRodCenterOfPercussion()` - For uniform rod: (2/3)L

### Universal Gravitation (`gravitation.hpp`)

Newton's synthesis of terrestrial and celestial mechanics:

#### Universal Gravitation
- `universalGravitationForce()` - F = Gm₁m₂/r²
- `gravitationalFieldStrength()` - g = GM/r²

#### Moon's Motions Connected with Fall of Apple
- `calculateSurfaceGravity()` - g at Earth's surface
- `calculateMoonCentripetalAccel()` - Moon's orbital acceleration
- `verifyInverseSquareLaw()` - g_surface/g_moon = (r_moon/R_earth)²
- `calculateDistanceRatio()` - r_moon/R_earth ≈ 60

#### Mass of the Earth
- `calculateEarthMass()` - From surface gravity: M = gR²/G
- `calculateEarthMassFromMoon()` - From Moon's orbit
- `calculateCentralMassFromOrbit()` - General formula for any satellite

#### Significance of Kepler's Third Law
- `calculatePeriodKeplerThird()` - T² = (4π²/GM)r³
- `calculateKeplerConstant()` - K = 4π²/(GM), same for all satellites
- `verifyKeplerThirdLaw()` - Check if orbit obeys law
- `compareSatellitePeriods()` - T₁/T₂ = (r₁/r₂)^(3/2)

#### Gravitational Potential Energy
- `gravitationalPotentialEnergy()` - U = -GMm/r
- `calculateEscapeVelocity()` - v_esc = √(2GM/R)

#### Binary Systems
- `calculateReducedMass()` - μ = m₁m₂/(m₁ + m₂)
- `calculateCMDistance()` - Distance from center of mass

### Thermodynamics (`thermodynamics.hpp`)

Gas laws, pressure, and molecular velocities:

#### Boyle's Law and Variations
- `boylesLaw()` - P₁V₁ = P₂V₂ (isothermal)
- `charlesLaw()` - V₁/T₁ = V₂/T₂ (isobaric)
- `gayLussacsLaw()` - P₁/T₁ = P₂/T₂ (isochoric)
- `combinedGasLaw()` - P₁V₁/T₁ = P₂V₂/T₂

#### Ideal Gas Law
- `idealGasLawPressure()` / `idealGasLawVolume()` / `idealGasLawTemperature()` - PV = nRT
- `idealGasLawMoles()` - Calculate number of moles

#### Pressure Measurements
- `barometerPressure()` - Pressure from mercury column height
- `altitudeFromPressure()` - Altitude from atmospheric pressure

#### Kinetic Theory of Gases
- `rmsVelocity()` - Root-mean-square molecular velocity v_rms = √(3kT/m)
- `averageKineticEnergy()` - Average KE of gas molecules
- `gasPressureFromKineticTheory()` - P = (1/3)ρv²_rms

#### Elasticity of Gases
- `isothermalBulkModulus()` - K_T = P (isothermal compression)
- `adiabaticBulkModulus()` - K_S = γP (adiabatic compression)

### Fluid Mechanics (`fluid_mechanics.hpp`)

Fluid statics, dynamics, and flow in pipes:

#### Fluid Statics
- `pressureAtDepth()` - P = P₀ + ρgh
- `gaugePressure()` - Pressure relative to atmospheric
- `absolutePressure()` - Total pressure including atmospheric

#### Continuity Equation
- `continuityEquation()` - A₁v₁ = A₂v₂ (mass conservation)
- `volumeFlowRate()` - Q = Av
- `massFlowRate()` - ṁ = ρAv

#### Bernoulli's Equation
- `bernoulliPressure()` - P + ½ρv² + ρgh = constant
- `bernoulliVelocity()` - Calculate velocity from pressure difference
- `dynamicPressure()` - q = ½ρv²

#### Torricelli's Theorem (Efflux)
- `effluxVelocity()` - v = √(2gh) (velocity of jet from orifice)
- `volumetricEffluxRate()` - Flow rate through orifice
- `effluxGasIsothermal()` / `effluxGasAdiabatic()` - Gas efflux velocities

#### Jet Motion
- `horizontalJetRange()` - Range of horizontal jet
- `jetVelocityFromRange()` - Calculate initial velocity from range

#### Pipe Friction
- `pipeFrictionPressureDrop()` - Darcy-Weisbach equation ΔP = f(L/D)(½ρv²)
- `reynoldsNumber()` - Re = vD/ν (flow characterization)
- `headLoss()` - Energy loss due to friction
- `powerLoss()` - Power dissipated by friction

### Elasticity (`elasticity.hpp`)

Material deformation under stress:

#### Stress and Strain
- `calculateStress()` - σ = F/A
- `calculateStrain()` - ε = ΔL/L₀

#### Young's Modulus (Tensile/Compressive)
- `calculateYoungsModulus()` - E = σ/ε
- `calculateElongation()` - ΔL = FL/(AE)
- `calculateStressFromStrain()` - σ = Eε (Hooke's Law)

#### Bulk Modulus (Volume Elasticity)
- `calculateBulkModulus()` - K = -ΔP/(ΔV/V₀)
- `calculateVolumeChange()` - ΔV from applied pressure
- `calculateCompressibility()` - β = 1/K
- `calculateVolumetricStrain()` - ε_v = ΔV/V₀

#### Shear Modulus (Rigidity)
- `calculateShearModulus()` - G = τ/γ
- `calculateShearStress()` - τ = F/A
- `calculateShearStrain()` - γ = Δx/h

#### Poisson's Ratio
- `calculatePoissonsRatio()` - ν = -ε_transverse/ε_axial
- `calculateTransverseStrain()` - Lateral strain from axial strain
- `calculateShearModulusFromYoungs()` - G = E/(2(1+ν))
- `calculateBulkModulusFromYoungs()` - K = E/(3(1-2ν))

#### Beam Bending and Deflection
- `calculateBeamDeflectionCenterLoad()` - δ = FL³/(48EI) (simply supported)
- `calculateCantileverDeflection()` - δ = FL³/(3EI)
- `calculateBeamDeflectionUniformLoad()` - Distributed load deflection
- `calculateRectangularMomentOfInertia()` - I = bh³/12
- `calculateCircularMomentOfInertia()` - I = πd⁴/64
- `calculateBendingStress()` - σ = My/I

#### Elastic Energy
- `calculateElasticEnergy()` - U = ½k(ΔL)²
- `calculateEnergyDensity()` - u = σ²/(2E)

### Surface Tension (`surface_tension.hpp`)

Surface phenomena in liquids:

#### Surface Tension Fundamentals
- `calculateSurfaceTensionForce()` - F = γL
- `calculateSurfaceEnergy()` - E = γA
- `calculateWorkInStretching()` - Work to stretch liquid film

#### Pressure Due to Curved Surfaces
- `calculateDropletPressure()` - ΔP = 2γ/r (spherical droplet)
- `calculateBubblePressure()` - ΔP = 4γ/r (soap bubble, two surfaces)
- `calculateCylindricalPressure()` - ΔP = γ/r (liquid jet)
- `calculateYoungLaplacePressure()` - ΔP = γ(1/R₁ + 1/R₂)

#### Capillary Rise
- `calculateCapillaryRise()` - h = 2γcosθ/(ρgr)
- `calculateWaterCapillaryRise()` - Rise in glass tube (perfect wetting)
- `calculateMercuryCapillaryDepression()` - Depression in glass tube
- `calculateTubeRadiusFromRise()` - Determine radius from observed rise
- `calculateCapillaryRiseBetweenPlates()` - Rise between parallel plates

#### Meniscus and Contact Angle
- `calculateCapillaryLength()` - a = √(γ/(ρg))
- `calculateContactAngle()` - Young's equation cosθ = (γ_SV - γ_SL)/γ_LV

#### Droplet and Bubble Dynamics
- `calculateDropletRadius()` - Radius from volume
- `calculateDropletSplittingEnergy()` - Energy to split droplet
- `calculateSplitDropletRadius()` - Radius after splitting
- `calculateBubblePressureDifference()` - ΔP between connected bubbles

#### Applications
- `calculateMaxSupportedWeight()` - Force supported by surface tension
- `calculateRingDetachmentForce()` - Force to pull ring from surface
- `calculateSurfaceTensionFromRise()` - Experimental determination

### Wave Mechanics (`wave_mechanics.hpp`)

Wave phenomena, sound, and vibrations:

#### Wave Fundamentals
- `calculateWavelength()` - λ = v/f
- `calculateFrequency()` - f = v/λ
- `calculateWaveVelocity()` - v = fλ
- `calculateAngularFrequency()` - ω = 2πf
- `calculateWaveNumber()` - k = 2π/λ
- `calculatePeriod()` - T = 1/f

#### Velocity of Sound
- `newtonsFormulaSound()` - v = √(E/ρ) (for solids/liquids)
- `laplaceFormulaSound()` - v = √(γP/ρ) (for gases, Laplace correction)
- `soundVelocityFromTemperature()` - v = √(γRT/M)
- `soundVelocityInAir()` - v ≈ 331.3 + 0.606T (empirical)
- `soundVelocityInSolid()` / `soundVelocityInLiquid()` - Speed in different media

#### Sound Intensity
- `calculateSoundIntensity()` - I = P/A
- `calculateSoundLevelDecibels()` - β = 10log₁₀(I/I₀)
- `intensityFromDecibels()` - Convert dB to W/m²
- `intensityFromPressure()` - I = ΔP²/(2ρv)
- `sphericalWaveIntensity()` - I = P/(4πr²)

#### Doppler Effect
- `dopplerFrequency()` - f' = f(v + v_o)/(v - v_s)
- `dopplerSourceApproaching()` / `dopplerSourceReceding()` - Moving source
- `dopplerObserverApproaching()` - Moving observer
- `calculateBeatFrequency()` - f_beat = |f₁ - f₂|

#### String Vibrations
- `stringWaveVelocity()` - v = √(T/μ)
- `calculateLinearDensity()` - μ = m/L
- `stringFundamentalFrequency()` - f₁ = v/(2L)
- `stringHarmonicFrequency()` - f_n = nf₁
- `stringHarmonicWavelength()` - λ_n = 2L/n
- `calculateRequiredTension()` - Tension for desired frequency

#### Standing Waves and Resonance
- `openTubeResonance()` - f_n = nv/(2L) (both ends open)
- `closedTubeResonance()` - f_n = nv/(4L) (one end closed, n odd)

#### Wave Energy
- `waveEnergyDensity()` - u = ½ρω²A²
- `wavePower()` - Power transmitted by wave

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

# Build scientific demo
g++ -std=c++11 -I./include examples/scientific_demo.cpp -o scientific_demo
./scientific_demo
```

### Option 2: Using Make

```bash
# Build all demos
make

# Run basic demo (Newton's Laws, Kinematics, Dynamics)
make run

# Run advanced demo (Units, Inclined Plane, Energy, Projectiles, Orbits, Gravitation)
make run-advanced

# Run scientific demo (Thermodynamics, Fluids, Elasticity, Surface Tension, Waves)
make run-scientific

# Run all demos sequentially
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

# Run scientific demo
./scientific_demo
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
   - `physics::harmonic_motion` - Simple harmonic motion
   - `physics::rotational_dynamics` - Rotational mechanics
   - `physics::orbital` - Orbital mechanics
   - `physics::gravitation` - Universal gravitation
   - `physics::thermodynamics` - Gas laws and thermodynamics
   - `physics::fluid_mechanics` - Fluid dynamics
   - `physics::elasticity` - Material deformation
   - `physics::surface_tension` - Surface phenomena
   - `physics::wave_mechanics` - Waves and acoustics
3. **Comprehensive Documentation**: Every function fully documented with units and constraints
4. **Error Handling**: Input validation with meaningful exception messages
5. **SI Units**: Consistent use of SI units throughout (with conversion utilities)
6. **Pure Functions**: No side effects, easy to test and reason about

## Current Features

This showcase currently implements:

**Basic Mechanics:**
- ✅ Newton's Three Laws of Motion
- ✅ Kinematics (constant acceleration)
- ✅ Dynamics (force causing motion)
- ✅ Unit Conversions (Dyne-Gram relationship, CGS, SI, Imperial)

**Motion and Energy:**
- ✅ Inclined Plane Mechanics (velocity at foot)
- ✅ Energy and Momentum (KE/momentum compared)
- ✅ 2D Projectile Motion

**Circular and Oscillatory Motion:**
- ✅ Uniform Circular Motion (centripetal force expressions)
- ✅ Simple Harmonic Motion (SHM, pendulums, vibrating masses)

**Rotational Mechanics:**
- ✅ Rotational Dynamics (torque, angular momentum, moment of inertia)

**Gravitation:**
- ✅ Orbital Mechanics (Earth satellites)
- ✅ Universal Gravitation (Apple-Moon connection, Kepler's Laws)

**Thermodynamics:**
- ✅ Gas Laws (Boyle's, Charles's, Gay-Lussac's, Ideal Gas Law)
- ✅ Pressure Measurements (barometer, altitude)
- ✅ Kinetic Theory (RMS velocity, molecular KE)
- ✅ Elasticity of Gases

**Fluid Mechanics:**
- ✅ Fluid Statics (pressure at depth)
- ✅ Continuity Equation and Flow Rates
- ✅ Bernoulli's Equation
- ✅ Torricelli's Theorem (efflux velocity)
- ✅ Pipe Friction (Darcy-Weisbach)

**Elasticity:**
- ✅ Young's Modulus (stress-strain)
- ✅ Bulk Modulus (volume elasticity)
- ✅ Shear Modulus and Poisson's Ratio
- ✅ Beam Bending and Deflection

**Surface Tension:**
- ✅ Capillary Rise and Depression
- ✅ Droplet and Bubble Pressure
- ✅ Young-Laplace Equation
- ✅ Surface Energy

**Wave Mechanics:**
- ✅ Sound Waves (Newton's formula, Laplace correction)
- ✅ Sound Intensity and Decibels
- ✅ Doppler Effect
- ✅ String Vibrations and Harmonics
- ✅ Standing Waves and Resonance

## Future Extensions

This showcase can be further extended with:
- 3D motion with vectors
- Non-uniform circular motion
- Heat transfer and calorimetry
- Electromagnetic waves
- Optics (reflection, refraction, interference)
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
