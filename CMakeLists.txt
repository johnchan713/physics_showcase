cmake_minimum_required(VERSION 3.10)

# Project name and version
project(PhysicsShowcase VERSION 1.0 LANGUAGES CXX)

# Specify C++ standard (updated to C++17 for novel theories)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(USE_FFTW "Enable FFTW library for production-quality FFT in Navier-Stokes solver" OFF)
option(BUILD_TESTS "Build test executables for novel theories" ON)

# Compiler warnings
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/new_theory)

# Find FFTW if requested
if(USE_FFTW)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(FFTW fftw3)
    endif()

    if(FFTW_FOUND)
        include_directories(${FFTW_INCLUDE_DIRS})
        link_directories(${FFTW_LIBRARY_DIRS})
        add_compile_definitions(USE_FFTW)
        message(STATUS "FFTW3 found - enabling production FFT mode for Navier-Stokes solver")
        message(STATUS "  Include dirs: ${FFTW_INCLUDE_DIRS}")
        message(STATUS "  Libraries: ${FFTW_LIBRARIES}")
        set(FFTW_AVAILABLE TRUE)
    else()
        message(WARNING "FFTW3 not found - Navier-Stokes solver will use conceptual mode")
        message(STATUS "  To enable FFTW: install libfftw3-dev and reconfigure with -DUSE_FFTW=ON")
        set(FFTW_AVAILABLE FALSE)
    endif()
endif()

# Add executables
add_executable(physics_demo
    examples/main.cpp
)

# Novel Theory Tests
if(BUILD_TESTS)
    add_executable(test_quantum_gravity tests/new_theory/test_quantum_gravity_phenomenology.cpp)
    add_executable(test_adscft tests/new_theory/test_ads_cft_applications.cpp)
    add_executable(test_emergent tests/new_theory/test_emergent_spacetime.cpp)
    add_executable(test_navier_stokes tests/new_theory/test_navier_stokes.cpp)
    add_executable(test_ns_solver tests/new_theory/test_navier_stokes_solver.cpp)

    # Link FFTW to Navier-Stokes solver if available
    if(USE_FFTW AND FFTW_AVAILABLE)
        target_link_libraries(test_ns_solver ${FFTW_LIBRARIES})
    endif()

    # Link math library
    target_link_libraries(test_quantum_gravity m)
    target_link_libraries(test_adscft m)
    target_link_libraries(test_emergent m)
    target_link_libraries(test_navier_stokes m)
    target_link_libraries(test_ns_solver m)

    # Enable CTest
    enable_testing()
    add_test(NAME QuantumGravity COMMAND test_quantum_gravity)
    add_test(NAME AdSCFT COMMAND test_adscft)
    add_test(NAME EmergentSpacetime COMMAND test_emergent)
    add_test(NAME NavierStokesRegularity COMMAND test_navier_stokes)
    add_test(NAME NavierStokesSolver COMMAND test_ns_solver)
endif()

# Installation rules (optional)
install(TARGETS physics_demo DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

# Print build information
message(STATUS "")
message(STATUS "========================================")
message(STATUS "Physics Showcase - Build Configuration")
message(STATUS "========================================")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Build tests: ${BUILD_TESTS}")
if(USE_FFTW)
    if(FFTW_AVAILABLE)
        message(STATUS "FFTW mode: ENABLED (Production)")
    else()
        message(STATUS "FFTW mode: REQUESTED but not found")
    endif()
else()
    message(STATUS "FFTW mode: Disabled (Conceptual)")
endif()
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "========================================")
message(STATUS "")
