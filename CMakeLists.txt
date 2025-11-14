cmake_minimum_required(VERSION 3.10)

# Project name and version
project(PhysicsShowcase VERSION 1.0 LANGUAGES CXX)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add executables
add_executable(physics_demo
    examples/main.cpp
)

add_executable(advanced_demo
    examples/advanced_demo.cpp
)

add_executable(scientific_demo
    examples/scientific_demo.cpp
)

# Installation rules (optional)
install(TARGETS physics_demo advanced_demo scientific_demo DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

# Print build information
message(STATUS "")
message(STATUS "Physics Showcase - Build Configuration")
message(STATUS "========================================")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
