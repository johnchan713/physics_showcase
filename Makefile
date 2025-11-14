# Makefile for Physics Showcase

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -I./include

# Target executables
TARGET1 = physics_demo
TARGET2 = advanced_demo
TARGET3 = scientific_demo
TARGETS = $(TARGET1) $(TARGET2) $(TARGET3)

# Source files
SOURCES1 = examples/main.cpp
SOURCES2 = examples/advanced_demo.cpp
SOURCES3 = examples/scientific_demo.cpp

# Object files
OBJECTS1 = $(SOURCES1:.cpp=.o)
OBJECTS2 = $(SOURCES2:.cpp=.o)
OBJECTS3 = $(SOURCES3:.cpp=.o)

# Default target - build all three
all: $(TARGETS)

# Build the basic demo
$(TARGET1): $(OBJECTS1)
	$(CXX) $(CXXFLAGS) -o $(TARGET1) $(OBJECTS1)
	@echo "Basic demo build complete! Run with: ./$(TARGET1)"

# Build the advanced demo
$(TARGET2): $(OBJECTS2)
	$(CXX) $(CXXFLAGS) -o $(TARGET2) $(OBJECTS2)
	@echo "Advanced demo build complete! Run with: ./$(TARGET2)"

# Build the scientific demo
$(TARGET3): $(OBJECTS3)
	$(CXX) $(CXXFLAGS) -o $(TARGET3) $(OBJECTS3)
	@echo "Scientific demo build complete! Run with: ./$(TARGET3)"

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run the basic demo
run: $(TARGET1)
	./$(TARGET1)

# Run the advanced demo
run-advanced: $(TARGET2)
	./$(TARGET2)

# Run the scientific demo
run-scientific: $(TARGET3)
	./$(TARGET3)

# Run all demos
run-all: $(TARGETS)
	@echo "Running basic demo..."
	./$(TARGET1)
	@echo "\nRunning advanced demo..."
	./$(TARGET2)
	@echo "\nRunning scientific demo..."
	./$(TARGET3)

# Clean build artifacts
clean:
	rm -f $(OBJECTS1) $(OBJECTS2) $(OBJECTS3) $(TARGETS)
	@echo "Clean complete!"

# Rebuild everything
rebuild: clean all

# Help message
help:
	@echo "Physics Showcase - Build System"
	@echo ""
	@echo "Available targets:"
	@echo "  make                 - Build all demos (default)"
	@echo "  make run             - Build and run basic demo"
	@echo "  make run-advanced    - Build and run advanced demo"
	@echo "  make run-scientific  - Build and run scientific demo"
	@echo "  make run-all         - Build and run all demos"
	@echo "  make clean           - Remove build artifacts"
	@echo "  make rebuild         - Clean and rebuild everything"
	@echo "  make help            - Show this help message"

.PHONY: all run run-advanced run-scientific run-all clean rebuild help
