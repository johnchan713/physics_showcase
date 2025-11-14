# Makefile for Physics Showcase

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -I./include

# Target executables
TARGET1 = physics_demo
TARGET2 = advanced_demo
TARGETS = $(TARGET1) $(TARGET2)

# Source files
SOURCES1 = examples/main.cpp
SOURCES2 = examples/advanced_demo.cpp

# Object files
OBJECTS1 = $(SOURCES1:.cpp=.o)
OBJECTS2 = $(SOURCES2:.cpp=.o)

# Default target - build both
all: $(TARGETS)

# Build the basic demo
$(TARGET1): $(OBJECTS1)
	$(CXX) $(CXXFLAGS) -o $(TARGET1) $(OBJECTS1)
	@echo "Basic demo build complete! Run with: ./$(TARGET1)"

# Build the advanced demo
$(TARGET2): $(OBJECTS2)
	$(CXX) $(CXXFLAGS) -o $(TARGET2) $(OBJECTS2)
	@echo "Advanced demo build complete! Run with: ./$(TARGET2)"

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run the basic demo
run: $(TARGET1)
	./$(TARGET1)

# Run the advanced demo
run-advanced: $(TARGET2)
	./$(TARGET2)

# Run both demos
run-all: $(TARGETS)
	@echo "Running basic demo..."
	./$(TARGET1)
	@echo "\nRunning advanced demo..."
	./$(TARGET2)

# Clean build artifacts
clean:
	rm -f $(OBJECTS1) $(OBJECTS2) $(TARGETS)
	@echo "Clean complete!"

# Rebuild everything
rebuild: clean all

# Help message
help:
	@echo "Physics Showcase - Build System"
	@echo ""
	@echo "Available targets:"
	@echo "  make               - Build both demos (default)"
	@echo "  make run           - Build and run basic demo"
	@echo "  make run-advanced  - Build and run advanced demo"
	@echo "  make run-all       - Build and run both demos"
	@echo "  make clean         - Remove build artifacts"
	@echo "  make rebuild       - Clean and rebuild everything"
	@echo "  make help          - Show this help message"

.PHONY: all run run-advanced run-all clean rebuild help
