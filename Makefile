# Makefile for Physics Showcase

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -I./include

# Target executable
TARGET = physics_demo

# Source files
SOURCES = examples/main.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Build the demo
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)
	@echo "Physics demo build complete! Run with: ./$(TARGET)"

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run the demo
run: $(TARGET)
	./$(TARGET)

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)
	@echo "Clean complete!"

# Rebuild everything
rebuild: clean all

# Help message
help:
	@echo "Physics Showcase - Build System"
	@echo ""
	@echo "Available targets:"
	@echo "  make         - Build physics demo (default)"
	@echo "  make run     - Build and run physics demo"
	@echo "  make clean   - Remove build artifacts"
	@echo "  make rebuild - Clean and rebuild everything"
	@echo "  make help    - Show this help message"

.PHONY: all run clean rebuild help
