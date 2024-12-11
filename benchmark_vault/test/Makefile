# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 -Itest/catch2 -Itest/cpp-httplib
LIBS = -lssl -lcrypto -lpthread

# Project files and directories
TARGET = vault_benchmark
SRCS = main.cc
DEPS_DIR = test

# URLs for dependencies
CATCH2_REPO = https://github.com/catchorg/Catch2.git
HTTPLIB_REPO = https://github.com/yhirose/cpp-httplib.git

# Default rule: compile the program
all: $(TARGET)

# Download and set up Catch2
$(DEPS_DIR)/catch2:
	git clone --depth 1 $(CATCH2_REPO) $(DEPS_DIR)/catch2

# Download and set up cpp-httplib
$(DEPS_DIR)/cpp-httplib:
	git clone --depth 1 $(HTTPLIB_REPO) $(DEPS_DIR)/cpp-httplib

# Compile the program
$(TARGET): $(SRCS) $(DEPS_DIR)/catch2 $(DEPS_DIR)/cpp-httplib
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS) $(LIBS)

# Clean up compiled files
clean:
	rm -f $(TARGET)

# Run the program
run: $(TARGET)
	./$(TARGET)
