cmake_minimum_required(VERSION 3.27)
project(ECE1724 CXX)

# Generate compile_commands.json for IntelliSense to work.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Cache third-party source code for reuse across CMake builds.
set(CPM_SOURCE_CACHE ${PROJECT_SOURCE_DIR}/third-party)

# CMake package manager for downloading and configuring dependencies.
include(cmake/CPM.cmake)

CPMAddPackage("gh:catchorg/Catch2@3.7.1")
CPMAddPackage("gh:fmtlib/fmt#11.0.2")
CPMAddPackage("gh:nlohmann/json@3.11.3")
CPMAddPackage("gh:yhirose/cpp-httplib@0.18.1")

add_subdirectory(source)
