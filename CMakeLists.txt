cmake_minimum_required(VERSION 3.15)

project(KMC)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directory that contains header/include files
include_directories(include)

# Read all source files from src directory
file(GLOB SOURCES "src/*.cpp")

# Find dependencies
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package (Boost REQUIRED)

add_compile_options(-O3)

# Create executable
add_executable(KMC ${SOURCES})
target_link_libraries (KMC Eigen3::Eigen Boost::boost)


