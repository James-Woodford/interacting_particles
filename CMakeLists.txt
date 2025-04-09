cmake_minimum_required(VERSION 3.10)
project(InteractingParticles LANGUAGES C CXX)
# add_subdirectory(src)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(${PROJECT_SOURCE_DIR}/src/interacting_particles/include)

find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)
include_directories(/usr/include/python3.10)

file(GLOB_RECURSE SIM_SRC_FILES ${PROJECT_SOURCE_DIR}/src/interacting_particles/*.cpp)
add_library(interacting_particles ${SIM_SRC_FILES})
# log source files
message(STATUS "Source files: ${SIM_SRC_FILES}")

target_link_libraries(interacting_particles PUBLIC Eigen3::Eigen)


pybind11_add_module(interacting_particles_py MODULE
        ${PROJECT_SOURCE_DIR}/src/interacting_particles/pybind.cpp
)

target_link_libraries(interacting_particles_py PRIVATE interacting_particles Eigen3::Eigen)

set_target_properties(interacting_particles_py PROPERTIES
        OUTPUT_NAME "interacting_particles"
)