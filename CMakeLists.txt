cmake_minimum_required(VERSION 3.22)

project(DiffusionQuantumMC VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(DQMCConfig.h.in DQMCConfig.h)

add_subdirectory(libs/DiffusionSolver)
add_subdirectory(libs/FileLoader)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(DiffusionQuantumMC DiffusionSolver FileLoader)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR})
