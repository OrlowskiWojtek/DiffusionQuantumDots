add_library(DiffusionSolver 
            src/results.cpp
            src/solver.cpp
            src/walkers.cpp
            src/blocking.cpp
            DiffusionParams/src/params.cpp
            DiffusionParams/src/harmonic_potential.cpp
            )

# Optimize this library

target_compile_options(DiffusionSolver PRIVATE
    "-O3;-Wall"
)

target_include_directories(DiffusionSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DiffusionSolver PUBLIC UnitHandler)
