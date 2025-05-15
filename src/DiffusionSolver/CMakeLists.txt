add_library(DiffusionSolver 
            Core/src/results.cpp
            Core/src/solver.cpp
            Core/src/walkers.cpp
            Core/src/blocking.cpp
            Core/src/electrons.cpp
            DiffusionParams/src/params.cpp
            DiffusionParams/src/harmonic_potential.cpp
            TrialFunctions/src/abstract_singlebody_orbital.cpp
            TrialFunctions/src/abstract_manybody_orbital.cpp
            TrialFunctions/src/harmonic_oscillator.cpp
            TrialFunctions/src/static_orbital.cpp
            TrialFunctions/src/jastrow_slater.cpp
            )

# Optimize this library

target_compile_options(DiffusionSolver PRIVATE
    "-O3;-Wall"
)

target_include_directories(DiffusionSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DiffusionSolver PUBLIC UnitHandler)
