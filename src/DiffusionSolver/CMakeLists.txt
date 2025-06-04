add_library(DiffusionSolver 
            Core/src/results.cpp
            Core/src/solver.cpp
            Core/src/walkers.cpp
            Core/src/blocking.cpp
            Core/src/electrons.cpp
            Core/src/solver_context.cpp
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

option(PURE_DIFFUSION "Use pure diffusion algorithm without trial wavefunction" OFF)

if (PURE_DIFFUSION)
    message("omitting trial wavefunction version")
    target_compile_definitions(DiffusionSolver PRIVATE "PURE_DIFFUSION")
endif()

find_package(OpenMP REQUIRED)

target_include_directories(DiffusionSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DiffusionSolver PUBLIC UnitHandler OpenMP::OpenMP_CXX WalkersVisualiser)
