add_library(DiffusionSolver src/results.cpp src/solver.cpp src/walkers.cpp src/blocking.cpp)

# Optimize this library
target_compile_options(DiffusionSolver PRIVATE
    "-O3;-Wall"
)

target_include_directories(DiffusionSolver PUBLIC include)
