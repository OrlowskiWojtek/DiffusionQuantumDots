#include "DQMCConfig.h"
#include "Core/include/solver.hpp"

#include <iomanip>
#include <iostream>
#include <chrono>

#define OUTPUT_TEXT_WIDTH 56

int main() {

    std::cout << "|--------------------------------------------------------|\n";
    std::cout << "| DIFFUSION QUANTUM MONTE CARLO PROGRAM FOR QUANTUM DOTS |\n";
    std::cout << "|" << std::setw(OUTPUT_TEXT_WIDTH / 2)
              << "VERSION: " << DiffusionQuantumMC_VERSION_MAJOR << "."
              << DiffusionQuantumMC_VERSION_MINOR << std::setw(OUTPUT_TEXT_WIDTH / 2 - 2) << "|\n";
    std::cout << "|--------------------------------------------------------|\n" << std::endl;


    auto start = std::chrono::high_resolution_clock::now();
    DiffusionQuantumSolver solver;
    solver.solve();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "[PAR] Time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms\n";

    return 0;
}
