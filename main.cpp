#include "DQMCConfig.h"
#include "solver.hpp"

#include <iomanip>
#include <iostream>

#define OUTPUT_TEXT_WIDTH 56

int main() {

    std::cout << "|--------------------------------------------------------|\n";
    std::cout << "| DIFFUSION QUANTUM MONTE CARLO PROGRAM FOR QUANTUM DOTS |\n";
    std::cout << "|" << std::setw(OUTPUT_TEXT_WIDTH / 2)
              << "VERSION: " << DiffusionQuantumMC_VERSION_MAJOR << "."
              << DiffusionQuantumMC_VERSION_MINOR
              << std::setw(OUTPUT_TEXT_WIDTH / 2 - 2) << "|" << std::endl;

    DiffusionQuantumSolver solver;
    solver.solve();

    return 0;
}
