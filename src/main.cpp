#include "DQMCConfig.h"
#include "Core/include/solver.hpp"
#include "DiffusionParams/include/params.hpp"
#include "include/UnitHandler.hpp"

#include <iomanip>
#include <iostream>
#include <chrono>

#define OUTPUT_TEXT_WIDTH 56

/**
 * Main entry point for the Diffusion Quantum Monte Carlo (DQMC) program
 * 
 * This program implements the Diffusion Quantum Monte Carlo algorithm with 
 * Fixed-Node Approximation for calculating ground and excited states of 
 * 2D quantum dots in anisotropic harmonic potentials.
 * 
 * The implementation follows the theoretical framework described in
 * many DQMC textbooks (mostly Foulkes 2001), where the Schrödinger equation is solved in 
 * imaginary time using a population of walkers that diffuse and branch
 * according to the trial wavefunction and local energy.
 */
int main() {

    std::cout << "|--------------------------------------------------------|\n";
    std::cout << "| DIFFUSION QUANTUM MONTE CARLO PROGRAM FOR QUANTUM DOTS |\n";
    std::cout << "|" << std::setw(OUTPUT_TEXT_WIDTH / 2)
              << "VERSION: " << DiffusionQuantumMC_VERSION_MAJOR << "."
              << DiffusionQuantumMC_VERSION_MINOR << std::setw(OUTPUT_TEXT_WIDTH / 2 - 1) << "|\n";
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
