#include "DQMCConfig.h"
#include "Core/include/solver.hpp"
#include "DiffusionParams/include/params.hpp"

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
    //DiffusionQuantumSolver solver;
    //solver.solve();
    
    DiffusionQuantumParams* p = DiffusionQuantumParams::getInstance();

    for(double a = 0; a < 0.1; a += 0.01){
        for(double b = 0; b < 4; b += 0.25){
            std::unique_ptr<DiffusionQuantumSolver> solver = std::make_unique<DiffusionQuantumSolver>();

            p->a = a;
            p->b = b;
            std::cout << a << "\t" << b << "\t";
            solver->solve();
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "[PAR] Time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms\n";

    return 0;
}
