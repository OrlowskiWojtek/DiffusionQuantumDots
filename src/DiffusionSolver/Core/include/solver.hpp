#ifndef DIFFUSION_QUANTUM_SOLVER_HPP
#define DIFFUSION_QUANTUM_SOLVER_HPP

#include "Core/include/electrons.hpp"
#include "DiffusionParams/include/params.hpp"
#include "Core/include/blocking.hpp"
#include <memory>

/**
 * Diffusion Quantum Solver class header
 * 
 * This class implements the core Diffusion Quantum Monte Carlo algorithm
 * for 2D quantum dots with Fixed-Node Approximation. It orchestrates the
 * entire simulation process including diffusion, branching, and data collection.
 * 
 * The implementation follows the DMC algorithm described in the documentation:
 * 1. Initialization of walkers
 * 2. Equilibration phase (without data collection)
 * 3. Production phase with diffusion, branching, and data accumulation
 * 4. Results analysis and output
 */
class DiffusionQuantumSolver{
public:
    DiffusionQuantumSolver();
    ~DiffusionQuantumSolver();
    
    void load_params();
    void solve();

private:
    std::unique_ptr<DiffusionQuantumElectrons> electrons;
    std::unique_ptr<EnergyBlockingAnalyzer> block_analyzer;
    DiffusionQuantumParams* params;

    DiffusionQuantumResults final_results;

    void init();
    void diffuse();
    void branch();
    void accumulate();

    void check_saving(int iter_idx);
    size_t save_counter;
};

#endif
