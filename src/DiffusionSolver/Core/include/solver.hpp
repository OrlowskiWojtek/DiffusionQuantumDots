#ifndef DIFFUSION_QUANTUM_SOLVER_HPP
#define DIFFUSION_QUANTUM_SOLVER_HPP

#include "Core/include/electrons.hpp"
#include "DiffusionParams/include/params.hpp"
#include "Core/include/blocking.hpp"
#include "include/visualiser.hpp"
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
 *
 * This class is just a container for actual algorithms used in \ref DiffusionQuantumElectrons class
 */
class DiffusionQuantumSolver{
public:
    DiffusionQuantumSolver();
    ~DiffusionQuantumSolver();
    
    void load_params();

    /**
     * Main solver method that executes the full DMC algorithm
     * This includes equilibration, post-equilibrium calculations, and results analysis
     */
    void solve();

private:
    // Manages the ensemble of electron walkers
    std::unique_ptr<DiffusionQuantumElectrons> electrons;

    // Handles statistical analysis of energy data using blocking technique
    // This is crucial for accurate error estimation in correlated Monte Carlo data
    std::unique_ptr<EnergyBlockingAnalyzer> block_analyzer;

    // Class performing plots at the end of simulation
    // Morphologica based
    std::unique_ptr<WalkersVisualiser> vis;

    // Singleton instance of simulation parameters
    DiffusionQuantumParams* params;

    void init();
    void diffuse();
    void branch();
    void accumulate();

    // Move walkers and check with metropolis sampling as in variational Monte Carlo
    // After initialization walkers should have distribution of |\psi|^2
    void initialize_distribution();

    // Few runs to sample VMC energy
    void sample_vmc_energy();

    // Finish initialization - calculate variational energy
    void finish_initialization();
};

#endif
