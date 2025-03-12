#ifndef DIFFUSION_QUANTUM_SOLVER_HPP
#define DIFFUSION_QUANTUM_SOLVER_HPP

#include "params.hpp"
#include "walkers.hpp"
#include "blocking.hpp"
#include <memory>

class DiffusionQuantumSolver{
public:
    DiffusionQuantumSolver();
    ~DiffusionQuantumSolver();
    
    void load_params();
    void solve();
private:
    
    std::shared_ptr<DiffusionWalkers> walkers;
    std::shared_ptr<EnergyBlockingAnalyzer> block_analyzer;

    DiffusionQuantumParams params;
    DiffusionQuantumResults final_results;

    void init();
    void diffuse();
    void branch();
    void accumulate();

    size_t save_counter;
};

#endif
