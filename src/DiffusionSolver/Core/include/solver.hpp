#ifndef DIFFUSION_QUANTUM_SOLVER_HPP
#define DIFFUSION_QUANTUM_SOLVER_HPP

#include "Core/include/electrons.hpp"
#include "DiffusionParams/include/params.hpp"
#include "Core/include/blocking.hpp"
#include <memory>

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
