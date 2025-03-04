#ifndef DIFFUSION_QUANTUM_SOLVER_HPP
#define DIFFUSION_QUANTUM_SOLVER_HPP

#include "params.hpp"
#include "walkers.hpp"
#include <memory>

class DiffusionQuantumSolver{
public:
    DiffusionQuantumSolver();
    ~DiffusionQuantumSolver();
    
    void load_params();
    void solve();
private:
    
    std::shared_ptr<DiffusionWalkers> walkers;
    DiffusionQuantumParams params;

    void init();
    void diffuse();
    void branch();
    void accumulate();
};

#endif
