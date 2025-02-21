#ifndef DIFFUSION_QUANTUM_SOLVER_HPP
#define DIFFUSION_QUANTUM_SOLVER_HPP

#include "params.hpp"
#include "walkers.hpp"

class DiffusionQuantumSolver{
public:
    DiffusionQuantumSolver();
    ~DiffusionQuantumSolver();
    
    void load();
    void solve();
private:
    
    DiffusionWalkers walkers;

    void init();
    void diffuse();
    void branch();
    void accumulate();
};

#endif
