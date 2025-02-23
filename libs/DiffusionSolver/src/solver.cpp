#include "solver.hpp"

DiffusionQuantumSolver::DiffusionQuantumSolver(){}
DiffusionQuantumSolver::~DiffusionQuantumSolver(){}

// TODO: init also filenames and other solve related systems
void DiffusionQuantumSolver::init(){
    walkers.init_walkers(params);
}

void DiffusionQuantumSolver::solve(){
    init();

    for(int i = 0; i < params.total_time_steps; i++){
        diffuse();
        branch();
        accumulate();
    }
}

void DiffusionQuantumSolver::diffuse(){
    walkers.diffuse();
}

void DiffusionQuantumSolver::branch(){
    walkers.eval_p();
    walkers.branch();
}
