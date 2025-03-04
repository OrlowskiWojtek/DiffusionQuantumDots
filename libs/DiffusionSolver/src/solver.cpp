#include "solver.hpp"

DiffusionQuantumSolver::DiffusionQuantumSolver() {}
DiffusionQuantumSolver::~DiffusionQuantumSolver() {}

// TODO: init also filenames and other solve related systems
void DiffusionQuantumSolver::init() { walkers.init_walkers(params); }

void DiffusionQuantumSolver::solve() {

    std::cout << "Starting initialization"
              << std::endl; // TODO: switch to formatted universal input / output

    init();

    std::cout << "Init finished" << std::endl;

    for (int i = 0; i < params.total_time_steps; i++) {
        diffuse();
        branch();
        accumulate();
    }
}

void DiffusionQuantumSolver::diffuse() { walkers.diffuse(); }

void DiffusionQuantumSolver::branch() {
    walkers.eval_p();
    walkers.branch();
}

void DiffusionQuantumSolver::accumulate() {}
