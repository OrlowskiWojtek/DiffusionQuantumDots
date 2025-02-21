#include "solver.hpp"

DiffusionQuantumSolver::DiffusionQuantumSolver(){}
DiffusionQuantumSolver::~DiffusionQuantumSolver(){}

void DiffusionQuantumSolver::init(){
    walkers.random_init();
}

