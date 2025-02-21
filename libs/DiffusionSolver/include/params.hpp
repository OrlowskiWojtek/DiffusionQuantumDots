#ifndef DIFFUSION_QUANTUM_PARAMS_HPP
#define DIFFUSION_QUANTUM_PARAMS_HPP

#include <functional>

struct DiffusionQuantumParams {
    double d_tau;         // time step value
    int total_time_steps; // total number of time steps valued d_tau
    int n0_walkers;       // beginning number of walkers alive
    int nmax_walkers;     // total number of walkers alive

    double xmin; // sampling minimum for visualisation
    double xmax; // sampling maximum for visualisation

    std::function<double(double)> pot; // potential in 1D quantum dot
};

#endif
