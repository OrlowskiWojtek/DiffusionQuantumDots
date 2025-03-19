#ifndef DIFFUSION_QUANTUM_PARAMS_HPP
#define DIFFUSION_QUANTUM_PARAMS_HPP

#include <functional>
#include <vector>
#include <cmath>

struct DiffusionQuantumParams {
    double d_tau = 0.0005;         // time step value
    int total_time_steps = 1e5; // total number of time steps valued d_tau
    int n0_walkers = 10000;       // beginning number of walkers alive
    int nmax_walkers = 11000;     // total number of walkers alive

    std::vector<int> save_hist_at = std::vector<int>({900000});
    double xmin = -5; // sampling minimum for visualisation
    double xmax = 5;  // sampling maximum for visualisation

    // std::function<double(double)> pot = [](double x){return (1. / 2. * 0.067 * std::pow(10. / 27211.6, 2) * std::pow(x, 2));}; // potential in 1D quantum dot
    std::function<double(double)> pot = [](double x){return (1. / 2. * std::pow(x, 2));}; // potential in 1D quantum dot
    
    int n_bins = 200; // number of bins used for generating wave function
    
    bool blocks_calibration = false;
    int n_block = pow(2, 10);

    void print_params();
};

#endif
