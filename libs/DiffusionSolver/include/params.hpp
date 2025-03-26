#ifndef DIFFUSION_QUANTUM_PARAMS_HPP
#define DIFFUSION_QUANTUM_PARAMS_HPP

#include <functional>
#include <vector>
#include <cmath>

struct DiffusionQuantumParams {
    double d_tau = 0.001;        // time step value
    int total_time_steps = 1e5;   // total number of time steps valued d_tau
    int eq_time_step = 1e4;       // time step to average from
    int n0_walkers = 1000;       // beginning number of walkers alive, also target number of walkers
    int nmax_walkers = 1100;     // maximal number of walkers alive - size of allocated vector

    std::vector<int> save_hist_at = std::vector<int>({});
    double xmin = -5; // sampling minimum for visualisation
    double xmax = 5;  // sampling maximum for visualisation

    // std::function<double(double)> pot = [](double x){return (1. / 2. * 0.067 * std::pow(10. / 27211.6, 2) * std::pow(x, 2));}; // potential in 1D quantum dot
    std::function<double(double)> pot = [](double x){return (1. / 2. * 0.16 * std::pow(x, 2));}; // potential in 1D quantum dot

    int n_bins = 200; // number of bins used for generating wave function

    bool blocks_calibration = false;
    int n_block = pow(2, 15);

    std::vector<double> nodes = std::vector<double>({});
};

#endif
