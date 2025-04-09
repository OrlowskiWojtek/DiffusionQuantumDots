#ifndef DIFFUSION_QUANTUM_PARAMS_HPP
#define DIFFUSION_QUANTUM_PARAMS_HPP

#include "include/walkers_struct.hpp"
#include <cmath>
#include <functional>
#include <vector>

// Class implements Singleton design pattern
class DiffusionQuantumParams {
private:
    static DiffusionQuantumParams *instance;
    DiffusionQuantumParams() {}

public:
    DiffusionQuantumParams(const DiffusionQuantumParams &) = delete;
    DiffusionQuantumParams &operator=(const DiffusionQuantumParams &) = delete;

    static DiffusionQuantumParams *getInstance() {
        if (!instance) {
            instance = new DiffusionQuantumParams();
        }

        return instance;
    }

    double d_tau = 0.001;       // time step value
    int total_time_steps = 1e5; // total number of time steps valued d_tau
    int eq_time_step = 1e4;     // time step to average from
    int n0_walkers = 1000;      // beginning number of walkers alive, also target number of walkers
    int nmax_walkers = 1100;    // maximal number of walkers alive - size of allocated vector

    std::vector<int> save_hist_at = std::vector<int>({});
    double xmin = -5; // sampling minimum for visualisation
    double xmax = 5;  // sampling maximum for visualisation

    int n_dims = 1; // number of dimensions
    std::function<double(walker)> pot = [](walker wlk) {
        return 0.5 * (0.64 * std::pow(wlk.x(), 2));
    }; // potential of 2D quantum dot

    // number of bins used for generating wave function // bins^D is used in D dimensions
    int n_bins = 200;

    bool blocks_calibration = false;
    int n_block = pow(2, 15);

    // std::vector<double> nodes = std::vector<double>({}); // This way of applying nodes doesnt
    // work in many dimensions
    std::function<double(walker)> trial_wavef = [](walker wlk) {
        return std::pow(0.8 * M_PI, 0.25) * std::exp(-0.8 * pow(wlk.x(), 2) / 2.) ;
    };
};

#endif
