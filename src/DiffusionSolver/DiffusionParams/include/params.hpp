#ifndef DIFFUSION_QUANTUM_PARAMS_HPP
#define DIFFUSION_QUANTUM_PARAMS_HPP

#include "DiffusionParams/include/harmonic_potential.hpp"
#include "TrialFunctions/include/abstract_singlebody_orbital.hpp"

#include "Core/include/walkers_struct.hpp"
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

// Class implements Singleton design pattern
// TODO: switch so it cointains ONLY params (no trial_wavef, only params to build wavef)
class DiffusionQuantumParams {
private:
    std::unique_ptr<HarmonicPotentialBuilder> pot_func_builder;

    static DiffusionQuantumParams *instance;
    DiffusionQuantumParams()
        : pot_func_builder(std::make_unique<HarmonicPotentialBuilder>()) {
        set_default_params();
    }

public:
    DiffusionQuantumParams(const DiffusionQuantumParams &) = delete;
    DiffusionQuantumParams &operator=(const DiffusionQuantumParams &) = delete;

    static DiffusionQuantumParams *getInstance() {
        if (!instance) {
            instance = new DiffusionQuantumParams();
        }

        return instance;
    }

    void set_default_params();
    double d_tau;         // time step value
    int total_time_steps; // total number of time steps valued d_tau
    int eq_time_step;      // time step to average from
    int n0_walkers;       // beginning number of walkers alive, also target number of walkers
    int nmax_walkers;     // maximal number of walkers alive - size of allocated vector

    std::vector<int> save_hist_at = std::vector<int>({}); // after equilibration phase
    double xmin; // sampling minimum for visualisation
    double xmax; // sampling maximum for visualisation

    int n_electrons; // number of electrons;
    int n_dims;      // number of dimensions

    HarmonicPotentialParams pot_params;

    double epsilon;   // relative permatibility
    std::vector<double> omegas; // omega in each direction

    std::function<double(const walker&)> pot;
    int n_bins;

    // TODO: revise blocks
    bool blocks_calibration;
    int n_block;

    std::unique_ptr<AbstractSinglebodyOrbital> trial_wavef;
};

#endif
