#ifndef DIFFUSION_QUANTUM_PARAMS_HPP
#define DIFFUSION_QUANTUM_PARAMS_HPP

#include "DiffusionParams/include/harmonic_potential.hpp"

#include <cmath>
#include <vector>

// Class implements Singleton design pattern
class DiffusionQuantumParams {
private:
    static DiffusionQuantumParams *instance;
    DiffusionQuantumParams() {
        set_default_params();
    }

    void check_params();

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
    double d_tau; // time step value

    int total_time_steps;   // total number of time steps valued d_tau
    int eq_time_step;       // time step to average from
    int initial_time_steps; // initial time steps to initilize walkers with |psi_trial| distribution
    int vmc_sampling_time_steps; // time steps for sampling vmc local energy
    int n0_walkers;              // beginning number of walkers alive, also target number of walkers
    int nmax_walkers;            // maximal number of walkers alive - size of allocated vector

    std::vector<int> save_hist_at = std::vector<int>({}); // after equilibration phase
    double xmin;                                          // sampling minimum for visualisation
    double xmax;                                          // sampling maximum for visualisation
    bool show_visualisation;                              // shows plots after simulation

    int n_electrons; // number of electrons;
    int n_dims;      // number of dimensions

    HarmonicPotentialParams pot_params;

    double effective_mass;           // effective mass of electron
    double epsilon;                  // relative permatibility
    std::vector<double> omegas;      // omega in each direction
    std::vector<ElectronSpin> spins; // spin for each electron
    int n_bins;

    std::tuple<int, int> vis_dim_idx_x; // used for visualisation - first dimension (first index is
                                        // electron idx, second is dimension idx)
    std::tuple<int, int> vis_dim_idx_y; // used for visualisation - second dimension (first index is
                                        // electron idx, second is dimension idx)

    int total_vis_idx_x; // used for visualisation of total walkers - sum over all electrons
    int total_vis_idx_y; // used for visualisation of total walkers - sum over all electrons

    // save energies to file every after equilibrium phase
    int save_every;

    // TODO: revise blocks
    bool blocks_calibration;
    int n_block;

    // TODO: remove, quick scan ofer a b parameters -> handle in minimum seeker
    double a = 0.005;
    double b = 2.;
};

#endif
