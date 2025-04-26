#include "DiffusionParams/include/params.hpp"
#include "include/UnitHandler.hpp"

DiffusionQuantumParams* DiffusionQuantumParams::instance = nullptr;

void DiffusionQuantumParams::set_default_params(){
    d_tau = 0.001;         // time step value
    total_time_steps = 10000; // total number of time steps valued d_tau
    eq_time_step = 8000;      // time step to average from
    n0_walkers = 10000;       // beginning number of walkers alive, also target number of walkers
    nmax_walkers = 11000;     // maximal number of walkers alive - size of allocated vector

    save_hist_at = std::vector<int>({}); // after equilibration phase
    xmin = -5;//UnitHandler::energy(UnitHandler::TO_AU, -5); // sampling minimum for visualisation
    xmax = 5;//UnitHandler::energy(UnitHandler::TO_AU, 5);  // sampling maximum for visualisation

    n_dims = 2; // number of dimensions
    epsilon = 13.6;   // relative permatibility

    pot_params.effective_mass = 1;//0.067;
    pot_params.dims = n_dims;
    pot_params.omegas = {0.8, 0.6, 0.};
    pot_func_builder->set_params(pot_params);

    pot = pot_func_builder->get_potential();
    n_bins = 200;

    // TODO: revise blocks
    blocks_calibration = false;
    n_block = pow(2, 15);

    trial_wavef = [](walker wlk) {
        wlk.x(); return 1; //std::pow(0.8 * M_PI, 0.25) * std::exp(-0.8 * pow(wlk.x(), 2) / 2.);
    };
}
