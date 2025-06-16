#include "DiffusionParams/include/params.hpp"
#include "include/UnitHandler.hpp"
#include <cassert>

DiffusionQuantumParams *DiffusionQuantumParams::instance = nullptr;

void DiffusionQuantumParams::set_default_params() {
    d_tau = 0.1;             // time step
    
    initial_time_steps = 200000;
    vmc_sampling_time_steps = 2000;

    eq_time_step = 20000;     
    total_time_steps = 80000;  

    n0_walkers = 10000;      // beginning number of walkers alive, also target number of walkers
    nmax_walkers = 20000;    // maximal number of walkers alive - size of allocated vector

    save_hist_at = std::vector<int>({});                 // after equilibration phase
    xmin = UnitHandler::length(UnitHandler::TO_AU, -50); // sampling minimum for visualisation
    xmax = UnitHandler::length(UnitHandler::TO_AU, 50);  // sampling maximum for visualisation

    n_electrons = 2;
    n_dims = 2;    // number of dimensions
    epsilon = 12.; // relative permatibility

    effective_mass = 0.067;
    omegas = {3., 5., 0.};

    max_drift_length = 50. / (std::sqrt(d_tau / effective_mass));

    std::transform(omegas.begin(), omegas.end(), omegas.begin(), [&](double om) {
        return UnitHandler::energy(UnitHandler::TO_AU, om);
    });

    pot_params.effective_mass = effective_mass;
    pot_params.dims = n_dims;
    pot_params.omegas = omegas;

    n_bins = 100;

    vis_dim_idx_x = std::make_tuple(0, 0);
    vis_dim_idx_y = std::make_tuple(1, 0);

    total_vis_idx_x = 0;
    total_vis_idx_y = 1;

    save_every = 50;

    // TODO: revise blocks
    blocks_calibration = true;
    n_block = pow(2, 15);

    show_visualisation = true;

    a = 0.25;
    b = 0.1575;

    check_params();
}

void DiffusionQuantumParams::check_params() {
    assert(std::get<0>(vis_dim_idx_x) < n_electrons);
    assert(std::get<1>(vis_dim_idx_x) < n_dims);

    assert(std::get<0>(vis_dim_idx_y) < n_electrons);
    assert(std::get<1>(vis_dim_idx_y) < n_dims);

    assert((total_vis_idx_x) < n_dims);
    assert((total_vis_idx_y) < n_dims);

    assert(save_every > 0);
}
