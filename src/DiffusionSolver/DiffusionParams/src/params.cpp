#include "DiffusionParams/include/params.hpp"
#include "include/UnitHandler.hpp"
#include <cassert>

DiffusionQuantumParams *DiffusionQuantumParams::instance = nullptr;

void DiffusionQuantumParams::set_default_params() {
    d_tau = 1.;              // time step value
    total_time_steps = 5000; // total number of time steps valued d_tau
    eq_time_step = 2000;     // time step to average from
    n0_walkers = 10000;      // beginning number of walkers alive, also target number of walkers
    nmax_walkers = 200000;    // maximal number of walkers alive - size of allocated vector

    save_hist_at = std::vector<int>({});                 // after equilibration phase
    xmin = UnitHandler::length(UnitHandler::TO_AU, -50); // sampling minimum for visualisation
    xmax = UnitHandler::length(UnitHandler::TO_AU, 50);  // sampling maximum for visualisation

    n_electrons = 1;
    n_dims = 2;    // number of dimensions
    epsilon = 12.; // relative permatibility

    effective_mass = 0.067;
    omegas = {3., 5., 0.};

    std::transform(omegas.begin(), omegas.end(), omegas.begin(), [&](double om) {
        return UnitHandler::energy(UnitHandler::TO_AU, om);
    });

    pot_params.effective_mass = effective_mass;
    pot_params.dims = n_dims;
    pot_params.omegas = omegas;

    n_bins = 100;

    vis_dim_idx_x = std::make_tuple(0, 0);
    vis_dim_idx_y = std::make_tuple(0, 1);

    // TODO: revise blocks
    blocks_calibration = true;
    n_block = pow(2, 15);

    show_visualisation = true;

    check_params();
}

void DiffusionQuantumParams::check_params() {
    assert(std::get<0>(vis_dim_idx_x) < n_electrons);
    assert(std::get<1>(vis_dim_idx_x) < n_dims);

    assert(std::get<0>(vis_dim_idx_y) < n_electrons);
    assert(std::get<1>(vis_dim_idx_y) < n_dims);
}
