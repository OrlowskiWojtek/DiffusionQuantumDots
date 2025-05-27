#include "DiffusionParams/include/params.hpp"
#include "include/UnitHandler.hpp"
#include <cassert>

DiffusionQuantumParams *DiffusionQuantumParams::instance = nullptr;

void DiffusionQuantumParams::set_default_params() {
    d_tau = 10;              // time step value
    total_time_steps = 300; // total number of time steps valued d_tau
    eq_time_step = 200;     // time step to average from
    n0_walkers = 10000;      // beginning number of walkers alive, also target number of walkers
    nmax_walkers = 15000;    // maximal number of walkers alive - size of allocated vector

    save_hist_at = std::vector<int>({});                 // after equilibration phase
    xmin = UnitHandler::length(UnitHandler::TO_AU, -10); // sampling minimum for visualisation
    xmax = UnitHandler::length(UnitHandler::TO_AU, 10);  // sampling maximum for visualisation

    n_electrons = 2;
    n_dims = 1;     // number of dimensions
    epsilon = 13.6; // relative permatibility

    effective_mass = 0.067;
    omegas = {5. / std::sqrt(effective_mass), 3. / std::sqrt(effective_mass), 0.};
    std::transform(omegas.begin(), omegas.end(), omegas.begin(), [&](double om) {
        return UnitHandler::energy(UnitHandler::TO_AU, om);
    });

    pot_params.effective_mass = effective_mass;
    pot_params.dims = n_dims;
    pot_params.omegas = omegas;

    n_bins = 200;

    vis_dim_idx_x = std::make_tuple(0, 0);
    vis_dim_idx_y = std::make_tuple(1, 0);

    // TODO: revise blocks
    blocks_calibration = true;
    n_block = pow(2, 15);

    check_params();
}

void DiffusionQuantumParams::check_params() {
    assert(std::get<0>(vis_dim_idx_x) < n_electrons);
    assert(std::get<1>(vis_dim_idx_x) < n_dims);

    assert(std::get<0>(vis_dim_idx_y) < n_electrons);
    assert(std::get<1>(vis_dim_idx_y) < n_dims);
}
