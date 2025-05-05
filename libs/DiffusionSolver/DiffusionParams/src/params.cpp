#include "DiffusionParams/include/params.hpp"
#include "TrialFunctions/include/harmonic_oscillator.hpp"
#include "include/UnitHandler.hpp"
#include <memory>

DiffusionQuantumParams *DiffusionQuantumParams::instance = nullptr;

void DiffusionQuantumParams::set_default_params() {
    d_tau = 15.;              // time step value
    total_time_steps = 1000; // total number of time steps valued d_tau
    eq_time_step = 700;      // time step to average from
    n0_walkers = 10000;       // beginning number of walkers alive, also target number of walkers
    nmax_walkers = 15000;     // maximal number of walkers alive - size of allocated vector

    save_hist_at = std::vector<int>({});                 // after equilibration phase
    xmin = UnitHandler::length(UnitHandler::TO_AU, -20); // sampling minimum for visualisation
    xmax = UnitHandler::length(UnitHandler::TO_AU, 20);  // sampling maximum for visualisation

    n_dims = 2;     // number of dimensions
    epsilon = 13.6; // relative permatibility

    std::vector<double> omegas({5., 3., 0.});
    std::transform(omegas.begin(), omegas.end(), omegas.begin(), [&](double om){ return UnitHandler::energy(UnitHandler::TO_AU, om);});

    double effective_mass = 1;

    pot_params.effective_mass = effective_mass;
    pot_params.dims = n_dims;
    pot_params.omegas = omegas;
    pot_func_builder->set_params(pot_params);

    pot = pot_func_builder->get_potential();
    n_bins = 200;

    // TODO: revise blocks
    blocks_calibration = false;
    n_block = pow(2, 15);

    HarmonicOscillatorOrbitalsParams p;
    p.dims = n_dims;
    p.effective_mass = effective_mass;
    p.omegas = omegas;
    p.excitations = std::vector<int>{1, 1, 0};
    trial_wavef = std::make_unique<HarmonicOscillatorOrbitals>(p);
}
