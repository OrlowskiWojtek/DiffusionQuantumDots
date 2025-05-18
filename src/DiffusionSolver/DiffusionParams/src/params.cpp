#include "DiffusionParams/include/params.hpp"
#include "TrialFunctions/include/harmonic_oscillator.hpp"
#include "TrialFunctions/include/harmonic_oscillator.hpp"
#include "TrialFunctions/include/jastrow_slater.hpp"
#include "include/UnitHandler.hpp"
#include <memory>

DiffusionQuantumParams *DiffusionQuantumParams::instance = nullptr;

void DiffusionQuantumParams::set_default_params() {
    d_tau = 10.;              // time step value
    total_time_steps = 30000; // total number of time steps valued d_tau
    eq_time_step = 20000;      // time step to average from
    n0_walkers = 1000;       // beginning number of walkers alive, also target number of walkers
    nmax_walkers = 1500;     // maximal number of walkers alive - size of allocated vector

    save_hist_at = std::vector<int>({});                 // after equilibration phase
    xmin = UnitHandler::length(UnitHandler::TO_AU, -10); // sampling minimum for visualisation
    xmax = UnitHandler::length(UnitHandler::TO_AU, 10);  // sampling maximum for visualisation

    n_electrons = 2;
    n_dims = 2;     // number of dimensions
    epsilon = 13.6; // relative permatibility

    double effective_mass = 0.067;
    std::vector<double> omegas({5. / std::sqrt(effective_mass), 3. / std::sqrt(effective_mass), 0.});
    std::transform(omegas.begin(), omegas.end(), omegas.begin(), [&](double om){ return UnitHandler::energy(UnitHandler::TO_AU, om);});

    pot_params.effective_mass = effective_mass;
    pot_params.dims = n_dims;
    pot_params.omegas = omegas;
    pot_func_builder->set_params(pot_params);

    pot = pot_func_builder->get_potential();
    n_bins = 200;

    // TODO: revise blocks
    blocks_calibration = false;
    n_block = pow(2, 15);

    JastrowSlaterOrbitalParams p;
    p.electron_number = n_electrons;
    p.omegas = omegas;
    p.effective_mass = effective_mass;
    p.dims = n_dims;
    p.a = 0.005;
    p.b = 2.;

    trial_wavef = std::make_unique<JastrowSlaterOrbital>(p);
}
