#include "Core/include/solver.hpp"
#include <armadillo>

DiffusionQuantumSolver::DiffusionQuantumSolver()
    : electrons(std::make_unique<DiffusionQuantumElectrons>())
    , block_analyzer(std::make_unique<EnergyBlockingAnalyzer>())
    , vis(std::make_unique<WalkersVisualiser>())
    , params(DiffusionQuantumParams::getInstance()) {
}

DiffusionQuantumSolver::~DiffusionQuantumSolver() {
}

// TODO: init also filenames and other solve related systems
void DiffusionQuantumSolver::init() {}

void DiffusionQuantumSolver::solve() {
    init();

    int equilibrate_loop = params->eq_time_step;
    int collect_loop = params->total_time_steps - params->eq_time_step;

#ifndef PURE_DIFFUSION
    for (int i = 0; i < params->initial_time_steps; i++) {
        initialize_distribution();
    }
    for (int i = 0; i < params->vmc_sampling_time_steps; i++) {
        sample_vmc_energy();
    }
    finish_initialization();
    std::cout << "===Finished initialization loop===" << std::endl;
#endif

    return;

    for (int i = 0; i < equilibrate_loop; i++) {
        diffuse();
        branch();
    }

    std::cout << "===Finished equilibration loop===" << std::endl;

    for (int i = 0; i < collect_loop; i++) {
        diffuse();
        branch();
        accumulate();
    }

    std::cout << "===Finished collection loop===" << std::endl;

    electrons->save_progress();
    electrons->get_results().save_to_file();

    if (params->show_visualisation) {
        vis->make_surf_plot(electrons->get_results().get_last_psi(),
                            electrons->get_results().get_last_total_psi(),
                            params->n_bins);
    }

    if (params->blocks_calibration) {
        block_analyzer->blocking_analysis(electrons->get_results().get_calib_mixed_energies(),
                                          electrons->get_results().get_calib_growth_energies());
    }
}

void DiffusionQuantumSolver::diffuse() {
    electrons->diffuse();
#ifndef PURE_DIFFUSION
    electrons->check_movement();
#endif
}

void DiffusionQuantumSolver::branch() {
    electrons->prepare_branch();
    electrons->branch();
}

void DiffusionQuantumSolver::accumulate() {
    electrons->count();
}

void DiffusionQuantumSolver::initialize_distribution() {
    electrons->initial_diffusion();
}

void DiffusionQuantumSolver::sample_vmc_energy() {
    electrons->sample_variational_energy();
}

void DiffusionQuantumSolver::finish_initialization() {
    electrons->finish_initial_diffusion();
}
