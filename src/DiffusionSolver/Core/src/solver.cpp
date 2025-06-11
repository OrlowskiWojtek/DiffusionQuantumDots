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
void DiffusionQuantumSolver::init() {
    std::sort(params->save_hist_at.begin(), params->save_hist_at.end());
    params->save_hist_at.push_back(params->total_time_steps - params->eq_time_step -
                                   1); // always save last one
    save_counter = 0;
}

void DiffusionQuantumSolver::solve() {
    init();

    int equilibrate_loop = params->eq_time_step;
    int collect_loop = params->total_time_steps - params->eq_time_step;

    for (int i = 0; i < 10000; i++) {
        initialize_distribution();
    }

    for (int i = 0; i < equilibrate_loop; i++) {
        diffuse();
        branch();
    }

    params->d_tau = params->equi_d_tau;
    for (int i = 0; i < collect_loop; i++) {
        diffuse();
        branch();
        accumulate();

        check_saving(i);
    }
 
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
}

void DiffusionQuantumSolver::branch() {
    electrons->prepare_branch();
    electrons->branch();
}

void DiffusionQuantumSolver::accumulate() {
    electrons->count();
}

void DiffusionQuantumSolver::check_saving(int iter_idx) {
    if (iter_idx == params->save_hist_at[save_counter]) {
        electrons->save_progress();
        if (save_counter < params->save_hist_at.size() - 1) {
            save_counter++;
        }
    }
}

void DiffusionQuantumSolver::initialize_distribution() {
    electrons->initial_diffusion();
}
