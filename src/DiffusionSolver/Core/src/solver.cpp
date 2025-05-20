#include "Core/include/solver.hpp"

DiffusionQuantumSolver::DiffusionQuantumSolver()
    : electrons(std::make_unique<DiffusionQuantumElectrons>())
    , block_analyzer(std::make_unique<EnergyBlockingAnalyzer>())
    , params(DiffusionQuantumParams::getInstance()) {}

DiffusionQuantumSolver::~DiffusionQuantumSolver() {}

// TODO: init also filenames and other solve related systems
void DiffusionQuantumSolver::init() {
    std::sort(params->save_hist_at.begin(), params->save_hist_at.end());
    params->save_hist_at.push_back(params->total_time_steps - params->eq_time_step - 1); // always save last one
    save_counter = 0;
}

void DiffusionQuantumSolver::solve() {
    init();

    int equilibrate_loop = params->eq_time_step;
    int collect_loop = params->total_time_steps - params->eq_time_step;

    for (int i = 0; i < equilibrate_loop; i++) {
        diffuse();
        branch();
    }


    for (int i = 0; i < collect_loop; i++) {
        diffuse();
        branch();
        accumulate();

        check_saving(i);
    }


    final_results = electrons->get_results();
    final_results.save_to_file();

    if (params->blocks_calibration) {
        std::vector<double> energies = final_results.get_energies();
        block_analyzer->blocking_analysis(energies);
    }
}

void DiffusionQuantumSolver::diffuse() {
    electrons->diffuse(); 
    electrons->accept_movement();
}

void DiffusionQuantumSolver::branch() {
    electrons->eval_p();
    electrons->branch();
}

void DiffusionQuantumSolver::accumulate() { electrons->count(); }

void DiffusionQuantumSolver::check_saving(int iter_idx) {
    if (iter_idx == params->save_hist_at[save_counter]) {
        electrons->save_progress();
        if (save_counter < params->save_hist_at.size() - 1) {
            save_counter++;
        }
    }
}
