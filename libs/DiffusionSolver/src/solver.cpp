#include "solver.hpp"

DiffusionQuantumSolver::DiffusionQuantumSolver()
    : walkers(std::make_unique<DiffusionWalkers>())
    , block_analyzer(std::make_unique<EnergyBlockingAnalyzer>()) {}

DiffusionQuantumSolver::~DiffusionQuantumSolver() {}

// TODO: init also filenames and other solve related systems
void DiffusionQuantumSolver::init() {
    walkers->init_walkers(params);

    std::sort(params.save_hist_at.begin(), params.save_hist_at.end());
    save_counter = 0;
}

void DiffusionQuantumSolver::solve() {
    init();

    for (int i = 0; i < params.total_time_steps; i++) {
        diffuse();
        branch();
        if (i > params.total_time_steps / 5) {
            accumulate();
        }

        if (i == params.save_hist_at[save_counter]) {
            walkers->save_progress();
            if (save_counter < params.save_hist_at.size() - 1) {
                save_counter++;
            }
        }
    }
    final_results = walkers->get_results();
    final_results.save_to_file();

    if (params.blocks_calibration) {
        std::vector<double> energies = final_results.get_energies();
        block_analyzer->blocking_analysis(energies);
    }
}

void DiffusionQuantumSolver::diffuse() { walkers->diffuse(); }

void DiffusionQuantumSolver::branch() {
    walkers->eval_p();
    walkers->branch();
}

void DiffusionQuantumSolver::accumulate() { walkers->count(); }
