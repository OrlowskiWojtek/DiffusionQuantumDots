#include "Core/include/electrons.hpp"
#include "Core/include/results.hpp"
#include "Core/include/walkers.hpp"
#include <numeric>
#include <omp.h>

DiffusionQuantumElectrons::DiffusionQuantumElectrons()
    : p(DiffusionQuantumParams::getInstance())
    , walkers_helper(std::make_unique<DiffusionWalkers>())
    , results(std::make_unique<DiffusionQuantumResults>())
    , general_context(std::make_unique<SolverContext>()) {

    current_it = 0;
    accu_it = 0;

    growth_estimator = 0;
    e_block = 0;

    num_alive = p->n0_walkers;
    target_alive = p->n0_walkers;

    init_rngs();
    init_context();
    init_containers();
}

void DiffusionQuantumElectrons::init_context() {
    for (int thr_idx = 0; thr_idx < omp_get_max_threads(); thr_idx++) {
        solver_contexts.emplace_back();
    }
}

void DiffusionQuantumElectrons::init_rngs() {
    uni_rng.seed(time(0));
    uniform_generator = boost::random::uniform_real_distribution<double>(0, 1);
}

void DiffusionQuantumElectrons::init_containers() {
    boost::random::mt19937 rng;
    rng.seed(time(0));
    boost::random::uniform_real_distribution<> initial_dist(p->xmin, p->xmax);

    electrons.resize(p->nmax_walkers);
    copy_electrons.resize(p->nmax_walkers);
    for (auto &ele : electrons) {
        ele.resize(p->n_electrons);
        std::for_each(ele.begin(), ele.end(), [&](walker &wlk) {
            wlk.cords.fill(0);
        });
        std::for_each(ele.begin(), ele.end(), [&](walker &wlk) {
            std::for_each(wlk.cords.begin(), wlk.cords.begin() + p->n_dims, [&](double &cord) {
                cord = initial_dist(rng);
            });
        });
    }

    std::copy(electrons.begin(), electrons.end(), copy_electrons.begin());

    diffusion_values.resize(p->nmax_walkers);
    for (auto &diff_value : diffusion_values) {
        diff_value.resize(p->n_electrons);
        std::for_each(diff_value.begin(), diff_value.end(), [&](walker &wlk) {
            wlk.cords.fill(0);
        });
    }

    p_values.resize(p->nmax_walkers);

    summed_walkers.resize(boost::extents[p->n_bins][p->n_bins][p->n_bins]);
    std::for_each(summed_walkers.data(), summed_walkers.data() + summed_walkers.num_elements(), [](int64_t &val) {
        val = 0;
    });
}

void DiffusionQuantumElectrons::diffuse() {
    current_it++;

#pragma omp parallel for
    for (int i = 0; i < num_alive; i++) {
        int tid = omp_get_thread_num();
        solver_contexts[tid].move_walkers(electrons[i], diffusion_values[i]);
    }
}

void DiffusionQuantumElectrons::eval_p() {
#pragma omp parallel for
    for (int i = 0; i < num_alive; i++) {
        int tid = omp_get_thread_num();

        if (solver_contexts[tid].apply_nodes(electrons[i], copy_electrons[i])) {
            p_values[i] = 0;
            continue;
        }

        p_values[i] = solver_contexts[tid].p_value(electrons[i], copy_electrons[i], growth_estimator);
    }
}

void DiffusionQuantumElectrons::branch() {
    new_alive = 0;

    for (int i = 0; i < num_alive; i++) {
        int m = static_cast<int>(p_values[i] + uniform_generator(uni_rng));
        set_alive(m, electrons[i]);
    }

    num_alive = new_alive;
    std::copy(copy_electrons.begin(), copy_electrons.begin() + num_alive, electrons.begin());

    update_growth_estimator();
}

void DiffusionQuantumElectrons::set_alive(int N, const electron_walker &wlk) {
    for (int i = new_alive; i < new_alive + N; i++) {
        if (i >= static_cast<int>(electrons.size())) {
            std::cout << "Error: The programme ran out of allocated memory. "
                         "Required resize."
                      << std::endl;
            break;
        }

        copy_electrons[i] = wlk;
    }

    new_alive += N;
}

// TODO: thing about it
void DiffusionQuantumElectrons::update_growth_estimator() {
    if (accu_it == 0) {
        int nblock = current_it % p->n_block;
        e_block = (growth_estimator + e_block * nblock) / (nblock + 1);
    } else {
        e_block = (growth_estimator + e_block * accu_it) / (accu_it + 1);
    }

    growth_estimator =
        e_block - 1. / p->d_tau * std::log(static_cast<double>(num_alive) / static_cast<double>(target_alive));
}

double DiffusionQuantumElectrons::trial_wavef(const electron_walker &wlk) { return (*general_context->orbital)(wlk); }

double DiffusionQuantumElectrons::local_energy_average() {
    double loc_ene_avg = std::accumulate(electrons.begin(),
                                         electrons.begin() + num_alive,
                                         0.,
                                         [this](double acc, const electron_walker &wlk) {
                                             return acc + general_context->local_energy(wlk);
                                         }) /
                         static_cast<double>(num_alive);

    return loc_ene_avg;
}

// need to test and reinvent error estimation
void DiffusionQuantumElectrons::count() {
    binning();

    mixed_estimator = local_energy_average();

    ground_state_estimator += mixed_estimator;
    accu_it++;

    if (p->blocks_calibration) {
        results->add_energies(mixed_estimator, e_block);
    }

    // TODO add overblocks estimation
}

// TODO smart binning for visualisation is needed - like class just for binning in all dimensions x1 vs x2, y1 vs y2
// etc. now it is summing all walkers into one bin
void DiffusionQuantumElectrons::binning() {
    std::for_each(electrons.begin(), electrons.begin() + num_alive, [&](const electron_walker &wlk) {
        std::for_each(wlk.begin(), wlk.end(), [&](const walker &wlk) {
            for (int dim = 0; dim < p->n_dims; dim++) {
                if (wlk.cords[dim] < p->xmin || wlk.cords[dim] > p->xmax) {
                    return;
                }
            }

            std::array<size_t, 3> walker_bin;
            walker_bin[0] = static_cast<int>((p->xmax - wlk.cords[0]) / (p->xmax - p->xmin) * p->n_bins);
            walker_bin[1] = static_cast<int>((p->xmax - wlk.cords[1]) / (p->xmax - p->xmin) * p->n_bins);
            walker_bin[2] = static_cast<int>((p->xmax - wlk.cords[2]) / (p->xmax - p->xmin) * p->n_bins);
            summed_walkers[walker_bin[0]][walker_bin[1]][walker_bin[2]]++; // TODO move hist to result class
        });
    });
}

void DiffusionQuantumElectrons::save_progress() {
    results->add_histogram(static_cast<double>(current_it) * p->d_tau,
                           current_it,
                           ground_state_estimator / accu_it,
                           e_block,
                           summed_walkers);
}

DiffusionQuantumResults &DiffusionQuantumElectrons::get_results() { return *results; }
