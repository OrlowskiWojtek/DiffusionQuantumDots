#include "Core/include/electrons.hpp"
#include "Core/include/results.hpp"
#include "Core/include/walkers.hpp"
#include "include/UnitHandler.hpp"
#include <numeric>
#include <omp.h>

DiffusionQuantumElectrons::DiffusionQuantumElectrons()
    : p(DiffusionQuantumParams::getInstance())
    , results(std::make_unique<DiffusionQuantumResults>())
    , general_context(std::make_unique<SolverContext>()) {

    current_it = 0;

    stats.reset();
    e_block = 0;

    num_alive = p->n0_walkers;
    target_alive = p->n0_walkers;

    eff_d_tau = p->d_tau;

    init_rngs();
    init_context();
    init_containers();
}

void DiffusionQuantumElectrons::init_context() {
    solver_contexts.clear();

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
        electron_walker &ele_wlk = ele.get_walker();
        ele_wlk.resize(p->n_electrons);
    
        std::for_each(ele_wlk.begin(), ele_wlk.end(), [&](walker &wlk) {
            wlk.cords.fill(0);
        });

        std::for_each(ele_wlk.begin(), ele_wlk.end(), [&](walker &wlk) {
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

    summed_walkers.resize(boost::extents[p->n_bins][p->n_bins]);
    std::for_each(summed_walkers.data(),
                  summed_walkers.data() + summed_walkers.num_elements(),
                  [](int64_t &val) {
                      val = 0;
                  });

    total_summed_walkers.resize(boost::extents[p->n_bins][p->n_bins]);
    std::for_each(total_summed_walkers.data(),
                  total_summed_walkers.data() + total_summed_walkers.num_elements(),
                  [](int64_t &val) {
                      val = 0;
                  });

    std::for_each(electrons.begin(), electrons.end(), [&](ElectronWalker &wlk) {
        general_context->calc_trial_wavef(wlk);
        general_context->calc_local_energy(wlk);
    });

    std::for_each(copy_electrons.begin(), copy_electrons.end(), [&](ElectronWalker &wlk) {
        general_context->calc_trial_wavef(wlk);
        general_context->calc_local_energy(wlk);
    });
}

void DiffusionQuantumElectrons::diffuse() {
    current_it++;

#pragma omp parallel for
    for (int i = 0; i < num_alive; i++) {
        int tid = omp_get_thread_num();
        solver_contexts[tid].move_walkers(electrons[i], diffusion_values[i]);
        solver_contexts[tid].calc_trial_wavef(electrons[i]);
    }
}

void DiffusionQuantumElectrons::check_movement() {
    int metropolis_accepted = 0;

#pragma omp parallel for reduction(+ : metropolis_accepted)
    for (int i = 0; i < num_alive; i++) {
        int tid = omp_get_thread_num();

        if (solver_contexts[tid].check_nodes(electrons[i], copy_electrons[i])) {
            solver_contexts[tid].reject_move(electrons[i], copy_electrons[i]);
            continue;
        }

        if (solver_contexts[tid].check_metropolis(
                electrons[i], copy_electrons[i], diffusion_values[i])) {
            solver_contexts[tid].calc_local_energy(electrons[i]);
            metropolis_accepted++;
        } else {
            solver_contexts[tid].reject_move(electrons[i], copy_electrons[i]);
        }
    }

    eff_d_tau =
        p->d_tau * static_cast<double>(metropolis_accepted) / static_cast<double>(num_alive);
}

void DiffusionQuantumElectrons::prepare_branch() {
#pragma omp parallel for
    for (int i = 0; i < num_alive; i++) {
        int tid = omp_get_thread_num();

#ifdef PURE_DIFFUSION
        if (solver_contexts[tid].check_nodes(electrons[i], copy_electrons[i])) {
            p_values[i] = 0;
            continue;
        }
#endif

        p_values[i] = solver_contexts[tid].p_value(
            electrons[i], copy_electrons[i], stats.growth_estimator, eff_d_tau);
    }
}

void DiffusionQuantumElectrons::branch() {
    new_alive = 0;

    for (int i = 0; i < num_alive; i++) {
        double random_number = uniform_generator(uni_rng);
        int m = static_cast<int>(p_values[i] + random_number);
        if (std::isnan(m) || m < 0) {
            m = 0;
        }
        if (m > 4) {
            m = 4;
        }
        set_alive(m, electrons[i]);
    }

    num_alive = new_alive;

    if (num_alive >= p->nmax_walkers) {
        num_alive = p->nmax_walkers / 2;
    }
    std::copy(copy_electrons.begin(), copy_electrons.begin() + num_alive, electrons.begin());

    update_growth_estimator();
}

void DiffusionQuantumElectrons::set_alive(int N, const ElectronWalker &wlk) {
    for (int i = new_alive; i < new_alive + N; i++) {
        if (i >= static_cast<int>(electrons.size())) {
            break;
        }

        copy_electrons[i] = wlk;
    }

    new_alive += N;
}

void DiffusionQuantumElectrons::update_growth_estimator() {
    if (stats.it == 0) {
        int nblock = current_it % p->n_block;
        e_block = (stats.growth_estimator + e_block * nblock) / (nblock + 1.);
    } else {
        e_block = stats.acc_growth_estimator / static_cast<double>(stats.it);
    }

    stats.growth_estimator =
        e_block - 1. / (p->d_tau) *
                      std::log(static_cast<double>(num_alive) / static_cast<double>(target_alive));
}

double DiffusionQuantumElectrons::local_energy_average() {
#ifndef PURE_DIFFUSION
    double loc_ene_avg = std::accumulate(electrons.begin(),
                                         electrons.begin() + num_alive,
                                         0.,
                                         [](double acc, const ElectronWalker &ele) {
                                             return acc + ele.local_energy;
                                         });
#else
    double loc_ene_avg = std::accumulate(electrons.begin(),
                                         electrons.begin() + num_alive,
                                         0.,
                                         [&](double acc, const ElectronWalker &ele) {
                                             return acc + general_context->get_potential(ele);
                                         });
#endif

    return loc_ene_avg / static_cast<double>(num_alive);
}

// need to test and reinvent error estimation
void DiffusionQuantumElectrons::count() {
    binning();
    total_binning();

    stats.mixed_estimator = local_energy_average();

    stats.acc_mixed_estimator += stats.mixed_estimator;

    stats.acc_growth_estimator += stats.growth_estimator;
    stats.acc_sq_growth_estimator += std::pow(stats.growth_estimator, 2);

    stats.it++;

    if (stats.it % p->save_every == 0) {
        results->save_energies(static_cast<double>(stats.it) * p->d_tau, num_alive, this->stats);
    }

    if (p->blocks_calibration) {
        results->add_energies(stats.mixed_estimator, stats.growth_estimator);
    }
}

void DiffusionQuantumElectrons::binning() {
    std::for_each(electrons.begin(), electrons.begin() + num_alive, [&](const ElectronWalker &wlk) {
        int x_ele = std::get<0>(p->vis_dim_idx_x);
        int x_dim = std::get<1>(p->vis_dim_idx_x);
        int y_ele = std::get<0>(p->vis_dim_idx_y);
        int y_dim = std::get<1>(p->vis_dim_idx_y);

        if (wlk.get_const_walker()[x_ele].cords[x_dim] < p->xmin ||
            wlk.get_const_walker()[x_ele].cords[x_dim] > p->xmax) {
            return;
        }

        if (wlk.get_const_walker()[y_ele].cords[y_dim] < p->xmin ||
            wlk.get_const_walker()[y_ele].cords[y_dim] > p->xmax) {
            return;
        }

        std::array<size_t, 2> walker_bin;

        walker_bin[0] = static_cast<int>((p->xmax - wlk.get_const_walker()[x_ele].cords[x_dim]) /
                                         (p->xmax - p->xmin) * p->n_bins);
        walker_bin[1] = static_cast<int>((p->xmax - wlk.get_const_walker()[y_ele].cords[y_dim]) /
                                         (p->xmax - p->xmin) * p->n_bins);

        if (p->n_dims == 1 && p->n_electrons == 1) {
            walker_bin[1] = 0;
        }

        summed_walkers[walker_bin[0]][walker_bin[1]]++;
    });
}

void DiffusionQuantumElectrons::total_binning() {
    std::for_each(electrons.begin(), electrons.begin() + num_alive, [&](const ElectronWalker &wlk) {
        std::for_each(
            wlk.get_const_walker().begin(), wlk.get_const_walker().end(), [&](const walker &wlk) {
                int x_dim = p->total_vis_idx_x;
                int y_dim = p->total_vis_idx_y;

                if (wlk.cords[x_dim] < p->xmin || wlk.cords[x_dim] > p->xmax) {
                    return;
                }

                if (wlk.cords[y_dim] < p->xmin || wlk.cords[y_dim] > p->xmax) {
                    return;
                }

                std::array<size_t, 2> walker_bin;

                walker_bin[0] = static_cast<int>((p->xmax - wlk.cords[x_dim]) /
                                                 (p->xmax - p->xmin) * p->n_bins);
                walker_bin[1] = static_cast<int>((p->xmax - wlk.cords[y_dim]) /
                                                 (p->xmax - p->xmin) * p->n_bins);

                if (p->n_dims == 1) {
                    walker_bin[1] = 0;
                }

                total_summed_walkers[walker_bin[0]][walker_bin[1]]++;
            });
    });
}

void DiffusionQuantumElectrons::save_progress() {
    results->add_histogram(static_cast<double>(current_it) * p->d_tau,
                           current_it,
                           stats,
                           summed_walkers,
                           total_summed_walkers);
}

DiffusionQuantumResults &DiffusionQuantumElectrons::get_results() {
    return *results;
}

void DiffusionQuantumElectrons::initial_diffusion() {
#pragma omp parallel for
    for (int i = 0; i < num_alive; i++) {
        int tid = omp_get_thread_num();
        copy_electrons[i] = electrons[i];

        solver_contexts[tid].apply_diffusion(electrons[i]);
        solver_contexts[tid].calc_trial_wavef(electrons[i]);
        if (!solver_contexts[tid].check_initial_metropolis(electrons[i], copy_electrons[i])) {
            solver_contexts[tid].reject_move(electrons[i], copy_electrons[i]);
        }
    }
}

void DiffusionQuantumElectrons::sample_variational_energy() {
#pragma omp parallel for
    for (int i = 0; i < num_alive; i++) {
        int tid = omp_get_thread_num();
        copy_electrons[i] = electrons[i];

        solver_contexts[tid].apply_diffusion(electrons[i]);
        solver_contexts[tid].calc_trial_wavef(electrons[i]);
        if (!solver_contexts[tid].check_initial_metropolis(electrons[i], copy_electrons[i])) {
            solver_contexts[tid].reject_move(electrons[i], copy_electrons[i]);
        } else {
            electrons[i].local_energy = solver_contexts[tid].local_energy(electrons[i]);
        }
    }

    acc_variational_energy += local_energy_average();
}

void DiffusionQuantumElectrons::finish_initial_diffusion() {
    variational_energy = acc_variational_energy / static_cast<double>(p->vmc_sampling_time_steps);

    std::cout << "===Received variational energy: "
              << UnitHandler::energy(UnitHandler::TO_DEFAULT, variational_energy)
              << "===" << std::endl;

    general_context->set_local_energy_cutoff(variational_energy);
    for (auto &cont : solver_contexts) {
        cont.set_local_energy_cutoff(variational_energy);
    }
}
