#include "Core/include/electrons.hpp"
#include "Core/include/results.hpp"
#include "Core/include/walkers.hpp"
#include <execution>
#include <numeric>

DiffusionQuantumElectrons::DiffusionQuantumElectrons()
    : p(DiffusionQuantumParams::getInstance())
    , walkers_helper(std::make_unique<DiffusionWalkers>())
    , results(std::make_unique<DiffusionQuantumResults>()) {

    current_it = 0;
    accu_it = 0;

    growth_estimator = 0;
    e_block = 0;

    num_alive = p->n0_walkers;
    target_alive = p->n0_walkers;

    init_rngs();
    init_containers();
}

void DiffusionQuantumElectrons::init_rngs() {
    movement_generator = boost::random::normal_distribution<double>(0, std::sqrt(p->d_tau));
    uniform_generator = boost::random::uniform_real_distribution<double>(0, 1);
}

void DiffusionQuantumElectrons::init_containers() {
    boost::random::mt19937 rng;
    rng.seed(std::time(0));
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

    drift_velocity.resize(p->n_electrons);
    m_front_walker_buffer.resize(p->n_electrons);
    m_back_walker_buffer.resize(p->n_electrons);
}

void DiffusionQuantumElectrons::diffuse() {
    current_it++;

    std::for_each(diffusion_values.begin(), diffusion_values.begin() + num_alive, [&](electron_walker &diff_value) {
        prepare_diffusion(diff_value);
    });

    std::for_each(electrons.begin(), electrons.begin() + num_alive, [&](electron_walker &ele) {
        prepare_drift(ele);
        apply_drift(ele);
    });

    for (int i = 0; i < num_alive; i++) {
        apply_diffusion(electrons[i], diffusion_values[i]);
    }
}

void DiffusionQuantumElectrons::eval_p() {
    std::transform(std::execution::seq, // smarter threading is needed
                   electrons.begin(),
                   electrons.begin() + num_alive,
                   copy_electrons.begin(),
                   p_values.begin(),
                   [&](electron_walker &wlk, electron_walker &prev_wlk) {
                       if (apply_nodes(wlk, prev_wlk)) {
                           return 0.;
                       }

                       return p_value(wlk, prev_wlk);
                   });
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

double DiffusionQuantumElectrons::trial_wavef(const electron_walker &wlk) { return (*p->trial_wavef)(wlk); }

bool DiffusionQuantumElectrons::apply_nodes(const electron_walker &wlk, const electron_walker &prev_wlk) {
    return (trial_wavef(wlk) * trial_wavef(prev_wlk)) <= 0;
}

double DiffusionQuantumElectrons::p_value(const electron_walker &wlk, const electron_walker &prev_wlk) {
    return std::exp(-p->d_tau * ((local_energy(wlk) + local_energy(prev_wlk)) / 2. - growth_estimator));
}

double DiffusionQuantumElectrons::local_energy(const electron_walker &wlk) {
    double dr = 1e-6;
    double dr2 = 1e-12;
    double kinetic_term = 0;

    double cent_value = trial_wavef(wlk);

    m_back_walker_buffer = wlk;
    m_front_walker_buffer = wlk;
    for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
        for (int d = 0; d < p->n_dims; d++) {

            m_back_walker_buffer[wlk_idx].cords[d] -= dr;
            m_front_walker_buffer[wlk_idx].cords[d] += dr;
            double back_value = trial_wavef(m_back_walker_buffer);
            double forw_value = trial_wavef(m_front_walker_buffer);

            kinetic_term += (back_value - 2 * cent_value + forw_value);

            m_back_walker_buffer[wlk_idx].cords[d] = wlk[wlk_idx].cords[d];
            m_front_walker_buffer[wlk_idx].cords[d] = wlk[wlk_idx].cords[d];
        }
    }

    kinetic_term = -0.5 * kinetic_term / (cent_value * dr2);

    double potential_term = std::accumulate(wlk.begin(), wlk.end(), 0., [this](double acc, const walker &single_wlk) {
        return acc + p->pot(single_wlk);
    });

    // TODO: hardcoded for now for coulomb interaction -> possibly change
    double interaction = 0;
    for (int i = 0; i < p->n_electrons; i++) {
        for (int j = i + 1; j < p->n_electrons; j++) {
            interaction += 1 / (p->epsilon * walkers_helper->distance(wlk[i], wlk[j], p->n_dims));
        }
    }

    return kinetic_term + interaction + potential_term;
}

double DiffusionQuantumElectrons::local_energy_average() {
    double loc_ene_avg = std::accumulate(electrons.begin(),
                                         electrons.begin() + num_alive,
                                         0.,
                                         [this](double acc, const electron_walker &wlk) {
                                             return acc + local_energy(wlk);
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

void DiffusionQuantumElectrons::prepare_drift(const electron_walker &ele_wlk) {
    double dr = 1e-6;
    double cent_value = trial_wavef(ele_wlk);

    m_front_walker_buffer = ele_wlk;

    for (int d = 0; d < p->n_dims; d++) {
        for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
            m_front_walker_buffer[wlk_idx].cords[d] += dr;

            double forw_value = trial_wavef(m_front_walker_buffer);
            drift_velocity[wlk_idx].cords[d] = (forw_value - cent_value) / (dr * cent_value);
            m_front_walker_buffer[wlk_idx].cords[d] = ele_wlk[wlk_idx].cords[d];
        }
    }
}

void DiffusionQuantumElectrons::apply_drift(electron_walker &ele_wlk) {
    for (int d = 0; d < p->n_dims; d++) {
        for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
            ele_wlk[wlk_idx].cords[d] += p->d_tau * drift_velocity[wlk_idx].cords[d];
        }
    }
}

void DiffusionQuantumElectrons::prepare_diffusion(electron_walker &diff_wlk) {
    std::for_each(diff_wlk.begin(), diff_wlk.end(), [&](walker &wlk) {
        std::for_each(wlk.cords.begin(), wlk.cords.begin() + p->n_dims, [&](double &cord) {
            cord = movement_generator(movement_rng);
        });
    });
}

void DiffusionQuantumElectrons::apply_diffusion(electron_walker &ele_wlk, const electron_walker &diff_wlk) {
    for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
        for (int d = 0; d < p->n_dims; d++) {
            ele_wlk[wlk_idx].cords[d] += diff_wlk[wlk_idx].cords[d];
        }
    }
}
