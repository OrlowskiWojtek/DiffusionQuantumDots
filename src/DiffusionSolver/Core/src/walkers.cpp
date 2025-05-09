#include "Core/include/walkers.hpp"
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <ctime>
#include <execution>
#include <numeric>

DiffusionWalkers::DiffusionWalkers()
    : p(DiffusionQuantumParams::getInstance())
    , results(std::make_unique<DiffusionQuantumResults>()) {}

DiffusionWalkers::~DiffusionWalkers() {}

void DiffusionWalkers::init_walkers() { // TODO: segmentize this function
    xmin = p->xmin;
    xmax = p->xmax;
    d_tau = p->d_tau;
    V = p->pot;
    num_alive = p->n0_walkers;
    target_alive = p->n0_walkers;
    n_bins = p->n_bins;
    calibrating = p->blocks_calibration;
    block_size = p->n_block;
    dims = p->n_dims;
    trial_wavef = p->trial_wavef->get_orbital();

    // initialize containers
    walkers.resize(p->nmax_walkers);
    copy_walkers.resize(p->nmax_walkers);
    p_values.resize(p->nmax_walkers);
    hist.resize(boost::extents[n_bins][n_bins][n_bins]);
    std::for_each(hist.data(), hist.data() + hist.num_elements(), [](int64_t &val) {
        val = 0;
    });

    boost::random::mt19937 rng;
    rng.seed(std::time(0));
    boost::random::uniform_real_distribution<> initial_dist(xmin, xmax);

    std::for_each(walkers.begin(), walkers.end(), [&](walker &wlk) {
        wlk.cords.fill(0);
    });

    std::for_each(walkers.begin(), walkers.end(), [&](walker &wlk) {
        std::for_each(wlk.cords.begin(), wlk.cords.begin() + dims, [&](double &cord) {
            cord = initial_dist(rng);
        });
    });

    std::copy(walkers.begin(), walkers.end(), copy_walkers.begin());

    movement_generator = boost::random::normal_distribution<double>(0, std::sqrt(p->d_tau));

    growth_estimator = 0;
    Eblock = 0;
    ground_state_estimator = 0.;
    ground_mean = 0.;

    current_it = 0;
    accumulation_it = 0;
    blocks_passed = 0;
}

void DiffusionWalkers::diffuse() {
    current_it++;

    std::for_each(walkers.begin(), walkers.begin() + num_alive, [&](walker &wlk) {
        apply_drift(wlk);
        apply_diffusion(wlk);
    });
}

void DiffusionWalkers::apply_drift(walker &wlk) {
    std::array<double, 3> drift_force = drift(wlk);
    std::transform(wlk.cords.begin(),
                   wlk.cords.begin() + dims,
                   drift_force.begin(),
                   wlk.cords.begin(),
                   [&](double &cord, double &drift_value) {
                       return cord + d_tau * drift_value;
                   });
}

void DiffusionWalkers::apply_diffusion(walker &wlk) {
    std::for_each(wlk.cords.begin(), wlk.cords.begin() + dims, [&](double &cord) {
        cord += movement_generator(rng);
    });
}

void DiffusionWalkers::eval_p() {
    std::transform(std::execution::par, // smarter threading is needed
                   walkers.begin(),
                   walkers.begin() + num_alive,
                   copy_walkers.begin(),
                   p_values.begin(),
                   [&](walker &wlk, walker &prev_wlk) {
                       if (apply_nodes(wlk, prev_wlk)) {
                           return 0.;
                       }

                       return p_value(wlk, prev_wlk);
                   });
}

double DiffusionWalkers::p_value(const walker &wlk, const walker &prev_wlk) {
    return std::exp(-d_tau * ((local_energy(wlk) + local_energy(prev_wlk)) / 2. - growth_estimator));
}

bool DiffusionWalkers::apply_nodes(const walker &wlk, const walker &prev_wlk) {
    return (trial_wavef(wlk) * trial_wavef(prev_wlk)) <= 0;
}

void DiffusionWalkers::reject_move(walker &wlk, walker &prev_wlk) { wlk = prev_wlk; }

void DiffusionWalkers::branch() {
    new_alive = 0;

    for (int i = 0; i < num_alive; i++) {
        int m = static_cast<int>(p_values[i] + uniform_generator(uni_rng));
        set_alive(m, walkers[i]);
    }

    num_alive = new_alive;
    std::copy(copy_walkers.begin(), copy_walkers.begin() + num_alive, walkers.begin());

    update_growth_estimator();
}

void DiffusionWalkers::set_alive(int N, const walker &wlk) {
    for (int i = new_alive; i < new_alive + N; i++) {
        if (i >= static_cast<int>(walkers.size())) {
            std::cout << "Error: The programme ran out of allocated memory. "
                         "Required resize."
                      << std::endl;
            break;
        }
        copy_walkers[i] = wlk;
    }

    new_alive += N;
}

void DiffusionWalkers::update_growth_estimator() {
    if (accumulation_it == 0) {
        int nblock = current_it % block_size;
        Eblock = (growth_estimator + Eblock * nblock) / (nblock + 1);
    } else {
        Eblock = (growth_estimator + Eblock * accumulation_it) / (accumulation_it + 1);
    }

    growth_estimator =
        Eblock - 1. / d_tau * std::log(static_cast<double>(num_alive) / static_cast<double>(target_alive));
}

void DiffusionWalkers::binning() {
    std::for_each(walkers.begin(), walkers.begin() + num_alive, [&](const walker &wlk) {
        for (int dim = 0; dim < dims; dim++) {
            if (wlk.cords[dim] < xmin || wlk.cords[dim] > xmax) {
                return;
            }
        }

        std::array<size_t, 3> walker_bin;
        walker_bin[0] = static_cast<int>((xmax - wlk.cords[0]) / (xmax - xmin) * n_bins);
        walker_bin[1] = static_cast<int>((xmax - wlk.cords[1]) / (xmax - xmin) * n_bins);
        walker_bin[2] = static_cast<int>((xmax - wlk.cords[2]) / (xmax - xmin) * n_bins);
        hist[walker_bin[0]][walker_bin[1]][walker_bin[2]]++;
    });
}

void DiffusionWalkers::count() {
    binning();

    current_Et = std::accumulate(walkers.begin(),
                                 walkers.begin() + num_alive,
                                 0.,
                                 [this](double acc, const walker &wlk) {
                                     return acc + local_energy(wlk);
                                 }) /
                 static_cast<double>(num_alive);

    ground_state_estimator += current_Et;
    accumulation_it++;

    if (calibrating) {
        results->add_energies(current_Et / static_cast<double>(num_alive), Eblock);
        return;
    }

    if (static_cast<size_t>(accumulation_it) % block_size == 0) {
        ground_mean += ground_state_estimator / static_cast<double>(num_alive) / block_size;
        ground_mean2 += std::pow(ground_state_estimator / static_cast<double>(num_alive) / block_size, 2);

        ground_state_estimator = 0;
        blocks_passed += 1;
    }
}

void DiffusionWalkers::save_progress() {
    results->add_histogram(
        static_cast<double>(current_it) * d_tau, current_it, ground_state_estimator / accumulation_it, Eblock, hist);
}

DiffusionQuantumResults &DiffusionWalkers::get_results() { return *results; }

double DiffusionWalkers::local_energy(const walker &wlk) {
    double dr = 1e-6;
    double dr2 = 1e-12;
    double kinetic_term = 0;
    double cent_value = trial_wavef(wlk);

    for (int d = 0; d < dims; d++) {
        walker back_walker(wlk.cords[0], wlk.cords[1], wlk.cords[2]);
        walker forw_walker(wlk.cords[0], wlk.cords[1], wlk.cords[2]);

        back_walker.cords[d] -= dr;
        forw_walker.cords[d] += dr;
        double back_value = trial_wavef(back_walker);
        double forw_value = trial_wavef(forw_walker);

        kinetic_term += (back_value - 2 * cent_value + forw_value);
    }

    kinetic_term = -0.5 * kinetic_term / (cent_value * dr2);

    return kinetic_term + V(wlk);
}

std::array<double, 3> DiffusionWalkers::drift(const walker &wlk) {
    std::array<double, 3> dr_force;

    double dr = 1e-6;
    double cent_value = trial_wavef(wlk);

    for (int d = 0; d < dims; d++) {
        walker forward_walker(wlk.cords[0], wlk.cords[1], wlk.cords[2]);
        forward_walker.cords[d] += dr;

        double forw_value = trial_wavef(forward_walker);
        dr_force[d] = (forw_value - cent_value) / (dr * cent_value);
    }

    return dr_force;
}

double DiffusionWalkers::trial_wf_value(const walker &wlk) { return trial_wavef(wlk); }

// TODO after upgrade libstdc++ update to views::zip with std::acumulate
double DiffusionWalkers::distance(const walker &wlk_a, const walker &wlk_b) {
    double s = 0;
    for(int i = 0; i < dims; i++){
        s += std::pow(wlk_a.cords[i] - wlk_b.cords[i], 2);
    }

    return std::sqrt(s);
}
