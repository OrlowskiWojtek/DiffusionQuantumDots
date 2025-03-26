#include "walkers.hpp"
#include "params.hpp"
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <numeric>

DiffusionWalkers::DiffusionWalkers()
    : results(std::make_unique<DiffusionQuantumResults>()) {}
DiffusionWalkers::~DiffusionWalkers() {}

void DiffusionWalkers::init_walkers(const DiffusionQuantumParams &params) {

    xmin = params.xmin;
    xmax = params.xmax;
    d_tau = params.d_tau;
    V = params.pot;
    num_alive = params.n0_walkers;
    target_alive = params.n0_walkers;
    n_bins = params.n_bins;
    calibrating = params.blocks_calibration;
    block_size = params.n_block;
    nodes = params.nodes;

    walkers.resize(params.nmax_walkers);
    copy_walkers.resize(params.nmax_walkers);
    p_values.resize(params.nmax_walkers);
    hist.resize(n_bins);

    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> initial_dist(xmin, xmax);

    std::for_each(walkers.begin(), walkers.end(), [&](walker &wlk) { wlk.x = initial_dist(rng); });
    std::copy(walkers.begin(), walkers.end(), copy_walkers.begin());

    movement_generator = boost::random::normal_distribution<double>(
        0, std::sqrt(params.d_tau)); // TODO remember about d factor in 3d space

    growth_estimator = 0;
    Eblock = 0;
    ground_state_estimator = 0.;

    current_it = 0;
    accumulation_it = 0;
    blocks_passed = 0;

    results->init_x(xmin, xmax, n_bins);
}

void DiffusionWalkers::diffuse() {
    current_it++;

    for (size_t i = 0; i < num_alive; i++) {
        walkers[i].x += movement_generator(rng);
    }
}

void DiffusionWalkers::eval_p() {
    for (size_t i = 0; i < num_alive; i++) {
        if (apply_nodes(i))
            continue;
        p_values[i] =
            std::exp(-d_tau * ((V(walkers[i].x) + V(copy_walkers[i].x)) / 2. - growth_estimator));
    }
}

bool DiffusionWalkers::apply_nodes(int i) {
    for (const double &node : nodes) {
        if ((walkers[i].x - node) * (copy_walkers[i].x - node) < 0) {
            p_values[i] = 0;
            return true;
        }
    }

    return false;
}

void DiffusionWalkers::branch() {
    new_alive = 0;

    for (size_t i = 0; i < num_alive; i++) {
        int m = static_cast<int>(p_values[i] + uniform_generator(uni_rng));
        set_alive(m, walkers[i].x);
    }

    num_alive = new_alive;
    std::copy(copy_walkers.begin(),
              copy_walkers.begin() + num_alive,
              walkers.begin());

    update_growth_estimator();
}

void DiffusionWalkers::set_alive(int N, double x) {
    for (size_t i = new_alive; i < new_alive + N; i++) {
        if (i >= walkers.size()) {
            std::cout << "Error: The programme ran out of allocated memory. "
                        "Required resize."
                      << std::endl;
            continue;
        }
        copy_walkers[i].x = x;
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
        Eblock -
        1. / d_tau * std::log(static_cast<double>(num_alive) / static_cast<double>(target_alive));
}

void DiffusionWalkers::count() {
    for (size_t i = 0; i < num_alive; i++) {
        if (walkers[i].x < xmin || walkers[i].x > xmax) {
            continue;
        }
        int bin = static_cast<int>((xmax - walkers[i].x) / (xmax - xmin) * n_bins);
        hist[bin]++;
    }

    current_Et = std::accumulate(walkers.begin(),
                                 walkers.begin() + num_alive,
                                 0.,
                                 [this](double acc, const walker &wlk) { return acc + V(wlk.x); });

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
    results->add_histogram(static_cast<double>(current_it) * d_tau, current_it, ground_mean / blocks_passed, Eblock, hist);
}

DiffusionQuantumResults &DiffusionWalkers::get_results() { return *results; }
