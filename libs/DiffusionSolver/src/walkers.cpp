#include "walkers.hpp"
#include "params.hpp"
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>

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

    walkers.resize(params.nmax_walkers);
    copy_walkers.resize(params.nmax_walkers);
    p_values.resize(params.nmax_walkers);
    hist.resize(n_bins);

    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> initial_dist(xmin, xmax);

    std::for_each(walkers.begin(), walkers.end(), [&](walker &wlk) { wlk.x = initial_dist(rng); });

    movement_generator = boost::random::normal_distribution<double>(
        0, std::sqrt(params.d_tau)); // TODO remember about d factor in 3d space

    Et = 0;
    current_it = 1;

    results->init_x(xmin, xmax, n_bins);
}

void DiffusionWalkers::diffuse() {
    std::copy(walkers.begin(), walkers.end(), this->copy_walkers.begin());

    for (size_t i = 0; i < num_alive; i++) {
        walkers[i].x += movement_generator(rng);
    }
}

void DiffusionWalkers::eval_p() {
    for (size_t i = 0; i < num_alive; i++) {
        p_values[i] =
            std::exp(-d_tau * ((V(walkers[i].x) + V(copy_walkers[i].x)) / 2. - growth_estimator));
    }
}

void DiffusionWalkers::branch() {
    new_alive = 0;

    for (size_t i = 0; i < num_alive; i++) {
        int m = static_cast<int>(p_values[i] + uniform_generator(uni_rng));
        set_alive(m, walkers[i].x);
    }

    num_alive = new_alive;
    std::copy(copy_walkers.begin(),
              copy_walkers.end(),
              walkers.begin()); // TODO : optimize so it needs no copy

    update_growth_estimator();
    current_it++;
}

void DiffusionWalkers::set_alive(int N, double x) {
    for (size_t i = new_alive; i < new_alive + N; i++) {
        if (i >= walkers.size()) {

            std::cout << "Error: The programme ran out of allocated memory. "
                         "Required resize."
                      << std::endl;
            walkers.emplace_back();
            walkers.back().x = x; // TODO -> resize also other containers like p_i and copy_walkers

            continue;
        }
        copy_walkers[i].x = x;
    }

    new_alive += N;
}

void DiffusionWalkers::update_growth_estimator() {
    double E_t = 0;

    for (size_t i = 0; i < num_alive; i++) {
        E_t += V(walkers[i].x);
    }

    Et += E_t / static_cast<double>(num_alive);

    growth_estimator =
        Et / static_cast<double>(current_it) -
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

    if (current_it == 100000) {
        save_progress();
    }
}

void DiffusionWalkers::save_progress() {
    results->add_histogram(static_cast<double>(current_it) * d_tau, current_it, hist);
}
