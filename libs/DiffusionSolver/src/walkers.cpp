#include "walkers.hpp"
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <numeric>

DiffusionWalkers::DiffusionWalkers() {}
DiffusionWalkers::~DiffusionWalkers() {}

void DiffusionWalkers::init_walkers(const DiffusionQuantumParams &params) {

    std::cout
        << "Starting initialization"
        << std::endl; // TODO: switch to formatted universal input / output

    walkers.resize(params.nmax_walkers);
    copy_walkers.resize(params.nmax_walkers);
    p_values.resize(params.nmax_walkers);

    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> initial_dist(params.xmin,
                                                            params.xmax);

    std::for_each(walkers.begin(), walkers.end(), [&](walker &wlk) {
        wlk.x = initial_dist(rng);
    });

    this->movement_generator = boost::random::normal_distribution<double>(
        0, std::sqrt(params.d_tau)); // TODO remember about d factor in 3d space

    d_tau = params.d_tau;
    V = params.pot;
}

void DiffusionWalkers::diffuse() {
    std::copy(
        this->walkers.begin(), this->walkers.end(), this->copy_walkers.begin());

    for (size_t i = 0; i < this->num_alive; i++) {
        this->walkers[i].x += this->movement_generator(rng);
    }
}

void DiffusionWalkers::eval_p() {
    for (size_t i = 0; i < this->num_alive; i++) {
        this->p_values[i] = std::exp(-d_tau * (V(walkers[i].x) + V(copy_walkers[i].x)) / 2. - growth_estimator);
    }
}

void DiffusionWalkers::branch() {
    for (size_t i = 0; i < this->num_alive; i++){
        int m = static_cast<int>(this->p_values[i] + this->uniform_generator(uni_rng));
        set_alive(m);
    }
}

void DiffusionWalkers::set_alive(int N) {
    
}
