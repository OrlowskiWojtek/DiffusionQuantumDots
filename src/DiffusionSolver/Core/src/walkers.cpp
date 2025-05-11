#include "Core/include/walkers.hpp"
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <ctime>

DiffusionWalkers::DiffusionWalkers()
    : p(DiffusionQuantumParams::getInstance()) {

    trial_wavef = p->trial_wavef->get_orbital();
    V = p->pot;
    movement_generator = boost::random::normal_distribution<double>(0, std::sqrt(p->d_tau));
}

DiffusionWalkers::~DiffusionWalkers() {}

void DiffusionWalkers::apply_diffusion(walker &wlk) {
    std::for_each(wlk.cords.begin(), wlk.cords.begin() + p->n_dims, [&](double &cord) {
        cord += movement_generator(rng);
    });
}

// TODO: revise if it is not necessery to move all walkers (from all electron_walkers) first and then calc derivative
double DiffusionWalkers::local_energy(const walker &wlk) {
    double dr = 1e-6;
    double dr2 = 1e-12;
    double kinetic_term = 0;
    double cent_value = trial_wavef(wlk);

    for (int d = 0; d < p->n_dims; d++) {
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

double DiffusionWalkers::trial_wf_value(const walker &wlk) { 
    return trial_wavef(wlk); 
}

// TODO after upgrade libstdc++ update to views::zip with std::acumulate
double DiffusionWalkers::distance(const walker &wlk_a, const walker &wlk_b) {
    double s = 0;
    for(int i = 0; i < p->n_dims; i++){
        s += std::pow(wlk_a.cords[i] - wlk_b.cords[i], 2);
    }

    return std::sqrt(s);
}
