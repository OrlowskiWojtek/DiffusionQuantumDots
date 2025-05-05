#include "DiffusionParams/include/harmonic_potential.hpp"
#include <cmath>

void HarmonicPotentialBuilder::fix_units() {
    omegas_squared.resize(p.omegas.size());
    std::transform(p.omegas.begin(), p.omegas.end(), omegas_squared.begin(), [](double om){ return std::pow(om, 2);});
}

void HarmonicPotentialBuilder::build_potential() {
    switch (p.dims) {
    case 1:
        potential = [&](const walker &wlk) {
            return 0.5 / p.effective_mass * (omegas_squared[0] * std::pow(wlk.cords[0], 2));
        };
        break;
    case 2:
        potential = [&](const walker &wlk) {
            return 0.5 / p.effective_mass *
                   (omegas_squared[0] * std::pow(wlk.cords[0], 2) + omegas_squared[1] * std::pow(wlk.cords[1], 2));
        };
        break;
    case 3:
        potential = [&](const walker &wlk) {
            return 0.5 / p.effective_mass *
                   (omegas_squared[0] * std::pow(wlk.cords[0], 2) + omegas_squared[1] * std::pow(wlk.cords[1], 2) +
                    omegas_squared[2] * std::pow(wlk.cords[2], 1));
        };
    default:
        potential = [&](const walker &wlk) {
            return 0.;
        };
        break;
    }
}

std::function<double(const walker &)> HarmonicPotentialBuilder::get_potential() { 
    fix_units();
    build_potential();
    return potential;
}

void HarmonicPotentialBuilder::set_params(const HarmonicPotentialParams& p){
    this->p = p;
}
