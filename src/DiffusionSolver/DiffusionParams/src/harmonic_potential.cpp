#include "DiffusionParams/include/harmonic_potential.hpp"
#include <cmath>

HarmonicPotentialFunctor::HarmonicPotentialFunctor(HarmonicPotentialParams p)
    : p(p) {
    precalculate();
    build_potential();
}

void HarmonicPotentialFunctor::precalculate() {
    omegas_squared.resize(p.omegas.size());
    std::transform(p.omegas.begin(), p.omegas.end(), omegas_squared.begin(), [](double om) {
        return std::pow(om, 2);
    });
}

void HarmonicPotentialFunctor::build_potential() {
    switch (p.dims) {
    case 1:
        potential = [&](const walker &wlk) {
            return 0.5 * p.effective_mass * (omegas_squared[0] * std::pow(wlk.cords[0], 2));
        };
        break;
    case 2:
        potential = [&](const walker &wlk) {
            return 0.5 * p.effective_mass *
                   (omegas_squared[0] * std::pow(wlk.cords[0], 2) +
                    omegas_squared[1] * std::pow(wlk.cords[1], 2));
        };
        break;
    case 3:
        potential = [&](const walker &wlk) {
            return 0.5 * p.effective_mass *
                   (omegas_squared[0] * std::pow(wlk.cords[0], 2) +
                    omegas_squared[1] * std::pow(wlk.cords[1], 2) +
                    omegas_squared[2] * std::pow(wlk.cords[2], 2));
        };
        break;
    default:
        potential = [&](const walker &wlk) {
            return 0.;
        };
        break;
    }
}

std::function<double(const walker &)> HarmonicPotentialFunctor::get_potential() {
    return potential;
}

double HarmonicPotentialFunctor::operator()(const walker &wlk) { return potential(wlk); }
