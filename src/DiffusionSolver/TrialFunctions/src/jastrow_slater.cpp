#include "TrialFunctions/include/jastrow_slater.hpp"
#include "Core/include/walkers_struct.hpp"
#include "TrialFunctions/include/harmonic_oscillator.hpp"
#include <memory>

JastrowSlaterOrbital::JastrowSlaterOrbital() { init_orbital(); }

void JastrowSlaterOrbital::add_single_orbital(std::unique_ptr<AbstractSinglebodyOrbital> single_orbital) {
    single_body_orbitals.push_back(std::move(single_orbital));
}

double JastrowSlaterOrbital::distance(const walker &wlk_a, const walker &wlk_b) {
    double s = 0;
    for (int i = 0; i < p->dims; i++) {
        s += std::pow(wlk_a.cords[i] - wlk_b.cords[i], 2);
    }

    return std::sqrt(s);
}

void JastrowSlaterOrbital::init_orbital() {
    slater_matrix.resize(p->electron_number, p->electron_number);
    for (int i = 0; i < p->electron_number; i++) {
        HarmonicOscillatorOrbitalsParams single_params; // TODO: try to use Aufbau principle in filling, now asume that
                                                        // in x direction is first excitement
        single_params.excitations = {i, 0, 0};
        single_params.omegas = p->omegas;
        single_params.effective_mass = p->effective_mass;
        single_params.dims = p->dims;

        add_single_orbital(std::make_unique<HarmonicOscillatorOrbitals>(single_params));
    }

    orbital = [&](const electron_walker &wlk) {
        return (*this)(wlk);
    };
}

double JastrowSlaterOrbital::operator()(const electron_walker& wlk){
    double jastrow_factor = 0;

    for (int i = 0; i < p->electron_number; i++) {
        for (int j = 0; j < p->electron_number; j++) {
            double r_ij = distance(wlk[i], wlk[j]);
            jastrow_factor += p->a * r_ij / (1. + p->b * r_ij);
        }
    }

    for (int i = 0; i < p->electron_number; i++) {
        for (int j = 0; j < p->electron_number; j++) {
            slater_matrix(i, j) = (*single_body_orbitals[i])(wlk[i]);
        }
    }

    return jastrow_factor * slater_matrix.determinant();
}
