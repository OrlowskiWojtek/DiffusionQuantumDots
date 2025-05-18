#include "TrialFunctions/include/jastrow_slater.hpp"
#include "Core/include/walkers.hpp"
#include "TrialFunctions/include/harmonic_oscillator.hpp"
#include <iostream>
#include <memory>

JastrowSlaterOrbital::JastrowSlaterOrbital() { init_orbital(); }

JastrowSlaterOrbital::JastrowSlaterOrbital(JastrowSlaterOrbitalParams p): p(p) {
    init_orbital(); 
}

void JastrowSlaterOrbital::add_single_orbital(std::unique_ptr<AbstractSinglebodyOrbital> single_orbital) {
    single_body_orbitals.push_back(std::move(single_orbital));
}

double JastrowSlaterOrbital::distance(const walker &wlk_a, const walker &wlk_b) {
    double s = 0;
    for (int i = 0; i < p.dims; i++) {
        s += std::pow(wlk_a.cords[i] - wlk_b.cords[i], 2);
    }

    return std::sqrt(s);
}

void JastrowSlaterOrbital::init_orbital() {
    slater_matrix.set_size(p.electron_number, p.electron_number);

    single_body_orbitals.clear();
    for (int i = 0; i < p.electron_number; i++) {
        HarmonicOscillatorOrbitalsParams single_params; // TODO: try to use Aufbau principle in filling, now asume that
                                                        // in x direction is first excitement
        single_params.excitations = {i, 0, 0};
        single_params.omegas = p.omegas;
        single_params.effective_mass = p.effective_mass;
        single_params.dims = p.dims;

        single_body_orbitals.push_back(std::make_unique<HarmonicOscillatorOrbitals>(single_params));
    }

    orbital = [&](const electron_walker &wlk) {
        return (*this)(wlk);
    };

    print_test_to_file();
}

double JastrowSlaterOrbital::operator()(const electron_walker& wlk){
    double jastrow_factor = 1;

    for (int i = 0; i < p.electron_number; i++) {
        for (int j = i + 1; j < p.electron_number; j++) {
            double r_ij = distance(wlk[i], wlk[j]);
            jastrow_factor *= std::exp(p.a * r_ij / (1. + p.b * r_ij));       
        }
    }

    for (int i = 0; i < p.electron_number; i++) {
        for (int j = 0; j < p.electron_number; j++) {
            slater_matrix(i, j) = (*single_body_orbitals[i])(wlk[j]);
        }
    }

    return jastrow_factor * arma::det(slater_matrix);
}

void JastrowSlaterOrbital::print(){
    std::cout << "Using Jastrow Slater based orbital" << std::endl;
}

std::function<double(const electron_walker& wlk)> JastrowSlaterOrbital::get_orbital(){
    return orbital;
}
