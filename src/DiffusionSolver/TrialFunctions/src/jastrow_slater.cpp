#include "TrialFunctions/include/jastrow_slater.hpp"
#include "Core/include/walkers.hpp"
#include "TrialFunctions/include/harmonic_oscillator.hpp"
#include <iostream>
#include <memory>
#include <numeric>

JastrowSlaterOrbital::JastrowSlaterOrbital() { init_orbital(); }

JastrowSlaterOrbital::JastrowSlaterOrbital(JastrowSlaterOrbitalParams p)
    : p(p) {
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
    spins_up = std::accumulate(p.spins.begin(), p.spins.end(), 0, [](int acc, ElectronSpin spin){ return acc + (spin == ElectronSpin::UP ? 1 : 0);});
    spins_down = std::accumulate(p.spins.begin(), p.spins.end(), 0, [](int acc, ElectronSpin spin){ return acc + (spin == ElectronSpin::DOWN ? 1 : 0);});

    slater_matrix_up.set_size(spins_up, spins_up);
    slater_matrix_down.set_size(spins_down, spins_down);

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

double JastrowSlaterOrbital::operator()(const electron_walker &wlk) {
    double J = 0;

    for (int i = 0; i < p.electron_number; i++) {
        for (int j = i + 1; j < p.electron_number; j++) {
            double r_ij = distance(wlk[i], wlk[j]);
            if(r_ij < 1e-5){
                r_ij = 1e-5;
            }
            J -= 1. / 2. * p.a / r_ij * (1. - std::exp(-r_ij / p.b)) ;
        }
    }

    double jastrow_factor = std::exp(J);

    for (int i = 0; i < spins_up; i++) {
        for (int j = 0; j < spins_up; j++) {
            slater_matrix_up(i, j) = (*single_body_orbitals[i])(wlk[j]);
        }
    }

    for (int i = 0; i < spins_down; i++) {
        for (int j = 0; j < spins_down; j++) {
            slater_matrix_down(i, j) = (*single_body_orbitals[i + spins_up])(wlk[j + spins_up]);
        }
    }

    if (p.electron_number == 1) {
        return slater_matrix_up(0, 0);
    }

    if (spins_up == 2) {
        return jastrow_factor * 1 / std::sqrt(2) *
               (slater_matrix_up(0, 0) * slater_matrix_up(1, 1) - slater_matrix_up(0, 1) * slater_matrix_up(1, 0));
    }

    if (spins_up == 1 && spins_down == 1) {
        return jastrow_factor * (slater_matrix_up(0, 0) * slater_matrix_down(0,0));
    }

    return jastrow_factor * arma::det(slater_matrix_up) * arma::det(slater_matrix_down);
}

void JastrowSlaterOrbital::print() { std::cout << "Using Jastrow Slater based orbital" << std::endl; }

std::function<double(const electron_walker &wlk)> JastrowSlaterOrbital::get_orbital() { return orbital; }
