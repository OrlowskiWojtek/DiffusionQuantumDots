#include "TrialFunctions/include/harmonic_oscillator.hpp"
#include "include/UnitHandler.hpp"
#include <boost/math/special_functions/hermite.hpp>
#include <fstream>
#include <functional>
#include <iostream>

HarmonicOscillatorOrbitals::HarmonicOscillatorOrbitals() {}

HarmonicOscillatorOrbitals::HarmonicOscillatorOrbitals(HarmonicOscillatorOrbitalsParams p)
    : p(p) {
    init_orbital();
}

// precalculating values
void HarmonicOscillatorOrbitals::init_orbital() {
    first_part_precalculated = 1;
    for (int i = 0; i < p.dims; i++) {
        first_part_precalculated *= 1. / sqrt(std::pow(2, p.excitations[i]) * factorial(p.excitations[i]));
    }

    second_part_precalculated = 1;
    for (int i = 0; i < p.dims; i++) {
        second_part_precalculated *= std::pow(p.effective_mass * p.omegas[i] / M_PI, 0.25);
    }

    orbital = [&](const walker &wlk) {
        double third_part = 1;
        for (int i = 0; i < p.dims; i++) {
            third_part *= std::exp(-p.effective_mass * p.omegas[i] * std::pow(wlk.cords[i], 2) / 2.);
        }

        double fourth_part = 1;
        for (int i = 0; i < p.dims; i++) {
            fourth_part *= boost::math::hermite(p.excitations[i], std::sqrt(p.effective_mass * p.omegas[i]) * wlk.cords[i]);
        }

        return first_part_precalculated * second_part_precalculated * third_part * fourth_part;
    };

    print_test_to_file();
}

int HarmonicOscillatorOrbitals::factorial(const int &n) {
    int f = 1;
    for (int i = 1; i <= n; ++i)
        f *= i;
    return f;
}

// TODO: optimization in precalculations of those values
double HarmonicOscillatorOrbitals::operator()(const walker &wlk) {
    double third_part = 1;
    for (int i = 0; i < p.dims; i++) {
        third_part *= std::exp(-p.effective_mass * p.omegas[i] * std::pow(wlk.cords[i], 2) / 2.);
    }

    double fourth_part = 1;
    for (int i = 0; i < p.dims; i++) {
        fourth_part *= boost::math::hermite(p.excitations[i], std::sqrt(p.effective_mass * p.omegas[i]) * wlk.cords[i]);
    }

    return first_part_precalculated * second_part_precalculated * third_part * fourth_part;
}

std::function<double(const walker &)> HarmonicOscillatorOrbitals::get_orbital() { return orbital; }

void HarmonicOscillatorOrbitals::print() {
    std::cout << "Using trial wavefunction: \n";
    std::cout << "Exact solution of harmonic oscillator at "
        << p.excitations[0]
        << p.excitations[1] 
        << p.excitations[2] << " state\n";
    std::cout << std::endl;
}

// Simple 2D test -> TODO
void HarmonicOscillatorOrbitals::print_test_to_file() {
    double xmin = UnitHandler::length(UnitHandler::TO_AU, -20.);
    double xmax = UnitHandler::length(UnitHandler::TO_AU, 20.);
    double dx = (xmax - xmin) / 100.;

    std::ofstream file("TrialWavefunctionTest");

    walker test_walker;
    test_walker.cords[2] = 0;

    for (double i = xmin; i < xmax; i += dx) {
        for (double j = xmin; j < xmax; j += dx) {
            test_walker.cords[0] = i;
            test_walker.cords[1] = j;
            file << (*this)(test_walker) << "\t";
        }
        file << "\n";
    }
    file << std::endl;
    file.close();
}
