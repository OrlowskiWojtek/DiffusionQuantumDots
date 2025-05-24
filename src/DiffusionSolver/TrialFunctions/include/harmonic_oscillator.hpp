#ifndef HARMONIC_OSCILATOR_SOLUTIONS_HPP
#define HARMONIC_OSCILATOR_SOLUTIONS_HPP

#include "TrialFunctions/include/abstract_singlebody_orbital.hpp"
#include <functional>

// exact solutions of schrodinger equation in harmonic oscillator
struct HarmonicOscillatorOrbitalsParams {
    int dims;
    double effective_mass;
    std::vector<double> omegas;
    std::vector<int> excitations;
};

class HarmonicOscillatorOrbitals : public AbstractSinglebodyOrbital {
public:
    HarmonicOscillatorOrbitals();
    HarmonicOscillatorOrbitals(HarmonicOscillatorOrbitalsParams p);

    double operator()(const walker &wlk) override;
    void print() override;

    std::function<double(const walker& wlk)> get_orbital() override;
private:

    HarmonicOscillatorOrbitalsParams p;

    double first_part_precalculated;
    double second_part_precalculated;
    double first_second_part_precalculated;

    std::vector<double> precalculated_eff_mass_omegas;
    std::vector<double> sqrt_precalculated_eff_mass_omegas;

    void init_orbital();
    std::function<double(const walker &)> orbital;

    int factorial(const int& n);
};

#endif
