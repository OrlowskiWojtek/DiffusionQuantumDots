#ifndef HARMONIC_OSCILATOR_SOLUTIONS_HPP
#define HARMONIC_OSCILATOR_SOLUTIONS_HPP

#include "TrialFunctions/include/abstract_wf.hpp"
#include <functional>

// exact solutions of schrodinger equation in harmonic oscillator
struct HarmonicOscillatorOrbitalsParams {
    int dims;
    double effective_mass;
    std::vector<double> omegas;
    std::vector<int> excitations;
};

class HarmonicOscillatorOrbitals : public AbstractOrbital {
public:
    HarmonicOscillatorOrbitals();
    HarmonicOscillatorOrbitals(HarmonicOscillatorOrbitalsParams p);

    double operator()(const walker &wlk) override;
    void print() override;

    std::function<double(const walker& wlk)> get_orbital() override;
private:
    void print_test_to_file() override;

    HarmonicOscillatorOrbitalsParams p;

    double first_part_precalculated;
    double second_part_precalculated;
    double third_part_precalculated;

    void init_orbital();
    std::function<double(const walker &)> orbital;

    int factorial(const int& n);
};

#endif
