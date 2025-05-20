#ifndef HARMONIC_POTENTIAL_HPP
#define HARMONIC_POTENTIAL_HPP

#include "Core/include/walkers.hpp"
#include <functional>

struct HarmonicPotentialParams{
    int dims;
    double effective_mass;
    std::vector<double> omegas;
};

class HarmonicPotentialFunctor {
public:
    HarmonicPotentialFunctor() = default;
    HarmonicPotentialFunctor(HarmonicPotentialParams p);

    std::function<double(const walker&)> get_potential();
    double operator()(const walker&);

private:
    std::function<double(const walker&)> potential;
    std::vector<double> omegas_squared;
    HarmonicPotentialParams p;

    void precalculate();
    void build_potential();
};

#endif
