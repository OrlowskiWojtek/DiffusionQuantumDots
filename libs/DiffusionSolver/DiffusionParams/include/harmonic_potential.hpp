#ifndef HARMONIC_POTENTIAL_HPP
#define HARMONIC_POTENTIAL_HPP

#include "include/walkers_struct.hpp"
#include <functional>

struct HarmonicPotentialParams{
    int dims;
    double effective_mass;
    std::vector<double> omegas;
};

class HarmonicPotentialBuilder {
public:
    HarmonicPotentialBuilder() = default;

    void set_params(const HarmonicPotentialParams&);
    std::function<double(const walker&)> get_potential();

private:
    std::function<double(const walker&)> potential;
    std::vector<double> omegas_squared;
    HarmonicPotentialParams p;

    void fix_units();
    void build_potential();
};

#endif
