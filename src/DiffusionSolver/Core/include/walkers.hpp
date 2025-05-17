#ifndef DIFFUSION_WALKERS_HPP
#define DIFFUSION_WALKERS_HPP

#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

#include "DiffusionParams/include/params.hpp"
#include "Core/include/walkers_struct.hpp"

#include <boost/random/uniform_real_distribution.hpp>

class DiffusionWalkers {
public:
    DiffusionWalkers();
    ~DiffusionWalkers();

    void apply_diffusion(walker &wlk);
    // TODO: void prepare_diffusion(walker &wlk);

    double trial_wf_value(const walker& wlk);
    double local_energy(const walker &wlk);
    double distance(const walker& wlk_a, const walker& wlk_b);

private:
    DiffusionQuantumParams *p;

    boost::random::mt19937 rng;
    boost::random::normal_distribution<double> movement_generator;

    std::function<double(const walker &)> V;
    std::function<double(const walker &)> trial_wavef;

    void reject_move(walker &wlk, walker &prev_wlk);
};

#endif
