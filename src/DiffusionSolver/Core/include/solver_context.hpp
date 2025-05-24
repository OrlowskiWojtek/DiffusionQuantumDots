#include "DiffusionParams/include/harmonic_potential.hpp"
#include "TrialFunctions/include/abstract_manybody_orbital.hpp"

#include "Core/include/walkers.hpp"
#include "DiffusionParams/include/params.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <memory>

// solver context - data for enabling multithreading in loops
// every context should have yet alone:
// - RNGs
// - orbitals
// - potentially potentials
// now only orbital is needed for eval_p multithreading;

class SolverContext {
public:
    SolverContext();

    double local_energy(ElectronWalker &wlk);
    double p_value(ElectronWalker &wlk, ElectronWalker &prev_wlk, double growth_estimator);

    bool apply_nodes(ElectronWalker &wlk, const ElectronWalker &prev_wlk);
    void move_walkers(ElectronWalker &wlk, electron_walker &diff_value);
    void check_movement(ElectronWalker &wlk, ElectronWalker &prev_wlk, electron_walker &diff_value);
    void calc_trial_wavef(ElectronWalker &wlk);

private:
    void init_potential();
    void init_orbital();
    void init_rng();

    electron_walker drift_velocity;
    electron_walker m_front_walker_buffer;
    electron_walker m_back_walker_buffer;
    double green_diffusion_norm;

    static boost::random::mt19937 s_seed_generator;

    std::unique_ptr<DiffusionWalkers> walkers_helper;

    boost::random::mt19937 uni_rng;
    boost::random::uniform_real_distribution<double> uniform_generator;

    boost::random::mt19937 movement_rng;
    boost::random::normal_distribution<double> movement_generator;

    std::unique_ptr<HarmonicPotentialFunctor> potential;
    std::unique_ptr<AbstractManybodyOrbital> orbital;

    void apply_drift(ElectronWalker &wlk);
    void apply_diffusion(ElectronWalker &wlk, electron_walker &diffusion_values);
    void prepare_drift(const ElectronWalker &wlk);

    double trial_wavef(const electron_walker &wlk);
    double green_diffusion_term(const ElectronWalker &wlk, const ElectronWalker &prev_wlk);

    DiffusionQuantumParams *p;
};
