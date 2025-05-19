#include "TrialFunctions/include/abstract_manybody_orbital.hpp"

#include "DiffusionParams/include/params.hpp"
#include "Core/include/walkers.hpp"
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

    std::unique_ptr<AbstractManybodyOrbital> orbital;

    double local_energy(const electron_walker& wlk);
    double p_value(const electron_walker& wlk, const electron_walker& prev_wlk, double growth_estimator);
private:

    electron_walker m_front_walker_buffer;
    electron_walker m_back_walker_buffer;
    
    std::unique_ptr<DiffusionWalkers> walkers_helper;

    double trial_wavef(const electron_walker& wlk);
    DiffusionQuantumParams* p;
};
