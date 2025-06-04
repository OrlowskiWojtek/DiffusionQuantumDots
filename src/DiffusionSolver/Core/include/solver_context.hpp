#include "DiffusionParams/include/harmonic_potential.hpp"
#include "TrialFunctions/include/abstract_manybody_orbital.hpp"

#include "Core/include/walkers.hpp"
#include "DiffusionParams/include/params.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <memory>

 /*! 
  * SolverContext is class containing basic calculations used in algorithm.
  * Local energy and movements are calculated here.
  * Main goal of class is to provide methods and elements for each thread
  * in order to avoid data racing.
  * Each context should have independent:
  * - random number generators with different initial seed
  * - trial wavefunction instance
  * - potential instance
  * - buffers used in calculations to avoid unnecessery allocation
*/
class SolverContext {
public:
    SolverContext();

    /** 
     * Calculates local energy and saves it inside @param wlk.
     * The local energy is defined as E_L(R) = (Hψᵀ(R))/ψᵀ(R).
    */
    void calc_local_energy(ElectronWalker &wlk);

    // Calculates p_value - probability of reproduction in branching step for 
    // @param wlk - current walker 
    // @param prev_wlk - walker in previous step
    // @param growth_estimator - current growth estimator
    double p_value(ElectronWalker &wlk, ElectronWalker &prev_wlk, double growth_estimator);

    // Checks if walker crossed node.
    //
    // @return true if walker crossed node
    bool check_nodes(ElectronWalker &wlk, const ElectronWalker &prev_wlk);

    // Manages whole movement step of single walker:
    // - random diffusion movement
    // - drift movement
    //
    // @param wlk - walker to move
    // @param diff_value - diffusion part of step - saved for quick estimation in acceptance step \see check_movement
    void move_walkers(ElectronWalker &wlk, electron_walker &diff_value);

    // Manages Metropolis acceptance
    //
    // @return true if movement has been accepted, otherwise returns false and moves walker back
    bool check_metropolis(ElectronWalker &wlk, ElectronWalker &prev_wlk, electron_walker &diff_value);

    // Calculates trial wavefunction value for walker
    // saves calculated value inside walker in order to avoid 
    // multiple calculations 
    void calc_trial_wavef(ElectronWalker &wlk);

    double get_potential(const ElectronWalker& wlk);

private:
    void init_potential();
    void init_orbital();
    void init_rng();

    // buffer used for storing calculated drift velocity
    electron_walker drift_velocity;

    // buffer used for derivative calculation using finite differences
    electron_walker m_front_walker_buffer;

    // buffer used for derivative calculation using finite differences
    electron_walker m_back_walker_buffer;

    // Static seed generator used to provide 
    // different contexts with different seed inside generators
    static boost::random::mt19937 s_seed_generator;

    // provides basic walkers helper functions
    std::unique_ptr<DiffusionWalkers> walkers_helper;

    boost::random::mt19937 uni_rng;
    boost::random::uniform_real_distribution<double> uniform_generator;

    boost::random::mt19937 movement_rng;
    boost::random::normal_distribution<double> movement_generator;

    std::unique_ptr<HarmonicPotentialFunctor> potential;
    std::unique_ptr<AbstractManybodyOrbital> orbital;

    // applies drift in movement step
    void apply_drift(ElectronWalker &wlk);

    // applies diffusion in movement step
    void apply_diffusion(ElectronWalker &wlk, electron_walker &diffusion_values);

    // preperes drift before diffusion
    void prepare_drift(const ElectronWalker &wlk);

    // calculates trial_wavefunction for wlk
    double trial_wavef(const electron_walker &wlk);

    // calculates local energy for wlk
    double local_energy(const ElectronWalker &wlk);

    // calculates Green function diffusion part (G_d(R'<-R))
    // used in metropolis back step propability evaluation
    double green_diffusion_term(const ElectronWalker &wlk, const ElectronWalker &prev_wlk);

    // parameters of simulation
    DiffusionQuantumParams *p;

    // Kinetic term in Hamiltonian
    double kinetic_term(const ElectronWalker& wlk);

    // Potential term in Hamiltonian
    double potential_term(const ElectronWalker& wlk);

    // Interaction term in Hamiltonian
    double interaction_term(const ElectronWalker& wlk);
};
