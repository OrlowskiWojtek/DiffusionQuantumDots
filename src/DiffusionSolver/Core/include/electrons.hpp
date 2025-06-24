#ifndef DIFFUSION_QUANTUM_ELECTRONS_HPP
#define DIFFUSION_QUANTUM_ELECTRONS_HPP

#include "Core/include/results.hpp"
#include "Core/include/solver_context.hpp"
#include "Core/include/walkers.hpp"
#include "DiffusionParams/include/params.hpp"
#include <Core/include/walkers.hpp>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

/**
 * Diffusion Quantum Electrons class header
 *
 * This class manages the ensemble of electron walkers in the DMC simulation.
 * It handles the core operations of diffusion, branching, and data collection
 * for the electron walkers representing the quantum system.
 *
 * The implementation follows the mathematical framework of DMC where walkers
 * evolve according to the imaginary-time Schrödinger equation:
 * ∂Ψ(R,τ)/∂τ = (1/2)∇²Ψ(R,τ) - (V(R) - E_ref)Ψ(R,τ)
 */
class DiffusionQuantumElectrons {
public:
    DiffusionQuantumElectrons();

    /**
     * Diffusion Quantum Electrons class header
     *
     * This class manages the ensemble of electron walkers in the DMC simulation.
     * It handles the core operations of diffusion, branching, and data collection
     * for the electron walkers representing the quantum system.
     *
     * The implementation follows the mathematical framework of DMC where walkers
     * evolve according to the imaginary-time Schrödinger equation:
     * ∂Ψ(R,τ)/∂τ = (1/2)∇²Ψ(R,τ) - (V(R) - E_ref)Ψ(R,τ)
     */
    void diffuse();

    /**
     * Performs the branching step for all electron walkers
     *
     * This implements the branching/death term -(V(R) - E_ref)Ψ(R,τ) by
     * replicating or removing walkers based on their weight:
     * w = exp[-(V(R) - E_ref)τ]
     */
    void branch();

    /**
     * Checks nodes crossing and metropolis acceptance
     *
     * Calculates acceptance rate that is needed for effective time calculations
     */
    void check_movement();

    /**
     * Prepares after movement walkers for branching step.
     *
     * Accepts or rejects proposed walker movements based on the trial wavefunction
     *
     * This implements the Fixed-Node Approximation by killing walkers that would
     * cross nodal surfaces where ψᵀ(R)ψᵀ(R') < 0. For non-node-crossing moves,
     * the acceptance probability is min[1, |Gd(R'<-R) ψᵀ(R')|²/(Gd(R<-R') |ψᵀ(R)|²)].
     *
     * Calculates the weight w = exp[-(V(R) - E_ref)τ] for each walker,
     * which determines its probability of replication or deletion.
     */
    void prepare_branch();

    /**
     * Accumulates statistics during the production phase
     *
     * Collects energy estimators and other observables for later analysis.
     * This includes the mixed estimator ⟨ψᵀ|H|ψDMC⟩/⟨ψᵀ|ψDMC⟩ and the
     * growth estimator based on population dynamics.
     */
    void count();

    /**
     * Saves the current state of the simulation
     * 
     * Moves required data to Results class.
     * Records walker distributions, energy estimators, and other relevant
     * data at specified checkpoints during the simulation.
     */
    void save_progress();

    /**
     * Saves the current state of the simulation
     *
     * Records walker distributions, energy estimators, and other relevant
     * data at specified checkpoints during the simulation.
     */
    DiffusionQuantumResults &get_results();

    /**
     * performs initial diffusion with metropolis sampling.
     *
     * used in case of importance sampling algorithm
     */
    void initial_diffusion();

    /**
     * After initial diffusion samples local_energy for obtaining approximate variational energy
     *
     */
    void sample_variational_energy();

    /**
     * After initial diffusion procedures.
     *
     * run few more iterations with logging local energy to variational_local_energy.
     */
    void finish_initial_diffusion();

private:
    // Vector of electron walkers representing the quantum system
    std::vector<ElectronWalker> electrons;

    // Temporary storage for walkers during diffusion and branching
    std::vector<ElectronWalker> copy_electrons;

    // Branching probabilities for each walker
    std::vector<double> p_values;

    // Diffusion values for walker movement
    std::vector<electron_walker> diffusion_values;

    // Histogram of walker positions for density estimation
    boost::multi_array<int64_t, 2> summed_walkers;

    // Histogram of walker positions for density estimation
    boost::multi_array<int64_t, 2> total_summed_walkers;

    // Effective time used in branching. See Reynolds 1982.
    double eff_d_tau;

    // Variational energy obtained after initial variational diffusion
    // Used in cutoff for local energy. See Runge 1993.
    double variational_energy;

    // Accumulator of variational energy
    double acc_variational_energy;

    // Singleton instance of simulation parameters
    DiffusionQuantumParams *p;

    // Walker population statistics
    int num_alive;    // Current number of active walkers
    int new_alive;    // Number of walkers after branching
    int target_alive; // Target walker population
    int current_it;   // Current iteration counter

    // Energy estimators
    AccumulatedStatistics stats;
    double e_block; // Block-averaged growth estimator

    /**
     * Calculates the average local energy across alive walkers
     *
     * The local energy is defined as E_L(R) = (Hψᵀ(R))/ψᵀ(R) and
     * provides an estimate of the energy when averaged over the
     * walker distribution.
     *
     * See \ref SolverContext::local_energy
     *
     * @return Average local energy across all walkers
     */
    double local_energy_average();

    /**
     * Sets the number of alive walkers and copies walker parameters if needed
     *
     * @param new_alive New number of alive walkers
     * @param wlk Reference walker to copy if population increases
     */
    void set_alive(int new_alive, const ElectronWalker &wlk);

    /**
     * Updates the growth estimator based on population dynamics
     *
     * The growth estimator is based on the rate of change of the
     * walker population and provides an alternative energy estimate.
     */
    void update_growth_estimator();

    /**
     * Performs binning of walker positions for density estimation
     *
     * Creates a histogram of walker positions to estimate the
     * electron density distribution.
     * Params of binning (which exact dimensions should be binned against each other)
     * are defined in \ref DiffusionQuantumParams singleton.
     */
    void binning();

    /**
     * Performs binning of walker positions for density estimation
     *
     * unlike \ref binning sums positions of walkers, so gives estimate
     * of total wavefunction.
     */
    void total_binning();

    // Container for accumulated results
    std::unique_ptr<DiffusionQuantumResults> results;

    // Random number generators for uniform distributions
    boost::random::mt19937 uni_rng;
    boost::random::uniform_real_distribution<double> uniform_generator;

    // Random number generators for Gaussian distributions (walker movement)
    boost::random::mt19937 movement_rng;
    boost::random::normal_distribution<double> movement_generator;

    // Context objects for managing trial wavefunctions and local energies (needed for
    // parallelization)
    std::unique_ptr<SolverContext> general_context;
    std::vector<SolverContext> solver_contexts;

    /**
     * Initializes random number generators
     *
     * Sets up the Mersenne Twister generators for uniform and
     * Gaussian distributions used in the simulation.
     */
    void init_rngs();

    /**
     * Initializes containers for walkers and other data structures
     *
     * Allocates memory and sets initial values for all containers
     * used in the simulation.
     */
    void init_containers();

    /**
     * Initializes the solver context
     *
     * Sets up the context objects that manage trial wavefunctions,
     * local energies, and other quantum mechanical operators.
     * Many contexts are needed because of parallelization
     */
    void init_context();
};

#endif
