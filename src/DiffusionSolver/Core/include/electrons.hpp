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

class DiffusionQuantumElectrons {
public:
    DiffusionQuantumElectrons();

    void diffuse();
    void accept_movement();

    void branch();
    void eval_p();

    void count();
    void save_progress();

    DiffusionQuantumResults &get_results();

private:
    std::vector<electron_walker> electrons;
    std::vector<electron_walker> copy_electrons;
    std::vector<double> p_values;

    std::vector<electron_walker> diffusion_values;
    boost::multi_array<int64_t, 3> summed_walkers;

    DiffusionQuantumParams *p;

    int num_alive;
    int new_alive;
    int target_alive;
    int current_it;
    int accu_it;

    double mixed_estimator; // average of local energy estimator
    double growth_estimator;
    double e_block;
    double ground_state_estimator;

    double local_energy_average();

    void set_alive(int new_alive, const electron_walker &wlk);
    void update_growth_estimator();
    void binning();

    std::unique_ptr<DiffusionQuantumResults> results;

    boost::random::mt19937 uni_rng;
    boost::random::uniform_real_distribution<double> uniform_generator;

    boost::random::mt19937 movement_rng;
    boost::random::normal_distribution<double> movement_generator;

    std::unique_ptr<SolverContext> general_context;
    std::vector<SolverContext> solver_contexts;

    void init_rngs();
    void init_containers();
    void init_context();
};

#endif
