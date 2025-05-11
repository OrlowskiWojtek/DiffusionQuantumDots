#ifndef DIFFUSION_QUANTUM_ELECTRONS_HPP
#define DIFFUSION_QUANTUM_ELECTRONS_HPP

#include <Core/include/walkers_struct.hpp>
#include "Core/include/walkers.hpp"
#include "DiffusionParams/include/params.hpp"
#include <vector>

class DiffusionQuantumElectrons{
public:
    DiffusionQuantumElectrons();

    void diffuse();
    void branch();
    void eval_p();

    void count();
    void save_progress();

    DiffusionQuantumResults &get_results();
private:
    // electron_walker representing simple walker in many dimensions systems
    // it is made in order to increase number of dimensions without 
    // increasing number of dimensions in walker_struct
    // for 1 electron_walker vector size is 1
    using electron_walker = std::vector<walker>;

    std::vector<electron_walker> electrons;
    std::vector<electron_walker> copy_electrons;
    std::vector<double> p_values;
    
    boost::multi_array<int64_t, 3> summed_walkers;

    DiffusionQuantumParams* p;
    std::unique_ptr<DiffusionWalkers> walkers_helper;

    int num_alive;
    int new_alive;  
    int target_alive;
    int current_it;
    int accu_it;

    double mixed_estimator; // average of local energy estimator 
    double growth_estimator;
    double e_block;
    double ground_state_estimator;

    bool apply_nodes(const electron_walker&, const electron_walker&);
    void apply_drift(electron_walker& wlk);
    double p_value(const electron_walker &, const electron_walker &);
    double trial_wavef(const electron_walker&);
    double local_energy(const electron_walker &);
    double local_energy_average();

    void set_alive(int new_alive, const electron_walker &wlk);
    void update_growth_estimator();
    void binning();

    std::unique_ptr<DiffusionQuantumResults> results;

    boost::random::mt19937 uni_rng;
    boost::random::uniform_real_distribution<double> uniform_generator;

    void init_rngs();
    void init_containers();

    std::array<double, 3> drift(const electron_walker &wlk);
};

#endif
