#ifndef DIFFUSION_WALKERS_HPP
#define DIFFUSION_WALKERS_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/multi_array.hpp>

#include "DiffusionParams/include/params.hpp"
#include "include/results.hpp"
#include "include/walkers_struct.hpp"

#include <boost/random/uniform_real_distribution.hpp>
#include <memory>
#include <vector>

class DiffusionWalkers {
public:
    DiffusionWalkers();
    ~DiffusionWalkers();

    void init_walkers();
    void diffuse();
    void branch();
    void eval_p();

    void count();
    void save_progress();

    DiffusionQuantumResults &get_results();

private:
    DiffusionQuantumParams* p;

    std::vector<walker> walkers;
    std::vector<walker> copy_walkers;
    std::vector<double> p_values;

    boost::multi_array<int64_t, 3> hist;

    boost::random::mt19937 rng;
    boost::random::normal_distribution<double> movement_generator;

    boost::random::mt19937 uni_rng;
    boost::random::uniform_real_distribution<double> uniform_generator =
        boost::random::uniform_real_distribution<double>(0, 1);

    double xmin, xmax;
    double growth_estimator;
    double ground_state_estimator;

    double ground_mean;
    double ground_mean2;
    int blocks_passed;

    double current_Et;
    double Eblock;
    bool calibrating;

    int current_it;
    int accumulation_it;

    double d_tau;
    std::function<double(walker)> V;
    std::function<double(walker)> trial_wavef;

    size_t num_alive;
    size_t new_alive;
    size_t target_alive;

    size_t n_bins;
    size_t block_size;

    int dims;

    std::unique_ptr<DiffusionQuantumResults> results;

    void set_alive(int new_alive, const walker& wlk);
    void update_growth_estimator();
    bool apply_nodes(int walker_idx);

    void binning();
};

#endif
