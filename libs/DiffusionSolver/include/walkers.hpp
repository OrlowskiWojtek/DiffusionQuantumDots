#ifndef DIFFUSION_WALKERS_HPP
#define DIFFUSION_WALKERS_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

#include "params.hpp"
#include "results.hpp"
#include <boost/random/uniform_real_distribution.hpp>
#include <memory>
#include <vector>

class DiffusionWalkers {
public:
    DiffusionWalkers();
    ~DiffusionWalkers();

    void init_walkers(const DiffusionQuantumParams &params);
    void diffuse();
    void branch();
    void eval_p();

    void count();
    void save_progress();

    DiffusionQuantumResults& get_results();

private:
    struct walker { // TODO -> add another dimensions
        double x;
    };

    std::vector<walker> walkers;
    std::vector<walker> copy_walkers;
    std::vector<double> p_values;
    std::vector<int64_t> hist;

    boost::random::mt19937 rng;
    boost::random::normal_distribution<double> movement_generator;

    boost::random::mt19937 uni_rng;
    boost::random::uniform_real_distribution<double> uniform_generator =
        boost::random::uniform_real_distribution<double>(0, 1);

    double xmin, xmax;
    double growth_estimator; // TODO switch to results class
    double ground_state_estimator; 
    double Et;
    double current_Et;
    bool calibrating;

    int current_it;
    int accumulation_it;

    double d_tau;
    std::function<double(double)> V;

    size_t num_alive;
    size_t new_alive;
    size_t target_alive;

    size_t n_bins;

    std::unique_ptr<DiffusionQuantumResults> results;

    void set_alive(int new_alive, double position);
    void update_growth_estimator();
};

#endif
