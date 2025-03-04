#ifndef DIFFUSION_WALKERS_HPP
#define DIFFUSION_WALKERS_HPP

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "params.hpp"
#include <boost/random/uniform_real_distribution.hpp>
#include <vector>

class DiffusionWalkers{
public:
    DiffusionWalkers();
    ~DiffusionWalkers();

    void init_walkers(const DiffusionQuantumParams& params);
    void diffuse();
    void branch();
    void eval_p();

private:
    struct walker{
        double x;
        bool alive;
    };

    std::vector<walker> walkers;
    std::vector<walker> copy_walkers;
    std::vector<double> p_values;
    
    boost::random::mt19937 rng;
    boost::random::normal_distribution<double> movement_generator;

    boost::random::mt19937 uni_rng;
    boost::random::uniform_real_distribution<double> uniform_generator = boost::random::uniform_real_distribution<double>(0,1);

    double growth_estimator; // TODO switch to results class
    double Et;
    int current_it;
    
    double d_tau;
    std::function<double(double)> V;

    size_t num_alive;
    size_t new_alive;
    size_t target_alive;

    void set_alive(int new_alive, double position);
    void update_growth_estimator();
    void generate_hist(int n_bins);
};

#endif
