#ifndef DIFFUSION_QUANTUM_RESULTS_HPP
#define DIFFUSION_QUANTUM_RESULTS_HPP

#include "DiffusionParams/include/params.hpp"
#include <boost/multi_array.hpp>
#include <cstdint>
#include <fstream>
#include <vector>

/**
 * Struct AccumulatedStatistics
 *
 * provides storage for basic results of calculations
 */
struct AccumulatedStatistics {
    // number of iterations of accumulation
    int it;

    // Estimator based on local energy / potential of walkers
    double mixed_estimator;

    // Based on population dynamics
    double growth_estimator;

    // Accumulated mixed estimator
    double acc_mixed_estimator;

    // Accumulated growth estimator
    double acc_growth_estimator;

    // Accumulated growth estimator used for error calcs
    double acc_sq_growth_estimator;

    // Estimate of state energy based on mixed estimator
    double ground_state_mixed_estimator;

    // Estimated of state energy based on growth estimator
    double ground_state_growth_estimator;

    // Standard error of growth estimator
    double growth_estimator_error;

    // calculates averages of estimators
    void finalise();

    // sets values to zero
    void reset();
};

/**
 * class DiffusionQuantumResults
 *
 * provides saving to files and final calculations utilities.
 */
class DiffusionQuantumResults {
public:
    DiffusionQuantumResults();
    ~DiffusionQuantumResults();

    void add_energies(double E, double g_est);
    void save_energies(double time, int num_alive, AccumulatedStatistics &stats);

    void add_histogram(double time,
                       int time_step,
                       AccumulatedStatistics &stats,
                       const boost::multi_array<int64_t, 2> &hist,
                       const boost::multi_array<int64_t, 2> &total_hist);
    void save_to_file();

    boost::multi_array<double, 2> &get_last_psi();
    boost::multi_array<double, 2> &get_last_total_psi();
    const std::vector<double> &get_calib_mixed_energies();
    const std::vector<double> &get_calib_growth_energies();

private:
    DiffusionQuantumParams *p;

    struct HistData {
        double time;
        int time_step;
        AccumulatedStatistics stats;
        boost::multi_array<double, 2> psi;
        boost::multi_array<double, 2> total_psi;

        HistData(double time,
                 int time_step,
                 AccumulatedStatistics stats,
                 boost::multi_array<double, 2> psi,
                 boost::multi_array<double, 2> total_psi)
            : time(time)
            , time_step(time_step)
            , stats(stats)
            , psi(psi)
            , total_psi(total_psi){};
    };

    std::vector<double> x;
    std::vector<HistData> time_evolution;

    std::vector<double> calibration_energies;
    std::vector<double> calibration_growth;

    std::ofstream averaged_energies_file;

    
    void init_x(double x_min, double x_max, int n);

    boost::multi_array<double, 2>
    normalize_pure_diffusion(const boost::multi_array<int64_t, 2> &hist);
    boost::multi_array<double, 2> normalize_trial_wavef(const boost::multi_array<int64_t, 2> &hist);
};

#endif
