#ifndef DIFFUSION_QUANTUM_RESULTS_HPP
#define DIFFUSION_QUANTUM_RESULTS_HPP

#include "DiffusionParams/include/params.hpp"
#include <boost/multi_array.hpp>
#include <cstdint>
#include <vector>

class DiffusionQuantumResults {
public:
    DiffusionQuantumResults();
    ~DiffusionQuantumResults();

    void add_energies(double E, double g_est);
    void add_histogram(double time,
                       int time_step,
                       double mean_energy,
                       double mean_growth_estimator,
                       const boost::multi_array<int64_t, 2> &hist);
    void save_to_file();

    boost::multi_array<double, 2> &get_last_psi();
    const std::vector<double> &get_energies();

private:
    DiffusionQuantumParams *p;

    struct HistData {
        double time;
        int time_step;
        double energy;
        double growth_estimator;
        boost::multi_array<double, 2> psi;

        HistData(double time, int time_step, double ene, double gest, boost::multi_array<double, 2> psi)
            : time(time)
            , time_step(time_step)
            , energy(ene)
            , growth_estimator(gest)
            , psi(psi) {};
    };

    std::vector<double> x;
    std::vector<HistData> time_evolution;

    std::vector<double> calibration_energies;
    std::vector<double> calibration_growth;

    void init_x(double x_min, double x_max, int n);
};

#endif
