#ifndef DIFFUSION_QUANTUM_RESULTS_HPP
#define DIFFUSION_QUANTUM_RESULTS_HPP

#include <cstdint>
#include <vector>
#include "walkers_struct.hpp"

class DiffusionQuantumResults {
public:
    DiffusionQuantumResults();
    ~DiffusionQuantumResults();

    void add_energies(double E, double g_est);
    void init_x(double x_min, double x_max, int n);
    void add_histogram(double time,
                       int time_step,
                       double mean_energy,
                       double mean_growth_estimator,
                       const std::vector<std::vector<int64_t>> &hist);
    void save_to_file();

    void set_dims(int ndims);
    const std::vector<double> &get_energies();

private:
    int dims;

    struct HistData {
        double time;
        int time_step;
        double energy;
        double growth_estimator;
        std::vector<walker> psi;

        HistData(double time, int time_step, double ene, double gest, std::vector<walker> psi)
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
};

#endif
