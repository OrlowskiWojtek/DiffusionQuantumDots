#ifndef DIFFUSION_QUANTUM_RESULTS_HPP
#define DIFFUSION_QUANTUM_RESULTS_HPP

#include <cstdint>
#include <vector>

class DiffusionQuantumResults{
public:
    DiffusionQuantumResults();
    ~DiffusionQuantumResults();

    void add_energies(double E, double g_est);
    void init_x(double x_min, double x_max, int n);
    void add_histogram(double time, int time_step, double energy, const std::vector<int64_t>& hist);
    void save_to_file();

    const std::vector<double>& get_energies();
private:  
    
    struct HistData{
        double time;
        int time_step;
        double energy;
        std::vector<double> psi;

        HistData(double time, int time_step, double ene, std::vector<double> psi):time(time), time_step(time_step), energy(ene), psi(psi){};
    };

    std::vector<double> x;
    std::vector<HistData> time_evolution;

    std::vector<double> calibration_energies;
    std::vector<double> calibration_growth;
};

#endif
