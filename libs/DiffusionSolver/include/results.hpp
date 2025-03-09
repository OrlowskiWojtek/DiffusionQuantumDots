#ifndef DIFFUSION_QUANTUM_RESULTS_HPP
#define DIFFUSION_QUANTUM_RESULTS_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

class DiffusionQuantumResults{
public:
    DiffusionQuantumResults();
    ~DiffusionQuantumResults();

    size_t num_alive;

    void add_energy(double E);
    void init_x(double x_min, double x_max, int n);
    void add_histogram(double time, int time_step, const std::vector<int64_t>& hist);
private:  
    
    struct HistData{
        double time;
        int time_step;
        std::vector<double> psi;

        HistData(double time, int time_step, std::vector<double> psi):time(time), time_step(time_step), psi(psi){};
    };

    std::vector<double> x;
    std::vector<HistData> time_evolution;

};

#endif
