#ifndef DIFFUSION_QUANTUM_RESULTS_HPP
#define DIFFUSION_QUANTUM_RESULTS_HPP

#include <cstddef>
#include <vector>

class DiffusionQuantumResults{
public:
    DiffusionQuantumResults();
    ~DiffusionQuantumResults();

    size_t num_alive;

    void add_histogram(const std::vector<double>& hist);
private:  
};

#endif
