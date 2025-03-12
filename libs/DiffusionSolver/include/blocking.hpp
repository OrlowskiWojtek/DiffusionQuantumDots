#ifndef ENERGY_BLOCKING_HPP
#define ENERGY_BLOCKING_HPP

#include <vector>

class EnergyBlockingAnalyzer{
public:
    EnergyBlockingAnalyzer();
    ~EnergyBlockingAnalyzer();

    void blocking_analysis(const std::vector<double>& energies);

private:
    std::vector<int> block_sizes;
    std::vector<double> reblocked;

    void save_to_file();
};

#endif
