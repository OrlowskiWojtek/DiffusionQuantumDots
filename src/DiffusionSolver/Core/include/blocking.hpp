#ifndef ENERGY_BLOCKING_HPP
#define ENERGY_BLOCKING_HPP

#include <vector>

class EnergyBlockingAnalyzer {
public:
    EnergyBlockingAnalyzer();
    ~EnergyBlockingAnalyzer();

    void blocking_analysis(const std::vector<double> &mixed_energies,
                           const std::vector<double> &growth_energies);

private:
    std::vector<int> block_sizes;
    std::vector<double> mixed_stderrs;
    std::vector<double> growth_stderrs;

    std::vector<double> reblock(const std::vector<double> &energies);

    void save_to_file();
};

#endif
