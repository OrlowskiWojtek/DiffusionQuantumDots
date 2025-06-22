#include "Core/include/blocking.hpp"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numeric>

EnergyBlockingAnalyzer::EnergyBlockingAnalyzer() {}
EnergyBlockingAnalyzer::~EnergyBlockingAnalyzer() {}

void EnergyBlockingAnalyzer::blocking_analysis(const std::vector<double> &mixed_energies,
                                               const std::vector<double> &growth_energies) {

    if (mixed_energies.size() != growth_energies.size()) {
        return;
    }

    size_t n = mixed_energies.size();

    block_sizes.clear();

    size_t block_size = 1;
    do {
        block_sizes.push_back(block_size);
        block_size *= 2;
    } while (block_size < n);

    mixed_stderrs = reblock(mixed_energies);
    growth_stderrs = reblock(growth_energies);

    save_to_file();
}

std::vector<double> EnergyBlockingAnalyzer::reblock(const std::vector<double> &energies) {
    std::vector<double> temp_energies{energies.begin(), energies.end()};
    std::vector<double> standard_reblocked_errors;

    for (int bs : block_sizes) {
        int num_blocks = energies.size() / bs;

        std::vector<double> reblocked(num_blocks, 0);
        for (size_t i = 0; i < reblocked.size(); i++) {
            reblocked[i] = std::accumulate(energies.begin() + bs * i,
                                           energies.begin() + bs * (i + 1),
                                           0.) /
                           static_cast<double>(bs);
        }

        double mean = std::accumulate(reblocked.begin(), reblocked.end(), 0.) /
                      static_cast<double>(reblocked.size());
        double mean2 = std::accumulate(reblocked.begin(),
                                       reblocked.end(),
                                       0.,
                                       [](double acc, double val) {
                                           return acc + std::pow(val, 2);
                                       }) /
                       static_cast<double>(reblocked.size());

        if (reblocked.size() - 1 != 0) {
            standard_reblocked_errors.push_back(
                std::sqrt((mean2 - std::pow(mean, 2)) / static_cast<double>(reblocked.size() - 1)));
        }
    }

    return standard_reblocked_errors;
}

void EnergyBlockingAnalyzer::save_to_file() {
    std::ofstream file("reblock_analysis.dqmc.dat");
    file << "#\tReblocking std errors\n";
    file << "#\tblock_size\tmixed_std_error\tgrowth_std_error\n";

    for (size_t i = 0; i < block_sizes.size(); i++) {
        file << block_sizes[i] << "\t" << mixed_stderrs[i] << "\t" << growth_stderrs[i] << "\n";
    }

    file.close();
}
