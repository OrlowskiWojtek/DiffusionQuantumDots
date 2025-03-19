#include "blocking.hpp"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numeric>

EnergyBlockingAnalyzer::EnergyBlockingAnalyzer() {}
EnergyBlockingAnalyzer::~EnergyBlockingAnalyzer() {}

void EnergyBlockingAnalyzer::blocking_analysis(const std::vector<double> &energies) {
    std::vector<double> temp_energies{energies.begin(), energies.end()};
    block_sizes.clear();
    stderrs.clear();

    size_t block_size = 1;
    do {
        block_sizes.push_back(block_size);
        block_size *= 2;
    } while (block_size < energies.size());

    temp_energies.resize(block_sizes.back());

    for (int bs : block_sizes) {
        std::vector<double> reblocked(temp_energies.size() / bs, 0);
        for (size_t i = 0; i < reblocked.size(); i++) {
            reblocked[i] = std::accumulate(temp_energies.begin() + bs * i,
                                           temp_energies.begin() + bs * (i + 1),
                                           0.) /
                           static_cast<double>(bs);
        }

        double mean = std::accumulate(reblocked.begin(), reblocked.end(), 0.) /
                      static_cast<double>(reblocked.size());
        double mean2 =
            std::accumulate(reblocked.begin(),
                            reblocked.end(),
                            0.,
                            [](double acc, double val) { return acc + std::pow(val, 2); }) /
            static_cast<double>(reblocked.size());

        stderrs.push_back(
            std::sqrt((mean2 - std::pow(mean, 2)) / static_cast<double>(reblocked.size() - 1)));
    }

    save_to_file();
}

void EnergyBlockingAnalyzer::save_to_file() {
    std::ofstream file("reblock_analysis.dqmc.dat");
    file << "#\tReblocking std errors\n";
    file << "#\tblock_size\tstd_error\n";

    for (size_t i = 0; i < block_sizes.size(); i++) {
        file << block_sizes[i] << "\t" << stderrs[i] << "\n";
    }

    file.close();
}
