#include "blocking.hpp"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>

EnergyBlockingAnalyzer::EnergyBlockingAnalyzer() {}
EnergyBlockingAnalyzer::~EnergyBlockingAnalyzer() {}

void EnergyBlockingAnalyzer::blocking_analysis(const std::vector<double> &energies) {
    std::vector<double> temp_energies{energies.begin(), energies.end()};
    block_sizes.clear();
    reblocked.clear();

    size_t i = 0;
    int block_size = 1;
    size_t numMeas = energies.size();

    while (numMeas >= 100) {
        double mean = 0;
        double mean2 = 0;
        for (size_t i = 0; i < numMeas / 2; i++) {
            mean += temp_energies[2 * i] + temp_energies[2 * i + 1];
            mean2 += std::pow(temp_energies[2 * i], 2) + std::pow(temp_energies[2 * i + 1], 2);
            temp_energies[i] = (temp_energies[2 * i] + temp_energies[2 * i + 1]) / 2.;
        }

        if (2 * i < numMeas) {
            mean += temp_energies[2 * i];
            mean2 += std::pow(temp_energies[2 * i], 2);
        }
        mean /= static_cast<double>(numMeas);
        mean2 /= static_cast<double>(numMeas);

        block_sizes.push_back(block_size);
        reblocked.push_back((mean2 - std::pow(mean, 2)) / (numMeas - 1));

        numMeas /= 2;
        block_size *= 2;
    }

    save_to_file();
}

void EnergyBlockingAnalyzer::save_to_file(){
    std::ofstream file("reblock_analysis.dqmc.dat");
    file << "#\tReblocking std errors\n";
    file << "#\tblock_size\tstd_error\n";

    for(size_t i = 0; i < block_sizes.size(); i++){
        file << block_sizes[i] << "\t" << reblocked[i] << "\n";
    }
    
    file.close();
}
