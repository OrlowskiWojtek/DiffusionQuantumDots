#include "include/results.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

DiffusionQuantumResults::DiffusionQuantumResults():
    p(DiffusionQuantumParams::getInstance()){
    
    init_x(p->xmin, p->xmax, p->n_bins);
}

DiffusionQuantumResults::~DiffusionQuantumResults() {}

void DiffusionQuantumResults::init_x(double xmin, double xmax, int n) {
    x.resize(n);
    for (int i = 0; i < n; i++) {
        x[i] =
            xmin + (xmax - xmin) / static_cast<double>(n) *
                       (i + 0.5);
    }
}

void DiffusionQuantumResults::add_histogram(double time,
                                            int time_step,
                                            double mean_energy,
                                            double mean_growth_estimator,
                                            const boost::multi_array<int64_t, 3> &hist) {

    boost::multi_array<double, 3> psi(boost::extents[x.size()][x.size()][x.size()]);
    double dx = (x[1] - x[0]);

    double int_val =
        std::accumulate(hist.data(),
                        hist.data() + hist.num_elements(),
                        0.,
                        [](double acc, const int64_t &val) { return acc + std::pow(val, 2); });
    
    int_val = std::sqrt(int_val * std::pow(dx, p->n_dims));

    std::transform(
        hist.data(), hist.data() + hist.num_elements(), psi.data(), [int_val](int64_t val) {
            return static_cast<double>(val) / int_val;
        });

    time_evolution.emplace_back(time, time_step, mean_energy, mean_growth_estimator, psi);

    std::cout << "|--------------------------------------------------------|\n";
    std::cout << "Adding histogram data" << std::endl;
    std::cout << "time = " << time << "\n";
    std::cout << "time step = " << time_step << "\n";
    std::cout << "energy = " << mean_energy << "\n";
    std::cout << "growth estimator = " << mean_growth_estimator << "\n";
    std::cout << "|--------------------------------------------------------|\n";
}

void DiffusionQuantumResults::save_to_file() {
    std::ofstream file("evolution.dqmc.dat");

    for (HistData &data : time_evolution) {
        file << "#\t" << "time:" << data.time << "\n";
        file << "#\t" << "time_step:" << data.time_step << "\n";
        file << "#\t" << "energy:" << data.energy << "\n";
        file << "#\t" << "xmin:" << x.front() << "\n";
        file << "#\t" << "xmax:" << x.back() << "\n";
        file << "#\t" << "nbins:" << x.size() << "\n";

        // TODO -> swtich to dimension version
        for (size_t i = 0; i < x.size(); i++) {
            for (size_t j = 0; j < x.size(); j++) {
                int j_idx = p->n_dims > 1 ? j : x.size() / 2; 
                file << time_evolution.back().psi[i][j_idx][x.size() / 2] << "\t";
            }
            file << "\n";
        }
    }

    file.close();

    std::ofstream file_energies("energies.dqmc.dat");
    file_energies << "#\t<V>\tg_est\n";

    if (calibration_energies.size() != calibration_growth.size()) {
        file_energies.close();
        return;
    }

    for (size_t i = 0; i < calibration_energies.size(); i++) {
        file_energies << calibration_energies[i] << "\t" << calibration_growth[i] << "\n";
    }

    file_energies.close();
}

void DiffusionQuantumResults::add_energies(double ene, double growth_estimator) {
    calibration_energies.push_back(ene);
    calibration_growth.push_back(growth_estimator);
}

const std::vector<double> &DiffusionQuantumResults::get_energies() { return calibration_energies; }
