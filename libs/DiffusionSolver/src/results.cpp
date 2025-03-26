#include "results.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

DiffusionQuantumResults::DiffusionQuantumResults() {}

DiffusionQuantumResults::~DiffusionQuantumResults() {}

void DiffusionQuantumResults::init_x(double xmin, double xmax, int n) {
    x.resize(n);
    for (int i = 0; i < n; i++) {
        x[i] =
            xmin + (xmax - xmin) / static_cast<double>(n) *
                       i; // TODO - put exact points (centers of blocks, so move half unit to right)
    }
}

void DiffusionQuantumResults::add_histogram(double time,
                                            int time_step,
                                            double mean_energy,
                                            double mean_growth_estimator,
                                            const std::vector<std::vector<int64_t>> &hist) {
    std::vector<walker> psi(hist.size());
    double dx = (x[1] - x[0]);

    double int_val =
        std::accumulate(hist.begin(),
                        hist.begin() + dims,
                        0.,
                        [](double acc, const std::vector<int64_t> &dim_val) {
                            double in_dim = std::accumulate(
                                dim_val.begin(), dim_val.end(), 0., [](double acc, double value) {
                                    return acc + std::pow(value, 2);
                                });
                            return acc + in_dim;
                        }) *
        std::pow(dx, dims);

    std::transform(hist.begin(), hist.end(), psi.begin(), [int_val](std::vector<int64_t> hist) {
        return walker(hist[0] / sqrt(int_val), hist[1] / sqrt(int_val), hist[2] / sqrt(int_val));
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

        for (size_t i = 0; i < x.size(); i++) {
            file << x[i] << "\t";
        }
        file << "\n";

        for (int dim; dim < dims; dim++) {
            
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

void DiffusionQuantumResults::set_dims(int ndims) { dims = ndims; }
