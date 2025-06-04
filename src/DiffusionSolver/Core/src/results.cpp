#include "Core/include/results.hpp"
#include "include/UnitHandler.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

DiffusionQuantumResults::DiffusionQuantumResults()
    : p(DiffusionQuantumParams::getInstance()) {

    averaged_mixed_estimator.reserve(p->total_time_steps);
    averaged_growth_estimator.reserve(p->total_time_steps);
    init_x(p->xmin, p->xmax, p->n_bins);
}

DiffusionQuantumResults::~DiffusionQuantumResults() {}

void DiffusionQuantumResults::init_x(double xmin, double xmax, int n) {
    x.resize(n);
    for (int i = 0; i < n; i++) {
        x[i] = xmin + (xmax - xmin) / static_cast<double>(n) * (i + 0.5);
    }
}

void DiffusionQuantumResults::add_histogram(double time,
                                            int time_step,
                                            double mean_energy,
                                            double mean_growth_estimator,
                                            const boost::multi_array<int64_t, 2> &hist) {

    boost::multi_array<double, 2> psi(boost::extents[x.size()][x.size()]);
    double dx = (x[1] - x[0]);

    double int_val = std::accumulate(
        hist.data(), hist.data() + hist.num_elements(), 0., [](double acc, const int64_t &val) {
            return acc + val;
        });

    int_val = int_val * std::pow(dx, p->n_dims);

    std::transform(
        hist.data(), hist.data() + hist.num_elements(), psi.data(), [int_val](int64_t val) {
            return static_cast<double>(val) / int_val;
        });

    time_evolution.emplace_back(time, time_step, mean_energy, mean_growth_estimator, psi);

    std::cout << "|--------------------------------------------------------|\n";
    std::cout << "Adding histogram data" << std::endl;
    std::cout << "time = " << time << "\n";
    std::cout << "time step = " << time_step << "\n";
    std::cout << "energy [meV] = "
              << UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, mean_energy) << "\n";
    std::cout << "growth estimator = "
              << UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, mean_growth_estimator)
              << "\n";
    std::cout << "|--------------------------------------------------------|\n";
}

// TODO: fragmentize
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
                file << data.psi[i][j] << "\t";
            }
            file << "\n";
        }
        file << std::endl;
    }

    file.close();

    std::ofstream file_energies("energies.dqmc.dat");
    file_energies << "#\tmixed_stimator[meV]\tgrowth_estimator[meV]\n";

    if (calibration_energies.size() != calibration_growth.size()) {
        file_energies.close();
        return;
    }

    for (size_t i = 0; i < calibration_energies.size(); i++) {
        file_energies << calibration_energies[i] << "\t" << calibration_growth[i] << "\n";
    }

    file_energies.close();

    std::ofstream file_avg_energies("averaged_energies.dqmc.dat");
    file_avg_energies << "#\tmixed_stimator[meV]\tgrowth_estimator[meV]\n";

    if (averaged_mixed_estimator.size() != averaged_growth_estimator.size()) {
        file_avg_energies.close();
        return;
    }

    for (size_t i = 0; i < averaged_mixed_estimator.size(); i++) {
        file_avg_energies << averaged_mixed_estimator[i] << "\t" << averaged_growth_estimator[i] << "\n";
    }

    file_avg_energies.close();
}

void DiffusionQuantumResults::add_energies(double mixed_estimator, double growth_estimator) {
    calibration_energies.push_back(
        UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, mixed_estimator));
    calibration_growth.push_back(
        UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, growth_estimator));
}

void DiffusionQuantumResults::add_avg_energies(double mixed_estimator, double growth_estimator) {
    averaged_mixed_estimator.push_back(
        UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, mixed_estimator));
    averaged_growth_estimator.push_back(
        UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, growth_estimator));
}

const std::vector<double> &DiffusionQuantumResults::get_energies() { return calibration_energies; }

// TODO: make checks
boost::multi_array<double, 2> &DiffusionQuantumResults::get_last_psi() {
    return time_evolution.back().psi;
};
