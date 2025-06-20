#include "Core/include/results.hpp"
#include "include/UnitHandler.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

DiffusionQuantumResults::DiffusionQuantumResults()
    : p(DiffusionQuantumParams::getInstance()) {

    init_x(p->xmin, p->xmax, p->n_bins);
    averaged_energies_file.open("averaged_energies.dqmc.dat");
    averaged_energies_file << "# ntarget:" << p->n0_walkers << "\n";
    averaged_energies_file
        << "#time[a.u.]\tpopulation\tavg_mixed_stimator[meV]\tavg_growth_estimator[meV]"
           "\tmixed_estimator[meV]\tgrowth_estimator[meV]\n";
}

DiffusionQuantumResults::~DiffusionQuantumResults() {
    if (averaged_energies_file.is_open()) {
        averaged_energies_file.close();
    }
}

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

#ifndef PURE_DIFFUSION
    boost::multi_array<double, 2> psi = normalize_trial_wavef(hist);
#else
    boost::multi_array<double, 2> psi = normalize_pure_diffusion(hist);
#endif

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
        file << "#\t" << "xmin:" << UnitHandler::length(UnitHandler::TO_DEFAULT, x.front()) << "\n";
        file << "#\t" << "xmax:" << UnitHandler::length(UnitHandler::TO_DEFAULT, x.back()) << "\n";
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
}

void DiffusionQuantumResults::add_energies(double mixed_estimator, double growth_estimator) {
    calibration_energies.push_back(
        UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, mixed_estimator));
    calibration_growth.push_back(
        UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, growth_estimator));
}

void DiffusionQuantumResults::save_energies(double time,
                                            int num_alive,
                                            double avg_mixed_estimator,
                                            double avg_growth_estimator,
                                            double mixed_estimator,
                                            double growth_estimator) {
    averaged_energies_file << time << "\t" << num_alive << "\t"
                           << UnitHandler::energy(UnitHandler::TO_DEFAULT, avg_mixed_estimator)
                           << "\t"
                           << UnitHandler::energy(UnitHandler::TO_DEFAULT, avg_growth_estimator)
                           << "\t" << UnitHandler::energy(UnitHandler::TO_DEFAULT, mixed_estimator)
                           << "\t" << UnitHandler::energy(UnitHandler::TO_DEFAULT, growth_estimator)
                           << "\n";
}

const std::vector<double> &DiffusionQuantumResults::get_calib_mixed_energies() {
    return calibration_energies;
}

const std::vector<double> &DiffusionQuantumResults::get_calib_growth_energies() {
    return calibration_growth;
}

// TODO: make checks
boost::multi_array<double, 2> &DiffusionQuantumResults::get_last_psi() {
    return time_evolution.back().psi;
};

boost::multi_array<double, 2>
DiffusionQuantumResults::normalize_trial_wavef(const boost::multi_array<int64_t, 2> &hist) {
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

    return psi;
}

boost::multi_array<double, 2>
DiffusionQuantumResults::normalize_pure_diffusion(const boost::multi_array<int64_t, 2> &hist) {
    boost::multi_array<double, 2> psi(boost::extents[x.size()][x.size()]);
    double dx = (x[1] - x[0]);

    double int_val = std::accumulate(
        hist.data(), hist.data() + hist.num_elements(), 0., [](double acc, const int64_t &val) {
            return acc + val * val;
        });

    int_val = int_val * std::pow(dx, p->n_dims);

    std::transform(
        hist.data(), hist.data() + hist.num_elements(), psi.data(), [int_val](int64_t val) {
            return static_cast<double>(val) / std::sqrt(int_val);
        });

    return psi;
}
