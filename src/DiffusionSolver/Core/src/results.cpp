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
    averaged_energies_file << "# time[a.u.]\t"
                           << "population\t"
                           << "avg_mixed_stimator[meV]\t"
                           << "avg_growth_estimator[meV]\t"
                           << "mixed_estimator[meV]\t"
                           << "growth_estimator[meV]\t"
                           << "growth_estimator_error_mean[meV]\n";
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
                                            AccumulatedStatistics &stats,
                                            const boost::multi_array<int64_t, 2> &hist) {

#ifndef PURE_DIFFUSION
    boost::multi_array<double, 2> psi = normalize_trial_wavef(hist);
#else
    boost::multi_array<double, 2> psi = normalize_pure_diffusion(hist);
#endif

    stats.finalise();
    time_evolution.emplace_back(time, time_step, stats, psi);

    std::cout << "|--------------------------------------------------------|\n";
    std::cout << "Adding histogram data" << std::endl;
    std::cout << "time = " << time << "\n";
    std::cout << "time step = " << time_step << "\n";
    std::cout << "mixed estimator [meV] = "
              << UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, stats.mixed_estimator)
              << "\n";
    std::cout << "mean mixed estimator [meV] = "
              << UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT,
                                     stats.ground_state_mixed_estimator)
              << "\n";
    std::cout << "growth estimator [meV] = "
              << UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT, stats.growth_estimator)
              << "\n";
    std::cout << "mean growth estimator [meV] = "
              << UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT,
                                     stats.ground_state_growth_estimator)
              << "\n";
    std::cout << "growth estimator error of mean [meV] ="
              << UnitHandler::energy(UnitHandler::ConvMode::TO_DEFAULT,
                                     stats.growth_estimator_error)
              << "\n";
    std::cout << "|--------------------------------------------------------|\n";
}

// TODO: fragmentize
void DiffusionQuantumResults::save_to_file() {
    std::ofstream file("evolution.dqmc.dat");

    for (HistData &data : time_evolution) {
        file << "#\t" << "time:" << data.time << "\n";
        file << "#\t" << "time_step:" << data.time_step << "\n";
        file << "#\t" << "energy [meV]:" << data.stats.mixed_estimator << "\n";
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
                                            AccumulatedStatistics &stats) {
    stats.finalise();

    averaged_energies_file
        << time << "\t" << num_alive << "\t"
        << UnitHandler::energy(UnitHandler::TO_DEFAULT, stats.ground_state_mixed_estimator) << "\t"
        << UnitHandler::energy(UnitHandler::TO_DEFAULT, stats.ground_state_growth_estimator) << "\t"
        << UnitHandler::energy(UnitHandler::TO_DEFAULT, stats.mixed_estimator) << "\t"
        << UnitHandler::energy(UnitHandler::TO_DEFAULT, stats.growth_estimator) << "\t"
        << UnitHandler::energy(UnitHandler::TO_DEFAULT, stats.growth_estimator_error) << "\n";
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

void AccumulatedStatistics::reset() {
    it = 0;

    mixed_estimator = 0;
    growth_estimator = 0;

    acc_growth_estimator = 0;
    acc_sq_growth_estimator = 0;

    acc_mixed_estimator = 0;
}

void AccumulatedStatistics::finalise() {
    if (it == 0) {
        return;
    }

    ground_state_growth_estimator = acc_growth_estimator / static_cast<double>(it);
    ground_state_mixed_estimator = acc_mixed_estimator / static_cast<double>(it);

    double growth_sq_mean = acc_sq_growth_estimator / static_cast<double>(it);
    growth_estimator_error =
        std::sqrt((growth_sq_mean - ground_state_growth_estimator * ground_state_growth_estimator) /
                  static_cast<double>(it));
}
