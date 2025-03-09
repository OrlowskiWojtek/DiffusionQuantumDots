#include "results.hpp"
#include <algorithm>
#include <cmath>
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
                                            const std::vector<int64_t> &hist) {
    std::vector<double> psi(x.size());
    double dx = (x[1] - x[0]);

    double int_val =
        std::accumulate(hist.begin(),
                        hist.end(),
                        0.,
                        [](double acc, double val) { return acc + std::pow(val, 2); }) *
        dx;

    std::transform(
        hist.begin(), hist.end(), psi.begin(), [int_val](double x) { return x / sqrt(int_val); });

    time_evolution.emplace_back(time, time_step, psi);
    std::cout << "---------------" << std::endl;
    for (size_t i = 0; i < x.size(); i++) {
        std::cout << x[i] << " " << psi[i] << "\n";
    }
}
