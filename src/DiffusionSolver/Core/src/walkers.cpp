#include "Core/include/walkers.hpp"
#include <cmath>

DiffusionWalkers::DiffusionWalkers() {}
DiffusionWalkers::~DiffusionWalkers() {}

double DiffusionWalkers::distance(const walker &wlk_a, const walker &wlk_b, int max_dims) {
    double s = 0;
    for (int i = 0; i < max_dims; i++) {
        s += std::pow(wlk_a.cords[i] - wlk_b.cords[i], 2);
    }

    return std::sqrt(s);
}
