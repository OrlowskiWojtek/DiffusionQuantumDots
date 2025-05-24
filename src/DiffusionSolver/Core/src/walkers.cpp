#include "Core/include/walkers.hpp"
#include <algorithm>
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

std::ostream& operator<<(std::ostream& os, const walker& wlk){
    std::for_each(wlk.cords.begin(), wlk.cords.end(), [&os](double cord){
        os << cord << " ";
    });

    return os;
}

std::ostream& operator<<(std::ostream& os, const electron_walker& ele_wlk){
    std::for_each(ele_wlk.begin(), ele_wlk.end(), [&os](const walker& wlk){
        os << "|" << wlk << "|"; 
    });

    return os;
}

electron_walker& ElectronWalker::get_walker(){
    return ele_wlk;
}

const electron_walker& ElectronWalker::get_const_walker() const{
    return ele_wlk;
}
