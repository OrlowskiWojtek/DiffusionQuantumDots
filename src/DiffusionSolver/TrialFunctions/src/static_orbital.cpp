#include "TrialFunctions/include/static_orbital.hpp"
#include "include/UnitHandler.hpp"
#include <functional>
#include <iostream>
#include <fstream>

StaticOrbital::StaticOrbital() {}

void StaticOrbital::print() { std::cout << "Using static orbital f(x) = 1"; }

std::function<double(const walker &wlk)> StaticOrbital::get_orbital() {
    return std::function<double(const walker &wlk)>([](const walker &) {
        return 1;
    });
}

double StaticOrbital::operator()(const walker& wlk){
    return 1;
}

