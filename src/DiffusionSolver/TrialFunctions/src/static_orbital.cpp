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


void StaticOrbital::print_test_to_file() {
    double xmin = UnitHandler::length(UnitHandler::TO_AU, -20.);
    double xmax = UnitHandler::length(UnitHandler::TO_AU, 20.);
    double dx = (xmax - xmin) / 100.;

    std::ofstream file("TrialWavefunctionTest");

    walker test_walker;
    test_walker.cords[2] = 0;

    for (double i = xmin; i < xmax; i += dx) {
        for (double j = xmin; j < xmax; j += dx) {
            test_walker.cords[0] = i;
            test_walker.cords[1] = j;
            file << (*this)(test_walker) << "\t";
        }
        file << "\n";
    }
    file << std::endl;
    file.close();
}
