#include "Core/include/walkers_struct.hpp"
#include "TrialFunctions/include/abstract_singlebody_orbital.hpp"
#include "include/UnitHandler.hpp"
#include <fstream>
#include <iostream>

void AbstractSinglebodyOrbital::print_test_to_file() {
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

void AbstractSinglebodyOrbital::print(){
    std::cout << "Using singlebody orbital" << std::endl;
}
