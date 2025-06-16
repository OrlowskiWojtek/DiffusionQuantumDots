#include "TrialFunctions/include/abstract_manybody_orbital.hpp"
#include "include/UnitHandler.hpp"
#include <fstream>
#include <iostream>

void AbstractManybodyOrbital::print_test_to_file() {
    double xmin = UnitHandler::length(UnitHandler::TO_AU, -50.);
    double xmax = UnitHandler::length(UnitHandler::TO_AU, 50.);
    double dx = (xmax - xmin) / 100.;

    std::ofstream file("TrialWavefunctionTest");

    electron_walker test_walker;
    test_walker.resize(2);
    test_walker[0].cords[2] = 0;

    for (double i = xmin; i < xmax; i += dx) {
        for (double j = xmin; j < xmax; j += dx) {
            test_walker[0].cords[0] = i;
            test_walker[1].cords[0] = j;
            file << (*this)(test_walker) << "\t";
        }
        file << "\n";
    }
    file << std::endl;
    file.close();
}

void AbstractManybodyOrbital::print(){
    std::cout << "Using manybody orbital" << std::endl;
}
