#ifndef STATIC_ORBITAL_HPP
#define STATIC_ORBITAL_HPP

#include "TrialFunctions/include/abstract_wf.hpp"
#include <functional>

class StaticOrbital : public AbstractOrbital {
public:
    StaticOrbital();

    double operator()(const walker &wlk) override;
    void print() override;

    std::function<double(const walker& wlk)> get_orbital() override;

private:
    void print_test_to_file() override;
};

#endif
