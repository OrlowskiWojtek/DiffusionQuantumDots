#ifndef ABSTRACT_WF_HPP
#define ABSTRACT_WF_HPP

#include "Core/include/walkers_struct.hpp"
#include "TrialFunctions/include/abstract_orbital.hpp"
#include <functional>

class AbstractSinglebodyOrbital: public AbstractOrbital{
public:
    virtual double operator()(const walker&) = 0; 
    void print() override;
    virtual std::function<double(const walker& wlk)> get_orbital() = 0;
    virtual ~AbstractSinglebodyOrbital() = default;

protected:
    void print_test_to_file() override;
};

#endif
