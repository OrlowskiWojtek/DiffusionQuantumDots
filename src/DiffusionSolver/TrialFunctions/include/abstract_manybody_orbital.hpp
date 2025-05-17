#ifndef ABSTRACT_MANYBODY_WF_HPP
#define ABSTRACT_MANYBODY_WF_HPP

#include "Core/include/walkers_struct.hpp"
#include "TrialFunctions/include/abstract_orbital.hpp"
#include <functional>

class AbstractManybodyOrbital : public AbstractOrbital{
public:
    virtual double operator()(const electron_walker&) = 0; 
    void print() override;
    virtual std::function<double(const electron_walker& wlk)> get_orbital() = 0;
    virtual ~AbstractManybodyOrbital() = default;

protected:
    void print_test_to_file() override;
};

#endif
