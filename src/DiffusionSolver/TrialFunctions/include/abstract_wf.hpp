#ifndef ABSTRACT_WF_HPP
#define ABSTRACT_WF_HPP

#include "Core/include/walkers_struct.hpp"
#include <functional>

class AbstractOrbital{
public:
    virtual double operator()(const walker&) = 0; 
    virtual void print() = 0;
    virtual std::function<double(const walker& wlk)> get_orbital() = 0;
    virtual ~AbstractOrbital() = default;

private:
    virtual void print_test_to_file() = 0;
};

#endif
