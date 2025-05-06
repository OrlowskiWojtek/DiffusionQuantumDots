#ifndef DIFFUSION_QUANTUM_ELECTRONS_HPP
#define DIFFUSION_QUANTUM_ELECTRONS_HPP

#include <Core/include/walkers_struct.hpp>
#include "DiffusionParams/include/params.hpp"
#include <vector>

class DiffusionQuantumElectrons{
public:
    DiffusionQuantumElectrons();

private:
    // electron_walker representing simple walker in many dimensions systems
    // it is made in order to increase number of dimensions without 
    // increasing number of dimensions in walker_struct
    // for 1 electron it vectors size is 1
    using electron_walker = std::vector<walker>;
    std::vector<electron_walker> electrons;

    DiffusionQuantumParams* p;
};

#endif
