#ifndef DIFFUSION_WALKERS_HPP
#define DIFFUSION_WALKERS_HPP

#include <vector>

class DiffusionWalkers{
public:
    DiffusionWalkers();
    ~DiffusionWalkers();

    void random_init();
private:
    struct walker{
        double x;
        bool alive;
    };

    std::vector<walker> walkers;
    std::vector<walker> copy_walkers;
};

#endif
