#ifndef DIFFUSION_WALKERS_HPP
#define DIFFUSION_WALKERS_HPP

#include <array>
#include <iostream>
#include <vector>

enum class ElectronSpin {UP, DOWN};

/**
 * Basic walker used in simulation
 * This is walker for single electron so it has only 3 dimensions
 * Default values for dimensions are equal to zero if are not used
 */
struct walker {
    ElectronSpin spin;
    std::array<double, 3> cords;

    const double &x() const { return cords[0]; };
    const double &y() const { return cords[1]; };
    const double &z() const { return cords[2]; };

    walker() = default;
    walker(double xval, double yval, double zval) {
        cords[0] = xval;
        cords[1] = yval;
        cords[2] = zval;
    }

    friend std::ostream &operator<<(std::ostream &, const walker &wlk);
};

/**
 * electron_walker representing simple walker in many dimensions systems
 * it is made in order to increase number of dimensions without
 * increasing number of dimensions in walker_struct
 * for 1 electron_walker vector size is 1
 *
 * Total number of dimensions is equal to n_dims x n_electrons
 */

using electron_walker = std::vector<walker>;

/**
 * ElectronWalker is class containing walker coordinates, walker local energy
 * and walker trial_wavef_value
 */
class ElectronWalker {
public:
    electron_walker &get_walker();
    const electron_walker &get_const_walker() const;

    double trial_wavef_value;
    double local_energy;

private:
    electron_walker ele_wlk;
};

std::ostream &operator<<(std::ostream &, const electron_walker &wlk);

class DiffusionWalkers {
public:
    DiffusionWalkers();
    ~DiffusionWalkers();

    double distance(const walker &wlk_a, const walker &wlk_b, int max_dims);
};

#endif
