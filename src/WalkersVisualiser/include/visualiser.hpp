#ifndef WALKERS_VISUALIZER_HPP
#define WALKERS_VISUALIZER_HPP

#include <boost/multi_array.hpp>

/*!
*   Class for visualisation of basic data after or during simulation.
*   It uses one graphical library and performs data plotting.
*
*   It should handle:
*   - plotting energies from results (mixed estimator & growth estimator)
*   - plotting energies average
*   - plotting walkers distribution in selected dimensions
*   - plotting blocking analysis
*/
class WalkersVisualiser{


public:
    void make_surf_plot(boost::multi_array<double, 2>& psi,
                        boost::multi_array<double, 2>& total_psi,
                        int n);

};

#endif
