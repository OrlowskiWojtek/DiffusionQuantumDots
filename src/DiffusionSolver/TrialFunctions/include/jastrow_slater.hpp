#ifndef JASTROW_SLATER_ORBITAL_HPP
#define JASTROW_SLATER_ORBITAL_HPP

#include "TrialFunctions/include/abstract_manybody_orbital.hpp"
#include "TrialFunctions/include/abstract_singlebody_orbital.hpp"
#include <memory>

#include <armadillo>

struct JastrowSlaterOrbitalParams{
    int electron_number;
    std::vector<ElectronSpin> spins;
    std::vector<double> omegas;
    double effective_mass;
    int dims;

    double a;
    double b;
};

class JastrowSlaterOrbital: public AbstractManybodyOrbital{
public:
    JastrowSlaterOrbital();
    JastrowSlaterOrbital(JastrowSlaterOrbitalParams p);

    double operator()(const electron_walker&) override; 
    void print() override;
    std::function<double(const electron_walker& wlk)> get_orbital() override;

    std::vector<std::unique_ptr<AbstractSinglebodyOrbital>> single_body_orbitals; // implementing composition pattern
private:
    std::function<double(const electron_walker& wlk)> orbital;
    JastrowSlaterOrbitalParams p;
    arma::Mat<double> slater_matrix_up;
    arma::Mat<double> slater_matrix_down;

    int spins_up;
    int spins_down;

    // take from walkers.hpp
    double distance(const walker&, const walker&);
    void init_orbital();

    void add_single_orbital(std::unique_ptr<AbstractSinglebodyOrbital>);
};

#endif
