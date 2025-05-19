#include "Core/include/solver_context.hpp"
#include "TrialFunctions/include/jastrow_slater.hpp"

SolverContext::SolverContext()
    : walkers_helper(std::make_unique<DiffusionWalkers>())
    , p(DiffusionQuantumParams::getInstance()) {

    JastrowSlaterOrbitalParams orbital_params;
    orbital_params.electron_number = p->n_electrons;
    orbital_params.omegas = p->omegas;
    orbital_params.effective_mass = p->effective_mass;
    orbital_params.dims = p->n_dims;
    orbital_params.a = 0.005;
    orbital_params.b = 2.;

    orbital = std::make_unique<JastrowSlaterOrbital>(orbital_params);
}

double SolverContext::p_value(const electron_walker &wlk, const electron_walker &prev_wlk, double growth_estimator) {
    return std::exp(-p->d_tau * ((local_energy(wlk) + local_energy(prev_wlk)) / 2. - growth_estimator));
}

double SolverContext::local_energy(const electron_walker &wlk) {
    double dr = 1e-6;
    double dr2 = 1e-12;
    double kinetic_term = 0;

    double cent_value = trial_wavef(wlk);

    m_back_walker_buffer = wlk;
    m_front_walker_buffer = wlk;
    for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
        for (int d = 0; d < p->n_dims; d++) {

            m_back_walker_buffer[wlk_idx].cords[d] -= dr;
            m_front_walker_buffer[wlk_idx].cords[d] += dr;
            double back_value = trial_wavef(m_back_walker_buffer);
            double forw_value = trial_wavef(m_front_walker_buffer);

            kinetic_term += (back_value - 2 * cent_value + forw_value);

            m_back_walker_buffer[wlk_idx].cords[d] = wlk[wlk_idx].cords[d];
            m_front_walker_buffer[wlk_idx].cords[d] = wlk[wlk_idx].cords[d];
        }
    }

    kinetic_term = -0.5 * kinetic_term / (cent_value * dr2);

    double potential_term = std::accumulate(wlk.begin(), wlk.end(), 0., [this](double acc, const walker &single_wlk) {
        return acc + p->pot(single_wlk);
    });

    // TODO: hardcoded for now for coulomb interaction -> possibly change
    double interaction = 0;
    for (int i = 0; i < p->n_electrons; i++) {
        for (int j = i + 1; j < p->n_electrons; j++) {
            interaction += 1 / (p->epsilon * walkers_helper->distance(wlk[i], wlk[j], p->n_dims));
        }
    }

    return kinetic_term + interaction + potential_term;
}

double SolverContext::trial_wavef(const electron_walker &wlk) { return (*this->orbital)(wlk); }
