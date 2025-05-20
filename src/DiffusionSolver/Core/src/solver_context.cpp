#include "Core/include/solver_context.hpp"
#include "TrialFunctions/include/jastrow_slater.hpp"
#include <cmath>

boost::random::mt19937 SolverContext::s_seed_generator(time(0));

SolverContext::SolverContext()
    : walkers_helper(std::make_unique<DiffusionWalkers>())
    , p(DiffusionQuantumParams::getInstance()) {

    drift_velocity.resize(p->n_electrons);
    green_diffusion_norm = std::pow(2 * M_PI * p->d_tau, -p->n_dims * p->n_electrons / 2.);

    init_potential();
    init_orbital();
    init_rng();
}

void SolverContext::init_potential() { potential = std::make_unique<HarmonicPotentialFunctor>(p->pot_params); }

void SolverContext::init_orbital() {
    JastrowSlaterOrbitalParams orbital_params;
    orbital_params.electron_number = p->n_electrons;
    orbital_params.omegas = p->omegas;
    orbital_params.effective_mass = p->effective_mass;
    orbital_params.dims = p->n_dims;
    orbital_params.a = 0.005;
    orbital_params.b = 2.;

    orbital = std::make_unique<JastrowSlaterOrbital>(orbital_params);
}

void SolverContext::init_rng() {
    uni_rng.seed(s_seed_generator());
    uniform_generator = boost::random::uniform_real_distribution<double>(0, 1);

    movement_rng.seed(s_seed_generator());
    movement_generator = boost::random::normal_distribution<double>(0, std::sqrt(p->d_tau));
}

double SolverContext::trial_wavef(const electron_walker &wlk) { return (*this->orbital)(wlk); }

void SolverContext::calc_trial_wavef(ElectronWalker &wlk) { wlk.trial_wavef_value = trial_wavef(wlk.get_walker()); }

double SolverContext::p_value(ElectronWalker &wlk, ElectronWalker &prev_wlk, double growth_estimator) {
    return std::exp(-p->d_tau * ((local_energy(wlk) + local_energy(prev_wlk)) / 2. - growth_estimator));
}

double SolverContext::local_energy(ElectronWalker &wlk) {
    double dr = 1e-5;
    double dr2 = 1e-10;
    double kinetic_term = 0;

    m_back_walker_buffer = wlk.get_const_walker();
    m_front_walker_buffer = wlk.get_const_walker();
    double cent_trial_wavef = wlk.trial_wavef_value;

    for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
        for (int d = 0; d < p->n_dims; d++) {

            m_back_walker_buffer[wlk_idx].cords[d] -= dr;
            m_front_walker_buffer[wlk_idx].cords[d] += dr;
            double back_value = trial_wavef(m_back_walker_buffer);
            double forw_value = trial_wavef(m_front_walker_buffer);

            kinetic_term += (back_value - 2 * cent_trial_wavef + forw_value);

            m_back_walker_buffer[wlk_idx].cords[d] = wlk.get_const_walker()[wlk_idx].cords[d];
            m_front_walker_buffer[wlk_idx].cords[d] = wlk.get_const_walker()[wlk_idx].cords[d];
        }
    }

    kinetic_term = -0.5 * kinetic_term / (cent_trial_wavef * dr2);

    double potential_term = std::accumulate(
        wlk.get_const_walker().begin(), wlk.get_const_walker().end(), 0., [this](double acc, const walker &single_wlk) {
            return acc + (*potential)(single_wlk);
        });

    // TODO: hardcoded for now for coulomb interaction -> possibly change
    double interaction = 0;
    for (int i = 0; i < p->n_electrons; i++) {
        for (int j = i + 1; j < p->n_electrons; j++) {
            interaction += 1 / (p->epsilon * walkers_helper->distance(
                                                 wlk.get_const_walker()[i], wlk.get_const_walker()[j], p->n_dims));
        }
    }

    double local_ene = kinetic_term + interaction + potential_term;
    wlk.local_energy = local_ene;

    return local_ene;
}

bool SolverContext::apply_nodes(ElectronWalker &wlk, const ElectronWalker &prev_wlk) {
    if ((wlk.trial_wavef_value * prev_wlk.trial_wavef_value) <= 0) {
        wlk = prev_wlk;
        return true;
    }

    return false;
}

void SolverContext::move_walkers(ElectronWalker &wlk, electron_walker &diff_value) {
    prepare_drift(wlk);
    apply_drift(wlk);
    apply_diffusion(wlk, diff_value);
}

void SolverContext::prepare_drift(const ElectronWalker &ele_wlk) {
    double dr = 1e-6;
    m_front_walker_buffer = ele_wlk.get_const_walker();
    double cent_trial_wavef = ele_wlk.trial_wavef_value;

    for (int d = 0; d < p->n_dims; d++) {
        for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
            m_front_walker_buffer[wlk_idx].cords[d] += dr;

            double forw_value = trial_wavef(m_front_walker_buffer);
            drift_velocity[wlk_idx].cords[d] = (forw_value - cent_trial_wavef) / (dr * cent_trial_wavef);
            m_front_walker_buffer[wlk_idx].cords[d] = ele_wlk.get_const_walker()[wlk_idx].cords[d];
        }
    }
}

void SolverContext::apply_drift(ElectronWalker &ele_wlk) {
    electron_walker &wlk = ele_wlk.get_walker();

    for (int d = 0; d < p->n_dims; d++) {
        for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
            wlk[wlk_idx].cords[d] += p->d_tau * drift_velocity[wlk_idx].cords[d];
        }
    }
}

void SolverContext::apply_diffusion(ElectronWalker &ele_wlk, electron_walker &diff_wlk) {
    electron_walker &wlk = ele_wlk.get_walker();

    for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
        for (int d = 0; d < p->n_dims; d++) {
            diff_wlk[wlk_idx].cords[d] = movement_generator(movement_rng);
            wlk[wlk_idx].cords[d] += diff_wlk[wlk_idx].cords[d];
        }
    }
}

double SolverContext::green_diffusion_term(const ElectronWalker &curr_wlk, const ElectronWalker &prev_wlk) {
    double movement_propability_nominator = 0;

    prepare_drift(prev_wlk);

    const electron_walker &curr_ele_wlk = curr_wlk.get_const_walker();
    const electron_walker &prev_ele_wlk = prev_wlk.get_const_walker();

    for (int d = 0; d < p->n_dims; d++) {
        for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
            movement_propability_nominator += std::pow(curr_ele_wlk[wlk_idx].cords[d] - prev_ele_wlk[wlk_idx].cords[d] -
                                                           p->d_tau * drift_velocity[wlk_idx].cords[d],
                                                       2);
        }
    }

    return green_diffusion_norm * std::exp(-movement_propability_nominator / (2 * p->d_tau));
}

void SolverContext::check_movement(ElectronWalker &wlk, ElectronWalker &prev_wlk, electron_walker &diff_value) {
    double movement_propability_nominator = 0;

    for (int d = 0; d < p->n_dims; d++) {
        for (int wlk_idx = 0; wlk_idx < p->n_electrons; wlk_idx++) {
            movement_propability_nominator += std::pow(diff_value[wlk_idx].cords[d], 2);
        }
    }

    double g_d_current = green_diffusion_norm * std::exp(-movement_propability_nominator / (2 * p->d_tau));
    double g_d_back = green_diffusion_term(prev_wlk, wlk);

    double trial_current = wlk.trial_wavef_value;
    double trial_previous = prev_wlk.trial_wavef_value;

    double p_acc = std::min(1., g_d_back * std::pow(trial_current, 2) / (g_d_current * std::pow(trial_previous, 2)));

    if (uniform_generator(uni_rng) > p_acc) {
        wlk = prev_wlk;
    }
}
