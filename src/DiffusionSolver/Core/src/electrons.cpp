#include "Core/include/electrons.hpp"
#include <execution>
#include <numeric>

DiffusionQuantumElectrons::DiffusionQuantumElectrons()
    : p(DiffusionQuantumParams::getInstance())
    , walkers_helper(std::make_unique<DiffusionWalkers>()) {

    electrons.resize(p->n_electrons);
    for (auto &ele : electrons) {
        ele.resize(p->nmax_walkers);
    }

    p_values.resize(p->nmax_walkers);

    current_it = 0;
    accu_it = 0;

    growth_estimator = 0;
    e_block = 0;

    init_rngs();
}

void DiffusionQuantumElectrons::diffuse() {
    current_it++;

    std::for_each(electrons.begin(), electrons.begin() + num_alive, [&](electron_walker &ele) {
        std::for_each(ele.begin(), ele.end(), [&](walker &wlk) {
            walkers_helper->apply_drift(wlk);
            walkers_helper->apply_diffusion(wlk);
        });
    });
}

void DiffusionQuantumElectrons::eval_p() {
    std::transform(std::execution::par, // smarter threading is needed
                   electrons.begin(),
                   electrons.begin() + num_alive,
                   copy_electrons.begin(),
                   p_values.begin(),
                   [&](electron_walker &wlk, electron_walker &prev_wlk) {
                       if (apply_nodes(wlk, prev_wlk)) {
                           return 0.;
                       }

                       return p_value(wlk, prev_wlk);
                   });
}

void DiffusionQuantumElectrons::branch() {
    new_alive = 0;

    for (int i = 0; i < num_alive; i++) {
        int m = static_cast<int>(p_values[i] + uniform_generator(uni_rng));
        set_alive(m, electrons[i]);
    }

    num_alive = new_alive;
    std::copy(electrons.begin(), copy_electrons.begin() + num_alive, electrons.begin());

    update_growth_estimator();
}

void DiffusionQuantumElectrons::set_alive(int N, const electron_walker &wlk) {
    for (int i = new_alive; i < new_alive + N; i++) {
        if (i >= static_cast<int>(electrons.size())) {
            std::cout << "Error: The programme ran out of allocated memory. "
                         "Required resize."
                      << std::endl;
            break;
        }
        copy_electrons[i] = wlk;
    }

    new_alive += N;
}

// TODO: thing about it
void DiffusionQuantumElectrons::update_growth_estimator(){
    if (accu_it == 0) {
        int nblock = current_it % p->n_block;
        e_block = (growth_estimator + e_block * nblock) / (nblock + 1);
    } else {
        e_block = (growth_estimator + e_block * accu_it) / (accu_it + 1);
    }

    growth_estimator =
        e_block - 1. / p->d_tau * std::log(static_cast<double>(num_alive) / static_cast<double>(target_alive));
}

double DiffusionQuantumElectrons::trial_wavef(const electron_walker &wlk) {
    return std::accumulate(wlk.begin(), wlk.end(), 0., [&](double sum, const walker &single_walker) {
        return sum + walkers_helper->trial_wf_value(single_walker);
    });
}

bool DiffusionQuantumElectrons::apply_nodes(const electron_walker &wlk, const electron_walker &prev_wlk) {
    return (trial_wavef(wlk) * trial_wavef(prev_wlk)) <= 0;
}

double DiffusionQuantumElectrons::p_value(const electron_walker & wlk, const electron_walker & prev_wlk) {
    return std::exp(-p->d_tau * ((local_energy(wlk) + local_energy(prev_wlk)) / 2. - growth_estimator));
}

double DiffusionQuantumElectrons::local_energy(const electron_walker& wlk){
    double local_ene = std::accumulate(wlk.begin(), wlk.end(), 0., [&](double sum, const walker& wlk){
        return sum + walkers_helper->local_energy(wlk);
    });

    // TODO: hardcoded for now for coulomb interaction -> possibly change
    double interaction = 0;
    for(int i = 0; i < p->n_electrons; i++){
        for(int j = i + 1; j < p->n_electrons; j++){
            interaction += 1 / (p->epsilon * walkers_helper->distance(wlk[i], wlk[j]));
        }
    }

    return local_ene + interaction;
}

void DiffusionQuantumElectrons::init_rngs(){
    uniform_generator = boost::random::uniform_real_distribution<double>(0, 1);
}
