#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_hamiltonian.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>

#include <bfpt_common/hamiltonian_12.hpp>
#include <bfpt_common/generic_dynamic_unique_kstate_hamiltonian.hpp>
#include <bfpt_common/do_common_recipie.hpp>

#include <armadillo>

#include <iostream>

#include <cassert>

// #######################################################################
// ## main...                                                           ##
// #######################################################################

double bfpt_gs(const size_t n_sites, const unsigned max_pt_order) {
    //using SiteStateT = model_monostar::MonostarSiteState;
    using KstateT = model_monostar::DynamicMonostarUniqueKstate;
    using BasisT = model_monostar::DynamicMonostarUniqueKstateBasis;
    bfpt_common::CommonRecipePrintFlags print_flags;
    //print_flags.print_populated_basis_flag = true; //TEMP
    //print_flags.print_unpopulated_basis_flag = true; //TEMP
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(model_monostar::classical_gs_kstate(n_sites)));
    const auto hamiltonian_12 = model_monostar::prepare_hamiltonian_12(1, 1);
    const bfpt_common::GenericKstateHamiltonian<KstateT> hamiltonian{n_sites, hamiltonian_12};
    return bfpt_common::do_common_recipe(hamiltonian, hamiltonian, basis,
                                         max_pt_order, 0, print_flags);
}

double bfpt_kn_es(const size_t n_sites, const unsigned max_pt_order, const unsigned k_n) {
    //using SiteStateT = model_monostar::MonostarSiteState;
    using KstateT = model_monostar::DynamicMonostarUniqueKstate;
    using BasisT = model_monostar::DynamicMonostarUniqueKstateBasis;
    bfpt_common::CommonRecipePrintFlags print_flags;
    //print_flags.print_populated_basis_flag = true; //TEMP
    //print_flags.print_unpopulated_basis_flag = true; //TEMP
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(model_monostar::classical_es_kstate(n_sites)));
    const auto hamiltonian_12 = model_monostar::prepare_hamiltonian_12(1, 1);
    const bfpt_common::GenericKstateHamiltonian<KstateT> hamiltonian{n_sites, hamiltonian_12};
    return bfpt_common::do_common_recipe(hamiltonian, hamiltonian, basis,
                                         max_pt_order, k_n, print_flags);
}

double bfpt_goldston(const size_t n_sites, const unsigned max_pt_order) {
    return bfpt_kn_es(n_sites, max_pt_order, 0);
}

int main() {
    const unsigned max_pt_order = 6;
    const size_t n_sites = 20;
    std::cout << "------------------------------------------" << std::endl;
    const double gs_energy = bfpt_gs(n_sites, max_pt_order);
    std::cout << "------------------------------------------" << std::endl;
    const double es_gold_energy = bfpt_goldston(n_sites, max_pt_order);
    std::cout << "------------------------------------------" << std::endl;
    std::vector<double> es_energies;
    for (unsigned k_n = 0; k_n < n_sites; k_n++) {
        double es_energy = bfpt_kn_es(n_sites, max_pt_order, k_n);
        es_energies.push_back(es_energy);
        std::cout << "------------------------------------------" << std::endl;
    }
    std::cout << "------------------------------------------" << std::endl;
    std::cout << " ├max_pt_order = " << max_pt_order << std::endl;
    std::cout << " ├state: gs " << max_pt_order << std::endl;
    std::cout << " ││enery = " << gs_energy << std::endl;
    std::cout << " ├state: es [goldston]" << es_gold_energy << std::endl;
    std::cout << " ││enery = " << es_gold_energy << std::endl;
    std::cout << " ││excitation enery = " << es_gold_energy - gs_energy << std::endl;
    for (unsigned k_n = 0; k_n < n_sites; k_n++) {
        const auto es_energy = es_energies[k_n];
        std::cout << " ├state: es [k_n = " << k_n << "]" << std::endl;
        std::cout << " ││enery = " << es_energy << std::endl;
        std::cout << " ││excitation enery = " << es_energy - gs_energy << std::endl;
        std::cout << " ││excitation enery (theory, infinite chain) = " << arma::datum::pi/2 * std::abs(std::sin(2 * arma::datum::pi * k_n / n_sites)) << std::endl;
    }
    std::cout << "------------------------------------------" << std::endl;
    //model_monostar::do_hardcoded_example_analyse();
    return 0;
}
