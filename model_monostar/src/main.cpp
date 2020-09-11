// #include <model_monostar/hardcoded_example.hpp>
#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_hamiltonian.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>

#include <bfpt_common/do_common_recipie.hpp>

#include <armadillo>

#include <iostream>
#include <iterator>

#include <cassert>
#include <complex>

using namespace std::complex_literals;

// #######################################################################
// ## main...                                                           ##
// #######################################################################

double bfpt_gs(const size_t n_sites, const unsigned max_pt_order) {
    bfpt_common::CommonRecipePrintFlags print_flags;
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(model_monostar::classical_gs_kstate(n_sites)));
    const model_monostar::DynamicMonostarHamiltonian hamiltonian{n_sites};
    return bfpt_common::do_common_recipe(hamiltonian, hamiltonian, basis,
                                         max_pt_order, 0, print_flags);
}

double bfpt_kn_es(const size_t n_sites, const unsigned max_pt_order, const unsigned k_n) {
    bfpt_common::CommonRecipePrintFlags print_flags;
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(model_monostar::classical_es_kstate(n_sites)));
    const model_monostar::DynamicMonostarHamiltonian hamiltonian{n_sites};
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
    }
    std::cout << "------------------------------------------" << std::endl;
    //model_monostar::do_hardcoded_example_analyse();
    return 0;
}
