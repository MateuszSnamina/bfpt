#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_hamiltonian_k0.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>

#include <extensions/adaptors.hpp>

// #include <armadillo>
#include <cassert>  // TODO:  move to hamiltonian sourcefile
#include <iostream>
#include <iterator>

// #######################################################################
// ## DynamicMonostarHamiltonian                                        ##
// #######################################################################

namespace model_monostar {

class DynamicMonostarHamiltonian {
   public:
    DynamicMonostarHamiltonian(const size_t n_sites);
    void push_back_conjugated_states_to_basis(
        const DynamicMonostarUniqueKstate& generator,
        DynamicMonostarUniqueKstateBasis& basis) const;
    // void fill_k0_hamiltonian_matrix_coll(
    //   size_t n_col,
    //   arma::mat& k0_hamiltonian_matrix) const;
   private:
    const size_t _n_sites;
};

inline DynamicMonostarHamiltonian::DynamicMonostarHamiltonian(const size_t n_sites)
    : _n_sites(n_sites) {}

inline void DynamicMonostarHamiltonian::push_back_conjugated_states_to_basis(
    const DynamicMonostarUniqueKstate& generator,
    DynamicMonostarUniqueKstateBasis& basis) const {
    assert(generator.n_sites() == _n_sites);
    const auto generator_range = generator.to_range();
    for (size_t i = 0, j = 1; i < _n_sites; i++, j = (i + 1) % _n_sites) {
        if (*std::next(std::begin(generator_range), i) == gs &&
            *std::next(std::begin(generator_range), j) == gs) {
            const auto conjugated_range = generator_range | extension::boost::adaptors::refined(i, es) | extension::boost::adaptors::refined(j, es);
            const auto conjugated_kstate_ptr = std::make_shared<DynamicMonostarUniqueKstate>(conjugated_range, kstate::ctr_from_range);
            basis.add_element(conjugated_kstate_ptr);
        }
    }
    for (size_t i = 0, j = 1; i < _n_sites; i++, j = (i + 1) % _n_sites) {
        if (*std::next(std::begin(generator_range), i) == es &&
            *std::next(std::begin(generator_range), j) == es) {
            const auto conjugated_range = generator_range | extension::boost::adaptors::refined(i, gs) | extension::boost::adaptors::refined(j, gs);
            const auto conjugated_kstate_ptr = std::make_shared<DynamicMonostarUniqueKstate>(conjugated_range, kstate::ctr_from_range);
            basis.add_element(conjugated_kstate_ptr);
        }
    }
}

}  // namespace model_monostar

// #######################################################################
// ## main...                                                           ##
// #######################################################################

void populate_pt_basis(
    const model_monostar::DynamicMonostarHamiltonian& hamiltonian,
    const unsigned max_pt_order,
    model_monostar::DynamicMonostarUniqueKstateBasis& basis) {
    // assert();
    unsigned last_chunk_size = basis.size();
    for (unsigned pt_order = 0; pt_order < max_pt_order; pt_order++) {
        const unsigned old_size = basis.size();
        assert(last_chunk_size <= old_size);
        for (unsigned idx = old_size - last_chunk_size; idx < old_size; idx++) {
            const auto el = basis.vec_index()[idx];
            assert(el);
            hamiltonian.push_back_conjugated_states_to_basis(*el, basis);
        }
        const unsigned new_size = basis.size();
        last_chunk_size = new_size - old_size;
    }
}

void bfpt_gs(const size_t n_sites, const unsigned max_pt_order) {
    // Define hamiltonian and basis:
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    model_monostar::DynamicMonostarHamiltonian hamiltonian{n_sites};
    // Define pt_orders subspace basis:
    std::vector<model_monostar::MonostarSiteState> generator_array(n_sites, model_monostar::gs);
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(generator_array, kstate::ctr_from_range));
    // ----
    std::cout << basis;
    // ----
    // Generate higher pt_orders subspace basis:
    populate_pt_basis(hamiltonian, max_pt_order, basis);
    // ----
    std::cout << basis;
    // ----
}

void bfpt_k0_es(const size_t n_sites, const unsigned max_pt_order) {
    // Define hamiltonian and basis:
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    model_monostar::DynamicMonostarHamiltonian hamiltonian{n_sites};
    // Define pt_orders subspace basis:
    std::vector<model_monostar::MonostarSiteState> generator_array(n_sites, model_monostar::gs);
    generator_array[0] = model_monostar::es;
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(generator_array, kstate::ctr_from_range));
    // ----
    std::cout << basis;
    // ----
    // Generate higher pt_orders subspace basis:
    populate_pt_basis(hamiltonian, max_pt_order, basis);
    // ----
    std::cout << basis;
    // ----
}

int main() {
    const unsigned max_pt_order = 8;
    const size_t n_sites = 20;
    //bfpt_gs(n_sites, max_pt_order);
    bfpt_k0_es(n_sites, max_pt_order);
    return 0;
}
