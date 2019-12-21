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

void bfpt_gs(const size_t n_sites) {
    model_monostar::DynamicMonostarUniqueKstateBasis basis{n_sites};
    model_monostar::DynamicMonostarHamiltonian hamiltonian{n_sites};
    /// ----
    std::vector<model_monostar::MonostarSiteState> generator_array(n_sites, model_monostar::gs);
    /// ----
    std::cout << basis;
    /// ----
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(generator_array, kstate::ctr_from_range));
    /// ----
    std::cout << basis;
    /// ----
    hamiltonian.push_back_conjugated_states_to_basis(*basis.vec_index()[0], basis);
    /// ----
    std::cout << basis;
    /// ----
    hamiltonian.push_back_conjugated_states_to_basis(*basis.vec_index()[1], basis);
    /// ----
    std::cout << basis;
    /// ----
}

int main() {
    const size_t n_sites = 6;
    bfpt_gs(n_sites);
    return 0;
}
