#ifndef MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP
#define MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP

#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_kstate.hpp>

#include <bfpt_common/i_dynamic_unique_kstate_hamiltonian.hpp>
#include <bfpt_common/i_dynamic_unique_kstate_populator.hpp>

#include <armadillo>

// #######################################################################
// ## DynamicMonostarHamiltonian                                        ##
// #######################################################################

namespace model_monostar {

class DynamicMonostarHamiltonian : public bfpt_common::IDynamicUniqueKstatePopulator<MonostarSiteState>,
                                   public bfpt_common::IDynamicUniqueKstateHamiltonian<MonostarSiteState> {
   public:
    DynamicMonostarHamiltonian(const size_t n_sites);
    // Generates all conjugated states:
    void push_back_coupled_states_to_basis(
        const DynamicMonostarUniqueKstate& generator,
        DynamicMonostarUniqueKstateBasis& basis) const override;
    // Generates hamiltonian matrix:
    arma::sp_cx_mat make_kn_hamiltonian_matrix(
        const DynamicMonostarUniqueKstateBasis& basis,
        const unsigned k_n) const override;

   private:
    void fill_kn_hamiltonian_matrix_coll(
        const DynamicMonostarUniqueKstateBasis& basis,
        size_t n_col,
        arma::sp_cx_mat& kn_hamiltonian_matrix,
        const unsigned k_n) const;
    void fill_kn_hamiltonian_matrix(
        const DynamicMonostarUniqueKstateBasis& basis,
        arma::sp_cx_mat& kn_hamiltonian_matrix,
        const unsigned k_n) const;

   private:
    const size_t _n_sites;
};

}  // namespace model_monostar

#endif