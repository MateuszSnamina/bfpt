#ifndef MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP
#define MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP

#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_kstate.hpp>

#include <bfpt_common/i_dynamic_unique_kstate_hamiltonian.hpp>
#include <bfpt_common/i_dynamic_unique_kstate_populator.hpp>

#include <armadillo>

#include <type_traits>
#include <map>
#include <tuple>

// #######################################################################
// ## DynamicMonostarHamiltonian                                        ##
// #######################################################################

namespace model_monostar {

template<typename SiteState>
struct SiteStatePair {
    static_assert(!std::is_reference_v<SiteState>, "SiteState must not be a reference type.");
    static_assert(!std::is_const_v<SiteState>, "SiteState must not be a cosnt type.");
    static_assert(!std::is_volatile_v<SiteState>, "SiteState must not be a volatile type.");
    SiteState state_1;
    SiteState state_2;
};

template<typename SiteState>
bool operator<(const SiteStatePair<SiteState>& lhs, const SiteStatePair<SiteState>& rhs) {
    return std::make_tuple(lhs.state_1, lhs.state_2) < std::make_tuple(rhs.state_1, rhs.state_2);
}

template<typename SiteState>
struct CoupleInfo {
    double coef;
    SiteStatePair<SiteState> state12;
};

class DynamicMonostarHamiltonian : public bfpt_common::IDynamicUniqueKstatePopulator<MonostarSiteState>,
        public bfpt_common::IDynamicUniqueKstateHamiltonian<MonostarSiteState> {
public:
    using SiteState = MonostarSiteState;
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
    /*const*/ std::map<SiteStatePair<SiteState>, double> _diag_info; //TODO restore const
    /*const*/ std::multimap<SiteStatePair<SiteState>, CoupleInfo<SiteState>> _full_off_diag_info; //TODO restore const
    /*const*/ std::multimap<SiteStatePair<SiteState>, CoupleInfo<SiteState>> _half_off_diag_info; //TODO restore const
};

}  // namespace model_monostar

#endif
