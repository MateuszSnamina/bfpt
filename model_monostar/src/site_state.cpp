
#include <model_monostar/monostar_hamiltonian_k0.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>

//#include <extensions/range_streamer.hpp>  // to delete?!
#include <extensions/adaptors.hpp>
#include <kstate/basis.hpp>
#include <kstate/basis_streamer.hpp>
#include <kstate/kstate_concrete.hpp>
#include <kstate/kstate_streamer.hpp>

#include <boost/operators.hpp>
#include <iostream>

// #######################################################################
// ## MonostarSiteState                                                 ##
// #######################################################################

namespace model_monostar {

struct MonostarSiteState : boost::totally_ordered<MonostarSiteState> {
    MonostarSiteState(bool is_excited) : _is_excited(is_excited){};
    bool _is_excited;
    bool operator<(const MonostarSiteState& other) const {
        return static_cast<int>(this->_is_excited) < static_cast<int>(other._is_excited);
    }
    bool operator==(const MonostarSiteState& other) const {
        return static_cast<int>(this->_is_excited) == static_cast<int>(other._is_excited);
    }
};  // namespace model_monostar

const MonostarSiteState gs{false};
const MonostarSiteState es{true};

std::ostream&
operator<<(std::ostream& stream, const MonostarSiteState& state) {
    stream << (state._is_excited ? '*' : '_');
    return stream;
}

}  // namespace model_monostar

// #######################################################################
// ## monostar_kstate_range_streamer_settings                           ##
// #######################################################################

namespace model_monostar {

const auto monostar_kstate_range_streamer_settings =
    extension::boost::RangeStreamerSettings()
        .set_stream_preparer([](::std::ostream& s) { s << "𝕊𝕥𝕒𝕣𝕂𝕤𝕥𝕒𝕥𝕖⦃"; })
        .set_stream_sustainer([](::std::ostream& s, size_t i) {})
        .set_stream_separer([](::std::ostream& s) {})
        .set_stream_finisher([](::std::ostream& s) { s << "⦄"; })
        .set_format_independence_flag(true);

// #######################################################################
// ## DynamicMonostarKstate nad DynamicMonostarUniqueKstate             ##
// #######################################################################

using DynamicMonostarKstate = kstate::DynamicKstate<MonostarSiteState>;
using DynamicMonostarUniqueKstate = kstate::DynamicUniqueKstate<MonostarSiteState>;

std::ostream&
operator<<(std::ostream& stream, const DynamicMonostarKstate& state) {
    auto kstate_streamer = kstate::KstateStreamer(stream).set_range_streamer_settings(monostar_kstate_range_streamer_settings);
    kstate_streamer.stream(state);
    return stream;
}

std::ostream&
operator<<(std::ostream& stream, const DynamicMonostarUniqueKstate& state) {
    auto kstate_streamer = kstate::KstateStreamer(stream).set_range_streamer_settings(monostar_kstate_range_streamer_settings);
    kstate_streamer.stream(state);
    return stream;
}

}  // namespace model_monostar

// #######################################################################
// ## DynamicMonostarUniqueKstateBasis                                  ##
// #######################################################################

namespace model_monostar {

using DynamicMonostarUniqueKstateBasis = kstate::Basis<DynamicMonostarUniqueKstate>;

std::ostream&
operator<<(std::ostream& stream, const DynamicMonostarUniqueKstateBasis& state) {
    auto basis_streamer = kstate::BasisStreamer(stream).set_range_streamer_settings_for_kstate(monostar_kstate_range_streamer_settings);
    basis_streamer.stream(state);
    return stream;
}

}  // namespace model_monostar

// #######################################################################
// ## DynamicMonostarHamiltonian                                        ##
// #######################################################################

namespace model_monostar {

class DynamicMonostarHamiltonian {
   public:
    void push_back_conjugated_states_to_basis(
        const DynamicMonostarUniqueKstate& generator,
        DynamicMonostarUniqueKstateBasis& basis) const;
    // void fill_k0_hamiltonian_matrix_coll(
    // size_t n_col,
    // arma::mat& k0_hamiltonian_matrix) const;
};

inline void DynamicMonostarHamiltonian::push_back_conjugated_states_to_basis(
    const DynamicMonostarUniqueKstate& generator,
    DynamicMonostarUniqueKstateBasis& basis) const {
    const auto xxx = generator.to_range();
    const auto yyy = xxx | extension::boost::adaptors::refined(0, es);
    const auto zzz = std::make_shared<DynamicMonostarUniqueKstate>(yyy, kstate::ctr_from_range);
    basis.add_element(zzz);
}

}  // namespace model_monostar

// #######################################################################
// ## main...                                                           ##
// #######################################################################

int main() {
    model_monostar::DynamicMonostarUniqueKstateBasis basis(6);
    /// ----
    std::array<model_monostar::MonostarSiteState, 6>
        generator_array = {model_monostar::gs, model_monostar::gs, model_monostar::gs, model_monostar::gs, model_monostar::gs, model_monostar::gs};
    // DynamicMonostarUniqueKstate generator(generator_array, kstate::ctr_from_range);
    /// ----
    std::cout << basis;
    /// ----
    basis.add_element(std::make_shared<model_monostar::DynamicMonostarUniqueKstate>(generator_array, kstate::ctr_from_range));
    /// ----
    std::cout << basis;
    /// ----
    model_monostar::DynamicMonostarHamiltonian hamiltonian;
    hamiltonian.push_back_conjugated_states_to_basis(*basis.vec_index()[0], basis);
    /// ----
    std::cout << basis;
    /// ----
    return 0;
}
