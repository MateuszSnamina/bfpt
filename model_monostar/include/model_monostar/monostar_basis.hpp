#ifndef MODEL_MONOSTAR_MONOSTAR_BASIS_HPP
#define MODEL_MONOSTAR_MONOSTAR_BASIS_HPP

#include <model_monostar/monostar_kstate.hpp>

#include <kstate/basis.hpp>
#include <kstate/basis_streamer.hpp>

// #######################################################################
// ## DynamicMonostarUniqueKstateBasis                                  ##
// #######################################################################

namespace model_monostar {

using DynamicMonostarUniqueKstateBasis = kstate::Basis<DynamicMonostarUniqueKstate>;

}  // namespace model_monostar

// #######################################################################
// ## DynamicMonostarUniqueKstateBasis - printing                       ##
// #######################################################################

namespace model_monostar {

inline std::ostream&
operator<<(std::ostream& stream, const DynamicMonostarUniqueKstateBasis& state) {
    using extension::boost::stream_pragma::RSS;
    const auto kstate_value_putter = [](std::ostream& os, DynamicMonostarUniqueKstate kstate) {
        os << kstate;
    };
    const auto range_streamer_settings_for_basis = RSS<DynamicMonostarUniqueKstate>().set_stream_value_putter(kstate_value_putter);
    auto basis_streamer = kstate::BasisStreamer<DynamicMonostarUniqueKstate>(stream);
    basis_streamer.set_range_streamer_settings(range_streamer_settings_for_basis);
    basis_streamer.stream(state);
    return stream;
}

}  // namespace model_monostar

#endif
