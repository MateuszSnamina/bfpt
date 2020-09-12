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

//TODO restore:

//inline std::ostream&
//operator<<(std::ostream& stream, const DynamicMonostarUniqueKstateBasis& state) {
//    auto basis_streamer = kstate::BasisStreamer(stream).set_range_streamer_settings_for_kstate(monostar_kstate_range_streamer_settings);
//    basis_streamer.stream(state);
//    return stream;
//}

}  // namespace model_monostar

#endif
