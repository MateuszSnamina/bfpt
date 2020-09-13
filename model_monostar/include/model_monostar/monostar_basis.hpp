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
operator<<(std::ostream& stream, const DynamicMonostarUniqueKstateBasis& basis) {
    using extension::boost::stream_pragma::RSS;
    using kstate::pramga::operator&&;
    using kstate::pramga::operator<<;
    const auto kstate_value_putter = [](std::ostream& os, DynamicMonostarUniqueKstate kstate) {
        os << kstate;
    };
    const auto range_streamer_settings = RSS<DynamicMonostarUniqueKstate>().set_stream_value_putter(kstate_value_putter);
    stream << (basis && range_streamer_settings);
    return stream;
}

}  // namespace model_monostar

#endif
