#pragma once

#include <monostar_system/monostar_kstate.hpp>

#include <kbasis/basis_streamer.hpp>
#include <kbasis/basis.hpp>

// #######################################################################
// ## DynamicMonostarKstateBasis                                        ##
// #######################################################################

namespace monostar_system {

using DynamicMonostarKstateBasis = kbasis::Basis<kstate_trait::TraitKstate<DynamicMonostarKstate>>;

}  // namespace monostar_system

// #######################################################################
// ## DynamicMonostarKstateBasis - printing                             ##
// #######################################################################

namespace monostar_system {

inline std::ostream&
operator<<(std::ostream& stream, const DynamicMonostarKstateBasis& basis) {
    using extension::boost::stream_pragma::RSS;
    using kbasis::pramga::operator&&;
    using kbasis::pramga::operator<<;
    const auto kstate_value_putter = [](std::ostream& os, DynamicMonostarKstate kstate) {
        os << kstate;
    };
    const auto range_streamer_settings = RSS<DynamicMonostarKstate>().set_stream_value_putter(kstate_value_putter);
    stream << (basis && range_streamer_settings);
    return stream;
}

}  // namespace monostar_system
