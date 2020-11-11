#pragma once

#include <monostar_system/monostar_kstate.hpp>

#include <kbasis/basis_streamer.hpp>
#include <kbasis/basis.hpp>

// #######################################################################
// ## MonostarKstateBasis                                               ##
// #######################################################################

namespace monostar_system {

using MonostarKstateBasis = kbasis::Basis<MonostarKstateTrait>;

}  // namespace monostar_system

// #######################################################################
// ## MonostarKstateBasis - printing                                    ##
// #######################################################################

namespace monostar_system {

inline std::ostream&
operator<<(std::ostream& stream, const MonostarKstateBasis& basis) {
    using extension::boost::stream_pragma::RSS;
    using kbasis::pramga::operator&&;
    using kbasis::pramga::operator<<;
    const auto kstate_value_putter = [](std::ostream& os, MonostarKstate kstate) {
        os << kstate;
    };
    const auto range_streamer_settings = RSS<MonostarKstate>().set_stream_value_putter(kstate_value_putter);
    stream << (basis && range_streamer_settings);
    return stream;
}

}  // namespace monostar_system
