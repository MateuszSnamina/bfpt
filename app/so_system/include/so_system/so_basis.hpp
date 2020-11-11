#pragma once

#include <so_system/so_kstate.hpp>

#include <kbasis/basis_streamer.hpp>
#include <kbasis/basis.hpp>

// #######################################################################
// ## SoKstateBasis                                                     ##
// #######################################################################

namespace so_system {

using SoKstateBasis = kbasis::Basis<SoKstateTrait>;

}  // namespace so_system

// #######################################################################
// ## SoKstateBasis - printing                                          ##
// #######################################################################

namespace so_system {

inline std::ostream&
operator<<(std::ostream& stream, const SoKstateBasis& basis) {
    using extension::boost::stream_pragma::RSS;
    using kbasis::pramga::operator&&;
    using kbasis::pramga::operator<<;
    const auto kstate_value_putter = [](std::ostream& os, SoKstate kstate) {
        os << kstate;
    };
    const auto range_streamer_settings = RSS<SoKstate>().set_stream_value_putter(kstate_value_putter);
    stream << (basis && range_streamer_settings);
    return stream;
}

}  // namespace so_system
