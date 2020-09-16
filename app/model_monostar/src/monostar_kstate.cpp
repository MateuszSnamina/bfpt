#include <model_monostar/monostar_kstate.hpp>

// #######################################################################
// ## monostar_kstate_range_streamer_settings                           ##
// #######################################################################

namespace model_monostar {

const extension::boost::RangeStreamerSettings<MonostarSiteState> monostar_kstate_range_streamer_settings =
    extension::boost::RangeStreamerSettings<MonostarSiteState>()
        .set_stream_preparer([](::std::ostream& s) { s << "𝕊𝕥𝕒𝕣𝕂𝕤𝕥𝕒𝕥𝕖⦃"; })
        .set_stream_sustainer([](::std::ostream& s, size_t i) {})
        .set_stream_separer([](::std::ostream& s) {})
        .set_stream_finisher([](::std::ostream& s) { s << "⦄"; })
        .set_format_independence_flag(true);
}
