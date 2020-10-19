#include <monostar_system/monostar_kstate.hpp>

#include <kstate_op_range/op_range_unique_shift.hpp>

#include <kstate_impl/kstate_streamer.hpp>

#include <vector>

// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

namespace monostar_system {

MonostarKstate classical_gs_kstate(const unsigned n_sites) {
    std::vector<MonostarSiteState> generator_array(n_sites, monostar_system::gs);
    return MonostarKstateTrait::from_range(generator_array);
}

MonostarKstate classical_es_kstate(const unsigned n_sites) {
    assert(n_sites > 0);
    std::vector<MonostarSiteState> generator_array(n_sites, monostar_system::gs);
    generator_array[0] = es;
    return MonostarKstateTrait::from_range(kstate_op_range::make_unique_shift(generator_array));
}

}  // namespace monostar_system


// #######################################################################
// ## DynamicMonostarKstate - printing                                  ##
// #######################################################################

namespace monostar_system {

std::ostream& operator<<(std::ostream& stream, const MonostarKstate& state) {
    using extension::boost::stream_pragma::RSS;
    using kstate_impl::pramga::operator||;
    using kstate_impl::pramga::operator<<;

    const extension::boost::RangeStreamerSettings<MonostarSiteState> monostar_kstate_range_streamer_settings =
            extension::boost::RangeStreamerSettings<MonostarSiteState>()
            .set_stream_preparer([](::std::ostream& s) { s << "𝕊𝕥𝕒𝕣𝕂𝕤𝕥𝕒𝕥𝕖⦃"; })
            .set_stream_sustainer([](::std::ostream&, size_t) {})
            .set_stream_separer([](::std::ostream&) {})
            .set_stream_finisher([](::std::ostream& s) { s << "⦄"; })
            .set_format_independence_flag(true);

    stream << (state || monostar_kstate_range_streamer_settings);
    return stream;
}

}  // namespace monostar_system
