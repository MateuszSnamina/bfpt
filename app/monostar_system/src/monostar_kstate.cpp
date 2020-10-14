#include <monostar_system/monostar_kstate.hpp>

#include <kstate_range_op/range_op_unique_shift.hpp>

#include <kstate_impl/kstate_streamer.hpp>

#include <vector>

// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

namespace monostar_system {

DynamicMonostarKstate classical_gs_kstate(const unsigned n_sites) {
    //std::vector<MonostarSiteState> generator_array(n_sites, monostar_system::gs);
    std::array<MonostarSiteState, 20> generator_array {
        monostar_system::gs, monostar_system::gs, monostar_system::gs, monostar_system::gs,
                monostar_system::gs, monostar_system::gs, monostar_system::gs, monostar_system::gs,
                monostar_system::gs, monostar_system::gs, monostar_system::gs, monostar_system::gs,
                monostar_system::gs, monostar_system::gs, monostar_system::gs, monostar_system::gs,
                monostar_system::gs, monostar_system::gs, monostar_system::gs, monostar_system::gs
    };
    return DynamicMonostarKstateTrait::from_range(generator_array);
    //return DynamicMonostarKstateTrait::from_range(generator_array);

}

DynamicMonostarKstate classical_es_kstate(const unsigned n_sites) {
    //std::vector<MonostarSiteState> generator_array(n_sites, monostar_system::gs);
    //generator_array[0] = es;
    std::array<MonostarSiteState, 20> generator_array {
        monostar_system::es, monostar_system::gs, monostar_system::gs, monostar_system::gs,
                monostar_system::gs, monostar_system::gs, monostar_system::gs, monostar_system::gs,
                monostar_system::gs, monostar_system::gs, monostar_system::gs, monostar_system::gs,
                monostar_system::gs, monostar_system::gs, monostar_system::gs, monostar_system::gs,
                monostar_system::gs, monostar_system::gs, monostar_system::gs, monostar_system::gs
    };
    return DynamicMonostarKstateTrait::from_range(kstate_range_op::make_unique_shift(generator_array));
}

}  // namespace monostar_system


// #######################################################################
// ## DynamicMonostarKstate - printing                                  ##
// #######################################################################

namespace monostar_system {

std::ostream& operator<<(std::ostream& stream, const DynamicMonostarKstate& state) {
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
