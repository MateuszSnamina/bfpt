#include <monostar_system/monostar_kstate.hpp>

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
    const auto generator = MonostarKstateTrait::from_range(generator_array);
    const auto generator_view = MonostarKstateTrait::to_view(generator);
    const size_t generator_n_unique_shift = MonostarKstateTrait::view_n_unique_shift(generator_view);
    const auto rotation_spec = kstate_view_amend_spec::rotated(generator_n_unique_shift);
    const auto generator_view_unique_shifted = MonostarKstateTrait::rotated_view(generator_view, rotation_spec);
    const auto uniquely_shifted_generator = MonostarKstateTrait::from_view(generator_view_unique_shifted);
    return uniquely_shifted_generator;
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
