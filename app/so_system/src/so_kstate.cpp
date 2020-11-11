#include <so_system/so_kstate.hpp>

#include <kstate_impl/kstate_streamer.hpp>

#include <vector>

// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

//namespace so_system {

//SoKstate classical_gs_kstate(const unsigned n_sites) {
//    std::vector<SoSiteState> generator_array(n_sites, so_system::gg);//TODO remove namespace?
//    return SoKstateTrait::from_range(generator_array);
//}

//SoKstate classical_es_kstate(const unsigned n_sites) {
//    assert(n_sites > 0);
//    std::vector<SoSiteState> generator_array(n_sites, so_system::gg);//TODO remove namespace?
//    generator_array[0] = so_system::eg;//TODO remove namespace?
//    const auto generator = SoKstateTrait::from_range(generator_array);
//    const auto generator_view = SoKstateTrait::to_view(generator);
//    const size_t generator_n_unique_shift = SoKstateTrait::view_n_unique_shift(generator_view);
//    const auto rotation_spec = kstate_view_amend_spec::rotated(generator_n_unique_shift);
//    const auto generator_view_unique_shifted = SoKstateTrait::rotated_view(generator_view, rotation_spec);
//    const auto uniquely_shifted_generator = SoKstateTrait::from_view(generator_view_unique_shifted);
//    return uniquely_shifted_generator;
//}

//}  // namespace so_system

// #######################################################################
// ## SoKstate - printing                                               ##
// #######################################################################

namespace so_system {

std::ostream& operator<<(std::ostream& stream, const SoKstate& state) {
    using extension::boost::stream_pragma::RSS;
    using kstate_impl::pramga::operator||;
    using kstate_impl::pramga::operator<<;

    const extension::boost::RangeStreamerSettings<SoSiteState> so_kstate_range_streamer_settings =
        extension::boost::RangeStreamerSettings<SoSiteState>()
            .set_stream_preparer([](::std::ostream& s) { s << "𝕊𝕥𝕒𝕣𝕂𝕤𝕥𝕒𝕥𝕖⦃"; })
            .set_stream_sustainer([](::std::ostream&, size_t) {})
            .set_stream_separer([](::std::ostream&) {})
            .set_stream_finisher([](::std::ostream& s) { s << "⦄"; })
            .set_format_independence_flag(true);

    stream << (state || so_kstate_range_streamer_settings);
    return stream;
}

}  // namespace so_system
