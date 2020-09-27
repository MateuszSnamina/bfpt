#ifndef MODEL_MONOSTAR_MONOSTAR_KSTATE_HPP
#define MODEL_MONOSTAR_MONOSTAR_KSTATE_HPP

#include <model_monostar/monostar_site_state.hpp>

#include <kstate/unique_shift.hpp>
#include <kstate/kstate_concrete.hpp>
#include <kstate/kstate_streamer.hpp>

#include <extensions/range_streamer.hpp>

#include <vector>

// #######################################################################
// ## monostar_kstate_range_streamer_settings                           ##
// #######################################################################

namespace model_monostar {

extern const extension::boost::RangeStreamerSettings<MonostarSiteState> monostar_kstate_range_streamer_settings;

}  // namespace model_monostar

// #######################################################################
// ## DynamicMonostarKstate                                             ##
// #######################################################################

namespace model_monostar {

using DynamicMonostarKstate = kstate::DynamicKstate<MonostarSiteState>;

}  // namespace model_monostar

// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

namespace model_monostar {

inline DynamicMonostarKstate classical_gs_kstate(const unsigned n_sites) {
    std::vector<MonostarSiteState> generator_array(n_sites, model_monostar::gs);
    return DynamicMonostarKstate(kstate::make_unique_shift(generator_array), kstate::ctr_from_range);
}

inline DynamicMonostarKstate classical_es_kstate(const unsigned n_sites) {
    std::vector<MonostarSiteState> generator_array(n_sites, model_monostar::gs);
    generator_array[0] = es;
    return DynamicMonostarKstate(kstate::make_unique_shift(generator_array), kstate::ctr_from_range);
}

}  // namespace model_monostar

// #######################################################################
// ## DynamicMonostarKstate - printing                                  ##
// #######################################################################

namespace model_monostar {

inline std::ostream&
operator<<(std::ostream& stream, const DynamicMonostarKstate& state) {
    using extension::boost::stream_pragma::RSS;
    using kstate::pramga::operator||;
    using kstate::pramga::operator<<;
    stream << (state || monostar_kstate_range_streamer_settings);
    return stream;
}

}  // namespace model_monostar

#endif
