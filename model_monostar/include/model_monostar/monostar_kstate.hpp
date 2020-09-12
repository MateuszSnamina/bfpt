#ifndef MODEL_MONOSTAR_MONOSTAR_KSTATE_HPP
#define MODEL_MONOSTAR_MONOSTAR_KSTATE_HPP

#include <model_monostar/monostar_site_state.hpp>

#include <kstate/kstate_concrete.hpp>
#include <kstate/kstate_streamer.hpp>

#include <extensions/range_streamer.hpp>

#include <vector>

// #######################################################################
// ## monostar_kstate_range_streamer_settings                           ##
// #######################################################################

namespace model_monostar {

extern const extension::boost::RangeStreamerSettings monostar_kstate_range_streamer_settings;

}  // namespace model_monostar

// #######################################################################
// ## DynamicMonostarKstate nad DynamicMonostarUniqueKstate             ##
// #######################################################################

namespace model_monostar {

using DynamicMonostarKstate = kstate::DynamicKstate<MonostarSiteState>;
using DynamicMonostarUniqueKstate = kstate::DynamicUniqueKstate<MonostarSiteState>;

}  // namespace model_monostar

// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

namespace model_monostar {

inline DynamicMonostarUniqueKstate classical_gs_kstate(const unsigned n_sites) {
    std::vector<MonostarSiteState> generator_array(n_sites, model_monostar::gs);
    return DynamicMonostarUniqueKstate(generator_array, kstate::ctr_from_range);
}

inline DynamicMonostarUniqueKstate classical_es_kstate(const unsigned n_sites) {
    std::vector<MonostarSiteState> generator_array(n_sites, model_monostar::gs);
    generator_array[0] = es;
    return DynamicMonostarUniqueKstate(generator_array, kstate::ctr_from_range);
}

}  // namespace model_monostar

// #######################################################################
// ## DynamicMonostarKstate nad DynamicMonostarUniqueKstate - printing  ##
// #######################################################################

namespace model_monostar {

inline std::ostream&
operator<<(std::ostream& stream, const DynamicMonostarKstate& state) {
    using namespace kstate::pramga;
    using namespace extension::boost::stream_pragma;
    kstate::pramga::operator<<(stream, state || RSS());
    //stream << (state || RSS()); //TODO refactor
    return stream;
}

inline std::ostream&
operator<<(std::ostream& stream, const DynamicMonostarUniqueKstate& state) {
    using namespace kstate::pramga;
    using namespace extension::boost::stream_pragma;
    kstate::pramga::operator<<(stream, state || RSS());
    //stream << (state || RSS()); //TODO refactor
    return stream;
}

}  // namespace model_monostar

#endif
