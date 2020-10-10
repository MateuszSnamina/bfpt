#include <monostar_system/monostar_site_state.hpp>

// #######################################################################
// ## MonostarSiteState                                                 ##
// #######################################################################

namespace monostar_system {

const MonostarSiteState gs{false};
const MonostarSiteState es{true};
const std::vector<MonostarSiteState> ordered_site_states{gs, es};//TODO remove

}  // namespace monostar_system
