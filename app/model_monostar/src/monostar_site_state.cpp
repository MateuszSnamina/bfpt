#include <model_monostar/monostar_site_state.hpp>

// #######################################################################
// ## MonostarSiteState                                                 ##
// #######################################################################

namespace model_monostar {

const MonostarSiteState gs{false};
const MonostarSiteState es{true};
const std::vector<MonostarSiteState> ordered_site_states{gs, es};//TODO remove

}  // namespace model_monostar
