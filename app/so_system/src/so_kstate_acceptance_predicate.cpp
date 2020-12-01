#include <so_system/so_kstate_acceptance_predicate.hpp>

// #######################################################################
// ## get_n_{spin,orbit}_site_excitations                               ##
// #######################################################################

namespace so_system {

unsigned get_n_spin_site_excitations(const SoKstate& kstate) {
    const unsigned n_sites = SoKstateTrait::n_sites(kstate);
    const auto kstate_view = SoKstateTrait::to_view(kstate);
    unsigned n_spin_site_excitations = 0;
    for (size_t idx = 0; idx < n_sites; idx++) {
        const auto site = SoKstateTrait::view_n_th_site_state(kstate_view, idx);
        if (site.is_spin_excited()) {
            n_spin_site_excitations++;
        }
    }
    return n_spin_site_excitations;
}

unsigned get_n_orbit_site_excitations(const SoKstate& kstate) {
    const unsigned n_sites = SoKstateTrait::n_sites(kstate);
    const auto kstate_view = SoKstateTrait::to_view(kstate);
    unsigned n_orbit_site_excitations = 0;
    for (size_t idx = 0; idx < n_sites; idx++) {
        const auto site = SoKstateTrait::view_n_th_site_state(kstate_view, idx);
        if (site.is_orbit_excited()) {
            n_orbit_site_excitations++;
        }
    }
    return n_orbit_site_excitations;
}

}  // namespace so_system

// #######################################################################
// ## {Spin,Orbit}SiteExcitationsAcceptancePredicate                    ##
// #######################################################################

namespace so_system {

SpinSiteExcitationsAcceptancePredicate::SpinSiteExcitationsAcceptancePredicate(
    std::optional<unsigned> n_max_site_spin_excitations)
    : _n_max_site_spin_excitations(n_max_site_spin_excitations) {
}

bool SpinSiteExcitationsAcceptancePredicate::operator()(const SoKstate& kstate) const {
    return (_n_max_site_spin_excitations
                ? get_n_spin_site_excitations(kstate) <= *_n_max_site_spin_excitations
                : true);
}

OrbitSiteExcitationsAcceptancePredicate::OrbitSiteExcitationsAcceptancePredicate(
    std::optional<unsigned> n_max_site_orbit_excitations)
    : _n_max_site_orbit_excitations(n_max_site_orbit_excitations) {
}

bool OrbitSiteExcitationsAcceptancePredicate::operator()(const SoKstate& kstate) const {
    return (_n_max_site_orbit_excitations
                ? get_n_orbit_site_excitations(kstate) <= *_n_max_site_orbit_excitations
                : true);
}

AcceptancePredicate::AcceptancePredicate(
    std::optional<unsigned> n_max_site_spin_excitations,
    std::optional<unsigned> n_max_site_orbit_excitations)
    : _spin_site_excitations_acceptance_predicate(n_max_site_spin_excitations),
      _orbit_site_excitations_acceptance_predicate(n_max_site_orbit_excitations) {
}

bool AcceptancePredicate::operator()(const SoKstate& kstate) const {
    return _spin_site_excitations_acceptance_predicate(kstate) &&
           _orbit_site_excitations_acceptance_predicate(kstate);
}

}  // namespace so_system
