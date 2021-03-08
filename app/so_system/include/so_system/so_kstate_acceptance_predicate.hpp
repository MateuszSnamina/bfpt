#pragma once

#include <so_system/so_kstate.hpp>

// #######################################################################
// ## get_n_{spin,orbit}_site_excitations                               ##
// #######################################################################

namespace so_system {

unsigned get_n_spin_site_excitations(const SoKstate& kstate);
unsigned get_n_orbit_site_excitations(const SoKstate& kstate);
bool are_all_orbit_site_excitations_near_domain_wall(const SoKstate& kstate);

}  // namespace so_system

// #######################################################################
// ## {Spin,Orbit}SiteExcitationsAcceptancePredicate                    ##
// #######################################################################

namespace so_system {

class SpinSiteExcitationsAcceptancePredicate {
   public:
    explicit SpinSiteExcitationsAcceptancePredicate(std::optional<unsigned> n_max_site_spin_excitations);
    bool operator()(const SoKstate& kstate) const;

   private:
    const std::optional<unsigned> _n_max_site_spin_excitations;
};

class OrbitSiteExcitationsAcceptancePredicate {
   public:
    explicit OrbitSiteExcitationsAcceptancePredicate(std::optional<unsigned> n_max_site_orbit_excitations);
    bool operator()(const SoKstate& kstate) const;

   private:
    const std::optional<unsigned> _n_max_site_orbit_excitations;
};

class AllOrbitSiteExcitationsNearDomainWallAcceptancePredicate {
   public:
    bool operator()(const SoKstate& kstate) const;
};

class AcceptancePredicate {
   public:
    AcceptancePredicate(
        std::optional<unsigned> n_max_site_spin_excitations,
        std::optional<unsigned> n_max_site_orbit_excitations,
        bool accept_orbit_site_excitations_only_if_near_domain_wall = false);
    bool operator()(const SoKstate& kstate) const;

   private:
    const SpinSiteExcitationsAcceptancePredicate _spin_site_excitations_acceptance_predicate;
    const OrbitSiteExcitationsAcceptancePredicate _orbit_site_excitations_acceptance_predicate;
    const AllOrbitSiteExcitationsNearDomainWallAcceptancePredicate _all_orbit_site_excitations_near_domain_wall_acceptance_predicate;
    const bool _accept_orbit_site_excitations_only_if_near_domain_wall;
};

}  // namespace so_system
