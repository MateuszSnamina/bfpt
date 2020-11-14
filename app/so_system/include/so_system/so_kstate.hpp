#pragma once

#include <so_system/so_site_state.hpp>

#include <kstate_impl/kstate_concrete_stl.hpp>
#include <kstate_impl/kstate_concrete_integral.hpp>

#include <iostream>

// #######################################################################
// ## SoKstate                                                          ##
// #######################################################################

namespace so_system {

using SoDynamicStlKstate = kstate_impl::DynamicStlKstate<SoSiteStateTrait>;

template <size_t N>
using StaticStlSoKstate = kstate_impl::StaticStlKstate<SoSiteStateTrait, N>;

using SoDynamicIntegral64Kstate = kstate_impl::DynamicIntegral64Kstate<SoSiteStateTrait, 2u>;

//template<size_t N>
//using StaticIntegralSoKstate = kstate_impl::StaticTwoLevelIntegral64Kstate<SoSiteStateTrait, N, 1u>;
// ^^^^ TODO: UNCOMMENT WHEN IT WILL BE IMPLEMENTED

/*
 *  Select the type you want to use:
 */

//using SoKstate = SoDynamicStlKstate;
//using SoKstate = StaticStlSoKstate<20>;
using SoKstate = SoDynamicIntegral64Kstate;

}  // namespace so_system

// #######################################################################
// ## SoKstateTrait - trait                                             ##
// #######################################################################

namespace so_system {

using SoKstateTrait = kstate_trait::TraitKstate<SoKstate>;

}

// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

namespace so_system {

SoKstate classical_gs_kstate(const unsigned n_sites);

SoKstate classical_es_kstate(const unsigned n_sites);

}  // namespace so_system

// #######################################################################
// ## SoKstate - printing                                               ##
// #######################################################################

namespace so_system {

std::ostream& operator<<(std::ostream& stream, const SoKstate& state);

}  // namespace so_system

// #######################################################################
// ## get_n_{spin,orbit}_site_excitations                               ##
// #######################################################################

namespace so_system {

inline unsigned get_n_spin_site_excitations(const SoKstate& kstate) {
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

inline unsigned get_n_orbit_site_excitations(const SoKstate& kstate) {
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

class SpinSiteExcitationsAcceptancePredicate {
   public:
    explicit SpinSiteExcitationsAcceptancePredicate(std::optional<unsigned> n_max_site_spin_excitations)
        : _n_max_site_spin_excitations(n_max_site_spin_excitations) {
    }
    bool operator()(const SoKstate& kstate) const {
        return (
            _n_max_site_spin_excitations
                ? get_n_spin_site_excitations(kstate) < *_n_max_site_spin_excitations
                : true);
    }

   private:
    const std::optional<unsigned> _n_max_site_spin_excitations;
};

class OrbitSiteExcitationsAcceptancePredicate {
   public:
    explicit OrbitSiteExcitationsAcceptancePredicate(std::optional<unsigned> n_max_site_orbit_excitations)
        : _n_max_site_orbit_excitations(n_max_site_orbit_excitations) {
    }
    bool operator()(const SoKstate& kstate) const {
        return (
            _n_max_site_orbit_excitations
                ? get_n_orbit_site_excitations(kstate) < *_n_max_site_orbit_excitations
                : true);
    }

   private:
    const std::optional<unsigned> _n_max_site_orbit_excitations;
};

class AcceptancePredicate {
   public:
    AcceptancePredicate(
        std::optional<unsigned> n_max_site_spin_excitations,
        std::optional<unsigned> n_max_site_orbit_excitations)
        : _spin_site_excitations_acceptance_predicate(n_max_site_spin_excitations),
          _orbit_site_excitations_acceptance_predicate(n_max_site_orbit_excitations) {
    }
    bool operator()(const SoKstate& kstate) const {
        return _spin_site_excitations_acceptance_predicate(kstate) &&
               _orbit_site_excitations_acceptance_predicate(kstate);
    }

   private:
    const SpinSiteExcitationsAcceptancePredicate _spin_site_excitations_acceptance_predicate;
    const OrbitSiteExcitationsAcceptancePredicate _orbit_site_excitations_acceptance_predicate;
};

}  // namespace so_system
