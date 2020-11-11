#pragma once

#include <kstate_trait/trait_site_state.hpp>

#include <boost/operators.hpp>

#include <iostream>
#include <vector>
#include <cassert>

// #######################################################################
// ## SoSiteState -- class                                              ##
// #######################################################################

namespace so_system {

class SoSiteState : boost::totally_ordered<SoSiteState> {
   public:
    constexpr SoSiteState(bool spin_excitation, bool orbit_excitation);
    constexpr explicit SoSiteState(unsigned char idx);
    constexpr bool operator<(const SoSiteState& other) const;
    constexpr bool operator==(const SoSiteState& other) const;
    friend std::ostream& operator<<(std::ostream&, const SoSiteState&);
    constexpr unsigned char get_index() const;
    constexpr bool is_spin_excited() const;
    constexpr bool is_orbit_excited() const;

   private:
    constexpr static unsigned char spin_bit_mask = 1u;
    constexpr static unsigned char orbit_bit_mask = 2u;
    unsigned char _idx;
    constexpr unsigned char excitation_bools_to_idx(bool spin_excitation, bool orbit_excitation);
};

// ***********************************************************************

inline constexpr SoSiteState::SoSiteState(bool spin_excitation, bool orbit_excitation)
    : SoSiteState(excitation_bools_to_idx(spin_excitation, orbit_excitation)) {
}

inline constexpr SoSiteState::SoSiteState(unsigned char idx)
    : _idx(idx) {
    assert(idx < 4u);
}

inline constexpr unsigned char SoSiteState::excitation_bools_to_idx(bool spin_excitation, bool orbit_excitation) {
    return (spin_excitation ? spin_bit_mask : 0u) | (orbit_excitation ? orbit_bit_mask : 0u);
}

constexpr unsigned char SoSiteState::get_index() const {
    return _idx;
}

inline constexpr bool SoSiteState::is_spin_excited() const {
    return _idx & spin_bit_mask;
}

inline constexpr bool SoSiteState::is_orbit_excited() const {
    return _idx & orbit_bit_mask;
}

inline constexpr bool SoSiteState::operator<(const SoSiteState& other) const {
    return (this->_idx < other._idx);
}

inline constexpr bool SoSiteState::operator==(const SoSiteState& other) const {
    return (this->_idx == other._idx);
}

inline std::ostream&
operator<<(std::ostream& stream, const SoSiteState& state) {
    if (!state.is_spin_excited() && !state.is_orbit_excited()) {
        stream << '_';
    } else if (!state.is_spin_excited() && state.is_orbit_excited()) {
        stream << '.';
    } else if (state.is_spin_excited() && !state.is_orbit_excited()) {
        stream << '*';
    } else if (state.is_spin_excited() && state.is_orbit_excited()) {
        stream << 'o';
    }
    return stream;
}

}  // namespace so_system

// #######################################################################
// ## SoSiteState -- global values                                      ##
// #######################################################################

namespace so_system {

constexpr SoSiteState gg{false, false};  // spin: ground + orbit: ground
constexpr SoSiteState ge{false, true};   // spin: ground + orbit: excited
constexpr SoSiteState eg{true, false};   // spin: excited + orbit: ground
constexpr SoSiteState ee{true, true};    // spin: excited + orbit: excited

}  // namespace so_system

// #######################################################################
// ## SoSiteState -- implement trait                                    ##
// #######################################################################

namespace kstate_trait {

template <>
struct TraitSiteState<so_system::SoSiteState> {
    static constexpr bool is_site_state_trait = true;
    using SiteStateT = so_system::SoSiteState;

    constexpr static unsigned site_basis_dim() {
        return 4u;
    }

    constexpr static unsigned get_index(const SiteStateT& state) {
        return state.get_index();
    }

    constexpr static SiteStateT from_index(unsigned idx) {
        if (idx < site_basis_dim()) {
            return SiteStateT(static_cast<unsigned char>(idx));
        } else {
            throw std::domain_error("Index out of range.");
        }
    }
};

}  // namespace kstate_trait

namespace so_system {

using SoSiteStateTrait = kstate_trait::TraitSiteState<so_system::SoSiteState>;

}  // end of namespace so_system
