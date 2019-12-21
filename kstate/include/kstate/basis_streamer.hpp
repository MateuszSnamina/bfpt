#ifndef KSTATE_BASIS_STREAMER_HPP
#define KSTATE_BASIS_STREAMER_HPP

#include <kstate/basis.hpp>
#include <kstate/kstate_streamer.hpp>

#include <iomanip>

// #######################################################################
// ##  basis_default_range_streamer_settings                            ##
// #######################################################################

namespace kstate {

static extension::boost::RangeStreamerSettings basis_default_range_streamer_settings{
    [](std::ostream& s) { s << "𝔹𝔸𝕊𝕀𝕊-BEGIN" << std::endl; },
    [](std::ostream& s, size_t i) { s << " - " << std::right << std::setw(6) << i << " : "; },
    [](std::ostream& s) { s << std::endl; },
    [](std::ostream& s) { s << std::endl
                            << "𝔹𝔸𝕊𝕀𝕊-END" << std::endl; }};
}  // namespace kstate

// #######################################################################
// ##  BasisStreamer                                                    ##
// #######################################################################

namespace kstate {

class BasisStreamer {
   public:
    BasisStreamer(std::ostream& os);
    template <typename Element>
    BasisStreamer& stream(const Basis<Element>&);
    std::ostream& ostream();
    const std::ostream& ostream() const;
    BasisStreamer& set_range_streamer_settings_for_kstate(extension::boost::RangeStreamerSettings);
    BasisStreamer& set_range_streamer_settings_for_basis(extension::boost::RangeStreamerSettings);

   private:
    extension::boost::RangeStreamerSettings _range_streamer_settings_for_kstate = kstate_default_range_streamer_settings;
    extension::boost::RangeStreamerSettings _range_streamer_settings_for_basis = basis_default_range_streamer_settings;
    std::ostream& _os;
};

// ***********************************************************************

inline BasisStreamer::BasisStreamer(std::ostream& os)
    : _os(os) {}

// ***********************************************************************

inline BasisStreamer& BasisStreamer::set_range_streamer_settings_for_kstate(extension::boost::RangeStreamerSettings _) {
    _range_streamer_settings_for_kstate = _;
    return *this;
}

inline BasisStreamer& BasisStreamer::set_range_streamer_settings_for_basis(extension::boost::RangeStreamerSettings _) {
    _range_streamer_settings_for_basis = _;
    return *this;
}

template <typename Element>
BasisStreamer& BasisStreamer::stream(const Basis<Element>& basis) {
    KstateStreamer kss(_os, _range_streamer_settings_for_kstate);
    _range_streamer_settings_for_basis._stream_preparer(_os);
    for (const auto& _ : basis.vec_index() | ::boost::adaptors::indexed(0)) {
        assert(_.value());
        const auto& index = _.index();
        const auto& kstate = *_.value();
        if (index != 0) {
            _range_streamer_settings_for_basis._stream_separer(_os);
        }
        _range_streamer_settings_for_basis._stream_sustainer(_os, index);
        kss.stream(kstate);
    }
    _range_streamer_settings_for_basis._stream_finisher(_os);
    return *this;
}

inline std::ostream& BasisStreamer::ostream() {
    return _os;
}

inline const std::ostream& BasisStreamer::ostream() const {
    return _os;
}

}  // namespace kstate

// #######################################################################
// ##  BasisStringStreamer                                              ##
// #######################################################################

namespace kstate {

class BasisStringStreamer {
   public:
    BasisStringStreamer();
    template <typename Element>
    BasisStringStreamer& stream(const Basis<Element>&);
    BasisStringStreamer& set_range_streamer_settings_for_kstate(extension::boost::RangeStreamerSettings);
    BasisStringStreamer& set_range_streamer_settings_for_basis(extension::boost::RangeStreamerSettings);

   public:  // Function to retreive the streaming result
    std::string str() const;

   private:
    std::ostringstream _oss;
    BasisStreamer _bs;
};

// ***********************************************************************

inline BasisStringStreamer::BasisStringStreamer()
    : _bs(_oss) {}

// ***********************************************************************

inline BasisStringStreamer& BasisStringStreamer::set_range_streamer_settings_for_kstate(extension::boost::RangeStreamerSettings _) {
    _bs.set_range_streamer_settings_for_kstate(_);
    return *this;
}

inline BasisStringStreamer& BasisStringStreamer::set_range_streamer_settings_for_basis(extension::boost::RangeStreamerSettings _) {
    _bs.set_range_streamer_settings_for_basis(_);
    return *this;
}

template <typename Element>
BasisStringStreamer& BasisStringStreamer::stream(const Basis<Element>& basis) {
    _bs.stream(basis);
    return *this;
}

inline std::string BasisStringStreamer::str() const {
    return _oss.str();
}

}  // namespace kstate

#endif