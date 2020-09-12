#ifndef KSTATE_BASIS_STREAMER_HPP
#define KSTATE_BASIS_STREAMER_HPP

#include <kstate/basis.hpp>
#include <kstate/kstate_streamer.hpp>

#include <iomanip>

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
    extension::boost::RangeStreamerSettings _range_streamer_settings_for_kstate;
    extension::boost::RangeStreamerSettings _range_streamer_settings_for_basis;
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

// TODO refactor the whole funciton!!!
template <typename Element>
BasisStreamer& BasisStreamer::stream(const Basis<Element>& basis) {
    using extension::boost::stream_pragma::operator<<;
    using extension::boost::stream_pragma::RSS;
    using kstate::pramga::operator||;
    // Defaults:
    const ::std::function<void(::std::ostream&)> default_stream_preparer =
            [](std::ostream& s) { s << "𝔹𝔸𝕊𝕀𝕊-BEGIN" << std::endl; };
    const ::std::function<void(::std::ostream&, size_t)> default_stream_sustainer =
            [](std::ostream& s, size_t i) { s << " - " << std::right << std::setw(6) << i << " : "; };
    const ::std::function<void(::std::ostream&)> default_stream_separer =
            [](std::ostream& s) { s << std::endl; };
    const ::std::function<void(::std::ostream&)> default_stream_finisher =
            [](std::ostream& s) { s << std::endl << "𝔹𝔸𝕊𝕀𝕊-END" << std::endl; };
    bool default_format_independence_flag = true;
    // Apply overrules:
    const auto basis_stream_preparer = (
                _range_streamer_settings_for_basis._stream_preparer ?
                *_range_streamer_settings_for_basis._stream_preparer :
                default_stream_preparer);
    const auto basis_stream_sustainer = (
                _range_streamer_settings_for_basis._stream_sustainer ?
                *_range_streamer_settings_for_basis._stream_sustainer :
                    default_stream_sustainer);
    const auto basis_stream_separer = (
                _range_streamer_settings_for_basis._stream_separer ?
                *_range_streamer_settings_for_basis._stream_separer :
                default_stream_separer);
    const auto basis_stream_finisher = (
                _range_streamer_settings_for_basis._stream_finisher ?
                *_range_streamer_settings_for_basis._stream_finisher :
                default_stream_finisher);
    const auto basis_format_independence_flag = (
                _range_streamer_settings_for_basis._format_independence_flag ?
                *_range_streamer_settings_for_basis._format_independence_flag :
                default_format_independence_flag);
    // Stream:
    const extension::std::StreamFromatStacker stream_format_stacker(
                _os, basis_format_independence_flag);
    basis_stream_preparer(_os);
    for (const auto& _ : basis.vec_index() | ::boost::adaptors::indexed(0)) {
        assert(_.value());
        const auto& index = _.index();
        const auto& kstate = *_.value();
        if (index != 0) {
            basis_stream_separer(_os);
        }
        basis_stream_sustainer(_os, index);
        _os << (kstate || _range_streamer_settings_for_kstate);
    }
    basis_stream_finisher(_os);
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
