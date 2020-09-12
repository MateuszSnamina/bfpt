#ifndef KSTATE_BASIS_STREAMER_HPP
#define KSTATE_BASIS_STREAMER_HPP

#include <kstate/basis.hpp>
#include <kstate/kstate_streamer.hpp>

#include <iomanip>

// #######################################################################
// ##  BasisStreamer                                                    ##
// #######################################################################

namespace kstate {

template<typename KstateT>
class BasisStreamer {
public:
    BasisStreamer(std::ostream& os);
    BasisStreamer<KstateT>& stream(const Basis<KstateT>&);
    std::ostream& ostream();
    const std::ostream& ostream() const;
    BasisStreamer<KstateT>& set_range_streamer_settings(extension::boost::RangeStreamerSettings<KstateT>);

private:
    extension::boost::RangeStreamerSettings<KstateT> _range_streamer_settings;
    std::ostream& _os;
};

// ***********************************************************************

template<typename KstateT>
BasisStreamer<KstateT>::BasisStreamer(std::ostream& os)
    : _os(os) {}

// ***********************************************************************

template<typename KstateT>
BasisStreamer<KstateT>&
BasisStreamer<KstateT>::set_range_streamer_settings(
        extension::boost::RangeStreamerSettings<KstateT> _) {
    _range_streamer_settings = _;
    return *this;
}

//// TODO refactor the whole funciton!!!
template<typename KstateT>
BasisStreamer<KstateT>&
BasisStreamer<KstateT>::stream(const Basis<KstateT>& basis) {
    using extension::boost::stream_pragma::RSS;
    using kstate::pramga::operator||;
    using kstate::pramga::operator<<;
    // Defaults for basis range streamer:
    const ::std::function<void(::std::ostream&)> default_stream_preparer =
            [](std::ostream& s) { s << "𝔹𝔸𝕊𝕀𝕊-BEGIN" << std::endl; };
    const ::std::function<void(::std::ostream&, size_t)> default_stream_sustainer =
            [](std::ostream& s, size_t i) { s << " - " << std::right << std::setw(6) << i << " : "; };
    //    const ::std::function<void(::std::ostream&, KstateT)> default_stream_value_putter =
    //            [](::std::ostream& s, KstateT ks) { s << ks; };
    const ::std::function<void(::std::ostream&)> default_stream_separer =
            [](std::ostream& s) { s << std::endl; };
    const ::std::function<void(::std::ostream&)> default_stream_finisher =
            [](std::ostream& s) { s << std::endl << "𝔹𝔸𝕊𝕀𝕊-END" << std::endl; };
    bool default_format_independence_flag = true;
    // Apply overrules:
    const auto basis_stream_preparer = (
                _range_streamer_settings._stream_preparer ?
                    *_range_streamer_settings._stream_preparer :
                    default_stream_preparer);
    const auto basis_stream_sustainer = (
                _range_streamer_settings._stream_sustainer ?
                    *_range_streamer_settings._stream_sustainer :
                    default_stream_sustainer);
    //    const auto basis_stream_value_putter = (
    //                _range_streamer_settings._stream_value_putter ?
    //                *_range_streamer_settings._stream_value_putter :
    //                    default_stream_value_putter);
    assert(_range_streamer_settings._stream_value_putter);
    const auto basis_stream_value_putter = *_range_streamer_settings._stream_value_putter ;
    const auto basis_stream_separer = (
                _range_streamer_settings._stream_separer ?
                    *_range_streamer_settings._stream_separer :
                    default_stream_separer);
    const auto basis_stream_finisher = (
                _range_streamer_settings._stream_finisher ?
                    *_range_streamer_settings._stream_finisher :
                    default_stream_finisher);
    const auto basis_format_independence_flag = (
                _range_streamer_settings._format_independence_flag ?
                    *_range_streamer_settings._format_independence_flag :
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
        basis_stream_value_putter(_os, kstate );
    }
    basis_stream_finisher(_os);
    return *this;
}

template<typename KstateT>
std::ostream&
BasisStreamer<KstateT>::ostream() {
    return _os;
}

template<typename KstateT>
const std::ostream&
BasisStreamer<KstateT>::ostream() const {
    return _os;
}

}  // namespace kstate

// #######################################################################
// ##  BasisStringStreamer                                              ##
// #######################################################################

namespace kstate {

template<typename KstateT>
class BasisStringStreamer {
public:
    BasisStringStreamer();
    BasisStringStreamer& stream(const Basis<KstateT>&);
    BasisStringStreamer<KstateT>& set_range_streamer_settings(extension::boost::RangeStreamerSettings<KstateT>);

public:  // Function to retreive the streaming result
    std::string str() const;

private:
    std::ostringstream _oss;
    BasisStreamer<KstateT> _bs;
};

// ***********************************************************************

template<typename KstateT>
BasisStringStreamer<KstateT>::BasisStringStreamer()
    : _bs(_oss) {}

// ***********************************************************************

template<typename KstateT>
BasisStringStreamer<KstateT>&
BasisStringStreamer<KstateT>::set_range_streamer_settings(
        extension::boost::RangeStreamerSettings<KstateT> _) {
    _bs.set_range_streamer_settings(_);
    return *this;
}

template<typename KstateT>
BasisStringStreamer<KstateT>&
BasisStringStreamer<KstateT>::stream(const Basis<KstateT>& basis) {
    _bs.stream(basis);
    return *this;
}

template<typename KstateT>
std::string
BasisStringStreamer<KstateT>::str() const {
    return _oss.str();
}

}  // namespace kstate

#endif
