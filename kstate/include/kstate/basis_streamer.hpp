#ifndef KSTATE_BASIS_STREAMER_HPP
#define KSTATE_BASIS_STREAMER_HPP

#include <kstate/basis.hpp>
#include <kstate/kstate_streamer.hpp>

#include <boost/range/adaptor/indirected.hpp>

#include <iomanip>

// #######################################################################
// ##  BasisStreamer                                                    ##
// #######################################################################

namespace kstate {

template<typename _KstateT>
class BasisStreamer {
public: // Helper types:
    using KstateT = _KstateT;
public:
    BasisStreamer(const Basis<KstateT>& basis, extension::boost::RangeStreamerSettings<KstateT> range_streamer_settings);
    ::std::ostream& stream(std::ostream&) const;
    std::string str() const;
private:
    const Basis<KstateT>& _basis;
    const extension::boost::RangeStreamerSettings<KstateT> _range_streamer_settings;
};

// ***********************************************************************

template<typename _KstateT>
BasisStreamer<_KstateT>::BasisStreamer(
        const Basis<KstateT>& basis,
        extension::boost::RangeStreamerSettings<KstateT> range_streamer_settings) :
    _basis(basis),
    _range_streamer_settings(range_streamer_settings) {
}

// ***********************************************************************

template<typename _KstateT>
::std::ostream&
BasisStreamer<_KstateT>::stream(std::ostream& os) const {
    // Defaults for basis range streamer:
    const ::std::function<void(::std::ostream&)> default_stream_preparer =
            [](std::ostream& s) { s << "𝔹𝔸𝕊𝕀𝕊-BEGIN" << std::endl; };
    const ::std::function<void(::std::ostream&, size_t)> default_stream_sustainer =
            [](std::ostream& s, size_t i) { s << " - " << std::right << std::setw(6) << i << " : "; };
    const ::std::function<void(::std::ostream&, KstateT)> default_stream_value_putter =
            [](::std::ostream& s, KstateT) { s << "[PLACE-FOR-KSTATE]"; };
    const ::std::function<void(::std::ostream&)> default_stream_separer =
            [](std::ostream& s) { s << std::endl; };
    const ::std::function<void(::std::ostream&)> default_stream_finisher =
            [](std::ostream& s) { s << std::endl << "𝔹𝔸𝕊𝕀𝕊-END" << std::endl; };
    bool default_format_independence_flag = true;
    // Apply overrules:
    const auto stream_preparer = (
                _range_streamer_settings._stream_preparer ?
                    *_range_streamer_settings._stream_preparer :
                    default_stream_preparer);
    const auto stream_sustainer = (
                _range_streamer_settings._stream_sustainer ?
                    *_range_streamer_settings._stream_sustainer :
                    default_stream_sustainer);
    const auto stream_value_putter = (
                _range_streamer_settings._stream_value_putter ?
                    *_range_streamer_settings._stream_value_putter :
                    default_stream_value_putter);
    const auto stream_separer = (
                _range_streamer_settings._stream_separer ?
                    *_range_streamer_settings._stream_separer :
                    default_stream_separer);
    const auto stream_finisher = (
                _range_streamer_settings._stream_finisher ?
                    *_range_streamer_settings._stream_finisher :
                    default_stream_finisher);
    const auto format_independence_flag = (
                _range_streamer_settings._format_independence_flag ?
                    *_range_streamer_settings._format_independence_flag :
                    default_format_independence_flag);
    // Stream:
    extension::boost::stream_range_impl<decltype(_basis.vec_index() | ::boost::adaptors::indirected)>(
                _basis.vec_index() | ::boost::adaptors::indirected,
                os,
                stream_preparer,
                stream_sustainer,
                stream_value_putter,
                stream_separer,
                stream_finisher,
                format_independence_flag);
    // Return:
    return os;
}

// ***********************************************************************

template<typename _KstateT>
std::string
BasisStreamer<_KstateT>::str() const {
    std::ostringstream oss;
    stream(oss);
    return oss.str();
}

}  // namespace kstate

#endif
