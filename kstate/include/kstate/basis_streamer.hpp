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
    BasisStreamer(const Basis<KstateT>& basis, extension::boost::RangeStreamerSettings<KstateT> range_streamer_settings);
    ::std::ostream& stream(std::ostream&) const;
    std::string str() const;
private:
    const Basis<KstateT>& _basis;
    const extension::boost::RangeStreamerSettings<KstateT> _range_streamer_settings;
};

// ***********************************************************************

template<typename KstateT>
BasisStreamer<KstateT>::BasisStreamer(
        const Basis<KstateT>& basis,
        extension::boost::RangeStreamerSettings<KstateT> range_streamer_settings) :
    _basis(basis),
    _range_streamer_settings(range_streamer_settings) {
}

// ***********************************************************************

//// TODO refactor the whole funciton!!!
template<typename KstateT>
::std::ostream&
BasisStreamer<KstateT>::stream(std::ostream& os) const {
    using extension::boost::stream_pragma::RSS;    //TODO: remove!
    using kstate::pramga::operator||;             //TODO: remove!
    using kstate::pramga::operator<<;             //TODO: remove!
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
                os, basis_format_independence_flag);
    basis_stream_preparer(os);
    for (const auto& _ : _basis.vec_index() | ::boost::adaptors::indexed(0)) {
        assert(_.value());
        const auto& index = _.index();
        const auto& kstate = *_.value();
        if (index != 0) {
            basis_stream_separer(os);
        }
        basis_stream_sustainer(os, index);
        basis_stream_value_putter(os, kstate );
    }
    basis_stream_finisher(os);

    return os;
}

template<typename KstateT>
std::string
BasisStreamer<KstateT>::str() const {
    std::ostringstream oss;
    stream(oss);
    return oss.str();
}

}  // namespace kstate

#endif
