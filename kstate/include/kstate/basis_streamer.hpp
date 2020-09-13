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

// ***********************************************************************

template<typename _BasisT>
class BasisStreamer;

// ***********************************************************************

template <typename BasisT>
BasisStreamer<BasisT>
make_basis_streamer(
        BasisT&&,
        const extension::boost::RangeStreamerSettings<typename remove_cvref_t<BasisT>::Element>);

// ***********************************************************************

template<typename _BasisT>
class BasisStreamer {
public: // Helper types:
    using BasisT = _BasisT;
    using ElementT = typename remove_cvref_t<BasisT>::Element;
public: // Factory function:
    friend
    BasisStreamer<BasisT> make_basis_streamer<BasisT>(BasisT&&, const extension::boost::RangeStreamerSettings<ElementT>);
public: // API:
    ::std::ostream& stream(std::ostream&) const;
    std::string str() const;
private:
    BasisStreamer(std::add_lvalue_reference_t<BasisT>, extension::boost::RangeStreamerSettings<ElementT>);
private:
    const _BasisT _basis;
    const extension::boost::RangeStreamerSettings<ElementT> _range_streamer_settings;
};

// ***********************************************************************

template<typename _BasisT>
BasisStreamer<_BasisT>::BasisStreamer(
        std::add_lvalue_reference_t<_BasisT> basis,
        extension::boost::RangeStreamerSettings<ElementT> range_streamer_settings) :
    _basis(std::forward<_BasisT>(basis)),
    _range_streamer_settings(range_streamer_settings) {
}

// ***********************************************************************

template<typename _BasisT>
::std::ostream&
BasisStreamer<_BasisT>::stream(std::ostream& os) const {
    // Defaults for basis range streamer:
    const ::std::function<void(::std::ostream&)> default_stream_preparer =
            [](std::ostream& s) { s << "𝔹𝔸𝕊𝕀𝕊-BEGIN" << std::endl; };
    const ::std::function<void(::std::ostream&, size_t)> default_stream_sustainer =
            [](std::ostream& s, size_t i) { s << " - " << std::right << std::setw(8) << i << " : "; };
    const ::std::function<void(::std::ostream&, ElementT)> default_stream_value_putter =
            [](::std::ostream& s, ElementT) { s << "[PLACE-FOR-KSTATE]"; };
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

template<typename _BasisT>
std::string
BasisStreamer<_BasisT>::str() const {
    std::ostringstream oss;
    stream(oss);
    return oss.str();
}

// ***********************************************************************

template<typename BasisT>
BasisStreamer<BasisT> make_basis_streamer(
        BasisT&& basis,
        const extension::boost::RangeStreamerSettings<typename remove_cvref_t<BasisT>::Element> range_streamer_settings) {
    return BasisStreamer<BasisT>(std::forward<BasisT>(basis), range_streamer_settings);
}

} // end of namespace kstate

// #######################################################################
// ##  pragma: operator||, operator<<                                   ##
// #######################################################################

namespace kstate::pramga {

template<typename BasisT>
BasisStreamer<BasisT>
operator&&(
        BasisT&& basis,
        extension::boost::RangeStreamerSettings<typename remove_cvref_t<BasisT>::Element> range_streamer_settings) {
    return make_basis_streamer<BasisT>(std::forward<BasisT>(basis), range_streamer_settings);
}

template<typename BasisT>
std::ostream& operator<<(
        std::ostream& os,
        const BasisStreamer<BasisT>& basis_streamer) {
    basis_streamer.stream(os);
    return os;
}

}  // namespace kstate::pramga

#endif
