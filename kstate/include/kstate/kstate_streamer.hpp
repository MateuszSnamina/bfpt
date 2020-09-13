#ifndef KSTATE_KSTATE_STREAMER_HPP
#define KSTATE_KSTATE_STREAMER_HPP

#include <kstate/kstate_abstract.hpp>
#include <kstate/remove_cvref.hpp>

#include <extensions/range_streamer.hpp>

#include <boost/range.hpp>
#include <boost/range/adaptor/indexed.hpp>

#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

// The doublestruck font: https://en.wikipedia.org/wiki/Blackboard_bold

// #######################################################################
// ##  KstateRangeStreamerSettings                                      ##
// #######################################################################

namespace kstate {

template<typename _T>
struct KstateRangeStreamerSettings {
    using T = _T;
    KstateRangeStreamerSettings(extension::boost::RangeStreamerSettings<_T> range_streamer_settings) :
        _range_streamer_settings(range_streamer_settings) {
    }
    const extension::boost::RangeStreamerSettings<_T> _range_streamer_settings;
};

}

// #######################################################################
// ##  RangeStreamer                                                    ##
// #######################################################################

namespace kstate {

// ***********************************************************************

template<typename _KstateT>
class KstateRangeStreamer;

// ***********************************************************************

template <typename KstateT>
KstateRangeStreamer<KstateT>
make_kstate_range_streamer(
        KstateT&& kstate,
        const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> kstate_range_streamer_settings);

// ***********************************************************************

// Template param `KstateT` should be `const T&`, `T&` or `T`.
// creating objects by means of `make_range_streamer` factory function ensure this.
template<typename _KstateT>
class KstateRangeStreamer {
    static_assert(! std::is_rvalue_reference<_KstateT>::value, "KstateT must not be a rvalue reference.");
    using KstateT = _KstateT;
    using SiteType = typename remove_cvref_t<KstateT>::SiteType;
    // static_assert: KstateT is subclass of KstateT<>;
public: // Factory function:
    friend
    KstateRangeStreamer<KstateT>
    make_kstate_range_streamer<KstateT>(
            KstateT&& kstate,
            const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> kstate_range_streamer_settings);
public: // API:
    std::ostream& stream(std::ostream& os) const;
    std::string str() const;
private:
    KstateRangeStreamer(
            std::add_rvalue_reference_t<KstateT> kstate,
            const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> kstate_range_streamer_settings);
private:
    const KstateT _kstate;
    const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> _kstate_range_streamer_settings;
};

// ***********************************************************************

template<typename _KstateT>
KstateRangeStreamer<_KstateT>::KstateRangeStreamer(
        std::add_rvalue_reference_t<KstateT> kstate,
        const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> kstate_range_streamer_settings) :
    _kstate(std::forward<KstateT>(kstate)),
    _kstate_range_streamer_settings(kstate_range_streamer_settings) {
}

// ***********************************************************************

template<typename _KstateT>
std::ostream&
KstateRangeStreamer<_KstateT>::stream(std::ostream& os) const{
    // Defaults:
    const auto default_stream_preparer = [](std::ostream& s) { s << "𝕂𝕤𝕥𝕒𝕥𝕖⦃"; };
    const auto default_stream_sustainer = [](std::ostream&, size_t) {};
    const auto default_stream_value_putter = [](std::ostream& s, SiteType t) { s << t; };
    const auto default_stream_separer = [](std::ostream& s) { s << "∙"; };
    const auto default_stream_finisher = [](std::ostream& s) { s << "⦄"; };
    const bool default_format_independence_flag = true;
    // Apply overrules:
    const auto stream_preparer = (
                _kstate_range_streamer_settings._range_streamer_settings._stream_preparer ?
                    *_kstate_range_streamer_settings._range_streamer_settings._stream_preparer :
                    default_stream_preparer);
    const auto stream_sustainer = (
                _kstate_range_streamer_settings._range_streamer_settings._stream_sustainer ?
                    *_kstate_range_streamer_settings._range_streamer_settings._stream_sustainer :
                    default_stream_sustainer);
    const auto stream_value_putter = (
                _kstate_range_streamer_settings._range_streamer_settings._stream_value_putter ?
                    *_kstate_range_streamer_settings._range_streamer_settings._stream_value_putter :
                    default_stream_value_putter);
    const auto stream_separer = (
                _kstate_range_streamer_settings._range_streamer_settings._stream_separer ?
                    *_kstate_range_streamer_settings._range_streamer_settings._stream_separer :
                    default_stream_separer);
    const auto stream_finisher = (
                _kstate_range_streamer_settings._range_streamer_settings._stream_finisher ?
                    *_kstate_range_streamer_settings._range_streamer_settings._stream_finisher :
                    default_stream_finisher);
    const auto format_independence_flag = (
                _kstate_range_streamer_settings._range_streamer_settings._format_independence_flag ?
                    *_kstate_range_streamer_settings._range_streamer_settings._format_independence_flag :
                    default_format_independence_flag);
    // Stream:
    extension::boost::stream_range_impl<typename remove_cvref_t<KstateT>::ConstAnyRangeType>(
                _kstate.to_any_range(),
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
KstateRangeStreamer<_KstateT>::str() const {
    std::ostringstream oss;
    stream(oss);
    return oss.str();
}

// ***********************************************************************

template<typename KstateT>
KstateRangeStreamer<KstateT> make_kstate_range_streamer(
        KstateT&& kstate,
        const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> kstate_range_streamer_settings) {
    return KstateRangeStreamer<KstateT>(std::forward<KstateT>(kstate), kstate_range_streamer_settings);
}

}  // namespace kstate

// #######################################################################
// ##  pragma: operator|, operator||, operator<<                        ##
// #######################################################################

namespace kstate::pramga {

template<typename KstateT>
KstateRangeStreamer<KstateT>
operator|(
        KstateT&& kstate,
        KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> kstate_range_streamer_settings) {
    return make_kstate_range_streamer<KstateT>(std::forward<KstateT>(kstate), kstate_range_streamer_settings);
}

template<typename KstateT>
KstateRangeStreamer<KstateT>
operator||(
        KstateT&& kstate,
        extension::boost::RangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> range_streamer_settings) {
    const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> kstate_range_streamer_settings{range_streamer_settings};
    return make_kstate_range_streamer<KstateT>(std::forward<KstateT>(kstate), kstate_range_streamer_settings);
}

template<typename KstateT>
std::ostream& operator<<(
        std::ostream& os,
        const KstateRangeStreamer<KstateT>& kstate_range_streamer) {
    kstate_range_streamer.stream(os);
    return os;
}

}  // namespace kstate::pramga

#endif
