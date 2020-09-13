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
// ##  RangeStreamer                                                    ##
// #######################################################################

namespace kstate {

// ***********************************************************************

template<typename _KstateT>
class KstateStreamer;

// ***********************************************************************

template <typename KstateT>
KstateStreamer<KstateT>
make_kstate_streamer(
        KstateT&&,
        const extension::boost::RangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType>);

// ***********************************************************************

template<typename _KstateT>
class KstateStreamer {
    static_assert(! ::std::is_rvalue_reference_v<_KstateT>,
    "KstateT must not be a rvalue reference.");
    static_assert(::std::is_lvalue_reference_v<_KstateT> || (!::std::is_reference_v<_KstateT> && !::std::is_const_v<_KstateT>),
    "KstateT must be of form: `T`, `T&` or `const T&`.");
    // static_assert: KstateT is subclass of KstateT<>; //TODO
public: // Helper types:
    using KstateT = _KstateT;
    using SiteType = typename remove_cvref_t<KstateT>::SiteType;
public: // Factory function:
    friend
    KstateStreamer<KstateT>
    make_kstate_streamer<KstateT>(KstateT&&, const extension::boost::RangeStreamerSettings<SiteType>);
public: // API:
    std::ostream& stream(std::ostream& os) const;
    std::string str() const;
private:
    KstateStreamer(std::add_rvalue_reference_t<KstateT>, const extension::boost::RangeStreamerSettings<SiteType>);
private:
    const KstateT _kstate;
    const extension::boost::RangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> _range_streamer_settings;
};

// ***********************************************************************

template<typename _KstateT>
KstateStreamer<_KstateT>::KstateStreamer(
        std::add_rvalue_reference_t<KstateT> kstate,
        const extension::boost::RangeStreamerSettings<SiteType> range_streamer_settings) :
    _kstate(std::forward<KstateT>(kstate)),
    _range_streamer_settings(range_streamer_settings) {
}

// ***********************************************************************

template<typename _KstateT>
std::ostream&
KstateStreamer<_KstateT>::stream(std::ostream& os) const{
    // Defaults:
    const auto default_stream_preparer = [](std::ostream& s) { s << "𝕂𝕤𝕥𝕒𝕥𝕖⦃"; };
    const auto default_stream_sustainer = [](std::ostream&, size_t) {};
    const auto default_stream_value_putter = [](std::ostream& s, SiteType t) { s << t; };
    const auto default_stream_separer = [](std::ostream& s) { s << "∙"; };
    const auto default_stream_finisher = [](std::ostream& s) { s << "⦄"; };
    const bool default_format_independence_flag = true;
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
KstateStreamer<_KstateT>::str() const {
    std::ostringstream oss;
    stream(oss);
    return oss.str();
}

// ***********************************************************************

template<typename KstateT>
KstateStreamer<KstateT> make_kstate_streamer(
        KstateT&& kstate,
        const extension::boost::RangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> range_streamer_settings) {
    return KstateStreamer<KstateT>(std::forward<KstateT>(kstate), range_streamer_settings);
}

}  // namespace kstate

// #######################################################################
// ##  pragma: operator||, operator<<                                   ##
// #######################################################################

namespace kstate::pramga {

template<typename KstateT>
KstateStreamer<KstateT>
operator||(
        KstateT&& kstate,
        const extension::boost::RangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteType> range_streamer_settings) {
    return make_kstate_streamer<KstateT>(std::forward<KstateT>(kstate), range_streamer_settings);
}

template<typename KstateT>
std::ostream& operator<<(
        std::ostream& os,
        const KstateStreamer<KstateT>& kstate_streamer) {
    kstate_streamer.stream(os);
    return os;
}

}  // namespace kstate::pramga

#endif
