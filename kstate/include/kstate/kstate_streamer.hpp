#ifndef KSTATE_KSTATE_STREAMER_HPP
#define KSTATE_KSTATE_STREAMER_HPP

#include <extensions/range_streamer.hpp>
#include <kstate/kstate_abstract.hpp>

#include <boost/range.hpp>
#include <boost/range/adaptor/indexed.hpp>

#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

// The doublestruck font: https://en.wikipedia.org/wiki/Blackboard_bold

// #######################################################################
// ##  remove_cvref (C++20)                                             ##
// #######################################################################

namespace kstate {

template< class T >
struct remove_cvref {
    typedef std::remove_cv_t<std::remove_reference_t<T>> type;
};

template< class T >
using remove_cvref_t = typename remove_cvref<T>::type;

}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@  VARIANT 1  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// #######################################################################
// ##  KstateRangeStreamerSettings                                      ##
// #######################################################################

namespace kstate {

template<typename T>
struct KstateRangeStreamerSettings {
    KstateRangeStreamerSettings(extension::boost::RangeStreamerSettings<T> range_streamer_settings) :
        _range_streamer_settings(range_streamer_settings) {
    }
    const extension::boost::RangeStreamerSettings<T> _range_streamer_settings;
};

}

// #######################################################################
// ##  RangeStreamer                                                    ##
// #######################################################################

namespace kstate {

// ***********************************************************************

template <typename KstateT>
class KstateRangeStreamer;

// ***********************************************************************

template <typename KstateT>
KstateRangeStreamer<KstateT>
make_kstate_range_streamer(
        KstateT&& kstate,
        const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteTypeXXX> kstate_range_streamer_settings);

// ***********************************************************************

// Template param `KstateT` should be `const T&`, `T&` or `T`.
// creating objects by means of `make_range_streamer` factory function ensure this.
template <typename KstateT>
class KstateRangeStreamer {
    static_assert(! ::std::is_rvalue_reference<KstateT>::value, "KstateT must not be a rvalue reference.");
    // static_assert: KstateT is subclass of KstateT<>;
private:
    KstateRangeStreamer(
            ::std::add_rvalue_reference_t<KstateT> kstate,
            const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteTypeXXX> kstate_range_streamer_settings) :
        _kstate(::std::forward<KstateT>(kstate)),
        _kstate_range_streamer_settings(kstate_range_streamer_settings) {
    }
public:
    friend
    KstateRangeStreamer<KstateT>
    make_kstate_range_streamer<KstateT>(
            KstateT&& kstate,
            const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteTypeXXX> kstate_range_streamer_settings);
    ::std::ostream& stream(::std::ostream& os) const {
        // Defaults:
        const ::std::function<void(::std::ostream&)> default_stream_preparer =
                [](std::ostream& s) { s << "𝕂𝕤𝕥𝕒𝕥𝕖⦃"; };
        const ::std::function<void(::std::ostream&, size_t)> default_stream_sustainer =
                [](std::ostream&, size_t) {};
        const ::std::function<void(::std::ostream&, typename remove_cvref_t<KstateT>::SiteTypeXXX)> default_stream_value_putter =
                [](::std::ostream& s, typename remove_cvref_t<KstateT>::SiteTypeXXX t) { s << t; };
        const ::std::function<void(::std::ostream&)> default_stream_separer =
                [](std::ostream& s) { s << "∙"; };
        const ::std::function<void(::std::ostream&)> default_stream_finisher =
                [](std::ostream& s) { s << "⦄"; };
        bool default_format_independence_flag = true;
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
    ::std::string str() const {
        ::std::ostringstream oss;
        stream(oss);
        return oss.str();
    }
    const KstateT _kstate;
    const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteTypeXXX> _kstate_range_streamer_settings;
};

// ***********************************************************************
template <typename KstateT>
KstateRangeStreamer<KstateT> make_kstate_range_streamer(
    KstateT&& kstate,
    const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteTypeXXX> kstate_range_streamer_settings) {
    return KstateRangeStreamer<KstateT>(::std::forward<KstateT>(kstate), kstate_range_streamer_settings);
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
        KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteTypeXXX> kstate_range_streamer_settings) {
    return make_kstate_range_streamer<KstateT>(std::forward<KstateT>(kstate), kstate_range_streamer_settings);
}

template<typename KstateT>
KstateRangeStreamer<KstateT>
operator||(
        KstateT&& kstate,
        extension::boost::RangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteTypeXXX> range_streamer_settings) {
    const KstateRangeStreamerSettings<typename remove_cvref_t<KstateT>::SiteTypeXXX> kstate_range_streamer_settings{range_streamer_settings};
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

//// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//// @@@@@@@@@@@@@@@@@  VARIANT 2  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//// #######################################################################
//// ##  make_range_streamer_from_kstate, operator||                      ##
//// #######################################################################

//namespace kstate {

//template <typename KstateT>
//extension::boost::RangeStreamer<typename KstateT::ConstAnyRangeType >
//make_range_streamer_from_kstate(
//    const KstateT& kstate,
//    extension::boost::RangeStreamerSettings range_streamer_settings) {
//    // Defaults:
//    const ::std::function<void(::std::ostream&)> default_stream_preparer =
//            [](std::ostream& s) { s << "𝕂𝕤𝕥𝕒𝕥𝕖⦃"; };
//    const ::std::function<void(::std::ostream&, size_t)> default_stream_sustainer =
//            [](std::ostream&, size_t) {};
//    const ::std::function<void(::std::ostream&)> default_stream_separer =
//            [](std::ostream& s) { s << "∙"; };
//    const ::std::function<void(::std::ostream&)> default_stream_finisher =
//            [](std::ostream& s) { s << "⦄"; };
//    bool default_format_independence_flag = true;
//    // Apply overrules:
//    const auto stream_preparer = (
//                range_streamer_settings._stream_preparer ?
//                *range_streamer_settings._stream_preparer :
//                default_stream_preparer);
//    const auto stream_sustainer = (
//                range_streamer_settings._stream_sustainer ?
//                *range_streamer_settings._stream_sustainer :
//                    default_stream_sustainer);
//    const auto stream_separer = (
//                range_streamer_settings._stream_separer ?
//                *range_streamer_settings._stream_separer :
//                default_stream_separer);
//    const auto stream_finisher = (
//                range_streamer_settings._stream_finisher ?
//                *range_streamer_settings._stream_finisher :
//                default_stream_finisher);
//    const auto format_independence_flag = (
//                range_streamer_settings._format_independence_flag ?
//                *range_streamer_settings._format_independence_flag :
//                default_format_independence_flag);
//    // Build range_streamer_settings:
//    const auto range_streamer_settings_override = extension::boost::RangeStreamerSettings()
//            .set_stream_preparer(stream_preparer)
//            .set_stream_sustainer(stream_sustainer)
//            .set_stream_separer(stream_separer)
//            .set_stream_finisher(stream_finisher)
//            .set_format_independence_flag(format_independence_flag);
//    // Build range_streamer:
//    const auto range_streamer = extension::boost::make_range_streamer(
//                kstate.to_any_range(),
//                range_streamer_settings_override);
//    return range_streamer;
//}

//}  // namespace kstate

//namespace kstate::pramga {

//template<typename KstateT>
//extension::boost::RangeStreamer<typename KstateT::ConstAnyRangeType >
//operator|| (
//        const KstateT& kstate,
//        extension::boost::RangeStreamerSettings range_streamer_settings) {
//    return make_range_streamer_from_kstate<KstateT>(kstate, range_streamer_settings);
//}

//}  // namespace kstate::pramga

#endif
