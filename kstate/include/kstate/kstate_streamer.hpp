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

namespace kstate {

static extension::boost::RangeStreamerSettings kstate_default_range_streamer_settings{
    [](std::ostream& s) { s << "𝕂𝕤𝕥𝕒𝕥𝕖⦃"; },
    [](std::ostream& s, size_t i) {},
    [](std::ostream& s) { s << "∙"; },
    [](std::ostream& s) { s << "⦄"; }};

}

//// #######################################################################
//// ##  KstateStreamer                                                   ##
//// #######################################################################

//namespace kstate {

//// The doublestruck font: https://en.wikipedia.org/wiki/Blackboard_bold

//static extension::boost::RangeStreamerSettings kstate_default_range_streamer_settings{
//    [](std::ostream& s) { s << "𝕂𝕤𝕥𝕒𝕥𝕖⦃"; },
//    [](std::ostream& s, size_t i) {},
//    [](std::ostream& s) { s << "∙"; },
//    [](std::ostream& s) { s << "⦄"; }};

//class KstateStreamer {
//public:  // Ctor:
//    KstateStreamer(std::ostream& os, extension::boost::RangeStreamerSettings range_streamer_settings = kstate_default_range_streamer_settings);

//public:  // Internal ostream assessor:
//    std::ostream& ostream();
//    const std::ostream& ostream() const;

//public:  // Setters for fine streaming settings:
//    KstateStreamer& set_range_streamer_settings(extension::boost::RangeStreamerSettings);
//    KstateStreamer& set_stream_preparer(std::function<void(std::ostream&)>);
//    KstateStreamer& set_stream_sustainer(std::function<void(std::ostream&, size_t)>);
//    KstateStreamer& set_stream_separer(std::function<void(std::ostream&)>);
//    KstateStreamer& set_stream_finisher(std::function<void(std::ostream&)>);
//    KstateStreamer& set_format_independence_flag(bool = true);

//public:  // The streaming function:
//    template <typename SiteType>
//    KstateStreamer& stream(const Kstate<SiteType>& rng);

//private:
//    std::ostream& _os;
//    extension::boost::RangeStreamerSettings _range_streamer_settings = kstate_default_range_streamer_settings;
//};

//// ***********************************************************************

//inline KstateStreamer::KstateStreamer(std::ostream& os, extension::boost::RangeStreamerSettings range_streamer_settings)
//    : _os(os), _range_streamer_settings(range_streamer_settings) {}

//// ***********************************************************************

//inline KstateStreamer& KstateStreamer::set_range_streamer_settings(extension::boost::RangeStreamerSettings _) {
//    _range_streamer_settings = _;
//    return *this;
//}

//inline KstateStreamer& KstateStreamer::set_stream_preparer(std::function<void(std::ostream&)> _) {
//    _range_streamer_settings.set_stream_preparer(_);
//    return *this;
//}

//inline KstateStreamer& KstateStreamer::set_stream_sustainer(std::function<void(std::ostream&, size_t)> _) {
//    _range_streamer_settings.set_stream_sustainer(_);
//    return *this;
//}

//inline KstateStreamer& KstateStreamer::set_stream_separer(std::function<void(std::ostream&)> _) {
//    _range_streamer_settings.set_stream_separer(_);
//    return *this;
//}

//inline KstateStreamer& KstateStreamer::set_stream_finisher(std::function<void(std::ostream&)> _) {
//    _range_streamer_settings.set_stream_finisher(_);
//    return *this;
//}

//inline KstateStreamer& KstateStreamer::set_format_independence_flag(bool _) {
//    _range_streamer_settings.set_format_independence_flag(_);
//    return *this;
//}

//template <typename SiteType>
//KstateStreamer& KstateStreamer::stream(const Kstate<SiteType>& kstate) {
//    using namespace extension::boost::stream_pragma;
//    const auto range_stream_settings = RSS()
//            .set_stream_preparer(*_range_streamer_settings._stream_preparer)
//            .set_stream_sustainer(*_range_streamer_settings._stream_sustainer)
//            .set_stream_separer(*_range_streamer_settings._stream_separer)
//            .set_stream_finisher(*_range_streamer_settings._stream_finisher); //TODO look at it!
//    _os << (kstate.to_any_range() | range_stream_settings);
//    return *this;
//}

//inline std::ostream& KstateStreamer::ostream() {
//    return _os;
//}

//inline const std::ostream& KstateStreamer::ostream() const {
//    return _os;
//}

//// #######################################################################
//// ##  KstateStringStreamer                                             ##
//// #######################################################################

//class KstateStringStreamer {
//public:  // Ctor:
//    KstateStringStreamer(extension::boost::RangeStreamerSettings range_streamer_settings = kstate_default_range_streamer_settings);

//public:  // Setters for fine streaming settings:
//    KstateStringStreamer& set_stream_preparer(std::function<void(std::ostream&)>);
//    KstateStringStreamer& set_stream_sustainer(std::function<void(std::ostream&, size_t)>);
//    KstateStringStreamer& set_stream_separer(std::function<void(std::ostream&)>);
//    KstateStringStreamer& set_stream_finisher(std::function<void(std::ostream&)>);
//    KstateStringStreamer& set_format_independence_flag(bool = true);

//public:  // The streaming function:
//    template <typename SiteType>
//    KstateStringStreamer& stream(const Kstate<SiteType>& rng);

//public:  // Function to retreive the streaming result
//    std::string str() const;

//private:
//    std::ostringstream _oss;
//    KstateStreamer _rs;
//};

//// ***********************************************************************

//inline KstateStringStreamer::KstateStringStreamer(extension::boost::RangeStreamerSettings range_streamer_settings)
//    : _rs(_oss, range_streamer_settings) {
//}

//// ***********************************************************************

//inline KstateStringStreamer& KstateStringStreamer::set_stream_preparer(std::function<void(std::ostream&)> _) {
//    _rs.set_stream_preparer(_);
//    return *this;
//}

//inline KstateStringStreamer& KstateStringStreamer::set_stream_sustainer(std::function<void(std::ostream&, size_t)> _) {
//    _rs.set_stream_sustainer(_);
//    return *this;
//}

//inline KstateStringStreamer& KstateStringStreamer::set_stream_separer(std::function<void(std::ostream&)> _) {
//    _rs.set_stream_separer(_);
//    return *this;
//}

//inline KstateStringStreamer& KstateStringStreamer::set_stream_finisher(std::function<void(std::ostream&)> _) {
//    _rs.set_stream_finisher(_);
//    return *this;
//}

//inline KstateStringStreamer& KstateStringStreamer::set_format_independence_flag(bool _) {
//    _rs.set_format_independence_flag(_);
//    return *this;
//}

//template <typename SiteType>
//KstateStringStreamer& KstateStringStreamer::stream(const Kstate<SiteType>& kstate) {
//    _rs.stream(kstate);
//    return *this;
//}

//inline std::string KstateStringStreamer::str() const {
//    return _oss.str();
//}

//}  // namespace kstate

// #######################################################################
// ##  RangeStreamer                                                    ##
// #######################################################################

namespace kstate {

class KstateRangeStreamerSettings {
public:
    KstateRangeStreamerSettings(extension::boost::RangeStreamerSettings range_streamer_settings) :
        _range_streamer_settings(range_streamer_settings) {
    }
public: //private: //TODO make private!!!
    const extension::boost::RangeStreamerSettings _range_streamer_settings;
};

}

// #######################################################################
// ##  RangeStreamer                                                    ##
// #######################################################################

namespace kstate {

template <typename KstateT>
class KstateRangeStreamer;

template <typename KstateT>
KstateRangeStreamer<KstateT> make_kstate_range_streamer(
        KstateT&& kstate,
        const KstateRangeStreamerSettings kstate_range_streamer_settings);

// Template param `KstateT` should be `const T&`, `T&` or `T`.
// creating objects by means of `make_range_streamer` factory function ensure this.
template <typename KstateT>
class KstateRangeStreamer {
    static_assert(! ::std::is_rvalue_reference<KstateT>::value, "KstateT must not be a rvalue reference.");
    // static_assert: KstateT is subclass of KstateT<>;
private:
    KstateRangeStreamer(::std::add_rvalue_reference_t<KstateT> kstate, const KstateRangeStreamerSettings kstate_range_streamer_settings) :
        _kstate(::std::forward<KstateT>(kstate)),
        _kstate_range_streamer_settings(kstate_range_streamer_settings) {
    }
public:

    friend KstateRangeStreamer<KstateT> make_kstate_range_streamer<KstateT>(KstateT&& kstate, const KstateRangeStreamerSettings kstate_range_streamer_settings);

    ::std::ostream& stream(::std::ostream& os) const {
        // Defaults:
        const ::std::function<void(::std::ostream&)> default_stream_preparer =
                [](std::ostream& s) { s << "𝕂𝕤𝕥𝕒𝕥𝕖⦃"; };
        const ::std::function<void(::std::ostream&, size_t)> default_stream_sustainer =
                [](std::ostream&, size_t) {};
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
        // Build range_streamer_settings:
        const auto range_streamer_settings = extension::boost::RangeStreamerSettings()
                .set_stream_preparer(stream_preparer)
                .set_stream_sustainer(stream_sustainer)
                .set_stream_separer(stream_separer)
                .set_stream_finisher(stream_finisher)
                .set_format_independence_flag(format_independence_flag);
        // Stream:
        const auto range_streamer = make_range_streamer(_kstate.to_any_range(), range_streamer_settings);
        range_streamer.stream(os);
        return os;
    }

    ::std::string str() const {
        ::std::ostringstream oss;
        stream(oss);
        return oss.str();
    }

    const KstateT _kstate;
    const KstateRangeStreamerSettings _kstate_range_streamer_settings;
};  // end of class KstateRangeStreamer<KstateT>


template <typename KstateT>
KstateRangeStreamer<KstateT> make_kstate_range_streamer(
    KstateT&& kstate,
    const KstateRangeStreamerSettings kstate_range_streamer_settings) {
    return KstateRangeStreamer<KstateT>(::std::forward<KstateT>(kstate), kstate_range_streamer_settings);
}

}  // namespace kstate

namespace kstate::pramga {

template<typename KstateT>
KstateRangeStreamer<KstateT>
operator|(
        KstateT&& kstate,
        const KstateRangeStreamerSettings& kstate_range_streamer_settings) {
    return make_kstate_range_streamer<KstateT>(kstate, kstate_range_streamer_settings);
}

template<typename KstateT>
KstateRangeStreamer<KstateT>
operator|| (
        KstateT&& kstate,
        const extension::boost::RangeStreamerSettings& range_streamer_settings) {
    const KstateRangeStreamerSettings kstate_range_streamer_settings{range_streamer_settings};
    return make_kstate_range_streamer<KstateT>(kstate, kstate_range_streamer_settings);
}

template<typename KstateT>
std::ostream& operator<<(std::ostream& os, const KstateRangeStreamer<KstateT>& kstate_range_streamer) {
    kstate_range_streamer.stream(os);
    return os;
}

}  // namespace kstate::pramga

#endif
