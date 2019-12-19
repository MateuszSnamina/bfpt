#ifndef KSTATE_KSTATE_STREAMER_HPP
#define KSTATE_KSTATE_STREAMER_HPP

#include <extensions/range_streamer.hpp>
#include <kstate/kstate.hpp>

#include <boost/range.hpp>
#include <boost/range/adaptor/indexed.hpp>

#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

// #######################################################################
// ##  KstateStreamer                                                   ##
// #######################################################################

namespace kstate {

static RangeStreamerSettings kstate_default_range_streamer_settings{
    [](std::ostream& s) { s << "𝕂𝕊⦃"; },
    [](std::ostream& s, size_t i) {},
    [](std::ostream& s) { s << "∙"; },
    [](std::ostream& s) { s << "⦄"; }};

class KstateStreamer {
   public:
    KstateStreamer(std::ostream& os, RangeStreamerSettings range_streamer_settings = kstate_default_range_streamer_settings);
    KstateStreamer& set_range_streamer_settings(RangeStreamerSettings);
    KstateStreamer& set_stream_preparer(std::function<void(std::ostream&)>);
    KstateStreamer& set_stream_sustainer(std::function<void(std::ostream&, size_t)>);
    KstateStreamer& set_stream_finisher(std::function<void(std::ostream&)>);
    KstateStreamer& set_stream_separer(std::function<void(std::ostream&)>);
    template <typename SiteType>
    KstateStreamer& stream(const Kstate<SiteType>& rng);
    std::ostream& ostream();
    const std::ostream& ostream() const;

   private:
    // The doublestruck font: https://en.wikipedia.org/wiki/Blackboard_bold
    std::ostream& _os;
    RangeStreamerSettings _range_streamer_settings = kstate_default_range_streamer_settings;
};

// ***********************************************************************

inline KstateStreamer::KstateStreamer(std::ostream& os, RangeStreamerSettings range_streamer_settings)
    : _os(os), _range_streamer_settings(range_streamer_settings) {}

// ***********************************************************************

inline KstateStreamer& KstateStreamer::set_range_streamer_settings(RangeStreamerSettings _) {
    _range_streamer_settings = _;
    return *this;
}

inline KstateStreamer& KstateStreamer::set_stream_preparer(std::function<void(std::ostream&)> _) {
    _range_streamer_settings.set_stream_preparer(_);
    return *this;
}

inline KstateStreamer& KstateStreamer::set_stream_sustainer(std::function<void(std::ostream&, size_t)> _) {
    _range_streamer_settings.set_stream_sustainer(_);
    return *this;
}

inline KstateStreamer& KstateStreamer::set_stream_finisher(std::function<void(std::ostream&)> _) {
    _range_streamer_settings.set_stream_finisher(_);
    return *this;
}

inline KstateStreamer& KstateStreamer::set_stream_separer(std::function<void(std::ostream&)> _) {
    _range_streamer_settings.set_stream_separer(_);
    return *this;
}

template <typename SiteType>
KstateStreamer& KstateStreamer::stream(const Kstate<SiteType>& kstate) {
    extension::boost::RangeStreamer(_os)
        .set_stream_preparer(_range_streamer_settings._stream_preparer)
        .set_stream_sustainer(_range_streamer_settings._stream_sustainer)
        .set_stream_separer(_range_streamer_settings._stream_separer)
        .set_stream_finisher(_range_streamer_settings._stream_finisher)  //TODO make is a single set_range_streamer_settings call.
        .stream(kstate.to_any_range());
    return *this;
}

inline std::ostream& KstateStreamer::ostream() {
    return _os;
}

inline const std::ostream& KstateStreamer::ostream() const {
    return _os;
}

// #######################################################################
// ##  KstateStringStreamer                                             ##
// #######################################################################

class KstateStringStreamer {
   public:
    KstateStringStreamer(RangeStreamerSettings range_streamer_settings = kstate_default_range_streamer_settings);
    KstateStringStreamer& set_stream_preparer(std::function<void(std::ostream&)>);
    KstateStringStreamer& set_stream_sustainer(std::function<void(std::ostream&, size_t)>);
    KstateStringStreamer& set_stream_finisher(std::function<void(std::ostream&)>);
    KstateStringStreamer& set_stream_separer(std::function<void(std::ostream&)>);
    template <typename SiteType>
    KstateStringStreamer& stream(const Kstate<SiteType>& rng);
    std::string str() const;

   private:
    std::ostringstream _oss;
    KstateStreamer _rs;
};

// ***********************************************************************

inline KstateStringStreamer::KstateStringStreamer(RangeStreamerSettings range_streamer_settings)
    : _rs(_oss, range_streamer_settings) {
}

// ***********************************************************************

inline KstateStringStreamer& KstateStringStreamer::set_stream_preparer(std::function<void(std::ostream&)> _) {
    _rs.set_stream_preparer(_);
    return *this;
}

inline KstateStringStreamer& KstateStringStreamer::set_stream_sustainer(std::function<void(std::ostream&, size_t)> _) {
    _rs.set_stream_sustainer(_);
    return *this;
}

inline KstateStringStreamer& KstateStringStreamer::set_stream_finisher(std::function<void(std::ostream&)> _) {
    _rs.set_stream_finisher(_);
    return *this;
}

inline KstateStringStreamer& KstateStringStreamer::set_stream_separer(std::function<void(std::ostream&)> _) {
    _rs.set_stream_separer(_);
    return *this;
}

template <typename SiteType>
KstateStringStreamer& KstateStringStreamer::stream(const Kstate<SiteType>& kstate) {
    _rs.stream(kstate);
    return *this;
}

inline std::string KstateStringStreamer::str() const {
    return _oss.str();
}

}  // namespace kstate

#endif