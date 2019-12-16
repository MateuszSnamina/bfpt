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
// ##  KstateStreamer                                                    ##
// #######################################################################

namespace kstate {

class KstateStreamer {
   public:
    KstateStreamer(std::ostream& os);
    KstateStreamer& set_stream_preparer(std::function<void(std::ostream&)>);
    KstateStreamer& set_stream_sustainer(std::function<void(std::ostream&, size_t)>);
    KstateStreamer& set_stream_finisher(std::function<void(std::ostream&)>);
    KstateStreamer& set_stream_separer(std::function<void(std::ostream&)>);
    template <typename SinglePassRange>
    KstateStreamer& stream(const Kstate<SinglePassRange>& rng);
    std::ostream& ostream();
    const std::ostream& ostream() const;

   private:
    // The doublestruck font: https://en.wikipedia.org/wiki/Blackboard_bold
    std::ostream& _os;
    std::function<void(std::ostream&)> _stream_preparer = [](std::ostream& s) {
        s << "𝕂𝕊⦃";
    };
    std::function<void(std::ostream&, size_t)> _stream_sustainer =
        [](std::ostream& s, size_t i) {};
    std::function<void(std::ostream&)> _stream_separer = [](std::ostream& s) {
        s << "∙";
    };
    std::function<void(std::ostream&)> _stream_finisher = [](std::ostream& s) {
        s << "⦄";
    };
};

// ***********************************************************************

inline KstateStreamer::KstateStreamer(std::ostream& os)
    : _os(os) {}

// ***********************************************************************

inline KstateStreamer& KstateStreamer::set_stream_preparer(std::function<void(std::ostream&)> _) {
    _stream_preparer = _;
    return *this;
}

inline KstateStreamer& KstateStreamer::set_stream_sustainer(std::function<void(std::ostream&, size_t)> _) {
    _stream_sustainer = _;
    return *this;
}

inline KstateStreamer& KstateStreamer::set_stream_finisher(std::function<void(std::ostream&)> _) {
    _stream_finisher = _;
    return *this;
}

inline KstateStreamer& KstateStreamer::set_stream_separer(std::function<void(std::ostream&)> _) {
    _stream_separer = _;
    return *this;
}

template <typename SinglePassRange>
KstateStreamer& KstateStreamer::stream(const Kstate<SinglePassRange>& kstate) {
    extension::boost::RangeStreamer(_os)
        .set_stream_preparer(_stream_preparer)
        .set_stream_sustainer(_stream_sustainer)
        .set_stream_separer(_stream_separer)
        .set_stream_finisher(_stream_finisher)
        .stream(kstate.to_range());
    return *this;
}

inline std::ostream& KstateStreamer::ostream() {
    return _os;
}

inline const std::ostream& KstateStreamer::ostream() const {
    return _os;
}

// #######################################################################
// ##  KstateStringStreamer                                              ##
// #######################################################################

class KstateStringStreamer {
   public:
    KstateStringStreamer();
    KstateStringStreamer& set_stream_preparer(std::function<void(std::ostream&)>);
    KstateStringStreamer& set_stream_sustainer(std::function<void(std::ostream&, size_t)>);
    KstateStringStreamer& set_stream_finisher(std::function<void(std::ostream&)>);
    KstateStringStreamer& set_stream_separer(std::function<void(std::ostream&)>);
    template <typename SinglePassRange>
    KstateStringStreamer& stream(const Kstate<SinglePassRange>& rng);
    std::string str() const;

   private:
    std::ostringstream _oss;
    KstateStreamer _rs;
};

// ***********************************************************************

inline KstateStringStreamer::KstateStringStreamer()
    : _rs(_oss) {
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

template <typename SinglePassRange>
KstateStringStreamer& KstateStringStreamer::stream(const Kstate<SinglePassRange>& kstate) {
    _rs.stream(kstate);
    return *this;
}

inline std::string KstateStringStreamer::str() const {
    return _oss.str();
}

}  // namespace kstate

#endif