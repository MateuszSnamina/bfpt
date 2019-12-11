#ifndef EXTENSION_RANGE_STREAMER_HPP
#define EXTENSION_RANGE_STREAMER_HPP

#include <boost/range.hpp>
#include <boost/range/adaptor/indexed.hpp>
// #include <boost/range/adaptor/sliced.hpp>

#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

// ########################################################
// ##  RangeStreamer                                     ##
// ########################################################

class RangeStreamer {
 public:
  RangeStreamer(std::ostream& os);
  RangeStreamer& set_stream_preparer(std::function<void(std::ostream&)>);
  RangeStreamer& set_stream_sustainer(
      std::function<void(std::ostream&, size_t)>);
  RangeStreamer& set_stream_finisher(std::function<void(std::ostream&)>);
  RangeStreamer& set_stream_separer(std::function<void(std::ostream&)>);
  template <typename SinglePassRange>
  RangeStreamer& stream(const SinglePassRange& rng);
  std::ostream& ostream();

 private:
  std::ostream& _os;
  std::function<void(std::ostream&)> _stream_preparer = [](std::ostream& s) {
    s << "{";
  };
  std::function<void(std::ostream&, size_t)> _stream_sustainer =
      [](std::ostream& s, size_t i) { s << i << ":"; };
  std::function<void(std::ostream&)> _stream_separer = [](std::ostream& s) {
    s << ",";
  };
  std::function<void(std::ostream&)> _stream_finisher = [](std::ostream& s) {
    s << "}";
  };
};

inline RangeStreamer::RangeStreamer(std::ostream& os) : _os(os) {}

inline RangeStreamer& RangeStreamer::set_stream_preparer(
    std::function<void(std::ostream&)> _) {
  _stream_preparer = _;
  return *this;
}

inline RangeStreamer& RangeStreamer::set_stream_sustainer(
    std::function<void(std::ostream&, size_t)> _) {
  _stream_sustainer = _;
  return *this;
}

inline RangeStreamer& RangeStreamer::set_stream_finisher(
    std::function<void(std::ostream&)> _) {
  _stream_finisher = _;
  return *this;
}

inline RangeStreamer& RangeStreamer::set_stream_separer(
    std::function<void(std::ostream&)> _) {
  _stream_separer = _;
  return *this;
}

template <typename SinglePassRange>
RangeStreamer& RangeStreamer::stream(const SinglePassRange& rng) {
  const auto d = std::distance(std::begin(rng), std::end(rng));
  assert(d >= 0);
  _stream_preparer(_os);
  for (const auto& _ : rng | boost::adaptors::indexed(0)) {
    _stream_sustainer(_os, _.index());
    _os << _.value();
    if (d > 1 && _.index() != d - 1) {
      _stream_separer(_os);
    }
  }
  _stream_finisher(_os);
  return *this;
}

inline std::ostream& RangeStreamer::ostream() { return _os; }

// ########################################################
// ##  RangeStreamStreamer                               ##
// ########################################################

class RangeStreamStreamer {
 public:
  RangeStreamStreamer();
  RangeStreamStreamer& set_stream_preparer(std::function<void(std::ostream&)>);
  RangeStreamStreamer& set_stream_sustainer(
      std::function<void(std::ostream&, size_t)>);
  RangeStreamStreamer& set_stream_finisher(std::function<void(std::ostream&)>);
  RangeStreamStreamer& set_stream_separer(std::function<void(std::ostream&)>);
  template <typename SinglePassRange>
  RangeStreamStreamer& stream(const SinglePassRange& rng);
  std::string str() const;

 private:
  std::ostringstream _oss;
  RangeStreamer _rs;
};

inline RangeStreamStreamer::RangeStreamStreamer() : _rs(_oss) {}

inline RangeStreamStreamer& RangeStreamStreamer::set_stream_preparer(
    std::function<void(std::ostream&)> _) {
  _rs.set_stream_preparer(_);
  return *this;
}

inline RangeStreamStreamer& RangeStreamStreamer::set_stream_sustainer(
    std::function<void(std::ostream&, size_t)> _) {
  _rs.set_stream_sustainer(_);
  return *this;
}

inline RangeStreamStreamer& RangeStreamStreamer::set_stream_finisher(
    std::function<void(std::ostream&)> _) {
  _rs.set_stream_finisher(_);
  return *this;
}

inline RangeStreamStreamer& RangeStreamStreamer::set_stream_separer(
    std::function<void(std::ostream&)> _) {
  _rs.set_stream_separer(_);
  return *this;
}

template <typename SinglePassRange>
RangeStreamStreamer& RangeStreamStreamer::stream(const SinglePassRange& rng) {
  _rs.stream(rng);
  return *this;
}

inline std::string RangeStreamStreamer::str() const { return _oss.str(); }

#endif