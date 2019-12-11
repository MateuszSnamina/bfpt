#ifndef EXTENSIONS_BOOST_ADAPTORS_HPP
#define EXTENSIONS_BOOST_ADAPTORS_HPP

#include <boost/range/join.hpp>
#include <iterator>

#ifndef NDEBUG
#include <boost/numeric/conversion/cast.hpp>
#include <cassert>
#endif

namespace extension_implementation::boost::adaptors {

// ########################################################
// ##  rotated                                           ##
// ########################################################

class RotateHolder {
 public:
  RotateHolder(size_t n) : n(n){};
  size_t n;
};

/*
 * Args domain: h.n has to be a positive number, less than the range size.
 * The performed rotation is a left rotation.
 *
 */
template <typename ForwardRange>
inline ::boost::joined_range<
    const typename ::boost::iterator_range<
        typename ::boost::range_iterator<const ForwardRange>::type>,
    const typename ::boost::iterator_range<
        typename ::boost::range_iterator<const ForwardRange>::type>>
operator|(const ForwardRange &rng, const RotateHolder &h) {
#ifndef NDEBUG
  const auto d = std::distance(std::begin(rng), std::end(rng));
  assert(d >= 0);
  assert(::boost::numeric_cast<decltype(d)>(h.n) <= d);
  assert(h.n >= 0);
#endif
  const auto mid = std::next(std::begin(rng), h.n);
  return ::boost::join(::boost::make_iterator_range(mid, std::end(rng)),
                       ::boost::make_iterator_range(std::begin(rng), mid));
}

RotateHolder rotated(size_t n) { return RotateHolder(n); }

// ########################################################
// ##  doubled                                           ##
// ########################################################

class Doubler {};

template <typename SinglePassRange>
inline ::boost::joined_range<
    const typename ::boost::iterator_range<
        typename ::boost::range_iterator<const SinglePassRange>::type>,
    const typename ::boost::iterator_range<
        typename ::boost::range_iterator<const SinglePassRange>::type>>
operator|(const SinglePassRange &rng, const Doubler &) {
  return ::boost::join(
      ::boost::make_iterator_range(std::begin(rng), std::end(rng)),
      ::boost::make_iterator_range(std::begin(rng), std::end(rng)));
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static Doubler doubled{};
#pragma GCC diagnostic pop

}  // namespace extension_implementation::boost::adaptors

namespace extension::boost::adaptors {
using extension_implementation::boost::adaptors::doubled;
using extension_implementation::boost::adaptors::rotated;
}  // namespace extension::boost::adaptors

#endif