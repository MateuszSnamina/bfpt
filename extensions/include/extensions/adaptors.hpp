#ifndef EXTENSIONS_BOOST_ADAPTORS_HPP
#define EXTENSIONS_BOOST_ADAPTORS_HPP

#include <boost/range/join.hpp>
#include <iterator>

#ifndef NDEBUG
#include <boost/numeric/conversion/cast.hpp>
#include <cassert>
#endif

namespace extension_implementation::boost::adaptors {

class RotateHolder {
 public:
  RotateHolder(size_t n) : n(n){};
  size_t n;
};

template <typename SinglePassRange>
inline ::boost::joined_range<
    const typename ::boost::iterator_range<
        typename ::boost::range_iterator<SinglePassRange>::type>,
    const typename ::boost::iterator_range<
        typename ::boost::range_iterator<SinglePassRange>::type>>
operator|(SinglePassRange &rng, const RotateHolder &f) {
  const auto d = std::distance(std::begin(rng), std::end(rng));
  assert(d >= 0);
  assert(::boost::numeric_cast<decltype(d)>(f.n) <= d);
  const auto mid = std::next(std::begin(rng), f.n);
  return ::boost::join(::boost::make_iterator_range(mid, std::end(rng)),
                       ::boost::make_iterator_range(std::begin(rng), mid));
}

RotateHolder rotated(size_t n) { return RotateHolder(n); }

class Doubler {};

template <typename SinglePassRange>
inline ::boost::joined_range<
    const typename ::boost::iterator_range<
        typename ::boost::range_iterator<SinglePassRange>::type>,
    const typename ::boost::iterator_range<
        typename ::boost::range_iterator<SinglePassRange>::type>>
operator|(SinglePassRange &rng, const Doubler &) {
  return ::boost::join(
      ::boost::make_iterator_range(std::begin(rng), std::end(rng)),
      ::boost::make_iterator_range(std::begin(rng), std::end(rng)));
}

static Doubler doubled{};

}  // namespace extension_implementation::boost::adaptors

namespace extension::boost::adaptors {
using extension_implementation::boost::adaptors::doubled;
using extension_implementation::boost::adaptors::rotated;
}  // namespace extension::boost::adaptors

#endif