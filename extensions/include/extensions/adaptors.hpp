#ifndef EXTENSIONS_BOOST_ADAPTORS_HPP
#define EXTENSIONS_BOOST_ADAPTORS_HPP

#include <boost/range/join.hpp>

#include <iterator>

#ifndef NDEBUG
#include <boost/numeric/conversion/cast.hpp>
#include <cassert>
#endif

// #######################################################################
// ##  rotated                                                          ##
// #######################################################################

namespace extension_implementation::boost::adaptors {

class RotateHolder {
   public:
    RotateHolder(size_t n) : n(n){};
    size_t n;
};

inline RotateHolder rotated(size_t n) {
    return RotateHolder(n);
}

/*
 * Args domain: h.n has to be a positive number, less than the range size.
 * The performed rotation is a left rotation.
 */

template <typename ForwardRange>
using RotatedRangeType = ::boost::joined_range<
    const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>,
    const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>>;

template <typename ForwardRange>
RotatedRangeType<ForwardRange>
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

}  // namespace extension_implementation::boost::adaptors

// #######################################################################
// ##  doubled                                                          ##
// #######################################################################

namespace extension_implementation::boost::adaptors {

class Doubler {};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static Doubler doubled{};
#pragma GCC diagnostic pop

template <typename ForwardRange>
using DoubledRangeType = ::boost::joined_range<const ForwardRange, const ForwardRange>;

template <typename ForwardRange>
DoubledRangeType<ForwardRange>
operator|(const ForwardRange &rng, const Doubler &) {
    return ::boost::join(rng, rng);
}

}  // namespace extension_implementation::boost::adaptors

// #######################################################################
// ##  refined                                                          ##
// #######################################################################

// namespace extension_implementation::boost::adaptors {

// template<typename T>
// class RefinedHolder {
//    public:
//     RefinedHolder(size_t n, T v) : n(n), v(v){};
//     size_t n;
//     T[1] v;
// };

// inline RefinedHolder refined(size_t n) {
//     return RefinedHolder(n);
// }

// /*
//  * Args domain: h.n has to be a positive number, less than the range size.
//  * The performed refinement is only on the h.n site,
//  * and relies on the original value replacement with the new value h.v[0];
//  */

// template <typename ForwardRange>
// using RefinedRangeType = ::boost::joined_range<
//     const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>,
//     const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>,
//     const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>>;

// template <typename ForwardRange>
// RefinedRangeType<ForwardRange>
// operator|(const ForwardRange &rng, const RefinedHolder &h) {
// #ifndef NDEBUG
//     const auto d = std::distance(std::begin(rng), std::end(rng));
//     assert(d >= 0);
//     assert(::boost::numeric_cast<decltype(d)>(h.n) <= d);
//     assert(h.n >= 0);
// #endif
//     const auto mid = std::next(std::begin(rng), h.n);
//     return ::boost::join(::boost::make_iterator_range(mid, std::end(rng)),
//                          ::boost::make_iterator_range(std::begin(rng), mid));
// }

// }  // namespace extension_implementation::boost::adaptors

// ########################################################
// ##  export API to extension::boost::adaptors          ##
// ########################################################

namespace extension::boost::adaptors {

using extension_implementation::boost::adaptors::doubled;
using extension_implementation::boost::adaptors::DoubledRangeType;
using extension_implementation::boost::adaptors::rotated;
using extension_implementation::boost::adaptors::RotatedRangeType;

}  // namespace extension::boost::adaptors

#endif