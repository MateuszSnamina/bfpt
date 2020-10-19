#ifndef EXTENSIONS_BOOST_ADAPTORS_HPP
#define EXTENSIONS_BOOST_ADAPTORS_HPP

#include <boost/range.hpp>
#include <boost/range/join.hpp>

#include <iterator>

#ifndef NDEBUG
#include <boost/numeric/conversion/cast.hpp>
#include <cassert>
#endif

// #######################################################################
// ##  rotated                                                          ##
// #######################################################################

namespace kstate_view_amend_spec {

class RotateHolder {
   public:
    RotateHolder(size_t n) : n(n){};
    const size_t n;
};

inline RotateHolder rotated(size_t n) {
    return RotateHolder(n);
}

} // namespace kstate_view_amend_spec

namespace kstate_op_range::raw::adaptors::impl {

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
operator|(const ForwardRange &rng, const kstate_view_amend_spec::RotateHolder &h) {
#ifndef NDEBUG
    const auto d = ::boost::size(rng);
    assert(d >= 0);
    assert(::boost::numeric_cast<decltype(d)>(h.n) <= d);
    assert(h.n >= 0);
#endif
    const auto mid = ::std::next(::std::begin(rng), h.n);
    return ::boost::join(::boost::make_iterator_range(mid, ::std::end(rng)),
                         ::boost::make_iterator_range(::std::begin(rng), mid));
}

}  // namespace kstate_op_range::raw::adaptors::impl

// #######################################################################
// ##  doubled                                                          ##
// #######################################################################

namespace kstate_view_amend_spec {

class Doubler {};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static Doubler doubled{};
#pragma GCC diagnostic pop

} // namespace kstate_view_amend_spec

namespace kstate_op_range::raw::adaptors::impl {

template <typename ForwardRange>
using DoubledRangeType = ::boost::joined_range<const ForwardRange, const ForwardRange>;

template <typename ForwardRange>
DoubledRangeType<ForwardRange>
operator|(const ForwardRange &rng, const kstate_view_amend_spec::Doubler &) {
    return ::boost::join(rng, rng);
}

}  // namespace kstate_op_range::raw::adaptors::impl

// #######################################################################
// ##  refined                                                          ##
// #######################################################################

namespace kstate_view_amend_spec {

template <typename T>
struct RefinedHolder {
    RefinedHolder(size_t n, T v)
        : n(n),
          v{v} {};
    const size_t n;
    const T v[1];
};

template <typename T>
RefinedHolder<T> refined(size_t n, T v) {
    return RefinedHolder<T>(n, v);
}

} // namespace kstate_view_amend_spec

namespace kstate_op_range::raw::adaptors::impl {

// /*
//  * Args domain: h.n has to be a positive number, less than the range size.
//  * The performed refinement is only on the h.n site,
//  * and relies on the original value replacement with the new value h.v[0];
//  */

template <typename ForwardRange, typename T>
using RefinedRangeType =
    ::boost::joined_range<
        const ::boost::joined_range<
            const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>,
            const typename ::boost::iterator_range<typename ::boost::range_iterator<const T[1]>::type>>,
        const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>>;

template <typename ForwardRange, typename T>
RefinedRangeType<ForwardRange, T>
operator|(const ForwardRange &rng, const kstate_view_amend_spec::RefinedHolder<T> &h) {
#ifndef NDEBUG
    const auto d = ::boost::size(rng);
    assert(d >= 0);
    assert(::boost::numeric_cast<decltype(d)>(h.n) <= d);
    assert(h.n >= 0);
#endif
    const auto &mid = ::std::next(::std::begin(rng), h.n);
    const auto &mid_1 = ::std::next(::std::begin(rng), h.n + 1);
    const auto &r1 = ::boost::make_iterator_range(::std::begin(rng), mid);
    const auto &r2 = ::boost::make_iterator_range(::std::begin(h.v), ::std::end(h.v));
    const auto &r3 = ::boost::make_iterator_range(mid_1, ::std::end(rng));
    return ::boost::join(::boost::join(r1, r2), r3);
}

}  // namespace kstate_op_range::raw::adaptors::impl

// ########################################################
// ##  export API to extension::boost::adaptors          ##
// ########################################################

namespace kstate_op_range::raw::adaptors {

using kstate_op_range::raw::adaptors::impl::DoubledRangeType;
using kstate_op_range::raw::adaptors::impl::RotatedRangeType;
using kstate_op_range::raw::adaptors::impl::RefinedRangeType;
using kstate_op_range::raw::adaptors::impl::operator|;

}  // namespace kstate_op_range::raw::adaptors

#endif
