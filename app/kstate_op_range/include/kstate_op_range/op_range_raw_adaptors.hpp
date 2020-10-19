#pragma once

#include <kstate_view_amend_spec/amend_spec.hpp>

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

/*
 * Args domain: h.n has to be a positive number, less than the range size.
 * The performed rotation is a left rotation.
 */

namespace kstate_op_range::raw {

template <typename ForwardRange>
using RotatedRangeType = ::boost::joined_range<
    const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>,
    const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>>;

template <typename ForwardRange>
RotatedRangeType<ForwardRange>
rotated(const ForwardRange &rng, const kstate_view_amend_spec::RotateHolder &h) {
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

}  // namespace kstate_op_range::raw

namespace kstate_op_range::raw::adaptors {

template <typename ForwardRange>
RotatedRangeType<ForwardRange>
operator|(const ForwardRange &rng, const kstate_view_amend_spec::RotateHolder &h) {
    return rotated<ForwardRange>(rng, h);
}

}  // namespace kstate_op_range::raw::adaptors

// #######################################################################
// ##  doubled                                                          ##
// #######################################################################

namespace kstate_op_range::raw {

template <typename ForwardRange>
using DoubledRangeType = ::boost::joined_range<const ForwardRange, const ForwardRange>;

template <typename ForwardRange>
DoubledRangeType<ForwardRange>
doubled(const ForwardRange &rng) {
    return ::boost::join(rng, rng);
}

}  // namespace kstate_op_range::raw

namespace kstate_op_range::raw::adaptors {

template <typename ForwardRange>
DoubledRangeType<ForwardRange>
operator|(const ForwardRange &rng, const kstate_view_amend_spec::Doubler &) {
    return doubled<ForwardRange>(rng);
}

}  // namespace kstate_op_range::raw::adaptors

// #######################################################################
// ##  refined                                                          ##
// #######################################################################

/*
 * Args domain: h.n has to be a positive number, less than the range size.
 * The performed refinement is only on the h.n site,
 * and relies on the original value replacement with the new value h.v[0];
 */

namespace kstate_op_range::raw {

template <typename ForwardRange, typename T>
using RefinedRangeType =
    ::boost::joined_range<
        const ::boost::joined_range<
            const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>,
            const typename ::boost::iterator_range<typename ::boost::range_iterator<const T[1]>::type>>,
        const typename ::boost::iterator_range<typename ::boost::range_iterator<const ForwardRange>::type>>;

template <typename ForwardRange, typename T>
RefinedRangeType<ForwardRange, T>
refined(const ForwardRange &rng, const kstate_view_amend_spec::RefinedHolder<T> &h) {
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

}  // namespace kstate_op_range::raw

namespace kstate_op_range::raw::adaptors {

template <typename ForwardRange, typename T>
RefinedRangeType<ForwardRange, T>
operator|(const ForwardRange &rng, const kstate_view_amend_spec::RefinedHolder<T> &h) {
    return refined<ForwardRange, T>(rng, h);
}

}  // namespace kstate_op_range::raw::adaptors
