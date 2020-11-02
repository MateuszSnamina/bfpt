#pragma once

#include <boost/range.hpp>

#include <iterator>
#include <type_traits>
#include <vector>
#include <array>
#include <cassert>

// #######################################################################
// ## init_vector_from_range (a helper function)                        ##
// #######################################################################

namespace kstate_impl::helpers {

template <typename Range>
std::vector<typename boost::range_value<Range>::type>
init_vector_from_range(
    const Range& rng) {
    using Element = typename boost::range_value<Range>::type;
    using Vector = std::vector<Element>;
    return Vector(std::begin(rng), std::end(rng));
}

}  // end of namespace kstate_impl::helpers

// #######################################################################
// ## init_array_from_range (a helper function)                         ##
// #######################################################################

// Solution coppied from:
// https://stackoverflow.com/questions/10929202/initialize-stdarray-with-a-range-pair-of-iterators

namespace kstate_impl::helpers {

template <std::size_t... Indices>
struct indices {
    using next = indices<Indices..., sizeof...(Indices)>;
};

template <std::size_t N>
struct build_indices {
    using type = typename build_indices<N - 1>::type::next;
};

template <>
struct build_indices<0> {
    using type = indices<>;
};

template <std::size_t N>
using BuildIndices = typename build_indices<N>::type;

template <typename Iterator>
using ValueType = typename std::iterator_traits<Iterator>::value_type;

// internal overload with indices tag
template <std::size_t... I, typename RandomAccessIterator,
          typename Array = std::array<ValueType<RandomAccessIterator>, sizeof...(I)>>
Array make_array_impl(RandomAccessIterator first, indices<I...>) {
    return Array{{first[I]...}};
}

// externally visible interface
template <std::size_t N, typename RandomAccessIterator>
std::array<ValueType<RandomAccessIterator>, N>
make_array(RandomAccessIterator first, RandomAccessIterator last) {
    // last is not relevant if we're assuming the size is N
    // I'll assert it is correct anyway
    assert(last - first == N);
    return make_array_impl(first, BuildIndices<N>{});
}

template <std::size_t N, typename Range>
std::array<typename boost::range_value<Range>::type, N>
init_array_from_range(
    const Range& rng) {
    return make_array<N>(std::begin(rng), std::end(rng));
}

}  // end of namespace kstate_impl::helpers
