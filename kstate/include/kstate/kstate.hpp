#ifndef KSTATE_KSTATE_HPP
#define KSTATE_KSTATE_HPP

#include <boost/algorithm/string/predicate.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm/search.hpp>
#include <extensions/adaptors.hpp>
#include <extensions/range_streamer.hpp>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <iostream>  //debug
#include <iterator>

namespace {

template <typename ForwardRange>
size_t n_unique_shift(const ForwardRange& rng) {
  const auto d = std::distance(std::begin(rng), std::end(rng));
  size_t i = 0;
  for (size_t _ = 1; _ < d; _++) {
    // std::cout << extension::boost::RangeStreamStreamer()
    //                  .stream(rng | extension::boost::adaptors::rotated(i))
    //                  .str()
    //           << " vs "
    //           << extension::boost::RangeStreamStreamer()
    //                  .stream(rng | extension::boost::adaptors::rotated(_))
    //                  .str()
    //           << std::endl;
    // std::cout << std::boolalpha;
    const bool result_cmp = boost::lexicographical_compare(
        rng | extension::boost::adaptors::rotated(i),
        rng | extension::boost::adaptors::rotated(_));
    // std::cout << result_cmp << std::endl;
    if (result_cmp) {
      i = _;
    }
  }
  //   std::cout << i << std::endl;
  return i;
}

}  // namespace

namespace kstate {

// #######################################################################
// ## SimpleKstate                                                      ##
// #######################################################################

template <typename SiteType>
class SimpleKstate {
  using BufferType = typename std::vector<SiteType>;
  using IteratorType = typename std::vector<SiteType>::iterator;
  using ConstIteratorType = typename std::vector<SiteType>::const_iterator;
  using RangeType = typename boost::iterator_range<IteratorType>;
  using ConstRangeType = typename boost::iterator_range<ConstIteratorType>;

 public:
  SimpleKstate(BufferType&& v);
  template <typename SomeRangeType>
  SimpleKstate(const SomeRangeType& v);
  ConstRangeType to_range() const;
  size_t n_sites() const;
  size_t n_least_replication_shift() const;
  bool is_prolific(int n_k) const;
  bool compare(const SimpleKstate<SiteType>& other) const;
  std::optional<size_t> tranlational_compare(
      const SimpleKstate<SiteType>& other) const;
  std::string to_str() const;

 private:
  const BufferType _v;
  const size_t _n_sites;
};

// #######################################################################

template <typename SiteType>
SimpleKstate<SiteType>::SimpleKstate(SimpleKstate<SiteType>::BufferType&& v)
    : _v(std::move(v)), _n_sites(_v.size()) {}

template <typename SiteType>
template <typename SomeRangeType>
SimpleKstate<SiteType>::SimpleKstate(const SomeRangeType& r)
    : _v(std::begin(r), std::end(r)), _n_sites(_v.size()) {}

template <typename SiteType>
typename SimpleKstate<SiteType>::ConstRangeType
SimpleKstate<SiteType>::to_range() const {
  return _v;
}

template <typename SiteType>
size_t SimpleKstate<SiteType>::n_sites() const {
  return _n_sites;
}

template <typename SiteType>
size_t SimpleKstate<SiteType>::n_least_replication_shift() const {
  assert(n_sites() > 0);
  const auto r = to_range();
  const auto rdr = r | extension::boost::adaptors::doubled |
                   extension::boost::adaptors::rotated(1);
  const auto it = boost::range::search(rdr, r);
  const auto _ = std::distance(std::begin(rdr), it);
  assert(_ >= 0);
  return static_cast<size_t>(_ + 1);
}

template <typename SiteType>
bool SimpleKstate<SiteType>::is_prolific(int n_k) const {
  return !((n_least_replication_shift() * n_k) % n_sites());
}

template <typename SiteType>
bool SimpleKstate<SiteType>::compare(
    const SimpleKstate<SiteType>& other) const {
  return boost::range::equal(to_range(), other.to_range());
}

template <typename SiteType>
std::optional<size_t> SimpleKstate<SiteType>::tranlational_compare(
    const SimpleKstate<SiteType>& other) const {
  if (n_sites() != other.n_sites()) {
    return std::nullopt;
  }
  const auto r1 = to_range();
  const auto r2 = other.to_range();
  const auto r2d = r2 | extension::boost::adaptors::doubled;
  const auto it = boost::range::search(r2d, r1);
  return it == std::end(r2d)
             ? std::optional<size_t>()
             : static_cast<size_t>(std::distance(std::begin(r2d), it));
}

template <typename SiteType>
std::string SimpleKstate<SiteType>::to_str() const {
  return extension::boost::RangeStreamStreamer()
      .set_stream_preparer([](std::ostream& s) { s << "⦃"; })
      .set_stream_sustainer([](std::ostream& s, size_t i) {})
      .set_stream_separer([](std::ostream& s) { s << "∙"; })
      .set_stream_finisher([](std::ostream& s) { s << "⦄"; })
      .stream(to_range())
      .str();
}

// #######################################################################
// ## UniqueSimpleKstate                                                ##
// #######################################################################

template <typename SiteType>
class SimpleUniqueKstate : public SimpleKstate<SiteType> {
  using BufferType = typename std::vector<SiteType>;
  using IteratorType = typename std::vector<SiteType>::iterator;
  using ConstIteratorType = typename std::vector<SiteType>::const_iterator;
  using RangeType = typename boost::iterator_range<IteratorType>;
  using ConstRangeType = typename boost::iterator_range<ConstIteratorType>;

 public:
  template <typename SomeRangeType>
  SimpleUniqueKstate(const SomeRangeType& v);
};

// #######################################################################

template <typename SiteType>
template <typename SomeRangeType>
SimpleUniqueKstate<SiteType>::SimpleUniqueKstate(const SomeRangeType& r)
    : SimpleKstate<SiteType>(
          r | extension::boost::adaptors::rotated(n_unique_shift(r))) {}

}  // namespace kstate

#endif