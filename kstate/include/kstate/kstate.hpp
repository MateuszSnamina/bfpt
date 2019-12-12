#ifndef KSTATE_KSTATE_HPP
#define KSTATE_KSTATE_HPP

#include <extensions/adaptors.hpp>
#include <extensions/range_streamer.hpp>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
// #include<algorithm>
#include <boost/range.hpp>
#include <boost/range/algorithm/search.hpp>

#include <boost/range/algorithm/copy.hpp>  //debug
#include <iostream>                        //debug
#include <iterator>                        //debug

// namespace {

// template <typename SinglePassRange>
// std::string range_to_str(const SinglePassRange& r) {
//   using ValueType = typename boost::range_value<const SinglePassRange>::type;
//   std::ostringstream oss;
//   std::ostream_iterator<ValueType> osi(oss, ",");
//   boost::copy(r, osi);
//   return oss.str();
// }

// template <typename SinglePassRange>
// void range_to_str(const SinglePassRange& r, std::ostream& os) {
//   using ValueType = typename boost::range_value<const SinglePassRange>::type;
//   std::ostream_iterator<ValueType> osi(os, ",");
//   boost::copy(r, osi);
// }

// }  // namespace

namespace kstate {

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
  // std::string SimpleKstate<SiteType>::to_repr() const;
  std::string to_str() const;

 private:
  const BufferType _v;
  const size_t _n_sites;
};

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
  return ! ((n_least_replication_shift() * n_k) % n_sites());
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

}  // namespace kstate

#endif