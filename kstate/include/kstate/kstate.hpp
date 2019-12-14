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
// ## TagClasses                                                        ##
// #######################################################################

struct FromRange {};
struct FromBuffer {};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static FromRange from_range{};
static FromBuffer from_buffer{};
#pragma GCC diagnostic pop

// #######################################################################
// ## Kstate                                                            ##
// #######################################################################

template <typename ConstRangeType>
class Kstate {
 public:  // helper types:
  using SiteType = typename ::boost::range_value<ConstRangeType>::type;

 public:  // base on the two functions:
  virtual ConstRangeType to_range() const = 0;
  virtual size_t n_sites() const = 0;

 public:  // Kstate implements:
  size_t n_least_replication_shift() const;
  bool is_prolific(int n_k) const;
  template <typename OtherConstRangeType>
  bool compare(const Kstate<OtherConstRangeType>& other) const;
  template <typename OtherConstRangeType>
  std::optional<size_t> tranlational_compare(
      const Kstate<OtherConstRangeType>& other) const;
  std::string to_str() const;

 public:
  virtual ~Kstate() = default;
};

// ***********************************************************************

template <typename ConstRangeType>
size_t Kstate<ConstRangeType>::n_least_replication_shift() const {
  assert(n_sites() > 0);
  const auto r = to_range();
  const auto rdr = r | extension::boost::adaptors::doubled |
                   extension::boost::adaptors::rotated(1);
  const auto it = boost::range::search(rdr, r);
  const auto _ = std::distance(std::begin(rdr), it);
  assert(_ >= 0);
  return static_cast<size_t>(_ + 1);
}

template <typename ConstRangeType>
bool Kstate<ConstRangeType>::is_prolific(int n_k) const {
  return !((n_least_replication_shift() * n_k) % n_sites());
}

template <typename ConstRangeType>
template <typename OtherConstRangeType>
bool Kstate<ConstRangeType>::compare(
    const Kstate<OtherConstRangeType>& other) const {
  return boost::range::equal(to_range(), other.to_range());
}

template <typename ConstRangeType>
template <typename OtherConstRangeType>
std::optional<size_t> Kstate<ConstRangeType>::tranlational_compare(
    const Kstate<OtherConstRangeType>& other) const {
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

template <typename ConstRangeType>
std::string Kstate<ConstRangeType>::to_str() const {
  return extension::boost::RangeStreamStreamer()
      .set_stream_preparer([](std::ostream& s) { s << "⦃"; })
      .set_stream_sustainer([](std::ostream& s, size_t i) {})
      .set_stream_separer([](std::ostream& s) { s << "∙"; })
      .set_stream_finisher([](std::ostream& s) { s << "⦄"; })
      .stream(to_range())
      .str();
}

// #######################################################################
// ## SimpleKstate                                                      ##
// #######################################################################

/*
 * SimpleKstate uses std::vect as a buffer.
 * This is to make the preliminary test of concept.
 */

template <typename SiteType>
struct SimpleKstateTypes {
  using BufferType = typename std::vector<SiteType>;
  using IteratorType = typename BufferType::iterator;
  using ConstIteratorType = typename BufferType::const_iterator;
  using RangeType = typename boost::iterator_range<IteratorType>;
  using ConstRangeType = typename boost::iterator_range<ConstIteratorType>;
};

template <typename SiteType>
class SimpleKstate
    : public Kstate<typename SimpleKstateTypes<SiteType>::ConstRangeType> {
  using BufferType = typename SimpleKstateTypes<SiteType>::BufferType;
  using IteratorType = typename SimpleKstateTypes<SiteType>::IteratorType;
  using ConstIteratorType =
      typename SimpleKstateTypes<SiteType>::ConstIteratorType;
  using RangeType = typename SimpleKstateTypes<SiteType>::RangeType;
  using ConstRangeType = typename SimpleKstateTypes<SiteType>::ConstRangeType;

 public:
  SimpleKstate(BufferType&&);
  // SimpleKstate(const SimpleKstate<SiteType>&);  // needed????? TODO check
  template <typename OtherRangeType>
  SimpleKstate(const OtherRangeType&, FromRange);

 public:
  ConstRangeType to_range() const override;
  size_t n_sites() const override;

 private:
  const BufferType _v;
  const size_t _n_sites;
};

// ***********************************************************************

template <typename SiteType>
SimpleKstate<SiteType>::SimpleKstate(SimpleKstate<SiteType>::BufferType&& v)
    : _v(std::move(v)), _n_sites(_v.size()) {}

template <typename SiteType>
template <typename OtherRangeType>
SimpleKstate<SiteType>::SimpleKstate(const OtherRangeType& r, FromRange)
    : _v(std::begin(r), std::end(r)), _n_sites(_v.size()) {}

// template <typename SiteType>
// SimpleKstate<SiteType>::SimpleKstate(const SimpleKstate<SiteType>& s)
//     : _v(s._v), _n_sites(s.n_sites()) {}

// ***********************************************************************

template <typename SiteType>
typename SimpleKstate<SiteType>::ConstRangeType
SimpleKstate<SiteType>::to_range() const {
  return _v;
}

template <typename SiteType>
size_t SimpleKstate<SiteType>::n_sites() const {
  return _n_sites;
}

// #######################################################################
// ## KstateUniqueView                                                  ##
// #######################################################################

template <typename ViewedRangeType>
struct KstateUniqueViewTypes {
  using RangeType =
      typename extension::boost::adaptors::RotatedRangeType<ViewedRangeType>;
};

template <typename ViewedRangeType>
class KstateUniqueView
    : public Kstate<
          typename KstateUniqueViewTypes<ViewedRangeType>::RangeType> {
  using RangeType = typename KstateUniqueViewTypes<ViewedRangeType>::RangeType;

 public:
  KstateUniqueView(const Kstate<ViewedRangeType>&);

 public:
  RangeType to_range() const override;
  size_t n_sites() const override;

 private:
  const Kstate<ViewedRangeType>& _r;
  const size_t _n_unique_shift;
};

// ***********************************************************************

template <typename ViewedRangeType>
KstateUniqueView<ViewedRangeType>::KstateUniqueView(
    const Kstate<ViewedRangeType>& r)
    : _r(r), _n_unique_shift(n_unique_shift(r)) {}

// ***********************************************************************

template <typename ViewedRangeType>
typename KstateUniqueView<ViewedRangeType>::RangeType
KstateUniqueView<ViewedRangeType>::to_range() const {
  return _r | _n_unique_shift;
}

template <typename ViewedRangeType>
size_t KstateUniqueView<ViewedRangeType>::n_sites() const {
  return _r.n_sites();
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
  SimpleUniqueKstate(const SomeRangeType& v, FromRange);
};

// ***********************************************************************

template <typename SiteType>
template <typename SomeRangeType>
SimpleUniqueKstate<SiteType>::SimpleUniqueKstate(const SomeRangeType& r,
                                                 FromRange)
    : SimpleKstate<SiteType>(
          r | extension::boost::adaptors::rotated(n_unique_shift(r)),
          from_range) {}

}  // namespace kstate

#endif