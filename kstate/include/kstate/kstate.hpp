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
    const bool result_cmp = boost::lexicographical_compare(
        rng | extension::boost::adaptors::rotated(i),
        rng | extension::boost::adaptors::rotated(_));
    if (result_cmp) {
      i = _;
    }
  }
  return i;
}

}  // namespace

namespace kstate {

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
// ## DynamicKstate                                                     ##
// #######################################################################

// Helper tag classes:
struct CtrFromRange {};
struct CtrFromBuffer {};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static CtrFromRange ctr_from_range{};
static CtrFromBuffer ctr_from_buffer{};
#pragma GCC diagnostic pop

// ***********************************************************************

// Helper functions:
template <typename Range>
std::vector<typename boost::range_value<Range>::type> init_vector_from_range(
    const Range& rng) {
  using Element = typename boost::range_value<Range>::type;
  using Vector = std::vector<Element>;
  return Vector(std::begin(rng), std::end(rng));
}

// ***********************************************************************

/*
 * DynamicKstate uses std::vect as am internall buffer.
 */

template <typename SiteType>
struct DynamicKstateTypes {
  using BufferType = typename std::vector<SiteType>;
  using IteratorType = typename BufferType::iterator;
  using ConstIteratorType = typename BufferType::const_iterator;
  using RangeType = typename boost::iterator_range<IteratorType>;
  using ConstRangeType = typename boost::iterator_range<ConstIteratorType>;
};

template <typename SiteType>
class DynamicKstate
    : public Kstate<typename DynamicKstateTypes<SiteType>::ConstRangeType> {
  using BufferType = typename DynamicKstateTypes<SiteType>::BufferType;
  using IteratorType = typename DynamicKstateTypes<SiteType>::IteratorType;
  using ConstIteratorType =
      typename DynamicKstateTypes<SiteType>::ConstIteratorType;
  using RangeType = typename DynamicKstateTypes<SiteType>::RangeType;
  using ConstRangeType = typename DynamicKstateTypes<SiteType>::ConstRangeType;

 public:
  DynamicKstate(BufferType&&, CtrFromBuffer);
  template <typename OtherRangeType>
  DynamicKstate(const OtherRangeType&, CtrFromRange);

 public:
  ConstRangeType to_range() const override;
  size_t n_sites() const override;

 protected:
  const BufferType _v;
};

// ***********************************************************************

template <typename SiteType>
DynamicKstate<SiteType>::DynamicKstate(DynamicKstate<SiteType>::BufferType&& v,
                                       CtrFromBuffer)
    : _v(std::move(v)) {}

template <typename SiteType>
template <typename OtherRangeType>
DynamicKstate<SiteType>::DynamicKstate(const OtherRangeType& r, CtrFromRange)
    : _v(init_vector_from_range(r)) {}
//    : _v(std::begin(r), std::end(r)) {}

// ***********************************************************************

template <typename SiteType>
typename DynamicKstate<SiteType>::ConstRangeType
DynamicKstate<SiteType>::to_range() const {
  return _v;
}

template <typename SiteType>
size_t DynamicKstate<SiteType>::n_sites() const {
  return _v.size();
}

// #######################################################################
// ## UniqueDynamicKstate                                               ##
// #######################################################################

template <typename SiteType>
class DynamicUniqueKstate
    : public Kstate<typename DynamicKstateTypes<SiteType>::ConstRangeType> {
  // public DynamicKstate<SiteType> {
  using BufferType = typename DynamicKstateTypes<SiteType>::BufferType;
  using IteratorType = typename DynamicKstateTypes<SiteType>::IteratorType;
  using ConstIteratorType =
      typename DynamicKstateTypes<SiteType>::ConstIteratorType;
  using RangeType = typename DynamicKstateTypes<SiteType>::RangeType;
  using ConstRangeType = typename DynamicKstateTypes<SiteType>::ConstRangeType;

 public:
  template <typename SomeRangeType>
  DynamicUniqueKstate(const SomeRangeType& v, CtrFromRange);

 public:
  ConstRangeType to_range() const override;
  size_t n_sites() const override;

 protected:
  const BufferType _v;
};

// ***********************************************************************

template <typename SiteType>
template <typename SomeRangeType>
DynamicUniqueKstate<SiteType>::DynamicUniqueKstate(const SomeRangeType& r,
                                                   CtrFromRange)
    : _v(init_vector_from_range(
          r | extension::boost::adaptors::rotated(n_unique_shift(r)))) {}

// ***********************************************************************

template <typename SiteType>
typename DynamicUniqueKstate<SiteType>::ConstRangeType
DynamicUniqueKstate<SiteType>::to_range() const {
  return _v;
}

template <typename SiteType>
size_t DynamicUniqueKstate<SiteType>::n_sites() const {
  return _v.size();
}

}  // namespace kstate

#endif