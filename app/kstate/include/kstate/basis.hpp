#ifndef KSTATE_BASIS_HPP
#define KSTATE_BASIS_HPP

#include <kstate/trait_kstate.hpp>
#include <kstate/vec_map.hpp>
#include <kstate/trait_kstate.hpp>
#include <kstate/kstate_comparator.hpp>

#include <boost/multi_index/mem_fun.hpp>

// #######################################################################
// ## Basis                                                             ##
// #######################################################################

namespace kstate {


template <typename _KstateTrait>
class Basis {
    static_assert(IsTraitKstate<_KstateTrait>::value);
    static_assert(_KstateTrait::is_kstate_trait);
public:
    // Helper types:
    using KstateTrait = _KstateTrait;
    using KstateT = typename _KstateTrait::KstateT;
    using KstatePtrT = std::shared_ptr<KstateT>;
private:
    using Key = decltype(std::declval<KstateT>().to_range());
    using KeyExtractorT = boost::multi_index::const_mem_fun<KstateT, Key, &KstateT::to_range>;
    using ComparisonPredicateT = RangeComparator;
public:
    using VecIndexT = typename VecMap<KstateT, KeyExtractorT, ComparisonPredicateT>::VecIndex;
    using MapIndexT = typename VecMap<KstateT, KeyExtractorT, ComparisonPredicateT>::MapIndex;

public:
    // Constructors:
    Basis(size_t n_sites);
    // The physical realm:
    size_t n_sites() const;
    // Container size:
    unsigned size() const;
    // Container index interfaces:
    VecIndexT& vec_index();
    MapIndexT& map_index();
    const VecIndexT& vec_index() const;
    const MapIndexT& map_index() const;
    // Container Member functions -- accesors:
    template <typename OtherRangeType>
    boost::optional<unsigned> find_element_and_get_its_ra_index(const OtherRangeType& v) const;
    // Container Member functions -- modifiers:
    void add_element(KstatePtrT c);
private:
    size_t _n_sites;
    VecMap<KstateT, KeyExtractorT, ComparisonPredicateT> _vec_map;
};

// *************************************************************************

template <typename _KstateTrait>
Basis<_KstateTrait>::Basis(size_t n_sites)
    : _n_sites(n_sites) {
}

// *************************************************************************

template <typename _KstateTrait>
size_t
Basis<_KstateTrait>::n_sites() const {
    return _n_sites;
}

template <typename _KstateTrait>
unsigned
Basis<_KstateTrait>::size() const {
    return _vec_map.size();
}

template <typename _KstateTrait>
typename Basis<_KstateTrait>::VecIndexT&
Basis<_KstateTrait>::vec_index() {
    return _vec_map.vec_index();
}

template <typename _KstateTrait>
typename Basis<_KstateTrait>::MapIndexT&
Basis<_KstateTrait>::map_index() {
    return _vec_map.map_index();
}

template <typename _KstateTrait>
const typename Basis<_KstateTrait>::VecIndexT&
Basis<_KstateTrait>::vec_index() const {
    return _vec_map.vec_index();
}

template <typename _KstateTrait>
const typename Basis<_KstateTrait>::MapIndexT&
Basis<_KstateTrait>::map_index() const {
    return _vec_map.map_index();
}

template <typename _KstateTrait>
template <typename OtherRangeType>
boost::optional<unsigned>
Basis<_KstateTrait>::find_element_and_get_its_ra_index(const OtherRangeType& v) const {
    return _vec_map.find_element_and_get_its_ra_index(v);
}

template <typename _KstateTrait>
void Basis<_KstateTrait>::add_element(KstatePtrT c) {
    assert(c);
    assert(c->n_sites() == n_sites());
    _vec_map.vec_index().push_back(c);
}

}  // namespace kstate

#endif
