#ifndef KSTATE_BASIS_HPP
#define KSTATE_BASIS_HPP

#include <kstate_trait/trait_kstate.hpp>

#include <kstate_trait/trait_kstate.hpp>
#include <kstate_trait/kstate_comparator.hpp>

#include <kstate/vec_map.hpp>

//#include <boost/multi_index/mem_fun.hpp>

// #######################################################################
// ## Basis                                                             ##
// #######################################################################

namespace kstate {

template <typename _KstateTraitT>
struct BasisKeyExtractor {
    // helper types:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename KstateTraitT::KstateT;
    using ConstRangeT = typename KstateTraitT::ConstRangeT;
    // api needed by boost::multiindex:
    using result_type = ConstRangeT;
    result_type operator()(const std::shared_ptr<KstateT>& kstate_ptr) const{
        return KstateTraitT::to_range(*kstate_ptr);
    }
};


template <typename _KstateTraitT>
class Basis {
    static_assert(kstate_trait::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
public:
    // Helper types:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;
    using KstatePtrT = std::shared_ptr<KstateT>;
private:
    //using Key = decltype(std::declval<KstateT>().to_range());
    //using KeyExtractorT = boost::multi_index::const_mem_fun<KstateT, Key, &KstateT::to_range>;
    using KeyExtractorT = BasisKeyExtractor<KstateTraitT>;
    using ComparisonPredicateT = kstate_trait::RangeComparator;
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

template <typename _KstateTraitT>
Basis<_KstateTraitT>::Basis(size_t n_sites)
    : _n_sites(n_sites) {
}

// *************************************************************************

template <typename _KstateTraitT>
size_t
Basis<_KstateTraitT>::n_sites() const {
    return _n_sites;
}

template <typename _KstateTraitT>
unsigned
Basis<_KstateTraitT>::size() const {
    return _vec_map.size();
}

template <typename _KstateTraitT>
typename Basis<_KstateTraitT>::VecIndexT&
Basis<_KstateTraitT>::vec_index() {
    return _vec_map.vec_index();
}

template <typename _KstateTraitT>
typename Basis<_KstateTraitT>::MapIndexT&
Basis<_KstateTraitT>::map_index() {
    return _vec_map.map_index();
}

template <typename _KstateTraitT>
const typename Basis<_KstateTraitT>::VecIndexT&
Basis<_KstateTraitT>::vec_index() const {
    return _vec_map.vec_index();
}

template <typename _KstateTraitT>
const typename Basis<_KstateTraitT>::MapIndexT&
Basis<_KstateTraitT>::map_index() const {
    return _vec_map.map_index();
}

template <typename _KstateTraitT>
template <typename OtherRangeType>
boost::optional<unsigned>
Basis<_KstateTraitT>::find_element_and_get_its_ra_index(const OtherRangeType& v) const {
    return _vec_map.find_element_and_get_its_ra_index(v);
}

template <typename _KstateTraitT>
void Basis<_KstateTraitT>::add_element(KstatePtrT c) {
    assert(c);
    assert(c->n_sites() == n_sites());
    _vec_map.vec_index().push_back(c);
}

}  // namespace kstate

#endif
