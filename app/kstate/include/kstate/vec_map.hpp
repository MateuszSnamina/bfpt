#ifndef KSTATE_VEC_MAP_HPP
#define KSTATE_VEC_MAP_HPP

#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/optional.hpp>

#include <memory>

// #######################################################################
// ## VecMap                                                            ##
// #######################################################################

namespace kstate {
/*
 * VecMap<Element> is a map-like container and is a vec-like container
 * (simultaneously) for objects of class Element.
 *
 * In the container objects are storred as shared_ptr to the objecs.
 */

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
class VecMap {
   public:
    using ElementT = _ElementT;
    using ElementPtrT = std::shared_ptr<ElementT>;
    using KeyExtractorT = _KeyExtractorT;
    using ComparisonPredicateT = _ComparisonPredicateT;


   private:
    // Tags for random-access-index and search-index;
    struct Vec;
    struct Map;
    // Container type definition -- helper typedefs:
    using VecTagDef = boost::multi_index::tag<Vec>;
    using MapTagDef = boost::multi_index::tag<Map>;
    // Container type definition -- index typedefs:
    using VecIndexDef = boost::multi_index::random_access<VecTagDef>;
    using MapIndexDef = boost::multi_index::ordered_unique<MapTagDef, KeyExtractorT, ComparisonPredicateT>;
    // Container type definition -- final container typedef:
    using Container = boost::multi_index::multi_index_container<ElementPtrT, boost::multi_index::indexed_by<VecIndexDef, MapIndexDef>>;
    // The container:
    Container container;

   public:
    // Container type definition -- index typedefs:
    using VecIndex = typename Container::template index<Vec>::type;
    using MapIndex = typename Container::template index<Map>::type;
    // Container size:
    unsigned size() const;
    // Container index interfaces:
    VecIndex& vec_index();
    MapIndex& map_index();
    const VecIndex& vec_index() const;
    const MapIndex& map_index() const;
    // Container Member functions -- accesors:
    template <typename KeyComparableT>
    boost::optional<unsigned> find_element_and_get_its_ra_index(const KeyComparableT& v) const;
    // Container Member functions -- modifiers:
    void add_element(ElementPtrT c);
};

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
struct VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::Vec{};

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
struct VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::Map{};

// *************************************************************************
// ********  Member functions definitions     ******************************
// *************************************************************************

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
unsigned
VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::size() const {
    return container.size();
}

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
typename VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::VecIndex&
VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::vec_index() {
    return container.template get<Vec>();
}

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
typename VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::MapIndex&
VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::map_index() {
    return container.template get<Map>();
}

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
const typename VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::VecIndex&
VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::vec_index() const {
    return container.template get<Vec>();
}

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
const typename VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::MapIndex&
VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::map_index() const {
    return container.template get<Map>();
}

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
template <typename KeyComparableT>
boost::optional<unsigned>
VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::find_element_and_get_its_ra_index(const KeyComparableT& v) const {
    auto search_iter = map_index().find(v);
    if (search_iter == map_index().end()) return boost::optional<unsigned>();
    auto ra_iter = container.template project<Vec>(search_iter);
    return std::distance(vec_index().begin(), ra_iter);
}

template <typename _ElementT, typename _KeyExtractorT, typename _ComparisonPredicateT>
void VecMap<_ElementT, _KeyExtractorT, _ComparisonPredicateT>::add_element(ElementPtrT c) {
    vec_index().push_back(c);
}

}  // namespace kstate

#endif
