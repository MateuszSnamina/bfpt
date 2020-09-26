#ifndef KSTATE_VEC_MAP_HPP
#define KSTATE_VEC_MAP_HPP

#include <kstate/kstate_comparator.hpp>

#include <boost/multi_index/mem_fun.hpp>
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
 * The class Element has to implament a member function:
 * SomeRangeType Element::to_range() const.
 * The VecMap containers uses above member function
 * as a key extractor for the map.
 *
 * In the container objects are storred as shared_ptr to the objecs.
 */

template <typename _Element>
class VecMap {
   public:
    using Element = _Element;
    using Key = decltype(std::declval<Element>().to_range());
    using ElementPtr = std::shared_ptr<Element>;

   private:
    // Tags for random-access-index and search-index;
    struct Vec;
    struct Map;
    // Container type definition -- helper typedefs:
    using VecTagDef = boost::multi_index::tag<Vec>;
    using MapTagDef = boost::multi_index::tag<Map>;
    using KayExtractorDef = boost::multi_index::const_mem_fun<Element, Key, &Element::to_range>;
    // Container type definition -- index typedefs:
    using VecIndexDef = boost::multi_index::random_access<VecTagDef>;
    using MapIndexDef = boost::multi_index::ordered_unique<MapTagDef, KayExtractorDef, RangeComparator>;
    // Container type definition -- final container typedef:
    using Container = boost::multi_index::multi_index_container<ElementPtr, boost::multi_index::indexed_by<VecIndexDef, MapIndexDef>>;
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
    template <typename OtherRangeType>
    boost::optional<unsigned> find_element_and_get_its_ra_index(const OtherRangeType& v) const;
    // Container Member functions -- modifiers:
    void add_element(ElementPtr c);
};

template <typename _Element>
struct VecMap<_Element>::Vec{};

template <typename _Element>
struct VecMap<_Element>::Map{};

// *************************************************************************
// ********  Member functions definitions     ******************************
// *************************************************************************

template <typename _Element>
unsigned
VecMap<_Element>::size() const {
    return container.size();
}

template <typename _Element>
typename VecMap<_Element>::VecIndex&
VecMap<_Element>::vec_index() {
    return container.template get<Vec>();
}

template <typename _Element>
typename VecMap<_Element>::MapIndex&
VecMap<_Element>::map_index() {
    return container.template get<Map>();
}

template <typename _Element>
const typename VecMap<_Element>::VecIndex&
VecMap<_Element>::vec_index() const {
    return container.template get<Vec>();
}

template <typename _Element>
const typename VecMap<_Element>::MapIndex&
VecMap<_Element>::map_index() const {
    return container.template get<Map>();
}

template <typename _Element>
template <typename OtherRangeType>
boost::optional<unsigned>
VecMap<_Element>::find_element_and_get_its_ra_index(const OtherRangeType& v) const {
    auto search_iter = map_index().find(v);
    if (search_iter == map_index().end()) return boost::optional<unsigned>();
    auto ra_iter = container.template project<Vec>(search_iter);
    return std::distance(vec_index().begin(), ra_iter);
}

template <typename _Element>
void VecMap<_Element>::add_element(ElementPtr c) {
    vec_index().push_back(c);
}

}  // namespace kstate

#endif
