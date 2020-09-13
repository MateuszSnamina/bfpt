#ifndef KSTATE_BASIS_HPP
#define KSTATE_BASIS_HPP

#include <kstate/vec_map.hpp>

// #######################################################################
// ## Basis                                                             ##
// #######################################################################

namespace kstate {

template <typename _Element>
class Basis {
   public:
    // Helper types:
    using Element = _Element;
    // Constructors:
    Basis(size_t n_sites);
    // The physical realm:
    size_t n_sites() const;
    // Container type definition -- index typedefs:
    using VecIndex = typename VecMap<Element>::VecIndex;
    using ElementPtr = typename VecMap<Element>::ElementPtr;
    // Container size:
    unsigned size() const;
    // Container index interfaces:
    VecIndex& vec_index();
    const VecIndex& vec_index() const;
    // Container Member functions -- accesors:
    template <typename OtherRangeType>
    boost::optional<unsigned> find_element_and_get_its_ra_index(const OtherRangeType& v) const;
    // Container Member functions -- modifiers:
    void add_element(ElementPtr c);

   private:
    size_t _n_sites;
    VecMap<Element> _vec_map;
};

// *************************************************************************

template <typename _Element>
Basis<_Element>::Basis(size_t n_sites)
    : _n_sites(n_sites) {
}

// *************************************************************************

template <typename _Element>
size_t
Basis<_Element>::n_sites() const {
    return _n_sites;
}

template <typename _Element>
unsigned
Basis<_Element>::size() const {
    return _vec_map.size();
}

template <typename _Element>
typename Basis<_Element>::VecIndex&
Basis<_Element>::vec_index() {
    return _vec_map.vec_index();
}

template <typename _Element>
const typename Basis<_Element>::VecIndex&
Basis<_Element>::vec_index() const {
    return _vec_map.vec_index();
}

template <typename _Element>
template <typename OtherRangeType>
boost::optional<unsigned>
Basis<_Element>::find_element_and_get_its_ra_index(const OtherRangeType& v) const {
    return _vec_map.find_element_and_get_its_ra_index(v);
}

template <typename _Element>
void Basis<_Element>::add_element(ElementPtr c) {
    assert(c);
    assert(c->n_sites() == n_sites());
    _vec_map.vec_index().push_back(c);
}

}  // namespace kstate

#endif
