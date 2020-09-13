#ifndef KSTATE_BASIS_HPP
#define KSTATE_BASIS_HPP

#include <kstate/vec_map.hpp>

// #######################################################################
// ## Basis                                                             ##
// #######################################################################

namespace kstate {

template <typename _ElementT>
class Basis {
   public:
    // Helper types:
    using ElementT = _ElementT;
    // Constructors:
    Basis(size_t n_sites);
    // The physical realm:
    size_t n_sites() const;
    // Container type definition -- index typedefs:
    using VecIndexT = typename VecMap<ElementT>::VecIndex;
    using ElementPtrT = typename VecMap<ElementT>::ElementPtr;
    // Container size:
    unsigned size() const;
    // Container index interfaces:
    VecIndexT& vec_index();
    const VecIndexT& vec_index() const;
    // Container Member functions -- accesors:
    template <typename OtherRangeType>
    boost::optional<unsigned> find_element_and_get_its_ra_index(const OtherRangeType& v) const;
    // Container Member functions -- modifiers:
    void add_element(ElementPtrT c);

   private:
    size_t _n_sites;
    VecMap<ElementT> _vec_map;
};

// *************************************************************************

template <typename _ElementT>
Basis<_ElementT>::Basis(size_t n_sites)
    : _n_sites(n_sites) {
}

// *************************************************************************

template <typename _ElementT>
size_t
Basis<_ElementT>::n_sites() const {
    return _n_sites;
}

template <typename _ElementT>
unsigned
Basis<_ElementT>::size() const {
    return _vec_map.size();
}

template <typename _ElementT>
typename Basis<_ElementT>::VecIndexT&
Basis<_ElementT>::vec_index() {
    return _vec_map.vec_index();
}

template <typename _ElementT>
const typename Basis<_ElementT>::VecIndexT&
Basis<_ElementT>::vec_index() const {
    return _vec_map.vec_index();
}

template <typename _ElementT>
template <typename OtherRangeType>
boost::optional<unsigned>
Basis<_ElementT>::find_element_and_get_its_ra_index(const OtherRangeType& v) const {
    return _vec_map.find_element_and_get_its_ra_index(v);
}

template <typename _ElementT>
void Basis<_ElementT>::add_element(ElementPtrT c) {
    assert(c);
    assert(c->n_sites() == n_sites());
    _vec_map.vec_index().push_back(c);
}

}  // namespace kstate

#endif
