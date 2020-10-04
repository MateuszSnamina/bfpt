#ifndef KSTATE_BASIS_HPP
#define KSTATE_BASIS_HPP

#include <kstate/vec_map.hpp>
#include <kstate/trait_kstate.hpp>

#include <kstate/is_base_of_template.hpp> // needed only for static assert.//TODO remove
#include <kstate/kstate_abstract.hpp> // needed only for static assert.//TODO remove

// #######################################################################
// ## Basis                                                             ##
// #######################################################################

namespace kstate {

template <typename _KstateTrait>
class Basis {
//    static_assert(!std::is_array_v<_KstateT>);
//    static_assert(!std::is_function_v<_KstateT>);
//    static_assert(!std::is_void_v<std::decay<_KstateT>>);
//    static_assert(!std::is_null_pointer_v<std::decay<_KstateT>>);
//    static_assert(!std::is_pointer_v<std::decay<_KstateT>>);
//    static_assert(!std::is_member_object_pointer_v<_KstateT>);
//    static_assert(!std::is_member_function_pointer_v<_KstateT>);
//    static_assert(!std::is_const_v<_KstateT>);
//    static_assert(!std::is_volatile_v<_KstateT>);
//    static_assert(!std::is_reference_v<_KstateT>);
//    static_assert(std::is_enum_v<std::decay<_KstateT>> || std::is_union_v<std::decay<_KstateT>> || std::is_class_v<std::decay<_KstateT>>);//TODO remove
    static_assert(_KstateTrait::is_kstate_trait);
public:
    // Helper types:
    using KstateTrait = _KstateTrait;
    using KstateT = typename _KstateTrait::KstateT;
    using KStatePtrT = typename VecMap<KstateT>::ElementPtr;
    using VecIndexT = typename VecMap<KstateT>::VecIndex;
public:
    // Constructors:
    Basis(size_t n_sites);
    // The physical realm:
    size_t n_sites() const;
    // Container size:
    unsigned size() const;
    // Container index interfaces:
    VecIndexT& vec_index();
    const VecIndexT& vec_index() const;
    // Container Member functions -- accesors:
    template <typename OtherRangeType>
    boost::optional<unsigned> find_element_and_get_its_ra_index(const OtherRangeType& v) const;
    // Container Member functions -- modifiers:
    void add_element(KStatePtrT c);
private:
    size_t _n_sites;
    VecMap<KstateT> _vec_map;
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
const typename Basis<_KstateTrait>::VecIndexT&
Basis<_KstateTrait>::vec_index() const {
    return _vec_map.vec_index();
}

template <typename _KstateTrait>
template <typename OtherRangeType>
boost::optional<unsigned>
Basis<_KstateTrait>::find_element_and_get_its_ra_index(const OtherRangeType& v) const {
    return _vec_map.find_element_and_get_its_ra_index(v);
}

template <typename _KstateTrait>
void Basis<_KstateTrait>::add_element(KStatePtrT c) {
    assert(c);
    assert(c->n_sites() == n_sites());
    _vec_map.vec_index().push_back(c);
}

}  // namespace kstate

#endif
