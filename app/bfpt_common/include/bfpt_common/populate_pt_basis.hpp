#ifndef BFPT_COMMON_POPULATE_PT_BASIS_HPP
#define BFPT_COMMON_POPULATE_PT_BASIS_HPP

#include <bfpt_common/i_kstate_populator.hpp>

#include <kstate/basis.hpp>
#include <kstate/kstate_abstract.hpp>

#include <cassert>

// #######################################################################
// ##  populate_pt_basis                                                ##
// #######################################################################

/*
 * The function generates pt-basis of given order.
 * The input for the function is:
 *  - the desider pt-order,
 *  - the 0'th pt-order basis,
 *  - the populator object.
 * 
 * Desired pt-order is given as `max_pt_order` argument;
 * the 0'th pt-order basis is read from `basis` argument
 * (the argument has intent "inout",
 *  as input it gives the 0'th pt-order basis
 *  as output it gives the function result);
 * Populator object is given as `populator` arggument.
 * 
 */

namespace bfpt_common {

//template <typename SiteType>
template<typename KstateT>
void populate_pt_basis(
    const IKstatePopulator<KstateT>& populator,
    const unsigned max_pt_order,
    kstate::Basis<KstateT>& basis,
    unsigned n_threads = 1) {
    static_assert(!std::is_array_v<KstateT>);
    static_assert(!std::is_function_v<KstateT>);
    static_assert(!std::is_void_v<std::decay<KstateT>>);
    static_assert(!std::is_null_pointer_v<std::decay<KstateT>>);
    static_assert(!std::is_enum_v<std::decay<KstateT>>);
    static_assert(!std::is_union_v<std::decay<KstateT>>);
    static_assert(std::is_class_v<std::decay<KstateT>>);
    static_assert(!std::is_pointer_v<std::decay<KstateT>>);
    static_assert(!std::is_member_object_pointer_v<KstateT>);
    static_assert(!std::is_member_function_pointer_v<KstateT>);
    static_assert(!std::is_const_v<KstateT>);
    static_assert(!std::is_volatile_v<KstateT>);
    static_assert(!std::is_reference_v<KstateT>);
    static_assert(kstate::is_base_of_template_v<KstateT, kstate::Kstate>);
    unsigned last_chunk_size = basis.size();
    for (unsigned pt_order = 0; pt_order < max_pt_order; pt_order++) {
        const unsigned old_size = basis.size();
        assert(last_chunk_size <= old_size);
#pragma omp parallel for num_threads(n_threads)
        for (unsigned idx = old_size - last_chunk_size; idx < old_size; idx++) {
            const auto el = basis.vec_index()[idx];
            assert(el);
            populator.push_back_coupled_states_to_basis(*el, basis);
        }
        const unsigned new_size = basis.size();
        last_chunk_size = new_size - old_size;
    }
}

}  // namespace bfpt_common

#endif
