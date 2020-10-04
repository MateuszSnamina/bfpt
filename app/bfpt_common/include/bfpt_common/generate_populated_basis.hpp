#ifndef BFPT_COMMON_GENERATE_POPULATED_BASIS_HPP
#define BFPT_COMMON_GENERATE_POPULATED_BASIS_HPP

#include <bfpt_common/i_kstate_basis_populator.hpp>

#include <kstate/basis.hpp>
#include <kstate/kstate_stl.hpp>
#include <kstate/kstate_abstract.hpp>

#include <omp.h>

#include <cassert>

#include <set>
//#include <chrono> // performance debug sake

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

template<typename KstateTrait>
void generate_populated_basis(
        const IKstateBasisPopulator<KstateTrait>& basis_populator,
        const unsigned max_pt_order,
        kstate::Basis<KstateTrait>& basis,
        unsigned n_threads = 1) {
//    static_assert(!std::is_array_v<KstateT>);
//    static_assert(!std::is_function_v<KstateT>);
//    static_assert(!std::is_void_v<std::decay<KstateT>>);
//    static_assert(!std::is_null_pointer_v<std::decay<KstateT>>);
//    static_assert(!std::is_enum_v<std::decay<KstateT>>);
//    static_assert(!std::is_union_v<std::decay<KstateT>>);
//    static_assert(std::is_class_v<std::decay<KstateT>>);
//    static_assert(!std::is_pointer_v<std::decay<KstateT>>);
//    static_assert(!std::is_member_object_pointer_v<KstateT>);
//    static_assert(!std::is_member_function_pointer_v<KstateT>);
//    static_assert(!std::is_const_v<KstateT>);
//    static_assert(!std::is_volatile_v<KstateT>);
//    static_assert(!std::is_reference_v<KstateT>);
//    static_assert(kstate::is_base_of_template_v<KstateT, kstate::Kstate>);//TODO remove.
    // *********** asserts ********************************************************************
    static_assert(KstateTrait::is_kstate_trait);
    // *********** using **********************************************************************
    using KstateT = typename KstateTrait::KstateT;
    // *********** pt orders loop  ************************************************************
    unsigned last_chunk_size = basis.size();
    for (unsigned pt_order = 0; pt_order < max_pt_order; pt_order++) {
        const unsigned old_size = basis.size();
        assert(last_chunk_size <= old_size);
        // *********** prepare ****************************************************************
        std::vector<kstate::KstateSet<KstateT>> kstate_set_all(n_threads);
        // *********** filling ****************************************************************
        //const auto tp_fill_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
#pragma omp parallel for num_threads(n_threads)
        for (unsigned idx = old_size - last_chunk_size; idx < old_size; idx++) {
            const auto tid = omp_get_thread_num();
            const auto el = basis.vec_index()[idx];
            assert(el);
            const kstate::KstateSet<KstateT> newely_generated_states = basis_populator.get_coupled_states(*el);
            kstate_set_all[tid].insert(std::begin(newely_generated_states), std::end(newely_generated_states));
        }
        //const auto tp_fill_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
        //std::cout << "fill took     :" << std::chrono::duration_cast<std::chrono::nanoseconds>(tp_fill_2 - tp_fill_1).count() / 1e6 << "ms" << std::endl; // performance debug sake
        // *********** reduction **************************************************************
        //const auto tp_reduce_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
        for (unsigned d = 1; d < n_threads; d *=2) {
#pragma omp parallel num_threads(n_threads)
            {
                const auto tid = omp_get_thread_num();
                if (tid % (2 * d) == 0) {
                    const auto idx1 = tid;
                    const auto idx2 = tid + d;
                    if (idx2 < n_threads) {
                        kstate_set_all[idx1].insert(std::cbegin(kstate_set_all[idx2]), std::cend(kstate_set_all[idx2]));
                    }
                }
            }
        }
        //const auto tp_reduce_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
        //std::cout << "reduce took   :" << std::chrono::duration_cast<std::chrono::nanoseconds>(tp_reduce_2 - tp_reduce_1).count() / 1e6 << "ms" << std::endl; // performance debug sake
        // *********** bassis append **********************************************************
        for (const auto& new_element : kstate_set_all[0]) {
            basis.add_element(new_element);
        }
        // *********** update last_chunk_size *************************************************
        const unsigned new_size = basis.size();
        last_chunk_size = new_size - old_size;
    }
}

}  // namespace bfpt_common

#endif
