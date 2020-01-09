// #ifndef KSTATE_BASIS_POPULATE_HPP
// #define KSTATE_BASIS_POPULATE_HPP

// #include <kstate/basis.hpp>
// #include <kstate/kstate_concrete.hpp>

// // #######################################################################
// // ##  IDynamicUniqueKstatePopulator                                    ##
// // #######################################################################

// namespace kstate {

// template <typename Element>
// class IDynamicUniqueKstatePopulator {
//    public:
//     virtual void push_back_coupled_states_to_basis(
//         const DynamicUniqueKstate<Element>& generator,
//         Basis<DynamicUniqueKstate<Element>>& basis) const = 0;
// };

// }  // namespace kstate

// // #######################################################################
// // ##  populate_pt_basis                                                ##
// // #######################################################################

// /*
//  * The function generates pt-basis of given order.
//  * The input for the function is:
//  *  - the desider pt-order,
//  *  - the 0'th pt-order basis,
//  *  - the populator object.
//  * 
//  * Desired pt-order is given as `max_pt_order` argument;
//  * the 0'th pt-order basis is read from `basis` argument
//  * (the argument has intent "inout",
//  *  as input it gives the 0'th pt-order basis
//  *  as output it gives the function result);
//  * Populator object is given as `populator` arggument.
//  * 
//  */

// namespace kstate {

// template <typename Element>
// void populate_pt_basis(
//     const IDynamicUniqueKstatePopulator<Element>& populator,
//     const unsigned max_pt_order,
//     Basis<DynamicUniqueKstate<Element>>& basis) {
//     unsigned last_chunk_size = basis.size();
//     for (unsigned pt_order = 0; pt_order < max_pt_order; pt_order++) {
//         const unsigned old_size = basis.size();
//         assert(last_chunk_size <= old_size);
//         for (unsigned idx = old_size - last_chunk_size; idx < old_size; idx++) {
//             const auto el = basis.vec_index()[idx];
//             assert(el);
//             populator.push_back_coupled_states_to_basis(*el, basis);
//         }
//         const unsigned new_size = basis.size();
//         last_chunk_size = new_size - old_size;
//     }
// }

// }  // namespace kstate

// #endif