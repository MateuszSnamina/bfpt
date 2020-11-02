#pragma once

#include <type_traits>

namespace koperator_trait {

template <typename _KoperatorT>
struct TraitKoperator {
    static constexpr bool is_koperator_trait = false;
    using KoperatorT = _KoperatorT;
};

template <typename T>
struct IsTraitKoperator : std::false_type {
};

template <typename T>
struct IsTraitKoperator<TraitKoperator<T>> : std::true_type {
};

// #######################################################################
// ## Instruction: registering a type as Koperator                      ##
// #######################################################################

/*
 *
 * In order to register a type Foo as a Koperator you need to
 * specialize TraitKoperator that will serve an API for the type
 * to be used in templates.
 *
 * The snippet below shows how to specialize TraitKoperator
 * providing all the necessary field/types.
 *
 * -----------------------------------------------------
 *
 * namespace koperator_trait {
 *
 * template<typename foo_namespace::_Foo>
 * struct TraitKoperator<foo_namespace::_Foo> {
 *     // the is_kstate_trait flag:
 *     static constexpr bool is_koperator_trait = true;
 *     // helper types:
 *     using KoperatorT = foo_namespace::_Foo;
 *     using KstateTraitT = <HERE FILL THE RIGHT TYPE>;
 *     using BasisT = kbasis::Basis<KstateTraitT>;
 *     // function being the public API:
 *     static void fill_kn_operator_builder_matrix_coll(
 *             const KoperatorT& koperator,
 *             const BasisT& basis,
 *             size_t n_col,
 *             arma::sp_cx_mat& kn_operator_builder_matrix,
 *             const unsigned k_n) {
 *             <HERE FILL THE IMPLEMENTATION TYPE>
 *             }
 * }; // end of struct TraitKoperator
 *
 * }  // end of namespace koperator_trait
 *
 * -----------------------------------------------------
 *
 * TraitKoperator is an inteface for kstate operator matrix builder.
 * Technically, the operator is build from matrix `G`
 * such as `H = G + G^\dagger`.
 *
 * The class provides `fill_kn_operator_builder_matrix_coll` function
 * for filling a single column in `G` matrix.
 * (The arg name `kn_operator_builder_matrix` is for the matrix G).
 *
 */

}  // namespace koperator_trait
