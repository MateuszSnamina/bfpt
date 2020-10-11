#pragma once

#include <type_traits>

namespace koperator_trait {

template<typename _KoperatorT>
struct TraitKoperator {
    static constexpr bool is_koperator_trait = false;
    using KoperatorT = _KoperatorT;
};

template<typename T>
struct IsTraitKoperator : std::false_type {
};

template<typename T>
struct IsTraitKoperator<TraitKoperator<T>> : std::true_type {
};

} // end of namespace kstate
