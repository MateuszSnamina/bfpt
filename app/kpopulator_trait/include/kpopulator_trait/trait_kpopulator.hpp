#pragma once

#include <type_traits>
#include <cassert>

namespace kpopulator_trait {

template <typename _KpopulatorT>
struct TraitKpopulator {
    static constexpr bool is_kpopulator_trait = false;
    using KpopulatorT = _KpopulatorT;
};

template <typename T>
struct IsTraitKpopulator : std::false_type {
};

template <typename T>
struct IsTraitKpopulator<TraitKpopulator<T>> : std::true_type {
};

}  // namespace kpopulator_trait

// #######################################################################
// ## Instruction: registering a type as Kpopulator                     ##
// #######################################################################

//[TODO]
