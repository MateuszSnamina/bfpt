#pragma once

#include<kstate_trait/kstate_comparator.hpp>

#include<set>
#include<map>
#include<memory>

// #######################################################################
// ## Kstate{Set,Map}                                                   ##
// #######################################################################

namespace kstate_trait {

template<class KstateTraitT>
using KstateSet = std::set<std::shared_ptr<typename KstateTraitT::KstateT>, KstateComparator<KstateTraitT>>;

template<class KstateTraitT, typename T>
using KstateMap = std::map<std::shared_ptr<typename KstateTraitT::KstateT>, T, KstateComparator<KstateTraitT>>;

}
