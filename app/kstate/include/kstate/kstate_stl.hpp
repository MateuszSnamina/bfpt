#ifndef KSTATE_KSTATE_STL_HPP
#define KSTATE_KSTATE_STL_HPP

#include<kstate/kstate_comparator.hpp>

#include<set>
#include<map>
#include<memory>

// #######################################################################
// ## Kstate{Set,Map}                                                   ##
// #######################################################################

namespace kstate {

template<class KstateTraitT>
using KstateSet = std::set<std::shared_ptr<typename KstateTraitT::KstateT>, KstateComparator<KstateTraitT>>;

template<class KstateTraitT, typename T>
using KstateMap = std::map<std::shared_ptr<typename KstateTraitT::KstateT>, T, KstateComparator<KstateTraitT>>;

}


#endif // KSTATE_KSTATE_STL_HPP
