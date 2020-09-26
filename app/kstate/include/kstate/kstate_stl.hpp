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

template<class KstateT>
using KstateSet = std::set<std::shared_ptr<KstateT>, KstateComparator<KstateT>>;

template<class KstateT, typename T>
using KstateMap = std::map<std::shared_ptr<KstateT>, T, KstateComparator<KstateT>>;

}


#endif // KSTATE_KSTATE_STL_HPP
