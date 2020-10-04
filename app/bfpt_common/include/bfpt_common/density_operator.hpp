#ifndef BFPT_COMMON_DENSITY_OPERATOR_HPP
#define BFPT_COMMON_DENSITY_OPERATOR_HPP

#include <bfpt_common/operator_kernel.hpp>

#include <map>
#include <complex>

// #######################################################################
// ## DensityOperator{1,12}                                             ##
// #######################################################################

namespace bfpt_common {

template<typename SiteStateTraitT>
using DensityOperator1 = std::map<std::pair<StateKernel1<SiteStateTraitT>, StateKernel1<SiteStateTraitT>>, std::complex<double>>;

template<typename SiteStateTraitT>
using DensityOperator12 = std::map<std::pair<StateKernel12<SiteStateTraitT>, StateKernel12<SiteStateTraitT>>, std::complex<double>>;

} // end of namespace bfpt_common

#endif
