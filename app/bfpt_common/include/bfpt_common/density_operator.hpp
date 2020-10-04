#ifndef BFPT_COMMON_DENSITY_OPERATOR_HPP
#define BFPT_COMMON_DENSITY_OPERATOR_HPP

#include <bfpt_common/operator_kernel.hpp>

#include <map>
#include <complex>

// #######################################################################
// ## DensityOperator{1,12}                                             ##
// #######################################################################

namespace bfpt_common {

template<typename SiteState>
using DensityOperator1 = std::map<std::pair<StateKernel1<SiteState>, StateKernel1<SiteState>>, std::complex<double>>;

template<typename SiteStateT>
using DensityOperator12 = std::map<std::pair<StateKernel12<SiteStateT>, StateKernel12<SiteStateT>>, std::complex<double>>;

} // end of namespace bfpt_common

#endif
