#pragma once

#include <chainkernel/operator_kernel.hpp>

#include <map>
#include <complex>

// #######################################################################
// ## DensityOperator{1,12}                                             ##
// #######################################################################


//TODO remove the file!!!!!

namespace bfpt_common {

template<typename SiteStateTraitT>
using DensityOperator1 = std::map<std::pair<chainkernel::StateKernel1<SiteStateTraitT>, chainkernel::StateKernel1<SiteStateTraitT>>, std::complex<double>>;

template<typename SiteStateTraitT>
using DensityOperator12 = std::map<std::pair<chainkernel::StateKernel12<SiteStateTraitT>, chainkernel::StateKernel12<SiteStateTraitT>>, std::complex<double>>;

} // end of namespace bfpt_common
