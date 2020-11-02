#pragma once

#include <kstate_op_range/op_range_compare.hpp>
#include <kstate_op_range/op_range_least_replication_shift.hpp>

#include <kstate_trait/trait_site_state.hpp>

#include <extensions/range_streamer.hpp>

#include <boost/range.hpp>

#include <cassert>
#include <iterator>
#include <optional>
#include <string>
#include <type_traits>

// #######################################################################
// ## Kstate                                                           ##
// #######################################################################

namespace kstate_impl {

template <typename _SiteStateTraitT, typename _ConstRangeT>
class Kstate {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(std::is_same_v<typename _SiteStateTraitT::SiteStateT, typename boost::range_value<_ConstRangeT>::type>);

   public:  // Helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using ConstRangeT = _ConstRangeT;
    using SiteStateT = typename boost::range_value<ConstRangeT>::type;
    using TraversalTagT = typename boost::range_traversal<ConstRangeT>::type;
    //static_assert(std::is_same<TraversalTagT, boost::random_access_traversal_tag>::value ||
    //std::is_same<TraversalTagT, boost::forward_traversal_tag>::value);
   public:  // Base on the two functions:
    virtual ConstRangeT to_range() const noexcept = 0;
    virtual size_t n_sites() const noexcept = 0;

   public:  // instance spec functions:
    virtual size_t n_least_replication_shift() const noexcept;
    virtual double norm_factor() const noexcept;
    virtual bool is_prolific(int n_k) const noexcept;

   public:  // compare_XXX functions:
    template <typename OtherConstRangeT>
    bool compare_equality_range(const OtherConstRangeT& other) const noexcept;
    template <typename OtherConstRangeT>
    std::optional<size_t> compare_translational_equality_range(const OtherConstRangeT& other) const noexcept;
    virtual std::string to_str() const noexcept;

   public:
    virtual ~Kstate() = default;
};

// ***********************************************************************

template <typename _SiteStateTraitT, typename _ConstRangeT>
size_t
Kstate<_SiteStateTraitT, _ConstRangeT>::n_least_replication_shift() const noexcept {
    return kstate_op_range::n_least_replication_shift(to_range());
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
double
Kstate<_SiteStateTraitT, _ConstRangeT>::norm_factor() const noexcept {
    return kstate_op_range::norm_factor(to_range());
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
bool Kstate<_SiteStateTraitT, _ConstRangeT>::is_prolific(int n_k) const noexcept {
    return kstate_op_range::is_prolific(to_range(), n_k);
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
template <typename OtherConstRangeT>
bool Kstate<_SiteStateTraitT, _ConstRangeT>::compare_equality_range(const OtherConstRangeT& other) const noexcept {
    assert(this->n_sites() == boost::size(other));
    return kstate_op_range::compare_equality(to_range(), other);
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
template <typename OtherConstRangeT>
std::optional<size_t>
Kstate<_SiteStateTraitT, _ConstRangeT>::compare_translational_equality_range(const OtherConstRangeT& other) const noexcept {
    //assert(this->n_sites() == boost::size(other));
    return kstate_op_range::compare_translational_equality(to_range(), other);
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
std::string
Kstate<_SiteStateTraitT, _ConstRangeT>::to_str() const noexcept {
    using namespace extension::boost::stream_pragma;
    const auto range_stream_settings = RSS<SiteStateT>()
                                           .set_string_preparer("⦃")
                                           .set_null_sustainer()
                                           .set_string_separer("∙")
                                           .set_string_finisher("⦄");
    return (to_range() | range_stream_settings).str();
}

}  // namespace kstate_impl
