#include <kstate_trait/trait_site_state.hpp>

#include <kstate_impl/kstate_concrete_integral.hpp>

#include <gtest/gtest.h>

#include <iostream>
#include <array>

// #######################################################################
// ## MonostarSiteState -- class                                        ##
// #######################################################################

namespace {

class MonostarForTestsSiteState {
public:
    MonostarForTestsSiteState(bool is_excited) : is_excited(is_excited){};
    const bool is_excited;
};

const MonostarForTestsSiteState gs{false};  // ground-state (there is no star)
const MonostarForTestsSiteState es{true};  // excited-state (there is a star)

}  // namespace

// #######################################################################
// ## MonostarSiteState -- global values                                ##
// #######################################################################

namespace  {


}  // namespace

// #######################################################################
// ## MonostarSiteState -- implement trait                              ##
// #######################################################################

namespace kstate_trait {

template<>
struct TraitSiteState<MonostarForTestsSiteState> {
    static constexpr bool is_site_state_trait = true;
    using SiteStateT = MonostarForTestsSiteState;

    static unsigned site_basis_dim() {
        return 2u;
    }

    static unsigned get_index(const SiteStateT& state) {
        return (state.is_excited ? 1u : 0u);
    }

    static SiteStateT from_index(unsigned idx) {
        if (idx == 0u) {
            return gs;
        } else if (idx == 1u) {
            return es;
        } else {
            throw std::domain_error("Index out of range.");
        }
    }
};

} // end of namespace kstate

namespace {

using MonostarSiteStateTrait = kstate_trait::TraitSiteState<MonostarForTestsSiteState>;

} // end of namespace

// #######################################################################
// ## Tests                                                             ##
// #######################################################################

TEST(DynamicIntegralKstate, IntegralToSiteStateRange) {
    const auto r = kstate_impl::integral_to_site_state_range<MonostarSiteStateTrait, uint64_t>(0b01101, 5u);
    EXPECT_EQ(boost::size(r), 5);
    EXPECT_TRUE((*std::next(std::begin(r), 0)).is_excited);
    EXPECT_FALSE((*std::next(std::begin(r), 1)).is_excited);
    EXPECT_TRUE((*std::next(std::begin(r), 2)).is_excited);
    EXPECT_TRUE((*std::next(std::begin(r), 3)).is_excited);
    EXPECT_FALSE((*std::next(std::begin(r), 4)).is_excited);
    EXPECT_TRUE(std::next(std::begin(r), 0)->is_excited);
    EXPECT_FALSE(std::next(std::begin(r), 1)->is_excited);
    EXPECT_TRUE(std::next(std::begin(r), 2)->is_excited);
    EXPECT_TRUE(std::next(std::begin(r), 3)->is_excited);
    EXPECT_FALSE(std::next(std::begin(r), 4)->is_excited);
}

TEST(DynamicIntegralKstate, SiteStateRangeToIntegral) {
    std::array<MonostarForTestsSiteState, 5> siste_state_array{es, gs, gs, es, es};
    kstate_op_integral::IntegralBitsDynamic<uint64_t> integra_bits =
            kstate_impl::site_state_range_to_integral<MonostarSiteStateTrait, uint64_t>(siste_state_array);
    EXPECT_EQ(integra_bits.get_n_all_bits() , 5);
    EXPECT_EQ(integra_bits.get_buffer() , 0b11001);
}
