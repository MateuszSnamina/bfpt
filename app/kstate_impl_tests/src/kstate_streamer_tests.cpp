#include <kstate_tests/site_state_trait_for_int.hpp>
#include <kstate_tests/site_state_trait_for_double.hpp>

#include <kstate_impl/kstate_streamer.hpp>
#include <kstate_impl/kstate_concrete_stl.hpp>

#include <gtest/gtest.h>

using kstate_impl::ctr_from_range;

// #######################################################################
// ## KstateStringStreamer                                              ##
// #######################################################################

TEST(KstateStringStreamer, SixElements) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using extension::boost::stream_pragma::RSS;
    using kstate_impl::pramga::operator||;
    using kstate_impl::pramga::operator<<;
    using kstate_impl::pramga::operator||;
    const int v1[6] = {11, 12, 13, 14, 15, 16};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ((k1 || RSS<int>()).str(), "𝕂𝕤𝕥𝕒𝕥𝕖⦃11∙12∙13∙14∙15∙16⦄");
}

TEST(KstateStringStreamer, OneElement) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using extension::boost::stream_pragma::RSS;
    using kstate_impl::pramga::operator||;
    using kstate_impl::pramga::operator<<;
    const int v1[1] = {11};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ((k1 || RSS<int>()).str(), "𝕂𝕤𝕥𝕒𝕥𝕖⦃11⦄");
}

TEST(KstateStringStreamer, Empty) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using extension::boost::stream_pragma::RSS;
    using kstate_impl::pramga::operator||;
    using kstate_impl::pramga::operator<<;
    std::vector<int> v1;
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ((k1 || RSS<int>()).str(), "𝕂𝕤𝕥𝕒𝕥𝕖⦃⦄");
}

TEST(KstateStringStreamer, Fancy) {
    using SiteStateTrait = kstate_trait::TraitSiteState<double>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using extension::boost::stream_pragma::RSS;
    using kstate_impl::pramga::operator||;
    using kstate_impl::pramga::operator<<;
    double v1[3] = {1.1, 1.2, 1.3};
    const Kstate k1(v1, ctr_from_range);
    const auto kstate_stream_settings = RSS<double>()
                       .set_stream_preparer([](std::ostream& s) {
                           s << "BEGIN***|" << std::fixed << std::setprecision(2)
                             << std::showpos;
                       })
                       .set_stream_sustainer([](std::ostream& s, size_t i) {
                           s << std::setw(2) << char('A' + i) << ":" << std::setw(6);
                       })
                       .set_stream_separer([](std::ostream& s) { s << "|"; })
                       .set_stream_finisher([](std::ostream& s) { s << "|***END"; });
    ASSERT_EQ((k1 || kstate_stream_settings).str(), "BEGIN***| A: +1.10| B: +1.20| C: +1.30|***END");
}
