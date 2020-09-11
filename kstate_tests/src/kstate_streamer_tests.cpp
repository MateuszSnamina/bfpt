#include <kstate/kstate_streamer.hpp>

#include <kstate/kstate_concrete.hpp>

#include <gtest/gtest.h>

using kstate::ctr_from_range;

// #######################################################################
// ## KstateStringStreamer                                              ##
// #######################################################################

TEST(KstateStringStreamer, SixElements) {
    using namespace kstate::pramga;
    const int v1[6] = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    ASSERT_EQ((k1 | KSSS()).str(), "𝕂𝕤𝕥𝕒𝕥𝕖⦃11∙12∙13∙14∙15∙16⦄");
}

TEST(KstateStringStreamer, OneElement) {
    using namespace kstate::pramga;
    const int v1[1] = {11};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    ASSERT_EQ((k1 | KSSS()).str(), "𝕂𝕤𝕥𝕒𝕥𝕖⦃11⦄");
}

TEST(KstateStringStreamer, Empty) {
    using namespace kstate::pramga;
    std::vector<int> v1;
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    ASSERT_EQ((k1 | KSSS()).str(), "𝕂𝕤𝕥𝕒𝕥𝕖⦃⦄");
}

TEST(KstateStringStreamer, Fancy) {
    using namespace kstate::pramga;
    double v1[3] = {1.1, 1.2, 1.3};
    const kstate::DynamicKstate<double> k1(v1, ctr_from_range);
    const auto kstate_stream_settings = KSSS()
                       .set_stream_preparer([](std::ostream& s) {
                           s << "BEGIN***|" << std::fixed << std::setprecision(2)
                             << std::showpos;
                       })
                       .set_stream_sustainer([](std::ostream& s, size_t i) {
                           s << std::setw(2) << char('A' + i) << ":" << std::setw(6);
                       })
                       .set_stream_separer([](std::ostream& s) { s << "|"; })
                       .set_stream_finisher([](std::ostream& s) { s << "|***END"; });
    ASSERT_EQ((k1 | kstate_stream_settings).str(), "BEGIN***| A: +1.10| B: +1.20| C: +1.30|***END");
}
