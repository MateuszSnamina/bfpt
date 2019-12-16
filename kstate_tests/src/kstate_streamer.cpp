#include <kstate/kstate.hpp>
#include <kstate/kstate_streamer.hpp>

#include <gtest/gtest.h>

using kstate::ctr_from_range;

TEST(KStateStreamer, SixElements) {
    const int v1[6] = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    ASSERT_EQ(kstate::KstateStringStreamer().stream(k1).str(), "𝕂𝕊⦃11∙12∙13∙14∙15∙16⦄");
}

TEST(KStateStreamer, OneElement) {
    const int v1[1] = {11};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_EQ(kstate::KstateStringStreamer().stream(k1).str(), "𝕂𝕊⦃11⦄");
}

TEST(KStateStreamer, Empty) {
    std::vector<int> v1;
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_EQ(kstate::KstateStringStreamer().stream(k1).str(), "𝕂𝕊⦃⦄");
}

TEST(KStateStreamer, Fancy) {
    double v1[3] = {1.1, 1.2, 1.3};
    const kstate::DynamicKstate<double> k1(v1, ctr_from_range);
    const auto s = kstate::KstateStringStreamer()
                       .set_stream_preparer([](std::ostream& s) {
                           s << "BEGIN***|" << std::fixed << std::setprecision(2)
                             << std::showpos;
                       })
                       .set_stream_sustainer([](std::ostream& s, size_t i) {
                           s << std::setw(2) << char('A' + i) << ":" << std::setw(6);
                       })
                       .set_stream_separer([](std::ostream& s) { s << "|"; })
                       .set_stream_finisher([](std::ostream& s) { s << "|***END"; })
                       .stream(k1)
                       .str();
    ASSERT_EQ(s, "BEGIN***| A: +1.10| B: +1.20| C: +1.30|***END");
}