#include <extensions/range_streamer.hpp>

#include <iomanip>
#include <vector>

#include <gtest/gtest.h>

TEST(RangeStreamer, RawRangeStreamerSixElementRange) {
  int v1[6] = {13, 14, 15, 16, 11, 12};
  ASSERT_EQ(extension::boost::RangeStringStreamer().stream(v1).str(),
            "{0:13,1:14,2:15,3:16,4:11,5:12}");
}

TEST(RangeStreamer, RawRangeStreamerOneElementRange) {
  int v2[1] = {13};
  ASSERT_EQ(extension::boost::RangeStringStreamer().stream(v2).str(), "{0:13}");
}

TEST(RangeStreamer, RawRangeStreamerEmptyRange) {
  std::vector<int> v3;
  ASSERT_EQ(extension::boost::RangeStringStreamer().stream(v3).str(), "{}");
}

TEST(RangeStreamer, Fancy) {
  double v1[3] = {1.1, 1.2, 1.3};
  const auto s =
      extension::boost::RangeStringStreamer()
          .set_stream_preparer([](std::ostream& s) {
            s << "BEGIN***|" << std::fixed << std::setprecision(2)
              << std::showpos;
          })
          .set_stream_sustainer([](std::ostream& s, size_t i) {
            s << std::setw(2) << char('A' + i) << ":" << std::setw(6);
          })
          .set_stream_separer([](std::ostream& s) { s << "|"; })
          .set_stream_finisher([](std::ostream& s) { s << "|***END"; })
          .stream(v1)
          .str();
  ASSERT_EQ(s, "BEGIN***| A: +1.10| B: +1.20| C: +1.30|***END");
}