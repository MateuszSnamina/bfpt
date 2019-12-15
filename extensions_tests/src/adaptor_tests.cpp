#include <extensions/adaptors.hpp>

#include <array>
#include <list>
#include <vector>

#include <boost/range/algorithm.hpp>

#include <gtest/gtest.h>

TEST(ExtensionBoostAdaptorsRotated, FilledRawArray) {
    int v1[6] = {11, 12, 13, 14, 15, 16};
    int v2[6] = {13, 14, 15, 16, 11, 12};
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(0)).size(), 6);
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(2)).size(), 6);
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(6)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(2), v2));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(6), v1));
}

TEST(ExtensionBoostAdaptorsRotated, FilledStdArray) {
    std::array<int, 6> v1{11, 12, 13, 14, 15, 16};
    std::array<int, 6> v2{13, 14, 15, 16, 11, 12};
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(0)).size(), 6);
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(2)).size(), 6);
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(6)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(2), v2));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(6), v1));
}

TEST(ExtensionBoostAdaptorsRotated, FilledStdVector) {
    std::vector<int> v1{11, 12, 13, 14, 15, 16};
    std::vector<int> v2{13, 14, 15, 16, 11, 12};
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(0)).size(), 6);
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(2)).size(), 6);
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(6)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(2), v2));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(6), v1));
}

TEST(ExtensionBoostAdaptorsRotated, FilledStdList) {
    std::list<int> v1{11, 12, 13, 14, 15, 16};
    std::list<int> v2{13, 14, 15, 16, 11, 12};
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(2), v2));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(6), v1));
}

TEST(ExtensionBoostAdaptorsRotated, EmptyStdList) {
    std::list<int> v1;
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
}

TEST(ExtensionBoostAdaptorsDoubled, FilledRawArray) {
    int v1[6] = {11, 12, 13, 14, 15, 16};
    int v2[12] = {11, 12, 13, 14, 15, 16, 11, 12, 13, 14, 15, 16};
    ASSERT_EQ((v1 | extension::boost::adaptors::doubled).size(), 12);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::doubled, v2));
}