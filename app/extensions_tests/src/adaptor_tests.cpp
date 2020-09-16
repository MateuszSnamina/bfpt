#include <extensions/adaptors.hpp>

#include <array>
#include <list>
#include <vector>

#include <boost/range/algorithm.hpp>

#include <gtest/gtest.h>

// #######################################################################
// ##  rotated                                                          ##
// #######################################################################

TEST(ExtensionBoostAdaptorsRotated, FilledRawArray) {
    int v1[6] = {11, 12, 13, 14, 15, 16};
    int v2[6] = {13, 14, 15, 16, 11, 12};
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(0)).size(), 6);
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(2)).size(), 6);
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(6)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(2), v2));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(6), v1));
    const auto r = v1 | extension::boost::adaptors::rotated(2);
    ASSERT_TRUE(boost::equal(r, v2));
    ASSERT_TRUE(boost::equal(r, v2));
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
    const auto r = v1 | extension::boost::adaptors::rotated(2);
    ASSERT_TRUE(boost::equal(r, v2));
    ASSERT_TRUE(boost::equal(r, v2));
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
    const auto r = v1 | extension::boost::adaptors::rotated(2);
    ASSERT_TRUE(boost::equal(r, v2));
    ASSERT_TRUE(boost::equal(r, v2));
}

TEST(ExtensionBoostAdaptorsRotated, FilledStdList) {
    std::list<int> v1{11, 12, 13, 14, 15, 16};
    std::list<int> v2{13, 14, 15, 16, 11, 12};
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::rotated(0)), 6);
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::rotated(2)), 6);
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::rotated(6)), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(2), v2));
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(6), v1));
    const auto r = v1 | extension::boost::adaptors::rotated(2);
    ASSERT_TRUE(boost::equal(r, v2));
    ASSERT_TRUE(boost::equal(r, v2));
}

TEST(ExtensionBoostAdaptorsRotated, EmptyStdArray) {
    std::array<int, 0> v1;
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(0)).size(), 0);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
    const auto r = v1 | extension::boost::adaptors::rotated(0);
    ASSERT_TRUE(boost::equal(r, v1));
    ASSERT_TRUE(boost::equal(r, v1));
}

TEST(ExtensionBoostAdaptorsRotated, EmptyStdVector) {
    std::vector<int> v1;
    ASSERT_EQ((v1 | extension::boost::adaptors::rotated(0)).size(), 0);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
    const auto r = v1 | extension::boost::adaptors::rotated(0);
    ASSERT_TRUE(boost::equal(r, v1));
    ASSERT_TRUE(boost::equal(r, v1));
}

TEST(ExtensionBoostAdaptorsRotated, EmptyStdList) {
    std::list<int> v1;
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::rotated(0)), 0);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::rotated(0), v1));
    const auto r = v1 | extension::boost::adaptors::rotated(0);
    ASSERT_TRUE(boost::equal(r, v1));
    ASSERT_TRUE(boost::equal(r, v1));
}
// #######################################################################
// ##  doubled                                                          ##
// #######################################################################

TEST(ExtensionBoostAdaptorsDoubled, FilledRawArray) {
    int v1[6] = {11, 12, 13, 14, 15, 16};
    int v2[12] = {11, 12, 13, 14, 15, 16, 11, 12, 13, 14, 15, 16};
    ASSERT_EQ((v1 | extension::boost::adaptors::doubled).size(), 12);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::doubled, v2));
    const auto r = v1 | extension::boost::adaptors::doubled;
    ASSERT_TRUE(boost::equal(r, v2));
    ASSERT_TRUE(boost::equal(r, v2));
}

TEST(ExtensionBoostAdaptorsDoubled, FilledStdArray) {
    std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    std::array<int, 12> v2 = {11, 12, 13, 14, 15, 16, 11, 12, 13, 14, 15, 16};
    ASSERT_EQ((v1 | extension::boost::adaptors::doubled).size(), 12);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::doubled, v2));
    const auto r = v1 | extension::boost::adaptors::doubled;
    ASSERT_TRUE(boost::equal(r, v2));
    ASSERT_TRUE(boost::equal(r, v2));
}

TEST(ExtensionBoostAdaptorsDoubled, FilledStdVector) {
    std::vector<int> v1 = {11, 12, 13, 14, 15, 16};
    std::vector<int> v2 = {11, 12, 13, 14, 15, 16, 11, 12, 13, 14, 15, 16};
    ASSERT_EQ((v1 | extension::boost::adaptors::doubled).size(), 12);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::doubled, v2));
    const auto r = v1 | extension::boost::adaptors::doubled;
    ASSERT_TRUE(boost::equal(r, v2));
    ASSERT_TRUE(boost::equal(r, v2));
}

TEST(ExtensionBoostAdaptorsDoubled, FilledStdList) {
    std::list<int> v1 = {11, 12, 13, 14, 15, 16};
    std::list<int> v2 = {11, 12, 13, 14, 15, 16, 11, 12, 13, 14, 15, 16};
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::doubled), 12);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::doubled, v2));
    const auto r = v1 | extension::boost::adaptors::doubled;
    ASSERT_TRUE(boost::equal(r, v2));
    ASSERT_TRUE(boost::equal(r, v2));
}

TEST(ExtensionBoostAdaptorsDoubled, EmptyStdArray) {
    std::array<int, 0> v1 = {};
    ASSERT_EQ((v1 | extension::boost::adaptors::doubled).size(), 0);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::doubled, v1));
    const auto r = v1 | extension::boost::adaptors::doubled;
    ASSERT_TRUE(boost::equal(r, v1));
    ASSERT_TRUE(boost::equal(r, v1));
}

TEST(ExtensionBoostAdaptorsDoubled, EmptyStdVector) {
    std::vector<int> v1 = {};
    ASSERT_EQ((v1 | extension::boost::adaptors::doubled).size(), 0);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::doubled, v1));
    const auto r = v1 | extension::boost::adaptors::doubled;
    ASSERT_TRUE(boost::equal(r, v1));
    ASSERT_TRUE(boost::equal(r, v1));
}

TEST(ExtensionBoostAdaptorsDoubled, EmptyStdList) {
    std::list<int> v1 = {};
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::doubled), 0);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::doubled, v1));
    const auto r = v1 | extension::boost::adaptors::doubled;
    ASSERT_TRUE(boost::equal(r, v1));
    ASSERT_TRUE(boost::equal(r, v1));
}

// #######################################################################
// ##  refined                                                          ##
// #######################################################################

TEST(ExtensionBoostAdaptorsRefined, FilledRawArray) {
    int v1[6] = {11, 12, 13, 14, 15, 16};
    int v2[6] = {11, 12, 47, 14, 15, 16};
    int v3[6] = {11, 12, 47, 14, 48, 16};
    int v4[6] = {47, 12, 13, 14, 15, 16};
    int v5[6] = {47, 12, 13, 14, 48, 16};
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(2, 47)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(2, 47), v2));
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(2, 47) | extension::boost::adaptors::refined(4, 48)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(2, 47) | extension::boost::adaptors::refined(4, 48), v3));
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(0, 47)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47), v4));
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48), v5));
    auto r = v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48);
    ASSERT_TRUE(boost::equal(r, v5));
    ASSERT_TRUE(boost::equal(r, v5));
}

TEST(ExtensionBoostAdaptorsRefined, FilledStdArray) {
    std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    std::array<int, 6> v2 = {11, 12, 47, 14, 15, 16};
    std::array<int, 6> v3 = {11, 12, 47, 14, 48, 16};
    std::array<int, 6> v4 = {47, 12, 13, 14, 15, 16};
    std::array<int, 6> v5 = {47, 12, 13, 14, 48, 16};
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(2, 47)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(2, 47), v2));
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(2, 47) | extension::boost::adaptors::refined(4, 48)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(2, 47) | extension::boost::adaptors::refined(4, 48), v3));
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(0, 47)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47), v4));
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48), v5));
    auto r = v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48);
    ASSERT_TRUE(boost::equal(r, v5));
    ASSERT_TRUE(boost::equal(r, v5));
}

TEST(ExtensionBoostAdaptorsRefined, FilledStdVector) {
    std::vector<int> v1 = {11, 12, 13, 14, 15, 16};
    std::vector<int> v2 = {11, 12, 47, 14, 15, 16};
    std::vector<int> v3 = {11, 12, 47, 14, 48, 16};
    std::vector<int> v4 = {47, 12, 13, 14, 15, 16};
    std::vector<int> v5 = {47, 12, 13, 14, 48, 16};
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(2, 47)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(2, 47), v2));
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(2, 47) | extension::boost::adaptors::refined(4, 48)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(2, 47) | extension::boost::adaptors::refined(4, 48), v3));
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(0, 47)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47), v4));
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48)).size(), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48), v5));
    auto r = v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48);
    ASSERT_TRUE(boost::equal(r, v5));
    ASSERT_TRUE(boost::equal(r, v5));
}

TEST(ExtensionBoostAdaptorsRefined, FilledStdList) {
    std::list<int> v1 = {11, 12, 13, 14, 15, 16};
    std::list<int> v2 = {11, 12, 47, 14, 15, 16};
    std::list<int> v3 = {11, 12, 47, 14, 48, 16};
    std::list<int> v4 = {47, 12, 13, 14, 15, 16};
    std::list<int> v5 = {47, 12, 13, 14, 48, 16};
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::refined(2, 47)), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(2, 47), v2));
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::refined(2, 47) | extension::boost::adaptors::refined(4, 48)), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(2, 47) | extension::boost::adaptors::refined(4, 48), v3));
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::refined(0, 47)), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47), v4));
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48)), 6);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48), v5));
    auto r = v1 | extension::boost::adaptors::refined(0, 47) | extension::boost::adaptors::refined(4, 48);
    ASSERT_TRUE(boost::equal(r, v5));
    ASSERT_TRUE(boost::equal(r, v5));
}

TEST(ExtensionBoostAdaptorsRefined, OneElementRawArray) {
    int v1[1] = {11};
    int v4[1] = {47};
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(0, 47)).size(), 1);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47), v4));
}

TEST(ExtensionBoostAdaptorsRefined, OneElementStdArray) {
    std::array<int, 1> v1 = {11};
    std::array<int, 1> v4 = {47};
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(0, 47)).size(), 1);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47), v4));
}

TEST(ExtensionBoostAdaptorsRefined, OneElementStdVector) {
    std::vector<int> v1 = {11};
    std::vector<int> v4 = {47};
    ASSERT_EQ((v1 | extension::boost::adaptors::refined(0, 47)).size(), 1);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47), v4));
}

TEST(ExtensionBoostAdaptorsRefined, OneElementStdList) {
    std::list<int> v1 = {11};
    std::list<int> v4 = {47};
    ASSERT_EQ(boost::size(v1 | extension::boost::adaptors::refined(0, 47)), 1);
    ASSERT_TRUE(boost::equal(v1 | extension::boost::adaptors::refined(0, 47), v4));
}
