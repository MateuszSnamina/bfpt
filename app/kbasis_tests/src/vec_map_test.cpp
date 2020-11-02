#include <kbasis/vec_map.hpp>

#include <memory>
#include <string>

#include <gtest/gtest.h>

namespace {
struct Person {
    Person(std::string county, std::string given_name, std::string family_name) : county(county),
                                                                                  given_name(given_name),
                                                                                  family_name(family_name) {
    }
    const std::string county;
    const std::string given_name;
    const std::string family_name;
    std::string to_str() const {
        return "[" + county + "][" + given_name + "][" + family_name + "]";
    }
};

struct PersonKeyExtractor {
    typedef std::string result_type;
    std::string operator()(const std::shared_ptr<Person>& p) const {
        return p->family_name + p->given_name;
    }
};

struct PersonComparisonPredicate {
    bool operator()(const std::string& l, const std::string& r) const {
        return l > r;  // reverse order
    }
};
}  // namespace

TEST(VecMap, Empty) {
    kbasis::VecMap<Person, PersonKeyExtractor, PersonComparisonPredicate> vec_map;
    // test vec_map size:
    ASSERT_EQ(vec_map.size(), 0);
    ASSERT_TRUE(vec_map.vec_index().begin() == vec_map.vec_index().end());
    ASSERT_TRUE(vec_map.map_index().begin() == vec_map.map_index().end());
    // test not-finding a not existing element:
    const auto k100 = std::make_shared<Person>("usa", "aaa", "bbb");
    EXPECT_FALSE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k100)));
}

TEST(VecMap, OneElement) {
    kbasis::VecMap<Person, PersonKeyExtractor, PersonComparisonPredicate> vec_map;
    const auto k1 = std::make_shared<Person>("usa", "aaa", "bbb");
    vec_map.add_element(k1);
    // test vec_map size:
    ASSERT_EQ(vec_map.size(), 1);
    ASSERT_TRUE(vec_map.vec_index().begin() + 1 == vec_map.vec_index().end());
    ASSERT_TRUE(std::next(vec_map.map_index().begin(), 1) == vec_map.map_index().end());
    // test element access:
    EXPECT_TRUE(*vec_map.vec_index().begin() == k1);
    EXPECT_TRUE(*vec_map.map_index().begin() == k1);
    // test finding the existing element:
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)));
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)), 0);
    const auto k10 = std::make_shared<Person>("usa", "aaa", "ccc");
    EXPECT_FALSE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k10)));
    // test not-finding a not existing element:
    const auto k100 = std::make_shared<Person>("usa", "bbb", "bbb");
    EXPECT_FALSE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k100)));
}

TEST(VecMap, TwoDifferentElementsTest0) {
    const auto k1 = std::make_shared<Person>("usa", "aaa", "ccc");
    const auto k2 = std::make_shared<Person>("usa", "aaa", "bbb");
    kbasis::VecMap<Person, PersonKeyExtractor, PersonComparisonPredicate> vec_map;
    vec_map.add_element(k1);
    vec_map.add_element(k2);
    // test vec_map size:
    ASSERT_EQ(vec_map.size(), 2);
    ASSERT_TRUE(vec_map.vec_index().begin() + 2 == vec_map.vec_index().end());
    ASSERT_TRUE(std::next(vec_map.map_index().begin(), 2) == vec_map.map_index().end());
    // test element access:
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 0) == k1);
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 1) == k2);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 0) == k1);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 1) == k2);
    // test finding the existing elements:
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k2)));
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)), 0);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k2)), 1);
    // test not-finding a not existing element:
    const auto k100 = std::make_shared<Person>("usa", "bbb", "ccc");
    EXPECT_FALSE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k100)));
}

TEST(VecMap, TwoDifferentElementsTest1) {
    const auto k1 = std::make_shared<Person>("usa", "aaa", "ccc");
    const auto k2 = std::make_shared<Person>("usa", "aaa", "bbb");
    kbasis::VecMap<Person, PersonKeyExtractor, PersonComparisonPredicate> vec_map;
    vec_map.add_element(k2);
    vec_map.add_element(k1);
    // test vec_map size:
    ASSERT_EQ(vec_map.size(), 2);
    ASSERT_TRUE(vec_map.vec_index().begin() + 2 == vec_map.vec_index().end());
    ASSERT_TRUE(std::next(vec_map.map_index().begin(), 2) == vec_map.map_index().end());
    // test element access:
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 0) == k2);
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 1) == k1);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 0) == k1);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 1) == k2);
    // test finding the existing elements:
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k2)));
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k2)), 0);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)), 1);
    // test not-finding a not existing element:
    const auto k100 = std::make_shared<Person>("usa", "bbb", "ccc");
    EXPECT_FALSE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k100)));
}

TEST(VecMap, TwoSameElementsTest0) {
    const auto k1 = std::make_shared<Person>("usa", "aaa", "ccc");
    const auto k2 = std::make_shared<Person>("usa", "aaa", "ccc");
    kbasis::VecMap<Person, PersonKeyExtractor, PersonComparisonPredicate> vec_map;
    vec_map.add_element(k1);
    vec_map.add_element(k2);
    // test vec_map size:
    ASSERT_EQ(vec_map.size(), 1);
    ASSERT_TRUE(vec_map.vec_index().begin() + 1 == vec_map.vec_index().end());
    ASSERT_TRUE(std::next(vec_map.map_index().begin(), 1) == vec_map.map_index().end());
    // test element access:
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 0) == k1);
    EXPECT_FALSE(*(vec_map.vec_index().begin() + 0) == k2);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 0) == k1);
    EXPECT_FALSE(*std::next(vec_map.map_index().begin(), 0) == k2);
    // test finding the existing elements:
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k2)));
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)), 0);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k2)), 0);
    // test not-finding a not existing element:
    const auto k100 = std::make_shared<Person>("usa", "bbb", "ccc");
    EXPECT_FALSE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k100)));
}

TEST(VecMap, BigTest) {
    const auto k1 = std::make_shared<Person>("usa", "aaa", "ccc");  // key: "cccaaa"
    const auto k2 = std::make_shared<Person>("usa", "bbb", "ccc");  // key: "cccbbb"
    const auto k3 = std::make_shared<Person>("usa", "aaa", "ddd");  // key: "dddaaa"
    const auto k4 = std::make_shared<Person>("usa", "bbb", "ccc");  // key: "cccbbb"// is same as for k2
    const auto k5 = std::make_shared<Person>("pl", "aaa", "nnn");   // key: "nnnaaa"
    const auto k6 = std::make_shared<Person>("usa", "aa", "cccn");  // key: "ccccnaa"
    const auto k7 = std::make_shared<Person>("usa", "bb", "cccb");  // key: "cccbbb" // is same as for k2
    const auto k8 = std::make_shared<Person>("pl", "aaa", "ddd");   // key: "dddaaa" // is same as for k3
    const auto k9 = std::make_shared<Person>("pl", "aaa", "mmm");   // key: "mmmaaa"

    kbasis::VecMap<Person, PersonKeyExtractor, PersonComparisonPredicate> vec_map;
    vec_map.add_element(k1);
    vec_map.add_element(k2);
    vec_map.add_element(k3);
    vec_map.add_element(k4);
    vec_map.add_element(k5);
    vec_map.add_element(k6);
    vec_map.add_element(k7);
    vec_map.add_element(k8);
    vec_map.add_element(k9);
    // test vec_map size:
    ASSERT_EQ(vec_map.size(), 6);
    ASSERT_TRUE(vec_map.vec_index().begin() + 6 == vec_map.vec_index().end());
    ASSERT_TRUE(std::next(vec_map.map_index().begin(), 6) == vec_map.map_index().end());
    // test element access:
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 0) == k1);
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 1) == k2);
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 2) == k3);
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 3) == k5);
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 4) == k6);
    EXPECT_TRUE(*(vec_map.vec_index().begin() + 5) == k9);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 0) == k5);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 1) == k9);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 2) == k3);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 3) == k6);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 4) == k2);
    EXPECT_TRUE(*std::next(vec_map.map_index().begin(), 5) == k1);
    // test finding the existing elements:
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k2)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k3)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k4)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k5)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k6)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k7)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k8)));
    ASSERT_TRUE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k9)));
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k1)), 0);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k2)), 1);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k3)), 2);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k4)), 1);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k5)), 3);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k6)), 4);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k7)), 1);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k8)), 2);
    EXPECT_EQ(*vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k9)), 5);
    // test not-finding a not existing element:
    const auto k100 = std::make_shared<Person>("usa", "www", "ccc");
    EXPECT_FALSE(vec_map.find_element_and_get_its_ra_index(PersonKeyExtractor()(k100)));
}
