#include <kstate/basis_streamer.hpp>

#include <kstate/kstate_concrete.hpp>

#include <gtest/gtest.h>

using kstate::ctr_from_range;

kstate::Basis<kstate::DynamicKstate<int>> get_testet_filled_basis() {
    const int v0[3] = {7, 12, 13};
    const int v1[3] = {11, 12, 13};
    const int v2[3] = {13, 14, 15};
    const int v3[3] = {13, 14, 15};  // replica of v2
    const int v4[3] = {13, 14, 15};  // replica of v2
    const int v5[3] = {1, 20, 15};
    const int v6[3] = {1, 14, 15};
    const int v7[3] = {1, 20, 15};  // replica of v5
    const int v8[3] = {3, 20, 15};
    const int v9[3] = {3, 20, 18};
    const int v10[3] = {3, 21, 15};
    const int v11[3] = {3, 21, 10};
    const auto k0 = std::make_shared<kstate::DynamicKstate<int>>(v0, ctr_from_range);
    const auto k1 = std::make_shared<kstate::DynamicKstate<int>>(v1, ctr_from_range);
    const auto k2 = std::make_shared<kstate::DynamicKstate<int>>(v2, ctr_from_range);
    const auto k3 = std::make_shared<kstate::DynamicKstate<int>>(v3, ctr_from_range);
    const auto k4 = std::make_shared<kstate::DynamicKstate<int>>(v4, ctr_from_range);
    const auto k5 = std::make_shared<kstate::DynamicKstate<int>>(v5, ctr_from_range);
    const auto k6 = std::make_shared<kstate::DynamicKstate<int>>(v6, ctr_from_range);
    const auto k7 = std::make_shared<kstate::DynamicKstate<int>>(v7, ctr_from_range);
    const auto k8 = std::make_shared<kstate::DynamicKstate<int>>(v8, ctr_from_range);
    const auto k9 = std::make_shared<kstate::DynamicKstate<int>>(v9, ctr_from_range);
    const auto k10 = std::make_shared<kstate::DynamicKstate<int>>(v10, ctr_from_range);
    const auto k11 = std::make_shared<kstate::DynamicKstate<int>>(v11, ctr_from_range);
    kstate::Basis<kstate::DynamicKstate<int>> basis(3);
    basis.add_element(k0);
    basis.add_element(k1);
    basis.add_element(k2);
    basis.add_element(k3);
    basis.add_element(k4);
    basis.add_element(k5);
    basis.add_element(k6);
    basis.add_element(k7);
    basis.add_element(k8);
    basis.add_element(k9);
    basis.add_element(k10);
    basis.add_element(k11);
    return basis;
}

TEST(BasisStreamer, BasicTest) {
    using extension::boost::stream_pragma::RSS;
    using kstate::pramga::operator&&;

    const auto& basis = get_testet_filled_basis();
    const auto kstate_value_putter = [](std::ostream& os, kstate::DynamicKstate<int> kstate) {
        using extension::boost::stream_pragma::RSS;
        using kstate::pramga::operator||;
        using kstate::pramga::operator<<;
        os << ( kstate || RSS<int>());
    };
    const auto range_streamer_settings_for_basis = RSS<kstate::DynamicKstate<int>>().set_stream_value_putter(kstate_value_putter);
    const auto basis_streamer = basis && range_streamer_settings_for_basis;
    std::string expected_string =
            "𝔹𝔸𝕊𝕀𝕊-BEGIN\n"
            " -        0 : 𝕂𝕤𝕥𝕒𝕥𝕖⦃7∙12∙13⦄\n"
            " -        1 : 𝕂𝕤𝕥𝕒𝕥𝕖⦃11∙12∙13⦄\n"
            " -        2 : 𝕂𝕤𝕥𝕒𝕥𝕖⦃13∙14∙15⦄\n"
            " -        3 : 𝕂𝕤𝕥𝕒𝕥𝕖⦃1∙20∙15⦄\n"
            " -        4 : 𝕂𝕤𝕥𝕒𝕥𝕖⦃1∙14∙15⦄\n"
            " -        5 : 𝕂𝕤𝕥𝕒𝕥𝕖⦃3∙20∙15⦄\n"
            " -        6 : 𝕂𝕤𝕥𝕒𝕥𝕖⦃3∙20∙18⦄\n"
            " -        7 : 𝕂𝕤𝕥𝕒𝕥𝕖⦃3∙21∙15⦄\n"
            " -        8 : 𝕂𝕤𝕥𝕒𝕥𝕖⦃3∙21∙10⦄\n"
            "𝔹𝔸𝕊𝕀𝕊-END\n";
    EXPECT_EQ(basis_streamer.str(), expected_string);
}
