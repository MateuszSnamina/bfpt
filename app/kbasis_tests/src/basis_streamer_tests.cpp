#include <kbasis_tests/site_state_trait_for_int.hpp>

#include <kbasis/basis_streamer.hpp>

#include <kstate_impl/kstate_concrete_stl.hpp>

#include <gtest/gtest.h>

using kstate_impl::ctr_from_range;

kbasis::Basis<kstate_trait::TraitKstate<kstate_impl::DynamicStlKstate<kstate_trait::TraitSiteState<int>>>>
get_testet_filled_basis() {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicStlKstate<SiteStateTrait>;
    using KStateTrait = kstate_trait::TraitKstate<Kstate>;
    kbasis::Basis<KStateTrait> basis(3);
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
    const auto k0 = std::make_shared<Kstate>(v0, ctr_from_range);
    const auto k1 = std::make_shared<Kstate>(v1, ctr_from_range);
    const auto k2 = std::make_shared<Kstate>(v2, ctr_from_range);
    const auto k3 = std::make_shared<Kstate>(v3, ctr_from_range);
    const auto k4 = std::make_shared<Kstate>(v4, ctr_from_range);
    const auto k5 = std::make_shared<Kstate>(v5, ctr_from_range);
    const auto k6 = std::make_shared<Kstate>(v6, ctr_from_range);
    const auto k7 = std::make_shared<Kstate>(v7, ctr_from_range);
    const auto k8 = std::make_shared<Kstate>(v8, ctr_from_range);
    const auto k9 = std::make_shared<Kstate>(v9, ctr_from_range);
    const auto k10 = std::make_shared<Kstate>(v10, ctr_from_range);
    const auto k11 = std::make_shared<Kstate>(v11, ctr_from_range);
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
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicStlKstate<SiteStateTrait>;

    using extension::boost::stream_pragma::RSS;
    using kbasis::pramga::operator&&;

    const auto& basis = get_testet_filled_basis();
    const auto kstate_value_putter = [](std::ostream& os, Kstate kstate) {
        using extension::boost::stream_pragma::RSS;
        using kstate_impl::pramga::operator||;
        using kstate_impl::pramga::operator<<;
        os << (kstate || RSS<int>());
    };
    const auto range_streamer_settings_for_basis = RSS<Kstate>().set_stream_value_putter(kstate_value_putter);
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
