#include <monostar_app/main_printing.hpp>

#include <monostar_hamiltonians/hamiltonian_params_agile_affo_dacay.hpp>
#include <monostar_hamiltonians/get_orbital_theta.hpp>

#include <extensions/range_streamer.hpp>
#include <extensions/stream_fromat_stacker.hpp>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>

#include <iostream>
#include <iomanip>

#include <cassert>
#include <cmath>

// #######################################################################
// ## print_foo                                                         ##
// #######################################################################

namespace monostar_app {

void print_input_data(const InterpretedProgramOptions& interpreted_program_options) {
    using namespace extension::boost::stream_pragma;
    // Stream RAII:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    // Print:
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_sites                            = " << interpreted_program_options.n_sites << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_pt                               = " << interpreted_program_options.n_pt << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] model_type                         = " << interpreted_program_options.model_type << std::endl;
    if (interpreted_program_options.model_type == ModelType::AF || interpreted_program_options.model_type == ModelType::FM) {
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::J_classical     = " << interpreted_program_options.hamiltonian_params_af_fm.get_J_classical() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::J_quantum       = " << interpreted_program_options.hamiltonian_params_af_fm.get_J_quantum() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::J_nnn_classical = " << interpreted_program_options.hamiltonian_params_af_fm.get_J_nnn_classical() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::J_nnn_quantum   = " << interpreted_program_options.hamiltonian_params_af_fm.get_J_nnn_quantum() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::B               = " << interpreted_program_options.hamiltonian_params_af_fm.get_B() << std::endl;
    }
    if (interpreted_program_options.model_type == ModelType::JKL01) {
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_jkl01::J               = " << interpreted_program_options.hamiltonian_params_jkl01.get_J() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_jkl01::J_0             = " << interpreted_program_options.hamiltonian_params_jkl01.get_J_0() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_jkl01::J_1             = " << interpreted_program_options.hamiltonian_params_jkl01.get_J_1() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_jkl01::K               = " << interpreted_program_options.hamiltonian_params_jkl01.get_K() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_jkl01::K_0             = " << interpreted_program_options.hamiltonian_params_jkl01.get_K_0() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_jkl01::K_1             = " << interpreted_program_options.hamiltonian_params_jkl01.get_K_1() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_jkl01::L               = " << interpreted_program_options.hamiltonian_params_jkl01.get_L() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_jkl01::L_1             = " << interpreted_program_options.hamiltonian_params_jkl01.get_L_1() << std::endl;
    }
    if (interpreted_program_options.model_type == ModelType::FO) {
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::tau_z_coef         = " << interpreted_program_options.hamiltonian_params_fo.get_tau_z_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::tau_minus_coef     = " << interpreted_program_options.hamiltonian_params_fo.get_tau_minus_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pzz_coef           = " << interpreted_program_options.hamiltonian_params_fo.get_Pzz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pxz_coef           = " << interpreted_program_options.hamiltonian_params_fo.get_Pxz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pxx_coef           = " << interpreted_program_options.hamiltonian_params_fo.get_Pxx_coef() << std::endl;
    }
    if (interpreted_program_options.model_type == ModelType::AgileAFFO) {
        //std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::s_coef          = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_s_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::ss_coef         = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_ss_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::tau_z_coef      = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_tau_z_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::tau_minus_coef  = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_tau_minus_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::Pzz_coef        = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_Pzz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::Pxz_coef        = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_Pxz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::Pxx_coef        = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_Pxx_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::ss_Pzz_coef     = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_ss_Pzz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::ss_Pxz_coef     = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_ss_Pxz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::ss_Pxx_coef     = " << interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian().get_ss_Pxx_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::eps             = " << interpreted_program_options.hamiltonian_params_agile_affo.get_aglie_params().get_eps() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_agile::phi             = " << interpreted_program_options.hamiltonian_params_agile_affo.get_aglie_params().get_phi() << std::endl;
    }
    if (interpreted_program_options.model_type == ModelType::AgileAFFO || interpreted_program_options.model_type == ModelType::FO) {
        if (interpreted_program_options.orbital_theta) {
            std::cout << "[INFO   ] [PROGRAM_OPTIONS] reference orbital theta            = " << *interpreted_program_options.orbital_theta << std::endl;
        } else {
            std::cout << "[INFO   ] [PROGRAM_OPTIONS] reference orbital theta            = "
                      << "<auto: let the program choose the optimal value>" << std::endl;
            if (interpreted_program_options.model_type == ModelType::AgileAFFO) {
                std::cout << "[INFO   ] [PROGRAM_OPTIONS] average_ss (for θ calculations)    = " << interpreted_program_options.average_ss << std::endl;
            }
        }
    }
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] run_type                           = " << interpreted_program_options.run_type << std::endl;
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] es_momentum_domain                 = " << interpreted_program_options.es_momentum_domain << std::endl;
    }
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::unpopulated_basis_flag      = " << interpreted_program_options.print_flags.print_unpopulated_basis_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::unpopulated_basis_size_flag = " << interpreted_program_options.print_flags.print_unpopulated_basis_size_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::populated_basis_flag        = " << interpreted_program_options.print_flags.print_populated_basis_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::populated_basis_size_flag   = " << interpreted_program_options.print_flags.print_populated_basis_size_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::hamiltonian_stats           = " << interpreted_program_options.print_flags.print_hamiltonian_stats << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::sp_hamiltonian_flag         = " << interpreted_program_options.print_flags.print_sp_hamiltonian_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::hamiltonian_flag            = " << interpreted_program_options.print_flags.print_hamiltonian_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::eigen_values_flag           = " << interpreted_program_options.print_flags.print_eigen_values_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::eigen_vectors_flag          = " << interpreted_program_options.print_flags.print_eigen_vectors_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::pretty_vectors_flag         = " << interpreted_program_options.print_flags.print_pretty_vectors_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::density_operator_flag       = " << interpreted_program_options.print_flags.print_density_operator_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::pretty_min_max_n_kstates    = "
              << "[" << interpreted_program_options.print_flags.print_pretty_min_max_n_kstates.first << ":" << interpreted_program_options.print_flags.print_pretty_min_max_n_kstates.second << ")" << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::pretty_probability_treshold = " << interpreted_program_options.print_flags.print_pretty_probability_treshold << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_threads                          = " << interpreted_program_options.n_threads << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] is_optimization_iteration          = " << interpreted_program_options.is_optimization_iteration << std::endl;
}

void print_results_tree(
    const InterpretedProgramOptions& interpreted_program_options,
    const std::shared_ptr<monostar_hamiltonians::HamiltonianReferenceEnergies> reference_energies,
    const std::optional<bfpt_common::CommonRecipeReceipt>& gs_receipt_optional,
    const std::optional<std::vector<bfpt_common::CommonRecipeReceipt>>& es_receipts_optional,
    const unsigned n_sites) {
    // Helpers:
    const auto es_momentum_range_sapn = es_momentum_domain_variant_to_momentum_range_sapn(
        interpreted_program_options.es_momentum_domain,
        interpreted_program_options.n_sites);
    // [[expect]]:
    if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
        assert(gs_receipt_optional);
    }
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        assert(es_receipts_optional);
        assert(es_receipts_optional->size() == es_momentum_range_sapn.second - es_momentum_range_sapn.first);
    }
    // Print:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    const unsigned width = 22;
    const unsigned precision = 14;
    std::cout << std::showpos << std::fixed << std::setprecision(precision) << std::left;
    // Print gs:
    if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
        const auto gs_receipt = gs_receipt_optional.value();
        std::cout << " ├state: gs " << std::endl;
        std::cout << " ││enery             = "
                  << std::setw(width) << gs_receipt.energy << " = "
                  << std::setw(width) << gs_receipt.energy / n_sites << " * n_sites"
                  << std::endl;
        if (reference_energies) {
            if (const auto& gs_energy = reference_energies->get_gs_energy()) {
                std::cout << " ││enery (reference) = "
                          << std::setw(width) << *gs_energy << " = "
                          << std::setw(width) << *gs_energy / n_sites << " * n_sites"
                          << std::endl;
            }
        }
    }
    // Print es:
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        const auto es_receipts = es_receipts_optional.value();
        for (unsigned k_n = es_momentum_range_sapn.first; k_n < es_momentum_range_sapn.second; k_n++) {
            const auto es_result = es_receipts[k_n - es_momentum_range_sapn.first];
            std::cout << " ├state: es [k_n = " << k_n << "]" << std::endl;
            std::cout << " ││abs. enery        = "
                      << std::setw(width) << es_result.energy << " = "
                      << std::setw(width) << es_result.energy / n_sites << " * n_sites"
                      << std::endl;
            if (interpreted_program_options.run_type == RunType::EG) {
                const auto gs_receipt = gs_receipt_optional.value();
                std::cout << " ││exc. enery        = "
                          << std::setw(width) << es_result.energy - gs_receipt.energy << " = "
                          << std::setw(width) << (es_result.energy - gs_receipt.energy) / n_sites << " * n_sites"
                          << std::endl;
            }
            if (reference_energies) {
                if (const auto& get_es_absolute_enery = reference_energies->get_es_absolute_enery(k_n)) {
                    std::cout << " ││abs. enery (ref.) = "
                              << std::setw(width) << *get_es_absolute_enery << " = "
                              << std::setw(width) << (*get_es_absolute_enery) / n_sites << " * n_sites"
                              << std::endl;
                }
                if (const auto& es_exciation_enery = reference_energies->get_es_exciation_enery(k_n)) {
                    std::cout << " ││exc. enery (ref.) = "
                              << std::setw(width) << *es_exciation_enery << " = "
                              << std::setw(width) << (*es_exciation_enery) / n_sites << " * n_sites"
                              << std::endl;
                }
            }
        }
    }
}

void print_post_data(
    const InterpretedProgramOptions& interpreted_program_options,
    /*const std::unique_ptr<monostar_app::ReferenceEnergies> reference_energies,*/
    const std::optional<bfpt_common::CommonRecipeReceipt>& gs_receipt_optional,
    const std::optional<std::vector<bfpt_common::CommonRecipeReceipt>>& es_receipts_optional) {
    // Helpers:
    const auto es_momentum_range_sapn = es_momentum_domain_variant_to_momentum_range_sapn(
        interpreted_program_options.es_momentum_domain,
        interpreted_program_options.n_sites);
    // [[expect]]:
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        assert(es_receipts_optional->size() == es_momentum_range_sapn.second - es_momentum_range_sapn.first);
    }
    // Stream RAII:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    // Helpers:
    using extension::boost::stream_pragma::RSS;
    using extension::boost::stream_pragma::operator|;
    using extension::boost::stream_pragma::operator<<;
    using namespace boost::adaptors;
    // Print gs:
    if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
        const auto gs_receipt = gs_receipt_optional.value();
        std::cout << "[RESULT] [POST] gs_energy: " << gs_receipt.energy << std::endl;
    }
    // Print es -- domain:
    if (interpreted_program_options.run_type == RunType::EG) {
        const auto& n_sites = interpreted_program_options.n_sites;
        const auto nk_to_k = [n_sites](int n_k) -> double { return (2 * arma::datum::pi * n_k) / n_sites; };
        const auto domain = boost::irange(es_momentum_range_sapn.first, es_momentum_range_sapn.second) | transformed(nk_to_k);
        std::cout << "[RESULT] [POST] domain: " << (domain | RSS<double>().like_python_list()) << std::endl;
    }
    // Print es -- helpers:
    const auto receipt_to_energy = [](bfpt_common::CommonRecipeReceipt result) -> double { return result.energy; };
    // Print es -- absolute energies:
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        const auto es_receipt = es_receipts_optional.value();
        std::cout << "[RESULT] [POST] es_absolute_energies: " << (es_receipt | transformed(receipt_to_energy) | RSS<double>().like_python_list()) << std::endl;
    }
    // Print es -- excitation energies:
    if (interpreted_program_options.run_type == RunType::EG) {
        const auto gs_receipt = gs_receipt_optional.value();
        const auto es_receipts = es_receipts_optional.value();
        const auto gs_energy = gs_receipt.energy;
        const auto absolute_energy_into_excitation_energy = [gs_energy](double es_energy) -> double { return es_energy - gs_energy; };
        const auto exciation_energies = es_receipts | transformed(receipt_to_energy) | transformed(absolute_energy_into_excitation_energy);
        std::cout << "[RESULT] [POST] es_excitation_energies: " << (exciation_energies | RSS<double>().like_python_list()) << std::endl;
    }
}

void print_optimization_data(
    const InterpretedProgramOptions& interpreted_program_options,
    const std::optional<bfpt_common::CommonRecipeReceipt>& gs_receipt_optional,
    const std::optional<std::vector<bfpt_common::CommonRecipeReceipt>>& es_receipts_optional) {
    // Expects:
    if (interpreted_program_options.run_type == RunType::G) {
        assert(gs_receipt_optional);
    }
    if (interpreted_program_options.run_type == RunType::E) {
        assert(es_receipts_optional);
        const auto es_receipts = es_receipts_optional.value();
        assert(es_receipts.size() == 1);
    }
    // Stream RAII:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    // Printing:
    if (interpreted_program_options.run_type == RunType::G) {
        const auto gs_receipt = gs_receipt_optional.value();
        std::cout << "[OPTIMIZATION] [OK] " << std::setprecision(15) << gs_receipt.energy << std::endl;
    } else if (interpreted_program_options.run_type == RunType::E && es_momentum_domain_variant_to_enum(interpreted_program_options.es_momentum_domain) == EsMomentumDomain::one) {
        assert(es_receipts_optional);
        const auto es_receipts = es_receipts_optional.value();
        const auto es_result = es_receipts[0];
        std::cout << "[OPTIMIZATION] [OK] " << std::setprecision(15) << es_result.energy << std::endl;
    } else {
        std::cout << "[OPTIMIZATION] [ERROR] "
                  << "Not in any of optimization setups." << std::endl;
    }
}

void print_theta_opt(
    const monostar_hamiltonians::HamiltonianParamsFo& hamiltonian_fo_params,
    const std::optional<double> user_defined_overrule) {
    // Using:
    using namespace extension::boost::stream_pragma;
    using extension::boost::stream_pragma::RSS;
    using extension::boost::stream_pragma::operator|;
    using extension::boost::stream_pragma::operator<<;
    // Print:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    std::cout << std::showpos;
    std::cout << "[INFO   ] [THETA_OPT] H                                        = " << hamiltonian_fo_params.string_repr_in_orbital_operators() << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] H                                        = " << hamiltonian_fo_params.string_repr_in_trigonometric_functions() << std::endl;
    if (const auto& _ = hamiltonian_fo_params.get_theta_opt_analytical()) {
        std::cout << "[INFO   ] [THETA_OPT] optimal orbital θ (analytical)           = " << (_.unwrap() | RSS<double>().like_python_set()) << std::endl;
    } else {
        std::cout << "[INFO   ] [THETA_OPT] optimal orbital θ (analytical)           = "
                  << "<no known analicycal solution solver>" << std::endl;
    }
    std::cout << "[INFO   ] [THETA_OPT] optimal orbital θ (numerical)            = "
              << (hamiltonian_fo_params.get_theta_opt_numerical() | RSS<double>().like_python_set()) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] optimal orbital θ                        = "
              << (hamiltonian_fo_params.get_theta_opt() | RSS<double>().like_python_set()) << std::endl;
    const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(hamiltonian_fo_params, user_defined_overrule);  //may thorw!
    std::cout << "[INFO   ] [THETA_OPT] used orbital θ                           = "
              << orbital_theta_to_use << " = "
              << orbital_theta_to_use / M_PI << " π" << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] cos, sin [used orbital θ]                = "
              << std::cos(orbital_theta_to_use) << ", "
              << std::sin(orbital_theta_to_use) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] cos², sin² [used orbital θ]              = "
              << std::pow(std::cos(orbital_theta_to_use), 2) << ", "
              << std::pow(std::sin(orbital_theta_to_use), 2) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] τᶻ, τˣ, τ⁻, τ⁺ [used orbital θ]          = "
              << +std::cos(orbital_theta_to_use) << ", "
              << -std::cos(orbital_theta_to_use) << ", "
              << -std::sin(orbital_theta_to_use) << ", "
              << +std::sin(orbital_theta_to_use)
              << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] Pᶻ, Pˣ [used orbital θ]                  = "
              << (0.5 + 0.5 * std::cos(orbital_theta_to_use)) << ", "
              << (0.5 - 0.5 * std::cos(orbital_theta_to_use)) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] Pᶻᶻ, Pᶻˣ+Pˣᶻ, Pˣˣ [used orbital θ]       = "
              << std::pow(0.5 + 0.5 * std::cos(orbital_theta_to_use), 2) << ", "
              << 2 * (0.5 + 0.5 * std::cos(orbital_theta_to_use)) * (0.5 - 0.5 * std::cos(orbital_theta_to_use)) << ", "
              << std::pow(0.5 - 0.5 * std::cos(orbital_theta_to_use), 2) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] H [used orbital θ]                       = "
              << hamiltonian_fo_params.get_site_energy(orbital_theta_to_use) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] dH/dθ [used orbital θ]                   = "
              << hamiltonian_fo_params.get_site_energy_derivative(orbital_theta_to_use) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] d²H/dθ² [used orbital θ]                 = "
              << hamiltonian_fo_params.get_site_energy_derivative2(orbital_theta_to_use) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] d³H/dθ³ [used orbital θ]                 = "
              << hamiltonian_fo_params.get_site_energy_derivative3(orbital_theta_to_use) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] d⁴H/dθ⁴ [used orbital θ]                 = "
              << hamiltonian_fo_params.get_site_energy_derivative4(orbital_theta_to_use) << std::endl;
}

void print_theta_opt(
    const monostar_hamiltonians::HamiltonianParamsAffo& hamiltonian_affo_params,
    const std::optional<double> user_defined_overrule,
    const double average_ss) {
    // Using:
    using namespace extension::boost::stream_pragma;
    using extension::boost::stream_pragma::RSS;
    using extension::boost::stream_pragma::operator|;
    using extension::boost::stream_pragma::operator<<;
    // Print:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    std::cout << std::showpos;
    std::cout << "[INFO   ] [THETA_OPT] H                                        = " << hamiltonian_affo_params.string_repr_in_orbital_operators() << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] H                                        = " << hamiltonian_affo_params.string_repr_in_trigonometric_functions() << std::endl;
    if (const auto& _ = hamiltonian_affo_params.get_theta_opt_analytical(average_ss)) {
        std::cout << "[INFO   ] [THETA_OPT] optimal orbital θ (analytical)           = " << (_.unwrap() | RSS<double>().like_python_set()) << std::endl;
    } else {
        std::cout << "[INFO   ] [THETA_OPT] optimal orbital θ (analytical)           = "
                  << "<no known analicycal solution solver>" << std::endl;
    }
    std::cout << "[INFO   ] [THETA_OPT] optimal orbital θ (numerical)            = "
              << (hamiltonian_affo_params.get_theta_opt_numerical(average_ss) | RSS<double>().like_python_set()) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] optimal orbital θ                        = "
              << (hamiltonian_affo_params.get_theta_opt(average_ss) | RSS<double>().like_python_set()) << std::endl;
    const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(hamiltonian_affo_params, user_defined_overrule, average_ss);  //may thorw!
    std::cout << "[INFO   ] [THETA_OPT] used orbital θ                           = "
              << orbital_theta_to_use << " = "
              << orbital_theta_to_use / M_PI << " π" << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] cos, sin [used orbital θ]                = "
              << std::cos(orbital_theta_to_use) << ", "
              << std::sin(orbital_theta_to_use) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] cos², sin² [used orbital θ]              = "
              << std::pow(std::cos(orbital_theta_to_use), 2) << ", "
              << std::pow(std::sin(orbital_theta_to_use), 2) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] τᶻ, τˣ, τ⁻, τ⁺ [used orbital θ]          = "
              << +std::cos(orbital_theta_to_use) << ", "
              << -std::cos(orbital_theta_to_use) << ", "
              << -std::sin(orbital_theta_to_use) << ", "
              << +std::sin(orbital_theta_to_use)
              << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] Pᶻ, Pˣ [used orbital θ]                  = "
              << (0.5 + 0.5 * std::cos(orbital_theta_to_use)) << ", "
              << (0.5 - 0.5 * std::cos(orbital_theta_to_use)) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] Pᶻᶻ, Pᶻˣ+Pˣᶻ, Pˣˣ [used orbital θ]       = "
              << std::pow(0.5 + 0.5 * std::cos(orbital_theta_to_use), 2) << ", "
              << 2 * (0.5 + 0.5 * std::cos(orbital_theta_to_use)) * (0.5 - 0.5 * std::cos(orbital_theta_to_use)) << ", "
              << std::pow(0.5 - 0.5 * std::cos(orbital_theta_to_use), 2) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] H [used orbital θ]                       = "
              << hamiltonian_affo_params.get_site_energy(orbital_theta_to_use, average_ss) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] dH/dθ [used orbital θ]                   = "
              << hamiltonian_affo_params.get_site_energy_derivative(orbital_theta_to_use, average_ss) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] d²H/dθ² [used orbital θ]                 = "
              << hamiltonian_affo_params.get_site_energy_derivative2(orbital_theta_to_use, average_ss) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] d³H/dθ³ [used orbital θ]                 = "
              << hamiltonian_affo_params.get_site_energy_derivative3(orbital_theta_to_use, average_ss) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] d⁴H/dθ⁴ [used orbital θ]                 = "
              << hamiltonian_affo_params.get_site_energy_derivative4(orbital_theta_to_use, average_ss) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] ((integrate out orbitals)H)::B           = "
              << hamiltonian_affo_params.average_out_orbitals_1(orbital_theta_to_use).get_B() << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] ((integrate out orbitals)H)::J_classical = "
              << hamiltonian_affo_params.average_out_orbitals_1(orbital_theta_to_use).get_J_classical() << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] ((integrate out orbitals)H)::J_quantum   = "
              << hamiltonian_affo_params.average_out_orbitals_1(orbital_theta_to_use).get_J_quantum() << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] ((integrate out orbitals)H)::free        = "
              << hamiltonian_affo_params.average_out_orbitals_1(orbital_theta_to_use).get_free() << std::endl;
}

void print_decayed_agile_affo_hamiltonian_params(
    const monostar_hamiltonians::HamiltonianParamsAgileAffo& hamiltonian_params_agile_affo,
    const std::optional<double> user_defined_overrule_for_theta_opt,
    const double average_ss) {
    const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(
        hamiltonian_params_agile_affo.get_so_hamiltonian(),
        user_defined_overrule_for_theta_opt,
        average_ss);
    const monostar_hamiltonians::HamiltonianParamsJkl01 decayed_hamiltonian =
        monostar_hamiltonians::dacay_hamiltonian_params_agile_affo(
            hamiltonian_params_agile_affo,
            orbital_theta_to_use);
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::J       = " << decayed_hamiltonian.get_J() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::J_0     = " << decayed_hamiltonian.get_J_0() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::J_1     = " << decayed_hamiltonian.get_J_1() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::K       = " << decayed_hamiltonian.get_K() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::K_0     = " << decayed_hamiltonian.get_K_0() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::K_1     = " << decayed_hamiltonian.get_K_1() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::L       = " << decayed_hamiltonian.get_L() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::L_1     = " << decayed_hamiltonian.get_L_1() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::(K+L)   = " << decayed_hamiltonian.get_K() + decayed_hamiltonian.get_L() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::(J+J_0) = " << decayed_hamiltonian.get_J() + decayed_hamiltonian.get_J_0() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::(J+J_1) = " << decayed_hamiltonian.get_J() + decayed_hamiltonian.get_J_1() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::(K+K_0) = " << decayed_hamiltonian.get_K() + decayed_hamiltonian.get_K_0() << std::endl;
    std::cout << "[INFO   ] [AGILE-TO-JKL01] ((decay)hamiltonian_agile)::(K+K_1) = " << decayed_hamiltonian.get_K() + decayed_hamiltonian.get_K_1() << std::endl;
}

}  // end of namespace monostar_app
