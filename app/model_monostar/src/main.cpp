#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/interpreted_program_options.hpp>
#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_hamiltonian_kernel.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>
#include <model_monostar/reference_energies.hpp>

#include <bfpt_common/hamiltonian_kernel.hpp>
#include <bfpt_common/generic_kstate_hamiltonian.hpp>
#include <bfpt_common/do_common_recipie.hpp>

#include <extensions/range_streamer.hpp>
#include <extensions/stream_fromat_stacker.hpp>

#include <armadillo>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>

#include <iostream>

#include <cassert>

// #######################################################################
// ## main...                                                           ##
// #######################################################################

double bfpt_gs(
        const bfpt_common::HamiltonianKernel1<model_monostar::MonostarSiteState>& hamiltonian_kernel_1,
        const bfpt_common::HamiltonianKernel12<model_monostar::MonostarSiteState>& hamiltonian_kernel_12,
        const size_t n_sites, const unsigned max_pt_order,
        const bfpt_common::CommonRecipePrintFlags& print_flags,
        unsigned n_threads) {
    //using SiteStateT = model_monostar::MonostarSiteState;
    using KstateT = model_monostar::DynamicMonostarKstate;
    using BasisT = model_monostar::DynamicMonostarKstateBasis;
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(model_monostar::classical_gs_kstate(n_sites)));
    const bfpt_common::GenericKstateHamiltonian<KstateT> hamiltonian{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12};
    return bfpt_common::do_common_recipe(hamiltonian, hamiltonian, basis,
                                         max_pt_order, 0,
                                         print_flags, "[gs] ",
                                         n_threads);
}

double bfpt_kn_es(
        const bfpt_common::HamiltonianKernel1<model_monostar::MonostarSiteState>& hamiltonian_kernel_1,
        const bfpt_common::HamiltonianKernel12<model_monostar::MonostarSiteState>& hamiltonian_kernel_12,
        const size_t n_sites, const unsigned max_pt_order, const unsigned k_n,
        const bfpt_common::CommonRecipePrintFlags& print_flags,
        unsigned n_threads) {
    //using SiteStateT = model_monostar::MonostarSiteState;
    using KstateT = model_monostar::DynamicMonostarKstate;
    using BasisT = model_monostar::DynamicMonostarKstateBasis;
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(model_monostar::classical_es_kstate(n_sites)));
    const bfpt_common::GenericKstateHamiltonian<KstateT> hamiltonian{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12};
    return bfpt_common::do_common_recipe(hamiltonian, hamiltonian, basis,
                                         max_pt_order, k_n,
                                         print_flags, "[es (" + std::to_string(k_n) + ")] ",
                                         n_threads);
}

void print_input_data(const InterpretedProgramOptions& interpreted_program_options) {
    using namespace extension::boost::stream_pragma;
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_sites                            = " << interpreted_program_options.n_sites << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_pt                               = " << interpreted_program_options.n_pt << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] model_type                         = " << interpreted_program_options.model_type << std::endl;
    if (interpreted_program_options.model_type == ModelType::AF || interpreted_program_options.model_type == ModelType::FM){
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::J_classical     = " << interpreted_program_options.hamiltonian_af_fm_params.get_J_classical() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::J_quantum       = " << interpreted_program_options.hamiltonian_af_fm_params.get_J_quantum() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::B               = " << interpreted_program_options.hamiltonian_af_fm_params.get_B() << std::endl;
    }
    if (interpreted_program_options.model_type == ModelType::FO) {
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pdelta_coef        = " << interpreted_program_options.hamiltonian_fo_params.get_Pdelta_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pxx_coef           = " << interpreted_program_options.hamiltonian_fo_params.get_Pxx_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pxz_coef           = " << interpreted_program_options.hamiltonian_fo_params.get_Pxz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pzz_coef           = " << interpreted_program_options.hamiltonian_fo_params.get_Pzz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] reference orbital theta            = " << interpreted_program_options.theta_opt << std::endl;
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
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::pretty_min_max_n_kstates    = " << "[" << interpreted_program_options.print_flags.print_pretty_min_max_n_kstates.first << ":" << interpreted_program_options.print_flags.print_pretty_min_max_n_kstates.second << ")" << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::pretty_probability_treshold = " << interpreted_program_options.print_flags.print_pretty_probability_treshold << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_threads                         = " << interpreted_program_options.n_threads << std::endl;
}


void print_results_tree(
        const InterpretedProgramOptions& interpreted_program_options,
        const std::shared_ptr<model_monostar::ReferenceEnergies> reference_energies,
        const std::optional<double>& gs_energy,
        const std::optional<std::vector<double>>& es_energies) {
    const auto es_momentum_range_sapn = es_momentum_domain_variant_to_momentum_range_sapn(
                interpreted_program_options.es_momentum_domain,
                interpreted_program_options.n_sites);
    assert(es_energies->size() == es_momentum_range_sapn.second - es_momentum_range_sapn.first);
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
        std::cout << " ├state: gs "  << std::endl;
        std::cout << " ││enery = " << *gs_energy << std::endl;
        if (reference_energies) {
            std::cout << " ││enery (reference) = " << reference_energies->get_gs_energy() << std::endl;
        }
    }
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        for (unsigned k_n =  es_momentum_range_sapn.first; k_n < es_momentum_range_sapn.second; k_n++) {
            const auto es_energy = (*es_energies)[k_n];
            std::cout << " ├state: es [k_n = " << k_n << "]" << std::endl;
            std::cout << " ││abs. enery        = " << es_energy << std::endl;
            if (gs_energy) {
                std::cout << " ││exc. enery        = " << es_energy - *gs_energy << std::endl;
                if (reference_energies) {
                    std::cout << " ││exc. enery (ref.) = " << reference_energies->get_es_exciation_enery(k_n) << std::endl;
                }
            }
        }
    }
}

void print_post_data(
        const InterpretedProgramOptions& interpreted_program_options,
        /*const std::unique_ptr<model_monostar::ReferenceEnergies> reference_energies,*/
        const std::optional<double>& gs_energy,
        const std::optional<std::vector<double>>& es_energies) {
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    using  extension::boost::stream_pragma::RSS;
    using extension::boost::stream_pragma::operator|;
    using extension::boost::stream_pragma::operator<<;
    using namespace boost::adaptors;
    const auto es_momentum_range_sapn = es_momentum_domain_variant_to_momentum_range_sapn(
                interpreted_program_options.es_momentum_domain,
                interpreted_program_options.n_sites);
    assert(es_energies->size() == es_momentum_range_sapn.second - es_momentum_range_sapn.first);
    if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
        std::cout << "[RESULT] [POST] gs_energy: " << *gs_energy << std::endl;
    }
    if (interpreted_program_options.run_type == RunType::EG) {
        const auto& n_sites = interpreted_program_options.n_sites;
        const auto nk_to_k =
                [n_sites](int n_k)->double{return (2 * arma::datum::pi * n_k) / n_sites;};
        const auto domain = boost::irange(es_momentum_range_sapn.first, es_momentum_range_sapn.second) | transformed(nk_to_k);
        std::cout << "[RESULT] [POST] domain: " << (domain | RSS<double>().like_python_list()) << std::endl;
    }
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        std::cout << "[RESULT] [POST] es_absolute_energies: " << ((*es_energies) | RSS<double>().like_python_list()) << std::endl;
    }
    if (interpreted_program_options.run_type == RunType::EG) {
        const auto absolute_energy_into_excitation_energy =
                [gs_energy](double es_energy)->double{return es_energy - *gs_energy;};
        const auto exciation_energies = (*es_energies) | transformed(absolute_energy_into_excitation_energy);
        std::cout << "[RESULT] [POST] es_excitation_energies: " << (exciation_energies | RSS<double>().like_python_list()) << std::endl;
    }
}

int main(int argc, char** argv) {
    try {
        // ******************************************************************
        const RawProgramOptions raw_program_options = grep_program_options(argc, argv);
        const InterpretedProgramOptions interpreted_program_options = interpret_program_options(raw_program_options);
        // ******************************************************************
        print_input_data(interpreted_program_options);
        // ******************************************************************
        const auto hamiltonian_kernel_1 = [&interpreted_program_options]() {
            const double B = interpreted_program_options.hamiltonian_af_fm_params.get_B();
            return model_monostar::prepare_hamiltonian_kernel_1_af_fm(B);
        }();
        const auto hamiltonian_kernel_12 = [&interpreted_program_options]() {
            const auto J_classical = interpreted_program_options.hamiltonian_af_fm_params.get_J_classical();
            const auto J_quantum = interpreted_program_options.hamiltonian_af_fm_params.get_J_quantum();
            switch (interpreted_program_options.model_type) {
            case ModelType::AF:
                return model_monostar::prepare_hamiltonian_kernel_12_af(J_classical, J_quantum);
            case ModelType::FM:
                return model_monostar::prepare_hamiltonian_kernel_12_fm(J_classical, J_quantum);
            default:
                assert(false);
                return model_monostar::prepare_hamiltonian_kernel_12_af(J_classical, J_quantum);
            };
            //            if (interpreted_program_options.model_type == ModelType::AF) {
            //                return model_monostar::prepare_hamiltonian_kernel_12_af(J_classical, J_quantum);
            //            }
            //            if (interpreted_program_options.model_type == ModelType::FM) {
            //                return model_monostar::prepare_hamiltonian_kernel_12_fm(J_classical, J_quantum);
            //            }
            //            assert(false);
        }();
        // ******************************************************************
        const std::shared_ptr<model_monostar::ReferenceEnergies> reference_energies =
                [&interpreted_program_options]() {
            const auto J_classical = interpreted_program_options.hamiltonian_af_fm_params.get_J_classical();
            const auto J_quantum = interpreted_program_options.hamiltonian_af_fm_params.get_J_quantum();
            const auto B = interpreted_program_options.hamiltonian_af_fm_params.get_B();
            switch (interpreted_program_options.model_type) {
            case ModelType::AF:
                if (interpreted_program_options.model_type == ModelType::AF && J_classical == J_quantum && B == 0) {
                    return std::dynamic_pointer_cast<model_monostar::ReferenceEnergies>(
                                std::make_shared<model_monostar::ReferenceEnergiesAf>(
                                    interpreted_program_options.n_sites, J_classical));
                } else {
                    return std::shared_ptr<model_monostar::ReferenceEnergies>(nullptr);
                }
                return std::shared_ptr<model_monostar::ReferenceEnergies>(nullptr);
            case ModelType::FM:
                return std::dynamic_pointer_cast<model_monostar::ReferenceEnergies>(
                            std::make_shared<model_monostar::ReferenceEnergiesFm>(
                                interpreted_program_options.n_sites, J_classical, J_quantum, B));
            default:
                assert(false);
                return std::shared_ptr<model_monostar::ReferenceEnergies>(nullptr);
            }
        }();
        // ******************************************************************
        const std::optional<double> gs_energy =
                [&interpreted_program_options, &hamiltonian_kernel_1, &hamiltonian_kernel_12]() -> std::optional<double>{
            if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
                std::cout << "------------------------------------------" << std::endl;
                return bfpt_gs(
                            hamiltonian_kernel_1,
                            hamiltonian_kernel_12,
                            interpreted_program_options.n_sites, interpreted_program_options.n_pt,
                            interpreted_program_options.print_flags,
                            interpreted_program_options.n_threads);
                std::cout << "------------------------------------------" << std::endl;
            }
            return std::nullopt;
        }();
        // ******************************************************************
        const std::optional<std::vector<double>> es_energies =
                [&interpreted_program_options, &hamiltonian_kernel_1, &hamiltonian_kernel_12]() -> std::optional<std::vector<double>> {
            std::cout << "------------------------------------------" << std::endl;
            if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
                const auto es_momentum_range_sapn = es_momentum_domain_variant_to_momentum_range_sapn(
                            interpreted_program_options.es_momentum_domain,
                            interpreted_program_options.n_sites);
                std::vector<double> es_energies;
                for (unsigned k_n = es_momentum_range_sapn.first; k_n < es_momentum_range_sapn.second; k_n++) {
                    std::cout << "[PROGRESS] " << "solving n_k: " << k_n << std::endl;
                    const double es_energy = bfpt_kn_es(
                                hamiltonian_kernel_1,
                                hamiltonian_kernel_12,
                                interpreted_program_options.n_sites, interpreted_program_options.n_pt, k_n,
                                interpreted_program_options.print_flags,
                                interpreted_program_options.n_threads);
                    es_energies.push_back(es_energy);
                    std::cout << "------------------------------------------" << std::endl;
                }
                return es_energies;
            }
            return std::nullopt;
        }();
        // ******************************************************************
        print_results_tree(
                    interpreted_program_options,
                    reference_energies,
                    gs_energy,
                    es_energies);
        print_post_data(
                    interpreted_program_options,
                    /*reference_energies,*/
                    gs_energy,
                    es_energies);
        // ******************************************************************
    } catch (std::exception& e) {
        std::cerr << "[ERROR  ] Abnormal termination!" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
}

