#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/interpreted_program_options.hpp>
#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_hamiltonian.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>
#include <model_monostar/reference_energies.hpp>

#include <bfpt_common/hamiltonian_12.hpp>
#include <bfpt_common/generic_kstate_hamiltonian.hpp>
#include <bfpt_common/do_common_recipie.hpp>

#include <armadillo>

#include <iostream>

#include <cassert>

// #######################################################################
// ## main...                                                           ##
// #######################################################################

double bfpt_gs(
        const bfpt_common::Hamiltonian12<model_monostar::MonostarSiteState>& hamiltonian_12,
        const size_t n_sites, const unsigned max_pt_order,
        const bfpt_common::CommonRecipePrintFlags& print_flags) {
    //using SiteStateT = model_monostar::MonostarSiteState;
    using KstateT = model_monostar::DynamicMonostarUniqueKstate;
    using BasisT = model_monostar::DynamicMonostarUniqueKstateBasis;
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(model_monostar::classical_gs_kstate(n_sites)));
    const bfpt_common::GenericKstateHamiltonian<KstateT> hamiltonian{n_sites, hamiltonian_12};
    return bfpt_common::do_common_recipe(hamiltonian, hamiltonian, basis,
                                         max_pt_order, 0, print_flags);
}

double bfpt_kn_es(
        const bfpt_common::Hamiltonian12<model_monostar::MonostarSiteState>& hamiltonian_12,
        const size_t n_sites, const unsigned max_pt_order, const unsigned k_n,
        const bfpt_common::CommonRecipePrintFlags& print_flags) {
    //using SiteStateT = model_monostar::MonostarSiteState;
    using KstateT = model_monostar::DynamicMonostarUniqueKstate;
    using BasisT = model_monostar::DynamicMonostarUniqueKstateBasis;
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(model_monostar::classical_es_kstate(n_sites)));
    const bfpt_common::GenericKstateHamiltonian<KstateT> hamiltonian{n_sites, hamiltonian_12};
    return bfpt_common::do_common_recipe(hamiltonian, hamiltonian, basis,
                                         max_pt_order, k_n, print_flags);
}

void print_input_data(const InterpretedProgramOptions& interpreted_program_options) {
    using namespace extension::boost::stream_pragma;
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_sites                           = " << interpreted_program_options.n_sites << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_pt                              = " << interpreted_program_options.n_pt << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] model_type                        = " << interpreted_program_options.model_type << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] J_classical                       = " << interpreted_program_options.J_classical << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] J_quantum                         = " << interpreted_program_options.J_quantum << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] run_type                          = " << interpreted_program_options.run_type << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print_unpopulated_basis_flag      = " << interpreted_program_options.print_flags.print_unpopulated_basis_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print_unpopulated_basis_size_flag = " << interpreted_program_options.print_flags.print_unpopulated_basis_size_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print_populated_basis_flag        = " << interpreted_program_options.print_flags.print_populated_basis_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print_populated_basis_size_flag   = " << interpreted_program_options.print_flags.print_populated_basis_size_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print_sp_hamiltonian_flag         = " << interpreted_program_options.print_flags.print_sp_hamiltonian_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print_hamiltonian_flag            = " << interpreted_program_options.print_flags.print_hamiltonian_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print_eigen_values_flag           = " << interpreted_program_options.print_flags.print_eigen_values_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print_eigen_vectors_flag          = " << interpreted_program_options.print_flags.print_eigen_vectors_flag << std::endl;
}


int main(int argc, char** argv) {
    try {
        // ******************************************************************
        const RawProgramOptions raw_program_options = grep_program_options(argc, argv);
        const InterpretedProgramOptions interpreted_program_options = interpret_program_options(raw_program_options);
        // ******************************************************************
        print_input_data(interpreted_program_options);
        // ******************************************************************
        const auto hamiltonian_12 = [&interpreted_program_options](){
            const auto J_classical = interpreted_program_options.J_classical;
            const auto J_quantum = interpreted_program_options.J_quantum      ;
            if (interpreted_program_options.model_type == ModelType::AF) {
                return model_monostar::prepare_hamiltonian_12_af(J_classical, J_quantum);
            }
            if (interpreted_program_options.model_type == ModelType::FM) {
                return model_monostar::prepare_hamiltonian_12_fm(J_classical, J_quantum);
            }
            assert(false);
        }();
        // ******************************************************************
        const std::unique_ptr<model_monostar::ReferenceEnergies> reference_energies =
                [&interpreted_program_options]() -> std::unique_ptr<model_monostar::ReferenceEnergies> {
                if (interpreted_program_options.model_type == ModelType::AF &&
                    interpreted_program_options.J_classical == interpreted_program_options.J_quantum) {
                return std::make_unique<model_monostar::ReferenceEnergiesAf>(
                    interpreted_program_options.n_sites,
                    interpreted_program_options.J_classical);
    }
                if (interpreted_program_options.model_type == ModelType::FM) {
                return std::make_unique<model_monostar::ReferenceEnergiesFm>(
                    interpreted_program_options.n_sites,
                    interpreted_program_options.J_classical,
                    interpreted_program_options.J_quantum);
    }
                return nullptr;
    }();
        // ******************************************************************
        const std::optional<double> gs_energy =
                [&interpreted_program_options, &hamiltonian_12]() -> std::optional<double>{
            if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
                std::cout << "------------------------------------------" << std::endl;
                return bfpt_gs(
                            hamiltonian_12,
                            interpreted_program_options.n_sites, interpreted_program_options.n_pt,
                            interpreted_program_options.print_flags);
                std::cout << "------------------------------------------" << std::endl;
            }
            return std::nullopt;
        }();
        // ******************************************************************
        const std::optional<std::vector<double>> es_energies =
                [&interpreted_program_options, &hamiltonian_12]() -> std::optional<std::vector<double>> {
            std::cout << "------------------------------------------" << std::endl;
            if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
                std::vector<double> es_energies;
                for (unsigned k_n = 0; k_n < interpreted_program_options.n_sites; k_n++) {
                    const double es_energy = bfpt_kn_es(
                                hamiltonian_12,
                                interpreted_program_options.n_sites, interpreted_program_options.n_pt, k_n,
                                interpreted_program_options.print_flags);
                    es_energies.push_back(es_energy);
                    std::cout << "------------------------------------------" << std::endl;
                }
                return es_energies;
            }
            return std::nullopt;
        }();
        // ******************************************************************
        if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
            std::cout << " ├state: gs "  << std::endl;
            std::cout << " ││enery = " << *gs_energy << std::endl;
            if (reference_energies) {
                std::cout << " ││enery (reference) = " << reference_energies->get_gs_energy() << std::endl;
            }
        }
        if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
            for (unsigned k_n = 0; k_n < interpreted_program_options.n_sites; k_n++) {
                const auto es_energy = (*es_energies)[k_n];
                std::cout << " ├state: es [k_n = " << k_n << "]" << std::endl;
                std::cout << " ││abs. enery        = " << es_energy << std::endl;
                if (gs_energy) {
                    std::cout << " ││exc. enery        = " << es_energy - *gs_energy << std::endl;
                    if (reference_energies) {
                        std::cout << " ││exc. enery (ref.) = " << reference_energies->get_es_energy(k_n) << std::endl;
                    }
                }
            }
        }
        std::cout << "------------------------------------------" << std::endl;
    } catch (std::exception& e) {
        std::cerr << "[ERROR  ] Abnormal termination!" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
