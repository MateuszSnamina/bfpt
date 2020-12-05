#include <monostar_app/get_model_bits.hpp>

#include <monostar_hamiltonians/hamiltonian_kernel_af_fm.hpp>
#include <monostar_hamiltonians/hamiltonian_reference_energies_af_fm.hpp>
#include <monostar_hamiltonians/hamiltonian_kernel_fo.hpp>
#include <monostar_hamiltonians/hamiltonian_reference_energies_fo.hpp>
#include <monostar_hamiltonians/hamiltonian_kernel_jkl01.hpp>
#include <monostar_hamiltonians/hamiltonian_reference_energies_jkl01.hpp>

#include <monostar_hamiltonians/hamiltonian_params_agile_affo_dacay.hpp>
#include <monostar_hamiltonians/get_orbital_theta.hpp>

namespace monostar_app {

chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
get_prepare_hamiltonian_kernel_1(const InterpretedProgramOptions& interpreted_program_options) {
    switch (interpreted_program_options.model_type) {
        case ModelType::AF:
        case ModelType::FM:
            return monostar_hamiltonians::prepare_hamiltonian_kernel_1_af_fm(interpreted_program_options.hamiltonian_params_af_fm);
        case ModelType::FO: {
            const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(
                interpreted_program_options.hamiltonian_params_fo,
                interpreted_program_options.orbital_theta);
            return monostar_hamiltonians::prepare_hamiltonian_kernel_1_fo(interpreted_program_options.hamiltonian_params_fo, orbital_theta_to_use);
        }
        case ModelType::JKL01:
        case ModelType::AgileAFFO:
            return chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>{};
        default:
            throw std::domain_error("Invalid model_type enum value.");
    }
}

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
get_prepare_hamiltonian_kernel_12(const InterpretedProgramOptions& interpreted_program_options) {
    switch (interpreted_program_options.model_type) {
        case ModelType::AF:
            return monostar_hamiltonians::prepare_hamiltonian_kernel_12_af(interpreted_program_options.hamiltonian_params_af_fm);
        case ModelType::FM:
            return monostar_hamiltonians::prepare_hamiltonian_kernel_12_fm(interpreted_program_options.hamiltonian_params_af_fm);
        case ModelType::FO: {
            const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
            return monostar_hamiltonians::prepare_hamiltonian_kernel_12_fo(interpreted_program_options.hamiltonian_params_fo, orbital_theta_to_use);
        }
        case ModelType::JKL01:
        case ModelType::AgileAFFO:
            return chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>{};
        default:
            throw std::domain_error("Invalid model_type enum value.");
    };
}

chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>
get_prepare_hamiltonian_kernel_123(const InterpretedProgramOptions& interpreted_program_options) {
    switch (interpreted_program_options.model_type) {
        case ModelType::AF:
            return chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>{};
        case ModelType::FM:
            return chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>{};
        case ModelType::FO:
            return chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>{};
        case ModelType::JKL01:
            return monostar_hamiltonians::prepare_hamiltonian_kernel_123_jkl01(interpreted_program_options.hamiltonian_params_jkl01);
        case ModelType::AgileAFFO: {
            const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(
                interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian(),
                interpreted_program_options.orbital_theta,
                interpreted_program_options.average_ss);
            const monostar_hamiltonians::HamiltonianParamsJkl01 decayed_hamiltonian =
                monostar_hamiltonians::dacay_hamiltonian_params_agile_affo(
                    interpreted_program_options.hamiltonian_params_agile_affo,
                    orbital_theta_to_use);
            return monostar_hamiltonians::prepare_hamiltonian_kernel_123_jkl01(decayed_hamiltonian);
        }
        default:
            throw std::domain_error("Invalid model_type enum value.");
    };
}

chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>
get_prepare_hamiltonian_kernel_1234(const InterpretedProgramOptions& interpreted_program_options) {
    switch (interpreted_program_options.model_type) {
        case ModelType::AF:
            return chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>{};
        case ModelType::FM:
            return chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>{};
        case ModelType::FO:
            return chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>{};
        case ModelType::JKL01:
            return monostar_hamiltonians::prepare_hamiltonian_kernel_1234_jkl01(interpreted_program_options.hamiltonian_params_jkl01);
        case ModelType::AgileAFFO: {
            const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(
                interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian(),
                interpreted_program_options.orbital_theta,
                interpreted_program_options.average_ss);
            const monostar_hamiltonians::HamiltonianParamsJkl01 decayed_hamiltonian =
                monostar_hamiltonians::dacay_hamiltonian_params_agile_affo(
                    interpreted_program_options.hamiltonian_params_agile_affo,
                    orbital_theta_to_use);
            return monostar_hamiltonians::prepare_hamiltonian_kernel_1234_jkl01(decayed_hamiltonian);
        }
        default:
            throw std::domain_error("Invalid model_type enum value.");
    };
}

std::vector<utility::Named<arma::cx_mat22>>
get_one_site_metrices_for_average_calculation(const InterpretedProgramOptions& interpreted_program_options) {
    switch (interpreted_program_options.model_type) {
        case ModelType::AF:
        case ModelType::FM:
        case ModelType::JKL01:
        case ModelType::AgileAFFO:
            return monostar_hamiltonians::OneSiteSpinNamedMatrices::site_matrices_for_average_calculations_af_fm();
        case ModelType::FO: {
            const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
            return monostar_hamiltonians::OneSiteOrbitalNamedMatrices::matrices_for_average_calculations(orbital_theta_to_use);
        }
        default:
            throw std::domain_error("Invalid model_type enum value.");
    };
}

std::vector<utility::Named<arma::cx_mat44>>
get_two_site_metrices_for_average_calculation(const InterpretedProgramOptions& interpreted_program_options) {
    switch (interpreted_program_options.model_type) {
        case ModelType::AF:
        case ModelType::JKL01:
        case ModelType::AgileAFFO:
            return monostar_hamiltonians::TwoSitesSpinNamedMatrices::matrices_for_average_calculations_af();
        case ModelType::FM:
            return monostar_hamiltonians::TwoSitesSpinNamedMatrices::matrices_for_average_calculations_fm();
        case ModelType::FO: {
            const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
            return monostar_hamiltonians::TwoSitesOrbitalNamedMatrices::matrices_for_average_calculations(orbital_theta_to_use);
        }
        default:
            throw std::domain_error("Invalid model_type enum value.");
    };
}

std::shared_ptr<monostar_hamiltonians::HamiltonianReferenceEnergies>
get_hamiltonian_reference_energies(const InterpretedProgramOptions& interpreted_program_options) {
    switch (interpreted_program_options.model_type) {
        case ModelType::AF:
            return std::dynamic_pointer_cast<monostar_hamiltonians::HamiltonianReferenceEnergies>(
                std::make_shared<monostar_hamiltonians::HamiltonianReferenceEnergiesAf>(
                    interpreted_program_options.n_sites, interpreted_program_options.hamiltonian_params_af_fm));
        case ModelType::FM:
            return std::dynamic_pointer_cast<monostar_hamiltonians::HamiltonianReferenceEnergies>(
                std::make_shared<monostar_hamiltonians::HamiltonianReferenceEnergiesFm>(
                    interpreted_program_options.n_sites, interpreted_program_options.hamiltonian_params_af_fm));
        case ModelType::FO: {
            const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
            return std::dynamic_pointer_cast<monostar_hamiltonians::HamiltonianReferenceEnergies>(
                std::make_shared<monostar_hamiltonians::HamiltonianReferenceEnergiesFo>(
                    interpreted_program_options.n_sites, interpreted_program_options.hamiltonian_params_fo, orbital_theta_to_use));
        }
        case ModelType::JKL01: {
            return std::dynamic_pointer_cast<monostar_hamiltonians::HamiltonianReferenceEnergies>(
                std::make_shared<monostar_hamiltonians::HamiltonianReferenceEnergiesJkl01>(
                    interpreted_program_options.n_sites, interpreted_program_options.hamiltonian_params_jkl01));
        }
        case ModelType::AgileAFFO: {
            const double orbital_theta_to_use = monostar_hamiltonians::get_orbital_theta(
                interpreted_program_options.hamiltonian_params_agile_affo.get_so_hamiltonian(),
                interpreted_program_options.orbital_theta,
                interpreted_program_options.average_ss);
            const monostar_hamiltonians::HamiltonianParamsJkl01 decayed_hamiltonian =
                monostar_hamiltonians::dacay_hamiltonian_params_agile_affo(
                    interpreted_program_options.hamiltonian_params_agile_affo,
                    orbital_theta_to_use);
            std::make_shared<monostar_hamiltonians::HamiltonianReferenceEnergiesJkl01>(
                interpreted_program_options.n_sites, decayed_hamiltonian);  //TODO check whether decayed_hamiltonian os passed by ref and the ref dangles!
        }
        default:
            throw std::domain_error("Invalid model_type enum value.");
    }
}

}  // end of namespace monostar_app
