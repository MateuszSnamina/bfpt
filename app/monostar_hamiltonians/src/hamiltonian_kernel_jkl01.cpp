#include <monostar_hamiltonians/hamiltonian_kernel_jkl01.hpp>

#include <cmath>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12,123,1234}_jkl01                  ##
// #######################################################################

namespace monostar_hamiltonians {

chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_jkl01() {
    using namespace monostar_system;
    using OnDiagInfoType = std::map<chainkernel::StateKernel1<monostar_system::MonostarSiteStateTrait>, double>;
    using OffDiagInfoType = std::multimap<chainkernel::StateKernel1<monostar_system::MonostarSiteStateTrait>, chainkernel::CoupleInfoKernel1<MonostarSiteStateTrait>>;
    OnDiagInfoType on_diag_info;
    OffDiagInfoType half_off_diag_info;
    return chainkernel::OperatorKernel1<MonostarSiteStateTrait>{on_diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_jkl01([[maybe_unused]] const HamiltonianParamsJkl01& params){
    return prepare_hamiltonian_kernel_1_jkl01();
}

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_jkl01() {
    using namespace monostar_system;
    using OnDiagInfoType = std::map<chainkernel::StateKernel12<monostar_system::MonostarSiteStateTrait>, double>;
    using OffDiagInfoType = std::multimap<chainkernel::StateKernel12<monostar_system::MonostarSiteStateTrait>, chainkernel::CoupleInfoKernel12<MonostarSiteStateTrait>>;
    OnDiagInfoType on_diag_info;
    OffDiagInfoType half_off_diag_info;
    return chainkernel::OperatorKernel12<MonostarSiteStateTrait>{on_diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_jkl01([[maybe_unused]] const HamiltonianParamsJkl01& params){
    return prepare_hamiltonian_kernel_12_jkl01();
}

chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_123_jkl01(double L, double L_1) {
    using namespace monostar_system;
    using OnDiagInfoType = std::map<chainkernel::StateKernel123<monostar_system::MonostarSiteStateTrait>, double>;
    using OffDiagInfoType = std::multimap<chainkernel::StateKernel123<monostar_system::MonostarSiteStateTrait>, chainkernel::CoupleInfoKernel123<MonostarSiteStateTrait>>;
    OnDiagInfoType on_diag_info{
        {{gs, gs, gs}, L},
        {{es, es, es}, L},
        {{gs, es, gs}, L+2*L_1},
        {{es, gs, es}, L+2*L_1},
        {{gs, gs, es}, L+L_1},
        {{es, es, gs}, L+L_1},
        {{gs, es, es}, L+L_1},
        {{es, gs, gs}, L+L_1}
    };
    OffDiagInfoType half_off_diag_info;
    return chainkernel::OperatorKernel123<MonostarSiteStateTrait>{on_diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_123_jkl01(const HamiltonianParamsJkl01& params) {
    return prepare_hamiltonian_kernel_123_jkl01(params.get_L(), params.get_L_1());
}

chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1234_jkl01(
        double J, double J_0, double J_1,
        double K, double K_0, double K_1){
    using namespace monostar_system;
    using OnDiagInfoType = std::map<chainkernel::StateKernel1234<monostar_system::MonostarSiteStateTrait>, double>;
    using OffDiagInfoType = std::multimap<chainkernel::StateKernel1234<monostar_system::MonostarSiteStateTrait>, chainkernel::CoupleInfoKernel1234<MonostarSiteStateTrait>>;
    OnDiagInfoType on_diag_info{
        {{gs, gs, gs, gs}, -0.25*J+K},
        {{es, es, es, es}, -0.25*J+K},
        {{gs, gs, es, gs}, +0.25*J+K +0.25*J_0+K_0 +0.25*J_1+K_1},
        {{es, es, gs, es}, +0.25*J+K +0.25*J_0+K_0 +0.25*J_1+K_1},
        {{gs, es, gs, gs}, +0.25*J+K +0.25*J_0+K_0 +0.25*J_1+K_1},
        {{es, gs, es, es}, +0.25*J+K +0.25*J_0+K_0 +0.25*J_1+K_1},
        {{gs, es, es, gs}, -0.25*J+K -0.25*2*J_1+2*K_1},
        {{es, gs, gs, es}, -0.25*J+K -0.25*2*J_1+2*K_1},
        {{gs, gs, gs, es}, -0.25*J+K -0.25*J_1+K_1},
        {{es, es, es, gs}, -0.25*J+K -0.25*J_1+K_1},
        {{gs, gs, es, es}, +0.25*J+K +0.25*J_0+K_0},
        {{es, es, gs, gs}, +0.25*J+K +0.25*J_0+K_0},
        {{gs, es, gs, es}, +0.25*J+K +0.25*J_0+K_0 +0.25*2*J_1+2*K_1},
        {{es, gs, es, gs}, +0.25*J+K +0.25*J_0+K_0 +0.25*2*J_1+2*K_1},
        {{gs, es, es, es}, -0.25*J+K -0.25*J_1+K_1},
        {{es, gs, gs, gs}, -0.25*J+K -0.25*J_1+K_1}
    };
    OffDiagInfoType half_off_diag_info{
        {{gs, gs, gs, gs}, {{gs, es, es, gs}, 0.5 * J + 0.5 * J_1}},
        {{gs, gs, gs, es}, {{gs, es, es, es}, 0.5 * J + 0.5 * J_1}},
        {{es, gs, gs, gs}, {{es, es, es, gs}, 0.5 * J + 0.5 * J_1}},
        {{es, gs, gs, es}, {{es, es, es, es}, 0.5 * J + 0.5 * J_1}}
    };
    return chainkernel::OperatorKernel1234<MonostarSiteStateTrait>{on_diag_info, half_off_diag_info};
}

chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1234_jkl01(const HamiltonianParamsJkl01& params){
    return prepare_hamiltonian_kernel_1234_jkl01(
                params.get_J(), params.get_J_0(), params.get_J_1(),
                params.get_K(), params.get_K_0(), params.get_K_1());
}

}  // end of namespace monostar_hamiltonians
