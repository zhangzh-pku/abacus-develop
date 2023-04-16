#ifdef __PEXSI
#include "diago_pexsi.h"

#include "c_pexsi_interface.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "pexsi/pexsi_solver.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{

void DiagoPexsi::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoPEXSI", "diag");
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);
    MPI_Comm COMM_DIAG = MPI_COMM_WORLD;
    extern MPI_Group GRID_GROUP;
    extern MPI_Group GRID_WORLD;
    extern MPI_Group DIAG_WORLD;
    simplePEXSI(GRID_WORLD,
                DIAG_WORLD,
                GRID_GROUP,
                this->ParaV->blacs_ctxt,
                GlobalV::NLOCAL,
                this->ParaV->nb,
                this->ParaV->nrow,
                this->ParaV->ncol,
                'C',
                h_mat.p,
                s_mat.p,
                GlobalV::nelec,
                "PEXSIOPTION",
                this->DM,
                this->EDM,
                this->totalEnergyH,
                this->totalEnergyS,
                this->totalFreeEnergy);

    ModuleBase::WARNING_QUIT("DiagoPexsi", "Pexsi is not completed");
}
void DiagoPexsi::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoPEXSI", "diag");
    ModuleBase::WARNING_QUIT("DiagoPexsi", "Pexsi is not completed");
}

} // namespace hsolver
#endif