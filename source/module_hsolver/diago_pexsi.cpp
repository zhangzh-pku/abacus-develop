#include "c_pexsi_interface.h"
#include "diago_pexsi.h"

namespace hsolver
{

void DiagoPexsi::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in)
{
    ModuleBase::WARNING_QUIT("DiagoPexsi", "Pexsi is not completed");
}
void DiagoPexsi::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in)
{
    ModuleBase::WARNING_QUIT("DiagoPexsi", "Pexsi is not completed");
}

}