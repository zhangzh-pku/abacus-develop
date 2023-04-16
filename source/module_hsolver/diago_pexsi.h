#ifndef DIGAOPEXSI_H
#define DIGAOPEXSI_H

#ifdef  __PEXSI


#include "module_basis/module_ao/parallel_orbitals.h"
#include "diagh.h"

namespace hsolver
{

class DiagoPexsi : public DiagH<double>
{
  public:
    DiagoPexsi(const Parallel_Orbitals* ParaV_in)
    {
      this->ParaV = ParaV_in;
    }
    void diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in) override;
    void diag(hamilt::Hamilt<double>* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in) override;
    const Parallel_Orbitals* ParaV;
    double* DM;
    double* EDM;
    double totalEnergyH;
    double totalEnergyS;
    double totalFreeEnergy;

};

}

#endif
#endif
