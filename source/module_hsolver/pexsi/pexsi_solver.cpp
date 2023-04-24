#include "pexsi_solver.h"

#include "module_base/global_variable.h"
#include "simplePEXSI.h"

#include <cstring>
#include <mpi.h>

PEXSI_Solver::PEXSI_Solver(const int blacs_text,
                           const int nb,
                           const int nrow,
                           const int ncol,
                           const double* h,
                           const double* s,
                           double* DM,
                           double* EDM,
                           double& totalEnergyH,
                           double& totalEnergyS,
                           double& totalFreeEnergy)
{
    this->blacs_text = blacs_text;
    this->nb = nb;
    this->nrow = nrow;
    this->ncol = ncol;
    this->h = new double[nrow * ncol];
    this->s = new double[nrow * ncol];
    std::memcpy(this->h, h, nrow * ncol * sizeof(double));
    std::memcpy(this->s, s, nrow * ncol * sizeof(double));
    this->DM = new double[nrow * ncol];
    this->EDM = new double[nrow * ncol];
    this->totalEnergyH = 0.0;
    this->totalEnergyS = 0.0;
    this->totalFreeEnergy = 0.0;
}

int PEXSI_Solver::solve()
{
    extern MPI_Comm DIAG_WORLD;
    extern MPI_Comm GRID_WORLD;
    extern MPI_Group GRID_GROUP;
    return simplePEXSI(MPI_COMM_WORLD,
                       MPI_COMM_WORLD,
                       GRID_GROUP,
                       this->blacs_text,
                       GlobalV::NLOCAL,
                       this->nb,
                       this->nrow,
                       this->ncol,
                       'C',
                       this->h,
                       this->s,
                       GlobalV::nelec,
                       "PEXSIOPTION",
                       this->DM,
                       this->EDM,
                       this->totalEnergyH,
                       this->totalEnergyS,
                       this->totalFreeEnergy);
}