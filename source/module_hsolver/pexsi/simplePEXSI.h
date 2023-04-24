#include <mpi.h>
// a simple interface for calling pexsi with 2D block cyclic distributed matrix
int simplePEXSI(MPI_Comm comm_PEXSI, MPI_Comm comm_2D, MPI_Group group_2D, const int blacs_ctxt,  // communicator parameters
                const int size, const int nblk, const int nrow, const int ncol, char LAYOUT, // input matrix parameters
                double* H, double* S,                 // input matrices
                const double nElectronExact, const std::string PexsiOptionFile,        // pexsi parameters file
                double*& DM, double*& EDM,      // output matrices
                double& totalEnergyH, double& totalEnergyS, double& totalFreeEnergy);