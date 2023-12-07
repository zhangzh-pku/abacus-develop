#ifndef PEXSI_Solver_H
#define PEXSI_Solver_H

class PEXSI_Solver
{
  public:
    PEXSI_Solver(const int blacs_text,
                 const int nb,
                 const int nrow,
                 const int ncol,
                 const double* h,
                 const double* s,
                 double* DM,
                 double* EDM,
                 double& totalEnergyH,
                 double& totalEnergyS,
                 double& totalFreeEnergy);
    int solve();
    int blacs_text;
    int nb;
    int nrow;
    int ncol;
    double* h;
    double* s;
    double* DM;
    double* EDM;
    double totalEnergyH;
    double totalEnergyS;
    double totalFreeEnergy;
};

#endif // PEXSI_Solver_H