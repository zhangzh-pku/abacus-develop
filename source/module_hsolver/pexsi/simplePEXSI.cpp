// use PEXSI to solve a Kohn-Sham equation
// the H and S matrices are given by 2D block cyclic distribution
// the Density Matrix and Energy Density Matrix calculated by PEXSI are transformed to 2D block cyclic distribution
// #include "mpi.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <memory>
#include <mpi.h>
#include "c_pexsi_interface.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "DistCCSMatrix.h"
#include "DistBCDMatrix.h"
#include "DistMatrixTransformer.h"

inline void strtolower(char *sa, char *sb)
{
    char c;
    int len = strlen(sa);
    for (int i = 0; i < len; i++)
    {
        c = sa[i];
        sb[i] = tolower(c);
    }
    sb[len] = '\0';
}

inline void setDefaultOption(int* int_para, double* double_para)
{
    // options.spin=2;
    double_para[0]=2;
    // options.gap=0;
    double_para[2]=0;
    // ZERO_Limit=DBL_MIN;
    double_para[11]=DBL_MIN;
    // options.matrixType=0;
    int_para[3]=0;
    // options.solver=1;
    int_para[6]=1;
    // options.ordering=0;
    int_para[8]=0;
    // options.rowOrdering=0;
    int_para[9]=0;
    // options.symmetric=0;
    int_para[11]=0;
    // options.transpose=0;
    int_para[12]=0;
    // options.nPoints=2;
    int_para[14]=2;
    // options.verbosity=1;
    int_para[15]=1;
}

int loadPEXSIOption(MPI_Comm comm, const std::string PexsiOptionFile, PPEXSIOptions& options, int& numProcessPerPole, double& ZERO_Limit)
{

    // temp variable arrays read from conf file and will be bcast to all processors

    // parameter array of type int,
    //  0: numPole
    //  1: isInertiaCount
    //  2: maxPEXSIIter
    //  3: matrixType
    //  4: isSymbolicFactorize
    //  5: isConstructCommPattern
    //  6: solver
    //  7: symmetricStorage
    //  8: ordering
    //  9: rowOrdering
    // 10: npSymbFact
    // 11: symmetric
    // 12: transpose
    // 13: method
    // 14: nPoints
    // 15: verbosity
    // 16: numProcessPerPole
    int int_para[17];

    // parameter array of type double
    //  0: spin
    //  1: temperature
    //  2: gap
    //  3: deltaE
    //  4: muMin0
    //  5: muMax0
    //  6: mu0
    //  7: muInertiaTolerance
    //  8: muInertiaExpansion
    //  9: muPEXSISafeGuard
    // 10: numElectronPEXSITolerance
    // 11: ZERO_Limit
    double double_para[12];
    int myid;
    MPI_Comm_rank(comm, &myid);
    if(myid==0)
    {
        std::ifstream ifs(PexsiOptionFile.c_str());
        if(! ifs)
        {
            return 1;
        }
        setDefaultOption(int_para, double_para);

        ifs.clear();
        ifs.seekg(0);

        char key[128];
        char lowercase_key[128];
        const int LINE_LINGTH=1024;
        char unused_string[LINE_LINGTH];

        while(ifs.good())
        {
            ifs >> key;
            //~ cout<<"readin word is: "<<key<<endl;
            strtolower(key, lowercase_key);
            if(strcmp("spin", lowercase_key)==0)
            {
                //~ ifs>>options.spin;
                ifs>>double_para[0];
                //~ cout<<"double_para[0]: "<<key<<" = "<<double_para[0]<<endl;
            }
            else if(strcmp("temperature", lowercase_key)==0)
            {
                //~ ifs>>options.temperature;
                ifs>>double_para[1];
                //~ cout<<"double_para[1]: "<<key<<" = "<<double_para[1]<<endl;
            }
            else if(strcmp("gap", lowercase_key)==0)
            {
                //~ ifs>>options.gap;
                ifs>>double_para[2];
                //~ cout<<"double_para[2]: "<<key<<" = "<<double_para[2]<<endl;
            }
            else if(strcmp("deltae", lowercase_key)==0)
            {
                //~ ifs>>options.deltaE;
                ifs>>double_para[3];
                //~ cout<<"double_para[3]: "<<key<<" = "<<double_para[3]<<endl;
            }
            else if(strcmp("numpole", lowercase_key)==0)
            {
                //~ ifs>>options.numPole;
                ifs>>int_para[0];
                //~ cout<<"int_para[0]: "<<key<<" = "<<int_para[0]<<endl;
            }
            else if(strcmp("isinertiacount", lowercase_key)==0)
            {
                //~ ifs>>options.isInertiaCount;
                ifs>>int_para[1];
                //~ cout<<"int_para[1]: "<<key<<" = "<<int_para[1]<<endl;
            }
            else if(strcmp("maxpexsiiter", lowercase_key)==0)
            {
                //~ ifs>>options.maxPEXSIIter;
                ifs>>int_para[2];
                //~ cout<<"int_para[2]: "<<key<<" = "<<int_para[2]<<endl;
            }
            else if(strcmp("mumin0", lowercase_key)==0)
            {
                //~ ifs>>options.muMin0;
                ifs>>double_para[4];
                //~ cout<<"double_para[4]: "<<key<<" = "<<double_para[4]<<endl;
            }
            else if(strcmp("mumax0", lowercase_key)==0)
            {
                //~ ifs>>options.muMax0;
                ifs>>double_para[5];
                //~ cout<<"double_para[5]: "<<key<<" = "<<double_para[5]<<endl;
            }
            else if(strcmp("mu0", lowercase_key)==0)
            {
                //~ ifs>>options.mu0;
                ifs>>double_para[6];
                //~ cout<<"double_para[6]: "<<key<<" = "<<double_para[6]<<endl;
            }
            else if(strcmp("muinertiatolerance", lowercase_key)==0)
            {
                //~ ifs>>options.muInertiaTolerance;
                ifs>>double_para[7];
                //~ cout<<"double_para[7]: "<<key<<" = "<<double_para[7]<<endl;
            }
            else if(strcmp("muinertiaexpansion", lowercase_key)==0)
            {
                //~ ifs>>options.muInertiaExpansion;
                ifs>>double_para[8];
                //~ cout<<"double_para[8]: "<<key<<" = "<<double_para[8]<<endl;
            }
            else if(strcmp("mupexsisafeguard", lowercase_key)==0)
            {
                //~ ifs>>options.muPEXSISafeGuard;
                ifs>>double_para[9];
                //~ cout<<"double_para[9]: "<<key<<" = "<<double_para[9]<<endl;
            }
            else if(strcmp("numelectronpexsitolerance", lowercase_key)==0)
            {
                //~ ifs>>options.numElectronPEXSITolerance;
                ifs>>double_para[10];
                //~ cout<<"double_para[10]: "<<key<<" = "<<double_para[10]<<endl;
            }
            else if(strcmp("zero_limit", lowercase_key)==0)
            {
                ifs>>double_para[11];
            }
            else if(strcmp("matrixtype", lowercase_key)==0)
            {
                //~ ifs>>options.matrixType;
                ifs>>int_para[3];
                //~ cout<<"int_para[3]: "<<key<<" = "<<int_para[3]<<endl;
            }
            else if(strcmp("issymbolicfactorize", lowercase_key)==0)
            {
                //~ ifs>>options.isSymbolicFactorize;
                ifs>>int_para[4];
                //~ cout<<"int_para[4]: "<<key<<" = "<<int_para[4]<<endl;
            }
            else if(strcmp("isconstructcommpattern", lowercase_key)==0)
            {
                //~ ifs>>options.isConstructCommPattern;
                ifs>>int_para[5];
                //~ cout<<"int_para[5]: "<<key<<" = "<<int_para[5]<<endl;
            }
            else if(strcmp("solver", lowercase_key)==0)
            {
                //~ ifs>>options.solver;
                ifs>>int_para[6];
                //~ cout<<"int_para[6]: "<<key<<" = "<<int_para[6]<<endl;
            }
            else if(strcmp("symmetricstorage", lowercase_key)==0)
            {
                //~ ifs>>options.symmetricStorage;
                ifs>>int_para[7];
                //~ cout<<"int_para[7]: "<<key<<" = "<<int_para[7]<<endl;
            }
            else if(strcmp("ordering", lowercase_key)==0)
            {
                //~ ifs>>options.ordering;
                ifs>>int_para[8];
                //~ cout<<"int_para[8]: "<<key<<" = "<<int_para[8]<<endl;
            }
            else if(strcmp("rowordering", lowercase_key)==0)
            {
                //~ ifs>>options.rowOrdering;
                ifs>>int_para[9];
                //~ cout<<"int_para[9]: "<<key<<" = "<<int_para[9]<<endl;
            }
            else if(strcmp("npsymbfact", lowercase_key)==0)
            {
                //~ ifs>>options.npSymbFact;
                ifs>>int_para[10];
                //~ cout<<"int_para[10]: "<<key<<" = "<<int_para[10]<<endl;
            }
            else if(strcmp("symmetric", lowercase_key)==0)
            {
                //~ ifs>>options.symmetric;
                ifs>>int_para[11];
                //~ cout<<"int_para[11]: "<<key<<" = "<<int_para[11]<<endl;
            }
            else if(strcmp("transpose", lowercase_key)==0)
            {
                //~ ifs>>options.transpose;
                ifs>>int_para[12];
                //~ cout<<"int_para[12]: "<<key<<" = "<<int_para[12]<<endl;
            }
            else if(strcmp("method", lowercase_key)==0)
            {
                //~ ifs>>options.method;
                ifs>>int_para[13];
                //~ cout<<"int_para[13]: "<<key<<" = "<<int_para[13]<<endl;
            }
            else if(strcmp("npoints", lowercase_key)==0)
            {
                //~ ifs>>options.nPoints;
                ifs>>int_para[14];
                //~ cout<<"int_para[14]: "<<key<<" = "<<int_para[14]<<endl;
            }
            else if(strcmp("verbosity", lowercase_key)==0)
            {
                //~ ifs>>options.verbosity;
                ifs>>int_para[15];
                //~ cout<<"int_para[15]: "<<key<<" = "<<int_para[15]<<endl;
            }
            else if(strcmp("numprocessperpole", lowercase_key)==0)
            {
                //~ ifs>>options.verbosity;
                ifs>>int_para[16];
                //~ cout<<"int_para[16]: "<<key<<" = "<<int_para[16]<<endl;
            }
            else
            {
                if(key[0] == '#' || key[0] == '/')
                {
                    ifs.getline(unused_string, LINE_LINGTH);
                }
                else
                {
                    std::cout<<" THE PARAMETER NAME '" << key << "' IS NOT USED!" << std::endl;
                    return 1;
                }
            }
        }
    }

    // broadcast all options
    MPI_Bcast(int_para, 17, MPI_INT, 0, comm);
    MPI_Bcast(double_para, 12, MPI_DOUBLE, 0, comm);

    // setup PEXSI options from int_para and double_para
    options.numPole=int_para[0];
    options.isInertiaCount=int_para[1];
    options.maxPEXSIIter=int_para[2];
    options.matrixType=int_para[3];
    options.isSymbolicFactorize=int_para[4];
    options.isConstructCommPattern=int_para[5];
    options.solver=int_para[6];
    options.symmetricStorage=int_para[7];
    options.ordering=int_para[8];
    options.rowOrdering=int_para[9];
    options.npSymbFact=int_para[10];
    options.symmetric=int_para[11];
    options.transpose=int_para[12];
    options.method=int_para[13];
    options.nPoints=int_para[14];
    options.verbosity=int_para[15];
    numProcessPerPole=int_para[16];

    options.spin=double_para[0];
    options.temperature=double_para[1];
    options.gap=double_para[2];
    options.deltaE=double_para[3];
    options.muMin0=double_para[4];
    options.muMax0=double_para[5];
    options.mu0=double_para[6];
    options.muInertiaTolerance=double_para[7];
    options.muInertiaExpansion=double_para[8];
    options.muPEXSISafeGuard=double_para[9];
    options.numElectronPEXSITolerance=double_para[10];
    ZERO_Limit=double_para[11];

    return 0;
}

void splitNProc2NProwNPcol(const int NPROC, int& nprow, int& npcol)
{
    int integral_part = (int)sqrt(NPROC);
    if(NPROC%integral_part == 0)
    {
        nprow=integral_part;
        npcol=NPROC/integral_part;
    }
    else
    {
        int flag;
        int i;
        int low=pow(integral_part,2);
        int high=pow(integral_part+1,2);
        if( (NPROC-low) >= (high-NPROC))
        {
            flag=integral_part+1;
        }
        else
        {
            flag=integral_part;
        }
        for(i=flag; i>0; ++i)
        {
            if(NPROC%i == 0) break;
        }
        nprow=i;
        npcol=NPROC/i;
    }
}

int simplePEXSI(MPI_Comm comm_PEXSI, MPI_Comm comm_2D, MPI_Group group_2D, const int blacs_ctxt,  // communicator parameters
                const int size, const int nblk, const int nrow, const int ncol, char LAYOUT,  // matrix parameters
                double* H, double* S,                 // input matrices
                const double numElectronExact, const std::string PexsiOptionFile,         // pexsi parameters file
                double*& DM, double*& EDM,      // output matrices
                double& totalEnergyH, double& totalEnergyS, double& totalFreeEnergy)      // output energy
{
    int out_log=0;
    if(out_log == 1)
    {
        std::stringstream ss;
        int nproc_2D, nproc_PEXSI;
        int myid_2D, myid_PEXSI;
        if(comm_2D != MPI_COMM_NULL)
        {
            MPI_Comm_size(comm_2D, &nproc_2D);
            MPI_Comm_rank(comm_2D, &myid_2D);
            ss.str("");
            ss<<"\tIn 2D comm, myid = "<<myid_2D<<"; nproc = "<<nproc_2D<<std::endl;
            ss<<"H[0] = "<<H[0]<<", H["<<nrow*ncol-1<<"] = "<<H[nrow*ncol-1]<<std::endl;
            ss<<"S[0] = "<<S[0]<<", S["<<nrow*ncol-1<<"] = "<<S[nrow*ncol-1]<<std::endl;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss.str());
        }
        if(comm_PEXSI != MPI_COMM_NULL)
        {
            MPI_Comm_size(comm_PEXSI, &nproc_PEXSI);
            MPI_Comm_rank(comm_PEXSI, &myid_PEXSI);
            ss.str("");
            ss<<"\tIn PEXSI comm, myid = "<<myid_PEXSI<<"; nproc = "<<nproc_PEXSI;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss.str());
        }
    }
    if(comm_2D == MPI_COMM_NULL && comm_PEXSI == MPI_COMM_NULL) return 0;
    int myid;
    std::ofstream f_log;
    if(comm_PEXSI != MPI_COMM_NULL)
    {
        MPI_Comm_rank(comm_PEXSI, &myid);
        // for log
        #ifdef _DEBUG
        if(myid<100) log_openfile(myid, f_log);
        #endif
    }

    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
	//DONE(ofs_running,"set up PEXSI parameter, begin");
    // set up PEXSI parameter
    PPEXSIOptions options;
    PPEXSISetDefaultOptions(&options);
    int numProcessPerPole;
    double ZERO_Limit;
    loadPEXSIOption(comm_PEXSI, PexsiOptionFile, options, numProcessPerPole, ZERO_Limit);
    //OUT(ofs_running, "checkpoint01");
    // debug
    #ifdef _DEBUG
    if(comm_PEXSI != MPI_COMM_NULL)
    {
        if(myid<100) log_PEXSIOption(numElectronExact, f_log);
    }
    #endif
    // end debug
    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
	//DONE(ofs_running,"set up PEXSI parameter, finish");

    // set up PEXSI plan
    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "set_up_PEXSI_plan");
    //OUT(ofs_running, "checkpoint02");
    ModuleBase::timer::tick("DiagoPexsi","setup_PEXSI_plan");
    PPEXSIPlan plan;
    int info;
    int outputFileIndex;
    int pexsi_prow, pexsi_pcol;
    ModuleBase::timer::tick("DiagoPexsi","splitNProc2NProwNPcol");
    splitNProc2NProwNPcol(numProcessPerPole, pexsi_prow, pexsi_pcol);
    ModuleBase::timer::tick("DiagoPexsi","splitNProc2NProwNPcol");
    //OUT(ofs_running, "checkpoint03");
    #ifdef _DEBUG
    //if(comm_PEXSI != MPI_COMM_NULL)
    //{
        if(myid<100) log_PEXSIgrid(pexsi_prow, pexsi_pcol, f_log);
    //}
    #endif
    if(myid % (pexsi_prow * pexsi_pcol) == 0)
    {
        outputFileIndex=myid/(pexsi_prow*pexsi_pcol);
    }
    else
    {
        outputFileIndex=-1;
    }
    //OUT(ofs_running, "checkpoint04");
    ModuleBase::timer::tick("DiagoPexsi","PEXSIPlanInit");
    if(comm_PEXSI != MPI_COMM_NULL)
    {
        //OUT(ofs_running, "checkpoint05");
        plan = PPEXSIPlanInitialize(comm_PEXSI, pexsi_prow, pexsi_pcol, outputFileIndex, &info);
        #ifdef _DEBUG
        //OUT(ofs_running, "checkpoint06");
        if(myid<100) log_PEXSIinit(info, f_log);
        //OUT(ofs_running, "checkpoint07");
        #endif
    }
    ModuleBase::timer::tick("DiagoPexsi","PEXSIPlanInit");
    //OUT(ofs_running, "checkpoint08");
    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
    ModuleBase::timer::tick("DiagoPexsi","setup_PEXSI_plan");
	ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"set_up_PEXSI_plan finish");

    // create compressed column storage distribution matrix parameter
    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
	//DONE(ofs_running,"create compressed column storage distribution matrix parameter, begin");
    DistCCSMatrix DST_Matrix(comm_PEXSI, numProcessPerPole, size);
    //OUT(ofs_running, "checkpoint09");
    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
	//DONE(ofs_running,"create compressed column storage distribution matrix parameter, finish");

    #ifdef _DEBUG
    if(comm_PEXSI != MPI_COMM_NULL)
    {
        if(myid<100) log_DSTMatrix(DST_Matrix, f_log);
    }
    #endif

    // create block cyclic distribution matrix parameter
    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
	//DONE(ofs_running,"create block cyclic distribution matrix parameter, begin");
    //OUT(ofs_running, "checkpoint10");
    DistBCDMatrix SRC_Matrix(comm_2D, group_2D, blacs_ctxt, size, nblk, nrow, ncol, LAYOUT);
    //OUT(ofs_running, "checkpoint11");
    #ifdef _DEBUG
    if(comm_PEXSI != MPI_COMM_NULL)
    {
        if(myid<100) log_SRCMatrix(SRC_Matrix, f_log);
    }
    #endif
    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
	//DONE(ofs_running,"create block cyclic distribution matrix parameter, finish");

    double *HnzvalLocal=new double[1];
    double *SnzvalLocal=new double[1];
    double *DMnzvalLocal=new double[1];
    double *EDMnzvalLocal=new double[1];
    double *FDMnzvalLocal=new double[1];
    // transform H and S from 2D block cyclic distribution to compressed column sparse matrix
    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
	ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"transformBCDtoCCS, begin");
    //OUT(ofs_running, "checkpoint12");
    ModuleBase::timer::tick("Diago_LCAO_Matrix","TransMat2PEXSI");
    transformBCDtoCCS(SRC_Matrix, H, S, ZERO_Limit, DST_Matrix, HnzvalLocal, SnzvalLocal);
    ModuleBase::timer::tick("Diago_LCAO_Matrix","TransMat2PEXSI");
    //MPI_Barrier(MPI_COMM_WORLD);
    //LiuXh modify 2021-03-30, add DONE(ofs_running,"xx") for test
	ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"transformBCDtoCCS, finish");
    //OUT(ofs_running, "checkpoint13");
    if(comm_PEXSI != MPI_COMM_NULL)
    {
        // debug
        #ifdef _DEBUG
        if(myid<100) log_DSTparameter(DST_Matrix, HnzvalLocal, f_log);
        #endif
        // end debug

        // Load H and S to PEXSI
        int isSIdentity=0;
        //OUT(ofs_running, "checkpoint14");
        PPEXSILoadRealHSMatrix(plan, options,
                                size,
                                DST_Matrix.nnz, DST_Matrix.nnzLocal,
                                DST_Matrix.numColLocal, DST_Matrix.colptrLocal, DST_Matrix.rowindLocal,
                                HnzvalLocal,
                                isSIdentity,
                                SnzvalLocal,
                                &info);
        //OUT(ofs_running, "checkpoint15");
        #ifdef _DEBUG
        if(myid<100) log_HSload(f_log);
        #endif
        // call PEXSI to solve Kohn-Sham equation
        // PPEXSIDFTDriver2(plan, &options,
                         // numElectronExact,
                         // &muPEXSI,
                         // &numElectronPEXSI,
                         // &numTotalInertiaIter,
                         // &info);
        double mu;
        double nelec;
        double muMinInertia;
        double muMaxInertia;
        int numTotalPEXSIIter;
        int numTotalInertiaIter; // Number of total inertia[out]
        //OUT(ofs_running, "checkpoint16");
        //LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
	    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"PPEXSIDFTDriver, begin");
        ModuleBase::timer::tick("Diago_LCAO_Matrix","PEXSIDFT");
        PPEXSIDFTDriver(
        plan, // PEXSI plan[in]
        options, // PEXSI Options[in]
        numElectronExact, // exact electron number[in]
        &mu, // chemical potential[out]
        &nelec, // number of electrons[out]
        &muMinInertia, // Lower bound for mu after the last inertia[out]
        &muMaxInertia, // Upper bound for mu after the last inertia[out]
        &numTotalInertiaIter, // Number of total inertia[out]
        &numTotalPEXSIIter, // number of total pexsi evaluation procedure[out]
        &info ); // 0: successful; otherwise: unsuccessful
        //LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
        ModuleBase::timer::tick("Diago_LCAO_Matrix","PEXSIDFT");
	    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"PPEXSIDFTDriver, finish");
        //OUT(ofs_running, "checkpoint17");

        // debug
        #ifdef _DEBUG
        if(myid<100) log_PEXSIcalled(mu, nelec, muMinInertia, muMaxInertia, numTotalPEXSIIter, f_log);
        #endif
        // end debug

        // retrieve the results from the plan
        delete[] DMnzvalLocal;
        delete[] EDMnzvalLocal;
        delete[] FDMnzvalLocal;
        DMnzvalLocal=new double[DST_Matrix.nnzLocal];
        EDMnzvalLocal=new double[DST_Matrix.nnzLocal];
        FDMnzvalLocal=new double[DST_Matrix.nnzLocal];
        if(myid < numProcessPerPole)
        {
            PPEXSIRetrieveRealDFTMatrix(
                plan,
                DMnzvalLocal,
                EDMnzvalLocal,
                FDMnzvalLocal,
                &totalEnergyH,
                &totalEnergyS,
                &totalFreeEnergy,
                &info);
            #ifdef _DEBUG
            if(myid<100) log_DM(DST_Matrix, DMnzvalLocal, f_log);
            #endif
        }
        // clean PEXSI
        PPEXSIPlanFinalize(plan, &info);
        #ifdef _DEBUG
        if(myid<100) log_PEXSIFinalized(f_log);
        #endif
    }
    //OUT(ofs_running, "checkpoint18");

    // transform Density Matrix and Energy Density Matrix from compressed column sparse matrix
    // back to 2D block cyclic distribution if neccessary
    if(comm_2D != MPI_COMM_NULL)
    {
		delete[] DM;
		delete[] EDM;
		DM=new double[SRC_Matrix.nrow*SRC_Matrix.ncol];
		EDM=new double[SRC_Matrix.nrow*SRC_Matrix.ncol];
	}
    #ifdef _DEBUG
    //OUT(ofs_running, "checkpoint19");
    if(myid<100) log_DMEDM_in_BCD_allocated(f_log);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    //LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
	ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"transformCCStoBCD, begin");
    ModuleBase::timer::tick("Diago_LCAO_Matrix","TransMAT22D");
    transformCCStoBCD(DST_Matrix, DMnzvalLocal, EDMnzvalLocal, SRC_Matrix, DM, EDM);
    ModuleBase::timer::tick("Diago_LCAO_Matrix","TransMAT22D");
    //LiuXh modify 2021-04-29, add DONE(ofs_running,"xx") for test
	ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"transformCCStoBCD, finish");
    #ifdef _DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
    //OUT(ofs_running, "checkpoint20");
    if(comm_PEXSI != MPI_COMM_NULL)
    {
        if(myid<100) log_DMtransformed(f_log);
        if(myid<100) log_closefile(f_log);
        // output result
        // save local data of DMnzvalLocal
        /*
        ofstream f_DM;
        sprintf(fname,"DM_%2.2d.dat", myid);
        f_DM.open(fname, ios::out);
        for(int i=0; i<SRC_Matrix.nrow; ++i)
        {
            for(int j=0; j<SRC_Matrix.ncol; ++j)
            {
                f_DM<<DM[i*SRC_Matrix.ncol+j]<<"\t";
            }
            f_DM<<"\n";
        }
        f_DM.close();

        // save local data of EDMnzvalLocal
        ofstream f_EDM;
        sprintf(fname,"EDM_%2.2d.dat", myid);
        f_EDM.open(fname, ios::out);
        for(int i=0; i<SRC_Matrix.nrow; ++i)
        {
            for(int j=0; j<SRC_Matrix.ncol; ++j)
            {
                f_EDM<<EDM[i*SRC_Matrix.ncol+j]<<"\t";
            }
            f_EDM<<"\n";
        }
        f_EDM.close();
        */
    }
    #endif
	MPI_Barrier(MPI_COMM_WORLD);
    //OUT(ofs_running, "checkpoint21");
	MPI_Barrier(MPI_COMM_WORLD);
    delete[] DMnzvalLocal;
    delete[] EDMnzvalLocal;
    delete[] FDMnzvalLocal;
    delete[] HnzvalLocal;
    delete[] SnzvalLocal;
	MPI_Barrier(MPI_COMM_WORLD);
    //OUT(ofs_running, "checkpoint22");
	//MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}
