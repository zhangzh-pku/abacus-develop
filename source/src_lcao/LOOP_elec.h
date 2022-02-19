#ifndef LOOP_ELEC_H
#define LOOP_ELEC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "src_lcao/local_orbital_charge.h"

class LOOP_elec
{
	public:

	LOOP_elec(){};
	~LOOP_elec(){};

	// mohan add 2021-02-09
    void solve_elec_stru(const int& istep,
        std::vector<ModuleBase::matrix>& wfc_gamma,
        std::vector<ModuleBase::ComplexMatrix>& wfc_k,
        Local_Orbital_Charge &loc);

	private:

	// set matrix and grid integral
	void set_matrix_grid(void);

    void before_solver(const int& istep,
        std::vector<ModuleBase::matrix>& wfc_gamma,
        std::vector<ModuleBase::ComplexMatrix>& wfc_k,
        Local_Orbital_Charge &loc);

    void solver(const int& istep,
        std::vector<ModuleBase::matrix>& wfc_gamma,
        std::vector<ModuleBase::ComplexMatrix>& wfc_k,
        Local_Orbital_Charge &loc);

};

#endif
