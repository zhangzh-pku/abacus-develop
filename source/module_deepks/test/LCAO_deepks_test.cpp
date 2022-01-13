#include "LCAO_deepks_test.h"

namespace Test_Deepks
{
	Grid_Driver GridD(GlobalV::test_deconstructor,GlobalV::test_grid_driver,GlobalV::test_grid);
}

test_deepks::test_deepks()
{}

test_deepks::~test_deepks()
{}


void test_deepks::check_dstable(void)
{
	OGT.talpha.print_Table_DSR(ORB);
	bool same=this->compare_with_ref("S_I_mu_alpha.dat","S_I_mu_alpha_ref.dat");
	if(!same) std::cout << "S_I_mu_alpha.dat inconsistent" << std::endl;
}

void test_deepks::check_psialpha(void)
{
	std::vector<int> na;
	na.resize(ucell.ntype);
	for(int it=0;it<ucell.ntype;it++)
	{
		na[it] = ucell.atoms[it].na;
	}
	ld.init(ORB,
		ucell.nat,
		ucell.ntype,
		na);

	ld.build_psialpha(GlobalV::FORCE,
		ucell,
		ORB,
		Test_Deepks::GridD,
		ParaO.trace_loc_row,
		ParaO.trace_loc_col,
		OGT);

	ld.check_psialpha(GlobalV::FORCE,
		ucell,
		ORB,
		Test_Deepks::GridD,
		ParaO.trace_loc_row,
		ParaO.trace_loc_col,
		OGT);
	bool same=this->compare_with_ref("psialpha.dat","psialpha_ref.dat");
	if(!same) std::cout << "psialpha.dat inconsistent" << std::endl;
	same=this->compare_with_ref("dpsialpha_x.dat","dpsialpha_x_ref.dat");
	if(!same) std::cout << "dpsialpha_x.dat inconsistent" << std::endl;
	same=this->compare_with_ref("dpsialpha_y.dat","dpsialpha_y_ref.dat");
	if(!same) std::cout << "dpsialpha_y.dat inconsistent" << std::endl;	
	same=this->compare_with_ref("dpsialpha_z.dat","dpsialpha_z_ref.dat");
	if(!same) std::cout << "dpsialpha_z.dat inconsistent" << std::endl;		
}

void test_deepks::read_dm(void)
{
	ifstream ifs("dm");
	dm.create(GlobalV::NLOCAL,GlobalV::NLOCAL);

	for (int mu=0;mu<GlobalV::NLOCAL;mu++)
	{
		for (int nu=0;nu<GlobalV::NLOCAL;nu++)
		{
			double c;
			ifs >> c;
			dm(mu,nu)=c;
		}
	}
}

void test_deepks::read_dm_k(const int nks)
{
	dm_k.resize(nks);
	stringstream ss;
	for(int ik=0;ik<nks;ik++)
	{
        ss.str("");
        ss<<"dm_"<<ik;
        ifstream ifs(ss.str().c_str());
		dm_k[ik].create(GlobalV::NLOCAL,GlobalV::NLOCAL);

		for (int mu=0;mu<GlobalV::NLOCAL;mu++)
		{
			for (int nu=0;nu<GlobalV::NLOCAL;nu++)
			{
				std::complex<double> c;
				ifs >> c;
				dm_k[ik](mu,nu)=c;
			}
		}
	}
}

void test_deepks::check_pdm(void)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->read_dm();
		this->ld.cal_projected_DM(dm,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col);
	}
	else
	{
		this->read_dm_k(kv.nkstot);
		this->ld.cal_projected_DM_k(dm_k,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			kv.kvec_d);		
	}
	this->ld.check_projected_dm();
	bool same=this->compare_with_ref("pdm.dat","pdm_ref.dat");
	if(!same) std::cout << "pdm.dat inconsistent" << std::endl;	
}

void test_deepks::check_gdmx(void)
{
	this->ld.init_gdmx(ucell.nat);
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.cal_gdmx(dm,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col);
	}
	else
	{
		this->ld.cal_gdmx_k(dm_k,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			kv.kvec_d);			
	}
	this->ld.check_gdmx(ucell.nat);
	bool same;
	for(int ia=0;ia<ucell.nat;ia++)
	{
		stringstream ss;
		stringstream ss1;
		ss.str("");
        ss<<"gdmx_"<<ia<<".dat";
		ss1.str("");
        ss1<<"gdmx_"<<ia<<"_ref.dat";
		
		same=this->compare_with_ref(ss.str(),ss1.str());
		if(!same) std::cout << "gdmx_" << ia << ".dat inconsistent" << std::endl;

        ss.str("");
        ss<<"gdmy_"<<ia<<".dat";
		ss1.str("");
        ss1<<"gdmy_"<<ia<<"_ref.dat";
		same=this->compare_with_ref(ss.str(),ss1.str());
		if(!same) std::cout << "gdmy_" << ia << ".dat inconsistent" << std::endl;

        ss.str("");
        ss<<"gdmz_"<<ia<<".dat";
        ss1.str("");
        ss1<<"gdmz_"<<ia<<"_ref.dat";
		same=this->compare_with_ref(ss.str(),ss1.str());
		if(!same) std::cout << "gdmz_" << ia << ".dat inconsistent" << std::endl;
	}	
}

void test_deepks::check_descriptor(void)
{
	this->ld.cal_descriptor();
	this->ld.check_descriptor(ucell);
	bool same=this->compare_with_ref("descriptor.dat","descriptor_ref.dat");
	if(!same) std::cout << "descriptor.dat inconsistent" << std::endl;	
}

void test_deepks::check_gvx(void)
{
	this->ld.cal_gvx(ucell.nat);
	this->ld.check_gvx(ucell.nat);
	bool same;
	for(int ia=0;ia<ucell.nat;ia++)
	{
		stringstream ss;
		stringstream ss1;
		ss.str("");
        ss<<"gvx_"<<ia<<".dat";
		ss1.str("");
        ss1<<"gvx_"<<ia<<"_ref.dat";
		same=this->compare_with_ref(ss.str(),ss1.str());
		if(!same) std::cout << "gvx_" << ia << ".dat inconsistent" << std::endl;

        ss.str("");
        ss<<"gvy_"<<ia<<".dat";
		ss1.str("");
        ss1<<"gvy_"<<ia<<"_ref.dat";
		same=this->compare_with_ref(ss.str(),ss1.str());
		if(!same) std::cout << "gvy_" << ia << ".dat inconsistent" << std::endl;

        ss.str("");
        ss<<"gvz_"<<ia<<".dat";
        ss1.str("");
        ss1<<"gvz_"<<ia<<"_ref.dat";
		same=this->compare_with_ref(ss.str(),ss1.str());
		if(!same) std::cout << "gvz_" << ia << ".dat inconsistent" << std::endl;
	}
}

void test_deepks::check_edelta(void)
{
	this->ld.load_model("model.ptg");
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.allocate_V_delta(ucell.nat,ParaO.nloc);
	}
	else
	{
		this->ld.allocate_V_delta(ucell.nat,ParaO.nloc,kv.nkstot);
	}
	this->ld.cal_gedm(ucell.nat);

	ofstream ofs("E_delta.dat");
	ofs << std::setprecision(10) << this->ld.E_delta << std::endl;
	ofs.close();
	bool same=this->compare_with_ref("E_delta.dat","E_delta_ref.dat");
	if(!same) std::cout << "E_delta.dat inconsistent" << std::endl;
	
	this->ld.check_gedm();
	same=this->compare_with_ref("gedm.dat","gedm_ref.dat");
	if(!same) std::cout << "gedm.dat inconsistent" << std::endl;
}

void test_deepks::check_e_deltabands(void)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.cal_e_delta_band(dm,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			ParaO.nrow);
	}
	else
	{
		this->folding_nnr(kv);
		this->ld.cal_e_delta_band_k(dm_k,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			ParaO.nrow,
			ParaO.ncol);
	}

	ofstream ofs("E_delta_bands.dat");
	ofs << std::setprecision(10) << this->ld.e_delta_band << std::endl;
	ofs.close();
	bool same=this->compare_with_ref("E_delta_bands.dat","E_delta_bands_ref.dat");
	if(!same) std::cout << "E_delta_bands.dat inconsistent" << std::endl;
}

void test_deepks::check_v_delta()
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.add_v_delta(ucell,
            ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
            ParaO.trace_loc_col,
            ParaO.nrow,
            ParaO.ncol);
		this->ld.check_v_delta(ParaO.nrow,ParaO.ncol);
		bool same=this->compare_with_ref("H_V_delta.dat","H_V_delta_ref.dat");
		if(!same) std::cout << "H_V_delta.dat inconsistent" <<std::endl;
	}
	else
	{
		this->cal_nnr();
		this->ld.allocate_V_deltaR(nnr);
		this->ld.add_v_delta_k(ucell,
        	ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			nnr);
		this->ld.check_v_delta_k(nnr);
		bool same=this->compare_with_ref("H_V_deltaR.dat","H_V_deltaR_ref.dat");
		if(!same) std::cout << "H_V_deltaR.dat inconsistent" <<std::endl;
	}
}

void test_deepks::check_f_delta()
{
	ModuleBase::matrix svnl_dalpha;
	svnl_dalpha.create(3,3);
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		ld.cal_f_delta_gamma(dm,
            ucell,
            ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
			ParaO.trace_loc_col,
            1, svnl_dalpha);
	}
	else
	{
		ld.cal_f_delta_k(dm_k,
			ucell,
            ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			kv.kvec_d,
			1,svnl_dalpha);
	}
	ld.check_f_delta(ucell.nat, svnl_dalpha);

	bool same=this->compare_with_ref("F_delta.dat","F_delta_ref.dat");
	if(!same) std::cout << "F_delta.dat inconsistent" << std::endl;
}

bool test_deepks::compare_with_ref(
	const std::string f1,
	const std::string f2)
{

	ifstream file1(f1.c_str());
	ifstream file2(f2.c_str());
	double test_thr=1e-8;

	std::string word1;
	std::string word2;
	while(file1 >> word1)
	{
		file2 >> word2;
		if((word1[0]-'0'>=0 && word1[0]-'0'<10)||word1[0]=='-')
		{
			double num1 = std::stof(word1);
			double num2 = std::stof(word2);
			if(abs(num1-num2)>test_thr) return false;			
		}
		else
		{
			if(word1!=word2) return false;
		}
	}
	return true;
}