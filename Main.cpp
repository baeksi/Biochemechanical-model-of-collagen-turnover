// ******************* Problem description ********************************************
// Biochemomechanical approach to model subcellular synthesis of procollagen as well as
// its transition from an intermediate state of self-assembled microfibrils to mature 
// cross-linked fibers with mechano-regulated removal
//
// Theoretical reference: "A Biochemomechanical Model of Collagen Turnover in Arterial 
//			   Adaptations to Hemodynamic Loading",by H.Getachew, J.D.Humphrey, S.Baek		
//
//********************* Creators*******************************************************
// Originally written by Dr. Baek (2009?), modified by Hailu Getachew (2017) 
// at Michigan State University
// 
//********************* Current *******************************************************
// Modified by Vasilina Filonova (VF), on 12-06-2017, at University of Michigan
// Main changes:
//		- rewriten from C to C++
//		- organized to modulus and functions
//		- added or removed some comments
//*************************************************************************************

#include	"GrowthRemodeling_mod.h" 

//------------------- function declaration-----------------------------
	void OutputParamInFile (const double,const double); 
	
// ======================================= MAIN program =======================================
int  main()
{
// variables:
	// km2	- scale constant in vasoactive SMC stress
	// r_h	- homeostatic radius
	// ri_h - inner radius of the wall
	// dt	- time increment (days)
	// ----set local main variables--------------
	double		dt,km2,r_h,ri_h; 
	int			final_time = 500; // duration of time for computation (days)
	clock_t		startcputime, endcputime;
	double		cputime;
	
	ofstream outStreamT; 
	PointGR	point;
	
	// start clock for cpu and wall time
	startcputime = clock();

	//----------- get values needed in main()------------
	km2  = Kg_sh/3;
	ri_h = ro_h - h_h;  
	r_h  = 0.5 * (ri_h + ro_h); 
	dt = 1.0 /time_step ;
		
	//VF:-----output parameters into file------
	OutputParamInFile(km2,r_h);
	
	// print to screen
	cout << scientific;
	cout << "New code..."<< endl;
	cout << "Kg=" << Kg << ", Kg_mu=" << Kg_mu << endl;

	//----------- get homeostatic values ------------------------------------------------
	//update:CIm_p,mPm_h,mPm_p,Mc,Mc_h,Me,Mm,Mm_h,total_M, 
	//		CI_p[],M[],mP_h[],mP_p[],T1_c[],T2_c[],T2_m[],*pk1_p,*pk2_p,*pk3_p,*pm_p,*r_p
	point.ComputeHomeostaticValues(dt,r_h); 
	// ----------------------------------------------------------------------------------
	
	// open file ready for ouput
	string filename = "ProblemResults_f08.txt";
	outStreamT.open(filename.c_str());
	if (outStreamT.fail( )) {
      cerr << "error: open file for output failed!" << endl;
      abort();
	}

	// ---------------- start first loop --------------------------
	//use:   Mc_h,Mm_h,mPm_h,Me,mP_h,
	//update: CIm_p,mPm_p,Mc,Mm,total_M,
	//		  CI_p,mP_p,M,T1_c,T2_c,T2_m,pk1_p,pk2_p,pk3_p,pm_p,r_p
	point.FirstLoop_SolveGrowthRemodeling(outStreamT,dt,km2,r_h,ri_h,
										 final_time);
	// -------------------------------------------------------------
		
	// print to screen
	cout << "end\n";

	//close file
	outStreamT.close( );

	//get cput time
	endcputime = clock();
	cputime = (endcputime - startcputime)/CLOCKS_PER_SEC; 
	cout << "new code CPU time is " << cputime << " seconds" << endl;

	return 0; 
}	// end of main
	
//==================================================== FUNCTION ====================================
//**********************************************************
void OutputParamInFile (const double km2,const double r_h) 
//**********************************************************
// print homeoustatic values and problem parameters into a file
// Paremeters: 2
//		In (const): km2,r_h;
//
	{
		// declare ofstream class variable
		ofstream outStream; 	
		
		// open file for parameters output
		outStream.open("ProblemParameters.info"); 
		if (outStream.fail( ))
		{
			cout << "Output file opening failed.\n";
			abort(); 
		}
		
		// write paramaters
		outStream << scientific;
		outStream << "The homeostatic value for stress in media includes active one\n";
		outStream << "k_2=" << k_2 <<", mu_2=" << mu_2 << ", beta_1=" 
				  << beta_1 << ", zeta_c=" << zeta_c << ", K_act=" << K_act << endl;
		outStream << "Kg_mu=" << Kg_mu << ", Kg=" << Kg << ", km2=" << km2 << ", Kg_sh=" << Kg_sh << endl;
		outStream << "Homeostatics stress = " << stress_h << endl;
		outStream << "Gch=" << Gc_h << ", Gmh=" << Gm_h << ", Gez=" << Ge_1 << ", Get=" << Ge_2 << endl;
		outStream << "L0=" << Lambda_0 << ", Lm=" << Lambda_M << ", ro_h=" << ro_h<< ", r_h=" << r_h << ", h_h=" << h_h << endl;
		outStream << "ce=" << ce << ", cc1=" << cc1 << ", cc2=" << cc2 << endl;
		outStream << "cm1=" << cm1 << ", cm2=" << cm2 << ", cm3=" << cm3 << endl;
		outStream << "phi_k0=[" << phi_k0[0] << ", " << phi_k0[1] << ", " << phi_k0[2] << ", " << phi_k0[3] << "]\n";
		outStream << "P_h=" << P_h << ", S_basal=" << S_basal << ", alpha_h=" << alpha_h << endl;
		
		//close file
		outStream.close( );
		
	}

//*************** end of main program*********************************************************************