// Header file for global constants and for growth and remodeling class
// by Vasilina, 12-06-2017 
//
#ifndef _GrowthRemodeling_mod_H
#define _GrowthRemodeling_mod_H_
#include 	<cmath> 
#include	<cstdlib>
#include	<iomanip>
#include	<iostream>
#include 	<fstream>
#include	<string>
#include    <time.h>
using namespace std;

// UNITS: kg, meters,seconds (or days), 
//		pressure Pa=kg/m/s^2
//		J/kg=m^2/s^2
//		day = 86400 seconds
//		Warning: 
//			-K_act is written in (per days) while ce,cc1... use seconds (J/kg=m^2/s^2)
//				Need to choose one unit for time (either days or seconds) and use consistently!
// ----------------------------PARAMETERS-------------------------------------------------------
// ----Parameters from Table2 [Getachew,Humphrey,Baek] for mouse carotid artery:
//		ro_h (m)			- outer radius
//		h_h					- homeostatic wall thickness
//		Gc_h,Gm_h			- deposition stretch of the collagen and SMC constituent
//		Ge_1,Ge_2			- coeff. in tensor mapping from natural config. of elastin to the comput. ref. config. 
//		shear_h	 (Pa)		- homeostatic value of wall shear stress
//		stress_h (Pa)		- homeostatic value of normal stress
//		phi_e0, phi_m0,phi_k0[4] -  homeostatic mass fractions of elastin, SMC and collagen fibers
//		phi_f				- mass fractions of fluid in the vessel
//		P_h	(Pa)			- intramural homeostatic pressure
//		Lambda_0,Lambda_M	- the SMC stretch of maximum and zero contraction 
//		K_act (per day)		- coefficient for r_act ODE
//		S_basal	 (Pa)		- basal vasoactive tone of SMC
//		alpha_h (degree)	- homeostatic angle for collagne fiber orientation
//		ce (J/kg)			- parameter for strain energy of elastin (is ce1 in Table2) 
//		cc1 (J/kg)			- parameter for strain energy of collagen (is cc2 in Table2) 
//		cc2					- parameter for strain energy of collagen (is cc3 in Table2)
//		cm1, cm2 (J/kg)		- parameters for strain energy of smooth muscle
//		cm3					- parameter for strain energy of smooth muscle
// -----Other:
// VF: Wall parameters:
//		the same mechanical properties for multiple families of locally parallel collagen fibers
//		rho					- wall density (kg/m^3)
// VF: Model constant parameters:
//		beta_1 (nondim.)	- fraction of collagen released to extraceccular space
//		k_2 (per day)		- used to define beta_2 = k_2 /(k_2 + mu_2)
//		Kg_mu				-  multiplication coefficient in the functions of the collagen fiber removal rate (mu_F)
//		Kg					- gain-type parameter, sensitivity parameter for normal stress 
//		Kg_sh				- gain-type parameter, sensitivity parameter for shear stress  
//							//VF: was initialized as zero and then assign that value
//		L1					- principal stretch in axial direction is fixed 
//		mu_2 (per day)		- degradation rate in final state of extracellular collagen
//		MWH					- molecular weight of collagen //??Molecular Weight is taken to be 300 kDa = 4.981617e-22 kg
//		zeta_c				- relative stretch of collagen fibers
// VF: Intergation parameters:
//		a_max				- maximum age of mature collagen //VF: changed from double to int
//		Max_it				- maximum namber of iterations in 2nd and 3rd loops
//		time_step			- number of time steps per day
//		Tol					- tolerance criteria for errors in 2nd and 3rd iteration loops 
//
// --------------------------------VARIABLES---------------------------------------------------------------
// VF: model variables:
//		alpha[]				- angle between direction of collagen family k=1..4 and the basis vector in current configuration
//		CI[4]				- molar concentration on intermediate-state collagen
//		CI_p[4]				- concentration of procollagen (intermediate state of collagen)
//		CIm					- concentration of SMC
//		CIm_p				- concentration of intermediate SMC (as for collagen)
//		dwdL22_m, dwdL22_c,dwdL12_c	- derivatives of strain functions w.r.t stretches
//		h					- wall thickness
//		hc					- collagen thickness
//		hm					- SMC thickness
//		L2					- principal stretch in radial direction
//		Lk2_n				- stretches in collagen fiber used in constitutive relation
//		L2m_act				- circumferential active stretch of SMC (due to vascular smooth muscle tone)
//		Lm2_n2				- circumferential stretch in SMC fiber used in constitutive relation
//		mPm_p				- mass production rate of SMC
//		mP_h[4]				- homeostatic rate of mass production of collagen fibers
//		mP_p[4]				- rate of mass production of collagen fibers
//		mP[4]				- collagen mass rate productin (synthesis rate of intracellular procollagen and proliferation of SMC)
//		M[4]					- mass of idividual collagen fiber family
//		Mc					- mass of total collagen constitutient Mc=sum(M)
//		Mm					- total mass of SMC at current time 
//		Mm_h				- homeostatic total mass of SMC (at zero time)
//		pressure (Pa)		- intramural Pressure as a function of time
//		*pk1_a,*pk2_a,*pk3_a	- time array for (mass rate) * (survival fucntion), from th beginning 
//		*pk1_p,*pk2_p,*pk3_p	-  time array for (mass rate) * (survival fucntion) 
//		*pm_a				-  time array for (mass rate) * (survival fucntion), from th beginning
//		*pm_p				- time array for (mass rate) * (survival fucntion) 
//		r					- instanteneous wall inner radius (lumen radius)
//		r_act				- vessel radius for the active tone of SMC
//		*r_p				- time array of vessel radius
//		S					- function for the active smooth muscle tone
//		shear (Pa)			- shear stress
//		sigma_k				- normal stress for collagen fibers
//		total_M				- total mass of all the wall per unit area
//		t_act				- membrane normal stress due to vascular smooth muscle tone
//		T1_c[],T2_c[]		- Cauchy membrane stress(axial and circumferential) for collagen fibers
//		T2_m[]				- Cauchy membrane stress for SMC
//		zeta0,zeta1,zeta2,zetam	- relative stretch of collagen fibers
//		
//VF: integration variables:
//		final_time (days)	- time duration for computation 
//		num_it2				-
//
// --------------- set GLOBAL constants ---------------------------------------
	const double    alpha_h = 47.569;
	const double    beta_1 = 0.5;
	const double    cc1= 3091.6,	cc2 = 18.850,	ce = 657.89;
	const double    cm1 = 103.34,	cm2 = 11.526,	cm3 = 0.7936;
	const double	Gc_h = 1.07, Ge_1 = 2.19, Ge_2 = 1.64, Gm_h = 1.25;  
	const double	h_h = 0.000021175;
	const double	k_2 = 0.001, K_act = 0.1,Kg = 3,Kg_mu = 1,Kg_sh = 52;
	const double	L1 = 1.0,Lambda_0 = 0.65, Lambda_M = 1.65*0.9;  
	const double	mu_2 = 0.0;
	const double    MWH = 1;
	const double    P_h = 102.1445*133.32237;
	const double	PI = 3.14159265;
	const double    phi_e0 = 0.06, phi_m0 = 0.09, phi_f= 0.7, phi_k0[4] = {0.11137,0.034927,0.42685,0.42685};
	const double	rho = 1050.0,ro_h = 0.0003282;  
	const double    S_basal = 0.862055*1.75; //VF: why 1.75?
	const double    shear_h = 1.5, stress_h = 120000.0;
	const double	Tol = 1.e-11;
	const double	zeta_c = 1;
	const int       a_max = 350, Max_it = 2000, time_step = 10; 

class PointGR // 0D model for axisymmetrical cross-section
{ 
	private:
		int		num_pa;
		// homeostatic vars
		double	CIm_p,mPm_h,mPm_p,Mc,Mc_h,Me,Mm,Mm_h,total_M;
		double	CI_p[4],M[4],mP_h[4],mP_p[4],T1_c[3],T2_c[3],T2_m[3];
		double	*pk1_p,*pk2_p,*pk3_p,*pm_p,*r_p; 
		// loops vars
		double	h,hm,L2,L2m_act,Lk2_n,Lm2_n2,MWH1,pressure,r,r_act,S,shear,t_act;
		double	alpha[4],CI[4],mP[4],sigma_k[4];
		double	*pk1_a,*pk2_a,*pk3_a,*pm_a; 
		
	public:
		//constructor: allocate and initialize class variables
		PointGR();
		//destructor: release allocated memory
		~PointGR();

		// homeostatic
		void	ComputeHomeostaticValues(const double dt,const double r_h);
		double  mu_F(const double age, const double zeta);

		// loops
		void 	FirstLoop_SolveGrowthRemodeling(ofstream & outStreamT,const double dt, const double km2,const double r_h,const double ri_h,
					int & final_time);

		void 	SecondLoop_SolveMassProductionODE(ofstream & outStreamT, const double dt,const double km2,const double r_h,const double ri_h, const double t,
					double & CIm,double & zeta0,double & zeta1,double & zeta2,double & zetam,int & final_time);

		void 	ThirdLoop_SolveStressEquilibrium(ofstream & outStreamT,const double dt,const double km2,const double r_h,const double ri_h,const double t,
					double & hc,double & dwdL22_m,double & dwdL22_c,double & dwdL12_c,int & final_time,int & num_it2);
		
		void 	GetLZetas(const double dt,const double r_h,	      
					double & zeta0,double & zeta1,double & zeta2,double & zetam);

		void 	GetMandPointersP(const double dt,const double zeta0,const double zeta1,const double zeta2,const double zetam,
					double & mPm);
		
		double 	Given_P(const double time);
		
		double 	Given_flow(const double time);
		
		double 	f_mP(const double Mi_0, const double Mi_t,const double Sigma, const double mP_b, const double kkgg, const double tau);
		
		void 	OutputResultsEachTime(ofstream & outStreamT,const double t,const double zeta1,const double T1_c0,const double T2_c0,const double T2_m0,const double Mm);
};

#endif