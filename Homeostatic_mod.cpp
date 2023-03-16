// Functions of Homeostatic part
// by Vasilina, 12-06-2017 
#include	"GrowthRemodeling_mod.h"

//****************
PointGR::PointGR()
//****************
// costructor: allocate and initialize class variables
{
	int i;

	// get constant value
	num_pa = time_step * a_max + 1;

	// initialization to zero
	CIm_p=0.0; mPm_h=0.0; mPm_p=0.0; Mc=0.0; Mc_h=0.0; Me=0.0; Mm=0.0; Mm_h=0.0; total_M=0.0;
	h = 0.0; hm = 0.0; L2 = 0.0; L2m_act = 0.0; Lk2_n = 0.0; Lm2_n2 = 0.0; MWH1 = 0.0; 
	pressure = 0.0; r = 0.0; r_act = 0.0; S = 0.0; shear = 0.0;	t_act = 0.0;

	for (i = 0; i <4; ++i)
	{
		CI_p[i]=0.0; mP_h[i]=0.0; mP_p[i]=0.0; M[i]=0.0;
		alpha[i]=0.0; sigma_k[i]=0.0; CI[i]=0.0; mP[i]=0.0;
	}
	for (i = 0; i <3; ++i)
	{
		T1_c[i]=0.0;T2_c[i]=0.0;T2_m[i]=0.0;
	}
	//dynamic memory allocation for pointers
	pk1_p = new double[num_pa];
	pk2_p = new double[num_pa];
	pk3_p = new double[num_pa];
	pm_p  = new double[num_pa];
	r_p   = new double[num_pa];
	pk1_a = new double[num_pa]; 
	pk2_a = new double[num_pa];
	pk3_a = new double[num_pa];
	pm_a  = new double[num_pa]; 

	// initialize to zeros
	for (i=0; i<num_pa; i++)
	{
		*(pk1_p+i) = 0.0;
		*(pk2_p+i) = 0.0;
		*(pk3_p+i) = 0.0;
		*(pm_p+i) = 0.0;
		*(r_p+i) = 0.0;
		*(pk1_a+i) = 0.0;
		*(pk2_a+i) = 0.0;
		*(pk3_a+i) = 0.0;
		*(pm_a+i) = 0.0;
	}
}
//

//****************
PointGR::~PointGR()
//****************
// destructor: release allocated memory
{
	delete [] pk1_a;
	delete [] pk2_a;
	delete [] pk3_a;
	delete [] pm_a;
	delete [] pk1_p;
	delete [] pk2_p;
	delete [] pk3_p;
	delete [] pm_p;
	delete [] r_p;
}
//

//************************************************************************************************
void PointGR::ComputeHomeostaticValues(const double dt,const double r_h) 
//************************************************************************************************
// Compute some homeostatic values homeostatic values for pk, CI, a_mean, a_half
// Internal variables and parameters:
//		alpha_R[4]		- angle between direction of collagen family k=1..4 and the basis vector E1 in the referemce configuration
//		a_mean (day)	- the mean age of the mature collagen, mean life-time
//		beta_2 (1/day) - collagen conversion from intermediate state (self-assembled microfibrils) to final (fibrillar collagen with mature cross-links) parametrized state
//
// Parameters: 2
//		In (const): dt,r_h 
// Update implicitly within a class: CIm_p,mPm_h,mPm_p,Mc,Mc_h,Me,Mm,Mm_h,total_M,
//  				  CI_p[],M[],mP_h[],mP_p[],T1_c[],T2_c[],T2_m[],*pk1_p,*pk2_p,*pk3_p,*pm_p,*r_p
//
{ 
	// internal variables
	double	phi_c0;
	double  a_mean,a_half,beta_2;
	double 	dummy1, dummy2, dummy3, dummy4;
	//double 	mP[4],CI[4]; 	// VF: internal arrays
	double 	mPh[4],CIh[4]; 	// VF: rename  to differe from mP,CI class variables
	double  sigma_k_p[4],alpha_R[4]; //VF: won't be used
	double	CIm,mPm;
	int		flag,i ;
		
	//-----initialization of internal variables
	a_mean = 0;	flag = 1;
		  
	total_M = h_h * rho;
	phi_c0 = 1.0 - phi_e0 - phi_m0 - phi_f;
	Mc_h = total_M * phi_c0;
	Mm_h = total_M * phi_m0;
	
	// print to screen
	cout << scientific;
	cout << "Mc_h=" << Mc_h << ", total_M=" << total_M <<endl;  

	//VF:initialization of the first component 
	*pk1_p = 1.0;
	*pk2_p = 1.0;
	*pk3_p = 1.0;
	*pm_p = 1.0;
	*r_p = r_h;
	                                                                           
	for (i = 1; i < num_pa; i++)
	{
		*(pk1_p+i) = *(pk1_p+i-1) * exp(-0.5 * dt * (mu_F(i*dt,0) + mu_F((i-1)*dt,0))); 
		*(pk2_p+i) = *(pk2_p+i-1) * exp(-0.5 * dt * (mu_F(i*dt,0) + mu_F((i-1)*dt,0)));         
		*(pk3_p+i) = *(pk3_p+i-1) * exp(-0.5 * dt * (mu_F(i*dt,0) + mu_F((i-1)*dt,0)));         
		*(pm_p+i) = *(pm_p+i-1) * exp(-0.5 * dt * (mu_F(i*dt,0) + mu_F((i-1)*dt,0)));           
		
		// finding a half-life time
		if( flag && (*(pk1_p + i)> 0.5))
		{
			dummy1 = (i-1)*dt;
			dummy2 = i*dt;
			dummy3 = *(pk1_p+i-1);
			dummy4 = *(pk1_p+i);
			if(dummy4!=dummy3)	a_half = dummy1 + dt * (0.5 - dummy3)/(dummy4 - dummy3); //VF: what value is for "else"?
			flag = 0;
		}
		a_mean += 0.5 * dt * (*(pk2_p+i-1) + *(pk2_p+i)); 
		*(r_p+i) = r_h;
	}
  
	//print to screen
	cout << scientific;
	cout <<" "<< endl;
	cout << "half-life time =" << static_cast< float >(a_half) << endl;
	cout << "mean life time =" << static_cast< float >(a_mean) << endl;
	  
	beta_2 = k_2 /(k_2 + mu_2);
	Mc = 0.0;
	for(i=0;i<4;i++)
	{
		M[i]  = Mc_h * phi_k0[i]; 		
		Mc += M[i];
		//mP[i] = M[i]/(MWH * a_mean * beta_1 * beta_2); 
		// collagen mass production at t=0 is MWH * a_mean * beta_1 * beta_2
		mPh[i] = M[i]/(MWH * a_mean * beta_1 * beta_2); //VF, CI(t=0)=a_mean * beta_1 
		mP_p[i]= mPh[i]; //VF
		mP_h[i]= mPh[i]; //VF
		//CI[i] = M[i]/(MWH * k_2 * a_mean); 
		CIh[i] = M[i]/(MWH * k_2 * a_mean); //VF
		// concentration of procollagen at t=0
		CI_p[i] = CIh[i]; //VF
		sigma_k_p[i] = stress_h; //VF: sigma_k_p won't be used
	}
				
	Me = total_M*phi_e0;
	Mm = total_M*phi_m0;
	CIm = Mm / (MWH*k_2*a_mean);
	CIm_p = CIm;
	mPm = Mm /(MWH*a_mean*beta_1*beta_2);
	mPm_p = mPm;
	mPm_h = Mm / a_mean;   //added on october 14, 2015

	// VF: alpha_R won't be used 
	alpha_R[0] = 0.0;
	alpha_R[1] = PI/2.0;
	alpha_R[2] = alpha_h * PI/180.0;
	alpha_R[3] = 2.0 * PI - alpha_h * PI/180.0;

	for(i=0;i<num_pa;i++)
	{
		*(pk1_p+i) *= M[0]/a_mean;
		*(pk2_p+i) *= M[1]/a_mean;
		*(pk3_p+i) *= M[2]/a_mean;
		*(pm_p+i) *= Mm/a_mean;
	}
  
	for(i=0;i<3;i++)
	{
		T1_c[i] = stress_h; 
		T2_c[i] = stress_h; 
		T2_m[i] = stress_h;
	}
				
}	
	
//**********************************************************************
 double  PointGR::mu_F(const double age, const double zeta)
//**********************************************************************
// Compute rate of degradation
// WARNING: age variable is not used!
// Parameters:
//		In (const): age,zeta
	{
		double rate;

		if (zeta<zeta_c) rate = 1.0/100.0;
		else			 rate = 1.0/100.0 + pow(zeta-zeta_c,2.0) * Kg_mu/100.0;
		
		return  rate;	
	}	