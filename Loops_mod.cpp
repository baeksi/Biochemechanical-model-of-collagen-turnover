// Functions for Growth and Remodeling Loops part
// by Vasilina, 12-06-2017 

#include	"GrowthRemodeling_mod.h"

//****************************************************************************************************************************************
void PointGR::FirstLoop_SolveGrowthRemodeling(ofstream & outStreamT,const double dt, const double km2,const double r_h,const double ri_h,
				int & final_time)
//****************************************************************************************************************************************
// It is the First loop in the program. It contains Second Loop with Third Loop nested inside.
// It is main time-loop to solve ODE, in days
//
// Parameters: 	6=5+1
//		In (const):outStreamT,dt,km2,r_h,ri_h,
//		InOut (update): final_time
//
// Class vars used as const: Mc_h,Mm_h,mPm_h,Me,mP_h,
// Class vars updated:	CIm_p,mPm_p,Mc,Mm,total_M,final_time,
//					    CI,CI_p,mP,mP_p,M,sigma_k,T1_c,T2_c,T2_m,pk1_p,pk2_p,pk3_p,pm_p,r_p
// Use Functions:
// 		GetLZetas
//		GetMandPointersP
//		SecondLoop_SolveMassProductionODE
//		OutputResultsEachTime
//		Given_P
//		Given_flow
{				
	//internal variables
	double 	dummy1,dummy2,dummy3,dummy4,zeta0,zeta1,zeta2,zetam;
	double	CIm,mPm,sum,t;
	double  sigma_k_p[4]={}; //VF: sigma_k_p won't be used
	int		i;
	
	//VF: update class variable
	r_act = r_h;

	//VF: initialization of internal vars for a loop and for update in other functions
	mPm = 0.0;		sum = 0.0;		t = 0.0;
	
	// ---------------- start first loop ------------------------------
	while(t <= final_time)
	{
		t += dt; 

		// initialization inside a loop
		pressure = Given_P(t);
		zeta0=0; zeta1=0;	zeta2=0; zetam=0;
		
		dummy1 = fabs( (Given_P(t)-Given_P(t-dt))/Given_P(t) );                                  
		dummy2 = fabs( (Given_flow(t)-Given_flow(t-dt))/Given_flow(t) );                         
		dummy3 = fabs( (Given_P(t-dt)-Given_P(t-2.0*dt))/Given_P(t-dt) );                         
		dummy4 = fabs( (Given_flow(t-dt)-Given_flow(t-2.0*dt))/Given_flow(t-dt) );  
	
		// Get r and T1_c[0],T2_c[0],T2_m[0]
		/* if the change in the flow is significant then predictor is set to the 
		previous value, otherwise use 1st order predictor */
		if(dummy1 + dummy2 + dummy3 + dummy4 > 0.01)
		{                      
			r = *r_p;
			T1_c[0] = T1_c[1];
			T2_c[0] = T2_c[1];
			T2_m[0] = T2_m[1];
		}else
		{
			r = 0.5*(3.0* (*r_p) - *(r_p+1));
			T1_c[0]=0.5*(3.0*T1_c[1]-T1_c[2]);
			T2_c[0]=0.5*(3.0*T2_c[1]-T2_c[2]);
			T2_m[0]=0.5*(3.0*T2_m[1]-T2_m[2]);
		}      
		// Get L2, r_act
		L2 = r/r_h;	 																																								
		sum	+= K_act*L2*dt*exp(K_act*t);  // added for r_act, Lm_act, //VF: won't be used   	 	
		r_act += K_act*(r-r_act)*dt;
					
		// Get alpha[]
		alpha[0] = 0.0;
		alpha[1] = PI/2.0;
		alpha[2] = atan(L2/L1*tan(alpha_h*PI/180));																					 
		alpha[3] = 2.0*PI- alpha[2];
		
		//VF: Compute Lk2_n,Lm2_n2 and zeta0,zeta1,zeta2,zetam,
		//			Use: L2,r,alpha,pk1_a,pk2_a,pk3_a,pm_a,r_p,M
		//			Update:Lk2_n,Lm2_n2
		GetLZetas(dt,r_h,
			zeta0,zeta1,zeta2,zetam);
			
		//VF: moved out from the next loop to have shear defined in general
		if (t<=dt) shear = shear_h;
		
		//VF: Get CI[],mP[]
		for(i=0;i<4;i++)
		{
			sigma_k[i] = sqrt(pow(T1_c[0]*cos(alpha[i]),2.0)+pow(T2_c[0]*sin(alpha[i]),2.0));
			// mass production of collagen fiber
			mP[i] = f_mP(Mc_h, Mc, sigma_k[i], mP_h[i], Kg,shear);
			// concentration of collagen fiber
			CI[i] = CI_p[i] + 0.5*dt*((1-(k_2+mu_2)*dt)*beta_1*mP_p[i]
					+ (-2.0*(k_2+mu_2)+dt*(k_2+mu_2)*(k_2+mu_2))*CI_p[i] + beta_1*mP[i]);
		}

		//added after discussion with DR. Baek, it is evolution of reference coordinate
		if (MWH1==0)	MWH1 = MWH;				//VF: strange...
		else			MWH1 = L2*L1*h*MWH/h_h; //VF: note that h will be updated later, thus it was initialized as h=0 before a loop
	
		//VF: update mPm,Mc,Mm,total_M,M,pk1_a,pk2_a,pk3_a,pm_a
		//		Use: L2,mPm_h,Me,MWH1,Mc_h,Mm_h,shear,mP_h,T2_m,pk1_p,pk2_p,pk3_p,pm_p,CI		 
		//		Update:Mc,Mm,total_M,M,pk1_a,pk2_a,pk3_a,pm_a
		GetMandPointersP(dt,zeta0,zeta1,zeta2,zetam,
						mPm);

		//VF: Update thickness
		h = total_M / (rho*L1*L2);

		//Get CIm
		CIm = CIm_p + 0.5*dt*((1-(k_2+mu_2)*dt)*beta_1*mPm_p + (-2.0*(k_2+mu_2) + dt*(k_2+mu_2)*(k_2+mu_2))*CIm_p + beta_1*mPm);
 
		////---------- call second loop, which also contains third loop-----------------------------------
		//Use:    CIm_p,mPm_h,mPm_p,Mc_h,Me,Mm_h,MWH1,pressure, CI_p,mP_h,mP_p,pm_p,pk1_p,pk2_p,pk3_p,r_p,
		//Update: h,hm,L2,Lk2_n,L2m_act,Lm2_n2,Mc,Mm,r,r_act,
		//		  S,shear,t_act,total_M,pk1_a,pk2_a,pk3_a,pm_a,alpha,M,mP,sigma_k,T1_c,T2_c,T2_m
		SecondLoop_SolveMassProductionODE(outStreamT,dt,km2,r_h,ri_h,t,
										 CIm,zeta0,zeta1,zeta2,zetam,final_time);
		//------------------------------------------------------------------------------------------------
	
		for(i=0; i<4; i++)
		{
			mP_p[i] = mP[i];
			CI_p[i] = CI[i];		//VF: update CI_p
			sigma_k_p[i]=sigma_k[i]; //VF: won't be used
		}

		mPm_p = mPm;
		CIm_p = CIm; //VF: CIm is updated in SecondLoop
		for (i=0;i<num_pa;i++)
		{
			*(pk1_p+i) = *(pk1_a+i);
			*(pk2_p+i) = *(pk2_a+i);
			*(pk3_p+i) = *(pk3_a+i);
			*(pm_p+i)  = *(pm_a+i);

			//VF: update *r_p, from homeostatic it is *(r_p+i)=r_h
			if(i==num_pa-1)	*r_p=r;
			else			*(r_p+num_pa-i-1) = *(r_p+num_pa-i-2);
		}

		for (i=2;i>0;i--)
		{
			T1_c[i]=T1_c[i-1];
			T2_c[i]=T2_c[i-1];
			T2_m[i]=T2_m[i-1];
		}
		
		// output current results into file
		// use: r,h,shear,pressure,t_act,hm,S,L2m_act
		OutputResultsEachTime(outStreamT,t,zeta1,T1_c[0],T2_c[0],T2_m[0],Mm);
			
	} // ----------------------end first loop-------------------------------------

}

//**********************************************************************************************************************************************************
void PointGR::SecondLoop_SolveMassProductionODE(ofstream & outStreamT, const double dt,const double km2,const double r_h,const double ri_h, const double t,
					double & CIm,double & zeta0,double & zeta1,double & zeta2,double & zetam,int & final_time)
//**********************************************************************************************************************************************************
// Second loop: compute mass production.
// Parameters: 12=6+6 
//		In (const):outStreamT,dt,km2,r_h,ri_h,t
//		InOut (updtae): CIm,zeta0,zeta1,zeta2,zetam,final_time
// Internal class parameters: 
//		In: CIm_p,mPm_h,mPm_p,Mc_h,
//			 Me,Mm_h,MWH1,pressure, CI_p,mP_h,mP_p,pm_p,pk1_p,pk2_p,pk3_p,r_p,	  
//		Update:	h,hm,L2,Lk2_n,L2m_act,Lm2_n2,Mc,Mm,r,r_act,
//				S,shear,t_act,total_M,pk1_a,pk2_a,pk3_a,pm_a,alpha,CI,M,mP,sigma_k,T1_c,T2_c,T2_m;
// Use Funcstions:
//		ThirdLoop_SolveStressEquilibrium
//		GetLZetas
//		GetMandPointersP

{
	// internal variables
	double	dwdL22_m,dwdL22_c,dwdL12_c,error1,hc,mPm,T1_c_old,T2_c_old,T2_m_old;
	int 	i,num_it1,num_it2;

	// initialization
	error1 = 10.0; 		num_it1 = 0;		num_it2 = 0;	

	//--------------- call second loop, which also contains third loop---
	while(error1 > Tol && (num_it1 < Max_it && num_it2 < Max_it))
	{ 
		//VF: internal initialization for third loop or for the use in functions
		dwdL22_m = 0.0;			dwdL22_c = 0.0;			dwdL12_c = 0.0;
		hc = 0.0;				mPm=0.0;
		T1_c_old = T1_c[0];		T2_c_old = T2_c[0];		T2_m_old = T2_m[0];
		//VF: update to zero again:
		num_it2 = 0; zeta0=0; zeta1=0; zeta2=0; zetam=0;
			
		num_it1++;
		if(num_it1==Max_it)
		{
			outStreamT << "Iteration reachs maximum number !!\n";	
			final_time=0;
		}

		//------------------------third while loop---------------------------------
		//In:	outStreamT,dt,km2,Mc,Me,Mm,Mm_h,pressure,r_h,ri_h,t,total_M,num_pa,
		//		pk1_a,pk2_a,pk3_a,pm_a,r_p,
		//Update: h,hc,hm,dwdL22_m,dwdL22_c,dwdL12_c,L2,Lk2_n,Lm2_n2,L2m_act,
		//		  r,r_act,S,shear,t_act,final_time,num_it2,alpha
		ThirdLoop_SolveStressEquilibrium(outStreamT,dt,km2,r_h,ri_h,t,
								hc,dwdL22_m,dwdL22_c,dwdL12_c,final_time,num_it2);
		//-------------------------------------------------------------------------
	
		// update thickness of SMC constituent
		hm = Mm /((1-phi_f)*rho*L1*L2);
		// update thickness of collagen constituent
		hc = Mc / ((1-phi_f)*rho*L1*L2);
		// update axial and hoop Cauchy memebrane stress (passive) of collagen
		T1_c[0] = 2*L1*dwdL12_c/(L2*hc);
		T2_c[0] = 2*L2*dwdL22_c/(L1*hc);
		// update hoop Cauchy memebrane stress of SMC, with active component
		T2_m[0] = (2*L2*dwdL22_m/L1+t_act)/hm;                                                                                                
	
		// update radial stretch 
		L2=r/r_h;	
	
		//VF: update alpha
		alpha[2] = atan(L2/L1*tan(alpha_h*PI/180));                                                                                    
		alpha[3] = 2.0*PI-alpha[2];
	
		//VF: compute Lk2_n,Lm2_n2 and zeta0,zeta1,zeta2,zetam,
		//			Use: L2,r,alpha,pk1_a,pk2_a,pk3_a,pm_a,r_p,M
		//			Update:Lk2_n,Lm2_n2
		GetLZetas(dt,r_h,
				 zeta0,zeta1,zeta2,zetam);
	
		//Get mP[],sigma_k[],update CI[]
		for(i=0;i<4;i++)
		{
			sigma_k[i] = sqrt(pow(T1_c[0]*cos(alpha[i]),2.0)+pow(T2_c[0]*sin(alpha[i]),2.0));
			mP[i]= f_mP(Mc_h, Mc, sigma_k[i], mP_h[i], Kg,shear);
			CI[i]= (CI_p[i]+0.5*dt*(beta_1*mP_p[i]-(k_2+mu_2)*CI_p[i]+ beta_1*mP[i]))/(1.0+0.5*dt*(k_2+mu_2)); //VF: differs from 1st loop
		}
	
		//VF: update mPm,Mc,Mm,total_M,M,pk1_a,pk2_a,pk3_a,pm_a
		//		Use: L2,mPm_h,Me,MWH1,Mc_h,Mm_h,shear,mP_h,T2_m,pk1_p,pk2_p,pk3_p,pm_p,CI		 
		//		Update:Mc,Mm,total_M,M,pk1_a,pk2_a,pk3_a,pm_a
		GetMandPointersP(dt,zeta0,zeta1,zeta2,zetam,
					    mPm);

		CIm = (CIm_p+0.5*dt*(beta_1*mPm_p-(k_2+mu_2)*CIm_p+beta_1*mPm))/(1.0+0.5*dt*(k_2+mu_2)); //VF: differs from 1st loop
	
		//VF: Update thickness
		h = total_M / (rho*L1*L2);

		if((T1_c_old!=0)&&(T2_c_old != 0) && (T2_m_old != 0))
		{
			error1 = fabs((T1_c_old-T1_c[0])/T1_c_old)
					+fabs((T2_c_old-T2_c[0])/T2_c_old)+fabs((T2_m_old-T2_m[0])/T2_m_old);	
		}

	} //------------------end of second loop----------------------------------------------
}

//*******************************************************************************************************************************************************
void PointGR::ThirdLoop_SolveStressEquilibrium(ofstream & outStreamT,const double dt,const double km2,const double r_h,const double ri_h,const double t,
					double & hc,double & dwdL22_m,double & dwdL22_c,double & dwdL12_c,int & final_time,int & num_it2)
//*******************************************************************************************************************************************************
// Third loop: Newton-Raphson iteration loop to determine the mean radius at each time step from stress equilibrium
// Reference: [Baek,Rajagopal,Humphrey(2006)], Eq.35
//
// Parameters: 12=6+6 
//		In(const):outStreamT,dt,km2,r_h,ri_h,t,
//		InOut(update): hc,dwdL22_m,dwdL22_c,dwdL12_c,final_time,num_it2;
// Internal class parameters: 
//		In: Mc,Me,Mm,Mm_h,pressure,
//			total_M,pk1_a,pk2_a,pk3_a,pm_a,r_p,
//		Update:h,hm,L2,Lk2_n,Lm2_n2,L2m_act,
//			r,r_act,S,shear,t_act,alpha;
{
	// internal variables
	double	alpha_a,ddwddL22,ddwddL22_c,ddwddL22_m,dwdL22,dfun_rdr,dsdr,dt_actdr,error2,fun_r;
	double	L2_a,Le2_n1,Le2_n2,Lk2_a,Lk2,Lm2_n1,pk_a,r_old,ri,sigma_m,wt;
	int		a,i;
	
	// initialization
	error2 = 10.0;

	while(error2 > Tol && num_it2<Max_it)
	{
		//VF: Initialization
		pk_a = 0.0;  dt_actdr =0.0; 

		num_it2 = num_it2+1;                                                                
		if(num_it2==Max_it)
		{
			outStreamT << "Iteration for N-R method reachs maximum number!!\n";	
			final_time = 0;
		}

		r_old = r;                                                                                               
		L2 = r/r_h;        
		h = total_M /(rho*L1*L2);																											
		ri = r - 0.5*h;                                                                                                                  
		shear = shear_h*Given_flow(t)*pow(ri_h/ri, 3.0);                                                                    

		// VF: intialize
		alpha[0] = 0.0;
		alpha[1] = PI/2.0;
		alpha[2] = atan(L2/L1*tan(alpha_h*PI/180));                                 
		alpha[3] = 2.0*PI- alpha[2];

		dwdL22_c = 0.0; 
		dwdL12_c = 0.0;
		dwdL22_m = 0.0;
		ddwddL22_c = 0.0;
		ddwddL22_m=0.0;
	      	      
		Le2_n1 = Ge_1*Ge_1*L1*L1;
		Le2_n2 = Ge_2*Ge_2*L2*L2;
		dwdL22 = Me*(ce/2.0)*Ge_2*Ge_2*(1.0-1.0/(Le2_n1*Le2_n2*Le2_n2));
		ddwddL22 = Me*ce*pow(Ge_2,4.0)/(Le2_n1*pow(Le2_n2,3.0));
	      
		for(a=0;a<num_pa;a++)
		{                     
			// get active stretch L2_a
			if (a==0)
			{	L2_a = r/r_h;
				wt = 0.5*dt;
			}else
			{   L2_a = *(r_p+a-1)/r_h;
				if (a==num_pa-1)	wt = 0.5*dt;
				else 				wt = dt;
			}   

			// for each k-th collagen-fiber family
			for(i=0; i<4; i++)
			{																					
				if (i==0 || i==1) 	alpha_a=alpha[i];
				if (i==2)			alpha_a=atan(L2_a/L1*tan(alpha_h*PI/180));
				if (i==3)			alpha_a=2.0*PI-atan(L2_a/L1*tan(alpha_h*PI/180));
				
				//VF: mapped stretches of SMC and of collagne fibers
				Lk2_a = pow(L1*cos(alpha_a),2.0)+pow(L2_a*sin(alpha_a),2.0);
				Lk2 = pow(L1*cos(alpha[i]),2.0)+pow(L2*sin(alpha[i]),2.0);
					
				if 	(Lk2_a>Tol)		Lk2_n = Gc_h*Gc_h*Lk2/Lk2_a;																		  
				if	(Lk2_n<1) 		Lk2_n=1;
				if	(i==0) 			pk_a = *(pk1_a+a);
				if	(i==1) 			pk_a = *(pk2_a+a);
				if	(i==2 || i==3) 	pk_a = *(pk3_a+a);			    										
					
				dwdL22_c += 0.5*cc1*pk_a*Gc_h*Gc_h/Lk2_a*pow(sin(alpha_a),2.0)
						*(Lk2_n-1)*exp(cc2*(Lk2_n-1)*(Lk2_n-1))*wt;
				dwdL12_c += 0.5*cc1*pk_a*Gc_h*Gc_h/Lk2_a*pow(cos(alpha_a),2.0)
						*(Lk2_n-1)*exp(cc2*(Lk2_n-1)*(Lk2_n-1))*wt;
				ddwddL22_c += 0.5*cc1* pk_a *Gc_h*Gc_h*Gc_h*Gc_h/(Lk2_a*Lk2_a)
						*pow(sin(alpha_a),4.0)*(1.0+2.0*cc2*(Lk2_n-1)
						*(Lk2_n-1))*exp(cc2*(Lk2_n-1)*(Lk2_n-1))*wt;
			}
	          
			Lm2_n1 = Gm_h *Gm_h;
			Lm2_n2 = Gm_h * Gm_h*L2*L2 / (L2_a*L2_a);
				
			if	(Lm2_n2<1) Lm2_n2=1;								     	
				
			dwdL22_m += *(pm_a+a)* Gm_h*Gm_h/(L2_a*L2_a)
					*(0.5*cm1*(1-1.0/(Lm2_n1*Lm2_n2*Lm2_n2))
					+0.5*cm2*(Lm2_n2-1)*exp(cm3*(Lm2_n2-1)*(Lm2_n2-1)))*wt;
			ddwddL22_m += *(pm_a+a)* pow(Gm_h/L2_a, 4.0)
					*(cm1/(Lm2_n1*pow(Lm2_n2,3.0))
					+0.5*cm2*(1+2*cm3*(Lm2_n2-1)*(Lm2_n2-1))
					*exp(cm3*(Lm2_n2-1)*(Lm2_n2-1)))*wt;                                                                                                                                                                               
		}

		/* ---------- active muscle tone---------- */
		// update total thickness
		h = total_M/rho;            
		// get thickness of smooth muscle constituent
		hm = Mm/((1-phi_f)*rho*L1*L2);                                                                               
		sigma_m = 2*L2*dwdL22_m/(L1*hm);// not used     
						
		// get S-function for t_act (it differes from the paper?)
		//S is changed to hyperbolic tanged to avoid  if condition 																								
		S=(1/(1-exp(-pow(km2/20,2.0))))*(Mm/Mm_h)*S_basal*(1-exp(-pow((km2/20-km2*(shear/shear_h-1)),2.0)));     																				  																										 
		
		// update radial active stretch of SMC
		L2m_act=r/r_act;  	

		dsdr=6*(1/(1-exp(-pow(km2/20,2.0))))*(Mm/Mm_h)*S_basal*(km2/20-km2*(shear/shear_h-1))*(km2*shear/(ri*shear_h))	
				* (1+0.5*total_M/(rho*L1*L2*L2*r_h))*exp(-pow(km2/20-km2*(shear/shear_h-1),2.0));   																																												

		if (pow((Lambda_M-L2m_act)/(Lambda_M-Lambda_0),2.0)<=1)
		{																				
			// compute membrane normal stress due to vascular smooth muscle tone - active stress of SMC
			t_act = S*L2m_act*(1-pow((Lambda_M-L2m_act)/(Lambda_M-Lambda_0),2.0));			         			   		    
			if (r!=0)	
			{	
				dt_actdr = dsdr*L2m_act*(1-pow((Lambda_M-L2m_act)/(Lambda_M-Lambda_0),2.0))
						 + t_act/r
						 + 2*L2m_act*S*(Lambda_M-L2m_act)*L2m_act/(r*pow(Lambda_M-Lambda_0,2.0));	
			}							 
		}
		else
		{
			t_act = 0;	
			dt_actdr = 0;
		}																																																																											
		dwdL22 += dwdL22_m + dwdL22_c;
		ddwddL22 += ddwddL22_m + ddwddL22_c;                       
		hc = Mc / ((1-phi_f)*rho*L1*L2);
		fun_r = 2.0*L2*dwdL22/L1 + t_act - pressure*r;                       																							
		dfun_rdr = (1.0/r_h)*(2*dwdL22/L1 + 4*L2*L2*ddwddL22/L1) + dt_actdr - pressure;                    

		if (fabs(dfun_rdr)<Tol)
		{
			if (dfun_rdr<0) dfun_rdr=-Tol; 
			else dfun_rdr = Tol;
		}

		//VF: update r
		r -= 0.2*fun_r/dfun_rdr;  	
		r_act = r_act + K_act*(r-r_old); 

		if (r_old !=0)  error2 = fabs((r_old-r)/r_old);

	}// ------------- end of third loop--------------
}

//**********************************************************************************
void PointGR::GetLZetas(const double dt,const double r_h,	      
					double & zeta0,double & zeta1,double & zeta2,double & zetam)
//**********************************************************************************
// Compute all zetas needed for the First and Second loops
// Get relative tension of collagen fibers and SMC
//		-4th fiber has the same relative tension as 3rd one
// 	Parameters: 7 = 3 + 4 
//		In (const):	dt,r_h,M,
//		InOut(update):	zeta0,zeta1,zeta2,zetam;
// Internal class parameters: 
//		In:	L2,r,alpha,pk1_a,pk2_a,pk3_a,pm_a,r_p,
//		Update:	Lk2_n,Lm2_n2;
//
	{
		// internal variables
		double 	alpha_a,L2_a,Lk2,Lk2_a,Lm2_n1,wt;
		int 	a,j;
		
		for(a=0;a<num_pa;a++)
		{
			if (a==0)
			{	L2_a = r/r_h;
				wt = 0.5*dt;
			}
			else
			{
				L2_a = *(r_p+a-1)/r_h;
				if(a==num_pa-1) wt = 0.5*dt;
				else wt = dt;
			}      

			for (j=0;j<3;j++) // VF: 4th fiber has the same relative tension as 3rd one
			{
				if (j==0 || j==1) alpha_a=alpha[j];  // angle
				if (j==2)	alpha_a=atan(L2_a/L1*tan(alpha_h*PI/180));
					
				Lk2_a = pow(L1*cos(alpha_a),2.0)+pow(L2_a*sin(alpha_a),2.0);
				Lk2 = pow(L1*cos(alpha[j]),2.0)+pow(L2*sin(alpha[j]),2.0);
					
				if (Lk2_a>Tol) Lk2_n = Gc_h*Gc_h*Lk2/Lk2_a;   
				if(Lk2_n<1) Lk2_n=1; 
				if (j==0)   zeta0 +=(*(pk1_a+a))*(sqrt(Lk2_n))*(Lk2_n-1)*exp(cc2*pow(Lk2_n-1,2.0))*wt;
				if (j==1)	zeta1 +=(*(pk2_a+a))*(sqrt(Lk2_n))*(Lk2_n-1)*exp(cc2*pow(Lk2_n-1,2.0))*wt;
				if (j==2)	zeta2 +=(*(pk3_a+a))*(sqrt(Lk2_n))*(Lk2_n-1)*exp(cc2*pow(Lk2_n-1,2.0))*wt;											         
			}	
		
			Lm2_n1 = Gm_h *Gm_h; // VF: the same as in Third loop
				
			if(L2_a>Tol)	Lm2_n2 = Gm_h * Gm_h*L2*L2 / (L2_a*L2_a); 
			if(Lm2_n2<1)	Lm2_n2=1;								  
			
			zetam += *(pm_a+a)*Gm_h*(L2/L2_a)*(cm1*(1-1/(Lm2_n1*pow(Lm2_n2,2.0)))+cm2*(Lm2_n2-1)*exp(cm3*pow(Lm2_n2-1,2.0)))*wt;																																					    		                                                                                                                                                                            
		}
											
		zeta0 = zeta0/(M[0]*Gc_h*(Gc_h*Gc_h-1)*exp(cc2*pow(Gc_h*Gc_h-1,2.0)));
		zeta1 = zeta1/(M[1]*Gc_h*(Gc_h*Gc_h-1)*exp(cc2*pow(Gc_h*Gc_h-1,2.0)));
		zeta2 = zeta2/(M[2]*Gc_h*(Gc_h*Gc_h-1)*exp(cc2*pow(Gc_h*Gc_h-1,2.0)));																							
		zetam = 1;    //added on March 31,2015, to make smooth muscle cells degradation rate constant
		
	}

//**************************************************************************************************************************
void PointGR::GetMandPointersP(const double dt,const double zeta0,const double zeta1,const double zeta2,const double zetam,
					double & mPm)
//**************************************************************************************************************************
//Compute some variables needed for First and Second loops
// Parameters: 6=5+1
// 		In (const): dt,zeta0,zeta1,zeta2,zetam
// 		InOut (update):mPm;
// Internal class parameters: 
// 		In: L2,mPm_h,Me,MWH1,Mc_h,Mm_h,
//			shear,mP_h,T2_m,pk1_p,pk2_p,pk3_p,pm_p,
// 		Update:Mc,Mm,total_M,M,pk1_a,pk2_a,pk3_a,pm_a;
// Use Functions:
//		f_mP
{
	// internal variables
	int i;
																							
	*pk1_a = MWH1*k_2*CI[0];
	*pk2_a = MWH1*k_2*CI[1];
	*pk3_a = MWH1*k_2*CI[2];
	mPm = f_mP(Mm_h, Mm, T2_m[0], mPm_h, Kg,shear);
	*pm_a = mPm;  //VF: in the first loop there was "comment oct 14, 2015: *pm_a = MWH1*k_2*CIm;"
	
	//VF:why M, Mm are updated again?
	M[0]=0;
	M[1]=0;
	M[2]=0;
	Mm = 0;          
	for(i=1;i<num_pa;i++)
	{	    
		*(pk1_a+i) = *(pk1_p+i-1)*exp(-0.5*dt*(mu_F(i*dt,zeta0)+mu_F((i-1)*dt, zeta0)));
		*(pk2_a+i) = *(pk2_p+i-1)*exp(-0.5*dt*(mu_F(i*dt,zeta1)+mu_F((i-1)*dt, zeta1)));
		*(pk3_a+i) = *(pk3_p+i-1)*exp(-0.5*dt*(mu_F(i*dt,zeta2)+mu_F((i-1)*dt, zeta2)));
		*(pm_a+i) = *(pm_p+i-1)*exp(-0.5*dt*(mu_F(i*dt,zetam)+mu_F((i-1)*dt, zetam)));
		// mass production of collagen fiber 1,2,3
		M[0] += 0.5*dt*(*(pk1_a+i)+*(pk1_a+i-1));                                                                 
		M[1] += 0.5*dt*(*(pk2_a+i)+*(pk2_a+i-1));                                                                 
		M[2] += 0.5*dt*(*(pk3_a+i)+*(pk3_a+i-1));      
		// get total mass of SMC
		Mm += 0.5*dt*(*(pm_a+i)+*(pm_a+i-1));                                                                   
	}
	M[3] =	M[2];
	Mc = M[0]+M[1]+M[2]+M[3];                                                                                   	  
	total_M = (Mc+Mm+Me)/(1-phi_f);

	//delete[] H;
}

//***************************************************
 double  PointGR::Given_P(const double time)
//***************************************************
// Intramural Pressure as a function of time (day)
// Parameters: 1
//		In (const): time
{
	if (time>=1)	return 1.0*P_h;
	else			return P_h;
}

//******************************************************
double  PointGR::Given_flow(const double time)
//******************************************************
// relative increase in flow; 
// flow(t)=1.0 at the homeostatic state
// Parameters: 1
//		In (const): time
{
	if (time>=1)	return 0.8;   
	else			return 1.0;
}
	
//*********************************************************************************************************************************************
double PointGR::f_mP(const double Mi_0, const double Mi_t,const double Sigma, const double mP_b, const double kkgg, const double tau)
//*********************************************************************************************************************************************
// Mass production of the constituent
// Biochemomechanical model, equation (17) in [Getachew,Humphrey,Baek]	
// Compute mP - synthesis rate of intracellular procollagen and proliferation of SMC
// It is given by a stress-dependent scalar function of time.
//
// Mi_0,Mi_t are not used
// Parameters: 4
// 		In (const):Mi_0,Mi_t,Sigma,mP_b,kkgg,tau;
{
	double mP;

	// mass procollagen production is taken to be 1.0 here
	mP = mP_b*(1.0 + kkgg*(Sigma/stress_h -1.0) - Kg_sh*(tau/shear_h-1));
	
	if (mP>0)	return mP;
	else		return 0.0;
}

//*****************************************************************************************************************************************************************************
void PointGR::OutputResultsEachTime(ofstream & outStreamT,const double t,const double zeta1,const double T1_c0,const double T2_c0,const double T2_m0,const double Mm) 
//*****************************************************************************************************************************************************************************
//  Print results in the end of first loop
//	Parameters: 8
//		In (const):outStreamT,t,zeta1,T1_c0,T2_c0,T2_m0,Mm;
// Internal class parameters: 
//		In:r,h,shear,pressure,t_act,hm,S,L2m_act;
	{			
		// print into screen
		cout << scientific;
		cout << "time=" <<fixed << setprecision(2)<< t 
			<< scientific << setprecision(6) << ", r=" << r << '\t' << zeta1 << endl;

		// print current results into file
		outStreamT	<< fixed << setprecision(3)<< t <<'\t'<< scientific << setprecision(6) 
					<< r <<'\t'<< h <<'\t'<<T1_c0<<'\t'<<T2_c0<<'\t'<<T2_m0<<'\t'<< shear<<'\t'<<pressure*r/h<<'\t'<<t_act/hm<<'\t'<<S<<'\t'<<L2m_act<<'\t'<<Mm<<'\t' 
					<< fixed << setprecision(3) << endl;
	}