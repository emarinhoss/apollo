/*
 * CR_BKE.cpp
 *
 *  Created on: Jan 6, 2016
 *      Author: richard
 */

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <string.h>
#include <omp.h>
#include <gsl/gsl_linalg.h>

#include "cArray.h"
#include "CRSolve.h"
#include "MathX.h"

using namespace std;

// Constructor
CR_BKE::CR_BKE(int Ns, bool intNe, bool fxdNe, bool Eh, bool Ee)
{
	nq		= 0;
	filenum = 0;
	Natoms 	= 0;
	Ngroups = 0;
	NTstep	= 0;
	Nexlvls = 0;

	Nelec	= 0.;
	TeV		= 0.;
	TrV		= 0.;
    Eexce	= Edexe	= Eione	= Erece	= 0.;
    Eexcp	= Eionp	= 0.;
    Eexcs	= Eions	= 0.;
    Edexp	= Erecp	= 0.;
    Efree	= 0.;
    Etot	= 0.;
    Eint	= 0.;

	Zq		= NULL;
	pRk 	= NULL;
//	pSpect	= NULL;

	Rad 	= false;
	Irrad	= false;
	Nstate 	= Ns;
	DensNe 	= intNe;
	EnerEh 	= Eh;
	EnerEe 	= Ee;
	FxdNe	= fxdNe;
	if (FxdNe)	DensNe = false;		// Conflict of bools!
	SpecOut = false;

    NCRdim = Ns;
    if (DensNe)	{nev = NCRdim; NCRdim++;}
    if (EnerEh)	{neh = NCRdim; NCRdim++;}
    if (EnerEe)	{nee = NCRdim; NCRdim++;}

    state_QSS 	= NULL;
    RadSEZ		= NULL;
    RadSEbb		= NULL;
    NcFull		= NULL;
    Nstar 		= new double [NCRdim];
    Ncand 		= new double [NCRdim];
    Nc 			= new double [NCRdim];
    Scr 		= new double [NCRdim];
    Jcr 		= new double [NCRdim];
    dNc 		= new double [NCRdim];
    Kcr 		= new double*[NCRdim]; for(int n=0;n<NCRdim;n++)  Kcr[n] = new double[NCRdim];
    Jac 		= new double*[NCRdim]; for(int n=0;n<NCRdim;n++)  Jac[n] = new double[NCRdim];

    printlog= false;
    tcurr 	= 0.;
    tbeg 	= 0.;
    tend	= 0.;
    dt		= 0.;
    multxdt = 1.;
    tprint	= 0.;
    dtprint = 0.;
    logtprt = 0.;
    logdtprt= 0.;
    rel_tol = 0.;
    abs_tol = 0.;

    fp		= NULL;
    fo		= NULL;
    fTg		= NULL;
    fRad	= NULL;
    fChrg	= NULL;
    fp_full	= NULL;
    ftBolt	= NULL;

}

// Destructor
CR_BKE::~CR_BKE()
{
	delete[] Zq;
	delete[] Nc;
	delete[] Nstar;
	delete[] Ncand;
    delete[] Scr;
    delete[] dNc;
    if (Jcr) delete[] Jcr;
    if (RadSEZ) delete[] RadSEZ;
    if (RadSEbb) delete[] RadSEbb;

    for(int n=0;n<NCRdim;n++)
    {
    	if (Kcr) delete[] Kcr[n];
    	if (Jac[n]) delete[] Jac[n];
    }
    if (Kcr) delete[] Kcr;
    if (Jac) delete[] Jac;

    if (NcFull) delete[] NcFull;

    if (fp) 	 fclose(fp);
    if (fo) 	 fclose(fo);
    if (fTg)	 fclose(fTg);
    if (fRad)	 fclose(fRad);
    if (fChrg)	 fclose(fChrg);
    if (fp_full) fclose(fp_full);
    if (ftBolt)	 fclose(ftBolt);

    if (state_QSS)
    {
		for (int n=0;n<Natoms;n++)
			delete[] state_QSS[n];
		delete[] state_QSS;
    }
}

// Enable/Disable radiation in rate calculations
void CR_BKE::setRad(bool is_radiating)
{
	Rad = is_radiating;
}

void CR_BKE::setIrrad(bool is_radiated)
{
	Irrad = is_radiated;

	if (Irrad && !EnerEe)
	{
		printf("Currently unable to irradiate system without electron temperature\nNow exiting\n");
		exit(1);
	}
}

// Include tabulated rates to solver
void CR_BKE::addRates(KRatesGrp* pK)
{
	pRk = pK;
	NTstep = pK->NTstep;
	Natoms = pK->Natoms;
	Nstate = pK->Nstate;

	// Compile charged states
	Zq = new int[Nstate];
	for(int n=0;n<Nstate;n++) Zq[n] = 0;
	for(int a=0;a<Natoms;a++)
	{
		AtomX& A = *(pRk->atom[a]);
		int Nlev = pRk->Ngroups[a];
		for(int n=0;n<Nlev;n++) Zq[pRk->state[a][n]] = A.Z;
	}

	RadSEZ	= new double [Natoms];
	for(int n=0;n<Natoms;n++) RadSEZ[n] = 0.;

	RadSEbb	= new double [pK->KRates::NBBlns];
}

// Set up the electron temperature for electron-dominated processes
void CR_BKE::setTemp(double Te, double Tr)
{
    TeV = Te;
    TrV = Tr;
}

// Set chemistry variables to a particular solution
void CR_BKE::setsol(double* Ncr, double Ne)
{
    for(int n=0;n<NCRdim;n++) Nc[n] = Ncr[n];

    Nelec = Ne;
}

// Set solver temporal variables
void CR_BKE::setTimes(double ti, double tf, double tinc)
{
	tbeg = ti;	tcurr = ti;
	tend = tf;
	dt = tinc;

	tprint = ti;
	logtprt = log10(tinc);

	printf("\nStart time: %3.3e secs\nIncrementing time: %3.3e secs\nEnd time: %3.3e secs\n\n",tbeg,dt,tend);
}

// Set solver output variables
void CR_BKE::setPrintVars(bool prtlog, double pdt, double logdt)
{
	printlog 	= prtlog;
	dtprint 	= pdt;
	logdtprt 	= logdt;
}

// Initialize files for recording solution
void CR_BKE::initRecord()
{
	KRatesGrp& crk = *pRk;

	const char* dirParent 	= pRk->outputDir;
	const char* densGrp		= "DensGrp.txt";
	const char* densFul 	= "DensFul.txt";
	const char* chrgInfo 	= "ChrgInfo.txt";
	const char* globDens 	= "GlobDens.txt";
	const char* radEmiss 	= "RadEmiss.txt";
	const char*	boltzTimes	= "BoltzTimes.txt";
	const char* groupTemp	= "GroupTemps.txt";

	// --------- Section for Densities Output
	char currDir[100];
	strcpy(currDir,dirParent);
	strcat(currDir,densGrp);
	fp = fopen(currDir,"w");	// File opened: Track grouped densities
	currDir[0] = '\0';

	fprintf(fp,"# Time");
	for (int l = 0; l < NCRdim; l++) {
		stringstream ss;
		ss << "y" << l;
		string str = ss.str();
		fprintf(fp,",\t%s",str.c_str());
	}
	fprintf(fp,"\n");

	fflush(fp);

	// --------- Section for Densities Output
	strcpy(currDir,dirParent);
	strcat(currDir,densFul);
	fp_full = fopen(currDir,"w");	// File opened: Track all densities in all groups
	currDir[0] = '\0';

	fprintf(fp_full,"# Time");
	for (int l = 0; l < crk.Nlvful; l++) {
		stringstream ss;
		ss << "y" << l;
		string str = ss.str();
		fprintf(fp_full,"\t %s",str.c_str());
	}

	fprintf(fp_full,"\n");

	fflush(fp_full);

	// --------- Section for Charge Information Output
	strcpy(currDir,dirParent);
	strcat(currDir,chrgInfo);
	fChrg = fopen(currDir,"w");	// File opened: Track all densities in all groups
	currDir[0] = '\0';

	fprintf(fChrg,"# Time");
	for (int l = 0; l < Natoms; l++) {
		stringstream ss;
		ss << "CD" << l << "\tCF" << l;
		string str = ss.str();
		fprintf(fChrg,"\t%s",str.c_str());
	}

	fprintf(fChrg,"\n");

	fflush(fChrg);

	// --------- Section for Radiation Emission Output
	strcpy(currDir,dirParent);
	strcat(currDir,radEmiss);
	fRad = fopen(currDir,"w");	// File opened: Track all densities in all groups
	currDir[0] = '\0';

	fprintf(fRad,"# Time");
	for (int l = 0; l < Natoms; l++) {
		stringstream ss;
		ss << "Atom+" << l << "\t";
		string str = ss.str();
		fprintf(fRad,"\t%s",str.c_str());
	}

	fprintf(fRad,"\n");

	fflush(fRad);

	// --------- Section for Macroscopic Characteristics
	strcpy(currDir,dirParent);
	strcat(currDir,globDens);
	fo = fopen(currDir,"w");	// File opened: Track macroscopic densities and temperatures
	currDir[0] = '\0';

	fprintf(fo,"Time,\t\tdt,\t\tNe,\t\tZbar,\t\tTe,\t\trateNe,\t\trateEe,");
	fprintf(fo,"\t\tEexce,\t\tEdexe,\t\tEione,\t\tErece,");
	fprintf(fo,"\t\tEexcp,\t\tEexcs,\t\tEdexp,");
	fprintf(fo,"\t\tEionp,\t\tErecs,\t\tErecp,\t\tEfree,\t\tEtotal,\t\tEint\n");

	fflush(fo);

	// --------- Section for Boltzmann Plot Times
	strcpy(currDir,dirParent);
	strcat(currDir,boltzTimes);
	ftBolt = fopen(currDir,"w");	// File opened: Group temperature evolution
	currDir[0] = '\0';

	fprintf(ftBolt,"# File Number, Time Stamp (secs) \n");

	fflush(ftBolt);

	// --------- Section for Group Temperatures
	strcpy(currDir,dirParent);
	strcat(currDir,groupTemp);
	fTg = fopen(currDir,"w");	// File opened: Group temperature evolution
	currDir[0] = '\0';

	fprintf(fTg,"# Time, Temperatures \n");

	fflush(fTg);

	// ====================== Files opened for writing

	printf("Initialized files to write...\n\n");
}

void CR_BKE::monitorCollNe(double dt)
{
	double sqrtTeV = sqrt(TeV);
	double logTeV = log(TeV);
	double loglamb = 23.5 - log(sqrt(Nelec)/(TeV*sqrt(sqrtTeV))) - \
			sqrt(1e-5 + (logTeV-2)*(logTeV-2)/16);

	double collrate = 2.91e-6*Nelec*loglamb/(1e6*TeV*sqrtTeV);

	if (dt > 1./collrate) printf("[> e-e dt]");
}

void CR_BKE::recalcNeTe(double* Ncr)
{
	KRatesGrp& crk = *pRk;

	if 		(DensNe)			Nelec = Ncr[nev];
	else if (!FxdNe)
	{
		Nelec = 0.;
		if (crk.GrpScheme == "QSS")
		{
			for (int na=0;na<Nstate;na++)
				Nelec += (double)Zq[na]*Nc[na];
		}
		else
			for(int n=0;n<Nstate;n++)  	Nelec += Zq[n]*Ncr[n];
	}

	if 		(EnerEe) TeV = Ncr[nee]*2./(Nelec*3.);
}

// Determine Planckian radiation
double CR_BKE::get_PlanckField(double dE, double Tr)
{
	double hInu = 2*dE*dE*dE/(MKS_hc*MKS_hc*(exp(dE/Tr) - 1));
	return hInu;
}

/**
 * Bremmstrahlung:
 * see (5.21) p. 259 of [1]
 * @param  Zi  specie charge
 * @param  Te  electron temperature (eV)
 * @return  emissivity coefficient [eV*m^3/sec]
 */
double ZeldBrem(double Zi, double Te)
{
//        double coeff = 1.42e-40*(Zi*Zi*sqrt(Te+1.)); // [MKS formulation]
	double coeff = 9.54808e-20*(Zi*Zi*sqrt(Te+1.e-2)); // [eV*m^3/sec]
    return coeff;
}

void CR_BKE::invert_SS(double Ntot)
{
	KRates& crk = *pRk;

	int msize = Nstate;
	double *m_data =  cArray<double>(msize*msize);
	double *b_data  = cArray<double>(msize);
	int *Zq  = cArray<int>(msize);

	int n,m,na,ia;
	int low,upp;

	// find interpolation point in tables
	double lnTeV = max(crk.lnT_min,min(.99*crk.lnT_max,log(TeV) ) );
	float xt = (lnTeV-crk.lnT_min)/crk.dlnT;
	int it = (int) xt; xt -= it;

	// add-up excitations and deexcitations
	for (int bb = 0 ; bb < crk.NBBeXD; bb++)
	{
		na = crk.bbexd[bb].low[0];	// Determine ionic-atom stage
		n  = crk.bbexd[bb].low[1]; // Lower level
		m  = crk.bbexd[bb].upp[1]; // Upper level
		low = crk.state[na][n]; // Find corresponding state for lower level
		upp = crk.state[na][m]; // Find corresponding state for upper level
		//if(m != n+1) continue;
		// excitation (n -> m)
		double a_nm = (1.-xt) * crk.tab_keX[bb][it] + xt * crk.tab_keX[bb][it+1];
//		if (a_nm <= 0)
//			printf("a(%d -> %d) = %e\n",low,upp,a_nm);
//		a_nm = max(1.e-70,a_nm);
		m_data[upp*msize+low] += a_nm*Nelec;
		b_data[low]           += a_nm*Nelec;

		// deexcitation (m -> n)
		double b_nm = (1.-xt) * crk.tab_keD[bb][it] + xt * crk.tab_keD[bb][it+1];
//		if (b_nm <= 0)
//			printf("b(%d <- %d) = %e\n",low,upp,b_nm);
//		b_nm = max(1.e-70,b_nm);
		m_data[low*msize+upp] += b_nm*Nelec;
		b_data[upp]           += b_nm*Nelec;
	}

	// add-up ionizations and recombinations
	for(int bf=0;bf<crk.NBFeIR;bf++)
	{
		na = crk.bfeir[bf].low[0]; // neutral atom index
		ia = crk.bfeir[bf].upp[0]; // ionized atom index
		n  = crk.bfeir[bf].low[1]; // lower level (initial)
		m  = crk.bfeir[bf].upp[1]; // upper level
		low = crk.state[na][n]; // find corresponding state
		upp = crk.state[ia][m]; // find corresponding state
		// ionization (n -> m)
		double a_nm = (1.-xt) * crk.tab_keI[bf][it] + xt * crk.tab_keI[bf][it+1];
//		if (a_nm <= 0)
//			printf("a(%d -> %d) = %e\n",low,upp,a_nm);
//		a_nm = max(1.e-70,a_nm);
		m_data[upp*msize+low] += a_nm*Nelec;
		b_data[low]           += a_nm*Nelec;

		// recombination (m -> n)
		double b_nm = (1.-xt) * crk.tab_keR[bf][it] + xt * crk.tab_keR[bf][it+1];
//		if (b_nm <= 0)
//			printf("b(%d <- %d) = %e\n",low,upp,b_nm);
//		b_nm = max(1.e-70,b_nm);
		m_data[low*msize+upp] += b_nm*Nelec*Nelec;
		b_data[upp]           += b_nm*Nelec*Nelec;
	}

	// add-up autoionization and electron capture
	for(int bf=0;bf<crk.NBFaIR;bf++)
	{
		na = crk.bfair[bf].low[0]; // neutral atom index
		ia = crk.bfair[bf].upp[0]; // ionized atom index
		n  = crk.bfair[bf].low[1]; // lower level (initial)
		m  = crk.bfair[bf].upp[1]; // upper level
		low = crk.state[na][n]; // find corresponding state
		upp = crk.state[ia][m]; // find corresponding state
		// ionization (n -> m)
		double a_nm = crk.kAI[bf];
		printf("AI %i: %3.2e,\t",bf,a_nm);
//		if (a_nm <= 0)
//			printf("a(%d -> %d) = %e\n",low,upp,a_nm);
//		a_nm = max(1.e-70,a_nm);
		m_data[upp*msize+low] += a_nm;
		b_data[low]           += a_nm;

		// capture (m -> n)
		double b_nm = (1.-xt) * crk.tab_kDR[bf][it] + xt * crk.tab_kDR[bf][it+1];
		printf("EC %i: %3.2e\n",bf,b_nm);
//		if (b_nm <= 0)
//			printf("b(%d <- %d) = %e\n",low,upp,b_nm);
//		b_nm = max(1.e-70,b_nm);
		m_data[low*msize+upp] += b_nm*Nelec;
		b_data[upp]           += b_nm*Nelec;
	}

	if (Rad) {
		// add-up radiative (BB) transitions
		for(int bb=0;bb<crk.NBBlns;bb++)
		{
			na = crk.bbrad[bb].low[0];
			n  = crk.bbrad[bb].low[1]; // lower level
			m  = crk.bbrad[bb].upp[1]; // upper level
			low =crk.state[na][n]; // find corresponding state
			upp =crk.state[na][m]; // find corresponding state
			// spontaneous emission (m -> n)
			double A_nm = crk.rad_Abb[bb];
//			if (A_nm <= 0)
//				printf("A(%d <- %d) = %e\n",low,upp,A_nm);
//			A_nm = max(1.e-70,A_nm);
			m_data[low*msize+upp] += A_nm;
			b_data[upp]           += A_nm;
		}
		// add-up radiative (BF) transitions
		for(int bf=0;bf<crk.NBFpIR;bf++)
		{
			na = crk.bfpir[bf].low[0]; // neutral atom index
			ia = crk.bfpir[bf].upp[0]; // ionized atom index
			n  = crk.bfpir[bf].low[1]; // lower level (initial)
			m  = crk.bfpir[bf].upp[1]; // upper level
			low = crk.state[na][n]; // find corresponding state
			upp = crk.state[ia][m]; // find corresponding state
			// recombination (m -> n)
			double b_nm = (1.-xt) * crk.tab_kpR[bf][it][0] + xt * crk.tab_kpR[bf][it+1][0];
//			if (b_nm <= 0)
//				printf("b(%d <- %d) = %e\n",low,upp,b_nm);
//			b_nm = max(1.e-70,b_nm);
			m_data[low*msize+upp] += b_nm*Nelec;
			b_data[upp]           += b_nm*Nelec;
		}
	}

	for(int i=0;i<msize;i++) m_data[i*msize+i] -= b_data[i];
	for(int i=0;i<msize;i++) b_data[i] = 0.;

	// replace the last row for ion
	for(int a=0;a<Natoms;a++)
	{
		int Nlev = crk.Nalevs[a];
		for(int n=0;n<Nlev;n++) Zq[crk.state[a][n]] = crk.atom[a]->Z;
//		printf("Atom Charge %i: %i\n",a,crk.atom[a]->Z);
	}
	int mmax=msize-1;
	for(int i=0;i<msize;i++) m_data[mmax*msize+i] = 1;
	b_data[mmax] = Ntot;
//	for(int i=0;i<msize;i++) m_data[mmax*msize+i] = (double) Zq[i];
//	b_data[mmax] = Nelec;

	gsl_matrix_view mat = gsl_matrix_view_array (m_data, msize, msize);
	gsl_vector_view rhs = gsl_vector_view_array (b_data, msize);

	gsl_vector *sol = gsl_vector_alloc (msize);
	int _s;

	gsl_permutation * per = gsl_permutation_alloc (msize);
	gsl_linalg_LU_decomp (&mat.matrix, per, &_s);
	gsl_linalg_LU_solve  (&mat.matrix, per, &rhs.vector, sol);

//	printf ("x = \n");
//    gsl_vector_fprintf (stdout, sol, "%g");

	for(int i=0;i<msize;i++) b_data[i] = gsl_vector_get (sol,i);

	gsl_permutation_free (per);
	gsl_vector_free (sol);
	del_cArray<double>(msize*msize,m_data);
	del_cArray<int>(msize,Zq);

	for(int i=0;i<msize;i++) Nc[i] = b_data[i];

	del_cArray<double>(msize,b_data);
}

void CR_BKE::SolveSS(double Ntot)
{
	inst_fullN();
	invert_SS(Ntot);
	extractN_Grp(Nc);
	calcERates(NcFull);
	SolutionOutput();
}

// Include spectral calculations
//void CR_BKE::addSpecFeat(Spectrum* pSpec)
//{
//	pSpect = pSpec;
//}

// =========== GROUPING SECTION

void CR_BKE::inst_fullN()
{
	KRatesGrp& crk = *pRk;

	Nexlvls = crk.KRates::Nstate-crk.KRates::Natoms;
	NcFull	= new double[crk.KRates::Nstate];
	for (int n=0;n<crk.KRates::Nstate;n++) NcFull[n] = 0.;

	state_QSS = new int*[crk.KRates::Natoms];
	for(int na=0;na<Natoms;na++) state_QSS[na]  = new int[crk.KRates::Nalevs[na]];
	for(int na=0;na<Natoms;na++)
		for(int n=0;n<crk.KRates::Nalevs[na];n++)
		{
			if (n==0) continue;			// Ignore ground levels of ions
			state_QSS[na][n] = nq++;
			//printf("State %i, %i: %i\n",na,n,state[na][n]);
		}
}

// ====== Boltzmann grouping functions
void CR_BKE::inst_Tgroup()
{
	KRatesGrp& crk = *pRk;

	for(int na=0;na<crk.Natoms;na++)
		for(int m=0;m<crk.Ngroups[na];m++)
			crk.group[na][m].Tgroup = TeV;
}

void CR_BKE::extractN_Grp(double* Narr)
{
	double* Nq = Narr;

	KRatesGrp& crk = *pRk;

	for(int na=0;na<crk.Natoms;na++)
		for(int m=0;m<crk.Ngroups[na];m++)
		{
			int ka 		= crk.state[na][m];
			int nlvls 	= crk.group[na][m].Nlvlgrp;

			for(int t=0;t<nlvls;t++)
			{
				int n = crk.group[na][m].paridx[t];
				int lowf = crk.KRates::state[na][n];
				double ge = (double)crk.atom[na]->ge[n];

				if (Nq[ka] <= 0.)
					NcFull[lowf] = 1e-10*ge;
				else if ((crk.GrpScheme == "Boltzmann" || \
						crk.GrpScheme == "Custom_Boltzmann") && crk.group[na][m].grplvl)
				{
					double TgV	= crk.group[na][m].Tgroup;
					NcFull[lowf] = Nq[ka]*crk.boltz_weight(na,m,n,TgV);
				}
				else
				{
					double ggrp	= (double)crk.group[na][m].ge;
					NcFull[lowf] = Nq[ka]*ge/ggrp;
				}
			}
		}
}

void CR_BKE::boltz_temp(double *y)
{
	KRatesGrp& crk = *pRk;
	double* Nq = y;

	for(int na=0;na<crk.Natoms;na++)
		for(int m=0;m<crk.Ngroups[na];m++)
		{
			if (crk.group[na][m].baselvl)
			{
				int mplus = m+1;
				int low	 = crk.state[na][m];
				int grp	 = crk.state[na][mplus];

				double Ee = (double)crk.group[na][mplus].Ee - (double)crk.group[na][m].Ee;
				double ge = (double)crk.group[na][m].ge;
				double ggrp = (double)crk.group[na][mplus].ge;

				double TgV = -Ee/log(Nq[grp]*ge/(Nq[low]*ggrp));;
				if ((Nq[grp]*ge/(Nq[low]*ggrp)) <= 0. || (Nq[grp]*ge/(Nq[low]*ggrp)) >= 1.)
				{
					crk.group[na][m].Tgroup = INFINITY; continue;
				}
				crk.group[na][m].Tgroup = iter_Tgroup(na,m,TgV,Nq);
			}
			else if (crk.group[na][m].grplvl)
			{
				int mmins = m-1;
				crk.group[na][m].Tgroup = crk.group[na][mmins].Tgroup;
			}
		}
}

double CR_BKE::iter_Tgroup(int na, int n, double Tstar, double *y)
{
	KRatesGrp& crk = *pRk;
	double* Nq = y;

	int nplus = n+1;
	int low	 = crk.state[na][n];
	int grp	 = crk.state[na][nplus];

	double Ttemp = Tstar;	double Tf = Tstar;
	for (int iter = 0; iter < 15; iter++)
	{
		// 1st term (non-logarithmic)
		double ge = (double)crk.group[na][n].ge;
		double log1 = ge*Nq[grp]/Nq[low];

		// 2nd term (non-logarithmic)
		int lown0 = crk.group[na][n].paridx[0];
		double Qprstar = crk.part_funct_renorm(na,nplus,lown0,Ttemp);
		double log2 = Qprstar;

		// Denominator (non-logarithmic)
		double QEprstar = crk.part_funct_E1mom(na,nplus,lown0,Ttemp);
		double denom = QEprstar;

		Ttemp += (log1 - log2)*Ttemp*Ttemp/denom;
		if (abs(Ttemp-Tf)/Ttemp >= 1.e-3)	Tf = Ttemp;
		else								break;
	}

	return Tf;
}

double CR_BKE::calc_xiEnergy(int na, int n, int nnorm, double Tgroup, bool grp)
{
	KRatesGrp& crk = *pRk;

	double Emom0 = crk.part_funct_renorm(na,n,nnorm,Tgroup);
	double Emom1 = crk.part_funct_E1mom(na,n,nnorm,Tgroup);
	double Emom2 = crk.part_funct_E2mom(na,n,nnorm,Tgroup);

	double xiterm = Emom2/(Emom1*Emom0);

	if (!grp)
		xiterm -= Emom1/Emom0;

	return xiterm;
}

double CR_BKE::calc_deltEE(int na, int n, int ia, int m, double* y, double Tgroup)
{
	KRatesGrp& crk = *pRk;
	double* Nq = y;

	double Elow = crk.group[na][n].Ee;
	if (crk.group[na][n].baselvl)
	{
		int nsubgp = n+1;
		double Nrat = Nq[nsubgp]/Nq[n];
		Elow -= Nrat*calc_xiEnergy(na,nsubgp,n,Tgroup,false);
	}
	else if (crk.group[na][n].grplvl)
	{
		int nbase = n-1;
		Elow += calc_xiEnergy(na,n,nbase,Tgroup,true);
	}

	double Eupp = crk.group[ia][m].Ee;
	if (crk.group[ia][m].baselvl)
	{
		int msubgp = m+1;
		double Nrat = Nq[msubgp]/Nq[m];
		Eupp += Nrat*calc_xiEnergy(ia,msubgp,m,Tgroup,false);
	}
	else if (crk.group[ia][m].grplvl)
	{
		int mbase = m-1;
		Eupp -= calc_xiEnergy(ia,m,mbase,Tgroup,true);
	}

	double Egap = Eupp - Elow;
	return Egap;
}

double CR_BKE::calc_deltIE(int na, int n, int ia, int m, double* y, double Tgroup)
{
	KRatesGrp& crk = *pRk;
	double* Nq = y;

	double IEgnd = 0.;
	bool IEthresh = false;

	for(int bf=0;bf<crk.KRates::NBFeIR;bf++)
	{
		int latom	= crk.KRates::bfeir[bf].low[0];
		int low		= crk.KRates::bfeir[bf].low[1];
		int uatom 	= crk.KRates::bfeir[bf].upp[0];
		int upp 	= crk.KRates::bfeir[bf].upp[1];

		if (latom == na && uatom == ia && low == 0 && upp  == 0)
		{
			IEgnd = crk.KRates::bfeir[bf].dE;
			IEthresh = true;
			break;
		}
	}

	if (!IEthresh)
	{
        std::cout << "Ionization energy for ground level not found:" << crk.KRates::atom[na][0].name << " ," << crk.KRates::atom[na][0].Z;
//		printf("Ionization energy for ground level not found: %s, %i",
//				crk.KRates::atom[na][0].name,crk.KRates::atom[na][0].Z);
		exit(1);
	}

	double Elow = crk.group[na][n].Ee;
	if (crk.group[na][n].baselvl)
	{
		int nsubgp = n+1;
		double Nrat = Nq[nsubgp]/Nq[n];
		Elow -= Nrat*calc_xiEnergy(na,nsubgp,n,Tgroup,false);
	}
	else if (crk.group[na][n].grplvl)
	{
		int nbase = n-1;
		Elow += calc_xiEnergy(na,n,nbase,Tgroup,true);
	}

	double Eupp = crk.group[ia][m].Ee;
	if (crk.group[ia][m].baselvl)
	{
		int msubgp = m+1;

		double Emom0 = crk.part_funct_renorm(ia,msubgp,m,Tgroup);
		double Emom1 = crk.part_funct_E1mom(ia,msubgp,m,Tgroup);
		Eupp += Emom1/Emom0;

		if (crk.group[na][n].baselvl || crk.group[na][n].grplvl)
		{
			int nsubgp = n+1;
			double Nsubrat = 0.;
			if (crk.group[na][n].baselvl)	Nsubrat = Nq[nsubgp]/Nq[msubgp];
			else							Nsubrat = Nq[n]/Nq[msubgp];

			double Nrat = Nq[msubgp]/Nq[m];
			Eupp -= Nsubrat*Nrat*calc_xiEnergy(ia,msubgp,m,Tgroup,false);
		}
	}
	else if (crk.group[ia][m].grplvl)
	{
		int mbase = m-1;

		double Emom0 = crk.part_funct_renorm(ia,m,mbase,Tgroup);
		double Emom1 = crk.part_funct_E1mom(ia,m,mbase,Tgroup);
		Eupp += Emom1/Emom0;

		if (crk.group[na][n].baselvl || crk.group[na][n].grplvl)
		{
			int nsubgp = n+1;
			double Nsubrat = 0.;
			if (crk.group[na][n].baselvl)	Nsubrat = Nq[nsubgp]/Nq[m];
			else							Nsubrat = Nq[n]/Nq[m];

			Eupp -= Nsubrat*calc_xiEnergy(ia,m,mbase,Tgroup,false);
		}
	}

	double Egap = Eupp - (IEgnd - Elow);
	return Egap;
}

// ====== Quasi-steady-state (QSS) functions
void CR_BKE::setGnds_QSS()
{
	KRatesGrp& crk = *pRk;

	int na,low,lowf;
	for (na=0;na<crk.Natoms;na++)
	{
		low = crk.state[na][0];
		lowf = crk.KRates::state[na][0]; 	// Find corresponding variable for lower level

		NcFull[lowf] = Nc[low];
	}
}

void CR_BKE::extractN_QSS(double* y)
{
	KRatesGrp& crk = *pRk;
	double* Nq = y;

	for (int a=0;a<Natoms;a++)
	{
		if (crk.KRates::Nalevs[a] == 1) continue;

		int msize = crk.KRates::Nalevs[a]-1;
		double *m_data =  cArray<double>(msize*msize);
		double *b_data  = cArray<double>(msize);
	//	int *Zq  = cArray<int>(msize);

		int n,m,na,ia;
		int low,upp,lowf,uppf;
		bool gnd,ion;
		gnd = false;
		ion = false;

		// find interpolation point in tables
		double lnTeV = max(crk.lnT_min,min(.99*crk.lnT_max,log(TeV) ) );
		float xt = (lnTeV-crk.lnT_min)/crk.dlnT;
		int it = (int) xt; xt -= it;

		// --------------------------------------
		// Add up excitations and de-excitations
		// --------------------------------------
		for (int bb=0;bb<crk.KRates::NBBeXD;bb++)
		{
			na 	= crk.KRates::bbexd[bb].low[0];
			n 	= crk.KRates::bbexd[bb].low[1];
			ia 	= crk.KRates::bbexd[bb].upp[0];
			m 	= crk.KRates::bbexd[bb].upp[1];

			if (a!=na)		continue;
			if (n==0)		gnd = true;

			// Excitation
			double a_nm = (1.-xt) * crk.KRates::tab_keX[bb][it] + xt * crk.KRates::tab_keX[bb][it+1];
			// De-excitation
			double b_nm = (1.-xt) * crk.KRates::tab_keD[bb][it] + xt * crk.KRates::tab_keD[bb][it+1];

			low = n-1;
			upp = m-1;

			lowf = crk.KRates::state[na][n]; 	// Find corresponding variable for lower level
			uppf = crk.KRates::state[na][m]; 	// Find corresponding variable for upper level

	//		printf("a(%d -> %d) = %e\n",low,upp,a_nm);
			if (gnd) 	b_data[upp] -= a_nm*Nelec*Nq[lowf];
			else
			{
				m_data[low*msize+low] -= a_nm*Nelec;
				m_data[upp*msize+low] += a_nm*Nelec;
			}

	//		printf("b(%d <- %d) = %e\n",low,upp,b_nm);
			if (gnd) 	m_data[upp*msize+upp] -= b_nm*Nelec;
			else
			{
				m_data[low*msize+upp] += b_nm*Nelec;
				m_data[upp*msize+upp] -= b_nm*Nelec;
			}

			gnd = false;
			ion = false;
		}

		// --------------------------------------
		// Add up ionizations and recombinations
		// --------------------------------------
		for(int bf=0;bf<crk.KRates::NBFeIR;bf++)
		{
			na 	= crk.KRates::bfeir[bf].low[0];
			n 	= crk.KRates::bfeir[bf].low[1];
			ia 	= crk.KRates::bfeir[bf].upp[0];
			m 	= crk.KRates::bfeir[bf].upp[1];

			if (a!=na)			continue;
			if (n!=0 && m==0)	ion = true;
			else 				continue;

			// Ionization
			double a_nm = (1.-xt) * crk.KRates::tab_keI[bf][it] + xt * crk.KRates::tab_keI[bf][it+1];
			// Recombination
			double b_nm = (1.-xt) * crk.KRates::tab_keR[bf][it] + xt * crk.KRates::tab_keR[bf][it+1];

			low = n-1;

			lowf = crk.KRates::state[na][n]; 	// Find corresponding variable for lower level
			uppf = crk.KRates::state[ia][m]; 	// Find corresponding variable for upper level

	//		printf("a(%d -> %d) = %e\n",low,upp,a_nm);
	//		a_nm = max(1.e-60,a_nm);
			if (ion) 	m_data[low*msize+low] -= a_nm*Nelec;

	//		printf("b(%d <- %d) = %e\n",low,upp,b_nm);
	//		b_nm = max(1.e-60,b_nm);
			if (ion) 	b_data[low] -= b_nm*Nelec*Nelec*Nq[uppf];

			gnd = false;
			ion = false;
		}

		if (Rad) {
			// Radiative (BB) transitions
			for(int bb=0;bb<crk.KRates::NBBlns;bb++)
			{
				na 	= crk.KRates::bbrad[bb].low[0];
				n 	= crk.KRates::bbrad[bb].low[1];
				ia 	= crk.KRates::bbrad[bb].upp[0];
				m 	= crk.KRates::bbrad[bb].upp[1];

				if (a!=na)		continue;
				if (n==0)		gnd = true;

				double A_nm = crk.KRates::rad_Abb[bb];

				low = n-1;
				upp = m-1;

				lowf = crk.KRates::state[na][n]; 	// Find corresponding variable for lower level
				uppf = crk.KRates::state[na][m]; 	// Find corresponding variable for upper level

				if (gnd) m_data[upp*msize+upp] -= A_nm;
				else
				{
					m_data[low*msize+upp] += A_nm;
					m_data[upp*msize+upp] -= A_nm;
				}

				gnd = false;
				ion = false;
			}

			// add-up radiative (BF) transitions
			for(int bf=0;bf<crk.KRates::NBFeIR;bf++)
			{
				na 	= crk.KRates::bfpir[bf].low[0];
				n 	= crk.KRates::bfpir[bf].low[1];
				ia 	= crk.KRates::bfpir[bf].upp[0];
				m 	= crk.KRates::bfpir[bf].upp[1];

				if (a!=na)		continue;
				if (n!=0 && m==0)	ion = true;
				else				continue;

				double b_nm = (1.-xt) * crk.KRates::tab_kpR[bf][it][0] + xt * crk.KRates::tab_kpR[bf][it+1][0];

				low = n-1;
				upp = m-1;

				lowf = crk.KRates::state[na][n]; 	// Find corresponding variable for lower level
				uppf = crk.KRates::state[ia][m]; 	// Find corresponding variable for upper level

				if (ion) b_data[low] -= b_nm*Nelec*Nq[uppf];

				gnd = false;
				ion = false;
			}
		}

//		for (int row=0;row<msize;row++)
//		{
//			for (int col=0;col<msize;col++)
//				printf("%2.1e ",m_data[row*msize+col]);
//			printf("\n");
//		}

		gsl_matrix_view mat = gsl_matrix_view_array (m_data, msize, msize);
		gsl_vector_view rhs = gsl_vector_view_array (b_data, msize);

		gsl_vector *sol = gsl_vector_alloc (msize);
		int _s;

		gsl_permutation * per = gsl_permutation_alloc (msize);
		gsl_linalg_LU_decomp (&mat.matrix, per, &_s);
		gsl_linalg_LU_solve  (&mat.matrix, per, &rhs.vector, sol);

	//	printf ("x = \n");
	//    gsl_vector_fprintf (stdout, sol, "%g");

		for(int i=0;i<msize;i++) b_data[i] = gsl_vector_get (sol,i);

		gsl_permutation_free (per);
		gsl_vector_free (sol);
		del_cArray<double>(msize*msize,m_data);

		for (int num=0;num<msize;num++)
		{
			int lvl = num+1;
			lowf = crk.KRates::state[a][lvl];

			NcFull[lowf] = b_data[num];
		}
		del_cArray<double>(msize,b_data);

	}
}

void CR_BKE::calcERates(double *y)
{
	int n,m,na,ia;
	int low,upp;
	double grat,dE;

	KRatesGrp& crk = *pRk;
	double* Nq = y;

	Eexce	= Edexe	= Eione	= Erece	= 0.;
	Eexcp	= Eionp	= 0.;
	Eexcs	= Eions	= 0.;
	Edexp	= Erecp	= 0.;
	Efree	= 0.;
	for (na=0;na<Natoms;na++) RadSEZ[na] = 0.;
	for (int bb=0;bb<crk.KRates::NBBlns;bb++) RadSEbb[bb] = 0.;

    // Find interpolation point in tables
    double lnTeV = max(crk.lnT_min,min(.99*crk.lnT_max,log(TeV) ) );
    float xt = (lnTeV-crk.lnT_min)/crk.dlnT;
    int it = (int) xt; xt -= it;

    // ----------------------------------------
    // Electron-Impact: Bound-Bound Transitions
    // ----------------------------------------
    for(int bb=0;bb<crk.KRates::NBBeXD;bb++)
    {
		na 	= crk.KRates::bbexd[bb].low[0];
		n 	= crk.KRates::bbexd[bb].low[1];
		ia 	= crk.KRates::bbexd[bb].upp[0];
		m 	= crk.KRates::bbexd[bb].upp[1];
		dE	= crk.KRates::bbexd[bb].dE;

		low = crk.KRates::state[na][n]; 	// Find corresponding group for lower level
		upp = crk.KRates::state[na][m]; 	// Find corresponding group for upper level

		// Excitation
		double a_nm = (1.-xt) * crk.KRates::tab_keX[bb][it] + xt * crk.KRates::tab_keX[bb][it+1];
		// De-excitation
		double b_nm = (1.-xt) * crk.KRates::tab_keD[bb][it] + xt * crk.KRates::tab_keD[bb][it+1];

        // Track energy rates
        Eexce += (a_nm*Nelec*Nq[low])*dE;
		Edexe += (b_nm*Nelec*Nq[upp])*dE;
    }

    // ---------------------------------------
    // Electron-Impact: Bound-Free Transitions
    // ---------------------------------------
    for(int bf=0;bf<crk.KRates::NBFeIR;bf++)
    {
		na 	= crk.KRates::bfeir[bf].low[0];
		n 	= crk.KRates::bfeir[bf].low[1];
		ia 	= crk.KRates::bfeir[bf].upp[0];
		m 	= crk.KRates::bfeir[bf].upp[1];
		dE	= crk.KRates::bfeir[bf].dE;

		if (crk.GrpScheme == "QSS" && m!=0) continue;

		low = crk.KRates::state[na][n]; 	// Find corresponding group for lower level
		upp = crk.KRates::state[ia][m]; 	// Find corresponding group for upper level

		// Ionization
		double a_nm = (1.-xt) * crk.KRates::tab_keI[bf][it] + xt * crk.KRates::tab_keI[bf][it+1];
		// Recombination
		double b_nm = (1.-xt) * crk.KRates::tab_keR[bf][it] + xt * crk.KRates::tab_keR[bf][it+1];

        // Track energy rates
        Eione += (a_nm*Nelec*Nq[low])*dE;
		Erece += (b_nm*Nelec*Nelec*Nq[upp])*dE;
    }

    if (Rad) {
       // ------------------------------------------
    	// Radiation-Induced: Bound-Bound Transitions
    	// ------------------------------------------
		for(int bb=0;bb<crk.KRates::NBBlns;bb++)
		{
			na 	= crk.KRates::bbrad[bb].low[0];
			n 	= crk.KRates::bbrad[bb].low[1];
			ia 	= crk.KRates::bbrad[bb].upp[0];
			m 	= crk.KRates::bbrad[bb].upp[1];
			dE	= crk.KRates::bbrad[bb].dE;

			low = crk.KRates::state[na][n]; 	// Find corresponding group for lower level
			upp = crk.KRates::state[ia][m]; 	// Find corresponding group for upper level

			grat = (double)crk.KRates::atom[ia]->ge[m]/(double)crk.KRates::atom[na]->ge[n];

			// Photoexcitation (n -> m)
			double B_nm = grat*MKS_hc*MKS_hc*crk.KRates::rad_Abb[bb]/(2*dE*dE*dE);
			// Stimulated emission (m -> n)
			double B_mns = B_nm/grat;
			// Spontaneous Emission (m -> n)
			double A_mn = crk.KRates::rad_Abb[bb];

			// Track energy rates
			if (Irrad) {
				Eexcp += (B_nm*Nq[low]*get_PlanckField(dE,TrV))*dE;
				Eexcs += (B_mns*Nq[upp]*get_PlanckField(dE,TrV))*dE;
			}
			Edexp += (A_mn*Nq[upp])*dE;

			RadSEZ[na] += (A_mn*Nq[upp])*dE;

			if (na != 2) continue;
			RadSEbb[bb] += (A_mn*Nq[upp])*dE;
		}

       // -----------------------------------------
    	// Radiation-Induced: Bound-Free Transitions
    	// -----------------------------------------
		for(int bf=0;bf<crk.KRates::NBFpIR;bf++)
		{
			na 	= crk.KRates::bfpir[bf].low[0];
			n 	= crk.KRates::bfpir[bf].low[1];
			ia 	= crk.KRates::bfpir[bf].upp[0];
			m 	= crk.KRates::bfpir[bf].upp[1];

			if (crk.GrpScheme == "QSS" && m!=0) continue;

			low = crk.KRates::state[na][n]; 	// Find corresponding group for lower level
			upp = crk.KRates::state[ia][m]; 	// Find corresponding group for upper level

			double a_nmE, s_nmE;
			a_nmE = s_nmE = 0.;
			if (Irrad) {
				// Photoionization (n -> m)
				a_nmE = crk.tab_kpI[bf][1];
				// Stimulated recombination (m -> n)
				s_nmE = (1.-xt) * crk.tab_kpS[bf][it][1] + xt * crk.tab_kpS[bf][it+1][1];
			}
			// Radiative Recombination (m -> n)
			double b_nmE = (1.-xt) * crk.KRates::tab_kpR[bf][it][1] + xt * crk.KRates::tab_kpR[bf][it+1][1];

			// Track energy rates
			if (Irrad) {
				Eionp += (a_nmE*Nq[low])*TrV;
				Eions += (s_nmE*Nelec*Nq[upp])*TeV;
			}
			Erecp += (b_nmE*Nelec*Nq[upp])*TeV;
		}
    }

    // ----------------------------------
    // Account for Bremmstrahlung losses
    // ----------------------------------
	for (int n=0;n<Nstate;n++)
	{
		double Z = (double) Zq[n];
		double kBrem = ZeldBrem(Z,TeV);

		Efree += kBrem*Nq[n]*Nelec;
	}
}

// Simultaneous calculation of the right-hand side vector and Jacobian
void CR_BKE::source(double *y)
{
    int n,m,na,ia,gpl,gpu;
    int low,upp;
    double grat,dE;

    KRatesGrp& crk = *pRk;
    double* Nq = y;

    for(n=0;n<NCRdim;n++)
    {
        Jcr[n] = 0.0;
        dNc[n] = 0.0;
        Scr[n] = 0.0;
        for(m=0;m<NCRdim;m++) Kcr[n][m] = 0.;
        for(m=0;m<NCRdim;m++) Jac[n][m] = 0.;
    }

    // Find interpolation point in tables
    double lnTeV = max(crk.lnT_min,min(.99*crk.lnT_max,log(TeV) ) );
    float xt = (lnTeV-crk.lnT_min)/crk.dlnT;
    int it = (int) xt; xt -= it;

    // ----------------------------------------
	// Electron-Impact: Bound-Bound Transitions
	// ----------------------------------------
    for(int bb=0;bb<crk.NBBeXD;bb++)
    {
  		na 	= crk.bbexd[bb].low[0];
  		n 	= crk.bbexd[bb].low[1];
  		ia 	= crk.bbexd[bb].upp[0];
  		m 	= crk.bbexd[bb].upp[1];
  		dE	= crk.bbexd[bb].dE;

  		gpl = crk.gstate[na][n];
  		gpu = crk.gstate[ia][m];

  		if (na == ia && gpl == gpu) continue;

		low = crk.state[na][gpl]; 	// Find corresponding group for lower level
		upp = crk.state[na][gpu]; 	// Find corresponding group for upper level

		// Excitation
		double a_nm = (1.-xt) * crk.tab_keX[bb][it] + xt * crk.tab_keX[bb][it+1];
		double dlnTe = ((crk.tab_lnT[it+1]) - (crk.tab_lnT[it]));
		double a_nm_der = (crk.tab_keX[bb][it+1] - crk.tab_keX[bb][it]) / (dlnTe * TeV);
		// De-excitation
		double b_nm = (1.-xt) * crk.tab_keD[bb][it] + xt * crk.tab_keD[bb][it+1];
		double b_nm_der = (crk.tab_keD[bb][it+1] - crk.tab_keD[bb][it]) / (dlnTe * TeV);

		double low_multx = 1.;
		double upp_multx = 1.;
		if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
		{
			if (crk.group[na][gpl].grplvl)
			{
				double TgV = crk.group[na][gpl].Tgroup;
				low_multx = crk.boltz_weight(na,gpl,n,TgV);
			}
			if (crk.group[na][gpu].grplvl)
			{
				double TgV = crk.group[na][gpu].Tgroup;
				upp_multx = crk.boltz_weight(na,gpu,m,TgV);
			}
		}
		else
		{
			low_multx = (double)crk.atom[na]->ge[n]/(double)crk.group[na][gpl].ge;
			upp_multx = (double)crk.atom[na]->ge[m]/(double)crk.group[na][gpu].ge;
		}

		// Excitation (n -> m)
        Jcr[low]      += (a_nm*Nelec*low_multx);
        Kcr[upp][low] += (a_nm*Nelec*low_multx);
        Jac[low][low] -= (a_nm*Nelec*low_multx);
        Jac[upp][low] += (a_nm*Nelec*low_multx);
        // De-excitation (m -> n)
        Jcr[upp]      += (b_nm*Nelec*upp_multx);
        Kcr[low][upp] += (b_nm*Nelec*upp_multx);
        Jac[low][upp] += (b_nm*Nelec*upp_multx);
        Jac[upp][upp] -= (b_nm*Nelec*upp_multx);

        if (EnerEe) {
        	// Excitation (n -> m)
			Kcr[nee][low] -= (a_nm*Nelec*low_multx)*dE;
			Jac[low][nee] -= (a_nm_der*Nq[low]*low_multx*2./3.);
			Jac[upp][nee] += (a_nm_der*Nq[low]*low_multx*2./3.);
			Jac[nee][low] -= (a_nm*Nelec*low_multx)*dE; // New term
			Jac[nee][nee] -= (a_nm_der*Nq[low]*low_multx*2./3.)*dE;
			// De-excitation (m -> n)
			Kcr[nee][upp] += (b_nm*Nelec*upp_multx)*dE;
			Jac[upp][nee] -= (b_nm_der*Nq[upp]*upp_multx*2./3.);
			Jac[low][nee] += (b_nm_der*Nq[upp]*upp_multx*2./3.);
			Jac[nee][upp] += (b_nm*Nelec*upp_multx)*dE; // New term
			Jac[nee][nee] += (b_nm_der*Nq[upp]*upp_multx*2./3.)*dE;
        }

        if (DensNe) {
        	// Excitation (n -> m)
        	Jac[low][nev] -= (a_nm*Nq[low]*low_multx);
			Jac[upp][nev] += (a_nm*Nq[low]*low_multx);
			// De-excitation (m -> n)
			Jac[upp][nev] -= (b_nm*Nq[upp]*upp_multx);
			Jac[low][nev] += (b_nm*Nq[upp]*upp_multx);

			if (EnerEe) {
				// Excitation (n -> m)
				Jac[nee][nev] -= (a_nm*Nq[low]*low_multx)*dE;
				// De-excitation (m -> n)
				Jac[nee][nev] += (b_nm*Nq[upp]*upp_multx)*dE;
			}
        }
    }

    // ---------------------------------------
	// Electron-Impact: Bound-Free Transitions
	// ---------------------------------------
    for(int bf=0;bf<crk.NBFeIR;bf++)
    {
    	na 	= crk.bfeir[bf].low[0];
		n 	= crk.bfeir[bf].low[1];
		ia 	= crk.bfeir[bf].upp[0];
		m 	= crk.bfeir[bf].upp[1];
		dE	= crk.bfeir[bf].dE;

  		gpl = crk.gstate[na][n];
  		gpu = crk.gstate[ia][m];

  		if (na == ia && gpl == gpu) continue;

		low = crk.state[na][gpl]; 	// Find corresponding group for lower level
		upp = crk.state[ia][gpu]; 	// Find corresponding group for upper level

		// Ionization
		double a_nm = (1.-xt) * crk.tab_keI[bf][it] + xt * crk.tab_keI[bf][it+1];
		double dlnTe = ((crk.tab_lnT[it+1]) - (crk.tab_lnT[it]));
		double a_nm_der = (crk.tab_keI[bf][it+1] - crk.tab_keI[bf][it]) / (dlnTe * TeV);
		// Recombination
		double b_nm = (1.-xt) * crk.tab_keR[bf][it] + xt * crk.tab_keR[bf][it+1];
		double b_nm_der = (crk.tab_keR[bf][it+1] - crk.tab_keR[bf][it]) / (dlnTe * TeV);

		double low_multx = 1.;
		double upp_multx = 1.;
		if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
		{
			if (crk.group[na][gpl].grplvl)
			{
				double TgV = crk.group[na][gpl].Tgroup;
				low_multx = crk.boltz_weight(na,gpl,n,TgV);
			}
			if (crk.group[ia][gpu].grplvl)
			{
				double TgV = crk.group[ia][gpu].Tgroup;
				upp_multx = crk.boltz_weight(ia,gpu,m,TgV);
			}
		}
		else
		{
			low_multx = (double)crk.atom[na]->ge[n]/(double)crk.group[na][gpl].ge;
			upp_multx = (double)crk.atom[ia]->ge[m]/(double)crk.group[ia][gpu].ge;
		}

		// Ionization (n -> m)
        Jcr[low]      += (a_nm*Nelec*low_multx);
        Kcr[upp][low] += (a_nm*Nelec*low_multx);
        Jac[low][low] -= (a_nm*Nelec*low_multx);
        Jac[upp][low] += (a_nm*Nelec*low_multx);
        // Recombination (m -> n)
        Jcr[upp]      += (b_nm*Nelec*Nelec*upp_multx);
        Kcr[low][upp] += (b_nm*Nelec*Nelec*upp_multx);
        Jac[low][upp] += (b_nm*Nelec*Nelec*upp_multx);
        Jac[upp][upp] -= (b_nm*Nelec*Nelec*upp_multx);

        if (EnerEe) {
        	// Ionization (n -> m)
			Kcr[nee][low] -= (a_nm*Nelec*low_multx)*dE;
			Jac[nee][low] -= (a_nm*Nelec*low_multx)*dE; // New term
			Jac[low][nee] -= (a_nm_der*Nq[low]*low_multx*2./3.);
			Jac[upp][nee] += (a_nm_der*Nq[low]*low_multx*2./3.);
			Jac[nee][nee] -= (a_nm_der*Nq[low]*low_multx*2./3.)*dE;
			// Recombination (m -> n)
			Kcr[nee][upp] += (b_nm*Nelec*Nelec*upp_multx)*dE;
			Jac[nee][upp] += (b_nm*Nelec*Nelec*upp_multx)*dE; // New term
			Jac[upp][nee] -= (b_nm_der*Nq[upp]*Nelec*upp_multx*2./3.);
			Jac[low][nee] += (b_nm_der*Nq[upp]*Nelec*upp_multx*2./3.);
			Jac[nee][nee] += (b_nm_der*Nq[upp]*Nelec*upp_multx*2./3.)*dE;
        }

        if (DensNe) {
        	// Ionization (n -> m)
        	Jcr[nev]  	  -= (a_nm*Nq[low]*low_multx);
        	Jac[low][nev] -= (a_nm*Nq[low]*low_multx);
			Jac[upp][nev] += (a_nm*Nq[low]*low_multx);
			Jac[nev][low] += (a_nm*Nelec*low_multx);
			Jac[nev][nev] += (a_nm*Nq[low]*low_multx);
			// Recombination (m -> n)
			Jcr[nev]  	  += (b_nm*Nelec*Nq[upp]*upp_multx);
			Jac[upp][nev] -= (b_nm*Nq[upp]*2.*Nelec*upp_multx);
			Jac[low][nev] += (b_nm*Nq[upp]*2.*Nelec*upp_multx);
			Jac[nev][upp] -= (b_nm*Nelec*Nelec*upp_multx);
			Jac[nev][nev] -= (b_nm*Nq[upp]*2.*Nelec*upp_multx);

			if (EnerEe) {
				// Ionization (n -> m)
				Jac[nev][nee] += (a_nm_der*Nq[low]*low_multx*2./3.);
				Jac[nee][nev] -= (a_nm*Nq[low]*low_multx)*dE;
				// Recombination (m -> n)
				Jac[nev][nee] -= (b_nm_der*Nq[upp]*Nelec*upp_multx*2./3.);
				Jac[nee][nev] += (b_nm*Nq[upp]*2*Nelec*upp_multx)*dE;
			}
        }
    }

	// -----------------------------------
	// Autoionization and Electron Capture
	// -----------------------------------
//    for(int bf=0;bf<crk.NBFaIR;bf++)
//	{
//  		na 	= crk.bfair[bf].low[0];
//  		n 	= crk.bfair[bf].low[1];
//  		ia 	= crk.bfair[bf].upp[0];
//  		m 	= crk.bfair[bf].upp[1];
//  		dE	= crk.bfair[bf].dE;
//
//  		gpl = crk.gstate[na][n];
//  		gpu = crk.gstate[ia][m];
//
//  		if (na == ia && gpl == gpu) continue;
//
//		low = crk.state[na][gpl]; 	// Find corresponding group for lower level
//		upp = crk.state[ia][gpu]; 	// Find corresponding group for upper level
//
//		// Autoionization
//		double a_nm = crk.kAI[bf];
//		// Electron Capture
//		double b_nm = (1.-xt) * crk.tab_kDR[bf][it] + xt * crk.tab_kDR[bf][it+1];
//		double dlnTe = ((crk.tab_lnT[it+1]) - (crk.tab_lnT[it]));
//		double b_nm_der = (crk.tab_kDR[bf][it+1] - crk.tab_kDR[bf][it]) / (dlnTe * TeV);
//
//		double low_multx = 1.;
//		double upp_multx = 1.;
//		if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
//		{
//			if (crk.group[na][gpl].grplvl)
//			{
//				double TgV = crk.group[na][gpl].Tgroup;
//				low_multx = crk.boltz_weight(na,gpl,n,TgV);
//			}
//			if (crk.group[ia][gpu].grplvl)
//			{
//				double TgV = crk.group[ia][gpu].Tgroup;
//				upp_multx = crk.boltz_weight(ia,gpu,m,TgV);
//			}
//		}
//		else
//		{
//			low_multx = (double)crk.atom[na]->ge[n]/(double)crk.group[na][gpl].ge;
//			upp_multx = (double)crk.atom[ia]->ge[m]/(double)crk.group[ia][gpu].ge;
//		}
//
//		// Autoionization (n -> m)
//        Jcr[low]      += (a_nm*low_multx);
//        Kcr[upp][low] += (a_nm*low_multx);
//        Jac[low][low] -= (a_nm*low_multx);
//        Jac[upp][low] += (a_nm*low_multx);
//        // Electron Capture (m -> n)
//        Jcr[upp]      += (b_nm*Nelec*upp_multx);
//        Kcr[low][upp] += (b_nm*Nelec*upp_multx);
//        Jac[low][upp] += (b_nm*Nelec*upp_multx);
//        Jac[upp][upp] -= (b_nm*Nelec*upp_multx);
//
//        if (EnerEe) {
//        	// Autoionization (n -> m)
//			Kcr[nee][low] += (a_nm*low_multx)*dE;
//			Jac[nee][low] += (a_nm*low_multx)*dE;
//			// Electron Capture (m -> n)
//			Kcr[nee][upp] -= (b_nm*Nelec*upp_multx)*dE;
//			Jac[nee][upp] -= (b_nm*Nelec*upp_multx)*dE;
//			Jac[upp][nee] -= (b_nm_der*Nq[upp]*upp_multx*2./3.);
//			Jac[low][nee] += (b_nm_der*Nq[upp]*upp_multx*2./3.);
//			Jac[nee][nee] -= (b_nm_der*Nq[upp]*upp_multx*2./3.)*dE;
//        }
//
//        if (DensNe) {
//        	// Autoionization (n -> m)
//        	Kcr[nev][low] += (a_nm*low_multx);
//			Jac[nev][low] += (a_nm*low_multx);
//			// Electron Capture (m -> n)
//			Jcr[nev]  	  	+= (b_nm*Nq[upp]*upp_multx);
//			Jac[upp][nev] -= (b_nm*Nq[upp]*upp_multx);
//			Jac[low][nev] += (b_nm*Nq[upp]*upp_multx);
//			Jac[nev][upp] -= (b_nm*Nelec*upp_multx);
//			Jac[nev][nev] -= (b_nm*Nq[upp]*upp_multx);
//
//			if (EnerEe) {
//				// Electron Capture (m -> n)
//				Jac[nev][nee] -= (b_nm_der*Nq[upp]*upp_multx*2./3.);
//				Jac[nee][nev] -= (b_nm*Nq[upp]*upp_multx)*dE;
//			}
//        }
//	}

    if (Rad) {
        // ------------------------------------------
    	// Radiation-Induced: Bound-Bound Transitions
    	// ------------------------------------------
    	for(int bb=0;bb<crk.NBBlns;bb++)
		{
			na 	= crk.bbrad[bb].low[0];
			n 	= crk.bbrad[bb].low[1];
			ia 	= crk.bbrad[bb].upp[0];
			m 	= crk.bbrad[bb].upp[1];
			dE	= crk.bbrad[bb].dE;

			gpl = crk.gstate[na][n];
			gpu = crk.gstate[ia][m];

			if (na == ia && gpl == gpu) continue;

			low = crk.state[na][gpl]; 	// Find corresponding group for lower level
			upp = crk.state[ia][gpu]; 	// Find corresponding group for upper level

			grat = (double)crk.atom[ia]->ge[m]/(double)crk.atom[na]->ge[n];

			// Photoexcitation (n -> m)
			double B_nm = grat*MKS_hc*MKS_hc*crk.rad_Abb[bb]/(2*dE*dE*dE);
			// Stimulated emission (m -> n)
			double B_mns = B_nm/grat;
			// Spontaneous Emission (m -> n)
			double A_mn = crk.rad_Abb[bb];

			double low_multx = 1.;
			double upp_multx = 1.;
			if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
			{
				if (crk.group[na][gpl].grplvl)
				{
					double TgV = crk.group[na][gpl].Tgroup;
					low_multx = crk.boltz_weight(na,gpl,n,TgV);
				}
				if (crk.group[na][gpu].grplvl)
				{
					double TgV = crk.group[na][gpu].Tgroup;
					upp_multx = crk.boltz_weight(na,gpu,m,TgV);
				}
			}
			else
			{
				low_multx = (double)crk.atom[na]->ge[n]/(double)crk.group[na][gpl].ge;
				upp_multx = (double)crk.atom[ia]->ge[m]/(double)crk.group[ia][gpu].ge;
			}

			if (Irrad) {
				// Photoexcitation (n -> m)
				Jcr[low]      += (B_nm*low_multx*get_PlanckField(dE,TrV));
				Kcr[upp][low] += (B_nm*low_multx*get_PlanckField(dE,TrV));
				Jac[upp][low] += (B_nm*low_multx*get_PlanckField(dE,TrV));
				Jac[low][low] -= (B_nm*low_multx*get_PlanckField(dE,TrV));
				// Stimulated emission (m -> n)
				Jcr[upp]      += (B_mns*upp_multx*get_PlanckField(dE,TrV));
				Kcr[low][upp] += (B_mns*upp_multx*get_PlanckField(dE,TrV));
				Jac[upp][upp] -= (B_mns*upp_multx*get_PlanckField(dE,TrV));
				Jac[low][upp] += (B_mns*upp_multx*get_PlanckField(dE,TrV));
			}
			// Spontaneous emission (m -> n)
			Jcr[upp]      += (A_mn*upp_multx);
			Kcr[low][upp] += (A_mn*upp_multx);
			Jac[upp][upp] -= (A_mn*upp_multx);
			Jac[low][upp] += (A_mn*upp_multx);

		}

        // -----------------------------------------
    	// Radiation-Induced: Bound-Free Transitions
    	// -----------------------------------------
        for(int bf=0;bf<crk.NBFpIR;bf++)
        {
      		na 	= crk.bfpir[bf].low[0];
      		n 	= crk.bfpir[bf].low[1];
      		ia 	= crk.bfpir[bf].upp[0];
      		m 	= crk.bfpir[bf].upp[1];
      		dE	= crk.bfpir[bf].dE;

      		gpl = crk.gstate[na][n];
      		gpu = crk.gstate[ia][m];

      		if (na == ia && gpl == gpu) continue;

    		low = crk.state[na][gpl]; 	// Find corresponding group for lower level
    		upp = crk.state[ia][gpu]; 	// Find corresponding group for upper level

			double a_nm, a_nmE, s_nm, s_nmE, s_nm_der, s_nmE_der;
			a_nm = a_nmE = s_nm = s_nmE = s_nm_der = s_nmE_der = 0.;
			double dlnTe = ((crk.tab_lnT[it+1]) - (crk.tab_lnT[it]));
			if (Irrad) {
				// Photoionization (n -> m)
				a_nm = crk.tab_kpI[bf][0];
				a_nmE = crk.tab_kpI[bf][1];
				// Stimulated recombination (m -> n)
				s_nm = (1.-xt) * crk.tab_kpS[bf][it][0] + xt * crk.tab_kpS[bf][it+1][0];
				s_nmE = (1.-xt) * crk.tab_kpS[bf][it][1] + xt * crk.tab_kpS[bf][it+1][1];
				s_nm_der = (crk.tab_kpS[bf][it+1][0] - crk.tab_kpS[bf][it][0]) / (dlnTe * TeV);
				s_nmE_der = (crk.tab_kpS[bf][it+1][1] - crk.tab_kpS[bf][it][1]) / (dlnTe * TeV);
			}
			// Radiative Recombination (m -> n)
			double b_nm = (1.-xt) * crk.tab_kpR[bf][it][0] + xt * crk.tab_kpR[bf][it+1][0];
			double b_nmE = (1.-xt) * crk.tab_kpR[bf][it][1] + xt * crk.tab_kpR[bf][it+1][1];
			double b_nm_der = (crk.tab_kpR[bf][it+1][0] - crk.tab_kpR[bf][it][0]) / (dlnTe * TeV);
			double b_nmE_der = (crk.tab_kpR[bf][it+1][1] - crk.tab_kpR[bf][it][1]) / (dlnTe * TeV);

			double low_multx = 1.;
			double upp_multx = 1.;
			if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
			{
				if (crk.group[na][gpl].grplvl)
				{
					double TgV = crk.group[na][gpl].Tgroup;
					low_multx = crk.boltz_weight(na,gpl,n,TgV);
				}
				if (crk.group[ia][gpu].grplvl)
				{
					double TgV = crk.group[ia][gpu].Tgroup;
					upp_multx = crk.boltz_weight(ia,gpu,m,TgV);
				}
			}
			else
			{
				low_multx = (double)crk.atom[na]->ge[n]/(double)crk.group[na][gpl].ge;
				upp_multx = (double)crk.atom[ia]->ge[m]/(double)crk.group[ia][gpu].ge;
			}

			if (Irrad) {
				// Photoionization (n -> m)
				Jcr[low]      += (a_nm*low_multx);
				Kcr[upp][low] += (a_nm*low_multx);
				Jac[low][low] -= (a_nm*low_multx);
				Jac[upp][low] += (a_nm*low_multx);
				// Stimulated Recombination (m -> n)
				Jcr[upp]      += (s_nm*Nelec*upp_multx);
				Kcr[low][upp] += (s_nm*Nelec*upp_multx);
				Jac[low][upp] += (s_nm*Nelec*upp_multx);
				Jac[upp][upp] -= (s_nm*Nelec*upp_multx);
			}
			// Radiative Recombination (m -> n)
			Jcr[upp]      += (b_nm*Nelec*upp_multx);
			Kcr[low][upp] += (b_nm*Nelec*upp_multx);
			Jac[low][upp] += (b_nm*Nelec*upp_multx);
			Jac[upp][upp] -= (b_nm*Nelec*upp_multx);

			if (EnerEe) {
				if (Irrad) {
					// Photoionization (n -> m)
					Kcr[nee][low] += (a_nmE*low_multx)*TrV;
					Jac[nee][low] += (a_nmE*low_multx)*TrV;
					// Stimulated Recombination (m -> n)
					Kcr[nee][upp] -= (s_nmE*Nelec*upp_multx)*TeV;
					Jac[upp][nee] -= (s_nm_der*Nq[upp]*upp_multx*2./3.); // New term
					Jac[low][nee] += (s_nm_der*Nq[upp]*upp_multx*2./3.); // New term
					Jac[nee][upp] -= (s_nmE*Nelec*upp_multx)*TeV;
					Jac[nee][nee] -= (s_nmE_der*Nq[upp]*upp_multx*2./3.)*TeV;
				}
				// Radiative Recombination (m -> n)
				Kcr[nee][upp] -= (b_nmE*Nelec*upp_multx)*TeV;
				Jac[upp][nee] -= (b_nm_der*Nq[upp]*upp_multx*2./3.); // New term
				Jac[low][nee] += (b_nm_der*Nq[upp]*upp_multx*2./3.); // New term
				Jac[nee][upp] -= (b_nmE*Nelec*upp_multx)*TeV;
				Jac[nee][nee] -= (b_nmE_der*Nq[upp]*upp_multx*2./3.)*TeV;
			}

			if (DensNe) {
				if (Irrad) {
					// Photoionization (n -> m)
					Kcr[nev][low] += (a_nm*low_multx);
					Jac[nev][low] += (a_nm*low_multx);
					// Stimulated Recombination (m -> n)
					Jcr[nev]      += (s_nm*Nq[upp]*upp_multx);
					Jac[low][nev] += (s_nm*Nq[upp]*upp_multx);
					Jac[upp][nev] -= (s_nm*Nq[upp]*upp_multx);
					Jac[nev][upp] -= (s_nm*Nelec*upp_multx);
					Jac[nev][nev] -= (s_nm*Nq[upp]*upp_multx);
				}
				// Radiative Recombination (m -> n)
				Jcr[nev]      += (b_nm*Nq[upp]*upp_multx);
				Jac[low][nev] += (b_nm*Nq[upp]*upp_multx);
				Jac[upp][nev] -= (b_nm*Nq[upp]*upp_multx);
				Jac[nev][upp] -= (b_nm*Nelec*upp_multx);
				Jac[nev][nev] -= (b_nm*Nq[upp]*upp_multx);

				if (EnerEe) {
					// Stimulated Recombination (m -> n)
					Jac[nev][nee] -= (s_nm_der*Nq[upp]*upp_multx*2./3.);
					Jac[nee][nev] -= (s_nmE*Nq[upp]*upp_multx)*TeV;
					// Radiative Recombination (m -> n)
					Jac[nev][nee] -= (b_nm_der*Nq[upp]*upp_multx*2./3.);
					Jac[nee][nev] -= (b_nmE*Nq[upp]*upp_multx)*TeV;
				}
			}
        }
    }

    // ----------------------------------
    // Account for Bremmstrahlung losses
    // ----------------------------------
//	for (int n=0;n<Nstate;n++)
//	{
//		double Z = (double) Zq[n];
//		double kBrem = ZeldBrem(Z,TeV);
//
//		if (EnerEe) {
//			Kcr[nee][n]	  -= kBrem*Nelec;
//			Jac[nee][n]   -= kBrem*Nelec;
//			Jac[nee][nee] -= kBrem*Nq[n]/(3*TeV);
//
//			if (DensNe)
//				Jac[nee][nev]   -= kBrem*Nq[n];
//		}
//
//		Efree += kBrem*Nq[n]*Nelec;
//	}

    int NCRrho = Nstate;
    if (DensNe) NCRrho++;
    // Combine terms on RHS
    for(n=0;n<NCRrho;n++)
    {
        Scr[n] = -Jcr[n]*Nq[n];
        for(m=0;m<NCRrho;m++)
            Scr[n] += Kcr[n][m]*Nq[m];
    }

    if (EnerEe)
    	for(n=0;n<Nstate;n++)
    		Scr[nee] += Kcr[nee][n]*Nq[n];

//    for(n=0;n<NCRdim;n++)
//	{
//		Scr[n] = -Jcr[n]*Nq[n];
//		for(m=0;m<NCRdim;m++)
//		{
//			Scr[n] += Kcr[n][m]*Nq[m];
//		}
//	}

    //printf("Jac[%d][%d] = %12.5e  :  Nelec=%12.5e\n",19,19,Jac[19][19],Nelec);
}

void CR_BKE::source_QSS(double *y)
{
	int n,m,na,ia;
	int low,upp;

	KRatesGrp& crk = *pRk;
	double* Nq = y;

	for(n=0;n<NCRdim;n++)
	{
		Jcr[n] = 0.0;
		dNc[n] = 0.0;
		Scr[n] = 0.0;
		for(m=0;m<NCRdim;m++) Kcr[n][m] = 0.;
		for(m=0;m<NCRdim;m++) Jac[n][m] = 0.;
	}

	// Find interpolation point in tables
	double lnTeV = max(crk.lnT_min,min(.99*crk.lnT_max,log(TeV) ) );
	float xt = (lnTeV-crk.lnT_min)/crk.dlnT;
	int it = (int) xt; xt -= it;
	double lnNe  = max(crk.lnN_min,min(.99*crk.lnN_max,log(Nelec) ) );
	float xn = (lnNe-crk.lnN_min)/crk.dlnN;
	int in = (int) xn; xn -= in;

	double** dat = cArray<double>(2,2);
	double xb[2] = { 0., 1. };
	double yb[2] = { 0., 1. };

	// add-up ionizations and recombinations
	for(int bf=0;bf<crk.NBFeIR_QSS;bf++)
	{
		na = crk.bfeir_QSS[bf].low[0]; // neutral atom index
		ia = crk.bfeir_QSS[bf].upp[0]; // ionized atom index
		low = na; // find corresponding state
		upp = ia; // find corresponding state

		// ionization (n -> m)
		dat[0][0] = crk.tab_keI_QSS[in][it][bf];
		dat[0][1] = crk.tab_keI_QSS[in][it+1][bf];
		dat[1][0] = crk.tab_keI_QSS[in+1][it][bf];
		dat[1][1] = crk.tab_keI_QSS[in+1][it+1][bf];
		double a_nm = MathX::interpolate2d<double>(xb,yb,xn,xt,dat);

		// recombination (m -> n)
		dat[0][0] = crk.tab_keR_QSS[in][it][bf];
		dat[0][1] = crk.tab_keR_QSS[in][it+1][bf];
		dat[1][0] = crk.tab_keR_QSS[in+1][it][bf];
		dat[1][1] = crk.tab_keR_QSS[in+1][it+1][bf];
		double b_nm = MathX::interpolate2d<double>(xb,yb,xn,xt,dat);

		// Ionization (n -> m)
		Scr[low]      -= (a_nm*Nelec*Nq[low]);
        Scr[upp]	  += (a_nm*Nelec*Nq[low]);
		Jac[low][low] -= (a_nm*Nelec);
		Jac[upp][low] += (a_nm*Nelec);
        // Recombination (m -> n)
        Scr[upp]      -= (b_nm*Nelec*Nq[upp]);
        Scr[low]	  += (b_nm*Nelec*Nq[upp]);
		Jac[low][upp] += (b_nm*Nelec);
		Jac[upp][upp] -= (b_nm*Nelec);


//		if (isnan(Scr[low]) || isnan(Scr[upp])) printf("a_nm is %2.1e, b_nm is %2.1e, Nelec is %2.1e, Nlow is %2.1e, Nupp is %2.1e\n",a_nm,b_nm,Nelec,Nq[low],Nq[upp]);
//		if (isinf(Scr[low]) || isinf(Scr[upp])) printf("a_nm is %2.1e, b_nm is %2.1e, Nelec is %2.1e, Nlow is %2.1e, Nupp is %2.1e\n",a_nm,b_nm,Nelec,Nq[low],Nq[upp]);
//
//		if (isnan(Scr[low]) || isnan(Scr[upp]) || isinf(Scr[low]) || isinf(Scr[upp])) exit(1);
	}
}

// Conducts a 1st order Implicit Euler solve of the system
void CR_BKE::solveBKE(double* Nbase, double* Nstage, double dt)
{
    int n,m;

    double* Ncr = Nbase;
    double* Nst = Nstage;

    //rescale by actual time-step, add diagonal to Jacobian
    for(n=0;n<NCRdim;n++)
    {
        Scr[n] *= dt;
        for(m=0;m<NCRdim;m++) Jac[n][m] *= (-dt);
        Jac[n][n] += 1.;
    }
    MathX::Gauss_Jordan_point(NCRdim, dNc, Jac, Scr);

    for(int n=0;n<NCRdim;n++) Nst[n] = Ncr[n] + dNc[n];
}

void CR_BKE::resolveStep()
{
	KRatesGrp& crk = *pRk;

	recalcNeTe(Nc);

	source(Nc);
	solveBKE(Nc,Nstar,dt/2.);
	extractN_Grp(Nstar);
	recalcNeTe(Nstar);
	if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
		boltz_temp(Nstar);

	// ==============

	source(Nstar);
	solveBKE(Nstar,Nstar,dt/2.);
	extractN_Grp(Nstar);
	recalcNeTe(Nstar);
	if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
		boltz_temp(Nstar);
}

void CR_BKE::candidateStep()
{
	KRatesGrp& crk = *pRk;

	recalcNeTe(Nc);
	source(Nc);
	solveBKE(Nc,Ncand,dt);
	extractN_Grp(Ncand);
	recalcNeTe(Ncand);
	if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
		boltz_temp(Ncand);
}

void CR_BKE::setTolerances(double reltol, double abstol)
{
	rel_tol = reltol;
	abs_tol = abstol;
}

void CR_BKE::Integrate()
{
	KRatesGrp& crk = *pRk;

	int steps = 0;
	inst_fullN();
	if (crk.GrpScheme == "QSS") setGnds_QSS();
	else						extractN_Grp(Nc);
	calcERates(NcFull);
	SolutionOutput();

	printf("Integrating:\nStep %i at time = %3.3e\n",steps,tcurr);
	bool maxed = false;
	while (tcurr < tend)
	{
		if (crk.GrpScheme == "QSS")
		{
			source_QSS(Nc);
//			for (int i=0;i<NCRdim;i++)	printf("%i is %2.1e\n",i,Nc[i]);
//			for (int i=0;i<NCRdim;i++)	printf("%i is %2.1e\n",i,Scr[i]);
//			for (int i=0;i<NCRdim;i++)	{{for (int j=0;j<NCRdim;j++)  printf("%i,%i is %2.1e\t",i,j,Jac[i][j]);} printf("\n");}
//			exit(1);
			//calcERates(NcFull);
			solveBKE(Nc,Nc,dt);
			setGnds_QSS();
			extractN_QSS(NcFull);
			recalcNeTe(Nc);
		}
		else
		{
			source(Nc);
			//calcERates(NcFull);
			solveBKE(Nc,Nc,dt);
			extractN_Grp(Nc);
			recalcNeTe(Nc);
			if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
				boltz_temp(Nc);
		}

//		int ier = 0;
//		while (ier == 0)
//		{
//			// Must set NcFull for QSS adaptive timestepping
//			resolveStep();
//			candidateStep();	// must be after resolve to track energy rate transfers
//
//			exit(1);
//
//			double error = ErrorEstimate();
//			if (error >= 1. && !maxed)
//				dt *= sqrt(0.9/error);
//			else
//			{
//				for (int n=0;n<NCRdim;n++)
//					Nc[n] = Ncand[n];
//
//				ier = 1;
//			}
//			printf("\tError: %3.3e, dt: %3.3e\n",error,dt);
//		}

		steps++;
		tcurr += dt;
		printf("Step %i at time = %3.3e with dt = %3.3e ",steps,tcurr,dt);

		calcERates(NcFull);
		SolutionOutput();
		dt *= multxdt;
//		if (dt > 1.e-6)
//		{
//			maxed = true;
//			dt = 1.e-6;
//		}

		if (tcurr + dt > tend) dt = tend - tcurr;

		monitorCollNe(dt);

		printf("\n");
	}

//	for(int n=0;n<crk.KRates::Nstate;n++) {printf("NcFull %i: %e; ", n, NcFull[n]);
//	printf("\n");}
}


void CR_BKE::newIntegrate(double *qSrc)
{
    KRatesGrp& crk = *pRk;

    int steps = 0;
    inst_fullN();
    if (crk.GrpScheme == "QSS") setGnds_QSS();
    else						extractN_Grp(Nc);
    calcERates(NcFull);
    SolutionOutput();

    printf("Integrating:\nStep %i at time = %3.3e\n",steps,tcurr);
    bool maxed = false;
    while (tcurr < tend)
    {
        if (crk.GrpScheme == "QSS")
        {
            source_QSS(Nc);
//			for (int i=0;i<NCRdim;i++)	printf("%i is %2.1e\n",i,Nc[i]);
//			for (int i=0;i<NCRdim;i++)	printf("%i is %2.1e\n",i,Scr[i]);
//			for (int i=0;i<NCRdim;i++)	{{for (int j=0;j<NCRdim;j++)  printf("%i,%i is %2.1e\t",i,j,Jac[i][j]);} printf("\n");}
//			exit(1);
            //calcERates(NcFull);
            solveBKE(Nc,Nc,dt);
            setGnds_QSS();
            extractN_QSS(NcFull);
            recalcNeTe(Nc);
        }
        else
        {
            source(Nc);
            //calcERates(NcFull);
            solveBKE(Nc,Nc,dt);
            extractN_Grp(Nc);
            recalcNeTe(Nc);
            if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
                boltz_temp(Nc);
        }

//		int ier = 0;
//		while (ier == 0)
//		{
//			// Must set NcFull for QSS adaptive timestepping
//			resolveStep();
//			candidateStep();	// must be after resolve to track energy rate transfers
//
//			exit(1);
//
//			double error = ErrorEstimate();
//			if (error >= 1. && !maxed)
//				dt *= sqrt(0.9/error);
//			else
//			{
//				for (int n=0;n<NCRdim;n++)
//					Nc[n] = Ncand[n];
//
//				ier = 1;
//			}
//			printf("\tError: %3.3e, dt: %3.3e\n",error,dt);
//		}

        steps++;
        tcurr += dt;
//        printf("Step %i at time = %3.3e with dt = %3.3e ",steps,tcurr,dt);

        calcERates(NcFull);
        SolutionOutput();
        dt *= multxdt;
//		if (dt > 1.e-6)
//		{
//			maxed = true;
//			dt = 1.e-6;
//		}

        if (tcurr + dt > tend) dt = tend - tcurr;

        monitorCollNe(dt);

//        printf("\n");
    }

    for(unsigned k = 0; k<NCRdim; k++)
        qSrc[k] = dNc[k];

//	for(int n=0;n<crk.KRates::Nstate;n++) {printf("NcFull %i: %e; ", n, NcFull[n]);
//	printf("\n");}
}

double CR_BKE::ErrorEstimate()
{
	double L2norm = 0.;
	for (int n=0;n<NCRdim;n++)
	{
		double diff = Ncand[n] - Nstar[n];
		double scerr = Ncand[n]*rel_tol + abs_tol;
		L2norm += diff*diff/(scerr*scerr);
	}

	double err = sqrt(L2norm/NCRdim);
	return err;
}

// Function that controls the output of the results.
// Modify this routine according to your needs
void CR_BKE::SolutionOutput()
{
	KRatesGrp& crk = *pRk;
//	Spectrum& spt = *pSpect;
	int tr,m,l,na,ia,gpu,upp;

	cout << setiosflags(ios::showpoint);// | ios::fixed);

	// ---------------------------------------------------------- //
	// ----------- Output to File (Plasma Properties) ----------- //
	// ---------------------------------------------------------- //

	fprintf(fp,"%3.3e",tcurr);
	fprintf(fo,"%3.3e\t%3.3e",tcurr,dt);

	for (l=0;l<Nstate;l++)
		fprintf(fp,"\t%5.4e",Nc[l]);

	double Nhtot,Nztot;
	Nhtot = Nztot = 0.;
//	for (l=0;l<Nstate;l++)
//	{
//		Nhtot += Nc[l];
//		Nztot += (double)Zq[l]*Nc[l];
//	}

	for (int na=0;na<Natoms;na++)
	{
		AtomX& A = *(pRk->atom[na]);
		double Zatm = (double)A.Z;
		for(int m=0;m<crk.KRates::Nalevs[na];m++)
		{
			int low = crk.KRates::state[na][m];
			Nhtot += NcFull[low];
			Nztot += Zatm*NcFull[low];
		}
	}
	fprintf(fp,"\t%5.4e",Nelec);
	fprintf(fo,"\t%5.4e\t%3.3e\t%3.3e",Nelec,Nztot/Nhtot,TeV);

	if (DensNe)
		fprintf(fo,"\t%5.4e",Scr[nev]);
	else
		fprintf(fo,"\t%5.4e",0.);

	if (EnerEe)
		fprintf(fo,"\t%5.4e",Scr[nee]);
	else
		fprintf(fo,"\t%5.4e",0.);

	fprintf(fp,"\n");
	fflush(fp);

	Etot = (Efree+Erecp+Edexp);
	Eint += (Efree+Erecp+Edexp)*dt;
	fprintf(fo,"\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n", \
			Eexce,Edexe,Eione,Erece,Eexcp,Eexcs,Edexp,Eionp,Eions,Erecp,Efree,Etot,Eint);
	fflush(fo);

	// ------------------------------------------------------- //
	// ----------- Output to File (Line Radiation) ----------- //
	// ------------------------------------------------------- //

		// This must be outputted before the "Boltzmann Profiles" section

	if (tcurr >= tprint)
	{
		for(na=0;na<crk.Natoms;na++)
		{
			char filename[100];
			std::ofstream fileBolt;   // file stream
			sprintf(filename,crk.outputDir);
			sprintf(filename + strlen(filename),"LineRad/Ar%02d/%dLineRad%08d.txt",(na+1),(na+1),filenum);
			fileBolt.open (filename);  // open the file in filenames at index i
			fileBolt << "Energy,\tAr" << na+1 << ",\tlow\t,upp\n";

			for(int bb=0;bb<crk.KRates::NBBlns;bb++)
			{
				int atm 	= crk.bbrad[bb].low[0];
					if (atm != na) continue;
				int lo 		= crk.bbrad[bb].low[1];
				int up  	= crk.bbrad[bb].upp[1];
				double dE	= crk.bbrad[bb].dE;

				fileBolt << dE << ",\t" << RadSEbb[bb] << ",\t" << lo << ",\t" << up << "\n";
			}

			fileBolt.close();
		}
	}

	fprintf(fRad,"%3.3e",tcurr);
	for(na=0;na<crk.Natoms;na++)
		fprintf(fRad,"\t%5.4e",RadSEZ[na]);
	fprintf(fRad,"\n");

	// -------------------------------------- //
	// --- Output to File (Spectral Data) --- //
	// -------------------------------------- //

	if (SpecOut && tcurr >= tprint)
	{
//		spt.resetEmisBins();
//		spt.apply_BBemission(& crk,Nelec,RadSEbb);
////		spt.apply_BFemission(& crk,NcFull,Nelec,TeV);
////		spt.apply_FFemission(& crk,NcFull,Zq,Nelec,TeV);
//		spt.printSpec(crk.outputDir,filenum);
	}

	// ---------------------------------------------------------- //
	// --- Output to File (Level densities and Charge States) --- //
	// ---------------------------------------------------------- //

	fprintf(fp_full,"%3.3e",tcurr);
	fprintf(fChrg,"%3.3e",tcurr);
	for(na=0;na<crk.Natoms;na++)
	{
		double Na = 0.;
		double Nexc = 0.;
		for(m=0;m<crk.KRates::Nalevs[na];m++)
		{
			int ka 		= crk.KRates::state[na][m];
			Na += NcFull[ka];
			if (m > 0) Nexc += NcFull[ka];
			fprintf(fp_full,"\t %5.4e",NcFull[ka]);
		}
		fprintf(fChrg,"\t%5.4e\t%5.4e",Na,Nexc);
	}

	fprintf(fChrg,"\t%5.4e\n",Nhtot);
	fprintf(fp_full,"\t%3.3e\n",Nelec); fflush(fp_full);

	// ------------------------------------------- //
	// --- Output to File (Boltzmann Profiles) --- //
	// ------------------------------------------- //

	if (tcurr >= tprint)
	{
		fprintf(ftBolt,"%i\t%3.3e\n",filenum,tcurr); fflush(ftBolt);
		for(na=0;na<crk.Natoms;na++)
		{
			char filename[150];
			std::ofstream fileBolt;   // file stream
			sprintf(filename,crk.outputDir);
			sprintf(filename + strlen(filename),"BoltzDist/Ar%02d/%dBoltzPlot%08d.txt",(na+1),(na+1),filenum);
			fileBolt.open (filename);  // open the file in filenames at index i
			fileBolt << "Energy,\tAr" << na+1 << "\n";

			for(m=0;m<crk.KRates::Nalevs[na];m++)
			{
				int n		 = crk.KRates::state[na][m];
				double Ediff = (double)crk.atom[na]->Ee[m];
				double ge = (double)crk.atom[na]->ge[m];
				double logNg = log10(NcFull[n]/ge);

				fileBolt << Ediff << ",\t" << logNg << "\n";
			}

			fileBolt.close();
		}
		filenum++;

		// Increment printing dt until appropriately out of print step range
		double tlolim = tcurr + dt;
		if (printlog)
			while (tprint < tlolim)
			{
				logtprt += logdtprt;
				tprint = pow(10,logtprt);
			}
		else
			while (tprint < tlolim)
				tprint += dtprint;
	}

	// ------------------------------------------- //
	// --- Output to File (Group Temperatures) --- //
	// ------------------------------------------- //

	if (crk.GrpScheme == "Boltzmann" || crk.GrpScheme == "Custom_Boltzmann")
	{
		fprintf(fTg,"%3.3e",tcurr);
		for(na=0;na<crk.Natoms;na++)
			for(m=0;m<crk.Ngroups[na];m++)
				if (crk.group[na][m].baselvl)
				{
//					int mplus = m+1;
//					int low	 = crk.state[na][m];
//					int grp	 = crk.state[na][mplus];
//
//					double Ee = (double)crk.group[na][mplus].Ee - (double)crk.group[na][m].Ee;
//					double ge = (double)crk.group[na][m].ge;
//					double ggrp = (double)crk.group[na][mplus].ge;
//
//					double Tg = -Ee/log(Nc[grp]*ge/(Nc[low]*ggrp));
//
//					fprintf(fTg,"\t(%4.3e\t%4.3e)",Tg,crk.group[na][m].Tgroup);

					int piv = crk.state[na][m];
					int grp = crk.state[na][m+1];

					if (Nc[piv] <= 1.e16 || Nc[grp] <= 1.e16)
						fprintf(fTg,"\t%4.3e",numeric_limits<double>::quiet_NaN());
//					if (crk.group[na][m].Tgroup > 100.*TeV)
//						fprintf(fTg,"\t%4.3e",numeric_limits<double>::quiet_NaN());
					else
						fprintf(fTg,"\t%4.3e",crk.group[na][m].Tgroup);
				}

		fprintf(fTg,"\n");
	}

	fflush(fTg);


}  // SolutionOutput

