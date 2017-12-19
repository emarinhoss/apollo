/*
 * CRSolve.h
 *
 *  Created on: Jan 6, 2016
 *      Author: richard
 */

#ifndef CRSOLVE_H_
#define CRSOLVE_H_

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

#include "Constants.h"
#include "KRates.h"
#include "KRatesGrp.h"
//#include "LSpecs.h"
//#include "Spectrum.h"
//#include "StiffIntegratorT.h"

class CR_BKE
{
private:

public:
	int nq;			// State count (for QSS)
    int NTstep;		// Number of temperature bins
    int Natoms;		// Number of atoms considered
    int Nstate;		// Number of atomic variables for current solver
    int NCRdim;		// Size of Numerical Solve (>= Nstate)
    int Ngroups;	// Number of atomic groupings
    int nev;		// Variable for electron density
    int neh;		// Variable for heavy species energy
    int nee;		// Variable for electron energy
    int Nexlvls;	// Number of excited levels (Used for QSS Solver)
    int filenum;

    bool Rad;		// Active Radiation
    bool Irrad;		// Active Radiation
    bool DensNe;	// Active Electron Density Eqn.
    bool FxdNe;		// Fixed Electron Density (DensNe = false)
    bool EnerEh;	// Active Heavy Species Energy Eqn.
    bool EnerEe;	// Active Electron Energy Eqn.
    bool SpecOut;	// Active Spectral Calculations
    bool printlog;	// Boltzmann plot output in log steps

    double	Nelec;
    double  TeV;	// Electron Temperature
    double  TrV;	// Radiation Temperature
    double	Eexce,Edexe;
    double	Eione,Erece;
    double	Eexcp,Eexcs,Edexp;
    double	Eionp,Eions,Erecp;
    double  Efree;
    double  Etot;
    double 	Eint;
    double	tbeg,tend;
	double	tcurr;
	double	dt;		// Initial timestep
	double	multxdt;
	double	tprint;	// Boltzmann plot output variables
	double	dtprint;
	double	logtprt;
	double	logdtprt;
	double  rel_tol;
	double  abs_tol;

    int* Zq;

    int** state_QSS;

    double*  RadSEZ;
    double*	 RadSEbb;
    double*  Nstar;
    double*  Ncand;
    double*	 Nc;
    double*  Scr;
    double*  Jcr;
    double*  dNc;
    double*  NcFull;

    double** Kcr;
    double** Jac;

    FILE* fp;
	FILE* fo;
	FILE* fTg;
	FILE* fRad;
	FILE* fChrg;
	FILE* fp_full;
	FILE* ftBolt;

    KRatesGrp* pRk;

//    Spectrum* pSpect;

    CR_BKE(int,bool,bool,bool,bool);
    virtual ~CR_BKE();

    void setRad(bool);
    void setIrrad(bool);
    void addRates(KRatesGrp*);
//    void addSpecFeat(Spectrum*);
    void setTemp(double,double);
    void setsol(double*,double);
    void getsol(double*);
    void setTimes(double,double,double);
    void setPrintVars(bool,double,double);
    void setGnds_QSS();
    void setExcs_QSS(double*);
    void recalcNeTe(double*);
    void monitorCollNe(double);
    void calcERates(double*);
    virtual void source(double*);
    virtual void source_QSS(double*);
    virtual void extractN_QSS(double*);
    virtual void solveBKE(double*,double*,double);
    double tstep();
    double stiff(double*);

    void initRecord();
    double get_PlanckField(double,double);
    void inst_fullN();
    void inst_Tgroup();
    void extractN_Grp(double*);
    void boltz_temp(double*);
    double iter_Tgroup(int,int,double,double*);
    double calc_xiEnergy(int,int,int,double,bool);
    double calc_deltEE(int,int,int,int,double*,double);
    double calc_deltIE(int,int,int,int,double*,double);

    void SolveSS(double);
    void invert_SS(double);

    void Integrate();
    void resolveStep();
    void candidateStep();
    void setTolerances(double,double);
    double ErrorEstimate();
    void SolutionOutput();
};


#endif /* CRSOLVE_H_ */
