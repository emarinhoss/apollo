/*
 * KRatesGrp.h
 *
 *  Created on: Jun 7, 2016
 *      Author: richard
 */

#ifndef KRATESGRP_H_
#define KRATESGRP_H_

#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "KRates.h"

using namespace std;


/*
 * Isolated transitions for group-to-group transitions
 */
struct transition_iso
{
	int glotridx;	// transition index from full list
	int low[4]; 	// 0: atom index, 1: level index, 2: charge for lower level, 3: group index
	int upp[4]; 	// 0: atom index, 1: level index, 2: charge for upper level, 3: group index
	float dE;   	// energy threshold or ionization potential
};

/*
 * Atomic level-grouping details
 */
struct atomgroup
{
	bool	baselvl;	// True indicates Boltzmann group is base/ground level
	bool	grplvl;		// True indicates Boltzmann group is upper levels
	bool	boltrev;	// Revert to explicit levels due to low level count for group
	int		identT;
	int		Nlvlgrp;	// Number of levels within group
	int		Z;			// Group charge
	int		ge;			// Group-averaged degeneracy
	int*	paridx;		// Particular index of level relative to ionic listing
	int		nS;			// Number of shells associated with the group
    int*	qn;			// Principal quantum number	 ------ first idx: shell
    int*	ql;			// Angular quantum number	 ------ first idx: shell
    int*	qw;			// Orbital Occupation Number ------ first idx: shell
    float	Ee;			// Group-average energy [eV]

    double	Tgroup;
};


class KRatesGrp : public KRates
{
private:
    int nq;				// Instantiate state index

public:
	/* Name of grouping scheme */
	string GrpScheme;	// 1) Full			: Use KRates to conduct a full solve
							// 2) QSS			: Assumes excited states are always at steady state
							// 3) Uniform		: Conserves number density within the group
							// 4) Boltzmann		: Conserves number density & energy within the group

	/* Level Grouping */
	int		Natoms;		// Number of atoms considered
	int 	Nstate;		// Total number of groups in reduced solver
	int 	Nlvful;		// Total number of levels from full solver
    int* 	Ngroups;	// Number of groups considered for atom
    int**	gstate;		// Indexing between "full" atomic levels and grouped atomic levels
    int** 	state;		// Indexing between grouped atomic levels and complexity-reduced solver listing
    atomgroup** group;	// first idx: atom index; second idx: group index
    //--- QSS section
    int 	Natoms_QSS;
	int 	NNstep;
	double 	lnN_min, lnN_max, dlnN;
	double* tab_lnN;

    /* Rate Grouping */
    int NBBlniso;		// Number of bound-bound transitions isolated from full for grouping
    int NBBeXiso;		// Number of electron-impact excitation transitions isolated from full for grouping
    int NBFaIiso;		// Number of autoionization transitions isolated from full for grouping
    int NBFeIiso;		// Number of electron-impact ionization transitions isolated from full for grouping
    int NBFpIiso;		// Number of electron-impact ionization transitions isolated from full for grouping
	transition_iso* bbexdiso;
	transition_iso* bbradiso;
	transition_iso* bfairiso;
	transition_iso* bfeiriso;
	transition_iso* bfpiriso;
	//--- QSS section
	int Nstate_QSS;
	int NBFeIR_QSS;
	transition* bfeir_QSS;
	double***  tab_keI_QSS;
	double***  tab_keR_QSS;
	double***  tab_Lz; // plasma radiative cooling rate
	double***  tab_Kz; // electron heating/cooling rate
	double****  tab_r0;
	double****  tab_rp;

	//KRatesGrp();
    KRatesGrp(const char* spcsName=NULL,int Na=0);
    virtual ~KRatesGrp();

    void AtomDist(double,double,double);
    void setGroupRange(int,float,float);

    void checkGroup();
    void unifyAtoms();
    void printGroup();

    void DistPrep(string);
    int full_copy();
    int QSS_group();
    int unif_shell_group();
    int bolt_shell_group();
    int unif_cust_group();
    int bolt_cust_group();

    void add_BBtrans(bool); // Excitation/De-excitations
    void full_BB_copy();
    void unif_BB_group();
    void bolt_BB_group();
    void printBBRates();

	void add_BFtrans(bool); // Ionization/Recombinations
	void full_BF_copy();
	void unif_BF_group();
	void bolt_BF_group();

	void printBFRates();

	void xfer_ECexchg(); // Elastic Collisions

	double boltz_weight(int,int,int,double);

	//---QSS section
	void setNRange(int,float,float);
	void countQSSAtom();
	void add_QSS_BFtrans();
	void compute_qss_rates(int, double, double, double*, double*, double*, double*, double&, double&, double&, double&);

	double part_funct(int,int,double);
	double part_funct(int,int,int,double);
	double part_funct_renorm(int,int,int,float);
	double part_funct_E1mom(int,int,int,float);
	double part_funct_E2mom(int,int,int,float);
	double mom1_Egrp(int,int,float,bool);
	double mom2_Egrp(int,int,float,bool);
	double bolt_Egap(int,int,int,float,float,float,float);
	double bolt_Igap(int,int,int,int,float,float,float,float,float);

};

#endif /* KRATESGRP_H_ */
