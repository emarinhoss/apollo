/*
 * KRatesGrp.cpp
 *
 *  Created on: Jun 7, 2016
 *      Author: richard
 */

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <string.h>
#include <gsl/gsl_linalg.h>

#include "cArray.h"
#include "KRatesGrp.h"

/**
 * Constructor
 */
KRatesGrp::KRatesGrp(const char* spcsName, int Na) : KRates(spcsName,Na)
{
	nq 	= 0;

	NBBlniso = 0;
	NBBeXiso = 0;
	NBFaIiso = 0;
	NBFeIiso = 0;
	NBFpIiso = 0;

	Natoms = Na;
	Nstate = 0;
	Nlvful = 0;
	group = new atomgroup*[Na];
	Ngroups = new int[Na];
	gstate = new int*[Na];
	state = new int*[Na];

	bbexdiso=NULL;
	bbradiso=NULL;
	bfairiso=NULL;
	bfeiriso=NULL;
	bfpiriso=NULL;

	// ------- QSS
	Nstate_QSS = 0;
	Natoms_QSS = 0;
	NNstep = 0;
	lnN_min = 0.;
	lnN_max = 0.;
	dlnN = 0.;
	tab_lnN=NULL;
	tab_lnN=NULL;
	tab_keI_QSS=NULL;
	tab_keR_QSS=NULL;
	tab_Lz=NULL;
	tab_Kz=NULL;
	tab_r0=NULL;
	tab_rp=NULL;
	NBFeIR_QSS = 0;
	bfeir_QSS=NULL;
}

/**
 * Destructor
 */
KRatesGrp::~KRatesGrp()
{
	for (int n=0;n<Natoms;n++) {
		for (int m=0;m<Ngroups[n];m++) {
			delete[] group[n][m].paridx;
			delete[] group[n][m].qn;
			delete[] group[n][m].ql;
			delete[] group[n][m].qw;
		}
	}

	for (int n=0;n<Natoms;n++) {
		delete[] group[n];
		delete[] gstate[n];
		delete[] state[n];
	}

	delete[] group;
	delete[] gstate;
	delete[] state;

	//if(tab_lnT) delete[] tab_lnT;
	delete [] Ngroups;

	delete [] bbexdiso;
	delete [] bbradiso;
	delete [] bfairiso;
	delete [] bfeiriso;
	delete [] bfpiriso;

	// QSS
	if(tab_lnN) delete[] tab_lnN;
	if(bfeir_QSS) delete [] bfeir_QSS;
	del_cArray<double>(NNstep+1,NTstep+1,Natoms_QSS,tab_keI_QSS);
	del_cArray<double>(NNstep+1,NTstep+1,Natoms_QSS,tab_keR_QSS);
	del_cArray<double>(NNstep+1,NTstep+1,Natoms_QSS,1,tab_r0);
	del_cArray<double>(NNstep+1,NTstep+1,Natoms_QSS,1,tab_rp);
	del_cArray<double>(NNstep+1,NTstep+1,Natoms,tab_Lz);
	del_cArray<double>(NNstep+1,NTstep+1,Natoms,tab_Kz);

}

// Add Atom for QSS case
void KRatesGrp::countQSSAtom()
{
    for (int na=0;na<Natoms;na++)
    {
    	if (atom[na]->ion) {
			Natoms_QSS++;
			printf("Config. QSS mode for atom %s \n",atom[na]->name.c_str());
		}
    	Nstate_QSS++;
    }

    printf("Number of QSS Atoms: %i\n",Natoms_QSS);
}

/*
 * Assign groups to specific levels
 * // Consider porting this function into AtomXGrp.cpp
 */
void KRatesGrp::DistPrep(string GrpType)
{
	printf("Preparing distribution bins based on %s grouping:\n",GrpType.c_str());

	int totStates = 0;

	GrpScheme = GrpType;

	for(int n=0;n<Natoms;n++)
		gstate[n] = new int[Nalevs[n]];

	if (GrpScheme == "Full" || GrpScheme == "SS")
		totStates = full_copy();
	else if (GrpScheme == "QSS")
		totStates = QSS_group();
	else if (GrpScheme == "Uniform")
		totStates = unif_shell_group();
	else if (GrpScheme == "Boltzmann")
		totStates = bolt_shell_group();
	else if (GrpScheme == "Custom_Uniform")
		totStates = unif_cust_group();
	else if (GrpScheme == "Custom_Boltzmann")
		totStates = bolt_cust_group();
	else {
		printf("Not a valid grouping option: %s\n",GrpScheme.c_str());
		exit(1);
	}

	checkGroup();
	printGroup();

	Nstate = totStates;
	Ncr_init = new double[totStates]();
}
// [Active]
int KRatesGrp::full_copy()
{
	int nStates = 0;

	for(int n=0;n<Natoms;n++)
	{
		// Descriptors for number of states per group, shell configurations
		Ngroups[n] = Nalevs[n];
		nStates += Nalevs[n];
		group[n] = new atomgroup[Nalevs[n]];

		for(int m=0;m<Nalevs[n];m++)
		{
			group[n][m].Nlvlgrp = 1;
			group[n][m].ge		= 0;
			group[n][m].Ee 		= 0;
			group[n][m].paridx 	= NULL;
			group[n][m].qn 		= NULL;
			group[n][m].ql 		= NULL;
			group[n][m].qw 		= NULL;
			group[n][m].baselvl = false;
			group[n][m].grplvl 	= false;
			group[n][m].boltrev = false;
			group[n][m].identT	= NULL;

			gstate[n][m] = m;
		}

		for(int m=0;m<Nalevs[n];m++)
			group[n][m].paridx = new int [1];

		for(int m=0;m<Nalevs[n];m++)
		{
			group[n][m].paridx[0] = m;
			group[n][m].Z = atom[n]->Z;
			group[n][m].ge += atom[n]->ge[m];
			group[n][m].Ee += atom[n]->Ee[m];
		}
	}

	return nStates;
}
// [Inactive] - Extended neglect
int KRatesGrp::unif_shell_group()
{
	int nStates = 0;

	// Determine the number of groups for the atom.
	for(int n=0;n<Natoms;n++)
	{
		// Initialize number of group count for atom 'n'
		Ngroups[n] = 0;

		int* qn_grp;
		int* ql_grp;
		int* qw_grp;

		qn_grp = (int*) malloc (atom[n]->nS[0]*sizeof(int));
		ql_grp = (int*) malloc (atom[n]->nS[0]*sizeof(int));
		qw_grp = (int*) malloc (atom[n]->nS[0]*sizeof(int));

		Ngroups[n]++;

		for(int m=0;m<Nalevs[n];m++) {
			bool diffShell = false;
			for (int ns=0; ns<atom[n]->nS[m]; ns++) {
				if (m == 0) {
					qn_grp[ns] = atom[n]->qn[m][ns];
					ql_grp[ns] = atom[n]->ql[m][ns];
					qw_grp[ns] = atom[n]->qw[m][ns];
				}
				else if (qn_grp[ns] == atom[n]->qn[m][ns] && \
						ql_grp[ns] == atom[n]->ql[m][ns] && \
						qw_grp[ns] == atom[n]->qw[m][ns])
					continue;
				else {
					diffShell = true;
					Ngroups[n]++;
					break;
				}
			}

			if (diffShell) {
				qn_grp = (int*) realloc (qn_grp,sizeof(int)*atom[n]->nS[m]);
				ql_grp = (int*) realloc (ql_grp,sizeof(int)*atom[n]->nS[m]);
				qw_grp = (int*) realloc (qw_grp,sizeof(int)*atom[n]->nS[m]);

				for (int ns=0; ns<atom[n]->nS[m]; ns++)
				{
					qn_grp[ns] = atom[n]->qn[m][ns];
					ql_grp[ns] = atom[n]->ql[m][ns];
					qw_grp[ns] = atom[n]->qw[m][ns];
				}
			}
		}

		free(qn_grp);
		free(ql_grp);
		free(qw_grp);

		nStates += Ngroups[n];
	}

	// Descriptors for number of states per group, shell configurations
	for(int n=0;n<Natoms;n++) {

		// Initialize atomic groups
		group[n] = new atomgroup[Ngroups[n]];
		int grpidx = 0;

		for(int m=0;m<Ngroups[n];m++)
		{
			group[n][m].Nlvlgrp = 0;
			group[n][m].ge		= 0;
			group[n][m].Ee 		= 0;
			group[n][m].paridx 	= NULL;
			group[n][m].qn 		= NULL;
			group[n][m].ql 		= NULL;
			group[n][m].qw 		= NULL;
		}

		// Determine group shell configurations
		for(int m=0;m<Nalevs[n];m++)
		{
			bool sameSt = true;
			for (int ns=0;ns<atom[n]->nS[m];ns++)
			{
				if (m == 0)
				{
					group[n][grpidx].Nlvlgrp++;
					sameSt = false;
					break;
				}
				else if (group[n][grpidx].qn[ns] == atom[n]->qn[m][ns] && \
						group[n][grpidx].ql[ns] == atom[n]->ql[m][ns] && \
						group[n][grpidx].qw[ns] == atom[n]->qw[m][ns])
					continue;
				else
				{
					sameSt = false;
					grpidx++;
					group[n][grpidx].nS = atom[n]->nS[m];
					group[n][grpidx].Nlvlgrp++;
					break;
				}
			}

			if (sameSt == true)
			{
				group[n][grpidx].Nlvlgrp++;
				group[n][grpidx].ge += atom[n]->ge[m];
			}
			else
			{
				group[n][grpidx].qn = new int [atom[n]->nS[m]];
				group[n][grpidx].ql = new int [atom[n]->nS[m]];
				group[n][grpidx].qw = new int [atom[n]->nS[m]];

				for (int ns=0; ns<atom[n]->nS[m]; ns++)
				{
					group[n][grpidx].qn[ns] = atom[n]->qn[m][ns];
					group[n][grpidx].ql[ns] = atom[n]->ql[m][ns];
					group[n][grpidx].qw[ns] = atom[n]->qw[m][ns];
				}

				group[n][grpidx].Z = atom[n]->Z;
				group[n][grpidx].ge += atom[n]->ge[m];
				group[n][grpidx].Ee += atom[n]->ge[m]*atom[n]->Ee[m];
			}
		}

		grpidx++;	// Increased due to 0 index
		if (grpidx == Ngroups[n])
			printf("Atom %s grouping shells (%i) have been determined...\n", \
					atom[n]->name.c_str(),Ngroups[n]);
		else
			printf("Grouping error detected\n");
	}
	printf("\n");

	for(int n=0;n<Natoms;n++)
		for(int m=0; m<Ngroups[n]; m++)
			group[n][m].Ee /= group[n][m].ge;

	// Descriptor for indexing between full levels list and levels grouped
	for(int n=0;n<Natoms;n++)
	{
		for(int m=0;m<Ngroups[n];m++)
			group[n][m].paridx = new int [group[n][m].Nlvlgrp];

		int grpidx = 0;
		int pidx = 0;

		// Match shell configurations
		for(int m=0;m<Nalevs[n];m++)
		{
			for (int ns=0; ns<atom[n]->nS[m]; ns++)
			{
				if (group[n][grpidx].qn[ns] == atom[n]->qn[m][ns] && \
						group[n][grpidx].ql[ns] == atom[n]->ql[m][ns] && \
						group[n][grpidx].qw[ns] == atom[n]->qw[m][ns])
					continue;
				else
				{
					pidx = 0;
					grpidx++;
					break;
				}
			}
			group[n][grpidx].paridx[pidx] = m;
			pidx++;
		}

		grpidx++;	// Increased due to 0 index
		if (grpidx == Ngroups[n])
			printf("Atom %s uniform grouping verified and applied...\n",atom[n]->name.c_str());
		else
			printf("Grouping error detected\n");
	}
	printf("\n");

	return nStates;
}
// [Inactive] - Extended neglect
int KRatesGrp::bolt_shell_group()
{
	int nStates = 0;

	// Determine the number of groups for the atom.
	for(int n=0;n<Natoms;n++)
	{
		// Initialize number of group count for atom 'n'
		Ngroups[n] = 0;

		int* qn_grp;
		int* ql_grp;
		int* qw_grp;

		qn_grp = (int*) malloc (atom[n]->nS[0]*sizeof(int));
		ql_grp = (int*) malloc (atom[n]->nS[0]*sizeof(int));
		qw_grp = (int*) malloc (atom[n]->nS[0]*sizeof(int));

		Ngroups[n]++;

		int numGrp = 0;
		for(int m=0;m<Nalevs[n];m++) {
			bool diffShell = false;
			for (int ns=0;ns<atom[n]->nS[m];ns++) {
				if (m == 0) {
					qn_grp[ns] = atom[n]->qn[m][ns];
					ql_grp[ns] = atom[n]->ql[m][ns];
					qw_grp[ns] = atom[n]->qw[m][ns];
				}
				else if (qn_grp[ns] == atom[n]->qn[m][ns] && \
						ql_grp[ns] == atom[n]->ql[m][ns] && \
						qw_grp[ns] == atom[n]->qw[m][ns])
					continue;
				else {
					diffShell = true;
					Ngroups[n]++;
					break;
				}
			}

			if (diffShell)
			{
				qn_grp = (int*) realloc (qn_grp,sizeof(int)*atom[n]->nS[m]);
				ql_grp = (int*) realloc (ql_grp,sizeof(int)*atom[n]->nS[m]);
				qw_grp = (int*) realloc (qw_grp,sizeof(int)*atom[n]->nS[m]);

				for (int ns=0; ns<atom[n]->nS[m]; ns++)
				{
					qn_grp[ns] = atom[n]->qn[m][ns];
					ql_grp[ns] = atom[n]->ql[m][ns];
					qw_grp[ns] = atom[n]->qw[m][ns];
				}

				numGrp = 1;
			}
			else
			{
				numGrp++;
				if (numGrp == 2)	Ngroups[n]++;
			}
		}

		free(qn_grp);
		free(ql_grp);
		free(qw_grp);

		nStates += Ngroups[n];
	}

	// Descriptors for number of states per group, shell configurations
	for(int n=0;n<Natoms;n++) {

		// Initialize atomic groups
		group[n] = new atomgroup[Ngroups[n]];
		int grpidx = 0;

		for(int m=0;m<Ngroups[n];m++)
		{
			group[n][m].Nlvlgrp = 0;
			group[n][m].ge		= 0;
			group[n][m].baselvl = false;
			group[n][m].grplvl 	= false;
			group[n][m].boltrev = false;
			group[n][m].identT	= NULL;
			group[n][m].Ee 		= 0.;
			group[n][m].paridx 	= NULL;
			group[n][m].qn 		= NULL;
			group[n][m].ql 		= NULL;
			group[n][m].qw 		= NULL;
		}

		// Determine group shell configurations
		for(int m=0;m<Nalevs[n];m++)
		{
			bool sameSt = true;
			for (int ns=0;ns<atom[n]->nS[m];ns++)
			{
				if (m == 0)
				{
					group[n][grpidx].Nlvlgrp++;
					sameSt = false;
					break;
				}
				else if (group[n][grpidx].qn[ns] == atom[n]->qn[m][ns] && \
						group[n][grpidx].ql[ns] == atom[n]->ql[m][ns] && \
						group[n][grpidx].qw[ns] == atom[n]->qw[m][ns])
					continue;
				else
				{
					grpidx++;
					group[n][grpidx].Nlvlgrp++;
					sameSt = false;
					break;
				}
			}

			if (sameSt == true)
			{
				if (group[n][grpidx].Nlvlgrp == 1 && group[n][grpidx].grplvl == false)
				{
					group[n][grpidx].baselvl = true;
					grpidx++;

					group[n][grpidx].grplvl = true;

					group[n][grpidx].qn = new int [atom[n]->nS[m]];
					group[n][grpidx].ql = new int [atom[n]->nS[m]];
					group[n][grpidx].qw = new int [atom[n]->nS[m]];

					for (int ns=0; ns<atom[n]->nS[m]; ns++)
					{
						group[n][grpidx].qn[ns] = atom[n]->qn[m][ns];
						group[n][grpidx].ql[ns] = atom[n]->ql[m][ns];
						group[n][grpidx].qw[ns] = atom[n]->qw[m][ns];
					}

					group[n][grpidx].Z = atom[n]->Z;
				}

				group[n][grpidx].Nlvlgrp++;
				group[n][grpidx].ge += atom[n]->ge[m];
				group[n][grpidx].Ee += atom[n]->ge[m]*atom[n]->Ee[m];
			}
			else
			{
				group[n][grpidx].nS = atom[n]->nS[m];

				group[n][grpidx].qn = new int [atom[n]->nS[m]];
				group[n][grpidx].ql = new int [atom[n]->nS[m]];
				group[n][grpidx].qw = new int [atom[n]->nS[m]];

				for (int ns=0; ns<atom[n]->nS[m]; ns++)
				{
					group[n][grpidx].qn[ns] = atom[n]->qn[m][ns];
					group[n][grpidx].ql[ns] = atom[n]->ql[m][ns];
					group[n][grpidx].qw[ns] = atom[n]->qw[m][ns];
				}

				group[n][grpidx].Z = atom[n]->Z;
				group[n][grpidx].ge += atom[n]->ge[m];
				group[n][grpidx].Ee += atom[n]->ge[m]*atom[n]->Ee[m];
			}
		}

		grpidx++;	// Increased due to 0 index
		if (grpidx == Ngroups[n])
			printf("Atom %s grouping shells (%i) have been determined...\n", \
					atom[n]->name.c_str(),Ngroups[n]);
		else
			printf("Grouping error detected: %i vs. %i\n",grpidx,Ngroups[n]);
	}
	printf("\n");

	for(int n=0;n<Natoms;n++)
		for(int m=0; m<Ngroups[n]; m++)
			group[n][m].Ee /= group[n][m].ge;

	// Descriptor for indexing between full levels list and levels grouped
	for(int n=0;n<Natoms;n++)
	{
		for(int m=0;m<Ngroups[n];m++)
			group[n][m].paridx = new int [group[n][m].Nlvlgrp];

		int grpidx = 0;
		int pidx = 0;

		// Match shell configurations
		for(int m=0;m<Nalevs[n];m++)
		{
			for (int ns=0; ns<atom[n]->nS[m]; ns++)
			{
				if (group[n][grpidx].qn[ns] == atom[n]->qn[m][ns] && \
						group[n][grpidx].ql[ns] == atom[n]->ql[m][ns] && \
						group[n][grpidx].qw[ns] == atom[n]->qw[m][ns] && \
						(group[n][grpidx].baselvl == false || m == 0))
					continue;
				else
				{
					pidx = 0;
					grpidx++;
					break;
				}
			}
			group[n][grpidx].paridx[pidx] = m;
			pidx++;
		}

		grpidx++;	// Increased due to 0 index
		if (grpidx == Ngroups[n])
			printf("Atom %s uniform grouping verified and applied...\n",atom[n]->name.c_str());
		else
			printf("Grouping error detected: %i vs. %i\n",grpidx,Ngroups[n]);
	}
	printf("\n");

	// Revert to explicit level variable if excited group, n', has only 1 level
	for(int n=0;n<Natoms;n++)
		for(int m=0; m<Ngroups[n]; m++)
			if (group[n][m].grplvl && group[n][m].Nlvlgrp == 1)
			{
				group[n][m-1].baselvl = false;	group[n][m-1].boltrev = true;
				group[n][m].grplvl = false;		group[n][m].boltrev = true;
			}

	return nStates;
}
// [Active]
int KRatesGrp::unif_cust_group()
{
	int nStates = 0;

	for(int n=0;n<Natoms;n++)
	{
		int nalevs = 0;

		// Determine the number of groups for the atom.
		stringstream itostr;
		itostr << atom[n]->Z+1;
		string filenum = itostr.str();

		int str_length = filenum.length();
		for (int i = 0; i < 2 - str_length; i++) filenum = "0" + filenum;

		string grpsfile	= "./Group_Data/Ar" + filenum + "groups.dat";

		int numlines;
		std::ifstream f(grpsfile.c_str());
		std::string line;
		for (numlines = 0; std::getline(f, line); ++numlines);
		printf("----Atom %s: Number of groups: %i\n",atom[n]->name.c_str(),numlines);
		f.close();

		Ngroups[n] = numlines;
		nStates += Ngroups[n];

		// Descriptors for number of states per group, shell configurations
		group[n] = new atomgroup[Ngroups[n]];
		int grpidx = 0;

		for(int m=0;m<Ngroups[n];m++)
		{
			group[n][m].Nlvlgrp = 0;
			group[n][m].ge		= 0;
			group[n][m].Ee 		= 0;
			group[n][m].paridx 	= NULL;
			group[n][m].qn 		= NULL;
			group[n][m].ql 		= NULL;
			group[n][m].qw 		= NULL;
		}

		int numch = 300;
		int fsval = 2;

		char fsstrg[numch];
		FILE* grpsFile;

		grpsFile = fopen(grpsfile.c_str(),"r"); // ---------- Allow for file selection

		stringstream gpstr;
		gpstr.str(fsstrg);

		string cb[fsval];

		int grp = 0;
		int numSh = 0;
		// Value Storage
		while (!feof(grpsFile))  {
			fgets(fsstrg,numch,grpsFile);
			gpstr.str(fsstrg);

			for (int i=0;i<fsval;i++) gpstr >> cb[i];

			numSh = atoi(cb[1].c_str());
			nalevs += numSh;
			if (nalevs>Nalevs[n])
			{
				if (numSh-(nalevs-Nalevs[n])>0) numSh = numSh-(nalevs-Nalevs[n]);
				else if (numSh-(nalevs-Nalevs[n])<=0) break;
			}
			group[n][grp].Nlvlgrp = numSh;		// Number of levels in group
			grp = atoi(cb[0].c_str());

			gpstr.str(std::string());
			gpstr.clear();

			grp++;
		}

		gpstr.str(std::string());
		gpstr.clear();
		fclose(grpsFile);

		for(int m=0;m<Ngroups[n];m++)
			if (group[n][m].Nlvlgrp < 1)
			{
				printf("%i and %i\n",n,m);
				printf("Invalid grouping data\n");
				exit(1);
			}

		for(int m=0;m<Ngroups[n];m++)
			group[n][m].paridx = new int [group[n][m].Nlvlgrp];

		grpidx = 0;
		int lvlcnt = 0;
		for(int m=0;m<Nalevs[n];m++)
		{
			group[n][grpidx].paridx[lvlcnt] = m;
			lvlcnt++;
			group[n][grpidx].Z = atom[n]->Z;
			group[n][grpidx].ge += atom[n]->ge[m];
			group[n][grpidx].Ee += atom[n]->ge[m]*atom[n]->Ee[m];

			gstate[n][m] = grpidx;

			if (lvlcnt >= group[n][grpidx].Nlvlgrp)
			{
				lvlcnt = 0;
				grpidx++;
			}
		}

		for(int m=0; m<Ngroups[n]; m++)
			group[n][m].Ee /= group[n][m].ge;
	}

	return nStates;
}
// [Active]
int KRatesGrp::bolt_cust_group()
{
	int m;
	int nStates = 0;

	for(int n=0;n<Natoms;n++)
	{
		int nalevs = 0;

		// Determine the number of groups for the atom.
		stringstream itostr;
		itostr << atom[n]->Z+1;
		string filenum = itostr.str();

		int str_length = filenum.length();
		for (int i = 0; i < 2 - str_length; i++) filenum = "0" + filenum;

		string grpsfile	= "./Group_Data/Ar" + filenum + "groups.dat";

		int numlines;
		std::ifstream f(grpsfile.c_str());
		std::string line;
		for (numlines = 0; std::getline(f, line); ++numlines);
		f.close();

		printf("---- Atom +%i: Declared groups: %i / ",atom[n]->Z,numlines);

		int numch = 400;
		int fsval = 2;

		char fsstrg[numch];
		FILE* grpsFile;

		grpsFile = fopen(grpsfile.c_str(),"r"); // ---------- Allow for file selection

		stringstream gpstr;
		gpstr.str(fsstrg);

		string cb[fsval];

		int grp = 0;
		int numSh = 0;
		// Value Storage
		while (!feof(grpsFile))  {
			fgets(fsstrg,numch,grpsFile);
			gpstr.str(fsstrg);

			for (int i=0;i<fsval;i++) gpstr >> cb[i];

			numSh = atoi(cb[1].c_str());
			if (numSh>1 && grp == atoi(cb[0].c_str()))	numlines++;

			gpstr.str(std::string());
			gpstr.clear();

			grp++;
		}

		gpstr.str(std::string());
		gpstr.clear();
		fclose(grpsFile);

		printf("Contributing variables: %i\n",numlines);

		// Descriptors for number of states per group, shell configurations
		Ngroups[n] = numlines;
		nStates += Ngroups[n];
		group[n] = new atomgroup[Ngroups[n]];

		for(m=0;m<Ngroups[n];m++)
		{
			group[n][m].Nlvlgrp = 0;
			group[n][m].ge		= 0;
			group[n][m].Ee 		= 0;
			group[n][m].paridx 	= NULL;
			group[n][m].qn 		= NULL;
			group[n][m].ql 		= NULL;
			group[n][m].qw 		= NULL;
		}

		grpsFile = fopen(grpsfile.c_str(),"r"); // ---------- Allow for file selection
		gpstr.str(fsstrg);

		m = 0;	grp = 0;
		// Value Storage
		while (!feof(grpsFile))  {
			fgets(fsstrg,numch,grpsFile);
			gpstr.str(fsstrg);

			for (int i=0;i<fsval;i++) gpstr >> cb[i];

			numSh = atoi(cb[1].c_str());

			if (numSh==1)
			{
				group[n][m].Nlvlgrp = 1;
				group[n][m].baselvl = false;
				group[n][m].grplvl 	= false;
				group[n][m].boltrev = false;
				group[n][m].identT	= NULL;
				nalevs++;
				//printf("Atom +%i, Group %i: %i levels\n",n,m,group[n][m].Nlvlgrp);
			}
			else if (numSh>1 && grp == atoi(cb[0].c_str()))
			{
				group[n][m].Nlvlgrp = 1;
				group[n][m].baselvl = true;
				group[n][m].grplvl 	= false;
				group[n][m].boltrev = false;
				group[n][m].identT	= NULL;
				//printf("Atom +%i, Group %i: %i levels\n",n,m,group[n][m].Nlvlgrp);
				m++;
				group[n][m].Nlvlgrp = numSh-1;
				group[n][m].baselvl = false;
				group[n][m].grplvl 	= true;
				group[n][m].boltrev = false;
				group[n][m].identT	= NULL;
				//printf("Atom +%i, Group %i: %i levels\n",n,m,group[n][m].Nlvlgrp);
				nalevs += numSh;
			}

			gpstr.str(std::string());
			gpstr.clear();

			if (nalevs == Nalevs[n])	break;
			m++;	grp++;
		}

		gpstr.str(std::string());
		gpstr.clear();
		fclose(grpsFile);

		for(int m=0;m<Ngroups[n];m++)
			if (group[n][m].Nlvlgrp < 1)
			{
				printf("Invalid grouping data: Atom +%i, Group %i: %i levels\n", \
						n,m,group[n][m].Nlvlgrp);
				exit(1);
			}

		for(int m=0;m<Ngroups[n];m++)
			group[n][m].paridx = new int [group[n][m].Nlvlgrp];

		int grpidx = 0;
		int lvlcnt = 0;
		for(int m=0;m<Nalevs[n];m++)
		{
			group[n][grpidx].paridx[lvlcnt] = m;
			group[n][grpidx].Z = atom[n]->Z;
			group[n][grpidx].ge += atom[n]->ge[m];
			group[n][grpidx].Ee += atom[n]->ge[m]*atom[n]->Ee[m];
			lvlcnt++;

			gstate[n][m] = grpidx;

//			printf("Total (%i,%i): ge = %i and Ee = %3.3e\n",n,grpidx,group[n][grpidx].ge,group[n][grpidx].Ee);
//			printf("\tAdded (%i,%i): ge = %i and Ee = %3.3e\n",n,m,atom[n]->ge[m],atom[n]->Ee[m]);

			if (lvlcnt >= group[n][grpidx].Nlvlgrp)
			{
				lvlcnt = 0;
				grpidx++;
			}
		}
		for(int m=0; m<Ngroups[n]; m++)
			group[n][m].Ee /= group[n][m].ge;
	}

	return nStates;
}
// [Active]
int KRatesGrp::QSS_group()
{
	int nStates = 0;

	for(int n=0;n<Natoms;n++)
	{
		// Descriptors for number of states per group, shell configurations
		Ngroups[n] = 1;
		nStates++;
		group[n] = new atomgroup[1];

		group[n][0].Nlvlgrp = 1;
		group[n][0].ge		= 0;
		group[n][0].Ee 		= 0;
		group[n][0].paridx 	= NULL;
		group[n][0].qn 		= NULL;
		group[n][0].ql 		= NULL;
		group[n][0].qw 		= NULL;
		group[n][0].baselvl = false;
		group[n][0].grplvl 	= false;
		group[n][0].boltrev = false;
		group[n][0].identT	= NULL;
		group[n][0].paridx = new int [1];
			group[n][0].paridx[0] = 0;
		group[n][0].Z = atom[n]->Z;
		group[n][0].ge += atom[n]->ge[0];
		group[n][0].Ee += atom[n]->Ee[0];
	}

	countQSSAtom();

	return nStates;
}

void KRatesGrp::printGroup()
{
	FILE* fp;

	// Create the output directory if it does not exist
    char levelsDir[100];
	sprintf(levelsDir,"mkdir -p ");
	sprintf(levelsDir + strlen(levelsDir),outputDir);
	sprintf(levelsDir + strlen(levelsDir),"Level_Grouping/");
	system(levelsDir);
	levelsDir[0] = '\0';

	// List grouping results on a per ion basis
	sprintf(levelsDir,outputDir);
	sprintf(levelsDir + strlen(levelsDir),"Level_Grouping/");
    for(int n=0;n<Natoms;n++)
    {
    	string lvlFile;
		stringstream convert;
		convert << levelsDir << "AtomIdx_" << n <<".dat";
		lvlFile = convert.str();
		fp = fopen(lvlFile.c_str(),"w");

    	for(int m=0;m<Ngroups[n];m++)
    	{
    		fprintf(fp,"Group %i:\t",m);
    		for(int p=0;p<group[n][m].Nlvlgrp;p++)
    			fprintf(fp,"%i, ",group[n][m].paridx[p]);

    		if (GrpScheme == "Boltzmann" || GrpScheme == "Custom_Boltzmann")
    		{
    			if (group[n][m].baselvl) fprintf(fp,"-Boltzmann base level");
    			if (group[n][m].grplvl) fprintf(fp,"-Boltzmann excited level(s)");
    		}
    		fprintf(fp,"\n");
    	}
    	fflush(fp);
    }
    fclose(fp);
}

/*
 * Distribute densities across atomic states
 */
void KRatesGrp::AtomDist(double ae, double Na, double TeV)
{
	int n,m,t,na;

	FILE* btFile;

	char currDir[100];
	const char* boltzInit = "BoltzInit.txt";

	strcpy(currDir,outputDir);
	strcat(currDir,boltzInit);

	btFile = fopen(currDir,"w"); // File opened: Boltzmann distributions are plotted
	currDir[0] = '\0';
	fprintf(btFile,"Energy Level\tTheo:ln(N/g)\tNum:ln(N/g)\tg\tQ\n");

	float* ndens;
	ndens = (float*) malloc (Nstate*sizeof(float));
	for(n=0;n<Nstate;n++) ndens[n] = 0.;
//	ndens = (float*) malloc (KRates::Nstate*sizeof(float));
//	for(n=0;n<KRates::Nstate;n++) ndens[n] = 0.;

	double NLv1tot = (1.-ae)*Na;
	double NLv2tot = ae*Na;

	for(int n=0;n<Natoms;n++) {
		// Evenly distribute electrons to ionized stages
		int numLvs = Ngroups[n];
		double Qpart = 0.;
		double Za = atom[n]->Z;

		printf("Atom %i: %i levels\n",n,numLvs);

		for (int m=0;m<numLvs;m++) {
			double gval = group[n][m].ge;
			double Eval = group[n][m].Ee;

			Qpart += gval*exp(-Eval/TeV);
		}

		// Density Assignment
		for (int m=0;m<numLvs;m++) {
			int ka   = state[n][m];
			double gval = group[n][m].ge;
			double Eval = group[n][m].Ee;

			double NLv = 1.;
			if (Za == 0.)		NLv = NLv1tot;
			else if (Za == 1.)	NLv = NLv2tot;

			double nAssign = NLv*gval*exp(-Eval/TeV)/Qpart;
			double logNg = log10(nAssign/gval);

			if (logNg < 0.) 	ndens[ka] = 5.e-1*gval; // Value chosen out of convenience for initial iterations
			else 				ndens[ka] += nAssign;

			fprintf(btFile,"%3.3e \t %3.3e \t %3.3e \t %3.3e \t %3.3e\n",Eval,logNg,log10(ndens[ka]/gval),gval,Qpart);
		}
		fprintf(btFile,"\n");

	}

	fflush(btFile);
	fclose(btFile);


	for (na=0;na<Natoms;na++)
		for (m=0;m<Ngroups[na];m++)
		{
			int ka = state[na][m];
			Ncr_init[ka] += ndens[ka];
		}

	free(ndens);

	printf("Applied atomic density distribution...\n\n");
}

/*
 * Check for duplicate groups
 */
void KRatesGrp::checkGroup()
{
	if (GrpScheme == "Uniform" || GrpScheme == "Boltzmann")
		for(int n=0;n<Natoms;n++)
			for(int m=0;m<Ngroups[n];m++) {
				if (m == Ngroups[n]-1)
					break;
				else
					for (int l=m+1;l<Ngroups[n];l++)
					{
						if (group[n][m].nS != group[n][l].nS)
							continue;
						if (l==m+1 && group[n][m].baselvl && group[n][l].grplvl)
							continue;
						int sameShells = 0;
						for (int ns=0;ns<group[n][m].nS;ns++)
						{
							if (group[n][m].qn[ns] == group[n][l].qn[ns] && \
									group[n][m].ql[ns] == group[n][l].ql[ns] && \
									group[n][m].qw[ns] == group[n][l].qw[ns])
							{
								sameShells++;
								continue;
							}
							else
								break;
						}
						if (sameShells == group[n][m].nS)
						{
							printf("Duplicate shells exist across groups!\n");
							//printf("%i %i is the same with %i %i\n",n,m,n,l);
							exit(1);
						}
					}
			}

	printf("No duplicate groups present within system!\n\n");
}

/*
 * Add grouped levels to solver's variable listing
 */
void KRatesGrp::unifyAtoms()
{
	Nlvful = KRates::Nstate;

	for(int na=0;na<Natoms;na++) state[na]  = new int[Ngroups[na]];

	for(int na=0;na<Natoms;na++)
		for(int n=0;n<Ngroups[na];n++)
		{
			state[na][n] = nq++;
			//printf("State %i, %i: %i\n",na,n,state[na][n]);
		}

	printf("Atomic solver variables unified!\n\n");
}

/*
 * Group bound-bound type transitions
 */
void KRatesGrp::add_BBtrans(bool irrad)
{
	KRates::add_BBtrans(irrad);

//	printf("Preparing bound-bound grouped rates based on %s grouping:\n",GrpScheme.c_str());
//
//	if (GrpScheme == "Full" || GrpScheme == "SS")
//		full_BB_copy();
//	else if (GrpScheme == "Uniform" || GrpScheme == "Custom_Uniform")
//		full_BB_copy();
////		unif_BB_group();
//	else if (GrpScheme == "Boltzmann" || GrpScheme == "Custom_Boltzmann")
//		full_BB_copy();
////		bolt_BB_group();
//	else if (GrpScheme == "QSS")
//		printf("Incorporating BB transitions into effective BF transitions");
//	else
//	{
//		printf("Not a valid grouping option: %s",GrpScheme.c_str());
//		exit(1);
//	}
//
//	printf("\n");

//	printBBRates();
}

void KRatesGrp::full_BB_copy()
{
	int n,m,t,l,u,g,na,ia,bb,tg,tr,gpl,gpu;
	bool lbool,ubool;

	// ----------------------------------------
	// Electron-Impact: Bound-Bound Transitions
	// ----------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBBeXD;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbexd[n].low[0];
		int Nlev = KRates::Nalevs[na];

		l = KRates::bbexd[n].low[1];
		u = KRates::bbexd[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		NBBeXiso++;
	}

	// --- Isolate reducible transitions
	bbexdiso = new transition_iso[NBBeXiso];

	for(bb=0;bb<NBBeXiso;bb++)
	{
		// Only count relevant level transitions
		na = KRates::bbexd[bb].low[0];
		int Nlev = KRates::Nalevs[na];

		l = KRates::bbexd[bb].low[1];
		u = KRates::bbexd[bb].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		bbexdiso[bb].glotridx = bb;
		bbexdiso[bb].low[0] = na;
		bbexdiso[bb].upp[0] = na;
		bbexdiso[bb].low[1] = l;
		bbexdiso[bb].upp[1] = u;
		bbexdiso[bb].low[2] = KRates::bbexd[bb].low[2];
		bbexdiso[bb].upp[2] = KRates::bbexd[bb].upp[2];
		bbexdiso[bb].low[3] = l;
		bbexdiso[bb].upp[3] = u;
		bbexdiso[bb].dE = KRates::bbexd[bb].dE;
	}

	if (KRates::NBBeXD != NBBeXiso)
	{
		printf("Number of electron-impact excitation transitions stored do not match (%i vs. %i)\n", \
				KRates::NBBeXD,NBBeXiso);
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Electron-impact transitions isolated from full list...\n", NBBeXiso);

	// ------------------------------------------
	// Radiation-Induced: Bound-Bound Transitions
	// ------------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBBlns;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbrad[n].low[0];
		int Nlev = Nalevs[na];

		l = KRates::bbrad[n].low[1];
		u = KRates::bbrad[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		NBBlniso++;
	}

	// --- Isolate reducible transitions
	bbradiso = new transition_iso[NBBlniso];

	for(bb=0;bb<NBBlniso;bb++)
	{
		// Only count relevant level transitions
		na = KRates::bbrad[bb].low[0];
		int Nlev = KRates::Nalevs[na];

		l = KRates::bbrad[bb].low[1];
		u = KRates::bbrad[bb].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		bbradiso[bb].glotridx = bb;
		bbradiso[bb].low[0] = na;
		bbradiso[bb].upp[0] = na;
		bbradiso[bb].low[1] = l;
		bbradiso[bb].upp[1] = u;
		bbradiso[bb].low[2] = KRates::bbrad[bb].low[2];
		bbradiso[bb].upp[2] = KRates::bbrad[bb].upp[2];
		bbradiso[bb].low[3] = l;
		bbradiso[bb].upp[3] = u;
		bbradiso[bb].dE = KRates::bbrad[bb].dE;
	}

	if (KRates::NBBlns != NBBlniso)
	{
		printf("Number of bound-bound lines stored do not match (%i vs. %i)\n", \
				KRates::NBBlns,NBBlniso);
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Photon-induced transitions isolated from full list...\n", NBBlniso);
}

void KRatesGrp::unif_BB_group()
{
	int n,m,t,l,u,g,na,ia,bb,tr,gpl,gpu;
	bool lbool,ubool;

    // ----------------------------------------
    // Electron-Impact: Bound-Bound Transitions
    // ----------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBBeXD;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbexd[n].low[0];
		int Nlev = KRates::Nalevs[na];

		l = KRates::bbexd[n].low[1];
		u = KRates::bbexd[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p] && !lbool)
				{
					lbool = true; gpl = g;
				}
				if (u == group[na][g].paridx[p] && !ubool)
				{
					ubool = true; gpu = g;
				}
			}
			if(lbool && ubool) break;
		}

		if (ubool == true && lbool == true && gpl != gpu) NBBeXiso++;
		else if (ubool == true && lbool == true && gpl == gpu) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bbexdiso = new transition_iso[NBBeXiso];

	bb = 0;
	for(n=0;n<KRates::NBBeXD;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbexd[n].low[0];
		int Nlev = KRates::Nalevs[na];

		l = KRates::bbexd[n].low[1];
		u = KRates::bbexd[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p] && !lbool)
				{
					lbool = true; gpl = g;
				}
				if (u == group[na][g].paridx[p] && !ubool)
				{
					ubool = true; gpu = g;
				}
			}
			if (lbool && ubool) break;
		}

		if (ubool == true && lbool == true && gpl != gpu)
		{
			bbexdiso[bb].glotridx = n;
			bbexdiso[bb].low[0] = na;
			bbexdiso[bb].upp[0] = na;
			bbexdiso[bb].low[1] = l;
			bbexdiso[bb].upp[1] = u;
			bbexdiso[bb].low[2] = KRates::bbexd[n].low[2];
			bbexdiso[bb].upp[2] = KRates::bbexd[n].upp[2];
			bbexdiso[bb].low[3] = gpl;
			bbexdiso[bb].upp[3] = gpu;
			bbexdiso[bb].dE = KRates::bbexd[n].dE;
			bb++;
		}
	}

	if (bb != NBBeXiso)
	{
		printf("Number of electron-impact excitation transitions stored do not match (%i vs. %i)\n",bb,NBBeXiso);
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Electron-impact transitions isolated from full list...\n",NBBeXiso);

	// ------------------------------------------
	// Radiation-Induced: Bound-Bound Transitions
	// ------------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBBlns;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbrad[n].low[0];
		int Nlev = Nalevs[na];

		l = KRates::bbrad[n].low[1];
		u = KRates::bbrad[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p] && !lbool)
				{
					lbool = true; gpl = g;
				}

				if (u == group[na][g].paridx[p] && !ubool)
				{
					ubool = true; gpu = g;
				}
			}
			if(lbool && ubool) break;
		}

		if (ubool && lbool && gpl != gpu) NBBlniso++;
		else if (ubool && lbool && gpl == gpu) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bbradiso = new transition_iso[NBBlniso];

	bb = 0;
	for(n=0;n<KRates::NBBlns;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbrad[n].low[0];
		int Nlev = KRates::Nalevs[na];

		l = KRates::bbrad[n].low[1];
		u = KRates::bbrad[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p] && !lbool)
				{
					lbool = true; gpl = g;
				}
				if (u == group[na][g].paridx[p] && !ubool)
				{
					ubool = true; gpu = g;
				}
			}
			if (lbool && ubool) break;
		}

		if (ubool == true && lbool == true && gpl != gpu)
		{
			bbradiso[bb].glotridx = n;
			bbradiso[bb].low[0] = na;
			bbradiso[bb].upp[0] = na;
			bbradiso[bb].low[1] = l;
			bbradiso[bb].upp[1] = u;
			bbradiso[bb].low[2] = KRates::bbrad[n].low[2];
			bbradiso[bb].upp[2] = KRates::bbrad[n].upp[2];
			bbradiso[bb].low[3] = gpl;
			bbradiso[bb].upp[3] = gpu;
			bbradiso[bb].dE = KRates::bbrad[n].dE;
			bb++;
		}
	}

	if (bb != NBBlniso)
	{
		printf("Number of bound-bound lines stored do not match (%i vs. %i)\n",bb,NBBlniso);
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Photon-induced transitions isolated from full list...\n", NBBlniso);
}

void KRatesGrp::bolt_BB_group()
{
	int n,m,t,l,u,g,na,ia,bb,tg,tr,gpl,gpu;
	bool lbool,ubool;

    // ----------------------------------------
    // Electron-Impact: Bound-Bound Transitions
    // ----------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBBeXD;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbexd[n].low[0];
		int Nlev = KRates::Nalevs[na];

		l = KRates::bbexd[n].low[1];
		u = KRates::bbexd[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++)
			{
				if (l == group[na][g].paridx[p] && !lbool)
				{
					lbool = true; gpl = g;
				}
				if (u == group[na][g].paridx[p] && !ubool)
				{
					ubool = true; gpu = g;
				}
			}
			if(lbool && ubool) break;
		}

		if (ubool == true && lbool == true && gpl != gpu) NBBeXiso++;
		else if (ubool == true && lbool == true && gpl == gpu) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bbexdiso = new transition_iso[NBBeXiso];

	bb = 0;
	for(n=0;n<KRates::NBBeXD;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbexd[n].low[0];
		int Nlev = KRates::Nalevs[na];

		l = KRates::bbexd[n].low[1];
		u = KRates::bbexd[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++)
			{
				if (l == group[na][g].paridx[p] && !lbool)
				{
					lbool = true; gpl = g;
				}
				if (u == group[na][g].paridx[p] && !ubool)
				{
					ubool = true; gpu = g;
				}
			}
			if (lbool && ubool) break;
		}

		if (ubool == true && lbool == true && gpl != gpu)
		{
			bbexdiso[bb].glotridx = n;
			bbexdiso[bb].low[0] = na;
			bbexdiso[bb].upp[0] = na;
			bbexdiso[bb].low[1] = l;
			bbexdiso[bb].upp[1] = u;
			bbexdiso[bb].low[2] = KRates::bbexd[n].low[2];
			bbexdiso[bb].upp[2] = KRates::bbexd[n].upp[2];
			bbexdiso[bb].low[3] = gpl;
			bbexdiso[bb].upp[3] = gpu;
			bbexdiso[bb].dE = KRates::bbexd[n].dE;
			bb++;
		}
	}

	if (bb != NBBeXiso)
	{
		printf("Number of electron-impact excitation transitions stored do not match (%i vs. %i)\n",bb,NBBeXiso);
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Electron-impact transitions isolated from full list...\n",NBBeXiso);

	// ------------------------------------------
	// Radiation-Induced: Bound-Bound Transitions
	// ------------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBBlns;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbrad[n].low[0];
		int Nlev = Nalevs[na];

		l = KRates::bbrad[n].low[1];
		u = KRates::bbrad[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p] && !lbool)
				{
					lbool = true; gpl = g;
				}
				if (u == group[na][g].paridx[p] && !ubool)
				{
					ubool = true; gpu = g;
				}
			}
			if (lbool && ubool) break;
		}

		if (ubool && lbool && gpl != gpu) NBBlniso++;
		else if (ubool && lbool && gpl == gpu) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("%i: %i and %i\n",na,gpl,gpu);
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bbradiso = new transition_iso[NBBlniso];

	bb = 0;
	for(n=0;n<KRates::NBBlns;n++)
	{
		// Only count relevant level transitions
		na = KRates::bbrad[n].low[0];
		int Nlev = KRates::Nalevs[na];

		l = KRates::bbrad[n].low[1];
		u = KRates::bbrad[n].upp[1];
		if (l >= Nlev || u >=Nlev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p] && !lbool)
				{
					lbool = true; gpl = g;
				}
				if (u == group[na][g].paridx[p] && !ubool)
				{
					ubool = true; gpu = g;
				}
			}
			if (lbool && ubool) break;
		}

		if (ubool == true && lbool == true && gpl != gpu)
		{
			bbradiso[bb].glotridx = n;
			bbradiso[bb].low[0] = na;
			bbradiso[bb].upp[0] = na;
			bbradiso[bb].low[1] = l;
			bbradiso[bb].upp[1] = u;
			bbradiso[bb].low[2] = KRates::bbrad[n].low[2];
			bbradiso[bb].upp[2] = KRates::bbrad[n].upp[2];
			bbradiso[bb].low[3] = gpl;
			bbradiso[bb].upp[3] = gpu;
			bbradiso[bb].dE = KRates::bbrad[n].dE;
			bb++;
		}
	}

	if (bb != NBBlniso)
	{
		printf("Number of bound-bound lines stored do not match (%i vs. %i)\n",bb,NBBlniso);
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Photon-induced transitions isolated from full list...\n", NBBlniso);
}

void KRatesGrp::printBBRates()
{
	int t,l,u,bb,tg;

//	FILE* fp;
//	string ratexfil;
//	string ratexdir = "./Results/Grouping_Results/Rates_Grouping/e-XD/CollEx";
//	system("exec rm -r ./Results/Grouping_Results/Rates_Grouping/e-XD/CollEx*");
//	for(bb=0;bb<NBBeXD;bb++)
//	{
//		int atom_index = bbexd[bb].low[0];
//		l = bbexd[bb].low[1];
//		u = bbexd[bb].upp[1];
//
//		stringstream convert;
//		convert << ratexdir << "_" << atom_index << "_" <<  l \
//				<< "_" << atom_index << "_" <<  u << ".dat";
//		ratexfil = convert.str();
//
//		fp = fopen(ratexfil.c_str(),"w");
//
//		for(t=0;t<=NTstep;t++)
//		{
//			double TeV = exp(tab_lnT[t]);
//			tg = NTstpgp;
//			if (GrpScheme == "Boltzmann" || GrpScheme == "Custom_Boltzmann")
//				fprintf(fp,"%6.3f     %12.5e     %12.5e\n", \
//						TeV,tab_keX_bg[bb][t][tg],tab_keD_bg[bb][t][tg]);
//			else
//				fprintf(fp,"%6.3f     %12.5e     %12.5e\n",TeV,tab_keX[bb][t],tab_keD[bb][t]);
//		}
//
//		fflush(fp);
//	}
//	fclose(fp);
}

/*
 * Group bound-free type transitions
 */
void KRatesGrp::add_BFtrans(bool irrad)
{
	KRates::add_BFtrans(irrad);
	KRates::add_AItrans();

//	printf("\nPreparing bound-free grouped rates based on %s grouping:\n",GrpScheme.c_str());
//
//	if (GrpScheme == "Full" || GrpScheme == "SS")
//		full_BF_copy();
//	else if (GrpScheme == "Uniform" || GrpScheme == "Custom_Uniform")
//		full_BF_copy();
////		unif_BF_group();
//	else if (GrpScheme == "Boltzmann" || GrpScheme == "Custom_Boltzmann")
//		full_BF_copy();
////		bolt_BF_group();
//	else if (GrpScheme == "QSS")
//		add_QSS_BFtrans();
//	else {
//		printf("Not a valid grouping option: %s",GrpScheme.c_str());
//		exit(1);
//	}

	if (GrpScheme == "QSS")
		add_QSS_BFtrans();

	printf("\n");

//	printBFRates();
}

void KRatesGrp::full_BF_copy()
{
	int n,m,t,tg,l,u,g,na,ia,bb,bf,tr,gpl,gpu;
	bool lbool,ubool;

	// ---------------------------------------
	// Electron-Impact: Bound-Free Transitions
	// ---------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBFeIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfeir[n].low[0];
		ia = KRates::bfeir[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfeir[n].low[1];
		u = KRates::bfeir[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		NBFeIiso++;
	}

	// --- Isolate reducible transitions
	bfeiriso = new transition_iso[NBFeIiso];

	for(bf=0;bf<NBFeIiso;bf++)
	{
		// Only count relevant level transitions
		na = KRates::bfeir[bf].low[0];
		ia = KRates::bfeir[bf].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfeir[bf].low[1];
		u = KRates::bfeir[bf].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		bfeiriso[bf].glotridx = bf;
		bfeiriso[bf].low[0] = na;
		bfeiriso[bf].upp[0] = ia;
		bfeiriso[bf].low[1] = l;
		bfeiriso[bf].upp[1] = u;
		bfeiriso[bf].low[2] = KRates::bfeir[bf].low[2];
		bfeiriso[bf].upp[2] = KRates::bfeir[bf].upp[2];
		bfeiriso[bf].low[3] = l;
		bfeiriso[bf].upp[3] = u;
		bfeiriso[bf].dE = KRates::bfeir[bf].dE;
	}

	if (KRates::NBFeIR != NBFeIiso)
	{
		printf("Number of electron-impact ionization transitions stored do not match\n");
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Electron-impact transitions isolated from full list...\n", NBFeIiso);

	// -------------------------------------------------------
	// Autoionization/Electron Capture: Bound-Free Transitions
	// -------------------------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBFaIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfair[n].low[0];
		ia = KRates::bfair[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfair[n].low[1];
		u = KRates::bfair[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		NBFaIiso++;
	}

	// --- Isolate reducible transitions
	bfairiso = new transition_iso[NBFaIiso];

	for(bf=0;bf<NBFaIiso;bf++)
	{
		// Only count relevant level transitions
		na = KRates::bfair[bf].low[0];
		ia = KRates::bfair[bf].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfair[bf].low[1];
		u = KRates::bfair[bf].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		bfairiso[bf].glotridx = bf;
		bfairiso[bf].low[0] = na;
		bfairiso[bf].upp[0] = ia;
		bfairiso[bf].low[1] = l;
		bfairiso[bf].upp[1] = u;
		bfairiso[bf].low[2] = KRates::bfair[bf].low[2];
		bfairiso[bf].upp[2] = KRates::bfair[bf].upp[2];
		bfairiso[bf].low[3] = l;
		bfairiso[bf].upp[3] = u;
		bfairiso[bf].dE = KRates::bfair[bf].dE;
	}

	if (KRates::NBFaIR != NBFaIiso)
	{
		printf("Number of autoionization transitions stored do not match\n");
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Autoionization/Electron-Capture transitions isolated from full list...\n", NBFaIiso);

	// -----------------------------------------
 	// Radiation-Induced: Bound-Free Transitions
 	// -----------------------------------------

	// --- Count number of transitions to be reduced
	for(bf=0;bf<KRates::NBFpIR;bf++)
	{
		// Only count relevant level transitions
		na = KRates::bfpir[bf].low[0];
		ia = KRates::bfpir[bf].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfpir[bf].low[1];
		u = KRates::bfpir[bf].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		NBFpIiso++;
	}

	// --- Isolate reducible transitions
	bfpiriso = new transition_iso[NBFpIiso];

	for(bf=0;bf<KRates::NBFpIR;bf++)
	{
		// Only count relevant level transitions
		na = KRates::bfpir[bf].low[0];
		ia = KRates::bfpir[bf].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfpir[bf].low[1];
		u = KRates::bfpir[bf].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		bfpiriso[bf].glotridx = bf;
		bfpiriso[bf].low[0] = na;
		bfpiriso[bf].upp[0] = ia;
		bfpiriso[bf].low[1] = l;
		bfpiriso[bf].upp[1] = u;
		bfpiriso[bf].low[2] = KRates::bfpir[bf].low[2];
		bfpiriso[bf].upp[2] = KRates::bfpir[bf].upp[2];
		bfpiriso[bf].low[3] = l;
		bfpiriso[bf].upp[3] = u;
		bfpiriso[bf].dE = KRates::bfpir[bf].dE;
	}

	if (KRates::NBFpIR != NBFpIiso)
	{
		printf("Number of photoionization transitions stored do not match\n");
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Photon-induced transitions isolated from full list...\n",NBFpIiso);
}

void KRatesGrp::unif_BF_group()
{
	int n,m,t,l,u,g,na,ia,bb,bf,tr,gpl,gpu;
	bool lbool,ubool;

	// ---------------------------------------
	// Electron-Impact: Bound-Free Transitions
	// ---------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBFeIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfeir[n].low[0];
		ia = KRates::bfeir[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfeir[n].low[1];
		u = KRates::bfeir[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if(l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if(lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia) NBFeIiso++;
		else if (ubool == true && lbool == true && na == ia) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bfeiriso = new transition_iso[NBFeIiso];

	bf = 0;
	for(n=0;n<KRates::NBFeIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfeir[n].low[0];
		ia = KRates::bfeir[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfeir[n].low[1];
		u = KRates::bfeir[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if (lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia)
		{
			bfeiriso[bf].glotridx = n;
			bfeiriso[bf].low[0] = na;
			bfeiriso[bf].upp[0] = ia;
			bfeiriso[bf].low[1] = l;
			bfeiriso[bf].upp[1] = u;
			bfeiriso[bf].low[2] = KRates::bfeir[n].low[2];
			bfeiriso[bf].upp[2] = KRates::bfeir[n].upp[2];
			bfeiriso[bf].low[3] = gpl;
			bfeiriso[bf].upp[3] = gpu;
			bfeiriso[bf].dE = KRates::bfeir[n].dE;
			bf++;
		}
	}

	if (bf != NBFeIiso)
	{
		printf("Number of electron-impact ionization transitions stored do not match\n");
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Electron-impact transitions isolated from full list...\n",NBFeIiso);

	// -------------------------------------------------------
	// Autoionization/Electron Capture: Bound-Free Transitions
	// -------------------------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBFaIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfair[n].low[0];
		ia = KRates::bfair[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfair[n].low[1];
		u = KRates::bfair[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if(l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if(lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia) NBFaIiso++;
		else if (ubool == true && lbool == true && na == ia) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bfairiso = new transition_iso[NBFaIiso];

	bf = 0;
	for(n=0;n<KRates::NBFaIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfair[n].low[0];
		ia = KRates::bfair[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfair[n].low[1];
		u = KRates::bfair[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if (lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia)
		{
			bfairiso[bf].glotridx = n;
			bfairiso[bf].low[0] = na;
			bfairiso[bf].upp[0] = ia;
			bfairiso[bf].low[1] = l;
			bfairiso[bf].upp[1] = u;
			bfairiso[bf].low[2] = KRates::bfair[n].low[2];
			bfairiso[bf].upp[2] = KRates::bfair[n].upp[2];
			bfairiso[bf].low[3] = gpl;
			bfairiso[bf].upp[3] = gpu;
			bfairiso[bf].dE = KRates::bfair[n].dE;
			bf++;
		}
	}

	if (bf != NBFaIiso)
	{
		printf("Number of autoionization transitions stored do not match\n");
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Autoionization/Electron Capture transitions isolated from full list...\n",NBFaIiso);

	// -----------------------------------------
	// Radiation-Induced: Bound-Free Transitions
	// -----------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBFpIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfpir[n].low[0];
		ia = KRates::bfpir[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfpir[n].low[1];
		u = KRates::bfpir[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if(l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if(lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia) NBFpIiso++;
		else if (ubool == true && lbool == true && na == ia) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bfpiriso = new transition_iso[NBFpIiso];

	bf = 0;
	for(n=0;n<KRates::NBFpIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfpir[n].low[0];
		ia = KRates::bfpir[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfpir[n].low[1];
		u = KRates::bfpir[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if (lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia)
		{
			bfpiriso[bf].glotridx = n;
			bfpiriso[bf].low[0] = na;
			bfpiriso[bf].upp[0] = ia;
			bfpiriso[bf].low[1] = l;
			bfpiriso[bf].upp[1] = u;
			bfpiriso[bf].low[2] = KRates::bfpir[n].low[2];
			bfpiriso[bf].upp[2] = KRates::bfpir[n].upp[2];
			bfpiriso[bf].low[3] = gpl;
			bfpiriso[bf].upp[3] = gpu;
			bfpiriso[bf].dE = KRates::bfpir[n].dE;
			bf++;
		}
	}

	if (bf != NBFpIiso)
	{
		printf("Number of photoionization transitions stored do not match\n");
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Photon-induced transitions isolated from full list...\n",NBFpIiso);
}

void KRatesGrp::bolt_BF_group()
{
	int n,m,t,tg,l,u,g,na,ia,bb,bf,tr,gpl,gpu;
	bool lbool,ubool;

	// ---------------------------------------
	// Electron-Impact: Bound-Free Transitions
	// ---------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBFeIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfeir[n].low[0];
		ia = KRates::bfeir[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfeir[n].low[1];
		u = KRates::bfeir[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if(l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if(lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia) NBFeIiso++;
		else if (ubool == true && lbool == true && na == ia) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bfeiriso = new transition_iso[NBFeIiso];

	bf = 0;
	for(n=0;n<KRates::NBFeIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfeir[n].low[0];
		ia = KRates::bfeir[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfeir[n].low[1];
		u = KRates::bfeir[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if (lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia)
		{
			bfeiriso[bf].glotridx = n;
			bfeiriso[bf].low[0] = na;
			bfeiriso[bf].upp[0] = ia;
			bfeiriso[bf].low[1] = l;
			bfeiriso[bf].upp[1] = u;
			bfeiriso[bf].low[2] = KRates::bfeir[n].low[2];
			bfeiriso[bf].upp[2] = KRates::bfeir[n].upp[2];
			bfeiriso[bf].low[3] = gpl;
			bfeiriso[bf].upp[3] = gpu;
			bfeiriso[bf].dE = KRates::bfeir[n].dE;
			bf++;
		}
	}

	if (bf != NBFeIiso)
	{
		printf("Number of electron-impact ionization transitions stored do not match\n");
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Electron-impact transitions isolated from full list...\n",NBFeIiso);

	// ---------------------------------------
	// Autoionization/Electron-Capture: Bound-Free Transitions
	// ---------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBFaIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfair[n].low[0];
		ia = KRates::bfair[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfair[n].low[1];
		u = KRates::bfair[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if(l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if(lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia) NBFaIiso++;
		else if (ubool == true && lbool == true && na == ia) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bfairiso = new transition_iso[NBFaIiso];

	bf = 0;
	for(n=0;n<KRates::NBFaIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfair[n].low[0];
		ia = KRates::bfair[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfair[n].low[1];
		u = KRates::bfair[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if (lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia)
		{
			bfairiso[bf].glotridx = n;
			bfairiso[bf].low[0] = na;
			bfairiso[bf].upp[0] = ia;
			bfairiso[bf].low[1] = l;
			bfairiso[bf].upp[1] = u;
			bfairiso[bf].low[2] = KRates::bfair[n].low[2];
			bfairiso[bf].upp[2] = KRates::bfair[n].upp[2];
			bfairiso[bf].low[3] = gpl;
			bfairiso[bf].upp[3] = gpu;
			bfairiso[bf].dE = KRates::bfair[n].dE;
			bf++;
		}
	}

	if (bf != NBFaIiso)
	{
		printf("Number of autoionization transitions stored do not match\n");
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Autoionization/Electron Capture transitions isolated from full list...\n",NBFaIiso);

	// -----------------------------------------
	// Radiation-Induced: Bound-Free Transitions
	// -----------------------------------------

	// --- Count number of transitions to be reduced
	for(n=0;n<KRates::NBFpIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfpir[n].low[0];
		ia = KRates::bfpir[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfpir[n].low[1];
		u = KRates::bfpir[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if(l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if(lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia) NBFpIiso++;
		else if (ubool == true && lbool == true && na == ia) continue;
		else
		{
			printf("Error during grouping: Unable to find particular level\n");
			printf("Exiting simulation...\n\n");
			exit(1);
		}
	}

	// --- Isolate reducible transitions
	bfpiriso = new transition_iso[NBFpIiso];

	bf = 0;
	for(n=0;n<KRates::NBFpIR;n++)
	{
		// Only count relevant level transitions
		na = KRates::bfpir[n].low[0];
		ia = KRates::bfpir[n].upp[0];
		int Nnlev = KRates::Nalevs[na];
		int Nilev = KRates::Nalevs[ia];

		l = KRates::bfpir[n].low[1];
		u = KRates::bfpir[n].upp[1];
		if (l >= Nnlev || u >= Nilev ) continue;

		lbool = false; ubool = false;
		for(g=0;g<Ngroups[na];g++) {
			for(int p=0;p<group[na][g].Nlvlgrp;p++) {

				if (l == group[na][g].paridx[p])
				{
					lbool = true; gpl = g;
					break;
				}
			}
			if (lbool) break;
		}

		for(g=0;g<Ngroups[ia];g++) {
			for(int p=0;p<group[ia][g].Nlvlgrp;p++) {

				if (u == group[ia][g].paridx[p])
				{
					ubool = true; gpu = g;
					break;
				}
			}
			if (ubool) break;
		}

		if (ubool == true && lbool == true && na != ia)
		{
			bfpiriso[bf].glotridx = n;
			bfpiriso[bf].low[0] = na;
			bfpiriso[bf].upp[0] = ia;
			bfpiriso[bf].low[1] = l;
			bfpiriso[bf].upp[1] = u;
			bfpiriso[bf].low[2] = KRates::bfpir[n].low[2];
			bfpiriso[bf].upp[2] = KRates::bfpir[n].upp[2];
			bfpiriso[bf].low[3] = gpl;
			bfpiriso[bf].upp[3] = gpu;
			bfpiriso[bf].dE = KRates::bfpir[n].dE;
			bf++;
		}
	}

	if (bf != NBFpIiso)
	{
		printf("Number of photoionization transitions stored do not match\n");
		printf("Exiting simulation...\n\n");
		exit(1);
	}
	else
		printf("%i Photon-induced transitions isolated from full list...\n",NBFpIiso);
}

void KRatesGrp::printBFRates()
{
	int t,tg,bf;

//	// Print grouped rate calculations
//	FILE* fp;
//	string ratiofil;
//	string ratiodir = "./Results/Grouping_Results/Rates_Grouping/e-IR/CollIo";
//	system("exec rm -r ./Results/Grouping_Results/Rates_Grouping/e-IR/CollIo*");
//	for(bf=0;bf<NBFeIR;bf++)
//	{
//		int neu_index = bfeir[bf].low[0];
//		int ion_index = bfeir[bf].upp[0];
//		int n = bfeir[bf].low[1];
//		int m = bfeir[bf].upp[1];
//
//		stringstream convert;
//		convert << ratiodir << "_" << neu_index << "_" <<  n \
//				<< "_" << ion_index << "_" <<  m << ".dat";
//		ratiofil = convert.str();
//
//		fp = fopen(ratiofil.c_str(),"w");
//
//		for(t=0;t<=NTstep;t++)
//		{
//			double TeV = exp(tab_lnT[t]);
//			tg = NTstpgp;
//			if (GrpScheme == "Boltzmann" || GrpScheme == "Custom_Boltzmann")
//				fprintf(fp,"%6.3f     %12.5e     %12.5e\n", \
//						TeV,tab_keI_bg[bf][t][tg],tab_keR_bg[bf][t][tg]);
//			else
//				fprintf(fp,"%6.3f     %12.5e     %12.5e\n",TeV,tab_keI[bf][t],tab_keR[bf][t]);
//		}
//
//		fflush(fp);
//	}
//	fclose(fp);
}

// Define a logarithmic range of electron densities (for QSS)
void KRatesGrp::setNRange(int N, float Nmin, float Nmax)
{
    NNstep = N;
    tab_lnN = new double[NNstep+1];
    lnN_min = log(Nmin);
    lnN_max = log(Nmax);
    dlnN = (lnN_max-lnN_min)/NNstep;

    // Linearly define densities in log space
    for(int n=0;n<=NNstep;n++)
        tab_lnN[n] = lnN_min + n * dlnN;
}

/*
 * Calculate effective QSS rates
 */
void KRatesGrp::add_QSS_BFtrans()
{
	printf("Tabulating QSS CIR BF rates\n");

	// count transition
	NBFeIR_QSS=0;
	for(int na=0;na<KRates::Natoms;na++)
	{
		AtomX* neu = KRates::atom[na];
		AtomX* ion = neu->ion;  // ion (upper ion stage)
		if(ion) NBFeIR_QSS++;
	}

    // Initialize indexing arrays
	bfeir_QSS = new transition[NBFeIR_QSS];
	tab_keI_QSS = cArray<double>(NNstep+1,NTstep+1,NBFeIR_QSS);
	tab_keR_QSS = cArray<double>(NNstep+1,NTstep+1,NBFeIR_QSS);
	tab_Lz = cArray<double>(NNstep+1,NTstep+1,KRates::Natoms);
	tab_Kz = cArray<double>(NNstep+1,NTstep+1,KRates::Natoms);
	tab_r0 = cArray<double*>(NNstep+1,NTstep+1,Natoms_QSS);
	tab_rp = cArray<double*>(NNstep+1,NTstep+1,Natoms_QSS);

	int bf = 0;
	for(int i=0;i<KRates::Natoms;i++)
	{
		AtomX* neu = KRates::atom[i];
		AtomX* ion = neu->ion;  // ion (upper ion stage)
		if(!ion) continue;
		printf("... Atom %s\n",neu->name.c_str());
		int na = i;
		int ia = -1;
		// if ion exists, which is it in the list?
		for(int k=0;k<KRates::Natoms;k++)
		{   // select entry based on name comparison
			if(ion->name == KRates::atom[k]->name) { ia = k; break; }
		}
		if(ia < 0)
		{
			printf("could not find ion (%s) in list\n",ion->name.c_str()); abort();
		}

		for (int n=0;n<=NNstep;n++)
		for (int t=0;t<=NTstep;t++) {
			tab_r0[n][t][na] = new double[KRates::Nalevs[na]-1];
			tab_rp[n][t][na] = new double[KRates::Nalevs[na]-1];
		}

		double dE = 0.;
		for(int t=0;t<neu->NBFeIR;t++)
		{
			int low = neu->bfeir[t].low[1];
			int upp = neu->bfeir[t].upp[1];
			if (low==0 && upp==0)
				dE = neu->bfeir[t].dE;
		}
		if (dE==0.) printf("could not find transition %s(0) -> %s(0)\n",
				neu->name.c_str(),ion->name.c_str());

		bfeir_QSS[bf].low[0] = na;
		bfeir_QSS[bf].upp[0] = ia;
		bfeir_QSS[bf].low[1] = 0;
		bfeir_QSS[bf].upp[1] = 0;
		bfeir_QSS[bf].low[2] = neu->Z;
		bfeir_QSS[bf].upp[2] = ion->Z;
		bfeir_QSS[bf].dE     = dE;

		for (int n=0;n<=NNstep;n++)
		{
//			printf("n = %i\n",n);
			for (int t=0;t<=NTstep;t++)
			{
//				printf("T = %i\n",t);
				double Ne = exp(tab_lnN[n]);
				double Te = exp(tab_lnT[t]);
				int Nxstate = KRates::Nalevs[na]-1;
				double* Bi  = cArray<double>(Nxstate+1);
				double* Si  = cArray<double>(Nxstate+1);
				double *X0 = tab_r0[n][t][bf];
				double *Xp = tab_rp[n][t][bf];

				double Scr = 0.;
				double Acr = 0.;
				double Lz  = 0.;
				double Kz  = 0.;

				compute_qss_rates(na,Ne,Te,X0,Xp,Bi,Si,Scr,Acr,Lz,Kz);

				tab_keI_QSS[n][t][na] = Scr;
				tab_keR_QSS[n][t][na] = Acr;
				tab_Lz[n][t][na]      = Lz;
				tab_Kz[n][t][na]      = Kz;

				del_cArray<double>(Nxstate+1,Bi);
				del_cArray<double>(Nxstate+1,Si);
			}
		}

		// check for nan rates
		for (int n=0;n<=NNstep;n++)
		for (int t=0;t<=NTstep;t++)
		{
			if (isnan(tab_keR_QSS[n][t][na]) || isinf(tab_keR_QSS[n][t][na])) {
				for (int k=t+1;k<=NTstep;k++)
					if (!isnan(tab_keR_QSS[n][k][na]) && !isinf(tab_keR_QSS[n][k][na])) {
						tab_keR_QSS[n][t][na] = tab_keR_QSS[n][k][na];
						tab_keI_QSS[n][t][na] = tab_keI_QSS[n][k][na];
						break;
					}
			}
		}

		bf++;
	}

	printf("%d QSS collisional I/R\n",NBFeIR_QSS);
}

void KRatesGrp::compute_qss_rates(int a, double Ne, double Te, double* X0, double* Xp, double* Bi, double* Si,  double& Seff, double &Aeff, double &Lz, double &Kz)
{
	int Nxstate = KRates::Nalevs[a]-1;
	double* DD0 = cArray<double>(Nxstate*Nxstate);
	double* DDp = cArray<double>(Nxstate*Nxstate);
	double* D0 = cArray<double>(Nxstate);
	double* Dp = cArray<double>(Nxstate);
	double** XD = cArray<double>(Nxstate+1,KRates::Nstate+1);
	double** ci  = cArray<double>(Nxstate+1,2);

	AtomX &A0 = *KRates::atom[a];
	AtomX &A1 = *KRates::atom[a]->ion;

	set_cArray<double>(Nxstate,X0,0.);
	set_cArray<double>(Nxstate,Xp,0.);
	set_cArray<double>(Nxstate+1,Si,1.);
	set_cArray<double>(Nxstate+1,Bi,1.);

	// find interpolation point in tables
	double lnTe = max(lnT_min,min(.99*lnT_max,log(Te) ) );
	float xt = (lnTe-lnT_min)/dlnT;
	int it = (int) xt; xt -= it;

	for(int i=0;i<=Nxstate;i++)
		Bi[i] = (double)A0.ge[i]/(double)A0.ge[0]*exp(-(A0.Ee[i]-A0.Ee[0])/Te);

	for(int k=0;k<KRates::NBBeXD;k++)
	{
		int alow = KRates::bbexd[k].low[0];
		int llow = KRates::bbexd[k].low[1];
		int lupp = KRates::bbexd[k].upp[1];
		if (alow != a) continue;
		double k_lu = (1.-xt)*KRates::tab_keX[k][it] + xt*KRates::tab_keX[k][it+1];
		double k_ul = (1.-xt)*KRates::tab_keD[k][it] + xt*KRates::tab_keD[k][it+1];
		XD[llow][lupp] += Ne*k_lu;
		XD[lupp][llow] += Ne*k_ul;
	}

	for(int k=0;k<KRates::NBBlns;k++)
	{
		int alow = KRates::bbrad[k].low[0];
		int llow = KRates::bbrad[k].low[1];
		int lupp = KRates::bbrad[k].upp[1];
		if (alow != a) continue;
		XD[lupp][llow] += KRates::rad_Abb[k];
	}

	for(int k=0;k<KRates::NBFeIR;k++)
	{
		int alow = KRates::bfeir[k].low[0];
		int llow = KRates::bfeir[k].low[1];
		int lupp = KRates::bfeir[k].upp[1];
		double dE = KRates::bfeir[k].dE;
		if (alow != a) continue;
		double k_lu = (1.-xt)*KRates::tab_keI[k][it]    + xt*KRates::tab_keI[k][it+1];
		double k_ul = (1.-xt)*KRates::tab_keR[k][it]    + xt*KRates::tab_keR[k][it+1];
		double r_ul = (1.-xt)*KRates::tab_kpR[k][it][0] + xt*KRates::tab_kpR[k][it+1][0];
		if (lupp==0)
		{
			double grat= (double)A0.ge[llow]/(double)A1.ge[lupp];
			double T3o2= Te*sqrt(Te);
			double Sf = grat/T3o2/6.e27*exp(dE/Te);
//			if (isinf(Sf)) printf("Inf Sf for k = %i/%i\n",k,KRates::NBFeIR);
//			if (isnan(Sf)) printf("NaN Sf for k = %i/%i\n",k,KRates::NBFeIR);
			ci[llow][0] += Ne*k_lu;
			ci[llow][1] += Ne*k_ul + r_ul;
			Si[llow]     = Sf;
		}
	}

	// ===============================================

	for(int l=0;l<=Nxstate;l++)
		if (Si[l] == 1.) Si[l] = Si[1]/Bi[1]*Bi[l];

	// Excited State Interactions
	for (int i=0;i<Nxstate;i++)
	for (int j=i+1;j<Nxstate;j++)
	{
		int llow = i+1;
		int lupp = j+1;
		double Buplo = (double)A0.ge[lupp]/(double)A0.ge[llow]*exp(-(A0.Ee[lupp]-A0.Ee[llow])/Te);
		double Bloup = (double)A0.ge[llow]/(double)A0.ge[lupp]*exp(-(A0.Ee[llow]-A0.Ee[lupp])/Te);
//		if (isinf(Bloup) || isinf(Buplo)) printf("Flag 2!");
		DD0[i*Nxstate+i] += XD[llow][lupp];
		DD0[j*Nxstate+j] += XD[lupp][llow];
//		DD0[i*Nxstate+j] -= XD[lupp][llow]*Bi[lupp]/Bi[llow];
		DD0[i*Nxstate+j] -= XD[lupp][llow]*Buplo;
//		DD0[j*Nxstate+i] -= XD[llow][lupp]*Bi[llow]/Bi[lupp];
		DD0[j*Nxstate+i] -= XD[llow][lupp]*Bloup;				// Bloup may be infinity
	}

	// Ground/Ionic State Interactions
	for (int i=0;i<Nxstate;i++)
	{
		int llow = 0;
		int lupp = i+1;
		double Bloup = (double)A0.ge[llow]/(double)A0.ge[lupp]*exp(-(A0.Ee[llow]-A0.Ee[lupp])/Te);
//		double Bf = Bi[lupp];
		DD0[i*Nxstate+i] += XD[lupp][llow];
//		if (isinf(Bloup)) printf("Flag 3!");
//		D0[i] += XD[llow][lupp]/Bf;
		D0[i] += XD[llow][lupp]*Bloup;				// Bloup may be infinity

		llow = i+1;
		lupp = 0;
		double Sf = Si[llow];
		DD0[i*Nxstate+i] += ci[llow][0];
//		if (isinf(Sf)) printf("Inf Sf for i = %i/%i\n",i,Nxstate);
//		if (isnan(Sf)) printf("NaN Sf for i = %i/%i\n",i,Nxstate);
//		if (Sf <= 0.) printf("Sf <= 0. for i = %i/%i\n",i,Nxstate);
		Dp[i] += ci[llow][1]/Sf;
	}

	// ================================================

//	printf("Entering solve for na = %i: Te = %2.1e and Ne = %2.1e\n",a,Te,Ne);

	cpy_cArray<double>(Nxstate*Nxstate,DDp,DD0);

	gsl_matrix_view mat0 	= gsl_matrix_view_array (DD0, Nxstate, Nxstate);
	gsl_matrix_view matp 	= gsl_matrix_view_array (DDp, Nxstate, Nxstate);
	gsl_vector_view rhs0 	= gsl_vector_view_array (D0, Nxstate);
	gsl_vector_view rhsp 	= gsl_vector_view_array (Dp, Nxstate);
	gsl_vector_view sol0 	= gsl_vector_view_array (X0, Nxstate);
	gsl_vector_view solp 	= gsl_vector_view_array (Xp, Nxstate);

	int _s;
	gsl_permutation * per = gsl_permutation_alloc (Nxstate);
	gsl_linalg_LU_decomp (&mat0.matrix, per, &_s);
	gsl_linalg_LU_solve  (&mat0.matrix, per, &rhs0.vector, &sol0.vector);
	gsl_linalg_LU_decomp (&matp.matrix, per, &_s);
	gsl_linalg_LU_solve  (&matp.matrix, per, &rhsp.vector, &solp.vector);

	Seff = ci[0][0]/Ne;
	Aeff = ci[0][1];

	/* option 1*/
	for (int i=0;i<Nxstate;i++)
	{
		int llow = 0;
		int lupp = i+1;
		double Bf = Bi[lupp];
		double Sf = Si[lupp];

//		Seff += XD[llow][lupp]/Ne;
//		Seff -= X0[i]*Bf*XD[lupp][llow]/Ne;

		if (isinf(Sf) || isnan(Xp[i]))	continue;
		Aeff += Xp[i]*Sf*XD[lupp][llow];			// May have infinity issues with Sf

//		if (isinf(Sf)) printf("Inf Sf for i = %i/%i with Xp = %2.1e and XD = %2.1e\n",i,Nxstate,Xp[i],XD[lupp][llow]);
	}

	/* option 2 */
	for (int i=0;i<Nxstate;i++)
	{
		int llow = i+1;
		double Bf = Bi[llow];
		double Sf = Si[llow];

		if (isinf(Bf) || isnan(X0[i]))	continue;
		Seff += X0[i]*Bf*ci[llow][0]/Ne;

//		Aeff += ci[llow][1];
//		Aeff -= Xp[i]*Sf*ci[llow][0]/Ne;
	}

//	if (isnan(Aeff) || isinf(Aeff))
//	{
//		Seff = ci[0][0]/Ne;
//		Aeff = ci[0][1];
//	}



	gsl_permutation_free (per);
	del_cArray<double>(Nxstate*Nxstate,DD0);
	del_cArray<double>(Nxstate*Nxstate,DDp);
	del_cArray<double>(Nxstate,D0);
	del_cArray<double>(Nxstate,Dp);
	del_cArray<double>(Nxstate+1,Nxstate+1,XD);
	del_cArray<double>(Nxstate+1,2,ci);
}

/*
 * Transfer elastic collision rates
 */
void KRatesGrp::xfer_ECexchg()
{
//	int t;
//
//	//=== Elastic Collisions ===//
//
//    // Initialize indexing arrays
//	tab_kEC = new double[NTstep+1];
//	for(t=0;t<=NTstep;t++)
//		tab_kEC[t] = 0.;
//
//	for(t=0;t<=NTstep;t++)
//		tab_kEC[t] = KRates::tab_kEC[t];
//
//    printf("Elastic collision rates successfully transferred to grouped rates!\n\n");
//
//    FILE* fp;
//    fp = fopen("./Results/Rates/ECRates.dat","w");
//	for(t=0;t<=NTstep;t++)
//	{
//		double TeV = exp(tab_lnT[t]);
//		fprintf(fp,"%3.3e",TeV);
//		fprintf(fp,"\t%3.3e",tab_kEC[t]);
//		fprintf(fp,"\n");
//	}
//	fclose(fp);
}

/*
 * Calculate the transition's Boltzmann weight
 */
double KRatesGrp::boltz_weight(int natm, int ngrp, int nlvl, double TgV)
{
	// Obtain lowest energy information
	int lv0 = group[natm][ngrp].paridx[0];
	double E0 = (double)atom[natm]->Ee[lv0];

	// Denominator
	int numLvs = group[natm][ngrp].Nlvlgrp;
	double den = 1.e-10;
	double invTgV=1./TgV;
	for(int m=0;m<numLvs;m++)
	{
		int lv = group[natm][ngrp].paridx[m];
		double ge = (double) atom[natm]->ge[lv];
		double Ee = (double) atom[natm]->Ee[lv];

		Ee -= E0;
		den += ge*exp(-Ee*invTgV);
	}

	// Numerator
	double ge = (double) atom[natm]->ge[nlvl];
	double Ee = (double) atom[natm]->Ee[nlvl];
	Ee -= E0;

	double num = ge*exp(-Ee*invTgV);

	// Calculate weighting factor
	double bgweight = num/den;

	return bgweight;
}

/*
 * Calculate the energy moments for the group
 */
double KRatesGrp::part_funct(int natm, int ngrp, double TgV)
{
	int numLvs = group[natm][ngrp].Nlvlgrp;

	double Qpart = 1.e-10;
	for(int m=0;m<numLvs;m++)
	{
		int lv = group[natm][ngrp].paridx[m];
		double ge = (double) atom[natm]->ge[lv];
		double Ee = (double) atom[natm]->Ee[lv];

		Qpart += ge/exp(Ee/TgV);
	}

	return Qpart;
}

double KRatesGrp::part_funct(int natm, int ngrp, int isolvl, double TgV)
{
	int numLvs = group[natm][ngrp].Nlvlgrp;
	double Qpart = 1.;
	double ge_iso = (double)atom[natm]->ge[isolvl];
	double Ee_iso = (double)atom[natm]->Ee[isolvl];

	for(int m=0;m<numLvs;m++)
	{
		int lvl = group[natm][ngrp].paridx[m];
		if (isolvl == lvl) continue;
		double ge = (double)atom[natm]->ge[lvl];
		double Ee = (double)atom[natm]->Ee[lvl];
		double grat = ge/ge_iso;
		double dE = Ee - Ee_iso;

		Qpart += grat*exp(-dE/TgV);
	}

	return Qpart;
}

double KRatesGrp::part_funct_renorm(int natm, int ngrp, int sublvl, float TgV)
{
	int numLvs = group[natm][ngrp].Nlvlgrp;

	double Qpart = 0.;
	for(int m=0;m<numLvs;m++)
	{
		int lvl = group[natm][ngrp].paridx[m];

		double ge = KRates::atom[natm]->ge[lvl];
		double Ee = KRates::atom[natm]->Ee[lvl];

		double Ee_norm = KRates::atom[natm]->Ee[sublvl];

		Qpart += ge*exp(-(Ee-Ee_norm)/TgV);
	}

	return Qpart;
}

double KRatesGrp::part_funct_E1mom(int natm, int ngrp, int sublvl, float TgV)
{
	int numLvs = group[natm][ngrp].Nlvlgrp;

	double Qpart = 0.;
	for(int m=0;m<numLvs;m++)
	{
		int lv = group[natm][ngrp].paridx[m];

		double ge = KRates::atom[natm]->ge[lv];
		double Ee = KRates::atom[natm]->Ee[lv];

		double Ee_norm = KRates::atom[natm]->Ee[sublvl];

		Qpart += ge*(Ee-Ee_norm)*exp(-(Ee-Ee_norm)/TgV);
	}

	return Qpart;
}

double KRatesGrp::part_funct_E2mom(int natm, int ngrp, int sublvl, float TgV)
{
	int numLvs = group[natm][ngrp].Nlvlgrp;

	double Qpart = 0.;
	for(int m=0;m<numLvs;m++)
	{
		int lv = group[natm][ngrp].paridx[m];

		double ge = KRates::atom[natm]->ge[lv];
		double Ee = KRates::atom[natm]->Ee[lv];

		double Ee_norm = KRates::atom[natm]->Ee[sublvl];

		Qpart += ge*(Ee-Ee_norm)*(Ee-Ee_norm)*exp(-(Ee-Ee_norm)/TgV);
	}

	return Qpart;
}

double KRatesGrp::mom1_Egrp(int natm, int ngrp, float TgV, bool full)
{
	int numLvs = group[natm][ngrp].Nlvlgrp;

	double Qpart = 0.;
	for(int m=0;m<numLvs;m++)
	{
		int lv = group[natm][ngrp].paridx[m];

		double ge = KRates::atom[natm]->ge[lv];
		double Ee = KRates::atom[natm]->Ee[lv];
		double E0 = group[natm][ngrp-1].Ee;

		Qpart += ge*(Ee-E0)*exp(-Ee/TgV);
	}

	if (full && numLvs > 1)
		Qpart /= part_funct(natm,ngrp,TgV) + part_funct(natm,ngrp-1,TgV);
	else
		Qpart /= part_funct(natm,ngrp,TgV);

	return Qpart;
}

double KRatesGrp::mom2_Egrp(int natm, int ngrp, float TgV, bool full)
{
	int numLvs = group[natm][ngrp].Nlvlgrp;

	double Qpart = 0.;
	for(int m=0;m<numLvs;m++)
	{
		int lv = group[natm][ngrp].paridx[m];

		double ge = KRates::atom[natm]->ge[lv];
		double Ee = KRates::atom[natm]->Ee[lv];
		double E0 = group[natm][ngrp-1].Ee;

		Qpart += ge*(Ee-E0)*(Ee-E0)*exp(-Ee/TgV);
	}

	if (full && numLvs > 1)
		Qpart /= part_funct(natm,ngrp,TgV) + part_funct(natm,ngrp-1,TgV);
	else
		Qpart /= part_funct(natm,ngrp,TgV);

	return Qpart;
}

/*
 * Determine the Boltzmann grouping energy gaps
 */
double KRatesGrp::bolt_Egap(int natm, int glow, int gupp,
		float Nratn, float Nratm, float TgVl, float TgVu)
{
	double dE 	= 0.;
	double ksi 	= 0.;

	// Lower group
	double E1 = 0.;
	ksi = mom2_Egrp(natm,glow,TgVl,true)/mom1_Egrp(natm,glow,TgVl,false) - mom1_Egrp(natm,glow,TgVl,true);
	if (group[natm][glow].grplvl)
		E1 = mom1_Egrp(natm,glow,TgVl,false) + ksi;
	else if (group[natm][glow].baselvl)
		E1 = group[natm][glow].Ee - ksi*Nratn;
	else
		E1 = group[natm][glow].Ee;

	// Upper group
	double E2 = 0.;
	ksi = mom2_Egrp(natm,gupp,TgVu,true)/mom1_Egrp(natm,gupp,TgVu,false) - mom1_Egrp(natm,gupp,TgVu,true);
	if (group[natm][gupp].grplvl)
		E2 = mom1_Egrp(natm,gupp,TgVu,false) + ksi;
	else if (group[natm][gupp].baselvl)
		E2 = group[natm][gupp].Ee - ksi*Nratm;
	else
		E2 = group[natm][gupp].Ee;

	dE = E2 - E1;
	return dE;
}

double KRatesGrp::bolt_Igap(int natm, int glow, int iatm, int gupp,
		float Nratn, float Nratm, float Nrat, float TgVl, float TgVu)
{
	double dE 	= 0.;
	double ksi	= 0.;
	double Igap = 0.;

	// Determined ionization gap between "ground" levels of each ion
	for (int n=0;n<NBFeIR;n++)
	{
		if (KRates::bfpir[n].low[0] == natm && KRates::bfpir[n].upp[0] == iatm && \
				KRates::bfpir[n].low[1] == 0 && KRates::bfpir[n].upp[1] == 0)
			Igap = KRates::bfpir[n].dE;
		else
			continue;

		if (n == NBFeIR-1 && Igap>0.)
		{
			printf("Unable to find an ionization gap from ion index %i\n",natm);
			exit(1);
		}
	}

	double E1 = 0.;
	ksi = mom2_Egrp(natm,glow,TgVl,true)/mom1_Egrp(natm,glow,TgVl,false) - mom1_Egrp(natm,glow,TgVl,true);
	// Lower group
	if (group[natm][glow].grplvl)
		E1 = Igap - mom1_Egrp(natm,glow,TgVl,false) - ksi;
	else if (group[natm][glow].baselvl)
		E1 = Igap - group[natm][glow].Ee + ksi*Nratn;
	else
		E1 = Igap - group[natm][glow].Ee;

	double E2 = 0.;
	ksi = mom2_Egrp(iatm,gupp,TgVu,true)/mom1_Egrp(iatm,gupp,TgVu,false) - mom1_Egrp(iatm,gupp,TgVu,true);
	// Upper group
	if (group[iatm][gupp].grplvl)
	{
		E1 += mom1_Egrp(iatm,gupp,TgVu,false);
		E2 = ksi*Nrat;
	}
	else if (group[iatm][gupp].baselvl)
	{
		E1 += group[iatm][gupp].Ee;
		E2 = -ksi*Nratm*Nrat;
	}
	else
		E1 += group[iatm][gupp].Ee;

	dE = E2 + E1;
	return dE;
}
