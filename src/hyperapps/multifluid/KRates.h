/*
 * File:   KRates.h
 * Author: hle
 *
 * Created on April 28, 2015, 6:23 PM
 */

#ifndef KRATES_H
#define	KRATES_H

#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "cubature.h"

#include "AtomX.h"
#include "constant.h"

using namespace std;

struct params_exc { transition tran; float T; float Tr; float lambda; AtomX** A; };

class KRates
{
protected:
    int ia;				// Index for Atom sweep
    int bb;
    int bf;
    int nq;				// Index for State
public:
    int NTstep;
    double Trad;
    double lnT_min, lnT_max, dlnT;
    double* tab_lnT;

    int NBBlns;
    int NBBeXD;
    int NBFeIR;
    int NBFpIR;
    int NBFaIR;

	transition* bbexd;
	transition* bbrad;
	transition* bfeir;
	transition* bfpir;
	transition* bfair;

    double**	tab_keX;
    double**	tab_keD;
    double**	tab_keI;
    double**  	tab_keR;
    double** 	tab_kpI;
    double** 	tab_kDR;
    double*** 	tab_kpR;
    double*** 	tab_kpS;
    double*   	rad_Abb;
    double*   	kAI;

    int Natoms;
    int* Nalevs;
    AtomX** atom;
    int** state;
    int Nstate;

    // Initialization variables
    double* Ncr_init;

    const char* spcsname;
    const char* outputDir;

    KRates();
    KRates(const char*,int);
    virtual ~KRates();

    virtual void DistPrep(string);
   	void AtomDist(double,double,double);

    void setCount(int);
    void setRange(int,float,float);
    void setTPlanck(float);
    void setMass(AtomX*);
    void setOutputDir(const char*);
    virtual void addAtom(AtomX*);
    virtual void addAtom(AtomX*,int);
    virtual void addAtom(AtomX*,double);
    void nextIon(int,int);

    virtual void add_BBtrans(bool); // Excitation/De-excitations
    virtual void add_BFtrans(bool); // Ionization/Recombinations
    virtual void add_AItrans(); // Autoionization/electron capture
    void add_ECexchg(); 		// Elastic Collisions
    void check();

    double* solve_steady_state(double Ne, double Te, bool is_radiating);
};

#endif	/* KRATES_H */



