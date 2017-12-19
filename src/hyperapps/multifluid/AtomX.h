/*
 * File:   AtomX.h
 * Author: hle
 *
 * Created on April 1, 2015, 6:00 PM
 */

#ifndef ATOMX_H
#define	ATOMX_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "constant.h"

using namespace std;

struct transition
{
   int low[3]; // 0: atom index, 1: level/group index, 2: charge for lower level
   int upp[3]; // 0: atom index, 1: level/group index, 2: charge for upper level
   float dE;   // energy threshold or ionization potential

   int NEpts ; // data for xs
   float lnEmin; // data for xs
   float lnEmax; // data for xs
   float dlnE; // data for xs
   double* xs_dat;
};

struct Aconfig
{
	float Ee;float Ei;
	int ge;
	vector<int> sublev_idx;
	int nS; // number of shell
	int qn[8],ql[8],qw[8];
};

class AtomX
{
private:
public:

	/* Atomic data */

	// Energy levels
    int    Nlevels; // Number of electronic levels
    int    Znuc;	// nuclear charge
    int    Z;       // charge
    int*   ge;      // level degeneracy
    int*   nS;      // Number of distinct shells considered
    int**  qn;      // Principal quantum number		// first idx: atomic state, second idx: shell
    int**  ql;		// Angular quantum number		// first idx: atomic state, second idx: shell
    int**  qw;		// Orbital Occupation Number	// first idx: atomic state, second idx: shell
    float* Ee;      // level energy [eV]
    float* Ei;      // /* ---- Ignored ---- */ ionization energy [eV]
    float  amu;     // mass in amu
    float* Asum;    // natural broadening HWHM

    int Nalevs;		// Number of active levels considered for full simulation

    // Atomic Lines
    int NBBlns;       // Number of lines
    transition* bbrad;
    double** bbradfg;  // first idx: line, second idx: wavelength in nm (0), dE (1), fosc(2), Arad(3)

    // e-impact BB
    int NBBeXD;        // Number of e-impact excitation/deexcitation
    int Nexdxs;        // Number of xs data points
    transition* bbexd;
    double*** bbeexxs;  // first idx: transition,  second idx: set, and
    				   // third idx: electron energy (0), xs (1) (excitation), and collision strength (2)
    double**  bbeexxs_intp;

    // e-impact BF
    int NBFeIR;        // Number of e-impact ionization/recombination
    int Neirxs;        // Number of xs data points
    transition* bfeir;
    double*** bfeioxs;  // first idx: transition,  second idx: set, and
    				   // third idx: electron energy (0), xs (1) (ionization), and collision strength (2)
    double**  bfeioxs_intp;

    // p-impact BF
    int NBFpIR;        // Number of photoionization/radiative recombination
    int Npirxs;        // Number of xs data points
    transition* bfpir;
    double*** bfpioxs;  // first idx: transition,  second idx: set, and third idx: photon energy (0) and xs (1) (photoionization)
    double**  bfpioxs_intp;

    // auto-ionization
    int NBFaIR;
    transition* bfair;
    double* bfaiort;

    // e-impact (DW) BF
    int NBFdIR;        // Number of distorted wave ionization/radiative recombination
	int Ndirxs;        // Number of xs data points
	transition* bfdir;
	double*** bfdioxs;  // first idx: transition,  second idx: set, and third idx: photon energy (0) and xs (1) (ionization)

	// e-impact EC
	int Neecxs;		    // Number of xs data points against eV
	double** bheecxs;   // first idx: data line,   second idx: electron energy (0) and xs (1) (elastic collisions)

	string name;	// name of atomic species (used for identification)
    AtomX* ion;			// pointer to next ion stage

    double mX;			// mass of atomic species


    AtomX();
    AtomX(string);
    virtual ~AtomX();

	void hydroIon();    
	void set_name(string);
    void set_ion(AtomX*);
    void AtomParams(string);
   	void setEi(float);
   	void setmass(float);
   	double memory_size();

    double eexc_rate(int bb, float TeV);
    double eion_rate(int bf, float TeV);
    double eexc_rate(int bb, float xlu, float TeV, bool fwd);
    double eion_rate(int bf, float xlu, float TeV, bool fwd);
    double* prec_rate(int bf, float TeV);
    double ecol_rate(float TeV);

    double ecxs_calc(double);

};

#endif /* ATOMX_H */
