/*
 * CRParams.h
 *
 *  Created on: Jan 9, 2017
 *      Author: richard
 */

#ifndef CRPARAMS_H_
#define CRPARAMS_H_

#include "myIO.h"
#include <sstream>
#include <string.h>

using namespace std;

class CRParams
{
public:
	// Kinetic Rates
	int Natoms;
	int temp_bins;
	int dens_bins;
    double temp_min;
	double temp_max;
	double dens_min;
	double dens_max;
	string atom_sym;
	// CR Solver Configuration
	bool irrad_bool;
	bool rad_bool;
	bool fxdNe_bool;
	bool intNe_bool;
	bool intTh_bool;
	bool intTe_bool;
	string grouping_type;
	string output_dir;
	// Temporal Parametrization
	double t_begin;
	double t_end;
	double t_dt;
	double multxdt;
	// Initial Conditions
	double temp_elec;
	double temp_heavy;
	double temp_rad;
	double pressure;
	double dens_elec;
	// Output Conditions
	bool   print_log;
	double printdt;
	double logdtprt;
	// Spectra Calculation and Parameters
	bool   print_spectra;
	double wave_min;
	double wave_max;
	double wave_resol;
	// Tolerances
	double rel_tol;
	double abs_tol;

	double* atom_cutoffs;

	int p2s; // print to screen
	vector<string> ops_name;

	CRParams();
	~CRParams();

	void init(myiolib::jFileSys* jSys);
};


#endif /* CRPARAMS_H_ */
