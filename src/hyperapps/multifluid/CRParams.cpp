/*
 * CRParams.cpp
 *
 *  Created on: Nov 16, 2017
 *      Author: richard
 */

#include "CRParams.h"

CRParams::CRParams()
{
	// Kinetic Rates
	Natoms 		= 0;
	temp_bins 	= 0;
	dens_bins 	= 0;
	temp_min 	= 0.;
	temp_max 	= 0.;
	dens_min 	= 0.;
	dens_max 	= 0.;
	// CR Solver Configuration
	irrad_bool 	= false;
	rad_bool 	= false;
	fxdNe_bool 	= false;
	intNe_bool 	= false;
	intTh_bool 	= false;
	intTe_bool 	= false;
	// Temporal Parametrization
	t_begin 	= 0.;
	t_end 		= 0.;
	t_dt 		= 0.;
	multxdt 	= 0.;
	// Initial Conditions
	temp_elec 	= 0.;
	temp_heavy 	= 0.;
	temp_rad 	= 0.;
	pressure 	= 0.;
	dens_elec 	= 0.;
	// Output Conditions
	print_log 	= false;
	printdt 	= 0.;
	logdtprt 	= 0.;
	// Spectra Calculation and Parameters
	print_spectra = false;
	wave_min 	= 0.;
	wave_max 	= 0.;
	wave_resol 	= 0.;
	// Tolerances
	rel_tol 	= 0.;
	abs_tol 	= 0.;

	atom_cutoffs = NULL;

	p2s 		= 0;
}

CRParams::~CRParams()
{
	if (atom_cutoffs) delete[] atom_cutoffs;
}

void CRParams::init(myiolib::jFileSys* jSys) {

	// Kinetic Rates
	bool found_atom_sym 		= false;
	bool found_Natoms 			= false;
	bool found_temp_bins 		= false;
	bool found_temp_min 		= false;
	bool found_temp_max 		= false;
	bool found_dens_bins 		= false;
	bool found_dens_min 		= false;
	bool found_dens_max 		= false;
	// CR Solver Configuration
	bool found_irrad_bool 		= false;
	bool found_rad_bool 		= false;
	bool found_intNe_bool 		= false;
	bool found_fxdNe_bool 		= false;
	bool found_intTh_bool 		= false;
	bool found_intTe_bool 		= false;
	bool found_grouping_type 	= false;
	bool found_output_dir		= false;
	// Temporal Parametrization
	bool found_t_begin 			= false;
	bool found_t_end 			= false;
	bool found_t_dt 			= false;
	bool found_multxdt 			= false;
	// Initial Conditions
	bool found_temp_elec 		= false;
	bool found_temp_heavy 		= false;
	bool found_temp_rad 		= false;
	bool found_pressure 		= false;
	bool found_dens_elec 		= false;
	// Output Conditions
	bool found_print_log 		= false;
	bool found_printdt 			= false;
	bool found_logdtprt 		= false;
	// Spectra Calculation and Parameters
	bool found_print_spectra	= false;
	bool found_wave_min 		= false;
	bool found_wave_max			= false;
	bool found_wave_resol 		= false;
	// Solver Tolerances
	bool found_rel_tol 			= false;
	bool found_abs_tol 			= false;

	string fname;
	fname = jSys->case_dir+jSys->casename+".pars";
//		printf("%s\n",fname.c_str());

	myiolib::Scan& scanner1 = *(new myiolib::Scan(fname));
	scanner1.readInputList("CRPARS");
	vector<myiolib::jString> toks;

	while (scanner1.next_line < (int) scanner1.input.size()) {
		toks = scanner1.nextLine().clean().parse();
		if (toks.size() > 1) {

			/* --------------------------------- */
			/* ------ KINETIC RATES SETUP ------ */
			/* --------------------------------- */

			if (toks[0].contains("atom_sym")) {
				found_atom_sym = true;
				atom_sym = toks[1].str();
			}
			if (toks[0].contains("Natoms")) {
				found_Natoms = true;
				Natoms = toks[1].iValue();
			}
			if (toks[0].contains("temp_bins")) {
				found_temp_bins = true;
				temp_bins = toks[1].iValue();
			}
			if (toks[0].contains("temp_min")) {
				found_temp_min = true;
				temp_min = toks[1].dValue();
			}
			if (toks[0].contains("temp_max")) {
				found_temp_max = true;
				temp_max = toks[1].dValue();
			}
			if (toks[0].contains("dens_bins")) {
				found_dens_bins = true;
				dens_bins = toks[1].iValue();
			}
			if (toks[0].contains("dens_min")) {
				found_dens_min = true;
				dens_min = toks[1].dValue();
			}
			if (toks[0].contains("dens_max")) {
				found_dens_max = true;
				dens_max = toks[1].dValue();
			}

			/* ------------------------------------- */
			/* ------ CR SOLVER CONFIGURATION ------ */
			/* ------------------------------------- */

			if (toks[0].contains("ird_bool")) {
				found_irrad_bool = true;
				irrad_bool = toks[1].equals("true");
			}
			if (toks[0].contains("rad_bool")) {
				found_rad_bool = true;
				rad_bool = toks[1].equals("true");
			}
			if (toks[0].contains("fxdNe_bool")) {
				found_fxdNe_bool = true;
				fxdNe_bool = toks[1].equals("true");
			}
			if (toks[0].contains("intNe_bool")) {
				found_intNe_bool = true;
				intNe_bool = toks[1].equals("true");
			}
			if (toks[0].contains("intTh_bool")) {
				found_intTh_bool = true;
				intTh_bool = toks[1].equals("true");
			}
			if (toks[0].contains("intTe_bool")) {
				found_intTe_bool = true;
				intTe_bool = toks[1].equals("true");
			}
			if (toks[0].contains("grouping_type")) {
				found_grouping_type = true;
				grouping_type = toks[1].str();
			}
			if (toks[0].contains("output_dir")) {
				found_output_dir = true;
				output_dir = toks[1].str();
			}

			/* -------------------------------------- */
			/* ------ TEMPORAL PARAMETRIZATION ------ */
			/* -------------------------------------- */

			if (toks[0].contains("t_begin")) {
				found_t_begin = true;
				t_begin = toks[1].dValue();
			}
			if (toks[0].contains("t_dt")) {
				found_t_dt = true;
				t_dt = toks[1].dValue();
			}
			if (toks[0].contains("t_end")) {
				found_t_end = true;
				t_end = toks[1].dValue();
			}
			if (toks[0].contains("multxdt")) {
				found_multxdt = true;
				multxdt = toks[1].dValue();
			}

			/* -------------------------------- */
			/* ------ INITIAL CONDITIONS ------ */
			/* -------------------------------- */

			if (toks[0].contains("temp_elec")) {
				found_temp_elec = true;
				temp_elec = toks[1].dValue();
			}
			if (toks[0].contains("temp_heavy")) {
				found_temp_heavy = true;
				temp_heavy = toks[1].dValue();
			}
			if (toks[0].contains("temp_rad")) {
				found_temp_rad = true;
				temp_rad = toks[1].dValue();
			}
			if (toks[0].contains("pressure")) {
				found_pressure = true;
				pressure = toks[1].dValue();
			}
			if (toks[0].contains("dens_elec")) {
				found_dens_elec = true;
				dens_elec = toks[1].dValue();
			}

			/* ------------------------------- */
			/* ------ OUTPUT CONDITIONS ------ */
			/* ------------------------------- */

			if (toks[0].contains("print_log")) {
				found_print_log = true;
				print_log = toks[1].equals("true");
			}
			if (toks[0].contains("printdt")) {
				found_printdt = true;
				printdt = toks[1].dValue();
			}
			if (toks[0].contains("logdtprt")) {
				found_logdtprt = true;
				logdtprt = toks[1].dValue();
			}

			/* -------------------------------------- */
			/* ------ SPECTRAL PARAMETRIZATION ------ */
			/* -------------------------------------- */

			if (toks[0].contains("print_spectra")) {
				found_print_spectra = true;
				print_spectra = toks[1].equals("true");
			}
			if (toks[0].contains("wave_min")) {
				found_wave_min = true;
				wave_min = toks[1].dValue();
			}
			if (toks[0].contains("wave_max")) {
				found_wave_max = true;
				wave_max = toks[1].dValue();
			}
			if (toks[0].contains("wave_resol")) {
				found_wave_resol = true;
				wave_resol = toks[1].dValue();
			}

			/* ------------------------------- */
			/* ------ OUTPUT CONDITIONS ------ */
			/* ------------------------------- */

			if (toks[0].contains("rel_tol")) {
				found_rel_tol = true;
				rel_tol = toks[1].dValue();
			}
			if (toks[0].contains("abs_tol")) {
				found_abs_tol = true;
				abs_tol = toks[1].dValue();
			}

		}
	}

	delete &scanner1;

	// Kinetics
	if (!found_atom_sym) 		{ printf("Element symbol needed! Exiting...\n"); exit(1);	}
	if (!found_Natoms) 			{ printf("Natoms not found! default to 2\n"); Natoms = 2;	}
	if (!found_temp_bins) 		{ printf("temp_bins not found! default to 100\n"); temp_bins = 100;	}
	if (!found_temp_min) 		{ printf("temp_min not found! default to 0.03\n"); temp_min = 0.03;	}
	if (!found_temp_max) 		{ printf("temp_max not found! default to 1500\n"); temp_max = 1500;	}
	if (!found_dens_bins) 		{ printf("dens_bins not found! default to 51\n"); dens_bins = 51.;	}
	if (!found_dens_min) 		{ printf("dens_min not found! default to 1.e8\n"); dens_min = 1.e8;	}
	if (!found_dens_max) 		{ printf("dens_max not found! default to 1.e28\n"); dens_max = 1.e28;	}
	// CR Solver Configuration
	if (!found_irrad_bool) 		{ printf("rad_bool not found! default to false\n"); irrad_bool = false;	}
	if (!found_rad_bool) 		{ printf("rad_bool not found! default to false\n"); rad_bool = false;	}
	if (!found_fxdNe_bool) 		{ printf("fxdNe_bool not found! default to false\n"); fxdNe_bool = false;	}
	if (!found_intNe_bool) 		{ printf("intNe_bool not found! default to false\n"); intNe_bool = false;	}
	if (!found_intTh_bool) 		{ printf("intTh_bool not found! default to false\n"); intTh_bool = false;	}
	if (!found_intTe_bool) 		{ printf("intTe_bool not found! default to false\n"); intTe_bool = false;	}
	if (!found_grouping_type) 	{ printf("grouping_type not found! default to \"Full\"\n"); grouping_type = "Full";	}
	if (!found_output_dir) 		{ printf("output_dir not found! default to \"./Output/\"\n"); output_dir = "./Output/";	}
	// Temporal Parameters
	if (!found_t_begin) 		{ printf("t_begin not found! default to 0.\n"); t_begin = 0.;	}
	if (!found_t_end) 			{ printf("t_end not found! default to 1.e-3\n"); t_end = 1.e-3;	}
	if (!found_t_dt) 			{ printf("t_dt not found! default to 1.e-9\n"); t_dt = 1.e-9;	}
	if (!found_multxdt) 		{ printf("multxdt not found! default to 1.e0\n"); multxdt = 1.e0;	}
	// Initial Conditions
	if (!found_temp_elec) 		{ printf("temp_elec not found! default to 50\n"); temp_elec = 50.;	}
	if (!found_temp_heavy) 		{ printf("temp_heavy not found! default to temp_min\n"); temp_heavy = temp_min;	}
	if (!found_temp_rad) 		{ printf("temp_rad not found! default to 100\n"); temp_heavy = 100;	}
	if (!found_pressure) 		{ printf("pressure not found! default to 1.e-3\n"); pressure = 1.e-3;	}
	if (!found_dens_elec) 		{ printf("dens_elec not found! default to 1.e3\n"); dens_elec = 1.e3;	}
	// Output Conditions
	if (!found_print_log) 		{ printf("print_log not found! default to false\n"); print_log = false;	}
	if (!found_printdt) 		{ printf("printdt not found! default to 1.e-8\n"); printdt = 1.e-8;	}
	if (!found_logdtprt) 		{ printf("logdtprt not found! default to 0.125\n"); logdtprt = 0.125;	}
	// Spectral Parameters
	if (!found_print_spectra) 	{ printf("print_spectra not found! default to 0.\n"); print_spectra = false;	}
	if (!found_wave_min) 		{ printf("wave_min not found! default to 1.e-3\n"); wave_min = 200.;	}
	if (!found_wave_max) 		{ printf("wave_max not found! default to 1.e-9\n"); wave_max = 800.;	}
	if (!found_wave_resol) 		{ printf("wave_resol not found! default to 1.e0\n"); wave_resol = 1.;	}
	// Output Conditions
	if (!found_rel_tol) 		{ printf("rel_tol not found! default to 1.e-3\n"); rel_tol = 1.e-3;	}
	if (!found_abs_tol) 		{ printf("abs_tol not found! default to 1.e-1\n"); abs_tol = 1.e-1;	}

	// ===================================================

	atom_cutoffs = new double [Natoms];
	for (int na=0;na<Natoms;na++) atom_cutoffs[na] = 0.;
	fname = jSys->case_dir+jSys->casename+".xoffs";
//		printf("%s\n",fname.c_str());

	myiolib::Scan& scanner2 = *(new myiolib::Scan(fname));
	scanner2.readInputList("ATOMXOFFS");

	while(scanner2.next_line < (int) scanner2.input.size()) {
		toks = scanner2.nextLine().clean().parse();
		if (toks.size() > 1) {

			for (int na=0;na<Natoms;na++)
			{
				char atomidx[5];
				sprintf(atomidx,"%02d",na);
				string atomstr = string(atomidx);
				if (toks[0].contains(atomstr))
					atom_cutoffs[na] = toks[1].dValue();
			}
		}
	}

	for (int na=0;na<Natoms;na++)
		if (atom_cutoffs[na] <= 0.)
		{
			printf("Invalid atom cutoff for ion +%i\n",na);
			exit(1);
		}
		else
			printf("Atom cutoff for ion +%i: %3.3e\n",na,atom_cutoffs[na]);

	delete &scanner2;
}
