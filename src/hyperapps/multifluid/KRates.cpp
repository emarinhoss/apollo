#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>
#include <gsl/gsl_linalg.h>
#include "cubature.h"

#include "cArray.h"
#include "MathX.h"
#include "KRates.h"

using namespace std;


int f_exc_ion_transform(unsigned ndim, const double *k, void *params, unsigned fdim, double *fval)
{
	struct params_exc * fparams = (struct params_exc *) params;
	double dE     = fparams->tran.dE;
	double lnEmin = fparams->tran.lnEmin;
	double lnEmax = fparams->tran.lnEmax;
	double dlnE   = fparams->tran.dlnE;
	double* xsel = fparams->tran.xs_dat;
	double T = fparams->T;

	double z = k[0];
	double omz = 1.-z;
	double xstar = dE/T;
	double x0 = xstar+z/omz;
	double lnE= log(x0*T);


	double xt = (max(lnEmin,min(.9999*lnEmax,lnE ))-lnEmin)/dlnE;
	int it = (int) xt; xt -= it;
	double sig  = (1.-xt)*xsel[it] + xt*xsel[it+1];
//	double sig  = (1.-xt)*log(xsel[it]) + xt*log(xsel[it+1]);
//	sig = exp(sig);
//	if (lnE < lnEmin) sig = 0.;
//	if (lnE > lnEmax) sig = 0.;

	double f = 0.;
	f += x0 * exp(-x0+xstar) * sig;		// exp(-dE/T) factor multiplied externally
	f /= omz*omz; // divided by jacobian
	fval[0] = f;

	return 0; // success
}

int f_pho_ion_transform(unsigned ndim, const double *k, void *params, unsigned fdim, double *fval)
{
	struct params_exc * fparams = (struct params_exc *) params;
//	AtomX** A = fparams->A;
//	int na= fparams->tran.low[0];
//	int ia= fparams->tran.upp[0];
//	int l = fparams->tran.low[1];
//	int u = fparams->tran.upp[1];
	double dE     = fparams->tran.dE;
	double lnEmin = fparams->tran.lnEmin;
	double lnEmax = fparams->tran.lnEmax;
	double dlnE   = fparams->tran.dlnE;
	double* xsel = fparams->tran.xs_dat;
	double Tr = fparams->Tr;

	double z = k[0];
	double omz = 1.-z;
	double xstar = dE/Tr;
	double x0 = z/omz;			// \varepsilon = 0.5*mv^2
	double x1 = x0+xstar;		// \varepsilon + dE = h\nu
	double lnE = log(x1*Tr);
	double effmult = 1 / (exp(x1) - 1);

	// PI xs
	double xt = (max(lnEmin,min(.9999*lnEmax,lnE ))-lnEmin)/dlnE;
	int it = (int) xt; xt -= it;
	double sigPI  = (1.-xt)*xsel[it] + xt*xsel[it+1];
//	if (lnE > lnEmax) sigPI = 0.;
//	if (lnE < lnEmin) sigPI = 0.;

	double f[2];

	f[0] = 0.;
	f[0] += x1 * x1 * effmult * sigPI;
	f[0] /= omz*omz; // divided by jacobian
	fval[0] = f[0];

	f[1] = 0.;
	f[1] += x1 * x1 * x0 * effmult * sigPI;
	f[1] /= omz*omz; // divided by jacobian
	fval[1] = f[1];

	return 0; // success
}

int f_rad_rec_transform(unsigned ndim, const double *k, void *params, unsigned fdim, double *fval)
{
	struct params_exc * fparams = (struct params_exc *) params;
	AtomX** A = fparams->A;
	int na= fparams->tran.low[0];
//	int ia= fparams->tran.upp[0];
	int l = fparams->tran.low[1];
	int u = fparams->tran.upp[1];
	double grat= (double)A[na]->ge[l]/2./(double)A[na]->ion->ge[u];
	double dE     = fparams->tran.dE;
	double lnEmin = fparams->tran.lnEmin;
	double lnEmax = fparams->tran.lnEmax;
	double dlnE   = fparams->tran.dlnE;
	double* xsel = fparams->tran.xs_dat;
	double T = fparams->T;

	double Eer = 0.510999e6/T; // normalized electron rest mass energy

	double z = k[0];
	double omz = 1.-z;
	double xstar = dE/T;
	double x0 = z/omz;			// \varepsilon = 0.5*mv^2
	double x1 = x0+xstar;		// \varepsilon + dE = h\nu
	double lnE= log(x1*T);

	// PI xs
	double xt = (max(lnEmin,min(.9999*lnEmax,lnE ))-lnEmin)/dlnE;
	int it = (int) xt; xt -= it;
	double sigPI  = (1.-xt)*xsel[it] + xt*xsel[it+1];
//	if (lnE > lnEmax) sigPI = 0.;

	// compute recombination xs from DB
	double x0xsig = grat*(x1)*(x1)/Eer*sigPI;		// KE_electron*xs_rr

	double f[2];

	f[0] = 0.;
	f[0] += exp(-x0) * x0xsig;
	f[0] /= omz*omz; // divided by jacobian
	fval[0] = f[0];

	f[1] = 0.;
	f[1] += x0 * exp(-x0) * x0xsig;
	f[1] /= omz*omz; // divided by jacobian
	fval[1] = f[1];

	return 0; // success
}

int f_stim_rec_transform(unsigned ndim, const double *k, void *params, unsigned fdim, double *fval)
{
	struct params_exc * fparams = (struct params_exc *) params;
	AtomX** A = fparams->A;
	int na= fparams->tran.low[0];
//	int ia= fparams->tran.upp[0];
	int l = fparams->tran.low[1];
	int u = fparams->tran.upp[1];
	double grat= (double)A[na]->ge[l]/2./(double)A[na]->ion->ge[u];
	double dE     = fparams->tran.dE;
	double lnEmin = fparams->tran.lnEmin;
	double lnEmax = fparams->tran.lnEmax;
	double dlnE   = fparams->tran.dlnE;
	double* xsel = fparams->tran.xs_dat;
	double T 	= fparams->T;
	double Tr 	= fparams->Tr;

	double Eer = 0.510999e6/T; // electron rest mass energy

	double z = k[0];
	double omz = 1.-z;
	double xstar = dE/T;
	double x0 = z/omz;
	double x1 = x0+xstar;
	double lnE= log(x1*T);
	double x1r = x1*T/Tr;
	double effmult = 1/(exp(x1r)-1);

	// PI xs
	double xt = (max(lnEmin,min(.9999*lnEmax,lnE ))-lnEmin)/dlnE;
	int it = (int) xt; xt -= it;
	double sigPI  = (1.-xt)*xsel[it] + xt*xsel[it+1];
//	if (lnE > lnEmax) sigPI = 0.;

	// compute recombination xs from DB
	double x0xsig = grat*(x1)*(x1)/Eer*sigPI;

	double f[2];

	f[0] = 0.;
	f[0] += exp(-x0) * x0xsig * effmult;
	f[0] /= omz*omz; // divided by jacobian
	fval[0] = f[0];

	f[1] = 0.;
	f[1] += x0 * exp(-x0) * x0xsig * effmult;
	f[1] /= omz*omz; // divided by jacobian
	fval[1] = f[1];

	return 0; // success
}

// Default Constructor
KRates::KRates():atom(NULL),NTstep(-1),tab_lnT(NULL),NBBlns(0),NBBeXD(0),NBFeIR(0),NBFpIR(0),NBFaIR(0)
{
	outputDir = "";
	spcsname = "";
	tab_keX = tab_keD = tab_keI = tab_keR = tab_kpI = tab_kDR = NULL;
	tab_kpR = tab_kpS = NULL;
	rad_Abb = kAI = NULL;
	bbexd = bfeir = bfpir = bfair = bbrad = NULL;

	Trad = 0.;
	lnT_min = 0;
	nq	= 0;
	bf	= 0;
	dlnT	= 0;
	bb	= 0;
	state = NULL;
	Ncr_init = NULL;
	Natoms = 0;
	lnT_max = 0;
	ia	= 0;
	Nstate = 0;
	Nalevs = 0;
}

// Constructor [Active]
KRates::KRates(const char* spcs,int Na):atom(NULL),NTstep(-1),tab_lnT(NULL),NBBlns(0),NBBeXD(0),NBFeIR(0),NBFpIR(0),NBFaIR(0)
{
	spcsname = spcs;
    Natoms = Na;			// Number of atoms considered
    atom = new AtomX*[Na];
    Nalevs = new int[Na];	// Number of levels for associated atom
    ia = 0;					// Initialize atom index
    nq = 0;					// Initialize state index
    state = new int*[Na];	//
    Nstate=0;
    Natoms = 0; // reset and recount when adding atoms
    tab_keX = tab_keD = tab_keI = tab_keR = tab_kpI = tab_kDR = NULL;
	tab_kpR = tab_kpS = NULL;
	rad_Abb = kAI = NULL;
	bbexd = bfeir = bfpir = bfair = bbrad = NULL;
	outputDir = NULL;

	Trad = 0.;
	bb = 0; lnT_max = 0; dlnT = 0; Ncr_init = NULL;
	lnT_min = 0; bf = 0;
}

// Destructor
KRates::~KRates()
{
    if(tab_lnT) {delete[] tab_lnT;tab_lnT=NULL;}
    delete [] Nalevs;
    delete [] bbexd;
    delete [] bfeir;
    delete [] bfpir;
    delete [] bfair;
    delete [] bbrad;

    delete [] Ncr_init;

    del_cArray<double>(NBBeXD,NTstep+1,tab_keX);
    del_cArray<double>(NBBeXD,NTstep+1,tab_keD);
    del_cArray<double>(NBFeIR,NTstep+1,tab_keI);
    del_cArray<double>(NBFeIR,NTstep+1,tab_keR);
    del_cArray<double>(NBFaIR,NTstep+1,tab_kDR);
    if (tab_kpR) del_cArray<double>(NBFpIR,NTstep+1,2,tab_kpR);
    if (tab_kpS) del_cArray<double>(NBFpIR,NTstep+1,2,tab_kpS);
    if (tab_kpI) del_cArray<double>(NBFpIR,2,tab_kpI);
    del_cArray<double>(NBBlns,rad_Abb);
    del_cArray<double>(NBFaIR,kAI);

//	if (tab_kEC)
//		delete [] tab_kEC;

	for(int n=0;n<Natoms;n++)
	{
		delete atom[n];
		delete [] state[n];
	}
	delete [] atom;
	delete [] state;
}

// Define a logarithmic range of temperatures
void KRates::setRange(int N, float Tmin, float Tmax)
{
    NTstep = N;
    tab_lnT = new double[NTstep+1];
    lnT_min = log(Tmin);
    lnT_max = log(Tmax);
    dlnT = (lnT_max-lnT_min)/NTstep;

    // Linearly define temperatures in log space
    for(int t=0;t<=NTstep;t++)
    {
        tab_lnT[t] = lnT_min + t * dlnT;
    }
}

// Define a logarithmic range of temperatures
void KRates::setTPlanck(float Tr)
{
    Trad = Tr;
}

// Set the mass and nuclear charge of the atom
void KRates::setMass(AtomX* spec)
{
//	char* nm = (char *)spec->name.c_str();
	char num[2] = {0};
	int nind = 0;
	for (int i = 0; i < 2; i++) {
		if(isalpha(spec->name[i])) {
			num[nind] = spec->name[i];
			nind++;
		}
	}

	if	(num[0] == 'A' && num[1] == 'r') {
		spec->mX	= 39.9/6.022e23/1000;
		spec->Znuc	= 18;
	} else if (num[0] == 'C' && num[1] == 'l') {
		spec->mX = 35.45/6.022e23/1000;
		spec->Znuc	= 17;
	} else if (num[0] == 'S' && num[1] == 'i') {
		spec->mX = 28.09/6.022e23/1000;
		spec->Znuc	= 14;
	}
}

namespace
{
    // The output directory comes from the input deck and used to be pasted
    // straight into a shell command that runs `rm -rf`. Anything the shell
    // treats specially is rejected rather than escaped, because getting the
    // escaping subtly wrong here deletes the user's files.
    bool crOutputDirIsSafe(const char* name, std::string& why)
    {
        if (name == NULL || *name == '\0')
        {
            why = "it is empty, which would make the cleanup command 'rm -rf *' "
                  "in the working directory";
            return false;
        }
        const std::string dir(name);
        if (dir == "/" || dir == "/." || dir == "." || dir == "./" || dir == "..")
        {
            why = "it names the root or the working directory itself";
            return false;
        }
        const std::string allowed =
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-/";
        if (dir.find_first_not_of(allowed) != std::string::npos)
        {
            why = "it contains a character outside [A-Za-z0-9._-/]";
            return false;
        }
        if (dir.find("..") != std::string::npos)
        {
            why = "it contains '..'";
            return false;
        }
        return true;
    }

    // mkdir -p, without a shell. Creates every missing component of `path`.
    bool crMakeDirs(const std::string& path)
    {
        std::string sofar;
        for (size_t i = 0; i <= path.size(); ++i)
        {
            if (i == path.size() || path[i] == '/')
            {
                if (!sofar.empty() && sofar != ".")
                {
                    if (mkdir(sofar.c_str(), 0755) != 0 && errno != EEXIST)
                        return false;
                }
            }
            if (i < path.size())
                sofar += path[i];
        }
        return true;
    }

    // Directory names are built by concatenation; join them with a single
    // separator regardless of whether the configured directory ends in one.
    std::string crJoin(const std::string& base, const std::string& leaf)
    {
        if (base.empty() || base[base.size()-1] == '/')
            return base + leaf;
        return base + "/" + leaf;
    }
}

// Prepare the output directories
void KRates::setOutputDir(const char* name)
{
	std::string why;
	if (!crOutputDirIsSafe(name, why))
	{
		fprintf(stderr,
		        "\nCR model: refusing to use output directory '%s': %s.\n"
		        "  Set output_dir in the input deck to a plain relative or "
		        "absolute path.\n\n",
		        name ? name : "(null)", why.c_str());
		exit(1);
	}

	outputDir = name;
	const std::string base(outputDir);

	DIR* dir = opendir(outputDir);
	if (dir)
	{
		closedir(dir);
		// The old code prompted on stdin and then ran
		//     system("exec rm -rf <outputDir>*")
		// on a 'Y'. `resp` was uninitialised and scanf's return value ignored,
		// so under a batch scheduler - where stdin is closed and scanf returns
		// EOF without writing anything - the branch taken was whatever byte
		// happened to be on the stack, one possibility being the delete. And
		// every MPI rank ran it at once. Deleting is now opt-in, interactive
		// only, and never the default.
		bool wipe = false;
		if (isatty(fileno(stdin)))
		{
			printf("\n>> Overwrite possible contents within %s? (Y/N)\n", outputDir);
			fflush(stdout);
			int resp = getchar();
			wipe = (resp == 'Y' || resp == 'y');
			if (!wipe && resp != 'N' && resp != 'n')
			{
				printf("Incomprehensible response\n");
				exit(1);
			}
		}
		else
		{
			printf("CR model: reusing existing output directory %s "
			       "(not deleting its contents; stdin is not a terminal).\n",
			       outputDir);
		}

		if (wipe)
		{
			printf("Emptying previously-made directory: %s\n\n", outputDir);
			// base passed crOutputDirIsSafe, so it holds no shell metacharacter;
			// it is still quoted so a future relaxation of that check cannot turn
			// into command injection.
			const std::string cmd = "rm -rf '" + base + "'/*";
			if (system(cmd.c_str()) != 0)
				fprintf(stderr, "CR model: warning: could not empty %s\n", outputDir);
		}
	}
	else if (ENOENT == errno)
	{
		printf("\nCreating output directory: %s\n\n", outputDir);
		if (!crMakeDirs(base))
		{
			fprintf(stderr, "CR model: could not create %s: %s\n",
			        outputDir, strerror(errno));
			exit(1);
		}
	}
	else
	{
		fprintf(stderr, "\nCR model: cannot access %s: %s\n",
		        outputDir, strerror(errno));
		exit(1);
	}

	// Subfolder Creation. These went through `system("mkdir -p ...")` with
	// unbounded sprintf into a 150-byte stack buffer; mkdir(2) needs neither.
	char ion[8];
	for (int na = 0; na < Natoms; na++)
	{
		snprintf(ion, sizeof(ion), "%02d", na+1);
		const std::string d = crJoin(base, std::string("BoltzDist/") + spcsname + ion);
		if (!crMakeDirs(d))
			fprintf(stderr, "CR model: could not create %s: %s\n",
			        d.c_str(), strerror(errno));
	}

	for (int na = 0; na < Natoms; na++)
	{
		snprintf(ion, sizeof(ion), "%02d", na+1);
		const std::string d = crJoin(base, std::string("LineRad/") + spcsname + ion);
		if (!crMakeDirs(d))
			fprintf(stderr, "CR model: could not create %s: %s\n",
			        d.c_str(), strerror(errno));
	}

	const std::string spectra = crJoin(base, "Spectra");
	if (!crMakeDirs(spectra))
		fprintf(stderr, "CR model: could not create %s: %s\n",
		        spectra.c_str(), strerror(errno));
}

// Add atom excited states based on cutoff energy
void KRates::addAtom(AtomX* spec, double lim)
{
	int cnt=0;
	for (int n=0;n<spec->Nlevels;n++) {
		if (spec->Ee[n]>lim) break;
		cnt++;
	}
	addAtom(spec,cnt);
}

// Add atom excited states based on previously defined levels
void KRates::addAtom(AtomX* spec)
{
	addAtom(spec,spec->Nlevels);
}

// Add atom excited states by manually defined number of levels
void KRates::addAtom(AtomX* spec, int Nlev)
{
    atom[ia]   = spec;				// Sets pointer to corresponding species
    Nlev = min(Nlev,spec->Nlevels);
    Nalevs[ia] = Nlev;
    atom[ia]->Nalevs = Nlev;
    state[ia]  = new int[Nlev];
    for(int n=0;n<Nlev;n++) state[ia][n] = nq++;
    ia++;			// Increment to work to next AtomX

    printf("Adding atom %s\n",spec->name.c_str());
    printf("Nlevels = %d Emax = %f\n",Nlev,spec->Ee[Nlev-1]);
//    for (int n=0;n<Nlev;n++)
//		printf("Ee[%d] = %f\n",n,spec->Ee[n]);
    Nstate += Nlev;
	Natoms++;

    setMass(spec);

    printf("Mass of %s is %e\n\n",spec->name.c_str(),spec->mX);
}

/**
 * Excitation/De-excitation
 */
void KRates::add_BBtrans(bool irrad)
{
	int n,t,l,u,na,bb;

	/*
	 * Bound-bound radiation
	 */
	printf("Tabulating BB rates:\n");

	// Count transitions
	NBBlns=0.;
	for(na=0;na<Natoms;na++)
	{
		AtomX& A = *atom[na];
		int Nlev = Nalevs[na];
		for(n=0;n<A.NBBlns;n++)
		{
			// Only count relevant level transitions
			l = A.bbrad[n].low[1];
			u = A.bbrad[n].upp[1];
			if (l >= Nlev || u >=Nlev ) continue;
			NBBlns++;
		}
	}

	// Initialize indexing arrays
	bbrad = new transition[NBBlns];
	rad_Abb = new double[NBBlns];

	bb = 0;
	for(na=0;na<Natoms;na++)
	{
		AtomX& A = *atom[na];
		int Nlev = Nalevs[na];
		for(n=0;n<A.NBBlns;n++)
		{
			// Only store relevant level transitions
			l = A.bbrad[n].low[1];
			u = A.bbrad[n].upp[1];
			if (l >= Nlev || u >=Nlev ) continue;
			bbrad[bb].low[0] = na;
			bbrad[bb].upp[0] = na;
			bbrad[bb].low[1] = l;
			bbrad[bb].upp[1] = u;
			bbrad[bb].low[2] = A.Z;
			bbrad[bb].upp[2] = A.Z;
			bbrad[bb].dE     = A.bbrad[n].dE;
			rad_Abb[bb]      = A.bbradfg[n][3];		// Einstein A coefficient for Spontaneous Emissions
			bb++;
		}
	}

	/*
	 * Electron-impact excitation/deexcitation
	 */

	// Count transitions
	NBBeXD=0;
	for(na=0;na<Natoms;na++)
	{
		AtomX& A = *atom[na];
		int Nlev = Nalevs[na];
		for(n=0;n<A.NBBeXD;n++)
		{
			// Only count relevant level transitions
			l = A.bbexd[n].low[1];
			u = A.bbexd[n].upp[1];
			if (l >= Nlev || u >=Nlev ) continue;
			NBBeXD++;
		}
	}

	// Initialize indexing arrays
	bbexd = new transition[NBBeXD];
	tab_keX = cArray<double>(NBBeXD,NTstep+1);
	tab_keD = cArray<double>(NBBeXD,NTstep+1);

	bb = 0;
	for(na=0;na<Natoms;na++)
	{
		AtomX& A = *atom[na];
		printf("... Atom %s\n",A.name.c_str());
		int Nlev = Nalevs[na];	// Number of levels explored for transitions
		for(n=0;n<A.NBBeXD;n++)
		{
			// Only store relevant level transitions
			l = A.bbexd[n].low[1];
			u = A.bbexd[n].upp[1];
			if (l >= Nlev || u >=Nlev ) continue;
			bbexd[bb].low[0] = na;
			bbexd[bb].upp[0] = na;
			bbexd[bb].low[1] = l;
			bbexd[bb].upp[1] = u;
			bbexd[bb].low[2] = A.Z;
			bbexd[bb].upp[2] = A.Z;
			bbexd[bb].dE     = A.bbexd[n].dE;

			bbexd[bb].NEpts = A.bbexd[n].NEpts;
			bbexd[bb].lnEmin= A.bbexd[n].lnEmin;
			bbexd[bb].lnEmax= A.bbexd[n].lnEmax;
			bbexd[bb].dlnE  = A.bbexd[n].dlnE;
			bbexd[bb].xs_dat= A.bbexd[n].xs_dat;

			for(t=0;t<=NTstep;t++)
			{
				double TeV = exp(tab_lnT[t]);
				double x_lu= bbexd[bb].dE/TeV;
				double grat= (double)A.ge[l]/(double)A.ge[u];

//				double kf = A.eexc_rate(n,TeV);
//				tab_keX[bb][t] = kf*exp(-x_lu);
//				tab_keD[bb][t] = grat*kf;
//				if (tab_keX[bb][t] < 1.e-200)
//				{
//					tab_keX[bb][t] = 0.;//1.e-70;
//					tab_keD[bb][t] = 0.;//1.e-70;
//				}

				double omega_H = 13.6/TeV;
				double ve_bar = 2.468e+06/sqrt(omega_H);
				const int N = 1;
				double res, err, tol = 1.e-3;
				double xl[N] = { 0. };
				double xu[N] = { 1. };
				params_exc par;
				par.T = TeV;
				par.tran = bbexd[bb];
				hcubature(1, f_exc_ion_transform, &par, N, xl, xu, 0, 0, tol, ERROR_INDIVIDUAL, &res, &err);
				tab_keX[bb][t] = ve_bar*res*exp(-x_lu);
				tab_keD[bb][t] = grat*ve_bar*res;
//				if (tab_keX[bb][t] < 1.e-200)
//				{
//					tab_keX[bb][t] = 0.;
//					tab_keD[bb][t] = 0.;
//				}

//				if (isnan(tab_keX[bb][t]))
//				{
//					printf("na = %d\n",na);
//					printf("n = %d\n",n);
//					printf("%d -> %d\n",l,u);
//					printf("TeV = %f\n",TeV);
//					printf("dE = %f\n",bbexd[bb].dE);
//					printf("xs\n");
//					for (int pp=0;pp<bbexd[bb].NEpts;pp++)
//					{
//						printf("%e \t %e\n",A.bbeexxs[n][pp][1],bbexd[bb].xs_dat[pp]);
//					}
//				}
			}

			bb++;
		}
	}

	// Print rate calculations
//	FILE* fp;

    printf("%d Collisional X/D\n",NBBeXD);
//    printf("  BB       lower   -    upper       dE[eV]     keX(%6.3f eV)\n",exp(tab_lnT[NTstep/2]));
//    for(bb=0;bb<NBBeXD;bb++)
//    {
//        int atom_index = bbexd[bb].low[0];
//        if (atom_index != 0) continue;
//        AtomX* sp = atom[atom_index];
//        l = bbexd[bb].low[1];
//        if (l != 0) continue;
//        u = bbexd[bb].upp[1];
//        if (u != 1 && u != 5 && u != 15 && u != 25 && u != 30) continue;
//        float dE = sp->Ee[u]-sp->Ee[l];//bbexd[bb].dE;
//        printf("%4d     (%2d,%4d) -  (%2d,%4d)     %6.3f     %12.5e     %12.5e\n",bb,atom_index,l,atom_index,u,dE,tab_keX[bb][NTstep/2],tab_keD[bb][NTstep/2]);
//    }

//    string ratexfil;
//    string ratexdir = "./Results/Rates/CollEx";
//    system("exec rm -r ./Results/Rates/CollEx*");
//    for(bb=0;bb<NBBeXD;bb++)
//    {
//        int atom_index = bbexd[bb].low[0];
//        if (atom_index != 0) continue;
//        AtomX* sp = atom[atom_index];
//        l = bbexd[bb].low[1];
//        if (l != 0) continue;
//        u = bbexd[bb].upp[1];
//        if (u != 1 && u != 5 && u != 15 && u != 25 && u != 30) continue;
//        float dE = sp->Ee[u]-sp->Ee[l];//bbexd[bb].dE;
//
//        stringstream convert;
//        convert << ratexdir << "_" << atom_index << "_" <<  l << "_" << atom_index << "_" <<  u << ".dat";
//        ratexfil = convert.str();
//
//		fp = fopen(ratexfil.c_str(),"w");
//
//		for(t=0;t<=NTstep;t++)
//		{
//			double TeV = exp(tab_lnT[t]);
//			fprintf(fp,"%6.3f     %12.5e     %12.5e\n",TeV,tab_keX[bb][t],tab_keD[bb][t]);
//		}
//
//		fflush(fp);
//    }

    printf("%d Atomic Lines\n\n",NBBlns);
//	printf("  BB       lower   -    upper       dE[eV]     Arad[s-1]\n",exp(tab_lnT[NTstep]));
//	for(bb=0;bb<NBBlns;bb++)
//	{
//		int atom_index = bbexd[bb].low[0];
//		AtomX* sp = atom[atom_index];
//		l = bbrad[bb].low[1];
//		u = bbrad[bb].upp[1];
//		float dE = sp->Ee[u]-sp->Ee[l];//bbrad[bb].dE;
//		printf("%4d     (%2d,%4d) -  (%2d,%4d)     %6.3f     %12.5e\n",bb,atom_index,l,atom_index,u,dE,rad_Abb[bb]);
//	}

//	fclose(fp);
}

/**
 * Ionization/Recombination
 */
void KRates::add_BFtrans(bool irrad)
{
	int n,t,l,u,na,ia,bf;

	/*
	 * Electron-impact ionization/recombination
	 */
	printf("Tabulating CIR BF rates:\n");

	// Count transitions
	NBFeIR=0;
	for(na=0;na<Natoms;na++)
	{
		AtomX* neu = atom[na];
		AtomX* ion = neu->ion;  // ion (upper ion stage)
		if(!ion) continue;
		ia = -1;
		// if ion exists, which is it in the list?
		for(int k=0;k<Natoms;k++)
		{   // select entry based on name comparison
			if(ion->name == atom[k]->name) { ia = k; break; }
		}
		if(ia < 0)
		{
			printf("could not find ion (%s) in list\n",ion->name.c_str()); abort();
		}
		int Nnlev = Nalevs[na];
		int Nilev = Nalevs[ia];

		for(n=0;n<neu->NBFeIR;n++)
		{
			l = neu->bfeir[n].low[1];
			u = neu->bfeir[n].upp[1];
			if (l >= Nnlev || u >=Nilev ) continue;
			NBFeIR++;
		}
	}

    // Initialize indexing arrays
	bfeir = new transition[NBFeIR];
	tab_keI = cArray<double>(NBFeIR,NTstep+1);
	tab_keR = cArray<double>(NBFeIR,NTstep+1);

	bf=0;
	for(na=0;na<Natoms;na++)
	{
		AtomX* neu = atom[na];
		printf("... Atom %s\n",neu->name.c_str());
		AtomX* ion = neu->ion;  // ion (upper ion stage)
		if(!ion) continue;
		ia = -1;
		// if ion exists, which is it in the list?
		for(int k=0;k<Natoms;k++)
		{   // select entry based on name comparison
			if(ion->name == atom[k]->name) { ia = k; break; }
		}
		if(ia < 0)
		{
			printf("could not find ion (%s) in list\n",ion->name.c_str()); abort();
		}
		int Nnlev = Nalevs[na];
		int Nilev = Nalevs[ia];

		for(n=0;n<neu->NBFeIR;n++)
		{
			// Only store relevant level transitions
			l = neu->bfeir[n].low[1];
			u = neu->bfeir[n].upp[1];
			if (l >= Nnlev || u >=Nilev ) continue;

			bfeir[bf].low[0] = na;
			bfeir[bf].upp[0] = ia;
			bfeir[bf].low[1] = l;
			bfeir[bf].upp[1] = u;
			bfeir[bf].low[2] = neu->Z;
			bfeir[bf].upp[2] = ion->Z;
			bfeir[bf].dE     = neu->bfeir[n].dE;

			bfeir[bf].NEpts = neu->bfeir[n].NEpts;
			bfeir[bf].lnEmin= neu->bfeir[n].lnEmin;
			bfeir[bf].lnEmax= neu->bfeir[n].lnEmax;
			bfeir[bf].dlnE  = neu->bfeir[n].dlnE;
			bfeir[bf].xs_dat= neu->bfeir[n].xs_dat;

			for(t=0;t<=NTstep;t++)
			{
				double TeV = exp(tab_lnT[t]);
				double x_lu= bfeir[bf].dE/TeV;
				double grat= (double)neu->ge[l]/(double)ion->ge[u];
				double T3o2= TeV*sqrt(TeV);

//				double kf = neu->eion_rate(n,TeV);
//				tab_keI[bf][t] = kf*exp(-x_lu);
//				tab_keR[bf][t] = grat/T3o2/6.e27*kf;
//				if (tab_keI[bf][t] < 1.e-200)
//				{
//					tab_keI[bf][t] = 0.;//1.e-70;
//					tab_keR[bf][t] = 0.;//1.e-70;
//				}

				double omega_H = 13.6/TeV;
				double ve_bar = 2.468e+06/sqrt(omega_H);
				const int N = 1;
				double res, err, tol = 1.e-3;
				double xl[N] = { 0. };
				double xu[N] = { 1. };
				params_exc par;
				par.T = TeV;
				par.tran = bfeir[bf];
				hcubature(1, f_exc_ion_transform, &par, N, xl, xu, 0, 0, tol, ERROR_INDIVIDUAL, &res, &err);
				tab_keI[bf][t] = ve_bar*res*exp(-x_lu);
				tab_keR[bf][t] = grat/T3o2/6.e27*ve_bar*res;
//				if (tab_keI[bf][t] < 1.e-200)
//				{
//					tab_keI[bf][t] = 0.;
//					tab_keR[bf][t] = 0.;
//				}
			}

			bf++;
		}
	}

	printf("Tabulating PIR BF rates:\n");

	// Count transitions
	NBFpIR=0;
	for(na=0;na<Natoms;na++)
	{
		AtomX* neu = atom[na];
		printf("... Atom %s\n",neu->name.c_str());
		AtomX* ion = neu->ion;  // ion (upper ion stage)
		if(!ion) continue;
		ia = -1;
		// if ion exists, which is it in the list?
		for(int k=0;k<Natoms;k++)
		{   // select entry based on name comparison
			if(ion->name == atom[k]->name) { ia = k; break; }
		}
		if(ia < 0)
		{
			printf("could not find ion (%s) in list\n",ion->name.c_str()); abort();
		}
		int Nnlev = Nalevs[na];
		int Nilev = Nalevs[ia];

		for(n=0;n<neu->NBFpIR;n++)
		{
			l = neu->bfpir[n].low[1];
			u = neu->bfpir[n].upp[1];
			if (l >= Nnlev || u >=Nilev ) continue;
			NBFpIR++;
		}
	}

    // Initialize indexing arrays
	bfpir = new transition[NBFpIR];
	if (irrad) {
		tab_kpI = cArray<double*>(NBFpIR);
		tab_kpS = cArray<double*>(NBFpIR,NTstep+1);
	}
	tab_kpR = cArray<double*>(NBFpIR,NTstep+1);

	bf=0;

//	char filename[150];
//	std::ofstream fileBolt;   // file stream
////	sprintf(filename,"./Rates/");
////	sprintf(filename + strlen(filename),"PhotoTransition.txt");
////	fileBolt.open (filename);  // open the file in filenames at index i
////	fileBolt << "Transition,\tRate\n";

	for(na=0;na<Natoms;na++)
	{
		AtomX* neu = atom[na];
		AtomX* ion = neu->ion;  // ion (upper ion stage)
		if(!ion) continue;
		ia = -1;
		// if ion exists, which is it in the list?
		for(int k=0;k<Natoms;k++)
		{   // select entry based on name comparison
			if(ion->name == atom[k]->name) { ia = k; break; }
		}
		if(ia < 0)
		{
			printf("could not find ion (%s) in list\n",ion->name.c_str()); abort();
		}
		int Nnlev = Nalevs[na];
		int Nilev = Nalevs[ia];

		for(n=0;n<neu->NBFpIR;n++)
		{
			// Only store relevant level transitions
			l = neu->bfpir[n].low[1];
			u = neu->bfpir[n].upp[1];
			if (l >= Nnlev || u >=Nilev ) continue;

			bfpir[bf].low[0] = na;
			bfpir[bf].upp[0] = ia;
			bfpir[bf].low[1] = l;
			bfpir[bf].upp[1] = u;
			bfpir[bf].low[2] = neu->Z;
			bfpir[bf].upp[2] = ion->Z;
			bfpir[bf].dE     = neu->bfpir[n].dE;

			bfpir[bf].NEpts = neu->bfpir[n].NEpts;
			bfpir[bf].lnEmin= neu->bfpir[n].lnEmin;
			bfpir[bf].lnEmax= neu->bfpir[n].lnEmax;
			bfpir[bf].dlnE  = neu->bfpir[n].dlnE;
			bfpir[bf].xs_dat= neu->bfpir[n].xs_dat;

			const int N = 1;
			double res[2], err[2], tol = 1.e-3;
			double xl[N] = { 0. };
			double xu[N] = { 1 };

			if (irrad)
			{
				tab_kpI[bf] = new double[2];
				double pionC = 3.9519911e27*Trad*Trad*Trad;		// 8*pi/(h^3*c^2)
				params_exc par;
				par.Tr = Trad;
				par.tran = bfpir[bf];
				par.A = atom;
				hcubature(2, f_pho_ion_transform, &par, N, xl, xu, 0, 0, tol, ERROR_INDIVIDUAL, res, err);
				tab_kpI[bf][0] = pionC*res[0];
				tab_kpI[bf][1] = pionC*res[1];

//				if (na == 0 && ia == 1 && l == 0 && u == 0)
//					printf("Photoionization rate: %3.3e and %3.3e\n",pionC*res[0],pionC*res[1]);

//				fileBolt << na << ",\t" << n << ",\t" << tab_kpI[bf][0] << ",\t" << tab_kpI[bf][1] <<"\n";

				if (isnan(tab_kpI[bf][0]))
				{
					printf("!!!!!!!!nan rate\n");
					printf("A(%d,%d) -> A(%d,%d)\n",ion->Z,u,neu->Z,l);
					printf("Tr = %e\n",Trad);

					for (int pp=0;pp<bfpir[bf].NEpts;pp++)
					{
						printf("%e \t %e\n",neu->bfpioxs[n][pp][0],bfpir[bf].xs_dat[pp]);
					}
				}
			}

//			filename[0] = 0;
//			sprintf(filename,"./Rates/");
//			sprintf(filename + strlen(filename),"%02dRadRecTransition%04d.txt",na,n);
//			fileBolt.open (filename);  // open the file in filenames at index i
//			fileBolt << "Temperature,\tRateRecN,\tRateRecE,\tRateStiN,\tRateStiE\n";

			for(t=0;t<=NTstep;t++)
			{
				double TeV = exp(tab_lnT[t]); 				// tab_kpR[bf][t] = neu->prec_rate(n,TeV);
				double omega_H = 13.6/TeV;
				double ve_bar = 2.468e+06/sqrt(omega_H);	// mean thermal velocity sqrt(8*kT/pi*mu)	[Refer to the formulary]
				params_exc par;
				par.T = TeV;
				par.tran = bfpir[bf];
				par.A = atom;

				if (irrad)
				{
					tab_kpS[bf][t] = new double[2];
					par.Tr = Trad;
					hcubature(2, f_stim_rec_transform, &par, N, xl, xu, 0, 0, tol, ERROR_INDIVIDUAL, res, err);
					tab_kpS[bf][t][0] = ve_bar*res[0];
					tab_kpS[bf][t][1] = ve_bar*res[1];
				}

				tab_kpR[bf][t] = new double[2];
				hcubature(2, f_rad_rec_transform, &par, N, xl, xu, 0, 0, tol, ERROR_INDIVIDUAL, res, err);
				tab_kpR[bf][t][0] = ve_bar*res[0];
				tab_kpR[bf][t][1] = ve_bar*res[1];

//				if (na == 0 && ia == 1 && l == 0 && u == 0)
//					printf("Recombination rates at Te = %3.3e: %3.3e and %3.3e\n",TeV,ve_bar*res[0],ve_bar*res[1]);

//				fileBolt << TeV << ",\t" << tab_kpR[bf][t][0] << ",\t" << tab_kpR[bf][t][1] << ",\t" << tab_kpS[bf][t][0] << ",\t" << tab_kpS[bf][t][1] << "\n";

				if (isnan(tab_kpR[bf][t][0]) || tab_kpR[bf][t][0] < 0.)
				{
					printf("!!!!!!!!nan rate\n");
					printf("A(%d,%d) -> A(%d,%d)\n",ion->Z,u,neu->Z,l);
					printf("T = %e\n",TeV);

					for (int pp=0;pp<bfpir[bf].NEpts;pp++)
					{
						printf("%e \t %e\n",neu->bfpioxs[n][pp][0],bfpir[bf].xs_dat[pp]);
					}
				}
			}

//			fileBolt.close();

			bf++;
		}
	}

////	fileBolt.close();
//	exit(1);

	// Print rate calculations
//	FILE* fp;

    printf("%d Collisional I/R\n",NBFeIR);
//    printf("  BF       lower   -    upper        I[eV]     keI(%6.3f eV)\n",exp(tab_lnT[NTstep/2]));
//    for(bf=0;bf<NBFeIR;bf++)
//    {
//        int neu_index = bfeir[bf].low[0];
//        if (neu_index != 0) continue;
//        int ion_index = bfeir[bf].upp[0];
//        if (ion_index != 1) continue;
//        AtomX* neu = atom[neu_index];
//        AtomX* ion = atom[ion_index];
//        int n = bfeir[bf].low[1];
//        if (n != 0) continue;
//        int m = bfeir[bf].upp[1];
//        if (m != 0 && m != 1 && m != 2 && m != 58) continue;
//        float dE = bfeir[bf].dE;
//        printf("%4d     (%2d,%4d) -  (%2d,%4d)     %6.3f     %12.5e     %12.5e\n",bf,neu_index,n,ion_index,m,dE,tab_keI[bf][NTstep/2],tab_keR[bf][NTstep/2]);
//    }

//    string ratiofil;
//	string ratiodir = "./Results/Rates/CollIo";
//	system("exec rm -r ./Results/Rates/CollIo*");
//	for(bf=0;bf<NBFeIR;bf++)
//	{
//		int neu_index = bfeir[bf].low[0];
//		if (neu_index != 0) continue;
//		int ion_index = bfeir[bf].upp[0];
//		if (ion_index != 1) continue;
//		AtomX* neu = atom[neu_index];
//		AtomX* ion = atom[ion_index];
//		int n = bfeir[bf].low[1];
//		if (n != 0) continue;
//		int m = bfeir[bf].upp[1];
//		if (m != 0 && m != 1 && m != 2 && m != 58) continue;
//		float dE = bfeir[bf].dE;
//
//		stringstream convert;
//		convert << ratiodir << "_" << neu_index << "_" <<  n << "_" << ion_index << "_" <<  m << ".dat";
//		ratiofil = convert.str();
//
//		fp = fopen(ratiofil.c_str(),"w");
//
//		for(t=0;t<=NTstep;t++)
//		{
//			double TeV = exp(tab_lnT[t]);
//			fprintf(fp,"%6.3f     %12.5e     %12.5e\n",TeV,tab_keI[bf][t],tab_keR[bf][t]);
//		}
//
//		fflush(fp);
//	}

    printf("%d Radiative I/R\n",NBFpIR);
//    printf("  BF       lower   -    upper        I[eV]     keI(%6.3f eV)\n",exp(tab_lnT[NTstep/2]));
//    for(bf=0;bf<NBFeIR;bf++)
//    {
//        int neu_index = bfeir[bf].low[0];
//        int ion_index = bfeir[bf].upp[0];
//        AtomX* neu = atom[neu_index];
//        AtomX* ion = atom[ion_index];
//        int n = bfeir[bf].low[1];
//        int m = bfeir[bf].upp[1];
//        float dE = bfeir[bf].dE;
//        printf("%4d     (%2d,%4d) -  (%2d,%4d)     %6.3f     %12.5e     %12.5e\n",bf,neu_index,n,ion_index,m,dE,tab_keI[bf][NTstep/2],tab_keR[bf][NTstep/2]);
//    }

//    fclose(fp);

    printf("\n");
}

/**
 * Autoionization/Electron capture
 */
void KRates::add_AItrans()
{
	int n,m,t,l,u,na,ia,bf;

	/*
	 * Autoionization rates
	 */
	printf("Tabulating AI/EC rates\n");

	// count transition
	NBFaIR=0;
	for(na=0;na<Natoms;na++)
	{
		AtomX* neu = atom[na];
		AtomX* ion = neu->ion;  // ion (upper ion stage)
		if(!ion) continue;
		ia = -1;
		// if ion exists, which is it in the list?
		for(int k=0;k<Natoms;k++)
		{   // select entry based on name comparison
			if(ion->name == atom[k]->name) { ia = k; break; }
		}
		if(ia < 0)
		{
			printf("could not find ion (%s) in list\n",ion->name.c_str()); abort();
		}
		int Nnlev = Nalevs[na];
		int Nilev = Nalevs[ia];

		for(n=0;n<neu->NBFaIR;n++)
		{
			l = neu->bfair[n].low[1];
			u = neu->bfair[n].upp[1];
			if (l >= Nnlev || u >=Nilev ) continue;
			NBFaIR++;
		}
	}

    // Initialize indexing arrays
	bfair = new transition[NBFaIR];
	kAI = cArray<double>(NBFaIR);
	tab_kDR = cArray<double>(NBFaIR,NTstep+1);

	bf=0;
	for(na=0;na<Natoms;na++)
	{
		AtomX* neu = atom[na];
		printf("... Atom %s\n",neu->name.c_str());
		AtomX* ion = neu->ion;  // ion (upper ion stage)
		if(!ion) continue;
		ia = -1;
		// if ion exists, which is it in the list?
		for(int k=0;k<Natoms;k++)
		{   // select entry based on name comparison
			if(ion->name == atom[k]->name) { ia = k; break; }
		}
		if(ia < 0)
		{
			printf("could not find ion (%s) in list\n",ion->name.c_str()); abort();
		}
		int Nnlev = Nalevs[na];
		int Nilev = Nalevs[ia];

		for(n=0;n<neu->NBFaIR;n++)
		{
			l = neu->bfair[n].low[1];
			u = neu->bfair[n].upp[1];
			if (l >= Nnlev || u >=Nilev ) continue;

			bfair[bf].low[0] = na;
			bfair[bf].upp[0] = ia;
			bfair[bf].low[1] = l;
			bfair[bf].upp[1] = u;
			bfair[bf].low[2] = neu->Z;
			bfair[bf].upp[2] = ion->Z;
			bfair[bf].dE     = neu->bfair[n].dE;

			kAI[bf] = neu->bfaiort[n];

			for(t=0;t<=NTstep;t++)
			{
				double TeV = exp(tab_lnT[t]);
				double x_lu= bfair[bf].dE/TeV;
				double grat= (double)neu->ge[l]/(double)ion->ge[u];
				double T3o2= TeV*sqrt(TeV);

				tab_kDR[bf][t] = neu->bfaiort[n]*grat/T3o2/6.e27*exp(-x_lu);
			}

			bf++;
		}
	}


    printf("%d Autoionization/Electron capture I/R\n\n",NBFaIR);
}

/**
 * Elastic Collisions
 */
//void KRates::add_ECexchg()
//{
//	int n,m,t,l,u,na,bf;
//
//	/*
//	 * Elastic Collisions
//	 */
//
//    // Initialize indexing arrays
//	tab_kEC = new double[NTstep+1];
//	for(t=0;t<=NTstep;t++)
//	{
//		tab_kEC[t] = 0.;
//	}
//
//	AtomX* neu = atom[0];
//
//	for(t=0;t<=NTstep;t++)
//	{
//		double TeV = exp(tab_lnT[t]);
//		tab_kEC[t] = neu->ecol_rate(TeV);
//	}
//
//    printf("Elastic Collision rates calculated\n\n");
//
////    FILE* fp;
////    fp = fopen("./Results/Rates/ECRates.dat","w");
////	for(t=0;t<=NTstep;t++)
////	{
////		double TeV = exp(tab_lnT[t]);
////		fprintf(fp,"%3.3e",TeV);
////		fprintf(fp,"\t%3.3e",tab_kEC[t]);
////		fprintf(fp,"\n");
////	}
////	fclose(fp);
//}

/**
 * Determine whether states are completely filled
 */
void KRates::check()
{
	int xcnt[Nstate];
	for(int n=0;n<Nstate;n++)
		xcnt[n]=0;
	for(int bb=0;bb<NBBeXD;bb++)
	{
		int na = bbexd[bb].low[0];
		int l = bbexd[bb].low[1];
		int u = bbexd[bb].upp[1];
		int ql = state[na][l];
		int qu = state[na][u];
		xcnt[ql]++;
		xcnt[qu]++;
	}
	for(int bf=0;bf<NBFeIR;bf++)
	{
		int na = bfeir[bf].low[0];
		int l  = bfeir[bf].low[1];
		int ia = bfeir[bf].upp[0];
		int u  = bfeir[bf].upp[1];
		int ql = state[na][l];
		int qu = state[ia][u];
		xcnt[ql]++;
		xcnt[qu]++;
	}
	for(int n=0;n<Nstate;n++)
	{
		if (xcnt[n]==0) printf("state %d inactive\n",n);
	}
}

/**
 * Determine number of states
 */
void KRates::DistPrep(string GrpType)
{
	printf("Preparing distribution bins...\n\n");

	int totStates = 0;

	if (GrpType == "Full") {
		for(int n=0;n<Natoms;n++)
			totStates += Nalevs[n];
	}
	else {
		printf("Not a valid grouping option: %s",GrpType.c_str());
		abort();
	}

	Ncr_init = new double[totStates];
}

/**
 * Distribute densities across all atomic states
 */
void KRates::AtomDist(double ae, double Na, double TeV)
{
	FILE* btFile;

	const char* dirString = "./Results/Current_Run/";
	const char* boltzInit = "BoltzInit.txt";

	char currDir[100];

	strcpy(currDir,dirString);
	strcat(currDir,boltzInit);
	btFile = fopen(currDir,"w"); // File opened: Boltzmann distributions are plotted
	currDir[0] = '\0';
	fprintf(btFile,"Energy Level \t ln(N/g)\n");

	double NLv1tot = (1.-ae)*Na;
	double NLv2tot = ae*Na;
	double Qpart = 0.;

	for(int n=0;n<Natoms;n++) {
		// Evenly distribute electrons to ionized stages
		int numLvs = Nalevs[n];
		double Za = atom[n]->Z;

		for (int m=0;m<numLvs;m++) {
			double gval = atom[n]->ge[m];
			double Eval = atom[n]->Ee[m];

			Qpart += gval*exp(-Eval/TeV);
		}

		// Theoretical and Density Assignment
		for (int m=0;m<numLvs;m++) {
			int ka   	= state[n][m];
			double gval = atom[n]->ge[m];
			double Eval = atom[n]->Ee[m];

			double NLv = 1.;
			if (Za == 0.)
				NLv = NLv1tot;
			else if (Za == 1.)
				NLv = NLv2tot;
			else
				NLv = 1.e-5;

			double nAssign = NLv*gval*exp(-Eval/TeV)/Qpart;
			double logNg = log10(nAssign/gval);

			if (logNg < 0.)		// Value chosen for 'loosened'/'unstiff' initial iteration
				Ncr_init[ka] = 5.e-1*gval;
			else
				Ncr_init[ka] += nAssign;

			fprintf(btFile,"%3.3e \t %3.3e \t %3.3e \t %3.3e \t %3.3e\n", \
				Eval,logNg,log10(Ncr_init[ka]/gval),gval,Qpart);
		}
		fprintf(btFile,"\n");

		Qpart = 0.;
	}
	fflush(btFile);
	fclose(btFile);

	printf("Applied atomic density distribution...\n\n");
}

/*
 * Solve for the steady state solution
 */
double* KRates::solve_steady_state(double Ne, double TeV, bool is_radiating)
{
	int msize = Nstate;
	double *m_data =  cArray<double>(msize*msize);
	double *b_data  = cArray<double>(msize);
	int *Zq  = cArray<int>(msize);

	int n,m,na;
	int low,upp;


	// find interpolation point in tables
	double lnTeV = max(lnT_min,min(.99*lnT_max,log(TeV) ) );
	float xt = (lnTeV-lnT_min)/dlnT;
	int it = (int) xt; xt -= it;

	// add-up excitations and deexcitations
	for (int bb = 0 ; bb < NBBeXD; bb++)
	{
		na = bbexd[bb].low[0];	// Determine ionic-atom stage
		n  = bbexd[bb].low[1]; // Lower level
		m  = bbexd[bb].upp[1]; // Upper level
		low = state[na][n]; // Find corresponding state for lower level
		upp = state[na][m]; // Find corresponding state for upper level
		//if(m != n+1) continue;
		// excitation (n -> m)
		double a_nm = (1.-xt) * tab_keX[bb][it] + xt * tab_keX[bb][it+1];
//		printf("a(%d -> %d) = %e\n",low,upp,a_nm);
		m_data[upp*msize+low] += a_nm*Ne;
		b_data[low]           += a_nm*Ne;

		// deexcitation (m -> n)
		double b_nm = (1.-xt) * tab_keD[bb][it] + xt * tab_keD[bb][it+1];
//		printf("b(%d <- %d) = %e\n",low,upp,b_nm);
		m_data[low*msize+upp] += b_nm*Ne;
		b_data[upp]           += b_nm*Ne;
	}

	// add-up ionizations and recombinations
	for(int bf=0;bf<NBFeIR;bf++)
	{
		na = bfeir[bf].low[0]; // neutral atom index
		ia = bfeir[bf].upp[0]; // ionized atom index
		n  = bfeir[bf].low[1]; // lower level (initial)
		m  = bfeir[bf].upp[1]; // upper level
		low = state[na][n]; // find corresponding state
		upp = state[ia][m]; // find corresponding state
		// ionization (n -> m)
		double a_nm = (1.-xt) * tab_keI[bf][it] + xt * tab_keI[bf][it+1];
//		printf("a(%d -> %d) = %e\n",low,upp,a_nm);
//		a_nm = max(1.e-60,a_nm);
		m_data[upp*msize+low] += a_nm*Ne;
		b_data[low]           += a_nm*Ne;

		// recombination (m -> n)
		double b_nm = (1.-xt) * tab_keR[bf][it] + xt * tab_keR[bf][it+1];
//		printf("b(%d <- %d) = %e\n",low,upp,b_nm);
//		b_nm = max(1.e-60,b_nm);
		m_data[low*msize+upp] += b_nm*Ne*Ne;
		b_data[upp]           += b_nm*Ne*Ne;
	}

	// add-up autoionization/electron capture
	for(int bf=0;bf<NBFaIR;bf++)
	{
		na = bfair[bf].low[0]; // neutral atom index
		ia = bfair[bf].upp[0]; // ionized atom index
		n  = bfair[bf].low[1]; // lower level (initial)
		m  = bfair[bf].upp[1]; // upper level
		low = state[na][n]; // find corresponding state
		upp = state[ia][m]; // find corresponding state
		// ionization (n -> m)
		double a_nm = kAI[bf];
		m_data[upp*msize+low] += a_nm;
		b_data[low]           += a_nm;

		// recombination (m -> n)
		double x_lu= bfair[bf].dE/TeV;
		double grat= (double)atom[na]->ge[n]/(double)atom[ia]->ge[m];
		double T3o2= TeV*sqrt(TeV);
		double b_nm = grat/T3o2/6.e27*a_nm*exp(-x_lu);
		m_data[low*msize+upp] += b_nm*Ne;
		b_data[upp]           += b_nm*Ne;
	}

	if (is_radiating) {
		for(int bb=0;bb<NBBlns;bb++)
		{
			na = bbrad[bb].low[0];
			n  = bbrad[bb].low[1]; // lower level
			m  = bbrad[bb].upp[1]; // upper level
			low =state[na][n]; // find corresponding state
			upp =state[na][m]; // find corresponding state
			// spontaneous emission (m -> n)
			double A_nm = rad_Abb[bb];
			m_data[low*msize+upp] += A_nm;
			b_data[upp]           += A_nm;
		}
		// add-up radiative (BF) transitions
		for(int bf=0;bf<NBFpIR;bf++)
		{
			na = bfpir[bf].low[0]; // neutral atom index
			ia = bfpir[bf].upp[0]; // ionized atom index
			n  = bfpir[bf].low[1]; // lower level (initial)
			m  = bfpir[bf].upp[1]; // upper level
			low = state[na][n]; // find corresponding state
			upp = state[ia][m]; // find corresponding state

			// recombination (m -> n)
			double b_nm = (1.-xt) * tab_kpR[bf][it][0] + xt * tab_kpR[bf][it+1][0];
			m_data[low*msize+upp] += b_nm*Ne;
			b_data[upp]           += b_nm*Ne;
		}
	}

	for(int i=0;i<msize;i++) m_data[i*msize+i] -= b_data[i];
	for(int i=0;i<msize;i++) b_data[i] = 0.;

	// replace the last row for ion
	for(int a=0;a<Natoms;a++)
	{
		int Nlev = Nalevs[a];
		for(int n=0;n<Nlev;n++) Zq[state[a][n]] = atom[a]->Z;
	}
	int mmax=msize-1;
	for(int i=0;i<msize;i++) m_data[mmax*msize+i] = Zq[i];
	b_data[mmax] = Ne;

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
	return b_data;
}
