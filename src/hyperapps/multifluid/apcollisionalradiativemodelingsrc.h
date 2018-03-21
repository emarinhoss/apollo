#ifndef APCOLLISIONALRADIATIVEMODELINGSRC_H
#define APCOLLISIONALRADIATIVEMODELINGSRC_H

// WarpX includes
#include <wxhyperbolicsrc.h>

// std includes
#include <omp.h>
#include <iostream>
#include <string>
#include <petsc.h>
#include <math.h>
#include <stdio.h>

// CR includes
#include "myIO.h"
#include "MathX.h"
#include "Constants.h"
#include "CRParams.h"
#include "KRatesGrp.h"
#include "CRSolve.h"
//#include "LSpecs.h"
//#include "Spectrum.h"
//#include "StiffIntegratorT.h"
//#include "ColRad.h"

template<class REAL>
class ApCollisionalRadiativeModelingSrc : public WxHyperbolicSrc<REAL>
{
  public:

    ApCollisionalRadiativeModelingSrc()
        : WxHyperbolicSrc<REAL>("CRmodel") {
    }

    void setup(const WxCryptSet& wxc)
    {
        // call base-class setup function
        WxHyperbolicSrc<REAL>::setup(wxc);

        // read charge and mass of particles
        _q = wxc.template get<REAL>("charge");
        _mi1 = wxc.template get<REAL>("ion_mass");
        _mi2 = wxc.template get<REAL>("neutral_mass");
        _me = wxc.template get<REAL>("electron_mass");

        _gas_gamma = wxc.template get<REAL>("gas_gamma");

        // CR folder to be used
        std::string fname = wxc.template get<std::string>("CR_location");

        jSys = new myiolib::jFileSys(fname);

        // Kinetic Rates
        pars.Natoms 	= 2;
        if (wxc.has("Natoms"))
          pars.Natoms = wxc.template get<int>("Natoms");

        pars.temp_bins 	= 90;
        if (wxc.has("temp_bins"))
          pars.temp_bins = wxc.template get<int>("temp_bins");

        pars.dens_bins 	= 51;
        if (wxc.has("dens_bins"))
          pars.temp_bins = wxc.template get<int>("dens_bins");

        pars.temp_min 	= 0.03;
        if (wxc.has("temp_min"))
          pars.temp_min = wxc.template get<REAL>("temp_min");

        pars.temp_max 	= 1000.;
        if (wxc.has("temp_max"))
          pars.temp_max = wxc.template get<REAL>("temp_max");

        pars.dens_min 	= 1.e9;
        if (wxc.has("dens_min"))
          pars.dens_min = wxc.template get<REAL>("dens_min");

        pars.dens_max 	= 1.e27;
        if (wxc.has("dens_max"))
          pars.dens_min = wxc.template get<REAL>("dens_max");

        // CR Solver Configuration
        pars.irrad_bool 	= false;
//        if (wxc.has("irrad_bool"))
//          pars.irrad_bool = wxc.template get<bool>("irrad_bool");

        pars.rad_bool       = false;
//        if (wxc.has("rad_bool"))
//          pars.rad_bool = wxc.template get<bool>("rad_bool");

        pars.fxdNe_bool 	= false;
//        if (wxc.has("fxdNe_bool"))
//          pars.fxdNe_bool = wxc.template get<bool>("fxdNe_bool");

        pars.intNe_bool 	= false;
//        if (wxc.has("intNe_bool"))
//          pars.intNe_bool = wxc.template get<bool>("intNe_bool");

        pars.intTh_bool 	= false;
//        if (wxc.has("intTh_bool"))
//          pars.intTh_bool = wxc.template get<bool>("intTh_bool");

        pars.intTe_bool 	= false;
//        if (wxc.has("intTe_bool"))
//          pars.intTe_bool = wxc.template get<bool>("intTe_bool");

        // Temporal Parametrization
        pars.multxdt 	= 1.;
        if (wxc.has("multxdt"))
          pars.multxdt = wxc.template get<REAL>("multxdt");

        // Initial Conditions
        pars.temp_rad 	= 300.;
        if (wxc.has("temp_rad"))
          pars.temp_rad = wxc.template get<REAL>("temp_rad");

        // Output Conditions
        pars.print_log 	= false;
        pars.printdt 	= 1.;
        pars.logdtprt 	= 1.;

        // Spectra Calculation and Parameters
        pars.print_spectra = false;
//        if (wxc.has("print_spectra"))
//          pars.print_spectra = wxc.template get<bool>("print_spectra");

        pars.wave_min 	= 20.e-9;
        if (wxc.has("wave_min"))
          pars.wave_min = wxc.template get<REAL>("wave_min");

        pars.wave_max 	= 600.e-9;
        if (wxc.has("wave_max"))
          pars.wave_max = wxc.template get<REAL>("wave_max");

        pars.wave_resol = 1.e-12;
        if (wxc.has("wave_resol"))
          pars.wave_resol = wxc.template get<REAL>("wave_resol");

        // Tolerances
        pars.rel_tol 	= 1.e-3;
        if (wxc.has("rel_tol"))
          pars.rel_tol = wxc.template get<REAL>("rel_tol");

        pars.abs_tol 	= 1.e9;
        if (wxc.has("abs_tol"))
          pars.abs_tol = wxc.template get<REAL>("abs_tol");

        pars.atom_cutoffs = NULL;

        pars.p2s 		= 0;

        pars.grouping_type = "QSS";
        pars.atom_sym = "Ar";
        pars.output_dir = "./CR_Run";

        jSys->display();

        // Define participating atomic models
        AtomX** Atoms = new AtomX*[pars.Natoms];
        for (int ida=0; ida<pars.Natoms; ida++) {
            stringstream itostr;
            itostr << ida;
            string filenum = itostr.str();
            string astr = pars.atom_sym + "+" + filenum;
            Atoms[ida] = (new AtomX(astr));
            Atoms[ida]->AtomParams(pars.atom_sym);
            printf("\n");
        }

        for (int ida=0; ida<pars.Natoms-1; ida++)
            Atoms[ida]->set_ion(Atoms[ida+1]);
        Atoms[pars.Natoms-1]->ion = NULL;

        // Define participating kinetic processes
        _kin =  KRatesGrp(pars.atom_sym.c_str(),pars.Natoms);
        for (int ida=0;ida<pars.Natoms;ida++)
            _kin.addAtom(Atoms[ida],pars.atom_cutoffs[ida]);
        delete [] Atoms;
        _kin.DistPrep(pars.grouping_type);
        _kin.setRange(pars.temp_bins,pars.temp_min,pars.temp_max);
        _kin.setNRange(pars.dens_bins,pars.dens_min,pars.dens_max);
        _kin.setTPlanck(pars.temp_rad);
        // Add all species to chemistry set and atomic grouping/distributions
        _kin.add_BBtrans(pars.irrad_bool);  // Add all bound-bound transitions ======> Output # of Transitions
        _kin.add_BFtrans(pars.irrad_bool);  // Add all bound-free  transitions ======> Output # of Transitions
        _kin.check();

        // Check if directory is present for following files
        _kin.setOutputDir(pars.output_dir.c_str());
        FILE* fRunInfo;
        char runInfoDir[100];
        strcpy(runInfoDir,_kin.outputDir);
        strcat(runInfoDir,"SimulationStats.txt");
        fRunInfo = fopen(runInfoDir,"w");
            fprintf(fRunInfo,"Run Type: %s\n",pars.grouping_type.c_str());

        // Count all states
        int Nstate = _kin.Nstate;
        int NCRdim = Nstate;
        printf("Total Number of Isolated Levels and Groups: %d \n",Nstate);
            fprintf(fRunInfo,"Total Number of Isolated Levels and Groups: %d \n",Nstate);
        int nev = 0, nee = 0;
        if (pars.intNe_bool) {
            nev = NCRdim++;
            printf("Electron density will be integrated\n");
            fprintf(fRunInfo,"Electron density will be integrated\n");
        }
        if (pars.intTe_bool) {
            nee = NCRdim++;
            printf("Electron temperature will be integrated\n");
            fprintf(fRunInfo,"Electron temperature will be integrated\n");
        }
        else
            fprintf(fRunInfo,"Electron temperature fixed: %e\n",pars.temp_elec);
        // Initial values
        _Nrc = new double [NCRdim];
        for(int n=0;n<NCRdim;n++) _Nrc[n] = 0.;
        // Density assignment
        double totNe = 0., totNh = 0.;
        for(int n=0;n<_kin.Natoms;n++) {
            double Za = _kin.atom[n]->Z;
            int numLvs = _kin.Ngroups[n];

            for (int m=0;m<numLvs;m++) {
                int ka   = _kin.state[n][m];
                _Nrc[ka] = _kin.Ncr_init[ka];
                if (pars.intNe_bool) _Nrc[nev] += Za*_Nrc[ka];
                if (pars.intTe_bool) _Nrc[nee] += 1.5*Za*_Nrc[ka]*pars.temp_elec;
                totNe += Za*_Nrc[ka];
                totNh += _Nrc[ka];
            }
        }

        // Initialize spectral output
    //	Spectrum spmat = Spectrum(pars.wave_min,pars.wave_max,pars.wave_resol);
    //	spmat.setSpecRange();

    //	if (pars.grouping_type != "QSS") Ne = totNe;	// What is this for?
    //		for(int n=0;n<NCRdim;n++) printf("Initialize %i: %e\n\n", n, _Nrc[n]); exit(1);
        printf("Total atoms: %e\nTotal electrons: %e\n",totNh,totNe);
            fprintf(fRunInfo,"Total atoms: %e\nTotal electrons: %e\n",totNh,totNe);
        printf("-------Ionization Fraction: %e\n\n",totNe/totNh);

        if (_kin.GrpScheme == "SS")
            printf("\n#### LU Decomposition Solve ####\n\n");
        else
            printf("\n#### Backward Euler Integration ####\n\n");

        // Initialize kinetics solver and associated structures
        _kso = CR_BKE::CR_BKE(Nstate,pars.intNe_bool,pars.fxdNe_bool,pars.intTh_bool,pars.intTe_bool);
        _kso->setRad(pars.rad_bool);
        _kso->setIrrad(pars.irrad_bool);
        _kso->addRates(&_kin);
        _kso->setTolerances(pars.rel_tol,pars.abs_tol);
        _kso->multxdt	 = pars.multxdt;
        _kso->setPrintVars(pars.print_log,pars.printdt,pars.logdtprt);
        _kso->SpecOut = pars.print_spectra;
    //	_kso.addSpecFeat(&spmat);
    }

  /**
   * Calculate value of Lorentz force source terms at given time/space
   * locations for given input conserved variables.
   *
   * @param n Length of q array
   * @param tx[0] is time and tx[1..3] are spatial location
   * @param q Conserved variables at tx
   * @param qaux Conserved variables at tx
   * @param s Source term
   */
    bool src(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s)
    {
        pars.t_begin = tx[0];
        pars.t_end = tx[0]+tx[4];

        // assumes that q is [rho_e, e_e, rho_i1, e_i1, rho_i2, e_i2]

        REAL rhoe  = q[0]; pars.dens_elec = rhoe;
        REAL ne = rhoe/_me;
        REAL ue = q[1]/rhoe;
        REAL ve = q[2]/rhoe;
        REAL we = q[3]/rhoe;
        REAL Ee = q[4];
        REAL Pe = (_gas_gamma-1)*(Ee-0.5*rhoe*(ue*ue+ve*ve+we*we));
        REAL Te = Pe/(MKS_kb*ne); pars.temp_elec = Te;

        REAL rho1  = q[5];
        REAL n1 = rho1/_mi1;
        REAL u1 = q[6]/rho1;
        REAL v1 = q[7]/rho1;
        REAL w1 = q[8]/rho1;
        REAL E1 = q[9];
        REAL P1 = (_gas_gamma-1)*(E1-0.5*rho1*(u1*u1+v1*v1+w1*w1)); pars.pressure 	= P1;
        REAL T1 = P1/(MKS_kb*n1); pars.temp_heavy = T1;

        REAL rho2  = q[10];
        REAL n2 = rho1/_mi2;
        REAL u2 = q[11]/rho2;
        REAL v2 = q[12]/rho2;
        REAL w2 = q[13]/rho2;
        REAL E2 = q[14];
        REAL P2 = (_gas_gamma-1)*(E2-0.5*rho2*(u2*u2+v2*v2+w2*w2));
        REAL T2 = P2/(MKS_kb*n2); pars.temp_heavy = T2;

        // Initialization conditions
        double Na = P1/(MKS_kb*T1);     // [m-3] - atomic number density
        double ae = pars.dens_elec/Na;  // approximate ionization fraction (Updated below)
        _kin.unifyAtoms();
        _kin.AtomDist(ae,Na,T1);        // =============> Output # of Levels

        _kso->setTemp(pars.temp_elec,pars.temp_rad);
        _kso->setsol(_Nrc,pars.dens_elec);
        _kso->inst_Tgroup();
        _kso->initRecord();
        _kso->setTimes(pars.t_begin,pars.t_end,pars.t_dt);


        // ------ Begin Simulation ------
        double *src;
        _kso->newIntegrate(src);
        // ------ End Simulation ------

        s[0] = (src[0]-n1)/tx[4];
        s[1] = (src[1]-n2)/tx[4];
//      s[2] = drho_it;
//      s[3] = de_it;
//      s[4] = drho_it;
//      s[5] = de_it;

      return true;
    }

  private:
    REAL _q, _me, _mi1, _mi2, _qbym, _gas_gamma; // charge, mass and charge/mass ratio

    unsigned _Natoms, _temp_bins, dens_min, dens_max;
    REAL temp_min, temp_max, dens_bins, t_dt;
    bool rad, ird, intNe, fxdNe, intTh, intTe;
    char _dirName;
    myiolib::jFileSys *jSys;
    CRParams pars;
    KRatesGrp _kin; // Kinetic processes
    CR_BKE *_kso; // Kinetic solver
    double *_Nrc;
};

#endif // APCOLLISIONALRADIATIVEMODELINGSRC_H
