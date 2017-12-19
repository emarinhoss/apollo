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

        // CR folder to be used
        std::string fname = wxc.template get<std::string>("CR_location");
        _dirName = fname[0];

        jSys = new myiolib::jFileSys(fname);
        pars.init(jSys);
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
        kin =  KRatesGrp(pars.atom_sym.c_str(),pars.Natoms);
        for (int ida=0;ida<pars.Natoms;ida++)
            kin.addAtom(Atoms[ida],pars.atom_cutoffs[ida]);
        delete [] Atoms;
        kin.DistPrep(pars.grouping_type);
        kin.setRange(pars.temp_bins,pars.temp_min,pars.temp_max);
        kin.setNRange(pars.dens_bins,pars.dens_min,pars.dens_max);
        kin.setTPlanck(pars.temp_rad);
        // Add all species to chemistry set and atomic grouping/distributions
        kin.add_BBtrans(pars.irrad_bool);  // Add all bound-bound transitions ======> Output # of Transitions
        kin.add_BFtrans(pars.irrad_bool);  // Add all bound-free  transitions ======> Output # of Transitions
        kin.check();

        // Check if directory is present for following files
        kin.setOutputDir(pars.output_dir.c_str());
        FILE* fRunInfo;
        char runInfoDir[100];
        strcpy(runInfoDir,kin.outputDir);
        strcat(runInfoDir,"SimulationStats.txt");
        fRunInfo = fopen(runInfoDir,"w");
            fprintf(fRunInfo,"Run Type: %s\n",pars.grouping_type.c_str());

        // Count all states
        int Nstate = kin.Nstate;
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
        REAL* Ncr = new REAL [NCRdim];
        for(int n=0;n<NCRdim;n++) Ncr[n] = 0.;
        // Density assignment
        REAL totNe = 0., totNh = 0.;
        for(int n=0;n<kin.Natoms;n++) {
            REAL Za = kin.atom[n]->Z;
            int numLvs = kin.Ngroups[n];

            for (int m=0;m<numLvs;m++) {
                int ka   = kin.state[n][m];
                Ncr[ka] = kin.Ncr_init[ka];
                if (pars.intNe_bool) Ncr[nev] += Za*Ncr[ka];
                if (pars.intTe_bool) Ncr[nee] += 1.5*Za*Ncr[ka]*pars.temp_elec;
                totNe += Za*Ncr[ka];
                totNh += Ncr[ka];
            }
        }

        // Initialize spectral output
    //	Spectrum spmat = Spectrum(pars.wave_min,pars.wave_max,pars.wave_resol);
    //	spmat.setSpecRange();

    //	if (pars.grouping_type != "QSS") Ne = totNe;	// What is this for?
    //		for(int n=0;n<NCRdim;n++) printf("Initialize %i: %e\n\n", n, Ncr[n]); exit(1);
        printf("Total atoms: %e\nTotal electrons: %e\n",totNh,totNe);
            fprintf(fRunInfo,"Total atoms: %e\nTotal electrons: %e\n",totNh,totNe);
        printf("-------Ionization Fraction: %e\n\n",totNe/totNh);

        if (kin.GrpScheme == "SS")
            printf("\n#### LU Decomposition Solve ####\n\n");
        else
            printf("\n#### Backward Euler Integration ####\n\n");

        // Initialize kinetics solver and associated structures
        kso = CR_BKE(Nstate,pars.intNe_bool,pars.fxdNe_bool,pars.intTh_bool,pars.intTe_bool);
        kso.setRad(pars.rad_bool);
        kso.setIrrad(pars.irrad_bool);
        kso.addRates(&kin);
        kso.setTemp(pars.temp_elec,pars.temp_rad);
        kso.setsol(Ncr,pars.dens_elec); delete [] Ncr;
        kso.inst_Tgroup();
        kso.initRecord();
        kso.setTimes(pars.t_begin,pars.t_end,pars.t_dt);
        kso.setTolerances(pars.rel_tol,pars.abs_tol);
        kso.multxdt	 = pars.multxdt;
        kso.setPrintVars(pars.print_log,pars.printdt,pars.logdtprt);
        kso.SpecOut = pars.print_spectra;
    //	kso.addSpecFeat(&spmat);
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
        REAL P1 = (_gas_gamma-1)*(E1-0.5*rho1*(u1*u1+v1*v1+w1*w1));
        REAL T1 = P1/(MKS_kb*n1); pars.temp_heavy = T1;




        // Initialization conditions
        REAL ThV = pars.temp_heavy;						// [eV]  - Heavy atom temperature
        REAL Pc = pars.pressure;							// Initial pressure of system
        REAL Na = Pc*MKS_atm/(MKS_kb*760*MKS_eV_K*T1);	// [m-3] - atomic number density
        REAL Ne = pars.dens_elec;							// [m-3] - approximate electron density
        REAL ae = Ne/Na;  								// approximate ionization fraction (Updated below)
        kin.unifyAtoms();
        kin.AtomDist(ae,Na,T1);		// =============> Output # of Levels


//      s[0] = drho_et;
//      s[1] = de_et;
//      s[2] = drho_it;
//      s[3] = de_it;
//      s[4] = drho_it;
//      s[5] = de_it;

      return true;
    }

  private:
    REAL _q, _me, _mi1, _mi_2, _qbym, _gas_gamma; // charge, mass and charge/mass ratio
    unsigned _Natoms, _temp_bins, dens_min, dens_max;
    REAL temp_min, temp_max, dens_bins, t_dt;
    bool rad, ird, intNe, fxdNe, intTh, intTe;
    char _dirName;
    myiolib::jFileSys *jSys;
    CRParams pars;
    KRatesGrp kin; // Kinetic processes
    CR_BKE kso; // Kinetic solver
};

#endif // APCOLLISIONALRADIATIVEMODELINGSRC_H
