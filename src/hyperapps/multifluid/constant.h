/* 
 * File:   constant.h
 * Author: cambier
 *
 * Created on July 21, 2011, 1:58 PM
 */

#ifndef CONSTANT_H
#define	CONSTANT_H
//
//#define M_PI 3.14159265

//typedef float real;

// device information
#define PAGE_LOCKED 0          // page-locked memory on host?
#define USE_CUBLAS  1          // use CUBLAS?
#define USE_ATOMIC  0          // use atomic operations?
#define BLOCK_SIZE 64          // default block size

// math information
#define FOR_ARRAY  1           // if true, arrays stired Fortran-style (1st-index varying first), else C-style - see hmatrix.h
#define MXN_RSEED  4096

// particle information
#define MXN_PDOM   200000      // maximum number of particles in domain
#define MXN_PCELL  100         // maximum number of particles in a cell

// species information
#define MXN_PSPEC      1024    // 2^10 = maximum number of particle species
#define MXN_LEVEL   4194304    // 2^22 = maximum number of levels

// domain information
#define MXN_CELLS  16777216    // 2^24
#define MXN_DOMAIN      256    // 2^8
#define DIM_FIELDS        3    // number field components
#define MXN_IJK        1024    // maximum index in any direction (i,j or k)

#define iunder  -999999
#define runder  -999999.9

// physical constants
#define MKS_a0      5.291772109e-11	// Bohr radius
#define MKS_Ry      1.097373157e+07
#define CGS_Ry		13.60569253
#define MKS_IH      13.6            // Hydrogen ionization potential (eV)

// electro-magnetism
#define MKS_c       2.997924580e+08
#define MKS_eps0    8.854187818e-12
#define MKS_mu0     1.256637061e-06
#define MKS_e       1.602189200e-19  // why different?
#define MKS_me      9.109534000e-31
#define MKS_mp      1.672510442e-27
#define MKS_amu     1.660565472e-27
#define MKS_mec     2.730924201e-22     // me*c   relativistic momentum scale
#define MKS_mec2    8.187104787e-14     // me*c2  relativistic energy scale
#define MKS_fstr	0.0072973525664		// Fine structure constant

// thermodynamics
#define MKS_Rb      8.314408694         // Boltzmann constant (mole)
#define MKS_Avo     6.022045e+23        // Avogadro number
#define MKS_kb      1.380662e-23        // Boltzmann constant
#define MKS_hPlanck 6.626176e-34        // Planck constant
#define MKS_hbar    1.054588664e-34     // Planck Cnst/(2*pi)
#define MKS_mole    6.023e+23           // mole (statistical weight unit)

// scaling
#define yotta       1.e+24
#define zetta       1.e+21
#define exa         1.e+18
#define peta        1.e+15
#define tera        1.e+12
#define giga        1.e+09
#define mega        1.e+06
#define kilo        1.e+03
#define milli       1.e-03
#define micro       1.e-06
#define nano        1.e-09
#define pico        1.e-12
#define femto       1.e-15
#define atto        1.e-18
#define zepto       1.e-21
#define yocto       1.e-24

// unit conversion
#define MKS_atm     101325.             // 1 bar in Pascal
#define MKS_psi     6894.75729          // 1 psi in Pascal
#define MKS_torr    133.322             // 1 Torr in Pascal
#define MKS_ETeV	1.6021766e-19		// 1 eV energy (Joules)
#define MKS_TeV     11604.5             // 1 eV temperature (Kelvin)
#define MKS_TRank   0.555555            // Rankine
#define MKS_Rcm     0.695038122         // convert 1/cm into Joule
#define MKS_KcalJ   4186.05             // kilo-cal into Joules
#define MKS_gram    1.e-03
#define MKS_erg     1.e-07

#define MKS_cm      1.e-02              // 1 cm = 0.01 m
#define MKS_mm      1.e-03              // 1 mm = 0.001 m
#define MKS_cm2     1.e-04
#define MKS_cm3     1.e-06
#define MKS_foot    0.3048
#define MKS_inch    2.54e-2
#define MKS_micron  1.e-06
#define MKS_km      1.e+03
#define MKS_mile    1562.

#define MKS_mn      60.                 // 1 minute = 60 sec
#define MKS_hr      3600.               // 1 hour = 3600 sec
#define MKS_ms      1.e-03              // 1 milli-sec
#define MKS_mus     1.e-06              // 1 micro-sec
#define MKS_ns      1.e-09              // 1 nao-sec
#define MKS_ps      1.e-12              // 1 pico-sec
#define MKS_fs      1.e-15              // 1 femto-sec

#define MKS_Gauss   1.e-04              // 1 Gauss in Tesla

#define MKS_Hart	4.35974417			// 1 Hartree in Joules

#define AU_Hart		27.21138602			// 1 Hartree in eV

#define iund  -999999
#define rund -999999.9

#endif	/* CONSTANT_H */

