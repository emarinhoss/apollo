/*
 * Constants.h
 *
 *  Created on: Feb 25, 2013
 *      Author: hle
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define M_PI 3.14159265358979323846
#define M_2PI 6.28318530718

// atomic and nuclear physics
#define MKS_a0      5.291772109e-11
#define MKS_Ry      1.097373157e+07
#define MKS_IH      13.6
#define RYD_eV		13.6056925330

// electro-magnetism
#define MKS_c       2.997924580e+08
#define MKS_eps0    8.854187818e-12
#define MKS_mu0     1.256637061e-06
#define MKS_e       1.602189200e-19
#define MKS_me      9.109534000e-31
#define MKS_mp      1.672510442e-27
#define MKS_amu     1.660565472e-27
#define MKS_mec     2.730924201e-22    	// me*c   relativistic momentum scale
#define MKS_mec2    8.187104787e-14    	// me*c2  relativistic energy scale

// thermodynamics
#define MKS_Rb      8.314408694        	// Boltzmann constant (mole)
#define MKS_Avo     6.022045e+23       	// Avogadro number
#define MKS_kb      1.380662e-23       	// Boltzmann constant
#define MKS_hPlanck 6.626176e-34       	// Planck constant
#define MKS_hbar    1.054588664e-34    	// Planck Cnst/(2*pi)
#define MKS_hc		1.98644568e-25     	// Planck constant * speed of light [J * m]
#define MKS_eV_hc	1.23984193e-6      	// Planck constant * speed of light [eV * m]
#define MKS_mole    6.023e+23          	// mole (statistical weight unit)

// unit conversion
#define MKS_atm     101325.            	// 1 bar in Pascal
#define MKS_eV_K    11605.5            	// 1 eV temperature (Kelvin)
#define MKS_eV_nm	1240			   	//
#define MKS_Rcm     0.695038122        	// convert 1/cm into Joule
#define MKS_KcalJ   4186.05            	// kilo-cal into Joules
#define MKS_cm      1.e-02             	// 1 cm = 0.01 m
#define MKS_cm2     1.e-04
#define MKS_cm3     1.e-06
#define MKS_us      1.e-06
#define MKS_ps      1.e-12
#define MKS_TeraHz  1.e+12


#endif /* CONSTANTS_H_ */
