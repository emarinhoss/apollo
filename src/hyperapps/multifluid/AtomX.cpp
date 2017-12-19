#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

// interpolation
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "myIO.h"
#include "constant.h"
#include "Constants.h"
#include "AtomX.h"
#include "cArray.h"
#include "MathX.h"

using namespace std;

// Default Constructor
AtomX::AtomX():ge(NULL),Ee(NULL),Ei(NULL),ion(NULL),Z(0),Znuc(0),amu(0.),mX(0.),name("Void")
{
	ge = NULL;
	nS = NULL;
	qn = NULL;
	ql = NULL;
	qw = NULL;
	Ee = NULL;
	Ei = NULL;
    Asum = NULL;
    Nalevs=Nlevels= 0;
    Neirxs=Neecxs=Ndirxs=Nexdxs=Npirxs= 0;
	NBBeXD=NBFeIR=NBFpIR=NBFaIR=NBFdIR=NBBlns=0;
	bbrad=bbexd=bfeir=bfpir=bfair=bfdir=NULL;
	bbeexxs=bfeioxs=bfpioxs=bfdioxs=NULL;
	bbradfg=bbeexxs_intp=bfeioxs_intp=bfpioxs_intp=bheecxs=NULL;
	bfaiort=NULL;
}

// Constructor [Active]
AtomX::AtomX(string s):ion(NULL),Z(0),Znuc(0),amu(0.),Ee(NULL),Ei(NULL)
{
	name = s;
	ge = NULL;
	nS = NULL;
	qn = NULL;
	ql = NULL;
	qw = NULL;
	Ee = NULL;
	Ei = NULL;
    Asum = NULL;
	NBBeXD=NBFeIR=NBFpIR=NBFaIR=NBBlns=0;
	bbrad=bbexd=bfeir=bfpir=bfair=NULL;
	bbeexxs=bfeioxs=bfpioxs=NULL;
	bbradfg=bbeexxs_intp=bfeioxs_intp=bfpioxs_intp=bheecxs=NULL;
	bfaiort=NULL;
    set_name(s);
}

// Destructor
AtomX::~AtomX()
{
    del_cArray<int>  (Nlevels,1,qn);
    del_cArray<int>  (Nlevels,1,ql);
    del_cArray<int>  (Nlevels,1,qw);
    del_cArray<int>  (Nlevels,ge);
    del_cArray<int>  (Nlevels,nS);
    del_cArray<float>(Nlevels,Ee);
    del_cArray<float>(Nlevels,Ei);
    del_cArray<float>(Nlevels,Asum);

    del_cArray<double>(NBBeXD,Nexdxs,2,bbeexxs);
    del_cArray<double>(NBFeIR,Neirxs,2,bfeioxs);
    del_cArray<double>(NBFpIR,Npirxs,2,bfpioxs);
    del_cArray<double>(NBBlns,4,bbradfg);
    del_cArray<double>(NBBeXD,Nexdxs,bbeexxs_intp);
	del_cArray<double>(NBFeIR,Neirxs,bfeioxs_intp);
	del_cArray<double>(NBFpIR,Npirxs,bfpioxs_intp);
	del_cArray<double>(NBFaIR,bfaiort);
	del_cArray<double>(Neecxs,2,bheecxs);

    del_cArray<transition>(NBBlns,bbrad);
    del_cArray<transition>(NBBeXD,bbexd);
    del_cArray<transition>(NBFeIR,bfeir);
    del_cArray<transition>(NBFpIR,bfpir);
    del_cArray<transition>(NBFaIR,bfair);
}

// Set the names of ground state atom and its ion stages
void AtomX::set_name(string s)
{
	// Accounts for higher stages so long
	// as the convention '__+(int)' is followed
    int ls = s.length();
    char* sc = (char *)s.data();
//    int ls = sizeof(s);

    char num[2] = {0};
    int nind = 0;
    for (int i = 0; i < ls; i++) {
    	if(sc[i]=='-') Z--;
    	if(sc[i]=='+') Z++;
		if(isdigit(sc[i])) {
			num[nind] = sc[i];
			nind++;
		}
    }
//    for (int i = 0; i < ls; i++) {
//    	if(s[i]=='-') Z--;
//    	if(s[i]=='+') Z++;
//		if(isdigit(s[i])) {
//			num[nind] = s[i];
//			nind++;
//		}
//    }
    int mult = 1;
    if (num[0] != 0) {
		sscanf(num,"%i",&mult);
    	Z = Z*mult;
    }

    printf("Atom construction: %s\n",name.c_str());
}

// Append ion stage to specified atom
void AtomX::set_ion(AtomX* inext)
{
    ion = inext;
}

void AtomX::hydroIon()
{
    Nlevels = 1;				// 1 energy state level for ionized Hydrogen
    ge = new int[Nlevels];
	Ee = new float[Nlevels];
	Ei = new float[Nlevels];
    Ee[0] = 0.;
    Ei[0] = 0.;
    ge[0] = 1 ;
}

// Store atomic and cross-sectional database into memory
void AtomX::AtomParams(string atomstr)
{
	// Provide file names for access
	printf("Reading data of Atom %s\n",name.c_str());
	ostringstream itostr;
	itostr << setw(2) << setfill('0') << Z+1;	// Files are numbered spectroscopically
	string filenum = itostr.str();

	string atname	= atomstr + "/" + atomstr;
	string fsname	= "./Atomic_Data/"+ atname + filenum + "fstruct_r.dat";
	string eiename 	= "./Atomic_Data/"+ atname + filenum + "eielev.dat";
	string eiiname 	= "./Atomic_Data/"+ atname + filenum + "eionlev.dat";
	string piname 	= "./Atomic_Data/"+ atname + filenum + "pionlev.dat";
	string osname 	= "./Atomic_Data/"+ atname + filenum + "osc_r.dat";
	string ainame 	= "./Atomic_Data/"+ atname + filenum + "ai.dat";
	string dwname 	= "./Atomic_Data/"+ atname + filenum + "dwilev.dat";
	string ecname	= "";
	if (Z == 0) ecname	= "./Atomic_Data/"+ atname + filenum + "eh_cc.dat";

	/* number of xs data points */
	Nexdxs = 9;
	Neirxs = 21;
	Npirxs = 21;
	Ndirxs = 21;

	/* ---------------------------------------------------- */
    /* -------------- Fine-Structure Section -------------- */
	/* ---------------------------------------------------- */

    int numlines;
    ifstream fi;
    std::string line;
    vector<myiolib::jString> tokens;
    int idx, lvl, lvl1, lvl2;
    float jfs;

    fi.open(fsname.c_str(),ios_base::in);
	if(fi.fail())
	{
		cout << "File:" << fsname <<": could not be opened" << endl;
		abort();
	}
	for (numlines = 0; std::getline(fi, line); ++numlines);

    Nlevels = numlines;
    ge = new int[Nlevels];
    nS = new int[Nlevels];
    qn = new int*[Nlevels];
    ql = new int*[Nlevels];
    qw = new int*[Nlevels];
    Ee = new float[Nlevels];
    Ei = new float[Nlevels];
    Asum = cArray<float>(Nlevels);
    char llabel[] = "spdfgh ";

    idx = 0;
	lvl = 0;
	fi.clear(); fi.seekg(0);
	while(!fi.eof())
	{
		if(getline(fi,line) )
		{
			myiolib::jString jline(line);
			if(jline.length() > 0) {
				tokens = jline.parse();
				if (idx != tokens[0].iValue()-1) { printf("Missing level %d\n", idx+1); abort(); }
				jfs = tokens[1].dValue();

				ge[idx] = 2*jfs + 1;			// Degeneracy
				nS[idx] = tokens[2].iValue();	// Shell
				Ee[idx] = tokens[3].dValue();	// Level Energy

				qn[idx] = new int[nS[idx]];
				ql[idx] = new int[nS[idx]];
				qw[idx] = new int[nS[idx]];

				if (tokens.size()-4 != 7*nS[idx] )
				{
					printf("No. of quantum numbers and shell don't match\n");
					abort();
				}

				for (int s=0;s<nS[idx];s++)
				{
					qn[idx][s] = tokens[s*7+4].iValue();
					ql[idx][s] = tokens[s*7+5].iValue();
					qw[idx][s] = tokens[s*7+6].iValue();
				}

				printf("%d: ",idx+1);
				for (int s=0;s<nS[idx];s++)
					printf("%d%c%d ",qn[idx][s],llabel[ql[idx][s]],qw[idx][s]);
				printf("\n");

				tokens.clear();
				idx++;
			}

		}
	}
	fi.close();
	tokens.clear();

	if (idx == Nlevels) {
		printf("Completed Fine-Structure Data\n");
		printf("----Number of atomic levels: %i\n",Nlevels);
	}
	else {
		printf("No. of levels and no. of txtlines don't match... check for empty txtlines\n");
		abort();
	}

	// renormalize energy level st ground state is 0
	double Eref = Ee[0];
	for(int i=0;i<Nlevels;i++) Ee[i] = Ee[i]-Eref;
	printf("----Normalized maximum energy is %3.3e\n",Ee[Nlevels-1]);

	if (Nlevels == 1) return void();

	/* ---------------------------------------------------------------- */
	/* -------------- Electron-Impact Excitation Section -------------- */
	/* ---------------------------------------------------------------- */

	fi.open(eiename.c_str(),ios_base::in);
	if(fi.fail())
	{
		cout << "File:" << eiename <<": could not be opened" << endl;
		abort();
	}
	for (numlines = 0; std::getline(fi, line); ++numlines);

	bbexd = new transition[numlines];
	bbeexxs = cArray<double>(numlines,Nexdxs,2);
	bbeexxs_intp = cArray<double>(numlines,Nexdxs);
	double** xs_dat = cArray<double>(2,30);
	// create interpolator
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();

	idx = 0;lvl1=lvl2=0;
	fi.clear(); fi.seekg(0);
	while(!fi.eof())
	{
		if(getline(fi,line) )
		{
			myiolib::jString jline(line);
			//jline = jline.clean();
			if(jline.length() > 0) {
				tokens = jline.clean().parse();

				if (tokens[1].iValue() == lvl1 && tokens[3].iValue() == lvl2)
					{ printf("repeated transition %d - %d\n",lvl1,lvl2); abort(); }
				lvl1 = tokens[1].iValue();
				lvl2 = tokens[3].iValue();

				int low = lvl1-1;int upp = lvl2-1;
				bbexd[idx].low[0] = 0;
				bbexd[idx].upp[0] = 0;
				bbexd[idx].low[1] = low;
				bbexd[idx].upp[1] = upp;
				bbexd[idx].low[2] = tokens[0].iValue();
				bbexd[idx].upp[2] = tokens[2].iValue();
				bbexd[idx].dE     = Ee[upp]-Ee[low];

				bool is_negative = false;

				for(int i=0; i<Nexdxs; i++) {
					bbeexxs[idx][i][0] = tokens[i*2+4].dValue();
					bbeexxs[idx][i][1] = tokens[i*2+5].dValue()*1.e-4; //cm2 to m2
					if (bbeexxs[idx][i][0]<0) is_negative = true;
					if (bbeexxs[idx][i][1]<0) is_negative = true;
					if (bbeexxs[idx][i][0]<0) bbeexxs[idx][i][0] = abs(bbeexxs[idx][i][0]);
					if (bbeexxs[idx][i][1]<0) bbeexxs[idx][i][1] = abs(bbeexxs[idx][i][1]);
					xs_dat[0][i] = log(bbeexxs[idx][i][0]);
					xs_dat[1][i] = log(bbeexxs[idx][i][1]);
		//			printf("xs_dat[0][%d] = %e\t",i,xs_dat[0][i]);
		//			printf("xs_dat[1][%d] = %e\n",i,xs_dat[1][i]);
				}

				if (is_negative)
				{
		//			printf("!!!!Negative excitation cross section values: Atom %s: %d -> %d\n",name.c_str(),low,upp);
		//			for(int i=0; i<Nexdxs; i++) {
		//				printf("E = %e\t",bbeexxs[idx][i][0]);
		//				printf("xs = %e\n",bbeexxs[idx][i][1]);
		//			}
		//			abort();
				}

				bbexd[idx].NEpts = Nexdxs;
				bbexd[idx].lnEmin = log(bbeexxs[idx][0][0]);
				bbexd[idx].lnEmax = log(.9999*bbeexxs[idx][Nexdxs-1][0]);
				bbexd[idx].dlnE = (bbexd[idx].lnEmax-bbexd[idx].lnEmin)/(float)(Nexdxs-1);
				bbexd[idx].xs_dat = bbeexxs_intp[idx];

		//		printf("idx = %d\n",idx);
		//		printf("NEpts = %d\n",bbexd[idx].NEpts);
		//		printf("lnEmin = %e\n",exp(bbexd[idx].lnEmin));
		//		printf("lnEmax = %e\n",exp(bbexd[idx].lnEmax));

				// create interpolator
                gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, Nexdxs);
                gsl_spline_init (spline, xs_dat[0], xs_dat[1], Nexdxs);

				for(int i=0; i<Nexdxs; i++) {
					double _lnE = bbexd[idx].lnEmin + (double)i*bbexd[idx].dlnE;
					_lnE = max(_lnE,xs_dat[0][0]);
                    bbeexxs_intp[idx][i] = exp((gsl_spline_eval (spline, _lnE, acc)));
				}
                gsl_spline_free (spline);

		//		printf("Here %i %i %e %e %e %e\n",bbexdlv[idx][0],bbexdlv[idx][1],bbeexxs[idx][0][0],bbeexxs[idx][0][1],bbeexxs[idx][8][0],bbeexxs[idx][8][1]);
				tokens.clear();
				idx++;
			}
		}
	}
	fi.close();
	tokens.clear();

	if (idx == numlines) {
		NBBeXD = idx;
		printf("Completed Electron-Impact Excitation Data\n");
		printf("----Number of electron-impact x/d: %i\n",NBBeXD);
	}
	else {
		printf("no of exc transitions and no of txtlines don't match... check for empty txtlines\n");
		abort();
	}

	/* ---------------------------------------------------------------- */
	/* -------------- Electron-Impact Ionization Section -------------- */
	/* ---------------------------------------------------------------- */
	fi.open(eiiname.c_str(),ios_base::in);
	if(fi.fail())
	{
		cout << "file:" << eiiname <<": could not be opened" << endl;
		abort();
	}
	for (numlines = 0; std::getline(fi, line); ++numlines);

	bfeir = new transition[numlines];
	bfeioxs      = cArray<double>(numlines,Neirxs,2);
	bfeioxs_intp = cArray<double>(numlines,Neirxs);

	idx = 0;lvl1=lvl2=0;
	fi.clear(); fi.seekg(0);
	while(!fi.eof())
	{
		if(getline(fi,line) )
		{
			myiolib::jString jline(line);
			//jline = jline.clean();
			if(jline.length() > 0) {
				tokens = jline.clean().parse();

				if (tokens[1].iValue() == lvl1 && tokens[3].iValue() == lvl2)
					{ printf("repeated transition %d - %d\n",lvl1,lvl2); abort(); }

				lvl1 = tokens[1].iValue();
				lvl2 = tokens[3].iValue();

				int low = lvl1-1;int upp = lvl2-1;
				bfeir[idx].low[0] = 0;
				bfeir[idx].upp[0] = 0;
				bfeir[idx].low[1] = low;
				bfeir[idx].upp[1] = upp;
				bfeir[idx].low[2] = tokens[0].iValue();
				bfeir[idx].upp[2] = tokens[2].iValue();
				bfeir[idx].dE     = tokens[4].dValue();

				bool is_negative = false;

				for(int i=0; i<Neirxs; i++) {
					bfeioxs[idx][i][0] = tokens[i*2+5].dValue();
					bfeioxs[idx][i][1] = tokens[i*2+6].dValue()*1.e-4; //cm2 to m2
					if (bfeioxs[idx][i][0]<0) is_negative = true;
					if (bfeioxs[idx][i][1]<0) is_negative = true;
					if (bfeioxs[idx][i][0]<0) bfeioxs[idx][i][0] = abs(bfeioxs[idx][i][0]);
					if (bfeioxs[idx][i][1]<0) bfeioxs[idx][i][1] = abs(bfeioxs[idx][i][1]);
					xs_dat[0][i] = log(bfeioxs[idx][i][0]);
					xs_dat[1][i] = log(bfeioxs[idx][i][1]);
		//			printf("xs_dat[0][%d] = %e\t",i,xs_dat[0][i]);
		//			printf("xs_dat[1][%d] = %e\n",i,xs_dat[1][i]);
				}

				if (is_negative)
				{
					printf("!!!!Negative ionization cross section values: Atom %s: %d -> %d\n",name.c_str(),bfeir[idx].low[1],bfeir[idx].upp[1]);
		//			for(int i=0; i<Neirxs; i++) {
		//				printf("E = %e\t",bfeioxs[idx][i][0]);
		//				printf("xs = %e\n",bfeioxs[idx][i][1]);
		//			}
		//			abort();
				}

				bfeir[idx].NEpts = Neirxs;
				bfeir[idx].lnEmin = log(bfeioxs[idx][0][0]);
				bfeir[idx].lnEmax = log(.9999*bfeioxs[idx][Neirxs-1][0]);
				bfeir[idx].dlnE = (bfeir[idx].lnEmax-bfeir[idx].lnEmin)/(float)(Neirxs-1);
				bfeir[idx].xs_dat = bfeioxs_intp[idx];

				// create interpolator
				gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, Neirxs);
				gsl_spline_init (spline, xs_dat[0], xs_dat[1], Neirxs);

				for(int i=0; i<Neirxs; i++) {
					double _lnE = bfeir[idx].lnEmin + (double)i*bfeir[idx].dlnE;
					_lnE = max(_lnE,xs_dat[0][0]);
					bfeioxs_intp[idx][i] = exp((gsl_spline_eval (spline, _lnE, acc)));
					if (bfeioxs_intp[idx][i]<0) { printf("negative xs after intepolation\n"); abort(); }
				}
				gsl_spline_free (spline);

		//		printf("Here %i %i %e %e %e %e\n",bfeirlv[idx][0],bfeirlv[idx][1],bfeioxs[idx][0][0],bfeioxs[idx][0][1],bfeioxs[idx][20][0],bfeioxs[idx][20][1]);
		//		printf("%i\n",idx);
				tokens.clear();
				idx++;
			}
		}
	}
	fi.close();
	tokens.clear();

	// Check if numlines is equal to idx
	if (idx == numlines) {
		NBFeIR = idx;
		printf("Completed Electron-Impact Ionization Data\n");
		printf("----Number of electron-impact i/r: %i\n",NBFeIR);
	}
	else {
		printf("no of eion transitions and no of txtlines don't match... check for empty txtlines\n");
		abort();
	}


	/* ----------------------------------------------------- */
	/* -------------- Photoionization Section -------------- */
	/* ----------------------------------------------------- */

	fi.open(piname.c_str(),ios_base::in);
	if(fi.fail())
	{
		cout << "file:" << piname <<": could not be opened" << endl;
		abort();
	}
	for (numlines = 0; std::getline(fi, line); ++numlines);

	bfpir = new transition[numlines];
	bfpioxs = cArray<double>(numlines,Npirxs,2);
	bfpioxs_intp = cArray<double>(numlines,Npirxs);

	idx = 0;lvl1=lvl2=0;
	fi.clear(); fi.seekg(0);
	while(!fi.eof())
	{
		if(getline(fi,line) )
		{
			myiolib::jString jline(line);
			//jline = jline.clean();
			if(jline.length() > 0) {
				tokens = jline.clean().parse();

				if (tokens[1].iValue() == lvl1 && tokens[3].iValue() == lvl2)
					{ printf("repeated transition %d - %d\n",lvl1,lvl2); abort(); }

				lvl1 = tokens[1].iValue();
				lvl2 = tokens[3].iValue();
				int low = lvl1-1;int upp = lvl2-1;

				bfpir[idx].low[0] = 0;
				bfpir[idx].upp[0] = 0;
				bfpir[idx].low[1] = low;
				bfpir[idx].upp[1] = upp;
				bfpir[idx].low[2] = tokens[0].iValue();
				bfpir[idx].upp[2] = tokens[2].iValue();
				bfpir[idx].dE     = tokens[4].dValue();

				bool is_negative = false;

				for(int i=0; i<Npirxs; i++) {
					bfpioxs[idx][i][0] = tokens[i*2+5].dValue();
					bfpioxs[idx][i][1] = tokens[i*2+6].dValue()*1.e-4; //cm2 to m2
					if (bfpioxs[idx][i][0]<0) is_negative = true;
					if (bfpioxs[idx][i][1]<0) is_negative = true;
					if (bfpioxs[idx][i][0]<0) bfpioxs[idx][i][0] = abs(bfpioxs[idx][i][0]);
					if (bfpioxs[idx][i][1]<0) bfpioxs[idx][i][1] = abs(bfpioxs[idx][i][1]);
					xs_dat[0][i] = log(bfpioxs[idx][i][0]);
					xs_dat[1][i] = log(bfpioxs[idx][i][1]);
		//			printf("xs_dat[0][%d] = %e\t",i,xs_dat[0][i]);
		//			printf("xs_dat[1][%d] = %e\n",i,xs_dat[1][i]);
				}

				if (is_negative)
				{
					printf("!!!!Negative photo-ionization cross section values: Atom %s: %d -> %d\n",name.c_str(),bfpir[idx].low[1],bfpir[idx].upp[1]);
		//			for(int i=0; i<Npirxs; i++) {
		//				printf("E = %e\t",bfpioxs[idx][i][0]);
		//				printf("xs = %e\n",bfpioxs[idx][i][1]);
		//			}
		//			abort();
				}

				bfpir[idx].NEpts = Npirxs;
				bfpir[idx].lnEmin = log(bfpioxs[idx][0][0]);
				bfpir[idx].lnEmax = log(.9999*bfpioxs[idx][Npirxs-1][0]);
				bfpir[idx].dlnE = (bfpir[idx].lnEmax-bfpir[idx].lnEmin)/(float)(Npirxs-1);
				bfpir[idx].xs_dat = bfpioxs_intp[idx];

				// create interpolator
				gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, Npirxs);
				gsl_spline_init (spline, xs_dat[0], xs_dat[1], Npirxs);

				for(int i=0; i<Npirxs; i++) {
					double _lnE = bfpir[idx].lnEmin + (double)i*bfpir[idx].dlnE;
					_lnE = max(_lnE,xs_dat[0][0]);
					bfpioxs_intp[idx][i] = exp((gsl_spline_eval (spline, _lnE, acc)));
				}
				gsl_spline_free (spline);

		//		printf("Here %i %i %e %e %e %e\n",bfpirlv[idx][0],bfpirlv[idx][1],bfpioxs[idx][0][0],bfpioxs[idx][0][1],bfpioxs[idx][20][0],bfpioxs[idx][20][1]);
				tokens.clear();
				idx++;
			}
		}
	}
	fi.close();
	tokens.clear();

	// Check if numlines is equal to idx
	if (idx == numlines) {
		NBFpIR = idx;
		printf("Completed Photo-Ionization Data\n");
		printf("----Number of photon-impact i/r: %i\n",NBFpIR);
	}
	else {
		printf("no of pion transitions and no of txtlines don't match... check for empty txtlines\n");
		abort();
	}

	/* --------------------------------------------------------- */
	/* -------------- Oscillator Strength Section -------------- */
	/* --------------------------------------------------------- */
	fi.open(osname.c_str(),ios_base::in);
	if(fi.fail())
	{
		cout << "file:" << osname <<": could not be opened" << endl;
		abort();
	}
	for (numlines = 0; std::getline(fi, line); ++numlines);

	bbrad = new transition[numlines];
	bbradfg = cArray<double>(numlines,4);

	idx = 0;lvl1=lvl2=0;
	fi.clear(); fi.seekg(0);
	while(!fi.eof())
	{
		if(getline(fi,line) )
		{
			myiolib::jString jline(line);
			//jline = jline.clean();
			if(jline.length() > 0) {
				tokens = jline.clean().parse();

				if (tokens[3].iValue() == lvl1 && tokens[4].iValue() == lvl2)
					{ printf("repeated transition %d - %d\n",lvl1,lvl2); abort(); }

				lvl1 = tokens[3].iValue();
				lvl2 = tokens[4].iValue();
				int low = lvl1-1;int upp = lvl2-1;

				bbrad[idx].low[0] = 0;
				bbrad[idx].upp[0] = 0;
				bbrad[idx].low[1] = low;
				bbrad[idx].upp[1] = upp;
				bbrad[idx].low[2] = Z;
				bbrad[idx].upp[2] = Z;
				bbrad[idx].dE     = Ee[upp]-Ee[low];

				float grat = (float)ge[low]/(float)ge[upp];
				bbradfg[idx][0] = tokens[0].dValue(); // wavelength in nm
				bbradfg[idx][1] = tokens[1].dValue(); // dE in eV
				bbradfg[idx][2] = tokens[2].dValue(); // oscillator strength
				bbradfg[idx][3] = 6.6702e+13/bbradfg[idx][0]/bbradfg[idx][0]*grat*bbradfg[idx][2]; // Einstein coef.
				tokens.clear();
				idx++;
			}
		}
	}
	fi.close();
	tokens.clear();
	del_cArray<double>(2,30,xs_dat);
	gsl_interp_accel_free (acc);

	// Check if numlines is equal to idx
	if (idx == numlines) {
		NBBlns = idx;
		printf("Completed Oscillator Strength Data\n");
		printf("----Number of atomic lines: %i\n",NBBlns);
	}
	else {
		printf("no of radlines and no of txtlines don't match... check for empty txtlines\n");
		abort();
	}

	// compute inverse radiative lifetimes
	for (int li=0;li<NBBlns;li++)
	{
		int upp    = bbrad[li].upp[1];
		float _Aul = bbradfg[li][3];
		Asum[upp] += _Aul;
	}

	/* -------------- Autoionization/electron capture Section -------------- */
	fi.open(ainame.c_str(),ios_base::in);
	if(fi.fail())
	{
		cout << "file:" << ainame <<": could not be opened" << endl;
		cout << "no autoionization data for " << name << endl;
	}
	for (numlines = 0; std::getline(fi, line); ++numlines);

	bfair   = new transition[numlines];
	bfaiort = cArray<double>(numlines);

	idx = 0;
	fi.clear(); fi.seekg(0);
	while(!fi.eof() && !fi.fail())
	{
		if(getline(fi,line) )
		{
			myiolib::jString jline(line);
			//jline = jline.clean();
			if(jline.length() > 0) {
				tokens = jline.clean().parse();

				int low = tokens[1].iValue()-1;
				int upp = tokens[3].iValue()-1;
//				printf("low = %d upp = %d\n",low,upp);

				bfair[idx].low[0] = 0;
				bfair[idx].upp[0] = 0;
				bfair[idx].low[1] = low;
				bfair[idx].upp[1] = upp;
				bfair[idx].low[2] = tokens[0].iValue();
				bfair[idx].upp[2] = tokens[2].iValue();
				bfair[idx].dE     = tokens[4].dValue();

				bfaiort[idx] = tokens[5].dValue();

				tokens.clear();
				idx++;
			}
		}
	}
	fi.close();
	tokens.clear();

	if (idx == numlines) {
		NBFaIR = idx;
		printf("Completed Auto-Ionization Data\n");
		printf("----Number of ai/dr: %i\n",NBFaIR);
	}
	else {
		printf("no of aion transitions and no of txtlines don't match... check for empty txtlines\n");
		abort();
	}

	/* --------------------------------------------------------- */
	/* -------------- Elastic Collisions Section --------------- */
	/* --------------------------------------------------------- */

	fi.open(ecname.c_str(),ios_base::in);
	if(fi.fail())
	{
		cout << "file:" << ecname <<": could not be opened" << endl;
		cout << "no elastic collision xs data for " << name << endl;
	}
	for (numlines = 0; std::getline(fi, line); ++numlines);

	bheecxs = new double*[numlines];
	for(int n=0;n<numlines;n++) bheecxs[n] = new double[2];

	idx = 0;
	fi.clear(); fi.seekg(0);
	while(!fi.eof() && !fi.fail())
	{
		if(getline(fi,line) )
		{
			myiolib::jString jline(line);
			//jline = jline.clean();
			if(jline.length() > 0) {
				tokens = jline.clean().parse();

				bheecxs[idx][0] = tokens[0].dValue();
				bheecxs[idx][1] = tokens[1].dValue();

				tokens.clear();
				idx++;
			}
		}
	}
	fi.close();
	tokens.clear();

	if (idx == numlines) {
		Neecxs = idx;
		printf("Completed Elastic Collision Data\n");
		printf("----Number of ec xs: %i\n",Neecxs);
	}
	else {
		printf("no of ec xs and no of txtlines don't match... check for empty txtlines\n");
		abort();
	}
}

// Set ionization potential
void AtomX::setEi(float Ip)
{
	for(int n=0;n<Nlevels;n++)
		Ei[n] = Ip-Ee[n];
}


// Set mass
void AtomX::setmass(float m)
{
	amu = m;
}

// compute memory size
double AtomX::memory_size()
{
	double _msize = 0.;
	_msize += (double)Nlevels*2*sizeof(int);
	for(int i=0;i<Nlevels;i++)
		_msize += (double)nS[i]*3*sizeof(int);
	_msize += (double)Nlevels*2*sizeof(float);

	_msize += (double)NBBlns*(4*sizeof(float)+7*sizeof(int));
	_msize += (double)NBBlns*(4*sizeof(float));

	_msize += (double)NBBeXD*(4*sizeof(float)+7*sizeof(int));
	_msize += (double)NBBeXD*3*Nexdxs*sizeof(float);

	_msize += (double)NBFeIR*(4*sizeof(float)+7*sizeof(int));
	_msize += (double)NBFeIR*3*Neirxs*sizeof(float);

	_msize += (double)NBFpIR*(4*sizeof(float)+7*sizeof(int));
	_msize += (double)NBFpIR*3*Npirxs*sizeof(float);
	return _msize*1.e-6;
}

// Rates
double AtomX::eexc_rate(int bb, float TeV)
{
	if (bb >= NBBeXD) {
		printf("!!!AtomX::eexc_rate: transition out of range, returning 0\n");
		return 0.;
	}

	// in MKS units
    double omega_H = 13.6/TeV;
    double ve_bar = 2.468e+06/sqrt(omega_H);
    int Npnts = Nexdxs;
    double estar = bbexd[bb].dE/TeV;
    double sum_int = 0.;
    for(int ke=1;ke<Npnts;ke++)
	{
		double e0  = bbeexxs[bb][ke-1][0]/TeV;
		double e1  = bbeexxs[bb][ke  ][0]/TeV;
		double xs0 = bbeexxs[bb][ke-1][1];
		double xs1 = bbeexxs[bb][ke  ][1];
		double de = e1-e0;
		double dfe = 0.;

		// option 1
		dfe += .5 * de * e0 * xs0 * exp(-e0+estar);
		dfe += .5 * de * e1 * xs1 * exp(-e1+estar);

		// option 2
//		double ec  = exp(0.5*(log(e0)+log(e1)));
//		double xsc = 0.5*(xs0+xs1);
//		dfe += de * xsc * ec * exp(-ec/TeV);

		sum_int = sum_int+dfe;
	}
    return ve_bar*sum_int;
}

double AtomX::eion_rate(int bf, float TeV)
{
	if (bf >= NBFeIR) {
		printf("!!!AtomX::eion_rate: transition out of range, returning 0\n");
		return 0.;
	}

	// in MKS units
    double omega_H = 13.6/TeV;
    double ve_bar = 2.468e+06/sqrt(omega_H);
    int Npnts = Neirxs;
    double estar = bfeir[bf].dE/TeV;
    double sum_int = 0.;
    for(int ke=1;ke<Npnts;ke++)
	{
		double e0  = bfeioxs[bf][ke-1][0]/TeV;
		double e1  = bfeioxs[bf][ke  ][0]/TeV;
		double xs0 = bfeioxs[bf][ke-1][1];
		double xs1 = bfeioxs[bf][ke  ][1];
		double de = e1-e0;
		double dfe = 0.;

		// option 1
		dfe += .5 * de * e0 * xs0 * exp(-e0+estar);
		dfe += .5 * de * e1 * xs1 * exp(-e1+estar);

		// option 2
//		double ec  = exp(0.5*(log(e0)+log(e1)));
//		double xsc = 0.5*(xs0+xs1);
//		dfe += de * xsc * ec * exp(-ec/TeV);

		sum_int = sum_int+dfe;
	}
    return ve_bar*sum_int;
}

double* AtomX::prec_rate(int bf, float TeV)
{
	if (bf >= NBFpIR) {
		printf("!!!AtomX::prec_rate: transition out of range, returning NULL\n");
		return NULL;
	}

	double* rate = new double[2];
	int l = bfpir[bf].low[1];
	int u = bfpir[bf].upp[1];
	double grat= (double)ge[l]/2./(double)ion->ge[u];
	double Ip = bfpir[bf].dE/TeV;
	double Ee = 0.510999e6/TeV; // electron rest mass energy

	// in MKS units
    double omega_H = 13.6/TeV;
    double ve_bar = 2.468e+06/sqrt(omega_H);
    int Npnts = Npirxs;
    double sum_int = 0.;
    double sum_int2 = 0.;

    double sigma0 = bfpioxs[bf][0][1];
    double dee= (bfpioxs[bf][0][0]-bfpir[bf].dE)/TeV;
    double ee = dee/2.;
    double xss = grat*(ee+Ip)*(ee+Ip)/(Ee*ee)*sigma0;
//    sum_int  += dee * ee * xss * exp(-ee);
//    sum_int2 += dee * ee * ee * xss * exp(-ee);

    for(int ke=1;ke<Npnts;ke++)
	{
		double e0  = bfpioxs[bf][ke-1][0]/TeV;
		double e1  = bfpioxs[bf][ke  ][0]/TeV;
		double xs0_pi = bfpioxs[bf][ke-1][1];
		double xs1_pi = bfpioxs[bf][ke  ][1];
		double xs0 = grat*e0*e0/(Ee*(e0-Ip))*xs0_pi;
		double xs1 = grat*e1*e1/(Ee*(e1-Ip))*xs1_pi;
		double de = e1-e0;
		double dfe = 0.;
		double dfe2 = 0.;

		// option 1
		dfe += .5 * de * (e0-Ip) * xs0 * exp(-(e0-Ip));
		dfe += .5 * de * (e1-Ip) * xs1 * exp(-(e1-Ip));

		dfe2 += .5 * de * (e0-Ip) * (e0-Ip) * xs0 * exp(-(e0-Ip));
		dfe2 += .5 * de * (e1-Ip) * (e1-Ip) * xs1 * exp(-(e1-Ip));

		sum_int  = sum_int+dfe;
		sum_int2 = sum_int2+dfe2;
	}
    rate[0] = 1.e+0*ve_bar*sum_int;
    rate[1] = ve_bar*sum_int2;
    return rate;
}

double AtomX::ecol_rate(float TeV)
{
//	if (bf >= NBFeIR) {
//		printf("!!!AtomX::eion_rate: transition out of range, returning 0\n");
//		return 0.;
//	}

	// in MKS units
    double omega_H = 13.6/TeV;
    double ve_bar = 2.468e+06/sqrt(omega_H);
//    int Npnts = Neecxs;
    int Npnts = 200;
    double interv = 6./((double) Npnts);
    double sum_int = 0.;
    for(int ke=1;ke<Npnts;ke++)
	{
//    	double e0  = bheecxs[ke-1][0]/TeV;
//		double e1  = bheecxs[ke  ][0]/TeV;
//		double xs0 = bheecxs[ke-1][1];
//		double xs1 = bheecxs[ke  ][1];

		double e0  = pow(10.,-3.+(ke-1)*interv);
		double e1  = pow(10.,-3.+(ke)*interv);
		double xs0 = ecxs_calc(e0);
		double xs1 = ecxs_calc(e1);
		e0 = e0/TeV; e1 = e1/TeV;

		double de = e1-e0;
		double dfe = 0.;

		// option 1
		dfe += .5 * de * e0 * e0 * xs0 * exp(-e0);
		dfe += .5 * de * e1 * e1 * xs1 * exp(-e1);

		// option 2
//		double ec  = exp(0.5*(log(e0)+log(e1)));
//		double xsc = 0.5*(xs0+xs1);
//		dfe += de * xsc * ec * exp(-ec/TeV);

		sum_int = sum_int+dfe;
	}
    return ve_bar*sum_int;
}

/**
 * Elastic collision cross section (Analytical Fit)
 * Argon cross sections
 *
 * Published for input data ranging from 1e-3 eV to 1e3 eV
 *
 * McEachran and Stauffer (2014)
 */
double AtomX::ecxs_calc(double E)
{
	double Ehart = E/AU_Hart;
	double k = sqrt(Ehart*(MKS_fstr*MKS_fstr*Ehart+2));

	double xs = 0.;
	if (k >= 0 && k < 0.907)
	{
		double a[6] = {6.53,-146.99,198.67,-325.50,265.56,552.07};
		double b[6] = {-10.68,-501.88,-187.37,76.11,-740.25,455.44};

		double num = a[0] + a[1]*k + a[2]*k*k + a[3]*k*k*log(k) + a[4]*k*k*k + \
				a[5]*k*k*k*log(k);
		double den = 1 + b[0]*k + b[1]*k*k + b[2]*k*k*log(k) + b[3]*k*k*k + \
				b[4]*k*k*k*log(k) + b[5]*k*k*k*k;

		xs = num/den;
	}
	else if (k >= 0.907 && k < 1.115)
	{
		double a[6] = {-321.1,1028.72,-1016.54,324.47,0.,0.};
		double b[6] = {0.,0.,0.,0.,0.,0.};

		double num = a[0] + a[1]*k + a[2]*k*k + a[3]*k*k*k + a[4]*k*k*k*k + a[5]*k*k*k*k*k;
		double den = 1 + b[0]*k + b[1]*k*k + b[2]*k*k*k + b[3]*k*k*k*k + b[4]*k*k*k*k*k + b[5]*k*k*k*k*k*k;

		xs = num/den;
	}
	else if (k >= 1.115 && k <= 8.577)
	{
		double a[6] = {433.19,-352.94,123.70,-4.97,0.,0.};
		double b[6] = {-26.76,57.14,-25.85,6.01,0.,0.};

		double num = a[0] + a[1]*k + a[2]*k*k + a[3]*k*k*k + a[4]*k*k*k*k + a[5]*k*k*k*k*k;
		double den = 1 + b[0]*k + b[1]*k*k + b[2]*k*k*k + b[3]*k*k*k*k + b[4]*k*k*k*k*k + b[5]*k*k*k*k*k*k;

		xs = num/den;
	}
	else
	{
		printf("Out of bounds for elastic collision rate integration: E = %3.2e!\n",E);
		exit(1);
	}

	return xs*1e-20;
}

