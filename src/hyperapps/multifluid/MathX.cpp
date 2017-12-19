#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <string.h>
#include <omp.h>
#include <iostream>

#include "MathX.h"

using namespace std;
using namespace Eigen;

double MathX::gamma(double a)
{
    double gam = 0.;
    double x = a-1.;
    double y = x+5.5;
    if(y > 0.)
    {
        double z = 1.+7.618009173e+01/(x+1.)
                     -8.650532033e+01/(x+2.)
                     +2.401409822e+01/(x+3.)
                     -1.231739516e-00/(x+4.)
                     -1.208580030e-03/(x+5.)
                     -0.536382000e-05/(x+6.)  ;
        if(z > 0.)
        {   gam = (x+0.5)*log(y)-y + log(2.50662827465 * z); }
        else
        {   cout << "invalid argument for MathX::gamma function, z<0., is " << z << endl; exit(0); }
    }
    else
    {   cout << "invalid argument for MathX::gamma function, y<0., is " << y << endl; exit(0); }
    return gam;
}

    /**
    * Returns the error function.
    * @param   arg
    * @return  erf
    */
double MathX::erf( double a )
{
    double val = 1.;
    if(a >=0.)
    {   double b = 1./(1.+0.3275911*a);
        val = 1. - b*( 0.254829592
                 + b*(-0.284496736
                 + b*( 1.421413741
                 + b*(-1.453152027 + b*1.061405429) ) ) ) * exp(-a*a);
    }
    else
    if(a < 0.)
    {   double b = 1./(1.-0.3275911*a);
        val =-1. + b*( 0.254829592
                 + b*(-0.284496736
                 + b*( 1.421413741
                 + b*(-1.453152027 + b*1.061405429) ) ) ) * exp(-a*a);
    }
    return (val);
}
    /**
    * Returns the complement of error function.
    * @param   a argument
    * @return  erfc
    */
double MathX::erfc( double a )
{
    return (1.0-MathX::erf(a));
}

double MathX::En(int n, double a, double b)
{
    int oQuad = 2;
    int Nquad = 0;
    double sum0 = 1;
    double error = 1.e+20;
    while(error > 1.e-08)
    {
        kMatrix2D<double>* _coef;
        _coef = MathX::gaussQuadNodes(oQuad++,a,b);
        kMatrix2D<double>& coef = *(_coef);
        Nquad = coef.length(1);
        double sum = 0.;
        for(int j=0;j<Nquad;j++)  sum += coef(1,j)*exp(-coef(0,j))/pow(coef(0,j),n);
        error = abs(sum/sum0-1.);
        sum0 = sum;
        delete _coef;
    }
    return sum0;
}

double MathX::Psi(double xn)
{
    int Npnts = 1000;
    double x0 = xn;
    double x1 = xn;
    double dx = (100.-xn)/Npnts;
    double sum = 0.;
    if(dx < 0.) return sum;
    for(int i=1;i<=Npnts;i++)
    {
        x1 = x0+dx;
        sum -= MathX::En(1,x0,x1);
        x0 = x1;
    }
    sum += (exp(-xn)/xn);
    return sum;
}

//double MathX::Psi(double xn)
//{
//    int Npnts = 1000;
//    double x0 = xn;
//    double x1 = xn;
//    double dx = (100.-xn)/Npnts;
//    double sum = 0.;
//    if(dx < 0.) return sum;
//
//    int i,j,k;
//    int oQuad = 4;
//    GaussQuadrature& gauss= *(new GaussQuadrature(oQuad));
//    int Nquad = gauss.Nq;
//    double* xi = new double[Nquad];
//    double* wg = new double[Nquad];
//    double xc,E1;
//
//    for(i=1;i<=Npnts;i++)
//    {
//        x1 = x0+dx;
//        xc = (x0+x1)/2.;
//        dx = (x1-x0)/2.;
//        for(k=0;k<Nquad;k++)
//        {
//            xi[k] = xc + gauss.qcoo[k] * dx;
//            wg[k] =      gauss.qwgt[k] * dx;
//        }
//        E1 = 0.; for(j=0;j<Nquad;j++)  E1 += wg[j]*exp(-xi[j])/xi[j];
//        sum -= E1;
//        x0 = x1;
//    }
//    sum += (exp(-xn)/xn);
//    delete[] xi;
//    delete[] wg;
//    delete &gauss;
//
//    return sum;
//}

double MathX::expIntegral(int n, double z)
{
	double epsilon = 0.001;
	double w = 1.;
	double dw, df;
	double xc = 1.;
	double sum = 0.;
	while(xc < 100.)
	{
		dw = epsilon*w;
		xc = w*z;
		df = exp(-xc)*dw/pow(w,n);
		sum += df;
		w += dw;
	}
	return sum;
}


kMatrix2D<double>* MathX::gaussQuadNodes(int order)
{   return MathX::gaussQuadNodes(order,-1.,1.);   }

kMatrix2D<double>* MathX::gaussQuadNodes(int order, double a, double b)
{
    int n = 1; for(int i=0;i<order;i++) n *= 2;
    double* qcoo = new double[n];
    double* qwgt = new double[n];
    double eps = 1.e-15;

    if(order==0)
    {
        qcoo[ 0] = 0.000000000000000;  qwgt[ 0] = 2.000000000000000;
    }
    else
    if(order==1)
    {
        qcoo[ 0] =-0.577350269189626;  qwgt[ 0] = 1.000000000000000;
        qcoo[ 1] = 0.577350269189626;  qwgt[ 1] = 1.000000000000000;
    }
    else
    if(order==2)
    {
        qcoo[ 0] =-0.339981043584856;  qwgt[ 0] = 0.652145154862546;
        qcoo[ 1] = 0.339981043584856;  qwgt[ 1] = 0.652145154862546;
        qcoo[ 2] =-0.861136311594053;  qwgt[ 2] = 0.347854845137454;
        qcoo[ 3] = 0.861136311594053;  qwgt[ 3] = 0.347854845137454;
    }
    else
    if(order==3)
    {
        qcoo[ 0] =-0.183434642495650;  qwgt[ 0] = 0.362683783378362;
        qcoo[ 1] = 0.183434642495650;  qwgt[ 1] = 0.362683783378362;
        qcoo[ 2] =-0.525532409916329;  qwgt[ 2] = 0.313706645877887;
        qcoo[ 3] = 0.525532409916329;  qwgt[ 3] = 0.313706645877887;
        qcoo[ 4] =-0.796666477413627;  qwgt[ 4] = 0.222381034453374;
        qcoo[ 5] = 0.796666477413627;  qwgt[ 5] = 0.222381034453374;
        qcoo[ 6] =-0.960289856497536;  qwgt[ 6] = 0.101228536290376;
        qcoo[ 7] = 0.960289856497536;  qwgt[ 7] = 0.101228536290376;
    }
    else
    if(order==4)
    {
        qcoo[ 0] =-0.09501250983764;  qwgt[ 0] = 0.18945061045507;
        qcoo[ 1] = 0.09501250983764;  qwgt[ 1] = 0.18945061045507;
        qcoo[ 2] =-0.28160355077926;  qwgt[ 2] = 0.18260341504492;
        qcoo[ 3] = 0.28160355077926;  qwgt[ 3] = 0.18260341504492;
        qcoo[ 4] =-0.45801677765723;  qwgt[ 4] = 0.16915651939500;
        qcoo[ 5] = 0.45801677765723;  qwgt[ 5] = 0.16915651939500;
        qcoo[ 6] =-0.61787624440264;  qwgt[ 6] = 0.14959598881658;
        qcoo[ 7] = 0.61787624440264;  qwgt[ 7] = 0.14959598881658;
        qcoo[ 8] =-0.75540440835500;  qwgt[ 8] = 0.12462897125553;
        qcoo[ 9] = 0.75540440835500;  qwgt[ 9] = 0.12462897125553;
        qcoo[10] =-0.86563120238783;  qwgt[10] = 0.09515851168249;
        qcoo[11] = 0.86563120238783;  qwgt[11] = 0.09515851168249;
        qcoo[12] =-0.94457502307323;  qwgt[12] = 0.06225352393865;
        qcoo[13] = 0.94457502307323;  qwgt[13] = 0.06225352393865;
        qcoo[14] =-0.98940093499165;  qwgt[14] = 0.02715245941175;
        qcoo[15] = 0.98940093499165;  qwgt[15] = 0.02715245941175;
    }
    else
    if(order==5)
    {
        qcoo[ 0] =-0.04830766568774;  qwgt[ 0] = 0.09654008851473;
        qcoo[ 1] = 0.04830766568774;  qwgt[ 1] = 0.09654008851473;
        qcoo[ 2] =-0.14447196158280;  qwgt[ 2] = 0.09563872007927;
        qcoo[ 3] = 0.14447196158280;  qwgt[ 3] = 0.09563872007927;
        qcoo[ 4] =-0.23928736225214;  qwgt[ 4] = 0.09384439908080;
        qcoo[ 5] = 0.23928736225214;  qwgt[ 5] = 0.09384439908080;
        qcoo[ 6] =-0.33186860228213;  qwgt[ 6] = 0.09117387869576;
        qcoo[ 7] = 0.33186860228213;  qwgt[ 7] = 0.09117387869576;
        qcoo[ 8] =-0.42135127613064;  qwgt[ 8] = 0.08765209300440;
        qcoo[ 9] = 0.42135127613064;  qwgt[ 9] = 0.08765209300440;
        qcoo[10] =-0.50689990893223;  qwgt[10] = 0.08331192422695;
        qcoo[11] = 0.50689990893223;  qwgt[11] = 0.08331192422695;
        qcoo[12] =-0.58771575724076;  qwgt[12] = 0.07819389578707;
        qcoo[13] = 0.58771575724076;  qwgt[13] = 0.07819389578707;
        qcoo[14] =-0.66304426693022;  qwgt[14] = 0.07234579410885;
        qcoo[15] = 0.66304426693022;  qwgt[15] = 0.07234579410885;
        qcoo[16] =-0.73218211874029;  qwgt[16] = 0.06582222277636;
        qcoo[17] = 0.73218211874029;  qwgt[17] = 0.06582222277636;
        qcoo[18] =-0.79448379596794;  qwgt[18] = 0.05868409347854;
        qcoo[19] = 0.79448379596794;  qwgt[19] = 0.05868409347854;
        qcoo[20] =-0.84936761373257;  qwgt[20] = 0.05099805926238;
        qcoo[21] = 0.84936761373257;  qwgt[21] = 0.05099805926238;
        qcoo[22] =-0.89632115576605;  qwgt[22] = 0.04283589802223;
        qcoo[23] = 0.89632115576605;  qwgt[23] = 0.04283589802223;
        qcoo[24] =-0.93490607593774;  qwgt[24] = 0.03427386291302;
        qcoo[25] = 0.93490607593774;  qwgt[25] = 0.03427386291302;
        qcoo[26] =-0.96476225558751;  qwgt[26] = 0.02539206530926;
        qcoo[27] = 0.96476225558751;  qwgt[27] = 0.02539206530926;
        qcoo[28] =-0.98561151154527;  qwgt[28] = 0.01627439473091;
        qcoo[29] = 0.98561151154527;  qwgt[29] = 0.01627439473091;
        qcoo[30] =-0.99726386184948;  qwgt[30] = 0.00701861000947;
        qcoo[31] = 0.99726386184948;  qwgt[31] = 0.00701861000947;
    }
    else
    if(order==6)
    {
        qcoo[ 0] =-0.02435029266342;  qwgt[ 0] = 0.04869095700914;
        qcoo[ 1] = 0.02435029266342;  qwgt[ 1] = 0.04869095700914;
        qcoo[ 2] =-0.07299312178780;  qwgt[ 2] = 0.04857546744150;
        qcoo[ 3] = 0.07299312178780;  qwgt[ 3] = 0.04857546744150;
        qcoo[ 4] =-0.12146281929612;  qwgt[ 4] = 0.04834476223480;
        qcoo[ 5] = 0.12146281929612;  qwgt[ 5] = 0.04834476223480;
        qcoo[ 6] =-0.16964442042399;  qwgt[ 6] = 0.04799938859646;
        qcoo[ 7] = 0.16964442042399;  qwgt[ 7] = 0.04799938859646;
        qcoo[ 8] =-0.21742364374001;  qwgt[ 8] = 0.04754016571483;
        qcoo[ 9] = 0.21742364374001;  qwgt[ 9] = 0.04754016571483;
        qcoo[10] =-0.26468716220877;  qwgt[10] = 0.04696818281621;
        qcoo[11] = 0.26468716220877;  qwgt[11] = 0.04696818281621;
        qcoo[12] =-0.31132287199021;  qwgt[12] = 0.04628479658131;
        qcoo[13] = 0.31132287199021;  qwgt[13] = 0.04628479658131;
        qcoo[14] =-0.35722015833767;  qwgt[14] = 0.04549162792742;
        qcoo[15] = 0.35722015833767;  qwgt[15] = 0.04549162792742;
        qcoo[16] =-0.40227015796399;  qwgt[16] = 0.04459055816376;
        qcoo[17] = 0.40227015796399;  qwgt[17] = 0.04459055816376;
        qcoo[18] =-0.44636601725346;  qwgt[18] = 0.04358372452932;
        qcoo[19] = 0.44636601725346;  qwgt[19] = 0.04358372452932;
        qcoo[20] =-0.48940314570705;  qwgt[20] = 0.04247351512365;
        qcoo[21] = 0.48940314570705;  qwgt[21] = 0.04247351512365;
        qcoo[22] =-0.53127946401989;  qwgt[22] = 0.04126256324262;
        qcoo[23] = 0.53127946401989;  qwgt[23] = 0.04126256324262;
        qcoo[24] =-0.57189564620263;  qwgt[24] = 0.03995374113272;
        qcoo[25] = 0.57189564620263;  qwgt[25] = 0.03995374113272;
        qcoo[26] =-0.61115535517239;  qwgt[26] = 0.03855015317862;
        qcoo[27] = 0.61115535517239;  qwgt[27] = 0.03855015317862;
        qcoo[28] =-0.64896547125466;  qwgt[28] = 0.03705512854024;
        qcoo[29] = 0.64896547125466;  qwgt[29] = 0.03705512854024;
        qcoo[30] =-0.68523631305423;  qwgt[30] = 0.03547221325688;
        qcoo[31] = 0.68523631305423;  qwgt[31] = 0.03547221325688;
        qcoo[32] =-0.71988185017161;  qwgt[32] = 0.03380516183714;
        qcoo[33] = 0.71988185017161;  qwgt[33] = 0.03380516183714;
        qcoo[34] =-0.75281990726053;  qwgt[34] = 0.03205792835485;
        qcoo[35] = 0.75281990726053;  qwgt[35] = 0.03205792835485;
        qcoo[36] =-0.78397235894334;  qwgt[36] = 0.03023465707240;
        qcoo[37] = 0.78397235894334;  qwgt[37] = 0.03023465707240;
        qcoo[38] =-0.81326531512280;  qwgt[38] = 0.02833967261426;
        qcoo[39] = 0.81326531512280;  qwgt[39] = 0.02833967261426;
        qcoo[40] =-0.84062929625258;  qwgt[40] = 0.02637746971505;
        qcoo[41] = 0.84062929625258;  qwgt[41] = 0.02637746971505;
        qcoo[42] =-0.86599939815409;  qwgt[42] = 0.02435270256871;
        qcoo[43] = 0.86599939815409;  qwgt[43] = 0.02435270256871;
        qcoo[44] =-0.88931544599511;  qwgt[44] = 0.02227017380838;
        qcoo[45] = 0.88931544599511;  qwgt[45] = 0.02227017380838;
        qcoo[46] =-0.91052213707850;  qwgt[46] = 0.02013482315353;
        qcoo[47] = 0.91052213707850;  qwgt[47] = 0.02013482315353;
        qcoo[48] =-0.92956917213194;  qwgt[48] = 0.01795171577570;
        qcoo[49] = 0.92956917213194;  qwgt[49] = 0.01795171577570;
        qcoo[50] =-0.94641137485840;  qwgt[50] = 0.01572603047602;
        qcoo[51] = 0.94641137485840;  qwgt[51] = 0.01572603047602;
        qcoo[52] =-0.96100879965205;  qwgt[52] = 0.01346304789672;
        qcoo[53] = 0.96100879965205;  qwgt[53] = 0.01346304789672;
        qcoo[54] =-0.97332682778991;  qwgt[54] = 0.01116813946013;
        qcoo[55] = 0.97332682778991;  qwgt[55] = 0.01116813946013;
        qcoo[56] =-0.98333625388463;  qwgt[56] = 0.00884675982636;
        qcoo[57] = 0.98333625388463;  qwgt[57] = 0.00884675982636;
        qcoo[58] =-0.99101337147674;  qwgt[58] = 0.00650445796898;
        qcoo[59] = 0.99101337147674;  qwgt[59] = 0.00650445796898;
        qcoo[60] =-0.99634011677196;  qwgt[60] = 0.00414703326056;
        qcoo[61] = 0.99634011677196;  qwgt[61] = 0.00414703326056;
        qcoo[62] =-0.99930504173577;  qwgt[62] = 0.00178328072170;
        qcoo[63] = 0.99930504173577;  qwgt[63] = 0.00178328072170;
    }

    else
    {   // general case
        double* u = new double[n+1];
        int m = (n+1)/2;
        int ia = 0;
        for(int k=m-1;k>=0;k--)
        {
            double* u = new double[n+1];
            double x = cos(M_PI*(k+0.75)/(n+0.5));
            double r = 0.;
            double w = 0.;
            double dpdx=0.;
            double err = 1.e+20;
            while(err > eps)
            {   u[0] = 1.;
                u[1] = x;
                for(int j=1;j<n;j++) u[j+1] = ( (2.*j+1.)*x*u[j] - j*u[j-1])/(j+1.);
                dpdx = n*(x*u[n]-u[n-1])/(x*x-1.);
                r = x;
                x = r-u[n]/dpdx;
                err = abs(x-r);
            }
            dpdx = n*(x*u[n]-u[n-1])/(x*x-1.);
            w = 2./((1.-x*x)*dpdx*dpdx);
            qcoo[ia] = - x;  qwgt[ia] = w;  ia++;
            qcoo[ia] = + x;  qwgt[ia] = w;  ia++;
        }
        delete[] u;
    }
    // sort in ascending order
    double* xi = new double[n];
    double* wg = new double[n];
    int j = 0;
    for(int k=n-2;k>=0;k-=2,j++)
    {
        xi[j] = qcoo[k];
        wg[j] = qwgt[k];
    }
    for(int k=1;k<n;k+=2,j++)
    {
        xi[j] = qcoo[k];
        wg[j] = qwgt[k];
    }

    // rescale nodes and weights
    kMatrix2D<double>* _coef = new kMatrix2D<double>(2,n);
    kMatrix2D<double>& coef = *(_coef);
    double xc = (b+a)/2.;
    double dx = (b-a)/2.;
    for(int k=0;k<n;k++)
    {
        coef(0,k) = xc + xi[k] * dx;
        coef(1,k) =      wg[k] * dx;
    }

    delete[] xi;
    delete[] wg;
    delete[] qcoo;
    delete[] qwgt;
    /**
     *  returns coef[0] = { array of nodes }
     *  returns coef[1] = { array of weights }
     */
    return _coef;
}

// Legendre Polynomials
double** MathX::LegendrePolyCoefs(int order)
{
    int Np = max(order,1);
    double** c = new double*[Np+1];
    // zero-th order
    c[0] = new double[1]; c[0][0] = 1.;
    // 1st-order
    c[1] = new double[2]; c[1][1] = 1.; c[1][0] = 0.;
    // higher orders
    for(int k=2;k<=Np;k++)
    {
        c[k] = new double[k+1];
        double a = (2.*k-1.)/k;
        double b = (1.*k-1.)/k;
        for(int i=0;i<k  ;i++) c[k][i+1] += a*c[k-1][i];
        for(int i=0;i<k-1;i++) c[k][i  ] -= b*c[k-2][i] ;
    }
    return c;
}

double MathX::LegendrePolynomial(int Np, double* cp, double z)
{
    double p1 = cp[Np];
    for(int n=Np-1;n>=0;n--) p1 = z*p1+cp[n];
    return p1;
}

double* MathX::LegendrePolynomial(int Np, double* cp, int Nz, double* z)
{
    double* p = new double[Nz];
    for(int i=0;i<Nz;i++) p[i] = cp[Np];
    for(int n=Np-1;n>=0;n--)
    {
        for(int i=0;i<Nz;i++) p[i] = z[i]*p[i]+cp[n];
    }
    return p;
}

// Natural splines
double* MathX::getSplineCoefs(int Npnts, double* x, double* y)
{
    double* a = new double[Npnts];
    double* b = new double[Npnts];
    double* d = new double[Npnts];
    double* r = new double[Npnts];

    Npnts--;
    int is = 0;
    int ie = Npnts;
    for(int i=0;i<=Npnts;i++)
    {
        a[i] = 1.;
        d[i] = 4.;
        b[i] = 1.;
        r[i] = y[i];
    }
    double dy_is = (y[is+1]-y[is]);
    double dy_ie = (y[ie]-y[ie-1]);
    a[is] = 0.; b[is] = 2.; r[is] += (dy_is/3.);
    a[ie] = 2.; b[ie] = 0.; r[ie] -= (dy_ie/3.);

    for(int i=is+1;i<=ie;i++)
    {    double frac = a[i]/d[i-1];
         d[i] -= frac*b[i-1];
         r[i] -= frac*r[i-1];
    }
    for(int i=is;i<=ie;i++)
    {    b[i] /= d[i];
         r[i] /= d[i];
    }
    for(int i=ie-1;i>=is;i--)
    {   r[i] -=(b[i] * r[i+1]);  }

    double* c = new double[Npnts+3];
    int ishft = 1;
    for(int i=is;i<=ie;i++) c[i+ishft] = r[i];
    c[ie+1+ishft] = c[ie-1+ishft] +(dy_ie/3.);
    c[is-1+ishft] = c[is+1+ishft] -(dy_is/3.);

    return c;
}

double* computeSpline(int Ndata, double* x, double* xp, double* cB)
{
    double* q = new double[Ndata];
    int ishft = 1;

    double xi = 0.;
    double dx = xp[1]-xp[0]; // assumes constant spacing !!
    for(int j=0;j<Ndata;j++)
    {
        double dz = x[j]/dx;
        int i = (int) dz;  dz -= i;
        double cm1 = cB[i-1+ishft];
        double cm0 = cB[i  +ishft];
        double cp1 = cB[i+1+ishft];
        double cp2 = cB[i+2+ishft];
        xi = 1.+ dz; double bm1 = (2.-xi)*(2.-xi)*(2.-xi);
        xi = 0.+ dz; double bm0 = 1. + 3.*(1.-xi)*( 1.+ (1.-xi) *(1.-(1.-xi) ) );
        xi =-1.+ dz; double bp1 = 1. + 3.*(1.+xi)*( 1.+ (1.+xi) *(1.-(1.+xi) ) );
        xi =-2.+ dz; double bp2 = (2.+xi)*(2.+xi)*(2.+xi);
        q[j] = cm1*bm1 + cm0*bm0 + cp1*bp1 + cp2*bp2;
    }
    return q;
}


void MathX::Gauss_Jordan_point(int Ndim, double* sol, double** jac, double* rhs)
{
    int n,np,ns,ms;
    int Nmax = Ndim-1;

    double errP;

    // add-up excitations and deexcitations
    for(ns=0;ns<=Nmax;ns++) sol[ns]=rhs[ns];

    for(np=0;np< Nmax;np++)
    {
        for(ms=np+1;ms<=Nmax;ms++)
        {   errP = jac[ms][np]/jac[np][np];
            sol[ms] -= errP*sol[np];
            for(ns=np;ns<=Nmax;ns++) jac[ms][ns] -= (errP*jac[np][ns]);
        }
    }
    sol[Nmax] = sol[Nmax]/jac[Nmax][Nmax];

    for(np=Nmax-1;np>=0;np--)
    {
        errP=0.;
        for(ns=np+1;ns<=Nmax;ns++) errP += (jac[np][ns]*sol[ns]);
        sol[np] = (sol[np]-errP)/jac[np][np] ;
    }
} // end Point-Implicit Solve



/**
 * LU decomposition of a square matrix using Gaussian elimination.
 * L and U are computed with the "daxpy"-based
 * elimination algorithm used in LINPACK and MATLAB.
 *
 * @param  X        matrix with dimensions <code>[m][m]</code>
 * @param  piv      vector to store pivot info with dimension <code>[m]</code>
 * @param  pivsign  wrapper for the pivot sign
 */
void MathX::LU_decomp(int N, double** X, int* piv)
{
    int i,j,k;
    for(i=0;i<N;i++) piv[i] = i;

    int pivsign = 1;
    // Main loop.
    for(k=0;k<N;k++)
    {   // Find pivot.
        int p = k;
        for (i=k+1;i<N;i++)
        {
            if(abs(X[i][k]) > abs(X[p][k]))  p = i;
        }
        // Exchange if necessary.
        if (p != k)
        {
            for(j=0;j<N;j++)
            {
                double t = X[p][j];
                X[p][j] = X[k][j];
                X[k][j] = t;
            }
            int t = piv[p];
            piv[p] = piv[k];
            piv[k] = t;
            pivsign = -pivsign;
        }

        // Compute multipliers and eliminate k-th column.
        if(X[k][k] != 0.)
        {
            for(i=k+1;i<N;i++)
            {
                double* Xi = X[i];
                double* Xk = X[k];
                Xi[k] /= Xk[k];
                for(j=k+1;j<N;j++)  Xi[j] -= Xi[k]*Xk[j];
            }
        }

    }
}

/**
 * Matrix inversion.
 *
 * @param  A   Matrix with dimensions <code>[n][n]</code>
 * @return     inv(A)
 */
double** MathX::invert(int N, double** A)
{
    int i,j,k;
    double** X = new double*[N]; for(i=0;i<N;i++) X[i] = new double[N];
    int* piv = new int[N];

    LU_decomp(N,A,piv);

    for(i=0;i<N;i++) X[i][piv[i]] = 1.;

    // Solve L*Y = B(piv,:)
    for(k=0;k<N;k++)
    {
        double* Xk = X[k];
        for(i=k+1;i<N;i++)
        {
            double* Xi = X[i];
            double* Ai = A[i];
            for(j=0;j<N;j++)  Xi[j] -= Xk[j]*Ai[k];
        }
    }
    // Solve U*X = Y;
    for(k=N-1;k>=0;k--)
    {
        double* Xk = X[k];
        double* Ak = A[k];
        for(j=0;j<N;j++) Xk[j] /= Ak[k];
        for(i=0;i<k;i++)
        {
            double* Xi = X[i];
            double* Ai = A[i];
            for(j=0;j<N;j++) Xi[j] -= Xk[j]*Ai[k];
        }
    }
    return X;
}

//void MathX::EigenDecomp(int N, double** A, double** Teig, double** Tinv, double* eVals)
//{
//    int i,j,k;
//
//    MatrixXd mat(N,N);
//    for(i=0;i<N;i++)
//    for(j=0;j<N;j++) mat(i,j) = A[i][j];
//
//    EigenSolver<MatrixXd> es;
//    es.compute(mat,/* compute eigenvectors ? */ true);
//    if (es.info() != Success) abort();
//    for(int i=0;i<N;i++)
//    {
//        complex<double> lambda = es.eigenvalues()(i);
//        eVals[i] = lambda.real();
//    }
//    MatrixXd evc(N,N);
//    MatrixXd inv(N,N);
//    for(int i=0;i<N;i++)
//    {
//        for(int j=0;j<N;j++)
//        {
//            complex<double> v = es.eigenvectors()(i,j);
//            evc(i,j) = v.real();
//        }
//    }
//    inv = evc.inverse();
//
//    for(int i=0;i<N;i++)
//    {
//        for(int j=0;j<N;j++)
//        {
//            Teig[i][j] = evc(i,j);
//            Tinv[i][j] = inv(i,j);
//        }
//    }
//}

//void MathX::EigenDecomp(int N, double** A, complex<double>** Teig, complex<double>** Tinv, complex<double>* eVals)
//{
//    int i,j,k;
//
//    MatrixXd mat(N,N);
//    for(i=0;i<N;i++)
//    for(j=0;j<N;j++) mat(i,j) = A[i][j];
//
//
//    EigenSolver<MatrixXd> es;
//    es.compute(mat,/* compute eigenvectors ? */ true);
//    if (es.info() != Success) abort();
//
//    for(int i=0;i<N;i++)
//    {
//        eVals[i] = es.eigenvalues()(i);
//    }
//    MatrixXcd evc(N,N);
//    MatrixXcd inv(N,N);
//    for(int i=0;i<N;i++)
//    {
//        for(int j=0;j<N;j++)
//        {
//            evc(i,j) = es.eigenvectors()(i,j);
//        }
//    }
//    inv = evc.inverse();
//
//
//    for(int i=0;i<N;i++)
//    {
//        for(int j=0;j<N;j++)
//        {
//            Teig[i][j] = evc(i,j);
//            Tinv[i][j] = inv(i,j);
//        }
//    }
//}

void LegendrePolynomials::setCoefficients(int maxOrder)
{
   Np = maxOrder+1;
   cp = new double*[Np];
   for(int k=0;k<=Np;k++)  cp[k] = new double[k+1];
   if(maxOrder>=0)
   {   cp[ 0][ 0]=+1.0000000000e+00;
   }
   if(maxOrder>=1)
   {   cp[ 1][ 0]=+0.0000000000e+00;  cp[ 1][ 1]=+1.0000000000e+00;
   }
   if(maxOrder>=2)
   {   cp[ 2][ 0]=-5.0000000000e-01;  cp[ 2][ 1]=+0.0000000000e+00;  cp[ 2][ 2]=+1.5000000000e+00;
   }

   if(maxOrder> 2)
   {   // higher orders
       for(int k=2;k<=Np;k++)
       {
           double a = (2.*k-1.)/k;
           double b = (1.*k-1.)/k;
           for(int i=0;i<k  ;i++) cp[k][i+1] += a*cp[k-1][i];
           for(int i=0;i<k-1;i++) cp[k][i  ] -= b*cp[k-2][i] ;
       }
   }
}

double LegendrePolynomials::evalPolynomial(int order, double z)
{
   double p1 = cp[order][order];
   for(int n=order-1;n>=0;n--) p1 = z*p1+cp[order][n];
   return p1;
}

GaussQuadrature::GaussQuadrature(int order)
{
    quadNodes(order, -1.,1.);
}
GaussQuadrature::GaussQuadrature(int order, double a, double b)
{
    quadNodes(order,a,b);
}
GaussQuadrature::~GaussQuadrature()
{
    if(qcoo) delete[] qcoo;
    if(qwgt) delete[] qwgt;
}

void GaussQuadrature::quadNodes(int order, double a, double b)
{
    Nq = 1; for(int i=0;i<order;i++) Nq *= 2;
    qcoo = new double[Nq];
    qwgt = new double[Nq];
    double eps = 1.e-15;

    if(order==0)
    {
        qcoo[ 0] = 0.000000000000000;  qwgt[ 0] = 2.000000000000000;
    }
    else
    if(order==1)
    {
        qcoo[ 0] =-0.577350269189626;  qwgt[ 0] = 1.000000000000000;
        qcoo[ 1] = 0.577350269189626;  qwgt[ 1] = 1.000000000000000;
    }
    else
    if(order==2)
    {
        qcoo[ 0] =-0.339981043584856;  qwgt[ 0] = 0.652145154862546;
        qcoo[ 1] = 0.339981043584856;  qwgt[ 1] = 0.652145154862546;
        qcoo[ 2] =-0.861136311594053;  qwgt[ 2] = 0.347854845137454;
        qcoo[ 3] = 0.861136311594053;  qwgt[ 3] = 0.347854845137454;
    }
    else
    if(order==3)
    {
        qcoo[ 0] =-0.183434642495650;  qwgt[ 0] = 0.362683783378362;
        qcoo[ 1] = 0.183434642495650;  qwgt[ 1] = 0.362683783378362;
        qcoo[ 2] =-0.525532409916329;  qwgt[ 2] = 0.313706645877887;
        qcoo[ 3] = 0.525532409916329;  qwgt[ 3] = 0.313706645877887;
        qcoo[ 4] =-0.796666477413627;  qwgt[ 4] = 0.222381034453374;
        qcoo[ 5] = 0.796666477413627;  qwgt[ 5] = 0.222381034453374;
        qcoo[ 6] =-0.960289856497536;  qwgt[ 6] = 0.101228536290376;
        qcoo[ 7] = 0.960289856497536;  qwgt[ 7] = 0.101228536290376;
    }
    else
    if(order==4)
    {
        qcoo[ 0] =-0.09501250983764;  qwgt[ 0] = 0.18945061045507;
        qcoo[ 1] = 0.09501250983764;  qwgt[ 1] = 0.18945061045507;
        qcoo[ 2] =-0.28160355077926;  qwgt[ 2] = 0.18260341504492;
        qcoo[ 3] = 0.28160355077926;  qwgt[ 3] = 0.18260341504492;
        qcoo[ 4] =-0.45801677765723;  qwgt[ 4] = 0.16915651939500;
        qcoo[ 5] = 0.45801677765723;  qwgt[ 5] = 0.16915651939500;
        qcoo[ 6] =-0.61787624440264;  qwgt[ 6] = 0.14959598881658;
        qcoo[ 7] = 0.61787624440264;  qwgt[ 7] = 0.14959598881658;
        qcoo[ 8] =-0.75540440835500;  qwgt[ 8] = 0.12462897125553;
        qcoo[ 9] = 0.75540440835500;  qwgt[ 9] = 0.12462897125553;
        qcoo[10] =-0.86563120238783;  qwgt[10] = 0.09515851168249;
        qcoo[11] = 0.86563120238783;  qwgt[11] = 0.09515851168249;
        qcoo[12] =-0.94457502307323;  qwgt[12] = 0.06225352393865;
        qcoo[13] = 0.94457502307323;  qwgt[13] = 0.06225352393865;
        qcoo[14] =-0.98940093499165;  qwgt[14] = 0.02715245941175;
        qcoo[15] = 0.98940093499165;  qwgt[15] = 0.02715245941175;
    }
    else
    if(order==5)
    {
        qcoo[ 0] =-0.04830766568774;  qwgt[ 0] = 0.09654008851473;
        qcoo[ 1] = 0.04830766568774;  qwgt[ 1] = 0.09654008851473;
        qcoo[ 2] =-0.14447196158280;  qwgt[ 2] = 0.09563872007927;
        qcoo[ 3] = 0.14447196158280;  qwgt[ 3] = 0.09563872007927;
        qcoo[ 4] =-0.23928736225214;  qwgt[ 4] = 0.09384439908080;
        qcoo[ 5] = 0.23928736225214;  qwgt[ 5] = 0.09384439908080;
        qcoo[ 6] =-0.33186860228213;  qwgt[ 6] = 0.09117387869576;
        qcoo[ 7] = 0.33186860228213;  qwgt[ 7] = 0.09117387869576;
        qcoo[ 8] =-0.42135127613064;  qwgt[ 8] = 0.08765209300440;
        qcoo[ 9] = 0.42135127613064;  qwgt[ 9] = 0.08765209300440;
        qcoo[10] =-0.50689990893223;  qwgt[10] = 0.08331192422695;
        qcoo[11] = 0.50689990893223;  qwgt[11] = 0.08331192422695;
        qcoo[12] =-0.58771575724076;  qwgt[12] = 0.07819389578707;
        qcoo[13] = 0.58771575724076;  qwgt[13] = 0.07819389578707;
        qcoo[14] =-0.66304426693022;  qwgt[14] = 0.07234579410885;
        qcoo[15] = 0.66304426693022;  qwgt[15] = 0.07234579410885;
        qcoo[16] =-0.73218211874029;  qwgt[16] = 0.06582222277636;
        qcoo[17] = 0.73218211874029;  qwgt[17] = 0.06582222277636;
        qcoo[18] =-0.79448379596794;  qwgt[18] = 0.05868409347854;
        qcoo[19] = 0.79448379596794;  qwgt[19] = 0.05868409347854;
        qcoo[20] =-0.84936761373257;  qwgt[20] = 0.05099805926238;
        qcoo[21] = 0.84936761373257;  qwgt[21] = 0.05099805926238;
        qcoo[22] =-0.89632115576605;  qwgt[22] = 0.04283589802223;
        qcoo[23] = 0.89632115576605;  qwgt[23] = 0.04283589802223;
        qcoo[24] =-0.93490607593774;  qwgt[24] = 0.03427386291302;
        qcoo[25] = 0.93490607593774;  qwgt[25] = 0.03427386291302;
        qcoo[26] =-0.96476225558751;  qwgt[26] = 0.02539206530926;
        qcoo[27] = 0.96476225558751;  qwgt[27] = 0.02539206530926;
        qcoo[28] =-0.98561151154527;  qwgt[28] = 0.01627439473091;
        qcoo[29] = 0.98561151154527;  qwgt[29] = 0.01627439473091;
        qcoo[30] =-0.99726386184948;  qwgt[30] = 0.00701861000947;
        qcoo[31] = 0.99726386184948;  qwgt[31] = 0.00701861000947;
    }
    else
    if(order==6)
    {
        qcoo[ 0] =-0.02435029266342;  qwgt[ 0] = 0.04869095700914;
        qcoo[ 1] = 0.02435029266342;  qwgt[ 1] = 0.04869095700914;
        qcoo[ 2] =-0.07299312178780;  qwgt[ 2] = 0.04857546744150;
        qcoo[ 3] = 0.07299312178780;  qwgt[ 3] = 0.04857546744150;
        qcoo[ 4] =-0.12146281929612;  qwgt[ 4] = 0.04834476223480;
        qcoo[ 5] = 0.12146281929612;  qwgt[ 5] = 0.04834476223480;
        qcoo[ 6] =-0.16964442042399;  qwgt[ 6] = 0.04799938859646;
        qcoo[ 7] = 0.16964442042399;  qwgt[ 7] = 0.04799938859646;
        qcoo[ 8] =-0.21742364374001;  qwgt[ 8] = 0.04754016571483;
        qcoo[ 9] = 0.21742364374001;  qwgt[ 9] = 0.04754016571483;
        qcoo[10] =-0.26468716220877;  qwgt[10] = 0.04696818281621;
        qcoo[11] = 0.26468716220877;  qwgt[11] = 0.04696818281621;
        qcoo[12] =-0.31132287199021;  qwgt[12] = 0.04628479658131;
        qcoo[13] = 0.31132287199021;  qwgt[13] = 0.04628479658131;
        qcoo[14] =-0.35722015833767;  qwgt[14] = 0.04549162792742;
        qcoo[15] = 0.35722015833767;  qwgt[15] = 0.04549162792742;
        qcoo[16] =-0.40227015796399;  qwgt[16] = 0.04459055816376;
        qcoo[17] = 0.40227015796399;  qwgt[17] = 0.04459055816376;
        qcoo[18] =-0.44636601725346;  qwgt[18] = 0.04358372452932;
        qcoo[19] = 0.44636601725346;  qwgt[19] = 0.04358372452932;
        qcoo[20] =-0.48940314570705;  qwgt[20] = 0.04247351512365;
        qcoo[21] = 0.48940314570705;  qwgt[21] = 0.04247351512365;
        qcoo[22] =-0.53127946401989;  qwgt[22] = 0.04126256324262;
        qcoo[23] = 0.53127946401989;  qwgt[23] = 0.04126256324262;
        qcoo[24] =-0.57189564620263;  qwgt[24] = 0.03995374113272;
        qcoo[25] = 0.57189564620263;  qwgt[25] = 0.03995374113272;
        qcoo[26] =-0.61115535517239;  qwgt[26] = 0.03855015317862;
        qcoo[27] = 0.61115535517239;  qwgt[27] = 0.03855015317862;
        qcoo[28] =-0.64896547125466;  qwgt[28] = 0.03705512854024;
        qcoo[29] = 0.64896547125466;  qwgt[29] = 0.03705512854024;
        qcoo[30] =-0.68523631305423;  qwgt[30] = 0.03547221325688;
        qcoo[31] = 0.68523631305423;  qwgt[31] = 0.03547221325688;
        qcoo[32] =-0.71988185017161;  qwgt[32] = 0.03380516183714;
        qcoo[33] = 0.71988185017161;  qwgt[33] = 0.03380516183714;
        qcoo[34] =-0.75281990726053;  qwgt[34] = 0.03205792835485;
        qcoo[35] = 0.75281990726053;  qwgt[35] = 0.03205792835485;
        qcoo[36] =-0.78397235894334;  qwgt[36] = 0.03023465707240;
        qcoo[37] = 0.78397235894334;  qwgt[37] = 0.03023465707240;
        qcoo[38] =-0.81326531512280;  qwgt[38] = 0.02833967261426;
        qcoo[39] = 0.81326531512280;  qwgt[39] = 0.02833967261426;
        qcoo[40] =-0.84062929625258;  qwgt[40] = 0.02637746971505;
        qcoo[41] = 0.84062929625258;  qwgt[41] = 0.02637746971505;
        qcoo[42] =-0.86599939815409;  qwgt[42] = 0.02435270256871;
        qcoo[43] = 0.86599939815409;  qwgt[43] = 0.02435270256871;
        qcoo[44] =-0.88931544599511;  qwgt[44] = 0.02227017380838;
        qcoo[45] = 0.88931544599511;  qwgt[45] = 0.02227017380838;
        qcoo[46] =-0.91052213707850;  qwgt[46] = 0.02013482315353;
        qcoo[47] = 0.91052213707850;  qwgt[47] = 0.02013482315353;
        qcoo[48] =-0.92956917213194;  qwgt[48] = 0.01795171577570;
        qcoo[49] = 0.92956917213194;  qwgt[49] = 0.01795171577570;
        qcoo[50] =-0.94641137485840;  qwgt[50] = 0.01572603047602;
        qcoo[51] = 0.94641137485840;  qwgt[51] = 0.01572603047602;
        qcoo[52] =-0.96100879965205;  qwgt[52] = 0.01346304789672;
        qcoo[53] = 0.96100879965205;  qwgt[53] = 0.01346304789672;
        qcoo[54] =-0.97332682778991;  qwgt[54] = 0.01116813946013;
        qcoo[55] = 0.97332682778991;  qwgt[55] = 0.01116813946013;
        qcoo[56] =-0.98333625388463;  qwgt[56] = 0.00884675982636;
        qcoo[57] = 0.98333625388463;  qwgt[57] = 0.00884675982636;
        qcoo[58] =-0.99101337147674;  qwgt[58] = 0.00650445796898;
        qcoo[59] = 0.99101337147674;  qwgt[59] = 0.00650445796898;
        qcoo[60] =-0.99634011677196;  qwgt[60] = 0.00414703326056;
        qcoo[61] = 0.99634011677196;  qwgt[61] = 0.00414703326056;
        qcoo[62] =-0.99930504173577;  qwgt[62] = 0.00178328072170;
        qcoo[63] = 0.99930504173577;  qwgt[63] = 0.00178328072170;
    }

    else
    {   // general case
        double* u = new double[Nq+1];
        int m = (Nq+1)/2;
        int ia = 0;
        for(int k=m-1;k>=0;k--)
        {
            double x = cos(M_PI*(k+0.75)/(Nq+0.5));
            double r = 0.;
            double w = 0.;
            double dpdx=0.;
            double err = 1.e+20;
            while(err > eps)
            {   u[0] = 1.;
                u[1] = x;
                for(int j=1;j<Nq;j++) u[j+1] = ( (2.*j+1.)*x*u[j] - j*u[j-1])/(j+1.);
                dpdx = Nq*(x*u[Nq]-u[Nq-1])/(x*x-1.);
                r = x;
                x = r-u[Nq]/dpdx;
                err = abs(x-r);
            }
            dpdx = Nq*(x*u[Nq]-u[Nq-1])/(x*x-1.);
            w = 2./((1.-x*x)*dpdx*dpdx);
            qcoo[ia] = - x;  qwgt[ia] = w;  ia++;
            qcoo[ia] = + x;  qwgt[ia] = w;  ia++;
        }
        delete[] u;
    }
    // sort in ascending order
    double* xi = new double[Nq];
    double* wg = new double[Nq];
    int j = 0;
    for(int k=Nq-2;k>=0;k-=2,j++)
    {
        xi[j] = qcoo[k];
        wg[j] = qwgt[k];
    }
    for(int k=1;k<Nq;k+=2,j++)
    {
        xi[j] = qcoo[k];
        wg[j] = qwgt[k];
    }

    // rescale nodes and weights
    double xc = (b+a)/2.;
    double dx = (b-a)/2.;
    for(int k=0;k<Nq;k++)
    {
        qcoo[k] = xc + xi[k] * dx;
        qwgt[k] =      wg[k] * dx;
    }
    delete[] xi;
    delete[] wg;
}

