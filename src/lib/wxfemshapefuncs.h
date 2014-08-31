#ifndef __wxfemshapefuncs__h__
#define __wxfemshapefuncs__h__

#include <cmath>

/**
   Values of shape function in
*/
template<typename REAL>
void
shapes1D(int nodes, int QuadPs, REAL *eta, REAL **phi, REAL **dphi)
{
    switch (nodes)
    {
    case 2:
        // for linear elements two nodes are needed
        for(unsigned k=0; k<QuadPs; ++k){
            phi[0][k] = 0.5*(1.-eta[k]);
            phi[1][k] = 0.5*(1.+eta[k]);

            dphi[0][k]= -0.5;
            dphi[1][k]= 0.5;
        }
        break;

    case 3:
        // for quadratic elements three nodes are needed
        for(unsigned k=0; k<QuadPs; ++k){
            phi[0][k] = -0.5*eta[k]*(1.-eta[k]);
            phi[1][k] = 1.-eta[k]*eta[k];
            phi[2][k] = 0.5*eta[k]*(1.+eta[k]);

            dphi[0][k]= -0.5*(1.-2.*eta[k]);
            dphi[1][k]= -2.*eta[k];
            dphi[2][k]= 0.5*(1.+2.*eta[k]);
        }
        break;
    case 4:
        // for cubic elements four nodes are needed
        for(unsigned k=0; k<QuadPs; ++k){
            phi[0][k] = -9./16.*(1.-eta[k])*(1./9.-eta[k]*eta[k]);
            phi[1][k] = 27./16.*(1.-eta[k]*eta[k])*(1./3.-eta[k]);
            phi[2][k] = 27./16.*(1.-eta[k]*eta[k])*(1./3.+eta[k]);
            phi[3][k] = -9./16.*(1.+eta[k])*(1./9.-eta[k]*eta[k]);

            dphi[0][k] = -9./16.*(-1.*(1./9.-eta[k]*eta[k])+(1-eta[k])*(-2.*eta[k]));
            dphi[1][k] = 27./16.*((-2.*eta[k])*(1./3.-eta[k])-(1.-eta[k]*eta[k]));
            dphi[2][k] = 27./16.*((-2.*eta[k])*(1./3.+eta[k])+(1.-eta[k]*eta[k]));
            dphi[3][k] = -9./16.*((1./9.-eta[k]*eta[k])+(1.+eta[k])*(-2.*eta[k]));
        }
        break;
    default:
        // cell centered value only
        phi[0][0] = 1.;

        dphi[0][0]= 0.;
    }
}

/**
   Values and weights of the Gauss-Lobatto quadrature
*/
template<typename REAL>
void
gaulob(int n, REAL *x, REAL *w)
{
    switch (n)
    {
    case 4:
        // for cubic elements four nodes are needed
        for(unsigned k=0; k<n; ++k){
            x[0] = -1.;
            x[1] = -sqrt(1./5.);
            x[2] = sqrt(1./5.);
            x[3] = 1.;

            w[0] = 1./6.;
            w[1] = 5./6.;
            w[2] = 5./6.;
            w[3] = 1./6.;
        }
        break;
    case 5:
        // for cubic elements four nodes are needed
        for(unsigned k=0; k<n; ++k){
            x[0] = -1.;
            x[1] = -sqrt(3./7.);
            x[2] = 0.0;
            x[3] = sqrt(3./7.);
            x[3] = 1.;

            w[0] = 1./10.;
            w[1] = 49./90.;
            w[2] = 32./45.;
            w[3] = 49./90.;
            w[3] = 1./10.;
        }
        break;
    default:
        // cell centered value only
        x[0] = -1.;
        x[1] = 0.;
        x[2] = 1.;

        w[0]= 1./3.;
        w[1]= 4./3.;
        w[2]= 1./3.;
    }
}


#endif // wxfemshapefuncs__h__
