#ifndef WXNODALDGGEOMETRY2D_H
#define WXNODALDGGEOMETRY2D_H


// WarpX lib includes
#include <wxindexer.h>
#include <petscdmplex.h>

// std includes
#include <vector>

template <typename REAL>
class wxNodalDGgeometry2D
{
  public:
/**
 * Create weights, absciae and Legendre polynomials of given spatial
 * order
 *
 * @param meqn Number of equations
 * @param spatialOrder Spatial order
 */
    wxNodalDGgeometry2D(DM dm, unsigned meqn, unsigned Spor);

/** Destroctor */
    virtual ~wxNodalDGgeometry2D();

/** Matrix Inverse */
    void invertMatrix(Mat A, Mat *invA);

/** return the inverse-Vandermonde Matrix */
    Mat inverseVandermonde(){
        return _IVand;
    }

/** return the number of nodes per element */
    unsigned NpElem(){
        return _NpE;
    }

/** return the number of nodes per element */
    unsigned NpFaces(){
        return _NpF;
    }

/** return the number of nodes per element */
    unsigned NfElem(){
        return _NfE;
    }

/** return Fmask, the nodes numbers along each face */
    void returnFmask(int *fmask){
        for(unsigned k=0; k<_NfE*_NpF; k++)
            fmask[k] = _Fmask[k];
    }

/** return the x-coordinate for element K node N */
    REAL Xcoordinate(unsigned K, unsigned N){
        return _xcoord[K][N];
    }

/** return the x-coordinate for element K node N */
    REAL Ycoordinate(unsigned K, unsigned N){
        return _ycoord[K][N];
    }

/** Return Face to Face connnectivity */
    void ElementTOElementANDFace(unsigned K, int *ftf){
        for(unsigned ff=0; ff<6; ff++)
            ftf[ff] = _ETETF[K][ff];
    }

///** calculate the  geometric factors for element k*/
    void GeometricFactors2d(int k, REAL geom[]);
//    void GeomFacs2d(int k, REAL *rx, REAL *sx, REAL *ry, REAL *sy, REAL *J);

///** calculate the face normals for a given element k*/
    void Normals2d(int k, REAL norms[]);
//    void FaceNodesNormals2d(int k, REAL *nx, REAL *ny, REAL *Fscale);

///** LIFT the flux: calculate the flux through the element boundaries */
//    void LIFT_flux(int k, REAL *nflux, REAL *nFrhs);

///** Calculate the weak derivatives */
//    void weakDericatives(unsigned K, REAL *DxnDy, REAL *Fflux, REAL *Gflux);

/** return minimun radius of inscribed circle */
    REAL dtscale2D(){
        return _dtscale;
    }
/** return minimun radius of inscribed circle */
    REAL rMin(){
        return _rmin;
   }

/** multi vector by the inverse mass matrix */
    void multiplyBYinverseMassMatrix(REAL* input, REAL* output);

///** return the evaluate the filter matrix */
//    void calculateFilter(REAL *filter, REAL *filterMatrix);

  private:
/**
 * Find which faces in one element connect to what
 * face in the neighboring element.
 */
    void FacePair2d(DM dm);

/**
 * Calculate the node coordinates for each element.
 */
    void CalculateNodeCoordinates2d(DM dm);

/** Transfer all the Matrices into row major arrays */
    void petscMatTOArray(Mat A, REAL *array);

/** No of equations */
    unsigned _meqn;
/** Polynomial order */
    unsigned _polyOr;
/** Data management */
    DM _dm;

    int _NpE; // number of points/nodes per element
    int _NpF; // number of points/nodes per face
    int _NfE;  // number of faces per element
    int _totNFace;  // total number of interior faces
    REAL *_r, *_s;  // (r,s) coordinates of reference nodes
    REAL *_Ds, *_Dr, *_Vand, *_VVT; // element matrices, VVT is the inverse of the mass matrix
    Mat _IVand;
    int *_Fmask;
    int **_EtoV; // element to verticies connectivity matrix
    int **_ETETF; // stores the element and face connections (Element K's face number F is connected to
                  // elmement K2's face number F2
    int _Klocal; // number of elements in this processor
    int _Vlocal; // number of nodes in this processor
    REAL **_xcoord; // node x-coordinates
    REAL **_ycoord; // node y-coordinates
    REAL _dtscale; // minimun value of element_radius/(element_perimeter/2)
    REAL _rmin; // minimun spacing between nodes

///** Label information to identify boundary faces */
//    std::vector<std::string> _bcLabels;

};

#endif // WXNODALDGGEOMETRY2D_H
