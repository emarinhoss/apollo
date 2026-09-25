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

/** return the total number of elements */
    unsigned totalElemNum(){
        return _Klocal;
    }

/** return the number of equation in the system */
    unsigned equationsNumber(){
        return _meqn;
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

/**
 * Is local cell K a real cell of the mesh, rather than a PETSc ghost cell?
 *
 * ApSolver::createMesh calls DMPlexConstructGhostCells (solvers/apsolver.cc),
 * which appends one ghost cell per boundary facet to the height-0 stratum. On
 * the shipped rmf_frc mesh that turns 7792 cells into 7984, and on the antenna
 * example's vacuumDisc.msh 3929 into 4086 - in both cases exactly the number of
 * boundary lines in the .msh. Ghost cells have no geometry of their own, so
 * anything that walks the stratum asking cells for normals or areas has to skip
 * them, or it reads zero-length edges.
 *
 * The nodal DG scheme tolerates them because its per-cell work is harmless on
 * one; the slope limiter did not, and stopped with "Edge length of 0 found in
 * element 7792". Note that a cell being a ghost is not the same as its data
 * being meaningless - ghost cells carry the exterior state and must still be
 * written when a pass writes the whole vector.
 *
 * Implemented by cone size rather than by DMPlexGetGhostCellStratum so that it
 * holds for any cell the DM reports without three vertices, whatever produced
 * it.
 */
    bool isRealCell(unsigned K) const {
        return K < _isRealCell.size() && _isRealCell[K] != 0;
    }

/** Return Face to Face connnectivity */
    void ElementTOElementANDFace(unsigned K, int *ftf){
        for(unsigned ff=0; ff<6; ff++)
            ftf[ff] = _ETETF[K][ff];
    }

/** Return Average operator for the limiter */
    void LimiterElementAVE(REAL *limAVG, REAL *dropAVG){
        for(unsigned k=0; k<_NpE; k++)
            limAVG[k] = 0.0;

        for(unsigned k=0; k<_NpE; k++)
            for(unsigned m=0; m<_NpE; m++)
                limAVG[k] += 0.5*_Mass[k+m*_NpE];

        int sk = 0;
        for(unsigned k=0; k<_NpE; k++)
            for(unsigned m=0; m<_NpE; m++)
                dropAVG[sk++] = -limAVG[m];

        for(unsigned m=0; m<_NpE; m++)
            dropAVG[m*_NpE+m] += 1.;
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

/** Calculate the area of a given element */
    REAL elementArea(int k);

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
    REAL *_Mass;
    Mat _IVand;
    int *_Fmask;
    int **_EtoV; // element to verticies connectivity matrix
    int **_ETETF; // stores the element and face connections (Element K's face number F is connected to
                  // elmement K2's face number F2
    int _Ktotal; // total number of element in the entire domain
    int _Klocal, _kLocalInt; // number of elements in this processor
    std::vector<char> _isRealCell; // per local cell, see isRealCell()
    int _Vlocal; // number of nodes in this processor
    REAL **_xcoord; // node x-coordinates
    REAL **_ycoord; // node y-coordinates
    REAL _dtscale; // minimun value of element_radius/(element_perimeter/2)
    REAL _rmin; // minimun spacing between nodes

///** Label information to identify boundary faces */
//    std::vector<std::string> _bcLabels;

};

#endif // WXNODALDGGEOMETRY2D_H
