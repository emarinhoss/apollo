/*
 * A stand-in for <petsc.h> that declares exactly what one PETSc release
 * declares, for the handful of symbols src/lib/petsc_compat.h shims.
 *
 * WHY THIS EXISTS. The shims are guarded by version, and a guard that names
 * the wrong version is invisible on the machine you happen to build on: Apollo
 * built fine against PETSc 3.19 for months while failing to compile at all
 * against 3.15, which is what Ubuntu 22.04 (and Google Colab) ships. Nothing
 * caught it because CI only ever ran the newer release.
 *
 * WHERE THE TABLE BELOW COMES FROM. Each boundary was read off upstream's own
 * headers rather than remembered, by fetching include/petscdmplex.h at tagged
 * releases from gitlab.com/petsc/petsc and grepping for the symbol:
 *
 *   DMPlexGetHybridBounds          declared through v3.12.5, absent in v3.13.6
 *   DMPlexSetAdjacencyUseCone      declared through v3.13.6, absent in v3.14.6
 *   DMPlexSetAdjacencyUseClosure   declared through v3.13.6, absent in v3.14.6
 *   DMPlexGetLabelValue            already absent in v3.9.4
 *
 * DMPlexGetHybridBounds is not in include/petsc/private/dmpleximpl.h at any of
 * those releases either, so there is no private declaration to collide with.
 *
 * This file is the "what PETSc has" half of the pair; petsc_compat.h is the
 * "when to substitute" half. test_petsc_compat.py compiles them together and
 * fails if they disagree in either direction - a shim that does not fire when
 * the symbol is gone, or one that fires when the symbol is still there.
 */
#ifndef APOLLO_TEST_PETSC_STUB_H
#define APOLLO_TEST_PETSC_STUB_H

#ifndef STUB_PETSC_MINOR
#error "compile with -DSTUB_PETSC_MINOR=<n> to select the PETSc release to emulate"
#endif

#define PETSC_VERSION_RELEASE  1
#define PETSC_VERSION_MAJOR    3
#define PETSC_VERSION_MINOR    STUB_PETSC_MINOR
#define PETSC_VERSION_SUBMINOR 0

/* Spelled the way petscversion.h spells it, so the guards are exercised as
 * written rather than against a simplified stand-in. */
#define PETSC_VERSION_LT(MAJOR, MINOR, SUBMINOR)                          \
  (PETSC_VERSION_RELEASE == 1 &&                                          \
   (PETSC_VERSION_MAJOR < (MAJOR) ||                                      \
    (PETSC_VERSION_MAJOR == (MAJOR) &&                                    \
     (PETSC_VERSION_MINOR < (MINOR) ||                                    \
      (PETSC_VERSION_MINOR == (MINOR) &&                                  \
       PETSC_VERSION_SUBMINOR < (SUBMINOR))))))
#define PETSC_VERSION_GE(MAJOR, MINOR, SUBMINOR) \
  (!PETSC_VERSION_LT(MAJOR, MINOR, SUBMINOR))

typedef int PetscErrorCode;
typedef int PetscInt;
typedef enum { PETSC_FALSE = 0, PETSC_TRUE = 1 } PetscBool;
typedef struct _p_DM *DM;

#define PETSC_NULLPTR nullptr

#define PetscFunctionBegin      do { } while (0)
#define PetscFunctionReturn(x)  return (x)
#define CHKERRQ(ierr)           do { if (ierr) return (ierr); } while (0)

/* Present at every release in range. */
extern PetscErrorCode DMPlexGetHeightStratum(DM, PetscInt, PetscInt *, PetscInt *);
extern PetscErrorCode DMGetLabelValue(DM, const char[], PetscInt, PetscInt *);

/* The legacy spellings, declared only where upstream still declares them.
 * Signatures copied verbatim from v3.12.5/v3.13.6 include/petscdmplex.h, so a
 * shim whose signature drifted is a compile error too. */
#if PETSC_VERSION_MINOR <= 12
extern PetscErrorCode DMPlexGetHybridBounds(DM, PetscInt *, PetscInt *, PetscInt *, PetscInt *);
#endif

#if PETSC_VERSION_MINOR <= 13
extern PetscErrorCode DMPlexSetAdjacencyUseCone(DM, PetscBool);
extern PetscErrorCode DMPlexSetAdjacencyUseClosure(DM, PetscBool);
#endif

#endif /* APOLLO_TEST_PETSC_STUB_H */
