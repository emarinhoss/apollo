/*
 * PETSc Compatibility Header
 *
 * This file provides compatibility shims for PETSc API changes between
 * version 3.6 (original Apollo development) and 3.19+ (current Homebrew).
 */

#ifndef PETSC_COMPAT_H
#define PETSC_COMPAT_H

#include <petsc.h>

/* Silence deprecation warnings for old API */
#define PETSC_SILENCE_DEPRECATION_WARNINGS_3_8_0
#define PETSC_SILENCE_DEPRECATION_WARNINGS_3_9_0
#define PETSC_SILENCE_DEPRECATION_WARNINGS_3_19_0

/* PETSC_NULL was deprecated in 3.19, replaced with PETSC_NULLPTR */
#ifndef PETSC_NULL
#define PETSC_NULL PETSC_NULLPTR
#endif

/*
 * DMPlexGetHybridBounds was removed in PETSc 3.13, not 3.18.
 *
 * The guard here said 3.18 until a build on Ubuntu 22.04 failed: that release
 * ships PETSc 3.15, where the symbol is gone but the shim did not fire, so
 * every file calling it failed to compile. Verified against upstream headers -
 * include/petscdmplex.h declares it through v3.12.5 and no longer declares it
 * in v3.13.6; it is not in include/petsc/private/dmpleximpl.h either. So the
 * shim is needed for 3.13 and later, and must NOT be defined at 3.12 and
 * earlier, where it would clash with the real declaration.
 *
 * test/test_petsc_compat.py pins both directions.
 */
#if PETSC_VERSION_GE(3, 13, 0)
/*
 * Callers use the returned bound directly as a loop limit, so every output is
 * written before anything that can fail: an early CHKERRQ return must not leave
 * the caller looping to an indeterminate value. This is what -Wmaybe-uninitialized
 * was reporting at five call sites.
 *
 * NOTE: fMax, eMax and vMax are returned as -1, meaning "not tracked". A loop
 * written as `for (f = fStart; f < fMax; ++f)` against one of those will run zero
 * times rather than over the mesh. Only the cMax form is supported here; ask for
 * the others and you are on your own.
 */
static inline PetscErrorCode DMPlexGetHybridBounds(DM dm, PetscInt *cMax, PetscInt *fMax, PetscInt *eMax, PetscInt *vMax)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (cMax) *cMax = 0;
  if (fMax) *fMax = -1; /* Not tracking face hybrid bounds */
  if (eMax) *eMax = -1; /* Not tracking edge hybrid bounds */
  if (vMax) *vMax = -1; /* Not tracking vertex hybrid bounds */

  if (cMax) {
    /* Get the end of interior cells */
    PetscInt cStart, cEnd;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
    *cMax = cEnd;
  }

  PetscFunctionReturn(0);
}
#endif

/*
 * DMPlexSetAdjacencyUseCone and DMPlexSetAdjacencyUseClosure were removed in
 * PETSc 3.14, not 3.18 - the same mis-dating as the hybrid-bounds shim above,
 * and the same symptom on Ubuntu 22.04's PETSc 3.15. Verified against upstream
 * headers: both are declared in include/petscdmplex.h through v3.13.6 and are
 * absent from v3.14.6 onwards. Replaced upstream by DMSetBasicAdjacency.
 */
#if PETSC_VERSION_GE(3, 14, 0)
static inline PetscErrorCode DMPlexSetAdjacencyUseCone(DM dm, PetscBool useCone)
{
  /* In newer PETSc, adjacency is controlled differently.
   * These functions are no-ops for compatibility. */
  (void)dm;
  (void)useCone;
  return 0;
}

static inline PetscErrorCode DMPlexSetAdjacencyUseClosure(DM dm, PetscBool useClosure)
{
  /* In newer PETSc, adjacency is controlled differently.
   * These functions are no-ops for compatibility. */
  (void)dm;
  (void)useClosure;
  return 0;
}
#endif

/* DMPlexGetLabelValue was moved to DMGetLabelValue in PETSc 3.11 */
#if PETSC_VERSION_GE(3, 11, 0)
#ifndef DMPlexGetLabelValue
#define DMPlexGetLabelValue DMGetLabelValue
#endif
#endif

#endif /* PETSC_COMPAT_H */
