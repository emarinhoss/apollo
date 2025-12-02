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
#define PETSC_SILENCE_DEPRECATION_WARNINGS_3_9_0
#define PETSC_SILENCE_DEPRECATION_WARNINGS_3_19_0

/* PETSC_NULL was deprecated in 3.19, replaced with PETSC_NULLPTR */
#ifndef PETSC_NULL
#define PETSC_NULL PETSC_NULLPTR
#endif

/* DMPlexGetHybridBounds was removed in PETSc 3.18 */
#if PETSC_VERSION_GE(3, 18, 0)
static inline PetscErrorCode DMPlexGetHybridBounds(DM dm, PetscInt *cMax, PetscInt *fMax, PetscInt *eMax, PetscInt *vMax)
{
  PetscInt pStart, pEnd;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetChart(dm, &pStart, &pEnd); CHKERRQ(ierr);

  if (cMax) {
    /* Get the end of interior cells */
    PetscInt cStart, cEnd;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
    *cMax = cEnd;
  }
  if (fMax) *fMax = -1; /* Not tracking face hybrid bounds */
  if (eMax) *eMax = -1; /* Not tracking edge hybrid bounds */
  if (vMax) *vMax = -1; /* Not tracking vertex hybrid bounds */

  PetscFunctionReturn(0);
}
#endif

/* DMPlexGetLabelValue was moved to DMGetLabelValue in PETSc 3.11 */
#if PETSC_VERSION_GE(3, 11, 0)
#ifndef DMPlexGetLabelValue
#define DMPlexGetLabelValue DMGetLabelValue
#endif
#endif

#endif /* PETSC_COMPAT_H */
