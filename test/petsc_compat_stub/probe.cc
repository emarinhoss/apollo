/*
 * Compiled once per emulated PETSc release by test_petsc_compat.py.
 *
 * Binding each name to a function pointer of the expected type is what makes
 * this a test rather than a smoke check: an undeclared name fails, a shim that
 * fired when the real declaration is still in scope fails as a static-after-
 * extern redeclaration, and a shim whose parameter list drifted fails to
 * convert. Nothing is linked - the extern declarations have no definitions.
 */
#include "petsc_compat.h"

void apollo_petsc_compat_probe(void)
{
  PetscErrorCode (*hybrid)(DM, PetscInt *, PetscInt *, PetscInt *, PetscInt *) =
      &DMPlexGetHybridBounds;
  PetscErrorCode (*use_cone)(DM, PetscBool)    = &DMPlexSetAdjacencyUseCone;
  PetscErrorCode (*use_closure)(DM, PetscBool) = &DMPlexSetAdjacencyUseClosure;
  PetscErrorCode (*label)(DM, const char[], PetscInt, PetscInt *) =
      &DMPlexGetLabelValue;

  (void)hybrid;
  (void)use_cone;
  (void)use_closure;
  (void)label;
}
