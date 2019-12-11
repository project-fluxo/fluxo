#ifndef OPTIM_H
#define OPTIM_H

#include <p8est_extended.h>
#include <p4est_to_p8est.h>

///////////////////////////#include <p4est_to_p8est.h>////////////////////////////////////////////////////
//
// optimize mortars: search for mortars being fully MPI_MINE and add them to innerMortars
//
//
///////////////////////////////////////////////////////////////////////////////
static void
MortarOptimisation(p4est_iter_face_info_t *info, void *user_data);

#endif /* OPTIM_H */