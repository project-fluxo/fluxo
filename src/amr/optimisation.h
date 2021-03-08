//!=================================================================================================================================
//! Copyright (c) 2018 - 2020 Alexander Astanin
//!
//! This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
//! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
//! of the License, or (at your option) any later version.
//!
//! FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
//!
//! You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
//!=================================================================================================================================
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
