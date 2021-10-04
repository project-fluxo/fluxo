//!==================================================================================================================================
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
//!==================================================================================================================================
#include "p4fluxo.h"
#include "optimisation.h"

static void
MortarOptimisation(p4est_iter_face_info_t *info, void *user_data) {
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;
    p4est_inner_data_t *ghost_data = (p4est_inner_data_t *) user_data;
    p4est_aux_data_t *aux_data = (p4est_aux_data_t *) p4est->user_pointer;
    int myrank = p4est->mpirank;
    int *AddToMortarInnerSide = (int *) aux_data->tmp[0]; //could be changed
    int *MortarSideNumeration = (int *) aux_data->tmp[1]; //could be changed
    int FirstMPIMINESide = *(int *) aux_data->tmp[2];     //couldn't be changed
    int LastMPIMINESide = *(int *) aux_data->tmp[3];      //couldn'tmake clean be changed
    int ghost = 0;
    int ghostor = 0;

    p4est_quadrant_t *quad;
    p4est_inner_data_t *dataquad;

    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);

    if (sides->elem_count == 1) {
        // There is nothing to do here with bounadry;

        return;
    }
    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);

    if (side[0]->is_hanging || side[1]->is_hanging) {
        //One Side is Mortar
        iBigSide = side[0]->is_hanging == 0 ? 0 : 1;
        iSmallSide = side[0]->is_hanging != 0 ? 0 : 1;
        int Bigface = side[iBigSide]->face;
        int SmallFace = side[iSmallSide]->face;
        if (side[iBigSide]->is.full.is_ghost == 0) //Big side is not MPI
        {
            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0) {
                    ghostor++; // not ghost side, but local
                    continue;
                } else {
                    int ghostid = side[iSmallSide]->is.hanging.quadid[j];
                    int SidesID = ghost_data[ghostid].SidesID[SmallFace];

                    if ((SidesID > FirstMPIMINESide) && (SidesID <= LastMPIMINESide))
                        ghost++;
                    // MPI MINE Side;
                }
            }
            if ((ghostor + ghost) > 4) exit(1);

            if ((ghost != 0) && ((ghostor + ghost) == 4)) {   // THat's It!
                (*AddToMortarInnerSide)++;
                quad = side[iBigSide]->is.full.quad;
                dataquad = (p4est_inner_data_t *) quad->p.user_data;
                int face = side[iBigSide]->face;
                (*MortarSideNumeration)++;
                dataquad->SidesID[face] = (*MortarSideNumeration);
                dataquad->flips[face] = 0;
                if (dataquad->SidesID[face] <= 0) {
                    printf("******************* ERROR MORTAR!!! \n");
                    fflush(stdout);
                    exit(1);
                }
            }

        } else { //We have only 4 small sides
            // Nothing to do
            return;

        }
    } else //Then it is a one side
    { //Nothing to do here. We optimize Mortar side
        return;
    }
}
