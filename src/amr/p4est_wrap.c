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
#include "optimisation.h"
#include "connectivity.h"
#include "p4fluxo.h"
// #include "ready.fc"
#include "optimisation.c"
#include "p8est_communication.h"
//#define NON_OPTIMIZED //Optimized Variant not ready yet

//#ifndef NON_OPTIMIZED

//Begin of ready.fc

void p4est_destroy_f(p4est_t **p4est) {
    p4est_destroy(*p4est);
}

void p4est_connectivity_destroy_f(p4est_connectivity_t **conn) {
    p4est_connectivity_destroy(*conn);
}


static void
ElementCounterSetOldToZero_iter(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    p4est_inner_data_t *dataquad = (p4est_inner_data_t *) q->p.user_data;
    // p4est_fortran_data_t *data = (p4est_fortran_data_t *)user_data;
    int i;
    //ElementID must be set in previous function (for example getData444)

    dataquad->weight = 0;
    for (i = 1; i < 8; i++) {
        dataquad->OldElementID[i] = 0;
    }
    dataquad->OldElementID[0] = dataquad->ElementID;

#ifndef NON_OPTIMIZED
    dataquad->IsChanged = 0;
    for (i = 0; i < 6; i++) {
        dataquad->OldSidesID[i] = dataquad->SidesID[i];

    }
#endif
    return;
}

static void
ElementCounter_iter(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    p4est_inner_data_t *dataquad = (p4est_inner_data_t *) q->p.user_data;
    p4est_fortran_data_t *data = (p4est_fortran_data_t *) user_data;
    int i;
    //ElementID must be set in previous function (for example getData444)
    dataquad->ElementID = ++data->nElems;
    dataquad->weight = 0;
#ifndef NON_OPTIMIZED
    dataquad->IsChanged = 0;
#endif
    for (i = 0; i < 6; i++) {
#ifndef NON_OPTIMIZED
        dataquad->OldSidesID[i] = dataquad->SidesID[i];
#endif
        dataquad->SidesID[i] = 0;

    }

    return;
}


static void
ElementCounterNew_iter(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    p4est_inner_data_t *dataquad = (p4est_inner_data_t *) q->p.user_data;
    int *nElems = (int *) user_data;
    int i;
    dataquad->ElementID = ++(*nElems);
    dataquad->weight = 0;
    for (i = 0; i < 8; i++) {
        dataquad->OldElementID[i] = 0;
    }
    dataquad->OldElementID[0] = dataquad->ElementID;
    for (i = 0; i < 6; i++) {
        dataquad->SidesID[i] = -1;
#ifndef NON_OPTIMIZED
        dataquad->OldSidesID[i] = -1;
#endif
        dataquad->flips[i] = -1;
        // dataquad->SidesRatio[i] = 0;
    }
    // printf("dataquad->ElementID = %d \n", dataquad->ElementID);
    // fflush(stdout);

    // data->nElems++;
    return;
}


int p4est_init_f(int mpicomm1) {
    // int mpiret = 0;

    /* Initialize MPI; see sc_mpi.h.
   * If configure --enable-mpi is given these are true MPI calls.
   * Else these are dummy functions that simulate a single-processor run. */
    // mpiret = sc_MPI_Init(NULL, NULL);
    // SC_CHECK_MPI(mpiret);
    // mpicomm = sc_MPI_COMM_WORLD;
#if MPI 
    mpicomm = MPI_Comm_f2c(mpicomm1); //Protect from different bits for Frotran and C    
#else  /*MPI*/
    mpicomm = sc_MPI_COMM_WORLD;
#endif  /*MPI*/
    // printf("mpicomm1 = %d, size = %d\n", mpicomm1, sizeof(mpicomm1));
    //  printf("mpicomm!!! %d \n",mpicomm);
    
    // printf("mpicomm = %d, size = %d \n", mpicomm, sizeof(mpicomm));
    // fflush(stdout);
    /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
    // sc_init(mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL); //SC_LP_SILENT , SC_LP_ERROR sc.h. */
    sc_init(mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL); //SC_LP_SILENT , SC_LP_ERROR
    // sc_init(MPI_Comm_f2c(mpicomm1), 1, 1, NULL, SC_LP_ESSENTIAL); //SC_LP_SILENT , C_LP_ERROR

    p4est_init(NULL, SC_LP_ERROR); //SC_LP_SILENT , SC_LP_ERROR
   
    return 0;
}


void p4est_finalize() {
    sc_finalize();

}

/** Broadcast a connectivity structure that exists only on one process to all.
!  *  On the other processors, it will be allocated using p8est_connectivity_new.
!  *  \param [in] conn_in For the root process the connectivity to be broadcast,
!  *                      for the other processes it must be NULL.
!  *  \param [in] root    The rank of the process that provides the connectivity.
!  *  \param [in] comm    The MPI communicator.
!  *  \return             For the root process this is a pointer to \a conn_in.
!  *                      Else, a pointer to a newly allocated connectivity
!  *                      structure with the same values as \a conn_in on the
!  *                      root process.
!  */
p8est_connectivity_t *p8est_conn_bcast(p8est_connectivity_t *conn_in,
                                       int root, /* The owner of CONN*/
                                       int comm) {

    int mpimyrank = 0;
#if MPI 
    MPI_Comm_rank(MPI_Comm_f2c(comm), &mpimyrank);
#endif  /*MPI*/

    if (mpimyrank != root) {
        p4est_connectivity_destroy(conn_in);
        conn_in = NULL;
    }

    fflush(stdout);
    sc_MPI_Comm mpicomm;
#if MPI 
    mpicomm = MPI_Comm_f2c(comm); //Protect from different bits for Frotran and C    
#else  /*MPI*/
    mpicomm = sc_MPI_COMM_WORLD;
#endif  /*MPI*/
    return p8est_connectivity_bcast(conn_in, root, mpicomm);
    // return NULL;
};

p4est_mpi_data_t *p4estGetMPIData_f(p4est_t *p4est) {
    p4est_mpi_data_t *data;
    data = (p4est_mpi_data_t *) malloc(sizeof(p4est_mpi_data_t));
    data->mpisize = p4est->mpisize;
    data->mpirank = p4est->mpirank;
    data->local_num_quad = p4est->local_num_quadrants;
    data->global_num_quad = p4est->global_num_quadrants;
    data->offsetMPI = p4est->global_first_quadrant;
    return data;
}

void p4estDelMPIData_f(p4est_mpi_data_t *data) {
    free(data);
    data = NULL;
}

void save_p4est(p4est_t *p4est, char in[]) {
    // printf("save p4est!! \n%s!\n", in);//p4est->local_num_quadrants);
    p4est_save_ext(in, p4est, 0, 0);
    return;
};



p4est_t *load_p4est(int mpicomm1, char *in) {

    p8est_connectivity_t *conn;
    p4est_t *p4est;
    sc_MPI_Comm mpicomm;
    // int myrank;
#if MPI 
    mpicomm = MPI_Comm_f2c(mpicomm1);
#else  /*MPI*/
    mpicomm = sc_MPI_COMM_WORLD;      
#endif  /*MPI*/
    // MPI_Comm_rank(mpicomm, &myrank);
    p4est = p8est_load_ext(in,      //const char *filename,
                           mpicomm, //sc_MPI_Comm mpicomm,
                           0,       //size_t data_size,
                           0,       //int load_data,
                           1,       //int autopartition,
                           0,       //int broadcasthead,
                           NULL,    //void *user_pointer,
                           &conn);   //p8est_connectivity_t * *connectivity);
    
    p4est->connectivity = conn;
    p4est_reset_data(p4est, sizeof(p4est_inner_data_t), NULL, NULL);
    int nEl = 0;
    p4est_iterate(p4est,                     /* the forest */
                  NULL,                   /* the ghost layer May be LAter!!! */
                  (void *) &nEl,        /* the synchronized ghost data */
                  ElementCounterNew_iter, /* callback to compute each quad's
                                             interior contribution to du/dt */
                  NULL,                   /* SidesCount_iter,            /* callback to compute each quads'
                                             faces' contributions to du/du */
                  NULL,                   /* there is no callback for the
                                             edges between quadrants */
                  NULL);
    // printf("Connectivity!! %d \n", p4est->connectivity->num_trees);
    // exit(0);
    return p4est;
};

p8est_connectivity_t *GetConnectivity(p4est_t *p4est) {
    // printf("p4est->connectivity = %p \n", p4est->connectivity);
    return p4est->connectivity;
};
p4est_connectivity_t *p8est_connectivity_new_periodic_f() {
    return p8est_connectivity_new_periodic();
}

// end of ready.fc

//#endif
// int KKK = 0; // Trigger for debugging. Delete after debugging 
//  * \return 1 if \a q should be refined,   0 otherwise.           * /
static int
refine_fn(p4est_t *p4est, p4est_topidx_t which_tree,
          p4est_quadrant_t *q) {
    int *ToRefine = (int *) p4est->user_pointer;

    p4est_inner_data_t *dataquad = (p4est_inner_data_t *) q->p.user_data;
    /*! if (ToRefine[dataquad->ElementID - 1] > 1)
    ! {
         exit(1);
     }*/
    // if (ToRefine[dataquad->ElementID - 1] > 0)
    // printf("ToRefine[dataquad->ElementID - 1] %d \n", ToRefine[dataquad->ElementID - 1]);
    // printf(" q->level %d \n",  q->level);
    if ((ToRefine[dataquad->ElementID - 1] > 0 && ToRefine[dataquad->ElementID - 1] > q->level)
        || (-ToRefine[dataquad->ElementID - 1] - 1 > q->level)) {
        // printf("REFINE!!!!!!!!!!!!!!! \n");
        return 1;
    } else {
        return 0;
    }
}

////////////////////////////////////////////////////////////
// \return 1 if \a children should be coarsened, 0 otherwise.       
// We use also the minimum level number with -1 (-1 means coarse to the root)
//
//                     
static int
coarse_fn(p4est_t *p4est,
          p4est_topidx_t which_tree,
          p4est_quadrant_t *children[]) {

    int *ToCoarse = (int *) p4est->user_pointer;
    p4est_inner_data_t *data;
    int Coarse8 = 0;
    int i;

    /* Check every 8 quads */
    for (i = 0; i < P8EST_CHILDREN; i++) {

        data = (p4est_inner_data_t *) children[i]->p.user_data;

        if (data->OldElementID[0] < 0 || data->ElementID == 0) //This Element was refined and we don't need to check it.
        {
            return 0;
        }
        if (ToCoarse[data->ElementID - 1] < 0) {
            if (children[i]->level >
                (-ToCoarse[data->ElementID - 1] - 1)) { Coarse8++; } // The level compared with -ToCoarsep[] -1
        } else {
            return 0;
        }
    }

    if (Coarse8 == 8) {
        return 1;
    } else {
        return 0;
    }
}

static void
replace_quads(p4est_t *p4est, p4est_topidx_t which_tree,
              int num_outgoing,
              p4est_quadrant_t *outgoing[],
              int num_incoming, p4est_quadrant_t *incoming[]) {

    int i, j;


    p4est_inner_data_t *parentquaddata, *childquaddata; // = (p4est_inner_data_t *)q->p.user_data;
    // return;

    if (num_outgoing > 1) {
        /* this is coarsening */
        // The quads are already initialized
        // This is a new quad
        parentquaddata = (p4est_inner_data_t *) incoming[0]->p.user_data;


        for (i = 0; i < P8EST_CHILDREN; i++) {
            childquaddata = (p4est_inner_data_t *) outgoing[i]->p.user_data;
            // parentquaddata->ElementID=0;
            parentquaddata->OldElementID[i] = childquaddata->OldElementID[0];

        }
    } else {

        /* this is refinement */
        parentquaddata = (p4est_inner_data_t *) outgoing[0]->p.user_data;

        if (parentquaddata->OldElementID[1] > 0)
            //This quad was Coarsed but Refined again by balance 2:1
        {
            // printf("refine Again\n");
            for (i = 0; i < P8EST_CHILDREN; i++) {
                childquaddata = (p4est_inner_data_t *) incoming[i]->p.user_data;
                childquaddata->OldElementID[0] = parentquaddata->OldElementID[i];

                int i1;
                for (i1 = 1; i1 < P8EST_CHILDREN; i1++) {
                    childquaddata->OldElementID[i1] = 0;
                }


#ifndef NON_OPTIMIZED
                // for (i1 = 0; i1 < 6; i1++)
                // // Safe for redefinition of Element. Insist on Metrics recalculation
                // {
                //     childquaddata->OldSidesID[i1]=-1;
                // }
#endif
            }
        } else {
            for (i = 0; i < P8EST_CHILDREN; i++) {
                childquaddata = (p4est_inner_data_t *) incoming[i]->p.user_data;
                childquaddata->OldElementID[0] = -parentquaddata->OldElementID[0];
                int j;
                for (j = 1; j < P8EST_CHILDREN; j++) {
                    childquaddata->OldElementID[j] = 0;
                }
                childquaddata->ElementID = 0;
            }
        }
    }
}

static void
SetSides_iter(p4est_iter_face_info_t *info, void *user_data) {
    int i, j, iBigSide = 0, iSmallSide = 0;

    p4est_t *p4est = info->p4est;
    p4est_SetSide_data_t *data = (p4est_SetSide_data_t *) user_data;
    p4est_quadrant_t *quad;
    p4est_inner_data_t *dataquad, *dataquad0, *dataquad1;
    int orientation = info->orientation;
    int GhostHere = 0;
    int which_face;
    int direction, minus = 0;

    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);


    if (sides->elem_count == 1) //Or sides->elem_num == 1
    {
        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        direction = side[0]->face / 2; /* 0 == x, 1 == y, 2 == z */
        quad = side[0]->is.full.quad;
        which_face = side[0]->face;
        dataquad = (p4est_inner_data_t *) side[0]->is.full.quad->p.user_data;
        dataquad->SidesID[which_face] = data->CurrentBCSide++;
        dataquad->flips[which_face] = 0;
        dataquad->weight++;
        return;
    }

    // side[i]->is.hanging.quad[j]->p.user_data;
    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    direction = side[0]->face / 2; /* 0 == x, 1 == y, 2 == z */
    minus = side[0]->face % 2;
    if (side[0]->is_hanging || side[1]->is_hanging) { //TODO: Complete Mortar Faces
        // ThereisMortar++;
        iBigSide = side[0]->is_hanging == 0 ? 0 : 1;
        iSmallSide = side[0]->is_hanging != 0 ? 0 : 1;
        int BigFace = side[iBigSide]->face;
        int SmallFace = side[iSmallSide]->face;
        if (side[iBigSide]->is.full.is_ghost == 0) {

            int flip = GetHFlip(BigFace, SmallFace, orientation);

            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0) {
                    int CurrentSide = data->CurrentInnerSide++;
                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = (p4est_inner_data_t *) quad->p.user_data;
                    dataquad->SidesID[SmallFace] = CurrentSide;
                    dataquad->weight++;
                    if (dataquad->SidesID[SmallFace] <= 0) {
                        printf("******************* ERROR!!! \n");
                        fflush(stdout);
                        exit(1);
                    }
                    dataquad->flips[SmallFace] = flip;
                    //ghost_data[ghostid] -> SideID[direction];
                }
                    //
                else {
                    GhostHere++;
                    //Ghost Small Side
                    // The ID has already set in ghost_data
                }
            }
            // Set the ID for BigSide when there is no any ghost side
            if (!GhostHere) { //SidesID must = 0
                int CurrentSide = data->CurrentMortrarInnerSide++;
                quad = side[iBigSide]->is.full.quad;
                dataquad = (p4est_inner_data_t *) quad->p.user_data;
                dataquad->SidesID[BigFace] = CurrentSide;
                dataquad->flips[BigFace] = 0;
                dataquad->weight += MORTAR_SIDE_WEIGHT;

                if (dataquad->SidesID[BigFace] <= 0) {
                    printf("******************* ERROR!!! \n");
                    fflush(stdout);
                    exit(1);
                }
            } else //One side is Ghost and not MPI MINE , so it is a
            {

                quad = side[iBigSide]->is.full.quad;
                dataquad = (p4est_inner_data_t *) quad->p.user_data;
                if (dataquad->SidesID[BigFace] <= 0) //The sideID not set yet
                {
                    int CurrentSide = data->CurrentMPIMortarSide++;
                    dataquad->SidesID[BigFace] = CurrentSide;
                    dataquad->flips[BigFace] = 0;
                    //dataquad->weight += MORTAR_SIDE_WEIGHT;
                }
                if (dataquad->SidesID[BigFace] <= 0) {
                    printf("******************* ERROR!!! \n");
                    fflush(stdout);
                    exit(1);
                }
            }
            // Then Number of the side has already set up
        } else {
            // We have Big Side MPI, so set up small MPI Sides Only
            return;

            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0) { // Set up when MPI MINE and MPI YOUR
                }
                    //
                else {
                    //Ghost Small Side
                    // Actually The Side is not on this Proc
                }
            }


        }

    } else { // No hanging sides
        if (side[0]->is.full.is_ghost || side[1]->is.full.is_ghost) { //Nothing to do HERe
        } else {
            int CurrentSide = data->CurrentInnerSide++;

            dataquad0 = (p4est_inner_data_t *) side[0]->is.full.quad->p.user_data;
            dataquad1 = (p4est_inner_data_t *) side[1]->is.full.quad->p.user_data;
            dataquad0->weight++;
            dataquad1->weight++;
            int opp_face0 = side[0]->face;
            int opp_face1 = side[1]->face;
            int flip = GetHFlip(opp_face0, opp_face1, orientation);
            if (minus) {
                dataquad0->flips[opp_face0] = 0; // I set it to 0 as a Master Side
                dataquad1->flips[opp_face1] = flip;
            } else {
                dataquad0->flips[opp_face0] = flip;
                dataquad1->flips[opp_face1] = 0; // I set it to 0 as a Master Side
            }
            dataquad0->SidesID[opp_face0] = CurrentSide;
            dataquad1->SidesID[opp_face1] = CurrentSide;

            if (dataquad0->SidesID[opp_face0] <= 0) {
                printf("******************* ERROR!!! \n");
                fflush(stdout);
                exit(1);
            }
        }
        return;
    }
}

//////////////////////////////////////////////////////////
// This is a simple function to use as quad iterator.
// It was used for debugging purpose
// 
// 
///////////////////////////////////////////
static void
Elm_iter(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    p4est_inner_data_t *dataquad = (p4est_inner_data_t *) q->p.user_data;
    p4est_t *p4est = info->p4est;
    // p4est_fortran_data_t *data = (p4est_fortran_data_t *)user_data;
    int i, j = 0;
    for (i = 0; i < 8; i++) {
        //         {
        if (dataquad->ElementID == 27 || dataquad->OldElementID[i] == 27) {
            printf("i = %d \n", i);
            printf("!!!!!!!! refine \n dataquad->OldElementID[0] = %d\n", dataquad->OldElementID[0]);
            printf("!!!!!!!!111 refine \n dataquad->OldElementID[%d] = %d\n", i, dataquad->OldElementID[i]);
            printf("!!!!!!1 refine \n dataquad->ElementID = %d\n", dataquad->ElementID);
            // exit(1);
        };
    }
    // if (dataquad->ElementID == 4 && p4est->mpirank == 1 && KKK>0)
    if (dataquad->OldElementID[0] == 0) {
        printf("Elm_iter Error in Numeration!!! ElemntId = %d \n", dataquad->ElementID);
        exit(1);
    };

}


static void
SetElementToSide_iter(p4est_iter_volume_info_t *info, void *user_data) {
    int i, j;
    p4est_quadrant_t *q = info->quad;
    p4est_t *p4est = info->p4est;
    p4est_inner_data_t *dataquad = (p4est_inner_data_t *) q->p.user_data;
    p4est_fortran_data_t *data = (p4est_fortran_data_t *) user_data;
    int *EtS = data->EtSPtr;
    int *StE = data->StEPtr;
    int ElemID = dataquad->ElementID;
    for (i = 0; i < 6; i++) {
        //ElemID begins with 1, but i and j with 0
        j = 0;
        int HIndex = P2H_side[i];
        EtS[(ElemID - 1) * 6 * 2 + (HIndex - 1) * 2 + (j)] = dataquad->SidesID[i];
        j = 1;
        EtS[(ElemID - 1) * 6 * 2 + (HIndex - 1) * 2 + (j)] = dataquad->flips[i];
    }
    for (i = 0; i < 6; i++) {
        int iSide = dataquad->SidesID[i];
        int Flip = dataquad->flips[i];
        if (Flip == 0) { //Master Side
            j = 1;
            StE[(iSide - 1) * 5 + (j - 1)] = dataquad->ElementID;
            j = 3;
            StE[(iSide - 1) * 5 + (j - 1)] = P2H_side[i];
        } else //Slave Side
        {

            j = 1;
            // fflush(stdout);
            j = 2;
            StE[(iSide - 1) * 5 + (j - 1)] = dataquad->ElementID;
            // i = 3;
            // SideToElement[(iSide - 1) * 5 + (i - 1)] = P2H_side[i];
            j = 4;
            StE[(iSide - 1) * 5 + (j - 1)] = P2H_side[i];
            j = 5;
            StE[(iSide - 1) * 5 + (j - 1)] = dataquad->flips[i];
            // }
        }
        // TODO!
        // IF(aSide%sideID .LE. nBCSides) BC(aSide%sideID)=aSide%BCIndex
    }
    return;
}

static void
SetSideToElement_iter(p4est_iter_face_info_t *info, void *user_data) {
    int i, j;
    p4est_t *p4est = info->p4est;
    p8est_connectivity_t *conn = p4est->connectivity;
    p4est_fortran_data_t *data = (p4est_fortran_data_t *) user_data;
    p4est_inner_data_t *ghost_data = (p4est_inner_data_t *) p4est->user_pointer;
    int *SideToElement = data->StEPtr;
    p4est_quadrant_t *quad;
    int orientation = info->orientation;
    int *MortarInfo = data->MIPtr;
    int *MortarType = data->MTPtr;

    p4est_inner_data_t *dataquad, *dataquad0, *dataquad1;
    int iBigSide = 0; //Number of Big side 0 or 1 in (mortar) in
    int iSmallSides = 0;
    int which_face, iSide;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);

    int AlreadyCount = 0; // This Flag set to 1 if there is a simple side and we just one number for this side

    if (sides->elem_count == 1) //Or sides->elem_num == 1
    {

        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        int treeid = side[0]->treeid;
        int face = side[0]->face;
        quad = side[0]->is.full.quad;
        dataquad = (p4est_inner_data_t *) side[0]->is.full.quad->p.user_data;
        int SideID = dataquad->SidesID[face];

        data->BCs[SideID - 1] = ((int32_t *) conn->tree_to_attr)[treeid * 6 + face];

        return;
    }
    // printf("Next Side!!\n");
    // side[i]->is.hanging.quad[j]->p.user_data;
    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    if (side[0]->is_hanging || side[1]->is_hanging) { //TODO: Complete Mortar Faces

        iBigSide = side[0]->is_hanging ? 1 : 0;
        iSmallSides = side[1]->is_hanging ? 1 : 0;
        int BigFace = side[iBigSide]->face;
        int SmallFace = side[iSmallSides]->face;
        p4est_inner_data_t *Bigdataquad;
        if (side[iBigSide]->is.full.is_ghost == 0) {
            quad = side[iBigSide]->is.full.quad;
            dataquad = (p4est_inner_data_t *) quad->p.user_data;
            int BigSideID = dataquad->SidesID[BigFace];
            // firstBCSide          = 1
            // firstMortarInnerSide = firstBCSide         +nBCSides
            //             firstInnerSide       = firstMortarInnerSide+nMortarInnerSides
            // firstMPISide_MINE    = firstInnerSide      +nInnerSides
            // firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE
            // firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR

            // lastBCSide           = firstMortarInnerSide-1
            // lastMortarInnerSide  = firstInnerSide    -1
            int ADD = 0;
            int firstMortarInnerSide = 1 + data->nBCSides; //1
            int firstInnerSide = firstMortarInnerSide + data->nMortarInnerSides; //1 + 0 = 1
            int firstMPISide = firstInnerSide + data->nInnerSides; // 1 + 8 = 9 
            int firstMortarMPISide = firstMPISide + data->nMPISides; // 9 + 14 = 21
            int lastMortarInnerSide = firstInnerSide - 1; // = 0
            if (BigSideID <= lastMortarInnerSide) { //Inner Mortar
                ADD = firstMortarInnerSide;
            } else {   //MPI Mortar
                ADD = firstMortarMPISide - data->nMortarInnerSides; // 21 - 0 = 21
            }


            int SideID = BigSideID + 1 - ADD; //  1 + 1 - 21 = -19 Index for MPI Sides


            i = 1;
            MortarType[(BigSideID - 1) * 2 + (i - 1)] = 1; //MortarType
            i = 2;
            MortarType[(BigSideID - 1) * 2 + (i - 1)] = SideID;

            p4est_inner_data_t *Smalldataquads[P8EST_HALF];
            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                if (side[iSmallSides]->is.hanging.is_ghost[j] == 0) {
                    // Smalldataquads[j] = (p4est_inner_data_t *)side[iSmallSides]->is.hanging.quad[j]->p.user_data;

                    quad = side[iSmallSides]->is.hanging.quad[j];
                    dataquad = (p4est_inner_data_t *) quad->p.user_data;
                    int iSide = dataquad->SidesID[SmallFace];
                    int flip = dataquad->flips[SmallFace];
                    // iSide = Smalldataquads[j]->SidesID[SmallFace];
                    // int GetHMortar(int PMortar, int PSide, int PnbSide, int PFlip)
                    int jIndex = 0; //f(j); //Convert index from P4 to H
                    int Pside = BigFace;
                    int PnbSide = SmallFace;
                    int PFlip = orientation;
                    // if (KKK==0){
                    jIndex = GetHMortar(j, Pside, PnbSide, PFlip);
                    i = 1;
                    MortarInfo[(SideID - 1) * 4 * 2 + (jIndex - 1) * 2 + (i - 1)] = iSide;
                    i = 2;
                    MortarInfo[(SideID - 1) * 4 * 2 + (jIndex - 1) * 2 + (i - 1)] = 0 * flip;
                    // }
                } else {
                    /// Use GHOST DATA
                    int ghostid = side[iSmallSides]->is.hanging.quadid[j];
                    int iSide = ghost_data[ghostid].SidesID[SmallFace];
                    int flip = ghost_data[ghostid].flips[SmallFace];
                    // iSide = Smalldataquads[j]->SidesID[SmallFace];
                    // int GetHMortar(int PMortar, int PSide, int PnbSide, int PFlip)
                    int jIndex = 0; //f(j); //Convert index from P4 to H
                    int Pside = BigFace;
                    int PnbSide = SmallFace;
                    int PFlip = orientation;
                    jIndex = GetHMortar(j, Pside, PnbSide, PFlip);

                    jIndex = GetHMortar(j, Pside, PnbSide, PFlip);
                    i = 1;
                    MortarInfo[(SideID - 1) * 4 * 2 + (jIndex - 1) * 2 + (i - 1)] = iSide;
                    i = 2;
                    MortarInfo[(SideID - 1) * 4 * 2 + (jIndex - 1) * 2 + (i - 1)] = flip;//Not Zero

                }
            }

        } else {

        }

        return;
    } else {
        return;
    }
}

static void
SidesCounter_iter(p4est_iter_face_info_t *info, void *user_data) {
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;
    p4est_fortran_data_t *data = (p4est_fortran_data_t *) user_data;
    p4est_quadrant_t *quad1, *quad2, *quad;
    p4est_inner_data_t *dataquad, *dataquad1, *dataquad2;
    int which_face;
    int ghost = 0;
    p4est_iter_face_side_t *side[2] = {NULL, NULL};
    sc_array_t *sides = &(info->sides);
    // int AlreadyCount = 0; // This Flag set to 1 if there is a simple side and we just one number for this side

    if (sides->elem_count == 1) {
        //printf("BCSide!!! info->tree_boundary = %d \n", info->tree_boundary );
        //  printf("BCSide!!! sides->elem_count = %d \n", sides->elem_count );

        data->nSides++;
        data->nBCSides++;
        return;
    }
    // side[i]->is.hanging.quad[j]->p.user_data;
    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);

    if (side[0]->is_hanging || side[1]->is_hanging) {
        //One Side is Mortar
        iBigSide = side[0]->is_hanging == 0 ? 0 : 1;
        iSmallSide = side[0]->is_hanging != 0 ? 0 : 1;

        if (side[iBigSide]->is.full.is_ghost == 0) //Big side is not MPI
        {

            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j]) {
                    ghost++;
                }
            }
            if (ghost == 0) //There is no small mpi sides. Add Mortar inner side
            {
                data->nMortarInnerSides++;
            } else //One or more small mpi sides. Add Big side as MPIMortar
            {
                data->nMortarMPISides++;
            }
            data->nMPISides += ghost;       //Increase number of small MPISides
            data->nSides += 4;              //Total nmber of sides
            data->nSides++;                 //Add also a MortarSide to nSides
            data->nInnerSides += 4 - ghost; //Number of inner sides
        } else { //We use only 4 small sides and don't take into account the big side

            for (j = 0; j < P8EST_HALF; j++) {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j]) {
                    ghost++;
                }
            }

            data->nSides += (4 - ghost);    //Small sides are MPI sides. Ghost is for another proc.
            data->nMPISides += (4 - ghost); //Small sides  are MPI Sides
        }
    } else //Then it is a one side
    {

        if (side[0]->is.full.is_ghost || side[1]->is.full.is_ghost) {
            //Debug
            int SideIn = 0;
            int SideGhost = 1;
            if (side[0]->is.full.is_ghost) {
                SideIn = 1;
                SideGhost = 0;
            }
            quad = side[SideIn]->is.full.quad;
            dataquad = (p4est_inner_data_t *) quad->p.user_data;
            //Debug
            data->nSides++;
            data->nMPISides++;

        } else //Not Ghost
        {
            // quad = side[i]->is.full.quad;
            data->nSides++;
            data->nInnerSides++;
        }
    }
}


static void
SidesCounter2_iter(p4est_iter_face_info_t *info, void *user_data) {
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;
    p4est_SetMPISide_data_t *data = (p4est_SetMPISide_data_t *) user_data;
    p4est_quadrant_t *quad;
    p4est_inner_data_t *dataquad;
    int which_face;
    int ghost = 0;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);
    int AlreadyCount = 0; // This Flag set to 1 if there is a simple side and we just one number for this side

    if (sides->elem_count == 1) {

        return;
    }
    // side[i]->is.hanging.quad[j]->p.user_data;
    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);

    if (side[0]->is_hanging || side[1]->is_hanging) {
        //One Side is Mortar
        iBigSide = side[0]->is_hanging == 0 ? 0 : 1;
        iSmallSide = side[0]->is_hanging != 0 ? 0 : 1;

        if (side[iBigSide]->is.full.is_ghost == 0) //Big side is not MPI
        {

            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j]) {
                    // ghost++;
                    int ghostid = side[iSmallSide]->is.hanging.quadid[j];
                    int proc = data->ghost_to_proc[ghostid];
                    data->nMPISidesCount[proc]++;
                }
            }

        } else { //We use only 4 small sides and don't take into account the big side

            for (j = 0; j < P8EST_HALF; j++) {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0) {
                    // ghost++;
                    int ghostid = side[iBigSide]->is.full.quadid;
                    int proc = data->ghost_to_proc[ghostid];
                    data->nMPISidesCount[proc]++;
                }
            }

        }
    } else //Then it is a one side
    {

        if (side[0]->is.full.is_ghost || side[1]->is.full.is_ghost) {
            int SideIn = 0;
            int SideGhost = 1;
            if (side[0]->is.full.is_ghost) {
                SideIn = 1;
                SideGhost = 0;
            }
            int ghostid = side[SideGhost]->is.full.quadid;
            int proc = data->ghost_to_proc[ghostid];
            data->nMPISidesCount[proc]++;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// This function is used to sed SideID from 1 to nMPISides_Proc(iNbProc) for every 
// NB Processor. Then it is used with shifting to setup the correct number
//
//
///////////////////////////////////////////////////////////////////////////////
static void
SetMPISidesAUXNumber(p4est_iter_face_info_t *info, void *user_data) {
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;
    p4est_inner_data_t *ghost_data = (p4est_inner_data_t *) user_data;
    p4est_aux_data_t *aux_data = (p4est_aux_data_t *) p4est->user_pointer;
    int myrank = p4est->mpirank;
    int *SideID = ((int *) (aux_data->tmp[0]));
    int *NbProc = (int *) aux_data->tmp[1];
    int *offsetMPISides_MINE = (int *) aux_data->tmp[2];
    int *offsetMPISides_YOUR = (int *) aux_data->tmp[3];
    int *nMPISides_Proc = (int *) aux_data->tmp[4];
    int *ghost_to_proc = (int *) aux_data->tmp[5];
    int iNBProc = *((int *) (aux_data->tmp[6]));
    int *nMPISides_MINE_Proc = aux_data->tmp[7];
    int *nMPISides_YOUR_Proc = aux_data->tmp[8];
    int direction, minus = 0;
    int orientation = info->orientation;
    p4est_quadrant_t *quad;
    p4est_inner_data_t *dataquad;

    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);

    /// <=== Here
    if (sides->elem_count == 1) {

        return;
    }

    // side[i]->is.hanging.quad[j]->p.user_data;
    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    direction = side[0]->face / 2; /* 0 == x, 1 == y, 2 == z */
    minus = side[0]->face % 2;

    if (side[0]->is_hanging || side[1]->is_hanging) {
        //One Side is Mortar
        iBigSide = side[0]->is_hanging == 0 ? 0 : 1;
        iSmallSide = side[0]->is_hanging != 0 ? 0 : 1;
        int Bigface = side[iBigSide]->face;
        int SmallFace = side[iSmallSide]->face;
        int flip = GetHFlip(Bigface, SmallFace, orientation);
        if (side[iBigSide]->is.full.is_ghost == 0) //Big side is not MPI
        {
            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j]) {
                    // Here Set UP
                    int ghostid = side[iSmallSide]->is.hanging.quadid[j];
                    int NBProc1 = ghost_to_proc[ghostid];
                    int NBProc2 = NbProc[iNBProc];
                    if (NBProc1 != NBProc2) continue;
                    // I need direction from
                    ghost_data[ghostid].SidesID[SmallFace] = (*SideID)++;

                    if (myrank < NBProc1) {
                        if (ghost_data[ghostid].SidesID[SmallFace] <= nMPISides_MINE_Proc[iNBProc]) {//MINE
                            ghost_data[ghostid].SidesID[SmallFace] += offsetMPISides_MINE[iNBProc];
                            ghost_data[ghostid].flips[SmallFace] = 0;
                        } else {//YOUR
                            ghost_data[ghostid].SidesID[SmallFace] =
                                    (ghost_data[ghostid].SidesID[SmallFace] - nMPISides_MINE_Proc[iNBProc]) +
                                    offsetMPISides_YOUR[iNBProc];
                            ghost_data[ghostid].flips[SmallFace] = flip;
                        }
                    } else // if (myrank < NBProc1)
                    {
                        if (ghost_data[ghostid].SidesID[SmallFace] <= nMPISides_YOUR_Proc[iNBProc]) {   //MINE
                            ghost_data[ghostid].SidesID[SmallFace] += offsetMPISides_YOUR[iNBProc];
                            ghost_data[ghostid].flips[SmallFace] = flip; //According to Fortran
                        } else {   //YOUR
                            ghost_data[ghostid].SidesID[SmallFace] =
                                    (ghost_data[ghostid].SidesID[SmallFace] - nMPISides_YOUR_Proc[iNBProc]) +
                                    offsetMPISides_MINE[iNBProc];
                            ghost_data[ghostid].flips[SmallFace] = 0;
                        }
                    }

                }
            }

        } else         //if (side[iBigSide]->is.full.is_ghost == 0)
        { //We use only 4 small sides and don't take into account the big side

            for (j = 0; j < P8EST_HALF; j++) {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0) {
                    // Here Set UP
                    int ghostid = side[iBigSide]->is.full.quadid;
                    int NBProc1 = ghost_to_proc[ghostid];
                    int NBProc2 = NbProc[iNBProc];
                    if (NBProc1 != NBProc2)
                        continue;
                    // I need direction from
                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = (p4est_inner_data_t *) quad->p.user_data;
                    dataquad->SidesID[SmallFace] = (*SideID)++;

                    if (myrank < NBProc1) {
                        if (dataquad->SidesID[SmallFace] <= nMPISides_MINE_Proc[iNBProc]) { //MINE
                            dataquad->SidesID[SmallFace] += offsetMPISides_MINE[iNBProc];
                            dataquad->flips[SmallFace] = 0;
                            if (dataquad->SidesID[SmallFace] <= 0) {
                                printf("******************* ERROR!!! \n");
                                fflush(stdout);
                                exit(1);
                            }
                        } else { //YOUR
                            dataquad->SidesID[SmallFace] =
                                    (dataquad->SidesID[SmallFace] - nMPISides_MINE_Proc[iNBProc]) +
                                    offsetMPISides_YOUR[iNBProc];
                            dataquad->flips[SmallFace] = flip;
                            if (dataquad->SidesID[SmallFace] <= 0) {
                                printf("******************* ERROR!!! \n");
                                fflush(stdout);
                                exit(1);
                            }
                        }
                    } else {
                        if (dataquad->SidesID[SmallFace] <= nMPISides_YOUR_Proc[iNBProc]) {//MINE
                            dataquad->SidesID[SmallFace] += offsetMPISides_YOUR[iNBProc];
                            dataquad->flips[SmallFace] = flip;
                            if (dataquad->SidesID[SmallFace] <= 0) {
                                printf("******************* ERROR!!! \n");
                                fflush(stdout);
                                exit(1);
                            }
                        } else { //YOUR
                            dataquad->SidesID[SmallFace] =
                                    (dataquad->SidesID[SmallFace] - nMPISides_YOUR_Proc[iNBProc]) +
                                    offsetMPISides_MINE[iNBProc];
                            dataquad->flips[SmallFace] = 0;
                            if (dataquad->SidesID[SmallFace] <= 0) {
                                printf("******************* ERROR!!! \n");
                                fflush(stdout);
                                exit(1);
                            }
                        }
                    }
                }
            }

        }
    } else //Then it is a one side
    {

        if (side[0]->is.full.is_ghost || side[1]->is.full.is_ghost) {
            int SideIn = 0;
            int SideGhost = 1;
            if (side[0]->is.full.is_ghost) {
                SideIn = 1;
                SideGhost = 0;
            }
            int face = side[SideIn]->face;

            int face1 = side[SideGhost]->face;
            int flip = GetHFlip(face, face1, orientation);
            int ghostid = side[SideGhost]->is.full.quadid;
            int NBProc1 = ghost_to_proc[ghostid];
            int NBProc2 = NbProc[iNBProc];

            if (NBProc1 != NBProc2)
                return;
            // I need direction from
            quad = side[SideIn]->is.full.quad;
            dataquad = (p4est_inner_data_t *) quad->p.user_data;
            dataquad->SidesID[face] = (*SideID)++;


            if (myrank < NBProc1) {
                if (dataquad->SidesID[face] <= nMPISides_MINE_Proc[iNBProc]) { //MINE
                    dataquad->SidesID[face] += offsetMPISides_MINE[iNBProc];
                    dataquad->flips[face] = 0;
                    if (dataquad->SidesID[face] <= 0) {
                        printf("******************* ERROR!!! \n");
                        fflush(stdout);
                        exit(1);
                    }
                } else { //YOUR
                    dataquad->SidesID[face] =
                            (dataquad->SidesID[face] - nMPISides_MINE_Proc[iNBProc]) +
                            offsetMPISides_YOUR[iNBProc];
                    if (dataquad->SidesID[face] <= 0) {
                        printf("******************* ERROR!!! \n");
                        fflush(stdout);
                        exit(1);
                    }

                    dataquad->flips[face] = flip;
                }
            } else //if (myrank < NBProc1)
            {
                if (dataquad->SidesID[face] <= nMPISides_YOUR_Proc[iNBProc]) {
                    dataquad->SidesID[face] += offsetMPISides_YOUR[iNBProc];
                    dataquad->flips[face] = flip;

                    // dataquad->flips[face] = 0;
                    if (dataquad->SidesID[face] <= 0) {
                        printf("******************* ERROR!!! \n");
                        fflush(stdout);
                        exit(1);
                    }
                } else {
                    dataquad->SidesID[face] =
                            (dataquad->SidesID[face] - nMPISides_YOUR_Proc[iNBProc]) +
                            offsetMPISides_MINE[iNBProc];
                    dataquad->flips[face] = 0;
                    if (dataquad->SidesID[face] <= 0) {
                        printf("******************* ERROR!!! \n");
                        fflush(stdout);
                        exit(1);
                    }
                }
            }
        }

    }
}


#ifndef NON_OPTIMIZED
///////////////////////////////////////////////////////////////////////////////
//
// Set the SidesRatio. = 0 is non hanging Side. -1 = Small Side, +1 = big Side
//  it is not used  the ghost layer.
//
///////////////////////////////////////////////////////////////////////////////
static void
SetSidesRatio(p4est_iter_face_info_t *info, void *user_data)
{
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;

    int Changes = 0;
   p4est_quadrant_t *quad;
    p4est_inner_data_t *dataquad;

    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);

    if (sides->elem_count == 1)
    {
        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        quad = side[0]->is.full.quad;
        int face = side[0]->face;
        dataquad = (p4est_inner_data_t *)quad->p.user_data;
        dataquad->SidesRatio[face] = 0;// For Boundary always 0.
        return;
    }
    // side[i]->is.hanging.quad[j]->p.user_data;
    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    if (side[0]->is_hanging || side[1]->is_hanging)
    {
        //One Side is Mortar
        iBigSide = side[0]->is_hanging == 0 ? 0 : 1;
        iSmallSide = side[0]->is_hanging != 0 ? 0 : 1;
        int Bigface = side[iBigSide]->face;
        int SmallFace = side[iSmallSide]->face;
        if (side[iBigSide]->is.full.is_ghost == 0) //Big side is not MPI
        {
            // int flip = GetHFlip(BigFace, SmallFace, orientation);
            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0)
                {
                    // not ghost side, but local
                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = (p4est_inner_data_t *)quad->p.user_data;
                    dataquad->SidesRatio[SmallFace] = -1;

                }
                else
                {
                }
            }
            quad = side[iBigSide]->is.full.quad;
            dataquad = (p4est_inner_data_t *)quad->p.user_data;
            dataquad->SidesRatio[Bigface] = 1;
        }
        else
        { //We have only 4 small sides
            // Nothing to do
            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0)
                {
                    // not ghost side, but local
                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = (p4est_inner_data_t *)quad->p.user_data;
                   dataquad->SidesRatio[SmallFace] = -1;
                    // continue;
                }
            }
        return;
        }
    }
    else //Then it is a one side
    {
        if (side[0]->is.full.is_ghost || side[1]->is.full.is_ghost)
        {
            int SideIn = 0;
            int SideGhost = 1;
            if (side[0]->is.full.is_ghost)
            {
                SideIn = 1;
                SideGhost = 0;
            }
            quad = side[SideIn]->is.full.quad;
            int face = side[SideIn]->face;
            dataquad = (p4est_inner_data_t *)quad->p.user_data;

            dataquad->SidesRatio[face] = 0;
        }
        else //No Ghosts
        {
            int face = side[0]->face;
            quad = side[0]->is.full.quad;
            dataquad = (p4est_inner_data_t *)quad->p.user_data;
            dataquad->SidesRatio[face] = 0;

            quad = side[1]->is.full.quad;
            face = side[1]->face;
            dataquad = (p4est_inner_data_t *)quad->p.user_data;
            dataquad->SidesRatio[face] = 0;
        }

            return;
    }
}
#endif //  NON_OPTIMIZED
#ifndef NON_OPTIMIZED
///////////////////////////////////////////////////////////////////////////////
//
// Check if two or five NB quads was changed - then set OldSidesID to Zero.
// We use also the ghost layer.
//
///////////////////////////////////////////////////////////////////////////////
static void
CheckChanges(p4est_iter_face_info_t *info, void *user_data)
{
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;
    int myrank = p4est->mpirank;
    int Changes = 0;
    p4est_quadrant_t *quad;
    p4est_inner_data_t *dataquad;

    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);

    if (sides->elem_count == 1)
    {
        // There is nothing to do here with bounadry;
        return;
    }
    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    if (side[0]->is_hanging || side[1]->is_hanging)
    {
        //One Side is Mortar
        iBigSide = side[0]->is_hanging == 0 ? 0 : 1;
        iSmallSide = side[0]->is_hanging != 0 ? 0 : 1;
        int Bigface = side[iBigSide]->face;
        int SmallFace = side[iSmallSide]->face;
        if (side[iBigSide]->is.full.is_ghost == 0) //Big side is not MPI
        {
            // int flip = GetHFlip(BigFace, SmallFace, orientation);
            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0)
                {
                     // not ghost side, but local
                     quad = side[iSmallSide]->is.hanging.quad[j];
                     dataquad = (p4est_inner_data_t *)quad->p.user_data;
                    if ((dataquad->SidesRatio[SmallFace] != -1) 
                        || (dataquad->OldElementID[0] < 0) 
                        || (dataquad->OldElementID[1] > 0)
                        || (dataquad->OldSidesID[SmallFace]<=0))
                        dataquad->IsChanged++;

                    // continue;
                }
                else
                {

                    // MPI MINE Side;
                }
            }
            quad = side[iBigSide]->is.full.quad;
            dataquad = (p4est_inner_data_t *)quad->p.user_data;
            // dataquad = (p4est_inner_data_t *)quad->p.user_data;
            // quad = side[iSmallSide]->is.hanging.quad[j];
            if ((dataquad->SidesRatio[Bigface] != 1) 
                || (dataquad->OldElementID[0] < 0) 
                || (dataquad->OldElementID[1] > 0) 
                || (dataquad->OldSidesID[Bigface] <= 0))
                dataquad->IsChanged++;

           

        }
        else
        { //We have only 4 small sides
            // Nothing to do
            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0)
                {
                    // not ghost side, but local
                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = (p4est_inner_data_t *)quad->p.user_data;

                    if ((dataquad->SidesRatio[SmallFace] != -1) 
                        || (dataquad->OldElementID[0] < 0) 
                        || (dataquad->OldElementID[1] > 0) 
                        || (dataquad->OldSidesID[SmallFace] <= 0))
                        dataquad->IsChanged++;
                    // continue;
                }
            }
            return;
        }
    }
    else //Then it is a one side
    {    //Nothing to do here. We optimize Mortar side

      if (side[0]->is.full.is_ghost || side[1]->is.full.is_ghost)
        {
            int SideIn = 0;
            int SideGhost = 1;
            if (side[0]->is.full.is_ghost)
            {
                SideIn = 1;
                SideGhost = 0;
            }
            quad = side[SideIn]->is.full.quad;
            int face = side[SideIn]->face;
            dataquad = (p4est_inner_data_t *)quad->p.user_data;
            if ((dataquad->SidesRatio[face] != 0) 
            || (dataquad->OldElementID[0] < 0) 
            || (dataquad->OldElementID[1] > 0) 
            || (dataquad->OldSidesID[face] <= 0))
                dataquad->IsChanged++;
        }
        else //No Ghosts
        {
            int face = side[0]->face;
            quad = side[0]->is.full.quad;
            dataquad = (p4est_inner_data_t *)quad->p.user_data;
            // dataquad->SidesRatio[face] = 0;
            if ((dataquad->SidesRatio[face] != 0) 
            || (dataquad->OldElementID[0] < 0) 
            || (dataquad->OldElementID[1] > 0) 
            || (dataquad->OldSidesID[face] <= 0))
                dataquad->IsChanged++;

            quad = side[1]->is.full.quad;
            face = side[1]->face;
            dataquad = (p4est_inner_data_t *)quad->p.user_data;
             if ((dataquad->SidesRatio[face] != 0) 
            || (dataquad->OldElementID[0] < 0) 
            || (dataquad->OldElementID[1] > 0) 
            || (dataquad->OldSidesID[face] <= 0))
                dataquad->IsChanged++;
            // dataquad->SidesRatio[face] = 0;
        }

        return;
    }
}

#endif

static void
ShiftMPISideNumeration(p4est_iter_face_info_t *info, void *user_data) {
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;
    p4est_inner_data_t *ghost_data = (p4est_inner_data_t *) user_data;
    p4est_aux_data_t *aux_data = (p4est_aux_data_t *) p4est->user_pointer;
    int myrank = p4est->mpirank;
    int AddToMortarInnerSide = *(int *) aux_data->tmp[0]; //could be changed

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
    // side[i]->is.hanging.quad[j]->p.user_data;
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
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0) { //Not Ghost


                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = quad->p.user_data;
                    if (dataquad->SidesID[SmallFace] > 0)
                        dataquad->SidesID[SmallFace] += AddToMortarInnerSide;
                    //It must be Error here, because SidesID =0 as default and not changed yet 
                    // continue;
                } else {
                    int ghostid = side[iSmallSide]->is.hanging.quadid[j];
                    if (ghost_data[ghostid].SidesID[SmallFace] > 0)
                        ghost_data[ghostid].SidesID[SmallFace] += AddToMortarInnerSide;
                    // MPI MINE Side;
                }
            }
        } else { //We have only 4 small sides

            for (j = 0; j < P8EST_HALF; j++) //Check if the side is MPI
            {
                if (side[iSmallSide]->is.hanging.is_ghost[j] == 0) { //Not Ghost

                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = quad->p.user_data;
                    if (dataquad->SidesID[SmallFace] > 0)
                        dataquad->SidesID[SmallFace] += AddToMortarInnerSide;
                    // continue;
                } else {
                    //Nothing to do with SmallMPI Sides in this case
                    // MPI MINE Side;
                }
            }
            return;
        }
    } else //Then it is a one side
    {
        int SideIn = 0;
        int SideGhost = 1;
        if (side[0]->is.full.is_ghost) {
            SideIn = 1;
            SideGhost = 0;
        }
        int face = side[SideIn]->face;
        // int ghostid = side[SideGhost]->is.full.quadid;

        quad = side[SideIn]->is.full.quad;
        dataquad = (p4est_inner_data_t *) quad->p.user_data;

        if (dataquad->SidesID[face] > 0)
            dataquad->SidesID[face] += AddToMortarInnerSide;
    }
}

p4est_t *p4est_new_f(p4est_connectivity_t **conn) {
    // printf("conn %p \n",conn);s
    fflush(stdout);
    p4est_t *p4;
    // p4 = p4est_new(mpicomm, *conn, 0, NULL, NULL);
    p4 = p4est_new_ext(mpicomm,                    /* communicator */
                       *conn,                      /* connectivity */
                       0,                          /* minimum quadrants per MPI process */
                       0,                          /* minimum level of refinement */
                       1,                          /* fill uniform */
                       sizeof(p4est_inner_data_t), /* data size */
                       NULL,                       /* initializes data function*/
                       NULL);                      /* context */
    int nElems = 0;
    // #ifndef NON_OPTIMIZED
    p4est_iterate(p4,                     /* the forest */
                  NULL,                   /* the ghost layer May be LAter!!! */
                  (void *) &nElems,         /* the synchronized ghost data */
                  ElementCounterNew_iter, /* callback to compute each quad's
                                             interior contribution to du/dt */
                  NULL, //SetSidesRatio,                   /* SidesCount_iter,            /* callback to compute each quads'
            //              faces' contributions to du/du */
                  NULL,                   /* there is no callback for the
                                             edges between quadrants */
                  NULL);                  /* there is no callback for the
                                             corners between quadrants */
    // free(p4est_fortran_data);
// #endif //  NON_OPTIMIZED
    return p4;
    // return NULL;
}

static void
ElementNumberChanges(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    p4est_fortran_data_t *data = (p4est_fortran_data_t *) user_data;
    p4est_inner_data_t *dataquad = (p4est_inner_data_t *) q->p.user_data;
    int i = 0;
    int Fside = 0;
    int iElem = dataquad->ElementID;
    int *ChangeElements = (int *) data->ChngElementPtr;
    for (i = 0; i < 8; ++i) {
        ChangeElements[(iElem - 1) * 8 + (i)] = dataquad->OldElementID[i];
    }
#ifndef NON_OPTIMIZED
    // if (dataquad->IsChanged == 0)
    // {
    //     // Not Changed
    //     i=1;
    //     ChangeElements[(iElem - 1) * 8 + (i)] = -1;
    //     for (i = 2; i < 8; i++)
    //     { /// Convertion of Sides!!!!!!!!!
    //         Fside = P2H_side[i-2]; //Hopr ID of Side 1..6
    //         ChangeElements[(iElem - 1) * 8 + (i)] = dataquad->OldSidesID[Fside-1];
    //     }
    // }
#endif //  NON_OPTIMIZED

    // printf("iElem = %d, i = %d\n", iElem, dataquad->OldElementID[0]);
    // fflush(stdout);
    return;
}

p4est_fortran_data_t *GetnNBProcs(p4est_t *p4est, void *FortranData) {
    p4est_fortran_data_t *p4est_fortran_data = (p4est_fortran_data_t *) FortranData;
    // p4est_fortran_data = (p4est_fortran_data_t *)malloc(sizeof(p4est_fortran_data_t));
    p8est_ghost_t *ghost = p4est_fortran_data->ghost_ptr = p8est_ghost_new(p4est, P8EST_CONNECT_FACE);
    p4est_fortran_data->ghost_data_ptr = (p4est_inner_data_t *) malloc(
            ghost->ghosts.elem_count * sizeof(p4est_inner_data_t));
    int rank = 0;
    p4est_locidx_t jl;
    p4est_locidx_t num_ghost = 0;
    int nNBProcs = 0;
    int iElem;
    num_ghost = (p4est_locidx_t) ghost->ghosts.elem_count;
    p4est_fortran_data->ghost_to_proc = (int *) malloc(num_ghost * sizeof(int));

    for (jl = 0; jl < num_ghost; ++jl) {

        while (ghost->proc_offsets[rank + 1] <= jl) {
            ++rank;
            P4EST_ASSERT(rank < p4est->mpisize);
        }
        p4est_fortran_data->ghost_to_proc[jl] = rank;
    }
    if (num_ghost > 0) nNBProcs = 1;
    int iGh;
    for (iGh = 1; iGh < num_ghost; ++iGh) {
        if (p4est_fortran_data->ghost_to_proc[iGh] != p4est_fortran_data->ghost_to_proc[iGh - 1])
            // if (ghost->proc_offsets[rank + 1] != ghost->proc_offsets[rank])
            nNBProcs++;
    }


    p4est_fortran_data->nNBProcs = nNBProcs;
    
    return p4est_fortran_data;
}

void FillElemsChange(p4est_t *p4est, void *FortranData){

    p4est_fortran_data_t *p4est_fortran_data = (p4est_fortran_data_t *) FortranData;
    int i, iElem;
    int nElems = p4est_fortran_data->nElems;
    int *ChangeElements = p4est_fortran_data->ChngElementPtr;// = (int*) malloc(8 * nElems * sizeof(int));
    for (iElem = 1; iElem <= nElems; ++iElem)
        for (i = 1; i <= 8; i++) {
            ChangeElements[(iElem - 1) * 8 + (i - 1)] = 0; 
        }
    //Set the ElementID and OldElement, OldElement[0]=ElementID, OldElement[1:7]=0
    // // Count  Elements and Sides
    p4est_fortran_data->nElems = 0;
    p4est_iterate(p4est,                     
                  NULL,                      
                  (void *) p4est_fortran_data, 
                  ElementCounter_iter,
                  NULL,          
                  NULL,          
                  NULL);         
    p4est_iterate(p4est,
                  NULL,
                  (void *) p4est_fortran_data,
                  ElementNumberChanges,
                  NULL,
                  NULL,
                  NULL);

   // p4est_fortran_data->nElems = p4est->local_num_quadrants;
    return;
}   

int GetNElems(p4est_t *p4est){
    return p4est->local_num_quadrants;
}

void SetEtSandStE(p4est_t *p4est, p4est_fortran_data_t *p4est_fortran_data) {
    // Allocate EtS and StE
    int i, j, iElem, iSide;
    i = 1;
    int nElems = p4est_fortran_data->nElems;
    int nSides = p4est_fortran_data->nSides;
    p4est_inner_data_t *ghost_data = p4est_fortran_data->ghost_data_ptr;
    p8est_ghost_t *ghost = p4est_fortran_data->ghost_ptr;

    int *EtS = p4est_fortran_data->EtSPtr;//= (int *)malloc(2 * 6 * nElems * sizeof(int));
    // for (iElem = 1; iElem <= nElems; ++iElem)
    //     for (j = 1; j <= 6; j++)
    //         for (i = 1; i <= 2; i++)
    //         {
    //             // EtS[i-1+(j-1)*6+(iElem-1)*6*2]=i+j*10+iElem*1000;
    //             EtS[(iElem - 1) * 6 * 2 + (j - 1) * 2 + (i - 1)] = -1; //iElem + j*1000 + i * 10000;
    //         }

    int *StE = p4est_fortran_data->StEPtr;//! = (int *)malloc(5 * nSides * sizeof(int));
    StE[1] = -1;
    // printf("Pointer  = %p \n!!", StE);

    for (iSide = 1; iSide <= nSides; ++iSide)
        for (i = 1; i <= 5; i++) {
            // EtS[i-1+(j-1)*6+(iElem-1)*6*2]=i+j*10+iElem*1000;
            StE[(iSide - 1) * 5 + (i - 1)] = -1; //iElem + j*1000 + i * 10000;
        }

    // Now Mortar Type
    int *MoTy = p4est_fortran_data->MTPtr;//! = (int *)malloc(2 * nSides * sizeof(int));
    for (iSide = 1; iSide <= nSides; ++iSide)
        for (i = 1; i <= 2; i++) {
            // EtS[i-1+(j-1)*6+(iElem-1)*6*2]=i+j*10+iElem*1000;
            MoTy[(iSide - 1) * 2 + (i - 1)] = -1; //iElem + j*1000 + i * 10000;
        }

    // Now MortarInfo
    int nMortarSides = p4est_fortran_data->nMortarInnerSides + p4est_fortran_data->nMortarMPISides;

    int *MoInf = p4est_fortran_data->MIPtr = (int *) malloc(2 * 4 * nMortarSides * sizeof(int));
    for (iSide = 1; iSide <= nMortarSides; ++iSide)
        for (i = 1; i <= 2; i++)
            for (j = 1; i <= 4; i++) {
                // EtS[i-1+(j-1)*6+(iElem-1)*6*2]=i+j*10+iElem*1000;
                // MoInf[(iSide - 1) * 5 + (i - 1)] = -1; //iElem + j*1000 + i * 10000;
                MoInf[(iSide - 1) * 4 * 2 + (j - 1) * 2 + (i - 1)] = -1;
            }

    int *BCs = p4est_fortran_data->BCs = (int *) malloc(p4est_fortran_data->nBCSides * sizeof(int32_t));
    for (iSide = 1; iSide <= p4est_fortran_data->nBCSides; ++iSide) {
        BCs[iSide - 1] = -1;
    }
    p4est->user_pointer = (void *) ghost_data;
    p4est_iterate(p4est,                      /* the forest */
                  ghost,                      /* the ghost layer May be LAter!!! */
                  (void *) p4est_fortran_data, /* the synchronized ghost data */
                  SetElementToSide_iter,      /* callback to compute each quad's
   //                           interior contribution to du/dt */
                  SetSideToElement_iter,      /* SidesCount_iter,            /* callback to compute each quads'
     //                         faces' contributions to du/du */
                  NULL,                       /* there is no callback for the
                              edges between quadrants */
                  NULL);
    p4est->user_pointer = NULL;
//     int *ChangeElements = p4est_fortran_data->ChngElementPtr;// = (int*) malloc(8 * nElems * sizeof(int));
//     for (iElem = 1; iElem <= nElems; ++iElem)
//         for (i = 1; i <= 8; i++) {
//             // EtS[i-1+(j-1)*6+(iElem-1)*6*2]=i+j*10+iElem*1000;
//             ChangeElements[(iElem - 1) * 8 + (i - 1)] = 0; //iElem + j*1000 + i * 10000;
//         }
// #ifndef NON_OPTIMIZED
//     // p4est_iterate(p4est,
//     //               NULL,
//     //               (void *)p4est_fortran_data,
//     //               NULL,

//     //               CheckChanges,

//     //               NULL,

//     //               NULL);
// #endif

//     p4est_iterate(p4est,
//                   NULL,
//                   (void *) p4est_fortran_data,
//                   ElementNumberChanges,
//                   NULL,
//                   NULL,
//                   NULL);

    pfree(ghost_data);
    p8est_ghost_destroy(ghost);
    return;

}

void GetData(p4est_t *p4est, p4est_fortran_data_t *p4est_fortran_data) {
    int i1, j1 = 0, k1 = 0;
    int rank = 0;
    p4est_locidx_t local_num_quad, num_ghost = 0;

    p8est_ghost_t *ghost = p4est_fortran_data->ghost_ptr;
    p4est_inner_data_t *ghost_data = p4est_fortran_data->ghost_data_ptr;
    int *ghost_to_proc = p4est_fortran_data->ghost_to_proc;
    int nNBProcs = p4est_fortran_data->nNBProcs;


    /* Calculate ghost information */

    int nElems;// = p4est_fortran_data->nElems = 0;
    int nSides = p4est_fortran_data->nSides = 0;
    p4est_fortran_data->nGlobalElems = p4est->global_num_quadrants;
    p4est_fortran_data->nElems = p4est->local_num_quadrants;
    p4est_fortran_data->nBCSides = 0;
    p4est_fortran_data->nMortarInnerSides = 0;
    p4est_fortran_data->nInnerSides = 0;
    p4est_fortran_data->nMPISides = 0;
    p4est_fortran_data->nMortarMPISides = 0;

    // if (p4est->mpisize > 1)
    {
        int i;
        for (i = 0; i < nNBProcs; i++) {
            p4est_fortran_data->nMPISides_Proc[i] = 0;
        }

        int iProc = 0;
        if (p4est->mpisize > 1) {
            for (rank = 0; rank < p4est->mpisize; rank++) {
                if (ghost->proc_offsets[rank + 1] != ghost->proc_offsets[rank]) {

                    p4est_fortran_data->nNbProc[iProc] = rank;
                    iProc++;
                }
            }
        }
    }


    //Set the ElementID and OldElement, OldElement[0]=ElementID, OldElement[1:7]=0
    // // Count  Elements and Sides
    p4est_iterate(p4est,                      /* the forest */
                  ghost,                      /* the ghost layer May be LAter!!! */
                  (void *) p4est_fortran_data, /* the synchronized ghost data */
                  NULL,//ElementCounter_iter,        /* callback to compute each quad's
                                          //   interior contribution to du/dt */
                  SidesCounter_iter,          /* callback to compute each quads'
                                             faces' contributions to du/du */
                  NULL,                       /* there is no callback for the
                                             edges between quadrants */
                  NULL);                      /* there is no callback for the
                                  // // Count  Elements and Sides*/


    nElems = p4est_fortran_data->nElems;
    nSides = p4est_fortran_data->nSides;


    // Set Sides Number in the Quad data
    // Use special struct called Counter - because every number begins not from the 0. e.g. nMortarSides begins with nBCSides

    p4est_SetMPISide_data_t SetMPISideData_t;
    p4est_SetMPISide_data_t *SetMPISideData = &SetMPISideData_t;
    SetMPISideData->nMPISidesCount = nullptr;

    p4est_aux_data_t aux_data_t;
    p4est_aux_data_t *aux_data = &aux_data_t;


    // if (p4est->mpisize > 1)
    {

        SetMPISideData->nMPISides_Proc = p4est_fortran_data->nMPISides_Proc;
        SetMPISideData->ghost_to_proc = ghost_to_proc;
        SetMPISideData->nMPISidesCount = (int *) malloc(p4est->mpisize * sizeof(int));
        // int INBProc = 0;
        for (rank = 0; rank < p4est->mpisize; rank++) {
            SetMPISideData->nMPISidesCount[rank] = 0;
        }
        //Size of 0:MPISIZE-1
        if (p4est->mpisize > 1) {
            //Set the ElementID and OldElement, OldElement[0]=ElementID, OldElement[1:7]=0
            // // Count  Elements and Sides
            p4est_iterate(p4est,                  /* the forest */
                          ghost,                  /* the ghost layer May be LAter!!! */
                          (void *) SetMPISideData, /* the synchronized ghost data */
                          NULL,                   /* callback to compute each quad's
                                             interior contribution to du/dt */
                          SidesCounter2_iter,     /* SidesCount_iter,            /* callback to compute each quads'
                                             faces' contributions to du/du */
                          NULL,                   /* there is no callback for the
                                             edges between quadrants */
                          NULL);                  /* there is no callback for the
                                             corners between quadrants */

        }
        int iNB;
        for (iNB = 0; iNB < nNBProcs; iNB++) {
            int rank = p4est_fortran_data->nNbProc[iNB];
            if (rank == p4est->mpirank)
                continue;

            if (SetMPISideData->nMPISidesCount[rank] > 0) {
                int iNBProc = p4est_fortran_data->nNbProc[iNB];

                SetMPISideData->nMPISides_Proc[iNB] = SetMPISideData->nMPISidesCount[rank];
                // printf("!==!==!= iNB = %d , SetMPISideData->nMPISides_Proc[iNB] = %d , iNBProc = %d \n", iNB, SetMPISideData->nMPISides_Proc[iNB], iNBProc);
            };
        }

        SetMPISideData->nMPISides_MINE_Proc = p4est_fortran_data->nMPISides_MINE_Proc;//(int *)malloc(nNBProcs * sizeof(int));
        SetMPISideData->nMPISides_YOUR_Proc = p4est_fortran_data->nMPISides_YOUR_Proc;//(int *)malloc(nNBProcs * sizeof(int));
        int i;
        for (i = 0; i < nNBProcs; i++) {
            SetMPISideData->nMPISides_MINE_Proc[i] = 0;
            SetMPISideData->nMPISides_YOUR_Proc[i] = 0;
        }
        // Now nMPISides_MINE_Proc and nMPISides_YOUR_Proc
        // nMPISides_MINE=SUM(nMPISides_MINE_Proc)
        int iNbProc = 0;
        for (iNbProc = 0; iNbProc < nNBProcs; iNbProc++) {
            if (p4est->mpirank < p4est_fortran_data->nNbProc[iNbProc]) {
                SetMPISideData->nMPISides_MINE_Proc[iNbProc] = SetMPISideData->nMPISides_Proc[iNbProc] / 2;
            } else {
                SetMPISideData->nMPISides_MINE_Proc[iNbProc] =
                        SetMPISideData->nMPISides_Proc[iNbProc] - SetMPISideData->nMPISides_Proc[iNbProc] / 2;
            }
            SetMPISideData->nMPISides_YOUR_Proc[iNbProc] =
                    SetMPISideData->nMPISides_Proc[iNbProc] - SetMPISideData->nMPISides_MINE_Proc[iNbProc];
        }
        SetMPISideData->nMPISides_MINE = 0;
        SetMPISideData->nMPISides_YOUR = 0;
        //nMPISides_MINE = SUM(nMPISides_MINE_Proc)
        // nMPISides_YOUR = SUM(nMPISides_YOUR_Proc)
        for (iNbProc = 0; iNbProc < nNBProcs; iNbProc++) {
            SetMPISideData->nMPISides_MINE += SetMPISideData->nMPISides_MINE_Proc[iNbProc];
            SetMPISideData->nMPISides_YOUR += SetMPISideData->nMPISides_YOUR_Proc[iNbProc];

        }


        // Set offsetMPISides_YOUR
        // and offsetMPISides_MINE
        SetMPISideData->offsetMPISides_MINE = p4est_fortran_data->offsetMPISides_MINE;
        SetMPISideData->offsetMPISides_YOUR = p4est_fortran_data->offsetMPISides_YOUR;
        if (nNBProcs > 0) {
            for (iNbProc = 0; iNbProc <= nNBProcs; iNbProc++) {
                SetMPISideData->offsetMPISides_MINE[iNbProc] = 0;
                SetMPISideData->offsetMPISides_YOUR[iNbProc] = 0;
            }

            SetMPISideData->offsetMPISides_MINE[0] = p4est_fortran_data->nInnerSides + p4est_fortran_data->nBCSides +
                                                     p4est_fortran_data->nMortarInnerSides;
            for (iNbProc = 0; iNbProc < nNBProcs; iNbProc++) {
                SetMPISideData->offsetMPISides_MINE[iNbProc + 1] =
                        SetMPISideData->offsetMPISides_MINE[iNbProc] + SetMPISideData->nMPISides_MINE_Proc[iNbProc];
            }

            SetMPISideData->offsetMPISides_YOUR[0] = SetMPISideData->offsetMPISides_MINE[nNBProcs];
            for (iNbProc = 0; iNbProc < nNBProcs; iNbProc++) {
                SetMPISideData->offsetMPISides_YOUR[iNbProc + 1] =
                        SetMPISideData->offsetMPISides_YOUR[iNbProc] + SetMPISideData->nMPISides_YOUR_Proc[iNbProc];
            }
        }

        for (i = 0; i < num_ghost; i++) {
            int Side = 0;
            for (Side = 0; Side < 6; Side++)
                ghost_data[i].SidesID[Side] = 0;
        }

        int iNBProc = 0;
        for (iNBProc = 0; iNBProc < nNBProcs; iNBProc++) // Main cycle over NB Processors
        {
            // First Step is to Set up SideID for MPI Sidebut Proc-by-Proc
            int SideID;
            SideID = 1;
            aux_data->tmp[0] = &SideID;
            aux_data->tmp[1] = p4est_fortran_data->nNbProc;
            aux_data->tmp[2] = SetMPISideData->offsetMPISides_MINE;
            aux_data->tmp[3] = SetMPISideData->offsetMPISides_YOUR;
            aux_data->tmp[4] = p4est_fortran_data->nMPISides_Proc;
            aux_data->tmp[5] = ghost_to_proc;
            aux_data->tmp[6] = &iNBProc;
            aux_data->tmp[7] = SetMPISideData->nMPISides_MINE_Proc;
            aux_data->tmp[8] = SetMPISideData->nMPISides_YOUR_Proc;

            p4est->user_pointer = aux_data;
            p4est_iterate(p4est,                /* the forest */
                          ghost,                /* Ghost*/
                          (void *) ghost_data,   /* the non-synchronized ghost data */
                          NULL,                 /* ElemIter */
                          SetMPISidesAUXNumber, /* SidesIter, */
                          NULL,                 /* EDGE Iter */
                          NULL);                /* CORNER Iter */
        }
        pfree(ghost_to_proc);
        // 1. Optimisation of MPIMortar side with MPI MINE small Side (shift to MPIInnerSide)
        int AddToMortarInnerSide = 0; //How many Mortar sides must be shifted to Innner.

        //if (OPi)
        {
            // Zero point for shifted MPIMortar Side
            int MortarSideNumeration = p4est_fortran_data->nBCSides + p4est_fortran_data->nMortarInnerSides;
            int FirstMPIMINESide = MortarSideNumeration + p4est_fortran_data->nInnerSides; //Greater >
            int LastMPIMINESide = FirstMPIMINESide + SetMPISideData->nMPISides_MINE;       // and Less <=

            aux_data->tmp[0] = &AddToMortarInnerSide;
            aux_data->tmp[1] = &MortarSideNumeration;
            aux_data->tmp[2] = &FirstMPIMINESide;
            aux_data->tmp[3] = &LastMPIMINESide;
            p4est->user_pointer = aux_data;


            p4est_iterate(p4est,              /* the forest */
                          ghost,              /* Ghost*/
                          (void *) ghost_data, /* the non-synchronized ghost data */
                          NULL,               /* ElemIter */
                          MortarOptimisation, /* SidesIter, */
                          NULL,               /* EDGE Iter */
                          NULL);              /* CORNER Iter */
        }
        // 2. Shift the Numeration of the rest of the sides

        if (AddToMortarInnerSide > 0) {
            aux_data->tmp[0] = &AddToMortarInnerSide;
            p4est->user_pointer = aux_data;
            p4est_iterate(p4est,                  /* the forest */
                          ghost,                  /* Ghost*/
                          (void *) ghost_data,     /* the non-synchronized ghost data */
                          NULL,                   /* ElemIter */
                          ShiftMPISideNumeration, /* SidesIter, */
                          NULL,                   /* EDGE Iter */
                          NULL);                  /* CORNER Iter */
        }


        // change Numeration for MPISide and arrays OffsetMPI and so on. see FLUXO
        p4est_fortran_data->nMortarInnerSides += AddToMortarInnerSide;
        p4est_fortran_data->nMortarMPISides -= AddToMortarInnerSide;
        for (iNbProc = 0; iNbProc < nNBProcs; iNbProc++) {

            SetMPISideData->offsetMPISides_MINE[iNbProc] += AddToMortarInnerSide;
            SetMPISideData->offsetMPISides_YOUR[iNbProc] += AddToMortarInnerSide;
        }
    }
    // Right now there is a Numeration for MPISides_MINE and MPISides_YOUR
    // Now we can set up the rest of Numeration

    p4est_SetSide_data_t SetSide_data_t;// = (p4est_SetSide_data_t *)malloc(sizeof(p4est_SetSide_data_t));

    p4est_SetSide_data_t *SetSide_data = &SetSide_data_t;

    SetSide_data->CurrentBCSide = 1;

    SetSide_data->CurrentMortrarInnerSide = SetSide_data->CurrentBCSide + p4est_fortran_data->nBCSides;

    SetSide_data->CurrentInnerSide = SetSide_data->CurrentMortrarInnerSide + p4est_fortran_data->nMortarInnerSides;
    SetSide_data->CurrentMPISide = SetSide_data->CurrentInnerSide + p4est_fortran_data->nInnerSides;
    SetSide_data->CurrentMPIMortarSide = SetSide_data->CurrentMPISide + p4est_fortran_data->nMPISides;

    // MPI Sides have been set up already

    if ((SetSide_data->CurrentMPIMortarSide + p4est_fortran_data->nMortarMPISides - 1) != nSides) {

        printf("SetSide_data->CurrentMortrarInnerSide = %d\n", SetSide_data->CurrentMortrarInnerSide);
        printf("SetSide_data->CurrentInnerSide = %d\n", SetSide_data->CurrentInnerSide);
        printf("SetSide_data->CurrentMPISide = %d\n", SetSide_data->CurrentMPISide);
        printf("SetSide_data->CurrentMPIMortarSide = %d\n", SetSide_data->CurrentMPIMortarSide);
        printf("p4est_fortran_data->nBCSides = %d\n", p4est_fortran_data->nBCSides);
        printf("p4est_fortran_data->nMortarInnerSides = %d\n", p4est_fortran_data->nMortarInnerSides);
        printf("p4est_fortran_data->nInnerSides; = %d\n", p4est_fortran_data->nInnerSides);
        printf("p4est_fortran_data->nMPISides; = %d\n", p4est_fortran_data->nMPISides);
        printf("p4est_fortran_data->nMortarMPISides = %d\n", p4est_fortran_data->nMortarMPISides);
        printf("Error with sides Numeration!!! nSides = %d \n", nSides);
        exit(1);
    }

    //... The other sides
    // The sum must be equal to nSides
    p4est_iterate(p4est,                /* the forest */
                  ghost,                 /* the ghost layer May be LAter!!! */
                  (void *) SetSide_data, /* the synchronized ghost data */
                  NULL,                 /* callback to compute each quad's
                                             interior contribution to du/dt */
                  SetSides_iter,
                  NULL,                 /* there is no callback for the
                                             edges between quadrants */
                  NULL);                /* there is no callback for the
    //                                          corners between quadrants */



    p4est->user_pointer = NULL;

    p4est_fortran_data->nMPISides_YOUR = SetMPISideData->nMPISides_YOUR;
    p4est_fortran_data->nMPISides_MINE = SetMPISideData->nMPISides_MINE;

    pfree(SetMPISideData->nMPISidesCount);

    return;
}


void RefineCoarse(p4est_t *p4est, void *ElemToRC) {
    int *ElemToRefineCoarse = (int *) ElemToRC;
    // p4est_fortran_data_t *back;
    // Refine And Coarse
    p4est->user_pointer = ElemToRC;
// #ifndef NON_OPTIMIZED
   

    p4est_iterate(p4est,                           /* the forest */
                  NULL,                            /* the ghost layer May be LAter!!! */
                  NULL,                            /* the synchronized ghost data */
                  ElementCounterSetOldToZero_iter, /* callback to compute each quad's
                                             interior contribution to du/dt */
                  NULL,// SetSidesRatio,                   /* SidesCount_iter,            /* callback to compute each quads'
            //                              faces' contributions to du/du */
                  NULL,                            /* there is no callback for the
                                             edges between quadrants */
                  NULL);

// #endif //  NON_OPTIMIZED
    int recursive = 0;
    int Callbackorphans = 0;
    int allowed_level = P4EST_QMAXLEVEL;
    p4est_refine_ext(p4est, recursive, allowed_level,
                     refine_fn, NULL,
                     replace_quads);
                     
    p4est_coarsen_ext(p4est, recursive, Callbackorphans,
                      coarse_fn, NULL, replace_quads);
    p4est_balance_ext(p4est, P4EST_CONNECT_FULL, NULL,
                      replace_quads);
    // p4est_balance_ext(p4est, P4EST_CONNECT_FACE, NULL,
    //                   replace_quads);
    p4est->user_pointer = NULL;
    fflush(stdout);

    return;//    return GetData(p4est);
}

// Build p4est data structures with existing connectivity from HDF5 mesh
void p4_connectivity_treevertex(p4est_topidx_t num_vertices,
                                p4est_topidx_t num_trees,
                                double *vertices,
                                p4est_topidx_t *tree_to_vertex,
                                p4est_topidx_t num_periodics,
                                p4est_topidx_t *join_faces,
                                p4est_connectivity_t **conn_out) {
    p4est_topidx_t tree;
    int face, i;
    p4est_connectivity_t *conn = NULL;

    P4EST_GLOBAL_PRODUCTIONF("creating connectivity from tree to vertex only...%d %d \n", num_vertices, num_trees);

    conn = p4est_connectivity_new(num_vertices, num_trees,
                                  0, 0,
                                  0, 0);
    for (i = 0; i < 3 * num_vertices; ++i) {
        conn->vertices[i] = vertices[i];
    }
    for (i = 0; i < 8 * num_trees; ++i) {
        conn->tree_to_vertex[i] = tree_to_vertex[i];
    }

    /*
   * Fill tree_to_tree and tree_to_face to make sure we have a valid
   * connectivity.
   */

    for (tree = 0; tree < conn->num_trees; ++tree) {
        for (face = 0; face < P4EST_FACES; ++face) {
            conn->tree_to_tree[P4EST_FACES * tree + face] = tree;
            conn->tree_to_face[P4EST_FACES * tree + face] = face;
        }
    }

    P4EST_ASSERT(p4est_connectivity_is_valid(conn));

    p4est_connectivity_complete(conn);

    P4EST_ASSERT(p4est_connectivity_is_valid(conn));
    // Join Faces
    if (num_periodics > 0) {
        for (i = 0; i < num_periodics; ++i) {
            p8est_connectivity_join_faces(conn, join_faces[5 * i], join_faces[5 * i + 1],
                                          join_faces[5 * i + 2], join_faces[5 * i + 3],
                                          join_faces[5 * i + 4]);
        }
    }

    P4EST_ASSERT(p4est_connectivity_is_valid(conn));

    P4EST_GLOBAL_PRODUCTIONF("New connectivity with %lld trees and %lld vertices\n",
                             (long long) conn->num_trees, (long long) conn->num_vertices);

    *conn_out = conn;
    // printf("Connectivity adress = %p\n", conn);
}


// int iElem=0;
static void
partition_init_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                  p4est_quadrant_t *q) {
    /* the data associated with a forest is accessible by user_pointer */
    p4est_balance_data_t *ctx = (p4est_balance_data_t *) p4est->user_pointer;
    /* the data associated with a quadrant is accessible by p.user_data */
    double *data = (double *) q->p.user_data;
    double *data1;
    int nElems = p4est->local_num_quadrants;
    int nVar = ctx->nVar;
    int PP_N = ctx->PP_N;
    int Datasize = ctx->DataSize;

    double *U = (double *) ctx->DataSetU;
    double *Elem_xGP = (double *) ctx->DataSetElem_xGP;

    int P1 = PP_N + 1;
    int P2 = P1 * P1;
    int P3 = P2 * P1;
    int D = 3;
    int i, j, k, iVar;
    int iElem = ctx->nElems;

    int Part = 0;
    i = j = k = 2;
    iVar = 1;
    // double a = U[(iElem - 1) * nVar * PP_Np13 + k * nVar * PP_Np12 + j * nVar * PP_Np1 + i*nVar + (iVar-1)];
    int dir = 1; //direction 1:3
    Part = 1;// This is U
    // double b = Elem_xGP[(iElem - 1) * D * P3 + k * D * P2 + j * D * P1 + i * D + (dir - 1)];
    double a = 0;
    for (iVar = 1; iVar <= nVar; iVar++)
        for (i = 0; i <= PP_N; i++)
            for (j = 0; j <= PP_N; j++)
                for (k = 0; k <= PP_N; k++) {
                    data[(Part - 1) * nVar * P3 + k * nVar * P2 + j * nVar * P1 + i * nVar + (iVar - 1)]
                            = U[(iElem - 1) * nVar * P3 + k * nVar * P2 + j * nVar * P1 + i * nVar + (iVar - 1)];
                }


    // data1 = &data[(Part - 1) * nVar * P3 + P3*nVar + P2*nVar + P1 * nVar + nVar - 1 ];// + PP_N * nVar * P2 + PP_N * nVar * P1 + PP_N * nVar + (nVar - 1)];
    data1 = &data[P3 *
                  nVar];// + P2 * nVar + P1 * nVar + nVar]; // + PP_N * nVar * P2 + PP_N * nVar * P1 + PP_N * nVar + (nVar - 1)];
    Part = 2; // This is Elem_xGP
    for (dir = 1; dir <= D; dir++)
        for (i = 0; i <= PP_N; i++)
            for (j = 0; j <= PP_N; j++)
                for (k = 0; k <= PP_N; k++) {
                    data1[k * D * P2 + j * D * P1 + i * D + (dir - 1)]
                            = Elem_xGP[(iElem - 1) * D * P3 + k * D * P2 + j * D * P1 + i * D + (dir - 1)];

                }


    ctx->nElems++;

}

static void
ReturnData(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    // p4est_inner_data_t *dataquad = (p4est_inner_data_t *)q->p.user_data;
    p4est_t *p4est = info->p4est;
    p4est_balance_data_t *ctx = (p4est_balance_data_t *) p4est->user_pointer;
    double *data = (double *) q->p.user_data;
    double *data1;
    int nElems = p4est->local_num_quadrants;
    int nVar = ctx->nVar;
    int PP_N = ctx->PP_N;
    int Datasize = ctx->DataSize;
    double *U = (double *) ctx->DataSetU;
    double *Elem_xGP = (double *) ctx->DataSetElem_xGP;

    int P1 = PP_N + 1;
    int P2 = P1 * P1;
    int P3 = P2 * P1;
    int D = 3;
    int i, j, k, iVar;
    int iElem = ctx->nElems;
    int Part = 0;
    i = j = k = 2;
    iVar = 1;
    int dir = 1; //direction 1:3
    Part = 1;    // This is U
    double a = 0;
    for (iVar = 1; iVar <= nVar; iVar++)
        for (i = 0; i <= PP_N; i++)
            for (j = 0; j <= PP_N; j++)
                for (k = 0; k <= PP_N; k++) {
                    U[(iElem - 1) * nVar * P3 + k * nVar * P2 + j * nVar * P1 + i * nVar + (iVar - 1)]
                            = data[(Part - 1) * nVar * P3 + k * nVar * P2 + j * nVar * P1 + i * nVar + (iVar - 1)];
                }

    data1 = &data[P3 *
                  nVar]; // + P2 * nVar + P1 * nVar + nVar]; // + PP_N * nVar * P2 + PP_N * nVar * P1 + PP_N * nVar + (nVar - 1)];
    Part = 2;                 // This is Elem_xGP
    for (dir = 1; dir <= D; dir++)
        for (i = 0; i <= PP_N; i++)
            for (j = 0; j <= PP_N; j++)
                for (k = 0; k <= PP_N; k++) {
                    // data[(Part - 1) * nVar * P3 + k * nVar * P2 + j * nVar * P1 + i * nVar + (dir - 1)]
                    //     = Elem_xGP[(iElem - 1) * D * P3 + k * D * P2 + j * D * P1 + i * D + (dir - 1)];
                    Elem_xGP[(iElem - 1) * D * P3 + k * D * P2 + j * D * P1 + i * D + (dir - 1)]
                            = data1[k * D * P2 + j * D * P1 + i * D + (dir - 1)];
                    // printf("a = %f \n",a);
                }
    ctx->nElems++;
}


void ResetElementNumber(p4est_t *p4est) {
    p4est_reset_data(p4est, sizeof(p4est_inner_data_t), NULL, NULL);
    int nEl = 0;
    p4est_iterate(p4est,
                  NULL,
                  (void *) &nEl,
                  ElementCounterNew_iter,
            // ElementCounterSetOldToZero_iter,
                  NULL,
                  NULL,
                  NULL);
    // ElementCounterNew_iter
    return;
}

// typedef void(* p4est_init_t)(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)

// typedef int(* p4est_weight_t)(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)

typedef struct p4est_Weights {
    int index;
    int8_t *Weights; // SideInfo in HDF5

} p4est_Weights_t;


static void
WeightsFunctionIter(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    p4est_inner_data_t *dataquad = (p4est_inner_data_t *) q->p.user_data;
    p4est_Weights_t *WeightStruct = (p4est_Weights_t *) user_data;
    WeightStruct->Weights[WeightStruct->index] = dataquad->weight;
    ++WeightStruct->index;
    return;
}


static int
weight_fn(p4est_t *p4est,
          p4est_topidx_t which_tree,
          p4est_quadrant_t *quad) {

    p4est_Weights_t *WeightStruct = (p4est_Weights_t *) p4est->user_pointer;
    int8_t Weight = WeightStruct->Weights[WeightStruct->index];
    ++WeightStruct->index;
  	return 1;
    return Weight;
}

// static int
// weight_fn(p4est_t *p4est,
//           p4est_topidx_t which_tree,
//           p4est_quadrant_t *quad)
// {
//         p4est_inner_data_t *dataquad = (p4est_inner_data_t *)quad->p.user_data;
//         return (int)dataquad->weight;
// }


void p4est_loadbalancing_init(p4est_t *p4est, void *user_pointer) {
    // p4set_gloidx_t  *dest; // before
    // p4set_gloidx_t  *src = (p4set_gloidx_t  *)desptr; // after
    p4est_balance_datav2_t *data = (p4est_balance_datav2_t *) user_pointer;
    p4est_gloidx_t *src = nullptr; //(p4set_gloidx_t  *)desptr; // after
    p4est_Weights_t WeightStruct;
    WeightStruct.index = 0;
    WeightStruct.Weights = (int8_t *) malloc(p4est->local_num_quadrants * sizeof(int8_t));

    p4est_iterate(p4est,
                  NULL,
                  (void *) &WeightStruct,
                  WeightsFunctionIter,
            // ElementCounterSetOldToZero_iter,
                  NULL,
                  NULL,
                  NULL);
    src = (p4est_gloidx_t *) malloc((p4est->mpisize + 1) * sizeof(p4est_gloidx_t));

    p4est_reset_data(p4est, 0, NULL, NULL);
    int i;
    for (i = 0; i <= p4est->mpisize; ++i) {
        src[i] = p4est->global_first_quadrant[i];

    }
    data->src_gfq = (p4est_gloidx_t *) src;
    p4est->user_pointer = (void *) &WeightStruct;
    WeightStruct.index = 0;
    const int allow_coarsening = 1;
    //p4est_partition(p4est, allow_coarsening, NULL);
    p4est_partition(p4est, allow_coarsening, weight_fn);
    p4est->user_pointer = nullptr;
    // printf(" data->src_gfq [1] = %d \n",  data->src_gfq[1]);
    data->nElems = p4est->local_num_quadrants;
    free(WeightStruct.Weights);
    // printf("mpirank = %d, nElems = %d \n", p4est->mpirank, data->nElems);
    return;
}

void p4est_loadbalancing_go(p4est_t *p4est, void *user_pointer) {

    p4est_balance_datav2_t *data = (p4est_balance_datav2_t *) user_pointer;

    p4est_transfer_fixed(p4est->global_first_quadrant, data->src_gfq, p4est->mpicomm, 8, data->Elem_xGP_new,
                         data->Elem_xGP_old, data->GPSize);
                         
    p4est_transfer_fixed(p4est->global_first_quadrant, data->src_gfq, p4est->mpicomm, 10, data->U_new, data->U_old,
                         data->DataSize);
    
    p4est_transfer_fixed(p4est->global_first_quadrant, data->src_gfq, p4est->mpicomm, 7, data->ElemWasCoarsened_new,
                         data->ElemWasCoarsened_old, data->CoarseSize);
    
#if SHOCK_NFVSE
    p4est_transfer_fixed(p4est->global_first_quadrant, data->src_gfq, p4est->mpicomm, 9, data->Alpha_new, data->Alpha_old,
                         sizeof(double));
#endif /*SHOCK_NFVSE*/

    //  printf("data->GPSize= %d \n", data->GPSize);
    // int a = MPI_Barrier(p4est->mpicomm);
    
    // int b = MPI_Barrier(p4est->mpicomm);
    free(data->src_gfq);
    // free(dest);

    return;

}


void p4est_loadbalancing(p4est_t *p4est, void *user_pointer) {
    p4est_balance_data_t *data = (p4est_balance_data_t *) user_pointer;
    int nVar = data->nVar;
    int PP_N = data->PP_N;
    int nElems = p4est->local_num_quadrants;
    data->nElems = 1;// We use this field as a counter
    int Datasize = data->DataSize;
    // Try to avoid copying the data

    p4est_t *p4est1 = p4est;//p4est_copy(p4est, 0); // 0 - no data to copy

    p4est_reset_data(p4est1, Datasize, partition_init_fn, user_pointer);
    int allow_coarsening = 1;
    p4est_partition(p4est1, allow_coarsening, NULL);
    // Return values to the two big array
    p4est1->user_pointer = user_pointer;

    nElems = p4est1->local_num_quadrants;
    // double *Elem_xGP = (double *)ctx->DataSetElem_xGP;
    // double *U = (double *)ctx->DataSetU;
    //Create to Array to return to the FLUXO
    int P1 = PP_N + 1;
    int P2 = P1 * P1;
    int P3 = P2 * P1;
    int Usize = nElems * nVar * P3;
    int GPsize = nElems * 3 * P3;
    data->DataSetU = (double *) malloc(Usize * sizeof(double));
    data->DataSetElem_xGP = (double *) malloc(GPsize * sizeof(double));
    data->nElems = 1;
    p4est_iterate(p4est1, /* the forest */
                  NULL,
                  NULL,
                  ReturnData,
                  NULL,
                  NULL,
                  NULL);
    // The return data is set up
    p4est_reset_data(p4est1, sizeof(p4est_inner_data_t), NULL, NULL);
    int nEl = 0;
    p4est_iterate(p4est1,                     /* the forest */
                  NULL,                   /* the ghost layer May be LAter!!! */
                  (void *) &nEl,        /* the synchronized ghost data */
                  ElementCounterNew_iter, /* callback to compute each quad's
                                             interior contribution to du/dt */
                  NULL,                   /* SidesCount_iter,            /* callback to compute each quads'
                                             faces' contributions to du/du */
                  NULL,                   /* there is no callback for the
                                             edges between quadrants */
                  NULL);

    // p4est_destroy(p4est);
    p4est = p4est1;
    // printf("Balance is DONE !!! \n");
    data->nElems = p4est1->local_num_quadrants;
    return;
    // exit(1);
}

//From Hopest, but was changed to connectivity,  not p4est
void p4_build_bcs(p8est_connectivity_t *connectivity,
                  p4est_topidx_t num_trees,
                  int32_t *bcelemmap) {
    int itree, iside;

    p8est_connectivity_t *conn = connectivity;
    P4EST_ASSERT(p4est_connectivity_is_valid(conn));
    // printf("Connectivity adress = %p\n", conn);
    P4EST_ASSERT(connectivity->num_trees == num_trees);
    p4est_connectivity_set_attr(conn, 6 * sizeof(int32_t));
    P4EST_ASSERT(p4est_connectivity_is_valid(conn));

    for (itree = 0; itree < num_trees; itree++) {
        for (iside = 0; iside < 6; iside++) {
            ((int32_t *) conn->tree_to_attr)[itree * 6 + iside] = bcelemmap[itree * 6 + iside];
        }
    }
}

//#include "save.fc"
