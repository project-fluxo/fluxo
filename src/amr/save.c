#include "connectivity.h"
#include "p4fluxo.h"


static void
SaveMeshSetElement_number(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    savemesh_inner_data_t *dataquad = (savemesh_inner_data_t *) q->p.user_data;
    p4est_t *p4est = info->p4est;
    int *iElem = (int *) user_data;
    int i, j;
    // p4est_fortran_data_t *data = (p4est_fortran_data_t *)user_data;
    (*iElem)++;
    dataquad->ElementID = (*iElem);
    for (i = 0; i < 6; i++) {
        dataquad->SidesID[i] = 0;
        dataquad->flips[i] = 100;
        for (j = 0; j < 4; j++) {
            dataquad->MortarSides[i][j] = 0;
            dataquad->nbElementID[i][j] = 0;
        }
    }
}

static void
SaveMeshSidesCounter_iter(p4est_iter_face_info_t *info, void *user_data) {
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;
    int *nSides = (int *) user_data;
    p4est_quadrant_t *quad1, *quad2, *quad;
    int *ghost_to_proc = (int *) p4est->user_pointer;
    int which_face;
    int ghost = 0;
    p4est_iter_face_side_t *side[2] = {NULL, NULL};
    sc_array_t *sides = &(info->sides);
    // int AlreadyCount = 0; // This Flag set to 1 if there is a simple side and we just one number for this side

    if (sides->elem_count == 1) {
        (*nSides)++;
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
                    // int ghostid = side[iSmallSide]->is.hanging.quadid[j];
                    // int proc=data->ghost_to_proc[ghostid];
                    // data->nMPISides_Proc[proc]++;
                } else {
                    //Not MPI Side
                    (*nSides)++;
                }
            }
            (*nSides) += 5; //Big Side + 4 Small Master Mortar

        } else { //We use only 4 small sides and don't take into account the big side

            for (j = 0; j < P8EST_HALF; j++) {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j]) {
                    // ghost++;
                } else {
                    (*nSides)++;
                }
            }

            // data->nSides += (4 - ghost);    //Small sides are MPI sides. Ghost is for another proc.
            // data->nMPISides += (4 - ghost); //Small sides  are MPI Sides
        }
    } else //Then it is a one side
    {

        if (side[0]->is.full.is_ghost || side[1]->is.full.is_ghost) {
            //One is Ghost Side
            int SideIn = 0;
            int SideGhost = 1;
            if (side[0]->is.full.is_ghost) {
                SideIn = 1;
                SideGhost = 0;
            }
            quad = side[SideIn]->is.full.quad;
            int id = side[SideGhost]->is.full.quadid;
            int proc = ghost_to_proc[id];
            if (proc < p4est->mpirank)
                (*nSides)++; // This side is  mine
        } else //Not Ghost
        {
            // quad = side[i]->is.full.quad;
            (*nSides)++;
            // nSides++;
            // data->nInnerSides++;
        }
    }
}

static void
SaveMeshRenumerateSideCounter_iter(p4est_iter_face_info_t *info, void *user_data) {
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;
    p8est_connectivity_t *conn = p4est->connectivity;
    int *nSides = (int *) user_data;
    p4est_quadrant_t *quad1, *quad2, *quad;
    int *ghost_to_proc = (int *) p4est->user_pointer;
    int which_face;
    int ghost = 0;
    int direction;
    int orientation = info->orientation;
    p4est_iter_face_side_t *side[2] = {NULL, NULL};
    sc_array_t *sides = &(info->sides);
    savemesh_inner_data_t *dataquad, *dataBigquad;

    // int AlreadyCount = 0; // This Flag set to 1 if there is a simple side and we just one number for this side

    // return; //DEBUG
    if (sides->elem_count == 1) {
        (*nSides)++;
        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        direction = side[0]->face / 2; /* 0 == x, 1 == y, 2 == z */
        quad = side[0]->is.full.quad;
        which_face = side[0]->face;
        int treeid = side[0]->treeid;
        int face = side[0]->face;
        dataquad = (savemesh_inner_data_t *) quad->p.user_data;
        int HFace = P2H_side[which_face] - 1;
        dataquad->SidesID[HFace] = (*nSides);
        dataquad->flips[HFace] = 0;
        dataquad->nbElementID[HFace][0] = ((int32_t *) conn->tree_to_attr)[treeid * 6 + face];
        return;
    }
    // return; //DEBUG
    // side[i]->is.hanging.quad[j]->p.user_data;
    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    int minus = side[0]->face % 2;
    if (side[0]->is_hanging || side[1]->is_hanging) {
        // return; //DEBUG
        //One Side is Mortar
        iBigSide = side[0]->is_hanging == 0 ? 0 : 1;
        iSmallSide = side[0]->is_hanging != 0 ? 0 : 1;
        int BigFace = side[iBigSide]->face;
        int SmallFace = side[iSmallSide]->face;
        int BigFaceH = P2H_side[BigFace] - 1;
        int SmallFaceH = P2H_side[SmallFace] - 1;
        int flip = GetHFlip(BigFace, SmallFace, orientation);
        if (side[iBigSide]->is.full.is_ghost == 0) //Big side is not MPI
        {
            // return; //DEBUG
            (*nSides)++;                     //Big Side
            quad = side[iBigSide]->is.full.quad;
            dataBigquad = (savemesh_inner_data_t *) quad->p.user_data;
            dataBigquad->SidesID[BigFaceH] = (*nSides);
            dataBigquad->flips[BigFaceH] = flip + 10 * (SmallFaceH + 1); //??????

            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {

                int jIndex = 0; //f(j); //Convert index from P4 to H
                int Pside = BigFace;
                int PnbSide = SmallFace;
                int PFlip = orientation;
                // if (KKK==0){
                jIndex = GetHMortar(j, Pside, PnbSide, PFlip) - 1; // HOPR index

                dataBigquad->MortarSides[BigFaceH][jIndex] = ++(*nSides);
            }


            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j]) {
                    // ghost++;
                    // int ghostid = side[iSmallSide]->is.hanging.quadid[j];
                    // int proc=data->ghost_to_proc[ghostid];
                    // data->nMPISides_Proc[proc]++;
                } else {
                    //Not MPI Side
                    // (*nSides)++;
                    int jIndex = 0; //f(j); //Convert index from P4 to H
                    int Pside = BigFace;
                    int PnbSide = SmallFace;
                    int PFlip = orientation;
                    // if (KKK==0){
                    jIndex = GetHMortar(j, Pside, PnbSide, PFlip) - 1; // HOPR index
                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = (savemesh_inner_data_t *) quad->p.user_data;
                    // dataquad->SidesID[SmallFaceH] = -(*nSides);
                    dataquad->SidesID[SmallFaceH] = -1 * (dataBigquad->MortarSides[BigFaceH][jIndex]);
                    dataquad->flips[SmallFaceH] = -flip - 10 * (BigFaceH + 1); //??????
                }
            }
            // if (ghost == 0) //There is no small mpi sides. Add Mortar inner side
            // {
            //     nSides++;
            // }
            // else //One or more small mpi sides. Add Big side as MPIMortar
            // {
            //     data->nMortarMPISides++;
            // }
            // data->nMPISides += ghost;       //Increase number of small MPISides
            // data->nSides += 4;              //Total nmber of sides
            // data->nSides++;                 //Add also a MortarSide to nSides
            // data->nInnerSides += 4 - ghost; //Number of inner sides
        } else { //We use only 4 small sides and don't take into account the big side
            // return; //DEBUG


            for (j = 0; j < P8EST_HALF; j++) {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j]) {
                    // ghost++;
                } else {
                    (*nSides)++;
                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = (savemesh_inner_data_t *) quad->p.user_data;
                    // dataquad->SidesID[SmallFaceH] = -(*nSides);
                    dataquad->flips[SmallFaceH] = -flip - 10 * (BigFaceH + 1); //??????
                }
            }


            // data->nSides += (4 - ghost);    //Small sides are MPI sides. Ghost is for another proc.
            // data->nMPISides += (4 - ghost); //Small sides  are MPI Sides
        }
    } else //Then it is a one side
    {

        if (side[0]->is.full.is_ghost || side[1]->is.full.is_ghost) {
            // return; //DEBUG
            //One is Ghost Side
            int SideIn = 0;
            int SideGhost = 1;
            if (side[0]->is.full.is_ghost) {
                SideIn = 1;
                SideGhost = 0;
            }
            quad = side[SideIn]->is.full.quad;
            int faceIn = side[SideIn]->face;
            int faceGhost = side[SideGhost]->face;

            int faceInH = P2H_side[faceIn] - 1;
            int faceGhostH = P2H_side[faceGhost] - 1;

            int id = side[SideGhost]->is.full.quadid;
            int proc = ghost_to_proc[id];
            if (proc < p4est->mpirank) {
                (*nSides)++; // This side is  mine
                dataquad = (savemesh_inner_data_t *) quad->p.user_data;
                dataquad->SidesID[faceInH] = (*nSides);
                dataquad->flips[faceInH] = GetHFlip(faceIn, faceGhost, orientation) + 10 * (faceGhostH + 1);
                // if (p4est->mpirank == 2)
                // {
                //     printf("dataquad->flips[%d] = %d \n", faceInH, dataquad->flips[faceInH]);
                // }
            } else {
                // This side is not mine
                dataquad = (savemesh_inner_data_t *) quad->p.user_data;
                dataquad->SidesID[faceInH] = 0;
                dataquad->flips[faceInH] = -GetHFlip(faceIn, faceGhost, orientation) - 10 * (faceGhostH + 1);
                // if (p4est->mpirank == 2)
                // {
                //     printf("dataquad->flips[%d] = %d \n", faceInH, dataquad->flips[faceInH]);
                // }
                ; //Will be Set up later
            }
        } else //Not Ghost
        {
            // quad = side[i]->is.full.quad;
            // return; //DEBUG
            (*nSides)++;
            quad = side[0]->is.full.quad;
            quad1 = side[1]->is.full.quad;
            int face = side[0]->face;
            int face1 = side[1]->face;

            int faceH = P2H_side[face] - 1;
            int face1H = P2H_side[face1] - 1;

            dataquad = (savemesh_inner_data_t *) quad->p.user_data;
            int flip = GetHFlip(face, face1, orientation);

            dataquad->SidesID[faceH] = (*nSides);
            if (minus) flip = -flip;
            if (flip > 0)
                dataquad->flips[faceH] = flip + 10 * (face1H + 1);
            else
                dataquad->flips[faceH] = flip - 10 * (face1H + 1);
            // printf("faceH = flip %d, face1H = %d , flip = %d\n",dataquad->flips[faceH], face1H, flip);

            dataquad = (savemesh_inner_data_t *) quad1->p.user_data;
            dataquad->SidesID[face1H] = (*nSides);
            if (flip > 0)
                dataquad->flips[face1H] = -flip - 10 * (faceH + 1);
            else
                dataquad->flips[face1H] = -flip + 10 * (faceH + 1);
            // printf("face1H flip = %d, \n", dataquad->flips[face1H]);
        }
    }
}

static void
SaveMeshGhostCounter_iter(p4est_iter_face_info_t *info, void *user_data) {
    int i, j, iBigSide = 0, iSmallSide = 0;
    p4est_t *p4est = info->p4est;
    p4est_savef_data_t *p4est_savef_data = (p4est_savef_data_t *) user_data;
    p4est_quadrant_t *quad1, *quad2, *quad;
    savemesh_inner_data_t *ghost_data = (savemesh_inner_data_t *) p4est->user_pointer;
    int which_face;
    int ghost = 0;
    int direction;
    int orientation = info->orientation;
    p4est_iter_face_side_t *side[2] = {NULL, NULL};
    sc_array_t *sides = &(info->sides);
    savemesh_inner_data_t *dataquad;

    if (sides->elem_count == 1) {
        // (*nSides)++;
        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        // direction = side[0]->face / 2; /* 0 == x, 1 == y, 2 == z */
        quad = side[0]->is.full.quad;
        which_face = side[0]->face;
        dataquad = (savemesh_inner_data_t *) quad->p.user_data;


        return;
    }

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);

    if (side[0]->is_hanging || side[1]->is_hanging) {
        ////We deal only with non hanging sides
        // //One Side is Mortar
        iBigSide = side[0]->is_hanging == 0 ? 0 : 1;
        iSmallSide = side[0]->is_hanging != 0 ? 0 : 1;
        int BigFace = side[iBigSide]->face;
        int SmallFace = side[iSmallSide]->face;
        int BigFaceH = P2H_side[BigFace] - 1;
        int SmallFaceH = P2H_side[SmallFace] - 1;

        if (side[iBigSide]->is.full.is_ghost == 0) //Big side is not MPI
        {

            p4est_quadrant_t *BigQuad = side[iBigSide]->is.full.quad;
            savemesh_inner_data_t *BigQuadData = BigQuad->p.user_data;

            for (j = 0; j < P8EST_HALF; j++) //Check if the other sides MPI
            {
                int jIndex = 0; //f(j); //Convert index from P4 to H
                int Pside = BigFace;
                int PnbSide = SmallFace;
                int PFlip = orientation;
                // if (KKK==0){
                jIndex = GetHMortar(j, Pside, PnbSide, PFlip) - 1;
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j]) {
                    // ghost++;
                    int ghostid = side[iSmallSide]->is.hanging.quadid[j];


                    // data->nMPISides_Proc[proc]++;
                    BigQuadData->nbElementID[BigFaceH][jIndex] = ghost_data[ghostid].ElementID;
                } else {
                    //Not MPI Side


                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = (savemesh_inner_data_t *) quad->p.user_data;
                    dataquad->nbElementID[SmallFaceH][0] = BigQuadData->ElementID;

                    BigQuadData->nbElementID[BigFaceH][jIndex] = dataquad->ElementID;

                }
            }

        } else { //We use only 4 small sides and don't take into account the big side
            int BigGhostid = side[iBigSide]->is.full.quadid;


            for (j = 0; j < P8EST_HALF; j++) {
                //quad=side[i]->is.hanging.quad[j];
                if (side[iSmallSide]->is.hanging.is_ghost[j]) {
                    // ghost++;
                } else {

                    int jIndex = 0; //f(j); //Convert index from P4 to H
                    int Pside = BigFace;
                    int PnbSide = SmallFace;
                    int PFlip = orientation;

                    jIndex = GetHMortar(j, Pside, PnbSide, PFlip) - 1;

                    quad = side[iSmallSide]->is.hanging.quad[j];
                    dataquad = (savemesh_inner_data_t *) quad->p.user_data;

                    dataquad->nbElementID[SmallFaceH][0] = ghost_data[BigGhostid].ElementID;
                    dataquad->SidesID[SmallFaceH] = -1 * (ghost_data[BigGhostid].MortarSides[BigFaceH][jIndex]);

                }
            }

        }
    } else //Then it is a one side
    {

        if (side[0]->is.full.is_ghost || side[1]->is.full.is_ghost) {
            //One is Ghost Side
            int SideIn = 0;
            int SideGhost = 1;
            if (side[0]->is.full.is_ghost) {
                SideIn = 1;
                SideGhost = 0;
            }
            p4est_quadrant_t *quadIn;
            quadIn = side[SideIn]->is.full.quad;
            // p4est_quadrant_t *quadIn = side[SideIn]->is.full.quad;
            int faceIn = side[SideIn]->face;
            int faceGhost = side[SideGhost]->face;

            int faceInH = P2H_side[faceIn] - 1;
            int faceGhostH = P2H_side[faceGhost] - 1;

            int GhostId = side[SideGhost]->is.full.quadid;
            savemesh_inner_data_t *DataIn = (savemesh_inner_data_t *) quadIn->p.user_data;
            savemesh_inner_data_t *DataGhost = &ghost_data[GhostId];
            if (DataIn->SidesID[faceInH] == 0) {
                DataIn->SidesID[faceInH] = DataGhost->SidesID[faceGhostH];
                DataIn->flips[faceInH] = -DataGhost->flips[faceGhostH];


            }
            DataIn->nbElementID[faceInH][0] = DataGhost->ElementID;

        } else //Not Ghost
        {
            quad1 = side[0]->is.full.quad;
            quad2 = side[1]->is.full.quad;
            savemesh_inner_data_t *Data1 = (savemesh_inner_data_t *) quad1->p.user_data;
            savemesh_inner_data_t *Data2 = (savemesh_inner_data_t *) quad2->p.user_data;
            int face1 = side[0]->face;
            int face2 = side[1]->face;

            int face1H = P2H_side[face1] - 1;
            int face2H = P2H_side[face2] - 1;

            Data1->nbElementID[face1H][0] = Data2->ElementID;
            Data2->nbElementID[face2H][0] = Data1->ElementID;

            return;

        }
    }
}


//////////////////////////////////////////////////////////
// 
// The Array ElemInfo is defined here. 
//
///////////////////////////////////////////
static void
SaveMeshFillArray_iter(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    savemesh_inner_data_t *dataquad = (savemesh_inner_data_t *) q->p.user_data;
    p4est_t *p4est = info->p4est;
    p4est_savef_data_t *data = (p4est_savef_data_t *) user_data;
    int *iElem = p4est->user_pointer;
    int *ElemInfo = data->ElemInfo;
    int *SideInfo = data->SideInfo;
    int iLocalSide = 0;
    int i, j = 0, nSidesHere = 0, A;
    static int SIndex = 0;
    int EIndex;
    // p4est->mpirank
    // p4est->global_first_quadrant[p4est->mpirank];
    int Offset = p4est->global_first_quadrant[p4est->mpirank];
    EIndex = dataquad->ElementID - Offset - 1;

    (*iElem)++;
    data->nLocalSides++;
    int ElemID = data->nLocalSides;
    // SideInfo[(iSide - 1) * 5 + (i - 1)] = ;
    // 1. Count the number of the sides from 6 to (4+1)*6 *every Big Side is 5 side (1 Big Side + 4 SmallMortar master)
    for (iLocalSide = 0; iLocalSide < 6; iLocalSide++) {
        if (dataquad->nbElementID[iLocalSide][1] != 0) {
            nSidesHere += 5;
        } else {
            nSidesHere++;
        }
    }
    A = data->nGlobalSides;
    i = 1;
    ElemInfo[(EIndex) * 6 + (i - 1)] = 208; //Element Type = HexaHedron non-linear
    i = 2;
    ElemInfo[(EIndex) * 6 + (i - 1)] = 1; // Zone =1

    i = 3;
    ElemInfo[(EIndex) * 6 + (i - 1)] = SIndex; // offsetIndSIDE

    i = 4;
    SIndex += nSidesHere;
    ElemInfo[(EIndex) * 6 + (i - 1)] = SIndex; // lastIndSIDE
    i = 5;
    // if (p4est->mpirank == 0) printf("EIndex = %d , dataquad->ElementID = %d \n", EIndex, dataquad->ElementID);
    // printf(" (EIndex - 1) * 6 + (i - 1) = %d , i = %d \n", (EIndex - 1) * 6 + (i - 1), i);
    // fflush(stdout);
    // ElemInfo[(8 - 1) * 6 + (6 - 1)] = 0;
    ElemInfo[(EIndex) * 6 + (i - 1)] = 0; // offsetIndNODE is filled up in FLUXO
    i = 6;
    ElemInfo[(EIndex) * 6 + (i - 1)] = 0; // lastIndNODE is filled up in FLUXO
    return;


}

//////////////////////////////////////////////////////////
//
// The Array SideInfo is defined here.
//
///////////////////////////////////////////
static void
SaveMeshFillArray2_iter(p4est_iter_volume_info_t *info, void *user_data) {
    p4est_quadrant_t *q = info->quad;
    savemesh_inner_data_t *dataquad = (savemesh_inner_data_t *) q->p.user_data;
    p4est_t *p4est = info->p4est;
    p4est_savef_data_t *data = (p4est_savef_data_t *) user_data;
    // int *iElem = p4est->user_pointer;
    int *ElemInfo = data->ElemInfo;
    int *SideInfo = data->SideInfo;
    int iLocalSide = 0;
    int i, j = 0, nSidesHere = 0, A;
    static int SIndex = 0;
    int EIndex;
    // p4est->mpirank
    // p4est->global_first_quadrant[p4est->mpirank];
    int Offset = p4est->global_first_quadrant[p4est->mpirank];
    EIndex = dataquad->ElementID - Offset - 1;

    for (iLocalSide = 0; iLocalSide < 6; iLocalSide++) // This is a HOPR Numeration
    {

        if (dataquad->nbElementID[iLocalSide][1] != 0) {
            // Mortar => 5 Sides.
            // printf("ElementID = %d , Mortar Side %d \n", dataquad->ElementID, iLocalSide);
            // First Big Mortar Side
            i = 1;
            SideInfo[(SIndex) * 5 + (i - 1)] = 4; //Side Type = HexaHedron non-linear
            i = 2;
            // if (dataquad->SidesID[iLocalSide] <0) Mortar Big Side
            SideInfo[(SIndex) * 5 + (i - 1)] = fabs(
                    dataquad->SidesID[iLocalSide]); // GlobalSideID positive (master) and negative (slave)

            i = 3;
            SideInfo[(SIndex) * 5 + (i - 1)] = -1; // Mortar Type * -1. Here is always 1 Type.
            i = 4;
            SideInfo[(SIndex) * 5 + (i - 1)] = 0; // 10*nbLocSide + FLIP
            i = 5;
            SideInfo[(SIndex) * 5 + (i - 1)] = 0; // BCID

            SIndex++;                           // 5 Times
            // Then 4 small master mortar sides
            /// First
            int iM = 0;
            for (iM = 0; iM < P8EST_HALF; iM++)// This is a HOPR Numeration - 1
            {

                i = 1;
                SideInfo[(SIndex) * 5 + (i - 1)] = 104; //Side Type = HexaHedron non-linear
                i = 2;
                // MortarSides[6][4]
                SideInfo[(SIndex) * 5 + (i - 1)] = fabs(
                        dataquad->MortarSides[iLocalSide][iM]); // GlobalSideID positive (master) and negative (slave)

                i = 3;
                SideInfo[(SIndex) * 5 + (i - 1)] = dataquad->nbElementID[iLocalSide][iM]; // Mortar Type * -1
                i = 4;
                SideInfo[(SIndex) * 5 + (i - 1)] = 0; // 10*nbLocSide + FLIP
                i = 5;
                SideInfo[(SIndex) * 5 + (i - 1)] = 0; // BCID

                SIndex++;                           //  4 Times
            }
        } else { // Non hanging side

            if (dataquad->flips[iLocalSide] == 0) //Boundary Side
            {
                i = 1;
                SideInfo[(SIndex) * 5 + (i - 1)] = 104; //Side Type = HexaHedron non-linear
                i = 2;
                // MortarSides[6][4]
                SideInfo[(SIndex) * 5 + (i -
                                         1)] = dataquad->SidesID[iLocalSide]; // GlobalSideID positive (master) and negative (slave)

                i = 3;
                SideInfo[(SIndex) * 5 + (i - 1)] = 0;//dataquad->nbElementID[iLocalSide][iM]; // Mortar Type * -1
                i = 4;
                SideInfo[(SIndex) * 5 + (i - 1)] = 0; // 10*nbLocSide + FLIP
                i = 5;
                SideInfo[(SIndex) * 5 + (i - 1)] = dataquad->nbElementID[iLocalSide][0]; // BCID

                SIndex++; //  Next Side
                continue;
            };

            if (dataquad->flips[iLocalSide] > 0 && dataquad->SidesID[iLocalSide] < 0)
                exit(1); // Error
            // Not mortar. => 1 Side/
            if (dataquad->flips[iLocalSide] < 0 && dataquad->SidesID[iLocalSide] < 0) //Small slave Mortar Side
            {
                i = 1;
                SideInfo[(SIndex) * 5 + (i - 1)] = -104; //Side Type = HexaHedron non-linear
                i = 4;
                SideInfo[(SIndex) * 5 + (i - 1)] = ((int) fabs(dataquad->flips[iLocalSide])) %
                                                   10; //It is a Flip!     //fabs(dataquad->flips[iLocalSide]); // 10*nbLocSide + FLIP dataquad->
            } else {
                i = 1; //Non Mortar Side
                SideInfo[(SIndex) * 5 + (i - 1)] = 4;
                i = 4;
                SideInfo[(SIndex) * 5 + (i - 1)] = fabs(
                        dataquad->flips[iLocalSide]); //fabs(dataquad->flips[iLocalSide]);

            }
            i = 2;
            // MortarSides[6][4]
            if (dataquad->flips[iLocalSide] < 0)                                           // Slave
                SideInfo[(SIndex) * 5 + (i - 1)] =
                        -1 * fabs(dataquad->SidesID[iLocalSide]); // GlobalSideID positive (master) and negative (slave)
            else
                SideInfo[(SIndex) * 5 + (i - 1)] = fabs(
                        dataquad->SidesID[iLocalSide]); // GlobalSideID positive (master) and negative (slave)

            i = 3;
            SideInfo[(SIndex) * 5 + (i - 1)] = dataquad->nbElementID[iLocalSide][0]; // Mortar Type * -1
            //Non Mortar Side

            i = 5;
            SideInfo[(SIndex) * 5 + (i - 1)] = 0; // BCID

            SIndex++; // 1 Times for each Side
        }


    }
}

void count_parallel_sides(p4est_t *p4est) {
    int mpiret;
    p4est_savef_data_t *p4est_savef_data = (p4est_savef_data_t *) p4est->user_pointer;

    int local = p4est_savef_data->nLocalSides;
    // p4est_gloidx_t *global_first_quadrant = p4est->global_first_quadrant;
    int i;
    const int num_procs = p4est->mpisize;

    p4est_savef_data->OffsetSideMPI[0] = 0;
    // global_first_quadrant[0] = 0;
    mpiret = sc_MPI_Allgather(&local, 1, MPI_INT,
                              p4est_savef_data->OffsetSideMPI + 1, 1, MPI_INT,
                              p4est->mpicomm);
    // mpiret = sc_MPI_Allgather(qlocal, 1, P4EST_MPI_GLOIDX,
    //                           global_first_quadrant + 1, 1, P4EST_MPI_GLOIDX,
    //                           p4est->mpicomm);
    SC_CHECK_MPI(mpiret);

    for (i = 0; i < num_procs; ++i) {
        p4est_savef_data->OffsetSideMPI[i + 1] += p4est_savef_data->OffsetSideMPI[i];
    }
    p4est_savef_data->nGlobalSides = p4est_savef_data->OffsetSideMPI[num_procs];
    // p4est->global_num_quadrants = global_first_quadrant[num_procs];
}

void count_parallel_index_sides(p4est_t *p4est) {
    int mpiret;
    p4est_savef_data_t *p4est_savef_data = (p4est_savef_data_t *) p4est->user_pointer;

    int local = p4est_savef_data->nSidesArrIndex;
    // p4est_gloidx_t *global_first_quadrant = p4est->global_first_quadrant;
    int i;
    const int num_procs = p4est->mpisize;

    p4est_savef_data->OffsetSideArrIndexMPI[0] = 0;
    // global_first_quadrant[0] = 0;
    mpiret = sc_MPI_Allgather(&local, 1, MPI_INT,
                              p4est_savef_data->OffsetSideArrIndexMPI + 1, 1, MPI_INT,
                              p4est->mpicomm);
    // mpiret = sc_MPI_Allgather(qlocal, 1, P4EST_MPI_GLOIDX,
    //                           global_first_quadrant + 1, 1, P4EST_MPI_GLOIDX,
    //                           p4est->mpicomm);
    SC_CHECK_MPI(mpiret);

    for (i = 0; i < num_procs; ++i) {
        p4est_savef_data->OffsetSideArrIndexMPI[i + 1] += p4est_savef_data->OffsetSideArrIndexMPI[i];
    }
    p4est_savef_data->nGlobalSides = p4est_savef_data->OffsetSideArrIndexMPI[num_procs];
    // p4est->global_num_quadrants = global_first_quadrant[num_procs];
}

p4est_savef_data_t *save_mesh(p4est_t *p4est1) {
    p8est_ghost_t *ghost = NULL;
    //  1. Count Local Sides (just number)
    p8est_t *p4est;
    p4est_savef_data_t *p4est_savef_data = (p4est_savef_data_t *) malloc(sizeof(p4est_savef_data_t));
    int nLocalSides = 0;
    // Iterator other sides
    int rank = 0;
    // int nNBProcs = 0;
    int *ghost_to_proc = NULL;
    // if (p4est->mpisize > 1)

    //Set New data
    p4est = p8est_copy(p4est1, 0);
    p4est->connectivity = p4est1->connectivity; // 0 - no data to copy

    p4est_reset_data(p4est, sizeof(savemesh_inner_data_t), NULL, NULL);
    ghost = p8est_ghost_new(p4est, P8EST_CONNECT_FACE);

    savemesh_inner_data_t *ghost_data = (savemesh_inner_data_t *) malloc(
            ghost->ghosts.elem_count * sizeof(savemesh_inner_data_t));

    int num_ghost = (p4est_locidx_t) ghost->ghosts.elem_count;
    ghost_to_proc = (int *) malloc(num_ghost * sizeof(int));

    int jl;
    for (jl = 0; jl < num_ghost; ++jl) {

        while (ghost->proc_offsets[rank + 1] <= jl) {
            ++rank;
            P4EST_ASSERT(rank < p4est->mpisize);
        }
        ghost_to_proc[jl] = rank;
    }

    p4est->user_pointer = ghost_to_proc;
    p4est_iterate(p4est,                     /* the forest */
                  ghost,                     /* the ghost layer May be LAter!!! */
                  (void *) &nLocalSides,      /* the synchronized ghost data */
                  NULL,                      //ElementCounter_iter,        /* callback to compute each quad's
            //                         interior contribution to du/dt */
                  SaveMeshSidesCounter_iter, /* callback to compute each quads'
                                        faces' contributions to du/du */
                  NULL,                      /* there is no callback for the
                                        edges between quadrants */
                  NULL);
    p4est->user_pointer = NULL;
    p4est_savef_data->nLocalSides = nLocalSides;
    //  2. Count Global Number of Sides (Array[mpisize])
    p4est_savef_data->OffsetSideMPI = (int *) malloc(p4est->mpisize * sizeof(int));
    p4est->user_pointer = p4est_savef_data;
    count_parallel_sides(p4est);
    p4est->user_pointer = NULL;
    int iElem = p4est->global_first_quadrant[p4est->mpirank];

    p4est_iterate(p4est,                     /* the forest */
                  ghost,                     /* the ghost layer May be LAter!!! */
                  (void *) &iElem,            /* the synchronized ghost data */
                  SaveMeshSetElement_number,    //ElementCounter_iter,        /* callback to compute each quad's
            //                         interior contribution to du/dt */
                  NULL, /* callback to compute each quads'
                                        faces' contributions to du/du */
                  NULL,                      /* there is no callback for the
                                        edges between quadrants */
                  NULL);
    // pfree(ghost_to_proc);
    //  3. Renumerate sides
    // 3a. Firstly set up local number in acoordance with OffsetSide

    nLocalSides = p4est_savef_data->OffsetSideMPI[p4est->mpirank];
    p4est->user_pointer = ghost_to_proc;
    p4est_iterate(p4est,                     /* the forest */
                  ghost,                     /* the ghost layer May be LAter!!! */
                  (void *) &nLocalSides,      /* the synchronized ghost data */
                  NULL,                      //ElementCounter_iter,        /* callback to compute each quad's
            //                         interior contribution to du/dt */
                  SaveMeshRenumerateSideCounter_iter, /* callback to compute each quads'
                                        faces' contributions to du/du */
                  NULL,                      /* there is no callback for the
                                        edges between quadrants */
                  NULL);
    p4est->user_pointer = NULL;
    // 3b. Set CALL GHOST EXCHANGE
    p4est_ghost_exchange_data(p4est, ghost, ghost_data);
    // 3c. Set Sides, that connected with other processors.
    // Create new data???
    p4est->user_pointer = ghost_data;
    p4est_iterate(p4est,                              /* the forest */
                  ghost,                              /* the ghost layer May be LAter!!! */
                  (void *) p4est_savef_data,           /* the synchronized ghost data */
                  NULL,                               //ElementCounter_iter,        /* callback to compute each quad's
            //                         interior contribution to du/dt */
                  SaveMeshGhostCounter_iter, /* callback to compute each quads'
                                        faces' contributions to du/du */
                  NULL,                               /* there is no callback for the
                                        edges between quadrants */
                  NULL);

    p4est->user_pointer = NULL;
    int OffsetElem = p4est->global_first_quadrant[p4est->mpirank];
    iElem = OffsetElem; //First set ElementID
    int nSave = p4est_savef_data->nLocalSides;
    int nGlobalSideSave = p4est_savef_data->nGlobalSides;

    p4est_savef_data->nGlobalSides = p4est_savef_data->OffsetSideMPI[p4est->mpirank];
    p4est_savef_data->nLocalSides = 0;//It is used as a counter in the array. So the value will be restored further.
    p4est->user_pointer = &iElem; // Number of the Element with offset.
    p4est_savef_data->ElemInfo = (int *) malloc(6 * p4est->local_num_quadrants * sizeof(int));

    p4est_iterate(p4est,                    /* the forest */
                  ghost,                    /* the ghost layer May be LAter!!! */
                  (void *) p4est_savef_data, /* the synchronized ghost data */
                  SaveMeshFillArray_iter,   /* callback to compute each quads'*/
                  NULL,                     //ElementCounter_iter,        /* callback to compute each quad's
            //                         interior contribution to du/dt */
            // faces' contributions to du/du */
                  NULL,                     /* there is no callback for the
                                    edges between quadrants */
                  NULL);


    p4est_savef_data->OffsetSideArrIndexMPI = (int *) malloc(p4est->mpisize * sizeof(int));
    p4est_savef_data->nSidesArrIndex = p4est_savef_data->ElemInfo[(p4est->local_num_quadrants - 1) * 6 + (4 - 1)];

    p4est->user_pointer = (void *) p4est_savef_data;
    count_parallel_index_sides(p4est);
    p4est_savef_data->SideInfo = (int *) malloc(5 * p4est_savef_data->nSidesArrIndex * sizeof(int));
    p4est_iterate(p4est,                    /* the forest */
                  NULL,                    /* the ghost layer May be LAter!!! */
                  (void *) p4est_savef_data, /* the synchronized ghost data */
                  SaveMeshFillArray2_iter,   /* callback to compute each quads'*/
                  NULL,                     //ElementCounter_iter,        /* callback to compute each quad's
            //                         interior contribution to du/dt */
            // faces' contributions to du/du */
                  NULL,                     /* there is no callback for the
                                    edges between quadrants */
                  NULL);

    p4est->user_pointer = NULL;


    // Set SideInfo Array
    // Add in cycle the the Index offset 
    // Test the Mortar case!
    // Pass to FLUXO
    // Save HDF5 Mesh
    // Add NODE Index
    // Check the Unique... Stuff
    // // p8est_ghost_destroy(ghost); // Should I create a new GHOST???
    // ghost = p8est_ghost_new(p4est, P8EST_CONNECT_FACE);
    p4est_savef_data->nLocalSides = nSave;
    p4est_savef_data->nGlobalSides = nGlobalSideSave;
    free(ghost_data);
    ghost_data = NULL;


    pfree(ghost_to_proc);
    p8est_ghost_destroy(ghost);
    p4est_destroy(p4est);

    return p4est_savef_data;
}

void
savef_data_destroy(p4est_savef_data_t *p4est_savef_data) {
    pfree(p4est_savef_data->OffsetSideMPI);
    pfree(p4est_savef_data->OffsetSideArrIndexMPI);
    pfree(p4est_savef_data->ElemInfo);
    pfree(p4est_savef_data->SideInfo);
    pfree(p4est_savef_data);
};

void save_p4est(p4est_t *p4est, char in[]) {
    printf("save p4est!! %d \n", p4est->local_num_quadrants);
    p4est_save_ext(in, p4est, 0, 0);

    return;
};

p4est_t *load_p4est(int mpicomm1, char in[]) {

    p8est_connectivity_t **conn;
    p4est_t *p4est;
    static sc_MPI_Comm mpicomm;
    mpicomm = MPI_Comm_f2c(mpicomm1);


    p4est = p8est_load_ext(in,      //const char *filename,
                           mpicomm, //sc_MPI_Comm mpicomm,
                           0,       //size_t data_size,
                           0,       //int load_data,
                           1,       //int autopartition,
                           1,       //int broadcasthead,
                           NULL,    //void *user_pointer,
                           conn);   //p8est_connectivity_t * *connectivity);
    p4est->connectivity = *conn;
    printf("Load p4est!! %d \n", p4est->connectivity->num_trees);

    return p4est;
};

p8est_connectivity_t *GetConnectivity(p4est_t *p4est) {

    printf("p4est->connectivity = %p \n", p4est->connectivity);
    return p4est->connectivity;
};