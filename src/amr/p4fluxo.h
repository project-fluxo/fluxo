#ifndef P4FLUXO_H
#define P4FLUXO_H

// #include "stdio.h"
// #include "malloc.h"
#if defined(__MACH__)
#include <stdlib.h>
#else
#include <malloc.h>
#endif

#include "mpi.h"
#include <sc_io.h>
#include "p8est.h"
#include <p8est_extended.h>
#include <p4est_to_p8est.h>

#define nullptr ((void*)0)
#define NON_OPTIMIZED
#define MORTAR_SIDE_WEIGHT 1
static sc_MPI_Comm mpicomm;


void pfree(void *A);

// The Data to save Mesh in HOPR file
typedef struct p4est_savef_data {
    int nGlobalSides;
    int nLocalSides;
    int nSidesArrIndex;
    int *OffsetSideMPI;
    int *OffsetSideArrIndexMPI;
    int *ElemInfo; // ElemInfo in HDF5
    int *SideInfo; // SideInfo in HDF5

} p4est_savef_data_t;

typedef struct savemesh_inner_data {

    // int OldElementID[8];
    int ElementID;
    int SidesID[6];
    int8_t flips[6];
    int MortarSides[6][4];
    int nbElementID[6][4];
    // int8_t CoarsedRefined;
} savemesh_inner_data_t;

typedef struct p4est_inner_data {

    int OldElementID[8];
    int ElementID;
    int SidesID[6];
    uint8_t flips[6];
    uint8_t weight;
#ifndef NON_OPTIMIZED
    int OldSidesID[6];
    int8_t SidesRatio[6];// 0 = Sides the same. -1 = Side is smaller as the NB. +1 = Side is Bigger as NB
    int8_t IsChanged;
#endif
} p4est_inner_data_t;

typedef struct p4est_balance_data {

    int nVar;
    int PP_N;
    int nElems;
    int DataSize;
    void *DataSetU;
    void *DataSetElem_xGP; //for                                                      each process and 1 beyond * /
    // int8_t CoarsedRefined;
} p4est_balance_data_t;


typedef struct p4est_balance_datav2 {
    int DataSize; //DataSize of 1 Elem
    int GPSize; //DataSize of 1 ElemxGP
    int nElems;// Number of new elements
    p4est_gloidx_t *src_gfq; // For transfer data. for p4est only
    void *Elem_xGP_new; // new Array ElemxGP
    void *Elem_xGP_old; // old Array ElemxGP
    void *U_new; // new U
    void *U_old; // old U 
    void *ElemWasCoarsened_new; // Was the elem coarsened?
    void *ElemWasCoarsened_old; // Was the elem coarsened?
} p4est_balance_datav2_t;

typedef struct p4est_mpi_data {

    int mpisize;
    int mpirank;
    p4est_locidx_t local_num_quad;
    p4est_gloidx_t global_num_quad;
    p4est_gloidx_t *offsetMPI; // size [mpisize+1] p4est_gloidx_t     *global_first_quadrant; /**< first global quadrant index
    //for                                                      each process and 1 beyond * /
    // int8_t CoarsedRefined;
} p4est_mpi_data_t;

typedef struct p4est_fortran_data {
    int nGlobalElems;
    int nElems;
    int nSides;
    int nBCSides;
    int nMortarInnerSides;
    int nInnerSides;
    int nMPISides;
    //int nMPISides_MINE;//int nMPISides_YOUR;
    int nMortarMPISides;
    //Pointers to Arrays ElemToSide, SideToElem usw.
    void *EtSPtr; // !
    void *StEPtr; // !
    void *MTPtr; // ! //Mortar Type Pointer
    void *MIPtr; // ! //Mortar Info Pointer
    void *ChngElementPtr; //! 
    void *ChngSidePtr;// !
    int *nNbProc; // ! //ALLOCATE(NbProc(nNbProcs),nMPISides_Proc(1:nNbProcs))
    int *nMPISides_Proc; // ! ///
    int nNBProcs;
    int *nMPISides_MINE_Proc; //!
    int *nMPISides_YOUR_Proc; //! 
    int nMPISides_YOUR;
    int nMPISides_MINE;
    int *offsetMPISides_MINE;// !
    int *offsetMPISides_YOUR; // !
    int *BCs;
    p8est_ghost_t *ghost_ptr;
    p4est_inner_data_t *ghost_data_ptr;
    int *ghost_to_proc;
} p4est_fortran_data_t;

//Auxilary data for Set MPI Sides
typedef struct p4est_SetMPISide_data {
    int nNBProcs;
    int *nNbProc; //ALLOCATE(NbProc(nNbProcs),nMPISides_Proc(1:nNbProcs))
    int *nMPISides_Proc;
    int *nMPISidesCount; //Size of mpisize
    int *ghost_to_proc;
    int *nMPISides_MINE_Proc;
    int *nMPISides_YOUR_Proc;
    int nMPISides_YOUR;
    int nMPISides_MINE;
    int *offsetMPISides_MINE;
    int *offsetMPISides_YOUR;

    // int nMPISides;
    //int nMPISides_MINE;//int nMPISides_YOUR;
    // int nMortarMPISides;
    //Pointers to Arrays ElemToSide, SideToElem usw.
} p4est_SetMPISide_data_t;

//Auxilary COMMON data 
typedef struct p4est_aux_data {

    void *tmp[9]; //Only void ptr are used to pass data to the function
} p4est_aux_data_t;

typedef struct p4est_SetSide_data {
    int CurrentBCSide;
    int CurrentMortrarInnerSide;
    int CurrentInnerSide;
    int CurrentMPISide;
    int CurrentMPIMortarSide;
    // int nInnerSides;
    // int nMPISides;
    //int nMPISides_MINE;//int nMPISides_YOUR;
    // int nMortarMPISides;
    //Pointers to Arrays ElemToSide, SideToElem usw.
} p4est_SetSide_data_t;

void free_data_memory(void *N);

void free_balance_memory(void *N);

#endif /* P4FLUXO_H */
