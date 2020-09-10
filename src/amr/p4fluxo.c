#include "p4fluxo.h"

void free_data_memory(void *N) {
    p4est_fortran_data_t *data;
    data = N;

    pfree(data->MIPtr);

    pfree(data->BCs);

    return;
}

void free_balance_memory(void *N) {
    p4est_balance_data_t *data = (p4est_balance_data_t *) N;

    pfree(data->DataSetU);

    pfree(data->DataSetElem_xGP);


    return;
}

void pfree(void *A) {

    if (A != NULL) {
        free(A);
        A = NULL;
    }
    return;
}