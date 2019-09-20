#include "p4fluxo.h"

void free_data_memory(void *N)
{
    p4est_fortran_data_t *data;
    data = N;

    pfree(data->nNbProc);

    pfree(data->nMPISides_Proc);

    pfree(data->nMPISides_MINE_Proc);

    pfree(data->nMPISides_YOUR_Proc);

    pfree(data->offsetMPISides_MINE);

    pfree(data->offsetMPISides_YOUR);

    pfree(data->EtSPtr);

    pfree(data->StEPtr);

    pfree(data->ChngElementPtr);

    // pfree(data->ChngSidePtr);

    pfree(data->MTPtr);

    pfree(data->MIPtr);

    pfree(data->BCs);
    // p4est_fortran_data_t *p4est_fortran_data;
    pfree(data);
    return;
}

void free_balance_memory(void *N)
{
    p4est_balance_data_t *data = (p4est_balance_data_t *)N;
    // p4est_fortran_data_t *data;
    // data = N;

    pfree(data->DataSetU);

    pfree(data->DataSetElem_xGP);

   
    return;
}

void pfree(void *A)
{
 
    if (A!=NULL)
    {
     free(A); 
     A=NULL;
    }
    return;
}