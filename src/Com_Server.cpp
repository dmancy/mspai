/*
    ======================================================================
    ======================================================================
    ==                                                                  ==
    ==  MSPAI:  Modified SPAI algorithm to comupte SParse Approximate   ==
    ==          Invers matrices.                                        ==
    ==                                                                  ==
    ==  Copyright (C)  2007, 2008, 2009 by                              ==
    ==                 Matous Sedlacek <sedlacek@in.tum.de>             ==
    ==                 Chair of Scientific Computing -- Informatics V   ==
    ==                 Technische Universität München                   ==
    ==                                                                  ==
    ==  This file is part of MSPAI.                                     ==
    ==                                                                  ==
    ==  MSPAI is free software: you can redistribute it and/or          ==
    ==  modify it under the terms of the GNU Lesser General Public      ==
    ==  License as published by the Free Software Foundation, either    ==
    ==  version 3 of the License, or (at your option) any later version.==
    ==                                                                  ==
    ==  MSPAI is distributed in the hope that it will be useful,        ==
    ==  but WITHOUT ANY WARRANTY; without even the implied warranty of  ==
    ==  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   ==
    ==  GNU Lesser General Public License for more details.             ==
    ==                                                                  ==
    ==  You should have received a copy of the GNU Lesser General       ==
    ==  Public License along with MSPAI.                                ==
    ==  If not, see <http://www.gnu.org/licenses/>.                     ==
    ==                                                                  ==
    ======================================================================
    ======================================================================
*/


//file includings
#include "Com_Server.h"


//============================================================================
//============================================================================
//================ Template specifications for double matrices ===============
//============================================================================
//============================================================================
template<> void 
Com_Server<double>::Initialize_nbr_done(void)
{
	nbr_done = 0;
}

template<>  void 
Com_Server<double>::Handle_Get_Col(const Matrix<double> *A, 
                                   const int     requestor)
{
    
    int             idx,
                    len_col = 0,
                    slen_col = 0,
                    len_row = 0,
                    next = 0,
                    index = 0,
                    *col_idcs_buf = NULL,
                    *row_idcs_buf = NULL;
                                
    double          *col_buf      = NULL,
                    *ptr_buf = NULL;
                    

    
    MPI_Status      status;
    MPI_Comm        world;
    
    Compressed_Lines<double>   
                    *lines = NULL;
    

    world = A->world;
    lines = A->c_lines;
    
    // Local pe is receiving column index of the
    // column which another pe is requesting remotely.
    if (MPI_Recv(static_cast<void*> (&idx), 1, MPI_INT, 
        requestor, get_rc_tag, world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Col\" by MPI_Recv.\n");


    // Getting the column data from local cache
    len_col         = lines->len_cols[idx];
    slen_col        = lines->len_scalar[idx] * A->block_sizes[idx + A->my_start_idx]; 
    col_idcs_buf    = lines->col_idcs[idx];
    col_buf         = lines->A[idx];
    if (A->c_lines->row_idcs)
    {
        row_idcs_buf    = lines->row_idcs[idx];
        len_row         = lines->len_rows[idx];
    }


    for (int i = 0; i < len_col; i++)
    {

      for (int j = 0; j < A->block_sizes[idx + A->my_start_idx] * A->block_sizes[col_idcs_buf[i]]; j++)
      {
        A->remote_col_buf[index++] = lines->A[idx][next + j];
      }
      next += A->block_sizes[idx + A->my_start_idx] * A->block_sizes[col_idcs_buf[i]];
      
    }
    

    // Send then synchronously the remotely 
    // requested column data
    
    // Sending column indices
    if (MPI_Send(static_cast<int*> (col_idcs_buf), len_col, MPI_INT, 
        requestor, send_cols_tag, world))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Col\" by MPI_Send.\n");        
    
    
    // If matrix is not symmetric, send row indices
    if (A->c_lines->row_idcs)
    {
        if (MPI_Send(static_cast<int*> (row_idcs_buf), len_row, MPI_INT, 
            requestor, send_rows_tag, world))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Handle_Get_Col\" by MPI_Send.\n");
    }

    // Sending column values 
    if (MPI_Send(static_cast<double*> (A->remote_col_buf), slen_col, MPI_DOUBLE, 
        requestor, send_vals_tag, world))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Col\" by MPI_Send.\n");

}



template<>  void 
Com_Server<double>::Handle_Get_Target_Col(const Matrix<double> *B, 
                                          const int     requestor)
{
    
    int             idx,
                    len_col = 0,
                    *col_idcs_buf = NULL;
                                                    
    double          *col_buf = NULL;
    
    MPI_Status      status;
    MPI_Comm        world;
    
    Compressed_Lines<double>   
                    *lines = NULL;
    

    world = B->world;
    lines = B->c_lines;
    
    // Local pe is receiving column index of the
    // column which another pe is requesting remotely.
    if (MPI_Recv(static_cast<void*> (&idx), 1, MPI_INT, 
        requestor, get_Bc_tag, world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Target_Col\" by MPI_Recv.\n");


    // Getting the column data from local cache
    len_col         = lines->len_cols[idx];
    col_idcs_buf    = lines->col_idcs[idx];
    col_buf         = lines->A[idx];


    // Send then synchronously the remotely 
    // requested column data
    
    // Sending column indices
    if (MPI_Send(static_cast<int*> (col_idcs_buf), len_col, 
        MPI_INT, requestor, send_Bcols_tag, world))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Target_Col\" by MPI_Send.\n");        
    

    // Sending column values 
    if (MPI_Send(static_cast<double*> (col_buf), len_col, 
        MPI_DOUBLE, requestor, send_Bvals_tag, world))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Target_Col\" by MPI_Send.\n");
}



template<>  void 
Com_Server<double>::Handle_Insert_Row_Solution(Matrix<double>   *M, 
                                               int              requestor)
{
    
    
    int         idx_and_len[3],
                idx, 
                nnz,
                *col_idcs_buf = NULL;
    
    double      *col_buf = NULL;
    
    MPI_Status  status;
    
    
    // Receiving data where to store the remotely computed 
    // column m_k. idx: where to store, nnz: length of data
    // to store
    if (MPI_Recv(static_cast<void *> (idx_and_len), 3, MPI_INT, 
                 requestor, put_Mcol_tag, M->world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Insert_Row_Solution\" by MPI_Recv.\n");
        
    idx = idx_and_len[0];
    nnz = idx_and_len[1];

    
    // Allocate memory and fill it with remote data
     
    col_idcs_buf = new int[nnz];
    
    if (MPI_Recv(static_cast<void *> (col_idcs_buf), nnz, MPI_INT, 
                requestor, put_Mcol_inds_tag, M->world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Insert_Row_Solution\" by MPI_Recv.\n");
    
        
    col_buf = new double[nnz];  
    
    if (MPI_Recv(static_cast<void *> (col_buf), nnz, MPI_DOUBLE, 
                requestor, put_Mcol_vals_tag, M->world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Insert_Row_Solution\" by MPI_Recv.\n");

    // Insert remotely computed data into local 
    // preconditioner solution    
    M->c_lines->col_idcs[idx]       = col_idcs_buf;
    M->c_lines->A[idx]              = col_buf;
    M->c_lines->len_cols[idx]       = nnz;
}



template<>  void 
Com_Server<double>::Get_Remote_Col( Matrix<double>      *A, 
                                    Matrix<double>      *&M,
                                    Matrix<double>      *B,
                                    Pattern             *P,
                                    Pattern             *UP,
                                    int                 col_idx,
                                    int&                col_len,
                                    int&                row_len, 
                                    int*&               col_idcs_buf, 
                                    int*&               row_idcs_buf,
                                    double*&            col_buf,
                                    Hash_Table<double> *&ht,
                                    int&                pe,
                                    int&                idx)
{       
    int                 send_flag,
                        cols_flag,
                        rows_flag,
                        vals_flag,
                        scalar_len;
        
    MPI_Status          send_status,
                        cols_status,
                        rows_status,
                        vals_status;
                
    MPI_Request         get_cols_request, 
                        get_rows_request,
                        get_vals_request, 
                        send_request;

    MPI_Comm            world;
        
    Timer               o_timer;


    world = A->world;   

    //Element is not in cache - get it from a remote pe
    if (MPI_Irecv(static_cast<void*> (A->remote_col_idcs_buf), 
        A->max_nnz, MPI_INT, pe, send_cols_tag, 
        world, &get_cols_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Irecv.\n");


    if (A->c_lines->row_idcs)
    {
        if (MPI_Irecv(static_cast<void*> (A->remote_row_idcs_buf), 
            A->max_nnz, MPI_INT, pe, send_rows_tag, 
            world, &get_rows_request))
                throw std::runtime_error(
                        "\n\tERROR in method \"Com_Server::"
                        "Get_Remote_Col\" by MPI_Irecv.\n");
    }


    if (MPI_Irecv(static_cast<void*> (A->remote_col_buf), 
        A->max_nnz, MPI_DOUBLE, pe, send_vals_tag, 
        world, &get_vals_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Irecv.\n");


    // Send the request
    if (MPI_Isend(static_cast<void*> (&idx), 1, MPI_INT, 
        pe, get_rc_tag, world, &send_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Isend.\n");


    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&send_request, &send_flag, &send_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Test.\n");
    }
    while (!send_flag);


    // Service requests until the data comes back. 
    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&get_cols_request, &cols_flag, &cols_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Test.\n");
    }
    while (!cols_flag);


    if (A->c_lines->row_idcs)
    {
        do 
        {
            Communicate(A, M, B, P, UP);
            if (MPI_Test(&get_rows_request, &rows_flag, &rows_status))
                throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Test.\n");
        }
        while (!rows_flag);
    }


    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&get_vals_request, &vals_flag, &vals_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Test.\n");
    }
    while (!vals_flag);
        
    col_len             = A->len_all_cols[col_idx];;
    col_idcs_buf        = A->remote_col_idcs_buf;   
    col_buf             = A->remote_col_buf;

    if (A->c_lines->row_idcs)   //unsymmetric matrix
    {
        row_len         = A->len_all_rows[col_idx];
        row_idcs_buf    = A->remote_row_idcs_buf;   
    }
    else                        //symmetric matrix
    {
        row_len         = col_len;
        row_idcs_buf    = A->remote_col_idcs_buf;
    }

    // Compute scalar_len

    for (int i = 0; i < col_len; i++)
        scalar_len += A->block_sizes[col_idx] * A->block_sizes[col_idcs_buf[i]];


    // Store the new data into cache 
    if (ht) //Does the hash table exist?
    {  
        // row structure ?
        if (A->c_lines->row_idcs) 
        {
            ht->Insert_Block( col_idx,
                        A->remote_col_idcs_buf,
                        A->remote_row_idcs_buf,
                        A->remote_col_buf,
                        col_len,
                        scalar_len,
                        row_len);
        }
        else 
        {
            ht->Insert_Block( col_idx,
                        A->remote_col_idcs_buf,
                        A->remote_col_idcs_buf,
                        A->remote_col_buf,
                        col_len,
                        scalar_len,
                        col_len);
        }
    }
}



template<>  void 
Com_Server<double>::Get_Remote_Target_Col(  Matrix<double>      *A, 
                                            Matrix<double>      *&M,
                                            Matrix<double>      *B, 
                                            Pattern             *P,
                                            Pattern             *UP,
                                            int                 col_idx,
                                            int&                col_len,
                                            int*&               col_idcs_buf, 
                                            double*&            col_buf,
                                            int&                pe,
                                            int&                idx)
{       
    int                 send_flag,
                        cols_flag,
                        vals_flag;
        
    MPI_Status          send_status,
                        cols_status,
                        vals_status;
                
    MPI_Request         get_cols_request, 
                        get_vals_request, 
                        send_request;

    MPI_Comm            world;
        
    Timer               o_timer;


    world = B->world;   

    //Element is not in cache - get it from a remote pe
    if (MPI_Irecv(static_cast<void*> (B->remote_col_idcs_buf), 
        B->max_nnz, MPI_INT, pe, send_Bcols_tag, 
        world, &get_cols_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Irecv.\n");


    if (MPI_Irecv(static_cast<void*> (B->remote_col_buf), 
        B->max_nnz, MPI_DOUBLE, pe, send_Bvals_tag, 
        world, &get_vals_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Irecv.\n");


    // Send the request
    if (MPI_Isend(static_cast<void*> (&idx), 1, MPI_INT, 
        pe, get_Bc_tag, world, &send_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Isend.\n");


    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&send_request, &send_flag, &send_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Test.\n");
    }
    while (!send_flag);


    // Service requests until the data comes back. 
    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&get_cols_request, &cols_flag, &cols_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Test.\n");
    }
    while (!cols_flag);


    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&get_vals_request, &vals_flag, &vals_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Test.\n");
    }
    while (!vals_flag);
        
    col_len             = B->len_all_cols[col_idx];;
    col_idcs_buf        = B->remote_col_idcs_buf;   
    col_buf             = B->remote_col_buf;
}


template<>  void 
Com_Server<double>::Insert_Row_Solution(Matrix<double>  *A, 
                                        Matrix<double>  *&M, 
                                        Matrix<double>  *B,
                                        Pattern         *P,
                                        Pattern         *UP,
                                        const int       col,
                                        double          *mk_Hat,
                                        const Index_Set *J) 
{
    int         pe,
                idx,
                nnz,
                *col_idcs_buf = NULL,
                flag,
                idx_and_len[3];
            
    double      *col_buf = NULL;
    
    MPI_Request requests[3];
    
    MPI_Status  statuses[3];

    
    pe = A->pe[col];
    idx = col - A->start_indices[pe];

    if (pe == A->my_id)         // The line is local.
    {
        // Just store it into local 
        // preconditioner solution
            
        nnz = J->len;
        M->c_lines->len_cols[idx] = nnz;
            
        col_buf         = new double[nnz];
        col_idcs_buf    = new int[nnz];
    
        M->c_lines->A[idx]          = col_buf;
        M->c_lines->col_idcs[idx]   = col_idcs_buf;
    
        for ( int i = 0; i < nnz; i++)
        {
            col_idcs_buf[i] = J->idcs[i];
            col_buf[i] = mk_Hat[i];
        }
    }
    else                        // The line is remote.
    {
        // This pe made the computation due to
        // load balancing. Send the results to 
        // remote pe which will store the solution
        // within it's own local memory.
        
        nnz = J->len;
        
        idx_and_len[0] = idx;
        idx_and_len[1] = nnz;
        idx_and_len[2] = nnz;

        // Prepare remote pe for index of column and length of data
        // which will be sent. 
        if (MPI_Isend(static_cast<void *> (idx_and_len), 3, MPI_INT, 
            pe, put_Mcol_tag, A->world, &requests[0]))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Insert_Row_Solution\" by MPI_Isend.\n");
    
        // Send the data.
        
        if (MPI_Isend(static_cast<void *> (J->idcs), nnz, MPI_INT, 
            pe, put_Mcol_inds_tag, A->world, &requests[1]))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Insert_Row_Solution\" by MPI_Isend.\n");
        
        if (MPI_Isend(static_cast<void *> (mk_Hat), nnz, MPI_DOUBLE, 
            pe, put_Mcol_vals_tag, A->world, &requests[2]))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Insert_Row_Solution\" by MPI_Isend.\n");
            
        // Check if remote pe got the data
        do 
        {
            Communicate(A, M, B, P, UP);
            if (MPI_Testall(3, requests, &flag, statuses))
                throw std::runtime_error(
                    "\n\tERROR in method \"Com_Server::"
                    "Insert_Row_Solution\" by MPI_Testall.\n"); 
        }
        while (! flag);
    }
}

template<>  void 
Com_Server<double>::Insert_Row_Solution_Block(Matrix<double>  *A, 
                                        Matrix<double>  *&M, 
                                        Matrix<double>  *B,
                                        Pattern         *P,
                                        Pattern         *UP,
                                        const int       col,
                                        double          *mk_Hat,
                                        const Index_Set *J) 
{
    int         pe,
                idx,
                nnz,
                scalar_nnz,   // number of scalar entries in the M column
                *col_idcs_buf = NULL,
                flag,
                idx_and_len[3];
            
    double      *col_buf = NULL;
    
    MPI_Request requests[3];
    
    MPI_Status  statuses[3];

    
    pe = A->pe[col];
    idx = col - A->start_indices[pe];

    if (pe == A->my_id)         // The line is local.
    {
        // Just store it into local 
        // preconditioner solution
            
        nnz = J->len; //number of non zeros blocks
        M->c_lines->len_cols[idx] = nnz;

        scalar_nnz = J->slen * A->block_sizes[col];
            
        col_buf         = new double[scalar_nnz];
        col_idcs_buf    = new int[nnz];
    
        M->c_lines->A[idx]          = col_buf;
        M->c_lines->col_idcs[idx]   = col_idcs_buf;

        idx = 0;
    
        /*
        for ( int i = 0; i < nnz; i++)
        {
            col_idcs_buf[i] = J->idcs[i];
            for (int k = 0; k < A->block_sizes[J->idcs[i]]*A->block_sizes[col]; k++)
            {
              col_buf[idx] = mk_Hat[idx];
              idx++;
            }
        }
        */
        memcpy(col_idcs_buf, J->idcs, nnz * sizeof(int));
        memcpy(col_buf, mk_Hat, J->slen * A->block_sizes[col] *  sizeof(double));
    }
    else                        // The line is remote.
    {
      printf("bad!!!!!\n");
        // This pe made the computation due to
        // load balancing. Send the results to 
        // remote pe which will store the solution
        // within it's own local memory.
        
        nnz = J->len;
        
        idx_and_len[0] = idx;
        idx_and_len[1] = nnz;
        idx_and_len[2] = nnz;

        // Prepare remote pe for index of column and length of data
        // which will be sent. 
        if (MPI_Isend(static_cast<void *> (idx_and_len), 3, MPI_INT, 
            pe, put_Mcol_tag, A->world, &requests[0]))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Insert_Row_Solution\" by MPI_Isend.\n");
    
        // Send the data.
        
        if (MPI_Isend(static_cast<void *> (J->idcs), nnz, MPI_INT, 
            pe, put_Mcol_inds_tag, A->world, &requests[1]))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Insert_Row_Solution\" by MPI_Isend.\n");
        
        if (MPI_Isend(static_cast<void *> (mk_Hat), nnz, MPI_DOUBLE, 
            pe, put_Mcol_vals_tag, A->world, &requests[2]))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Insert_Row_Solution\" by MPI_Isend.\n");
            
        // Check if remote pe got the data
        do 
        {
            Communicate(A, M, B, P, UP);
            if (MPI_Testall(3, requests, &flag, statuses))
                throw std::runtime_error(
                    "\n\tERROR in method \"Com_Server::"
                    "Insert_Row_Solution\" by MPI_Testall.\n"); 
        }
        while (! flag);
    }
}



//============================================================================
//============================================================================
//=============== Template specifications for COMPLEX matrices ===============
//============================================================================
//============================================================================


template<>  void 
Com_Server<COMPLEX>::Handle_Get_Col(const Matrix<COMPLEX> *A, 
                                    const int requestor)
{
    
    int                         idx,
                                len_col = 0,
                                len_row = 0,
                                *col_idcs_buf = NULL,
                                *row_idcs_buf = NULL;
                                
    COMPLEX                     *col_buf = NULL;
    
    MPI_Status                  status;
    MPI_Comm                    world;
    
    Compressed_Lines<COMPLEX>   *lines = NULL;
    
    
    MPI_Datatype                COMPLEX_MPI;
    
    
    // Build own MPI datatype. It contains 
    // two values, one for the real and one
    // for the imaginary part.
    MPI_Type_contiguous(2, MPI_DOUBLE, &COMPLEX_MPI);
    MPI_Type_commit(&COMPLEX_MPI);

    world = A->world;
    lines = A->c_lines;
    
    // Local pe is receiving column index of the
    // column which another pe is requesting remotely.
    if (MPI_Recv(static_cast<void*> (&idx), 1, MPI_INT, 
        requestor, get_rc_tag, world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Col\" by MPI_Recv.\n");

    
    // Getting the column data from local cache
    len_col         = lines->len_cols[idx];
    col_idcs_buf    = lines->col_idcs[idx];
    col_buf         = lines->A[idx];
    if (A->c_lines->row_idcs)
    {
        row_idcs_buf    = lines->row_idcs[idx];
        len_row         = lines->len_rows[idx];
    }
    
    // Send then synchronously the remotely 
    // requested column data
    
    // Sending column indices
    if (MPI_Send(static_cast<int*> (col_idcs_buf), len_col, MPI_INT, 
        requestor, send_cols_tag, world))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Col\" by MPI_Send.\n");
    
    
    // If matrix is not symmetric, sending row indices as well
    if (A->c_lines->row_idcs)
    {
        if (MPI_Send(static_cast<int*> (row_idcs_buf), len_row, MPI_INT, 
            requestor, send_rows_tag, world))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Handle_Get_Col\" by MPI_Send.\n");
    }
    
    
    // Sending column values
    if (MPI_Send(static_cast<COMPLEX*> (col_buf), len_col, COMPLEX_MPI, 
        requestor, send_vals_tag, world))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Col\" by MPI_Send.\n");
}



template<>  void 
Com_Server<COMPLEX>::Handle_Get_Target_Col( const Matrix<COMPLEX> *B, 
                                            const int requestor)
{
    
    int                         idx,
                                len_col = 0,
                                *col_idcs_buf = NULL;
                                                                
    COMPLEX                     *col_buf = NULL;
    
    MPI_Status                  status;
    MPI_Comm                    world;
    
    Compressed_Lines<COMPLEX>   *lines = NULL;
    
    
    MPI_Datatype                COMPLEX_MPI;
    
    
    // Build own MPI datatype. It contains 
    // two values, one for the real and one
    // for the imaginary part.
    MPI_Type_contiguous(2, MPI_DOUBLE, &COMPLEX_MPI);
    MPI_Type_commit(&COMPLEX_MPI);

    world = B->world;
    lines = B->c_lines;
    
    // Local pe is receiving column index of the
    // column which another pe is requesting remotely.
    if (MPI_Recv(static_cast<void*> (&idx), 1, MPI_INT, 
        requestor, get_Bc_tag, world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Target_Col\" by MPI_Recv.\n");

    
    // Getting the column data from local cache
    len_col         = lines->len_cols[idx];
    col_idcs_buf    = lines->col_idcs[idx];
    col_buf         = lines->A[idx];

    
    // Send then synchronously the remotely 
    // requested column data
    
    // Sending column indices
    if (MPI_Send(static_cast<int*> (col_idcs_buf), len_col, 
        MPI_INT, requestor, send_Bcols_tag, world))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Target_Col\" by MPI_Send.\n");
    
    
    // Sending column values
    if (MPI_Send(static_cast<COMPLEX*> (col_buf), len_col, 
        COMPLEX_MPI, requestor, send_Bvals_tag, world))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Get_Target_Col\" by MPI_Send.\n");
}



template<>  void 
Com_Server<COMPLEX>::Handle_Insert_Row_Solution(Matrix<COMPLEX> *M, 
                                               int requestor)
{
    int             idx_and_len[3],
                    idx, 
                    nnz,
                    *col_idcs_buf = NULL;
    
    COMPLEX         *col_buf = NULL;
    
    MPI_Status      status;
    
    MPI_Datatype    COMPLEX_MPI;
    
    
    // Build own MPI datatype. It contains 
    // two values, one for the real and one
    // for the imaginary part.
    MPI_Type_contiguous(2, MPI_DOUBLE, &COMPLEX_MPI);
    MPI_Type_commit(&COMPLEX_MPI);
    
    // Receiving data where to store the remotely computed 
    // column m_k. idx: where to store, nnz: length of data
    // to store
    if (MPI_Recv(static_cast<void *> (idx_and_len), 3, MPI_INT, 
                 requestor, put_Mcol_tag, M->world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Insert_Row_Solution\" by MPI_Recv.\n");
    
    idx = idx_and_len[0];
    nnz = idx_and_len[1];
    
    // Allocate memory and fill it with remote data
     
    col_idcs_buf = new int[nnz];
    
    if (MPI_Recv(static_cast<void *> (col_idcs_buf), nnz, MPI_INT, 
                requestor, put_Mcol_inds_tag, M->world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Insert_Row_Solution\" by MPI_Recv.\n");
        
    col_buf = new COMPLEX[nnz]; 
    
    if (MPI_Recv(static_cast<void *> (col_buf), nnz, COMPLEX_MPI, 
                requestor, put_Mcol_vals_tag, M->world, &status))
        throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Handle_Insert_Row_Solution\" by MPI_Recv.\n");

    
    // Insert remotely computed data into local 
    // preconditioner solution 
    M->c_lines->col_idcs[idx]       = col_idcs_buf;
    M->c_lines->A[idx]              = col_buf;
    M->c_lines->len_cols[idx]       = nnz;
}



template<>  void 
Com_Server<COMPLEX>::Get_Remote_Col(Matrix<COMPLEX>     *A, 
                                    Matrix<COMPLEX>     *&M,
                                    Matrix<COMPLEX>     *B,
                                    Pattern             *P,
                                    Pattern             *UP,
                                    int                 col_idx,
                                    int&                col_len,
                                    int&                row_len, 
                                    int*&               col_idcs_buf, 
                                    int*&               row_idcs_buf,
                                    COMPLEX*&           col_buf,
                                    Hash_Table<COMPLEX> *&ht,
                                    int&                pe,
                                    int&                idx)
{
    int                 send_flag,
                        cols_flag,
                        rows_flag,
                        vals_flag;
        
    MPI_Status          send_status,
                        cols_status,
                        rows_status,
                        vals_status;
                
    MPI_Request         get_cols_request, 
                        get_rows_request,
                        get_vals_request, 
                        send_request;
        
    MPI_Comm            world;
        

    world = A->world;   

    MPI_Datatype COMPLEX_MPI;
    MPI_Type_contiguous(2, MPI_DOUBLE, &COMPLEX_MPI);
    MPI_Type_commit(&COMPLEX_MPI);

    if (MPI_Irecv(static_cast<void*> (A->remote_col_idcs_buf), 
        A->max_nnz, MPI_INT, pe, send_cols_tag, 
        world, &get_cols_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Irecv.\n");


    if (A->c_lines->row_idcs)
    {
        if (MPI_Irecv(static_cast<void*> (A->remote_row_idcs_buf), 
            A->max_nnz, MPI_INT, pe, send_rows_tag, 
            world, &get_rows_request))
                throw std::runtime_error(
                    "\n\tERROR in method \"Com_Server::"
                    "Get_Remote_Col\" by MPI_Irecv.\n");
    }


    if (MPI_Irecv(static_cast<void*> (A->remote_col_buf), 
        A->max_nnz, COMPLEX_MPI, pe, send_vals_tag, 
        world, &get_vals_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Irecv.\n");
        

    // Send the request
    if (MPI_Isend(static_cast<void*> (&idx), 1, MPI_INT, 
        pe, get_rc_tag, world, &send_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Isend.\n");

    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&send_request, &send_flag, &send_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Test.\n");
    }
    while (!send_flag);


    // Service requests until the data comes back. 
    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&get_cols_request, &cols_flag, &cols_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Col\" by MPI_Test.\n");
    }
    while (!cols_flag);


    if (A->c_lines->row_idcs)
    {
        do 
        {
            Communicate(A, M, B, P, UP);
            if (MPI_Test(&get_rows_request, &rows_flag, &rows_status))
                throw std::runtime_error(
                    "\n\tERROR in method \"Com_Server::"
                    "Get_Remote_Col\" by MPI_Test.\n");
        }
        while (!rows_flag);
    }


    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&get_vals_request, &vals_flag, &vals_status))
            throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Get_Remote_Col\" by MPI_Test.\n");
    }
    while (!vals_flag);

    col_len         = A->len_all_cols[col_idx];
    col_idcs_buf    = A->remote_col_idcs_buf;   
    col_buf         = A->remote_col_buf;

    if (A->c_lines->row_idcs)   //unsymmetric matrix
    {
        row_len         = A->len_all_rows[col_idx];
        row_idcs_buf    = A->remote_row_idcs_buf;   
    }
    else                        //symmetric matrix
    {
        row_len         = col_len;
        row_idcs_buf    = A->remote_col_idcs_buf;
    }


    // Store the new gotten data into cache 
    if (ht) //Does the hash table exist?
    {  
        // row structure ?
        if (A->c_lines->row_idcs) 
        {
            ht->Insert( col_idx,
                        A->remote_col_idcs_buf,
                        A->remote_row_idcs_buf,
                        A->remote_col_buf,
                        col_len,
                        row_len);
        }
        else 
        {
            ht->Insert( col_idx,
                        A->remote_col_idcs_buf,
                        A->remote_col_idcs_buf,
                        A->remote_col_buf,
                        col_len,
                        col_len);
        }
    }
}



template<>  void 
Com_Server<COMPLEX>::Get_Remote_Target_Col( Matrix<COMPLEX>     *A, 
                                            Matrix<COMPLEX>     *&M,
                                            Matrix<COMPLEX>     *B,
                                            Pattern             *P,
                                            Pattern             *UP,
                                            int                 col_idx,
                                            int&                col_len,
                                            int*&               col_idcs_buf, 
                                            COMPLEX*&           col_buf,
                                            int&                pe,
                                            int&                idx)
{
    int                 send_flag,
                        cols_flag,
                        vals_flag;
        
    MPI_Status          send_status,
                        cols_status,
                        vals_status;
                
    MPI_Request         get_cols_request, 
                        get_vals_request, 
                        send_request;
        
    MPI_Comm            world;
        

    world = B->world;   

    MPI_Datatype COMPLEX_MPI;
    MPI_Type_contiguous(2, MPI_DOUBLE, &COMPLEX_MPI);
    MPI_Type_commit(&COMPLEX_MPI);

    if (MPI_Irecv(static_cast<void*> (B->remote_col_idcs_buf), 
        B->max_nnz, MPI_INT, pe, send_Bcols_tag, 
        world, &get_cols_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Irecv.\n");


    if (MPI_Irecv(static_cast<void*> (B->remote_col_buf), 
        B->max_nnz, COMPLEX_MPI, pe, send_Bvals_tag, 
        world, &get_vals_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Irecv.\n");
        

    // Send the request
    if (MPI_Isend(static_cast<void*> (&idx), 1, MPI_INT, 
        pe, get_Bc_tag, world, &send_request))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Isend.\n");

    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&send_request, &send_flag, &send_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Test.\n");
    }
    while (!send_flag);


    // Service requests until the data comes back. 
    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&get_cols_request, &cols_flag, &cols_status))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Get_Remote_Target_Col\" by MPI_Test.\n");
    }
    while (!cols_flag);


    do 
    {
        Communicate(A, M, B, P, UP);
        if (MPI_Test(&get_vals_request, &vals_flag, &vals_status))
            throw std::runtime_error(
            "\n\tERROR in method \"Com_Server::"
            "Get_Remote_Target_Col\" by MPI_Test.\n");
    }
    while (!vals_flag);

    col_len         = B->len_all_cols[col_idx];
    col_idcs_buf    = B->remote_col_idcs_buf;   
    col_buf         = B->remote_col_buf;
}



template<>  void 
Com_Server<COMPLEX>::Insert_Row_Solution(   Matrix<COMPLEX> *A, 
                                            Matrix<COMPLEX> *&M,
                                            Matrix<COMPLEX> *B,
                                            Pattern         *P,
                                            Pattern         *UP, 
                                            const int col,
                                            COMPLEX *mk_Hat,
                                            const Index_Set *J) 
{
    
    int         pe,
                idx,
                nnz,
                *col_idcs_buf = NULL,
                flag,
                idx_and_len[3];
            
    COMPLEX     *col_buf = NULL;
    
    MPI_Request requests[3];
    
    MPI_Status  statuses[3];

    
    pe = A->pe[col];
    idx = col - A->start_indices[pe];

    if (pe == A->my_id)         // The line is local. 
    {
        // Just store it into local 
        // preconditioner solution
            
        nnz = J->len;
        M->c_lines->len_cols[idx] = nnz;
            
        col_buf         = new COMPLEX[nnz];
        col_idcs_buf    = new int[nnz];
    
        M->c_lines->A[idx]          = col_buf;
        M->c_lines->col_idcs[idx]   = col_idcs_buf;
    
        for ( int i = 0; i < nnz; i++)
        {
            col_idcs_buf[i] = J->idcs[i];
            col_buf[i] = mk_Hat[i];
        }
    }
    else                        // The line is remote.
    {
        // This pe made the computation due to
        // load balancing. Send the results to 
        // remote pe which will store the solution
        // within it's own local memory.
            
        // Build own MPI datatype. It contains 
        // two values, one for the real and one
        // for the imaginary part.
        MPI_Datatype COMPLEX_MPI;
        MPI_Type_contiguous(2, MPI_DOUBLE, &COMPLEX_MPI);
        MPI_Type_commit(&COMPLEX_MPI);
        
        nnz = J->len;
        
        idx_and_len[0] = idx;
        idx_and_len[1] = nnz;
        idx_and_len[2] = nnz;

        // Prepare remote pe for index of column and length of data
        // which will be sent. 
        if (MPI_Isend(static_cast<void *> (idx_and_len), 2, MPI_INT, 
            pe, put_Mcol_tag, A->world, &requests[0]))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Insert_Row_Solution\" by MPI_Isend.\n");
    
        // Send the computed solution to remote pe 
        
        if (MPI_Isend(static_cast<void *> (J->idcs), nnz, MPI_INT, 
            pe, put_Mcol_inds_tag, A->world, &requests[1]))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Insert_Row_Solution\" by MPI_Isend.\n");
        
        if (MPI_Isend(static_cast<void *> (mk_Hat), nnz, COMPLEX_MPI, 
            pe, put_Mcol_vals_tag, A->world, &requests[2]))
            throw std::runtime_error(
                "\n\tERROR in method \"Com_Server::"
                "Insert_Row_Solution\" by MPI_Isend.\n");
            
        // Check if remote pe got the data
        do 
        {
            Communicate(A, M, B, P, UP);
            if (MPI_Testall(3, requests, &flag, statuses))
                throw std::runtime_error(
                    "\n\tERROR in method \"Com_Server::"
                    "Insert_Row_Solution\" by MPI_Testall.\n"); 
        }
        while (! flag);
    }
}
