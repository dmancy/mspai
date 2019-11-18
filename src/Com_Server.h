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


#ifndef GUARD_COM_SERVER_H
#define GUARD_COM_SERVER_H


//file includings
#include "Matrix.h"
#include "Index_Set.h"
#include "Hash_Table.h"
#include "Pattern.h"


//C++ includings
#include <mpi.h>
#include <iostream>
#include <stdexcept>


//Messages with these tags are 
// handled by com_server

/// communication id for sending
/// request to remote pe
const int get_rc_tag        = 1;

/// Communication id for requesting
/// a column of the target matrix
const int get_Bc_tag        = 2;

/// communication id for requesting
/// a column from a remote pe
const int request_Mcol_tag  = 3;

/// communication id for putting 
/// remotly computed preconditioner
/// column into own cache
const int put_Mcol_tag      = 4;

/// communication id that this pe
/// has finished
const int im_done_tag       = 5;

/// communication id that every 
/// pe has finished
const int done_signal_tag   = 6;

/// communication id if error 
/// occured within MPI communication
const int exit_tag          = 7;

/// Communication id for reqeusting
/// a column of the start pattern
const int get_p_tag         = 8;

/// Communication id for reqeusting
/// a column of the upper pattern
const int get_up_tag        = 9;


// Messages with these tags aren't 

/// communication id for requesting 
/// column indices
const int send_cols_tag     = 10;

/// communication id for requesting
/// row indices 
const int send_rows_tag     = 11;

/// communication id for requesting
/// column values 
const int send_vals_tag     = 12;

/// communication id for 
const int send_Mcol_tag     = 13;

/// communication id for sending 
/// column indices of remotly computed 
/// preconditioner column
const int put_Mcol_inds_tag = 14;

/// communication id for sending
/// column values of remotly computed 
/// preconditioner column
const int put_Mcol_vals_tag = 15;

/// communication id for sending
/// the start pattern indices
const int send_pcol_tag     = 16;

/// communication id for sending
/// the start pattern length
const int send_plen_tag     = 17;

/// communication id for sending
/// the upper pattern indices
const int send_upcol_tag    = 18;

/// communication id for sending
/// the upper pattern length
const int send_uplen_tag    = 19;

/// communication id for requesting
/// target column indices  
const int send_Bcols_tag    = 20;

/// communication id for requesting
/// target column values 
const int send_Bvals_tag    = 21;


/////////////////////////////////////////
///     \class Com_Server
///     \brief This class is responsible 
///            for the MPI communication 
///            between the remote cluster 
///            nodes
/////////////////////////////////////////
template <class T>
class Com_Server
{
    
    public:

        //============================================================
        //=========== Template methods - see Com_Server.imp ==========
        //============================================================

        ///////////////////////////////////////////////
        ///     \brief Initialize variables within
        ///            communication world. 
        ///
        ///     \todo Transform this into constructor
        /////////////////////////////////////////////// 
        void    Init();

        ///////////////////////////////////////////////
        ///     \brief Getter method to receive finished
        ///            status. Have all pes finished 
        ///            their work? 
        ///
        ///     \return  Whether all pe's have finished
        ///              or not
        ///////////////////////////////////////////////
        bool    Get_all_done();
        
        /////////////////////////////////////////////////////
        ///     \brief Communication method which maps the 
        ///            request to the specific handle method 
        ///            due to communication id. 
        ///     
        ///     Every time some MPI communication should 
        ///     be done between pes, this while loop is 
        ///     invoked. The MPI_Iprobe function is still
        ///     testing some requests and maps to the handle
        ///     method due to the communication id.
        ///     After a request was successful, check whether
        ///     all pes have finished their work.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        /////////////////////////////////////////////////////
        void    Communicate(Matrix<T> *A, 
                            Matrix<T> *&M,
                            Matrix<T> *B,
                            Pattern   *P,
                            Pattern   *UP);
        
        /////////////////////////////////////////////////////
        ///     \brief Tell Master pe = 0 that I've finished 
        ///            my work.
        ///     
        ///     Each pe will send a short message to the 
        ///     master pe, that he finished his work. This 
        ///     way it is possible to count the number of 
        ///     finished pes and receive the status that all
        ///     have finished at some time.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        /////////////////////////////////////////////////////
        void    Say_Im_Done(Matrix<T> *A,
                            Matrix<T> *&M,
                            Matrix<T> *B, 
                            Pattern   *P,
                            Pattern   *UP);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Some pe has finished his work and 
        ///            requests a new column of another pe.
        ///     
        ///     Each pe will send a short message to the 
        ///     master pe, that he finished his work. This 
        ///     way it is possible to count the number of 
        ///     finished pes and receive the status that all
        ///     have finished at some time.
        ///     This is part of the load balancing mechanism.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param pe The pe from which to request the
        ///               column
        ///     \param index The column index which is to be 
        ///                  set by the remote pe
        ///     \return Whether the remote pe has a  
        ///             column which has to be solved or not
        /////////////////////////////////////////////////////
        bool    Get_M_Col(  Matrix<T>   *A, 
                            Matrix<T>   *M,
                            Matrix<T>   *B,
                            Pattern     *P,
                            Pattern     *UP,
                            const int   pe, 
                            int&        index);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Some pe needs a start pattern column
        ///            from a remote pe.
        ///     
        ///     This is part of the load balancing mechanism.
        ///     If some pe has finished his work, he 
        ///     previously requested a preconditioner column
        ///     and needs now the specific start pattern 
        ///     column.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param col The start pattern column to request
        ///     \param J The index set which is to be set
        /////////////////////////////////////////////////////
        void    Get_P_Col(Matrix<T>   *A,
                          Matrix<T>   *M,
                          Matrix<T>   *B,
                          Pattern     *P,
                          Pattern     *UP,
                          const int   col,
                          Index_Set*& J);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Some pe needs an upper pattern column
        ///            from a remote pe.
        ///     
        ///     This is part of the load balancing mechanism.
        ///     If some pe has finished his work, he 
        ///     previously requested a preconditioner column
        ///     and needs now the specific upper pattern 
        ///     column.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param pe The pe to request the column from
        ///     \param idx The upper pattern column index 
        ///                to be requested
        ///     \param J The index set which is to be set
        /////////////////////////////////////////////////////
        void    Get_UP_Col(Matrix<T>   *A, 
                           Matrix<T>   *M,
                           Matrix<T>   *B,
                           Pattern     *P,
                           Pattern     *UP,
                           const int   pe,
                           int         idx,
                           Index_Set*& J);
        
                   
        ///////////////////////////////////////////////////////
        ///     \brief Requesting column data of the input 
        ///            matrix.
        ///     
        ///     If the requested column data is local, get
        ///     it from local cache. If not, look into the
        ///     hash table if this column was requested before.
        ///     If not, send the request to the specific pe and
        ///     receive the data. Finally store this request 
        ///     into hash table for later requests. Thus the
        ///     MPI traffic can be reduced.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param U_UP Union of upper pattern columns
        ///     \param col The comlumn which is requested
        ///     \param col_len The length of the requested column
        ///     \param row_len The length of the reqeusted row
        ///     \param col_idcs_buf Column indices of requested
        ///                         column
        ///     \param row_idcs_buf Row indices of requested row
        ///     \param col_buf  Values of requested column
        ///     \param ht Hash table to look for possible 
        ///               stored previously requested columns
        ///     \param pre_k_param Number of columns to prerequest
        ///                        while using upper pattern.
        ///     \param use_pre Whether to make the prerequeting code
        ///                    accessible or not.
        ///     \param pre_max_param Whether to prerequest at the 
        ///                          beginning only once or not.
        ///////////////////////////////////////////////////////
        void    Get_Col(Matrix<T>       *A, 
                        Matrix<T>       *&M,
                        Matrix<T>       *B,
                        Pattern         *P,
                        Pattern         *UP,
                        Index_Set*      U_UP,
                        int             col,
                        int&            col_len,
                        int&            row_len, 
                        int*&           col_idcs_buf, 
                        int*&           row_idcs_buf,
                        T*&             col_buf,
                        Hash_Table<T>   *&ht,
                        const int&      pre_k_param,
                        const bool      use_pre,
                        const int       pre_max_param);
        
        
        ///////////////////////////////////////////////////////
        ///     \brief Requesting column data of the target 
        ///            matrix.
        ///     
        ///     If the requested column data is local, get
        ///     it from local cache. 
        ///     If not, send the request to the specific pe and
        ///     receive the data. 
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param col The comlumn which is requested
        ///     \param col_len The length of the requested column
        ///     \param col_idcs_buf Column indices of requested
        ///                         column
        ///     \param col_buf  Values of requested column
        ///////////////////////////////////////////////////////
        void    Get_Target_Col( Matrix<T>   *A,
                                Matrix<T>   *&M,
                                Matrix<T>   *B, 
                                Pattern     *P,
                                Pattern     *UP,
                                int         col,
                                int&        col_len,
                                int*&       col_idcs_buf, 
                                T*&         col_buf);
        
        
        ///////////////////////////////////////////////////////
        ///     \brief Requesting column data of the input 
        ///            matrix.
        ///     
        ///     The requested column data is local, get
        ///     it from local cache. 
        ///
        ///     \param A Local matrix chunk.
        ///     \param col_len The length of the requested column.
        ///     \param row_len The length of the reqeusted row.
        ///     \param col_idcs_buf Column indices of requested
        ///                         column.
        ///     \param row_idcs_buf Row indices of requested row.
        ///     \param col_buf  Values of requested column.
        ///     \param idx Index of column which has to be extracted.
        ///////////////////////////////////////////////////////
        void    Get_Local_Col(  Matrix<T>   *A, 
                                int&        col_len,
                                int&        row_len, 
                                int*&       col_idcs_buf, 
                                int*&       row_idcs_buf,
                                T*&         col_buf,
                                const int&  idx);
        
        
        //============================================================
        //======== Template specifications for double matrices =======
        //============================================================
        
        
        ///////////////////////////////////////////////////////
        ///     \brief Requesting column data of the input 
        ///            matrix.
        ///     
        ///     A column has to be requested remotely from another 
        ///     pe.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param col_idx The comlumn which is requested
        ///     \param col_len The length of the requested column
        ///     \param row_len The length of the reqeusted row
        ///     \param col_idcs_buf Column indices of requested
        ///                         column
        ///     \param row_idcs_buf Row indices of requested row
        ///     \param col_buf  Values of requested column
        ///     \param ht Hash table to look for possible 
        ///               stored previously requested columns
        ///     \param pe Pe which stores this column locally.
        ///     \param idx The index of the column which has to 
        ///                be requested.
        ///////////////////////////////////////////////////////
        void    Get_Remote_Col( Matrix<double>      *A, 
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
                                int&                idx);
        
        
        ///////////////////////////////////////////////////////
        ///     \brief Requesting column data of the target 
        ///            matrix.
        ///     
        ///     A column has to be requested remotely from another 
        ///     pe.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param col_idx The comlumn which is requested
        ///     \param col_len The length of the requested column
        ///     \param col_idcs_buf Column indices of requested
        ///                         column
        ///     \param col_buf  Values of requested column
        ///     \param pe Pe which stores this column locally.
        ///     \param idx The index of the column which has to 
        ///                be requested.
        ///////////////////////////////////////////////////////
        void    Get_Remote_Target_Col(  Matrix<double>      *A,  
                                        Matrix<double>      *&M,
                                        Matrix<double>      *B, 
                                        Pattern             *P,
                                        Pattern             *UP,
                                        int                 col_idx,
                                        int&                col_len,
                                        int*&               col_idcs_buf, 
                                        double*&            col_buf,
                                        int&                pe,
                                        int&                idx);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Inserting computed preconditioner 
        ///            column into local cache if the 
        ///            index is local, else send remotly 
        ///            computed data to remote pe.
        ///
        ///     If the computed column is local, set it into
        ///     local cache. Else send the computed data
        ///     to the column-specific remote pe. This may 
        ///     occur when a pe already finished his own work
        ///     and requested a new column from a remote pe 
        ///     who was still working. This is part of the 
        ///     load balancing mechanism.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param col The computed column to insert
        ///     \param mk_Hat The solution vector containing
        ///                   results
        ///     \param J Index set which marks which indices
        ///              of mk_Hat have to be inserted into 
        ///              the preconditioner
        /////////////////////////////////////////////////////
        void    Insert_Row_Solution(Matrix<double>      *A,
                                    Matrix<double>      *&M,
                                    Matrix<double>      *B,
                                    Pattern             *P,
                                    Pattern             *UP,
                                    const int           col,
                                    double              *mk_Hat, 
                                    const Index_Set     *J);
        
        
        void    Insert_Row_Solution_Block(Matrix<double>      *A,
                                    Matrix<double>      *&M,
                                    Matrix<double>      *B,
                                    Pattern             *P,
                                    Pattern             *UP,
                                    const int           col,
                                    double              *mk_Hat, 
                                    const Index_Set     *J);
        
        //============================================================
        //======= Template specifications for COMPLEX matrices =======
        //============================================================
                
        
        ///////////////////////////////////////////////////////
        ///     \brief Requesting column data of the input 
        ///            matrix.
        ///     
        ///     A column has to be requested remotely from another 
        ///     pe.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param col_idx The comlumn which is requested
        ///     \param col_len The length of the requested column
        ///     \param row_len The length of the reqeusted row
        ///     \param col_idcs_buf Column indices of requested
        ///                         column
        ///     \param row_idcs_buf Row indices of requested row
        ///     \param col_buf  Values of requested column
        ///     \param ht Hash table to look for possible 
        ///               stored previously requested columns
        ///     \param pe Pe which stores this column locally.
        ///     \param idx The index of the column which has to 
        ///                be requested.
        ///////////////////////////////////////////////////////
        void    Get_Remote_Col( Matrix<COMPLEX>     *A, 
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
                                int&                idx);
        
        
        ///////////////////////////////////////////////////////
        ///     \brief Requesting column data of the target 
        ///            matrix.
        ///     
        ///     A column has to be requested remotely from another 
        ///     pe.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param col_idx The comlumn which is requested
        ///     \param col_len The length of the requested column
        ///     \param col_idcs_buf Column indices of requested
        ///                         column
        ///     \param col_buf  Values of requested column
        ///     \param pe Pe which stores this column locally.
        ///     \param idx The index of the column which has to 
        ///                be requested.
        ///////////////////////////////////////////////////////
        void    Get_Remote_Target_Col(  Matrix<COMPLEX>     *A,  
                                        Matrix<COMPLEX>     *&M,
                                        Matrix<COMPLEX>     *B, 
                                        Pattern             *P,
                                        Pattern             *UP,
                                        int                 col_idx,
                                        int&                col_len,
                                        int*&               col_idcs_buf, 
                                        COMPLEX*&           col_buf,
                                        int&                pe,
                                        int&                idx);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Inserting computed preconditioner
        ///            column into local cache if the index is
        ///            local, else send remotly computed data to 
        ///            remote pe.
        ///     
        ///     If the computed column is local, set it into
        ///     local cache. Else send the computed data
        ///     to the column-specific remote pe. This may 
        ///     occur when a pe already finished his own work
        ///     and requested a new column from a remote pe 
        ///     who was still working. This is part of the 
        ///     load balancing mechanism.
        ///
        ///     \param A Local matrix chunk
        ///     \param M Local preconditioner chunk
        ///     \param B Local target matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param UP Local upper pattern chunk
        ///     \param col The computed column to insert
        ///     \param mk_Hat The solution vector containing
        ///                   results
        ///     \param J Index set which marks which indices
        ///              of mk_Hat have to be inserted into 
        ///              the preconditioner
        /////////////////////////////////////////////////////
        void    Insert_Row_Solution(Matrix<COMPLEX>     *A,
                                    Matrix<COMPLEX>     *&M,
                                    Matrix<COMPLEX>     *B,
                                    Pattern             *P,
                                    Pattern             *UP,
                                    const int           col,
                                    COMPLEX             *mk_Hat, 
                                    const Index_Set     *J);
        
        
	void Initialize_nbr_done(void);

    private:
        
        //Variables
        
        /// Number of already finished pes
        static int  nbr_done;
        
        /// Flag whether everybody has finished his work
        bool        all_done;
        
        /// Flag to check if user wants to prerequest all columns
        /// at once at the beginning.
        static bool first;



        
        //============================================================
        //========== Template methods - see Com_Server.imp ===========
        //============================================================
        
        /////////////////////////////////////////////////////
        ///     \brief Master pe = 0 receives that another 
        ///            pe has finished his work
        ///     
        ///     Each pe will send a short message to the 
        ///     master pe, that he finished his work. This 
        ///     way it is possible to count the number of 
        ///     finished pes and receive the status that all
        ///     have finished at some time.
        ///
        ///     \param A Local matrix chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void        Handle_Im_Done(Matrix<T> *A, 
                                   int requestor);
        
        /////////////////////////////////////////////////////
        ///     \brief All pes have finished, status will be
        ///            set to all_done = true.
        ///     
        ///     This handle is invoked if all pes have 
        ///     finished their work. Master pe knows that
        ///     everybode has finished. 
        ///
        ///     \param A Local matrix chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Done_Signal(Matrix<T> *A, 
                                   int requestor);
        
        /////////////////////////////////////////////////////
        ///     \brief Checking if all pes have finished their
        ///            work.
        ///
        ///     This method is invoked every time a 
        ///     communication has finished successfully, to
        ///     check whether all pes have finished their 
        ///     work. 
        ///
        ///     \param P Local start pattern chunk
        /////////////////////////////////////////////////////
        void    Check_Done(Pattern *P);
        
        
        
        //============================================================
        //======= Template specifications for double matrices ========
        //============================================================
        
        /////////////////////////////////////////////////////
        ///     \brief Remote pe gets request to send column
        ///            data to the requestor.
        ///     
        ///     This handle is invoked if some pe needs column
        ///     data which is stored on a remote pe. First
        ///     the pe will receive the index of the column
        ///     and then it sends the requested data to the 
        ///     requestor.
        ///
        ///     \param A Local matrix chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Get_Col( const Matrix<double> *A, 
                                const int requestor);  
         
        
        /////////////////////////////////////////////////////
        ///     \brief Remote pe gets request to send column
        ///            data to the requestor.
        ///     
        ///     This handle is invoked if some pe needs column
        ///     data which is stored on a remote pe. First
        ///     the pe will receive the index of the column
        ///     and then it sends the requested data to the 
        ///     requestor.
        ///
        ///     \param B Local target matrix chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Get_Target_Col(const Matrix<double> *B, 
                                      const int    requestor);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Pe receives a computed solution column 
        ///            from a remote pe and must store it.
        ///     
        ///     This handle is invoked if some pe computed 
        ///     a preconditioner column after it has finished
        ///     its own work. This is part of the load balancing
        ///     mechanism. The remote pe receives the computed
        ///     data and must store it into its local cache.
        ///
        ///     \param M Local preconditioner chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Insert_Row_Solution(Matrix<double> *M, 
                                           int requestor);
        
        
        
        //============================================================
        //======= Template specifications for COMPLEX matrices =======
        //============================================================
        
        /////////////////////////////////////////////////////
        ///     \brief Remote pe gets request to send column
        ///            data to the requestor.
        ///     
        ///     This handle is invoked if some pe needs column
        ///     data which is stored on a remote pe. First
        ///     the pe will receive the index of the column
        ///     and then it sends the requested data to the 
        ///     requestor.
        ///
        ///     \param A Local matrix chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Get_Col( const Matrix<COMPLEX> *A, 
                                const int requestor);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Remote pe gets request to send column
        ///            data to the requestor.
        ///     
        ///     This handle is invoked if some pe needs column
        ///     data which is stored on a remote pe. First
        ///     the pe will receive the index of the column
        ///     and then it sends the requested data to the 
        ///     requestor.
        ///
        ///     \param B Local target matrix chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Get_Target_Col(const Matrix<COMPLEX> *B, 
                                      const int     requestor);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Pe receives a computed solution column
        ///            from a remote pe and must store it.
        ///     
        ///     This handle is invoked if some pe computed 
        ///     a preconditioner column after it has finished
        ///     its own work. This is part of the load balancing
        ///     mechanism. The remote pe receives the computed
        ///     data and must store it into its local cache.
        ///
        ///     \param M Local preconditioner chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Insert_Row_Solution(Matrix<COMPLEX> *M, 
                                           int requestor);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Pe is requested to send a unfinished
        ///            column to another remote pe who has
        ///            finished his work already.
        ///     
        ///     This handle is invoked if a remote pe has 
        ///     finished his work and is requesting a not 
        ///     solved column from a remote pe. If this pe
        ///     has some unfinished columns, he will send 
        ///     the index to the requestor. This is part of
        ///     the load balancing mechanism.
        ///
        ///     \param P Local start pattern chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Get_M_Col(Pattern   *P, 
                                 int         requestor);     
           
        
        /////////////////////////////////////////////////////
        ///     \brief Pe is requested to send a start
        ///            pattern column.
        ///     
        ///     This handle is invoked if a remote pe has 
        ///     finished his work and is requesting a 
        ///     start pattern column. This is part of the
        ///     load balancing mechanism. 
        ///
        ///     \param A Local matrix chunk
        ///     \param P Local start pattern chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Get_P_Col(const Matrix<T> *A, 
                                 const Pattern   *P,
                                 const int       requestor);
        
        
        /////////////////////////////////////////////////////
        ///     \brief Pe is requested to send an upper
        ///            pattern column.
        ///     
        ///     This handle is invoked if a remote pe has 
        ///     finished his work and is requesting an
        ///     upper pattern column. This is part of the
        ///     load balancing mechanism. 
        ///
        ///     \param A Local matrix chunk
        ///     \param UP Local upper pattern chunk
        ///     \param requestor Pe which invoked this handle
        /////////////////////////////////////////////////////
        void    Handle_Get_UP_Col(const Matrix<T> *A, 
                                  const Pattern   *UP,
                                  const int       requestor);
};

#include "Com_Server.imp"

#endif
