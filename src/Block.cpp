#include "Block.h"

Matrix<double> *Convert_To_Block_Matrix(Matrix<double> *A, int nblocks_local, int *block_sizes_local)
{

  int max_n = 0, count = 0, total_height = 0;
  int *block_bitvec = NULL, *block_start = NULL;
  double *fullcol;
  double infinity = 1.0E300;
  int *block_index_map = NULL;
  int max_block_columns;
  int idx, jcol, k, m, n, next;
  int jb, ib, start, start_idx;

  int col_buf_size, row_buf_size, A_buf_size;
  int max_col_buf_size, max_row_buf_size, max_A_buf_size;

  int *ptr_adr = NULL, *rptr_adr = NULL;
  double *A_adr = NULL;
  int next_ptr, next_rptr, next_Aptr;
  int num;
  int my_nnz = 0;
  

  Matrix<double> *B;
  B = new Matrix<double>(A->world);

  B->all_nbr_cols = new int[B->num_procs];
  B->start_indices = new int[B->num_procs];

  // Filling all_nbr_cols
  MPI_Barrier(B->world);
  MPI_Allgather(static_cast<void *>(&nblocks_local), 1, MPI_INT,
                static_cast<void *>(B->all_nbr_cols), 1, MPI_INT,
                B->world);

  max_block_columns = 0;
  for (int pe=0; pe<B->num_procs; pe++)
  {
    if (max_block_columns<B->all_nbr_cols[pe])
      max_block_columns = B->all_nbr_cols[pe];
  }

  B->my_nbr_cols = B->all_nbr_cols[B->my_id];

  B->n = 0;
  for (int i=0; i<B->num_procs; i++)
  {
    B->n += B->all_nbr_cols[i];
  }


  /* B is square */
  B->m = B->n;

  /*start indices : index of the block column*/
  B->start_indices[0] = 0;
  B->pe = new int[B->n];
  for (int pe=1; pe<B->num_procs; pe++)
  {
    B->start_indices[pe] = B->start_indices[pe-1] + B->all_nbr_cols[pe-1];
  }

  B->my_start_idx = B->start_indices[B->my_id];

  for (int pe=0; pe<(B->num_procs); pe++)
  {
    start_idx = B->start_indices[pe];
    for (int col=0; col<(B->all_nbr_cols[pe]); col++)
      B->pe[start_idx+col] = pe;
  }

  // Define the block size of all the block columns of B
  B->block_sizes = new int[B->n];

  MPI_Barrier(B->world);
  MPI_Allgatherv(static_cast<void *>(block_sizes_local), B->my_nbr_cols, MPI_INT, 
                 static_cast<void *>(B->block_sizes), B->all_nbr_cols, B->start_indices,
                 MPI_INT, B->world);

  B->block_size = 0;

  for (ib=0; ib<(B->n); ib++)
  {
    if (max_n < B->block_sizes[ib])
      max_n = B->block_sizes[ib];
  }

  B->max_block_size = max_n;

  /* Creation of B column sructure */
  B->c_lines = new Compressed_Lines<double>(B->my_nbr_cols, A->symmetric);


  fullcol = new double[A->n * max_n];
  block_bitvec = new int[B->n];

  for (int k=0; k< (A->n)*max_n; k++)
    fullcol[k] = infinity;

  for (ib=0; ib<B->n; ib++)
    block_bitvec[ib] = 0;

  /* Mapping of scalar indices to block indices */
    /// block_index_map[ col in A] = col in B

  int idx_map = 0;

  block_index_map = new int[A->n];
  for (ib=0; ib<B->n; ib++)
  {
    for (int k=0; k<B->block_sizes[ib]; k++)
      block_index_map[idx_map + k] = ib;

    idx_map += B->block_sizes[ib];
  }


  /* block_start[ib] gives the index of the first scalar column in block i */ 
  block_start = new int[B->n];
  block_start[0] = 0;
  for (ib=1; ib<B->n; ib++)
    block_start[ib] = block_start[ib-1] + B->block_sizes[ib-1];

  /* Fill the columns */
  //Sweep sub block columns in a processor
  //
  jb = 0;
  for (int startj=0; jb<max_block_columns; startj += block_sizes_local[jb], jb++)
  {
    if (jb < B->my_nbr_cols)
    {
      n = block_sizes_local[jb];
      count = 0; /* number of blocks in this block column jb */
      total_height = 0; /*number of entries in column block jb (only height not width)*/

      /*fill fullcol*/
      //Sweep scalar columns in block column jb
      for (int j=startj; j<startj+n; j++)
      {
        //Seep elements of A in scalar column j
        for (int k=0; k<A->c_lines->len_cols[j]; k++)
        {
          idx = A->c_lines->col_idcs[j][k]; //row index of element [j][k]
          jcol = j - startj;
          fullcol[jcol*A->n+idx] = A->c_lines->A[j][k];

          ib = block_index_map[idx];
          if (!block_bitvec[ib])
          {
            block_bitvec[ib] = 1;
            count++;
            total_height += B->block_sizes[ib];
          }
        }
      }


      /* Install block column jb */
      B->c_lines->col_idcs[jb] = new int[count];
      B->c_lines->A[jb] = new double[total_height * n];

      B->c_lines->len_cols[jb] = count;
      B->c_lines->len_scalar[jb] = total_height;

      my_nnz += total_height;

      /*Sweep Block rows of B in block column jb*/
      k = 0;
      next = 0;
      for (ib=0; ib < B->n; ib++)
      {
        /*if this block row is present */
        if (block_bitvec[ib])
        {
          B->c_lines->col_idcs[jb][k] = ib;
          start = block_start[ib];
          m = B->block_sizes[ib];
          for (int col=0; col<n; col++)
          {
            for (int row=0; row<m; row++)
            {
              if (fullcol[start + col*A->n + row] == infinity)
              {
                fullcol[start + col*A->n + row] = 0.0;
              }
              B->c_lines->A[jb][next + col*m + row] = fullcol[start + col*A->n + row];
              fullcol[start + col*A->n + row] = infinity;
            }
          }

          block_bitvec[ib] = 0;

          k++;
          next += m*n;
        }
      }
    }
  }

  /*i and j : column and row indices in the scalar matrix
   * ib and jb : column and row indices in the block matrix*/
  
  if (A->c_lines->row_idcs_buf)
  {
    B->c_lines->row_idcs = new int*[B->my_nbr_cols];
    jb = 0;
    for (int startj=0; jb<max_block_columns; startj += block_sizes_local[jb], jb++)
    {
      if (jb < B->my_nbr_cols)
      {
        n = block_sizes_local[jb];
        count = 0;  /*number of blocks in this block row jb*/

       /*fill fullcol (actually fullrow)*/
       for (int j=startj; j<startj+n; j++)
       {
        for (k=0; k<A->c_lines->len_rows[j]; k++)
        {

          idx = A->c_lines->row_idcs[j][k];
          jcol = j - startj;
          fullcol[jcol*A->n + idx] = 1.0;

          ib = block_index_map[idx];
          if (!block_bitvec[ib])
          {
            block_bitvec[ib] = 1;
            count++;
            block_start[ib] = idx;
          }
        }
       }

       /*Install block row structure*/
       B->c_lines->row_idcs[jb] = new int[count];
       B->c_lines->len_rows[jb] = count;
       k = 0;
       next = 0;
       for (ib=0; ib<B->n; ib++)
       {
         if (block_bitvec[ib])
         {
           B->c_lines->row_idcs[jb][k] = ib;
           k++;
           next += m*n;
           block_bitvec[ib] = 0;
         }
       }
      }
    }

  }

  /*Convert row and column data to consistent buffers*/

  col_buf_size = 0;
  row_buf_size = 0;
  A_buf_size = 0;

  idx = B->my_start_idx;
  for (int j=0; j<B->my_nbr_cols; j++, idx++)
  {
    col_buf_size += B->c_lines->len_cols[j];
    if (B->c_lines->row_idcs)
      row_buf_size += B->c_lines->len_rows[j];
    A_buf_size += (B->block_sizes[idx] * B->c_lines->len_scalar[j]);
  }

  MPI_Barrier(A->world);
  MPI_Allreduce(static_cast<void *>(&col_buf_size), static_cast<void *>(&max_col_buf_size), 1,
                MPI_INT, MPI_MAX, A->world);

  if (B->c_lines->row_idcs)
    MPI_Allreduce(static_cast<void *>(&row_buf_size), static_cast<void *>(&max_row_buf_size), 1,
                  MPI_INT, MPI_MAX, A->world);

  MPI_Allreduce(static_cast<void *>(&A_buf_size), static_cast<void *>(&max_A_buf_size), 1,
                MPI_INT, MPI_MAX, A->world);

  B->c_lines->col_idcs_buf = new int[max_col_buf_size];
  if (B->c_lines->row_idcs)
    B->c_lines->row_idcs_buf = new int[max_row_buf_size];
  B->c_lines->col_buf = new double[max_A_buf_size];

  ptr_adr = (B->c_lines->col_idcs_buf);
  rptr_adr = (B->c_lines->row_idcs_buf);
  A_adr = (B->c_lines->col_buf);
  idx = B->my_start_idx; 

  next_ptr  = 0;
  next_rptr = 0;
  next_Aptr = 0;

  for (int j=0;
      j<B->my_nbr_cols;
      j++, idx++)
  {

    num = B->c_lines->len_cols[j];
    memcpy(&(B->c_lines->col_idcs_buf[next_ptr]), B->c_lines->col_idcs[j], num*sizeof(int));
    next_ptr += B->c_lines->len_cols[j];
    delete [] B->c_lines->col_idcs[j];
    B->c_lines->col_idcs[j] = ptr_adr;
    ptr_adr += B->c_lines->len_cols[j];

    if (B->c_lines->row_idcs)
    {
      num = B->c_lines->len_rows[j];
      memcpy(&(B->c_lines->row_idcs_buf[next_rptr]), B->c_lines->row_idcs[j], num*sizeof(int));
      next_rptr += B->c_lines->len_rows[j];
      delete [] B->c_lines->row_idcs[j];
      B->c_lines->row_idcs[j] = rptr_adr;
      rptr_adr += B->c_lines->len_rows[j];
    }

    num = B->block_sizes[idx] * B->c_lines->len_scalar[j];
    memcpy(&(B->c_lines->col_buf[next_Aptr]), B->c_lines->A[j], num*sizeof(double));
    next_Aptr += B->block_sizes[idx] * B->c_lines->len_scalar[j];

    delete [] B->c_lines->A[j];
    B->c_lines->A[j] = A_adr;
    A_adr += B->block_sizes[idx] * B->c_lines->len_scalar[j];

  }

	/* Get length of all cols*/
	B->len_all_cols = new int[B->n];
	MPI_Barrier(B->world);
	MPI_Allgatherv(static_cast<void *>(B->c_lines->len_cols), 
	               B->my_nbr_cols, MPI_INT,
	               static_cast<void *>(B->len_all_cols), 
	               B->all_nbr_cols, B->start_indices, 
	               MPI_INT, B->world);
	
	/* Get length of all rows */
	B->len_all_rows = new int[B->n];
	MPI_Barrier(B->world);
	MPI_Allgatherv(static_cast<void *>(B->c_lines->len_rows), 
                    B->my_nbr_cols, MPI_INT,
                    static_cast<void *>(B->len_all_rows), 
                    B->all_nbr_cols, B->start_indices, 
                    MPI_INT, B->world);

  int max = 0;


	/* Get the maximum number of nnz per column/row of all pes */
	for (int i = 0; i < B->my_nbr_cols; i++)
		if (B->c_lines->len_rows[i] > max)
			max = B->c_lines->len_rows[i];

	for (int i=0; i< B->my_nbr_cols; i++)
		if (B->c_lines->len_cols[i] > max)
			max = B->c_lines->len_cols[i];

  B->my_nnz = my_nnz;

	MPI_Barrier(B->world);
	MPI_Allreduce(&max, &B->max_nnz, 1, MPI_INT, MPI_MAX, B->world);

  B->max_nnz *= B->max_block_size * B->max_block_size;


	/* Initialize the remote transfer buffer */
	B->remote_col_buf = new double[B->max_nnz];
	memset(B->remote_col_buf, 0, B->max_nnz * sizeof(double));

	B->remote_col_send = new double[B->max_nnz];
	memset(B->remote_col_send, 0, B->max_nnz * sizeof(double));

	B->remote_col_idcs_buf = new int[B->max_nnz];
	memset(B->remote_col_idcs_buf, 0, B->max_nnz * sizeof(int));

	B->remote_row_idcs_buf = new int[B->max_nnz];
	memset(B->remote_row_idcs_buf, 0, B->max_nnz * sizeof(int));

  delete [] fullcol;
  delete [] block_bitvec;
  delete [] block_start;
  delete [] block_index_map;

  return B;
}



void write_block(FILE *fptr, double *a, int m, int n)
{
  double val;
  static int count = 0;
  const char *str;
  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
    {
      val = a[i+m*j];
      if (val < 0.0) 
          str = " ";
      else
          str = "  ";

      fprintf(fptr, "%s%.13e", str, val);
      count += 1;
    }
    fprintf(fptr, "\n");
  }
}
