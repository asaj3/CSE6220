/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>

/*
 * TODO: Implement your solutions here
 */


void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
{
    // TODO
    // Initialize the variables for the Cartesian topology
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    
    // Initialize new communicator
    MPI_Comm comm_column;
    int col_dims[2] = {true, false};
    MPI_Cart_sub(comm, col_dims, &comm_column);
    
    if (coords[1] == 0) {
        
        //Determine the number of elements each processor will receive and the displacements
        int* ele_num = new int[dims[0]];
        int* displace = new int[dims[0]];
        
        for (int i = 0; i < dims[0]; i++) {
            ele_num[i] = block_decompose(n, dims[0], i);
        }
        
        displace[0] = 0;
        for (int i = 1; i < dims[0]; i++) {
            displace[i] = displace[i-1] + ele_num[i-1];
        }
        
        //Determine which processor is the root in the first column
        int rank_root;
        int coords_root[] = {0};
        MPI_Cart_rank(comm_column, coords_root, &rank_root);
        
        //Determine the local size for each processor in the first column
        int local_size = ele_num[coords[0]];
        (*local_vector) = new double[local_size];
        
        //Distribute the values
        MPI_Scatterv(input_vector, ele_num, displace, MPI_DOUBLE, *local_vector, local_size, MPI_DOUBLE, rank_root, comm_column);
        
        //Free up memories after the function finish the computation
        delete [] ele_num;
        delete [] displace;
    }
    return;
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO
    // Initialize the variables for the Cartesian topology
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    
    // Initialize new communicator
    MPI_Comm comm_column;
    int col_dims[2] = {true, false};
    MPI_Cart_sub(comm, col_dims, &comm_column);
    
    if (coords[1] == 0) {
        
        //Determine the number of elements each processor will receive and the displacements
        int* ele_num = new int[dims[0]];
        int* displace = new int[dims[0]];
        
        for (int i = 0; i < dims[0]; i++) {
            ele_num[i] = block_decompose(n, dims[0], i);
        }
        
        displace[0] = 0;
        for (int i = 1; i < dims[0]; i++) {
            displace[i] = displace[i-1] + ele_num[i-1];
        }
        
        //Determine which processor is the root in the first column
        int rank_root;
        int coords_root[] = {0};
        MPI_Cart_rank(comm_column, coords_root, &rank_root);
        
        //Determine the local size for each processor in the first column
        int local_size = ele_num[coords[0]];
        
        //Gather the values
        MPI_Gatherv(local_vector, local_size, MPI_DOUBLE, output_vector, ele_num, displace, MPI_DOUBLE, rank_root, comm_column);
        
        //Free up memories after the function finish the computation
        delete [] ele_num;
        delete [] displace;
    }
    return;
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // TODO
    // Initialize the variables for the Cartesian topology
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    
    // Initialize new communicators
    MPI_Comm comm_column;
    int col_dims[2] = {true, false};
    MPI_Cart_sub(comm, col_dims, &comm_column);
    
    //Initialize a temporary matrix to hole the intermediate values
    double* temp_mat = NULL;
    
    if (coords[1] == 0) {
        
        //Determine the number of elements each processor will receive and the displacements
        int* ele_num = new int[dims[0]];
        int* displace = new int[dims[0]];
        
        for (int i = 0; i < dims[0]; i++) {
            ele_num[i] = n * block_decompose(n, dims[0], i);
        }
        
        displace[0] = 0;
        for (int i = 1; i < dims[0]; i++) {
            displace[i] = displace[i-1] + ele_num[i-1];
        }
        
        //Determine which processor is the root in the first column
        int rank_root;
        int coords_root[] = {0};
        MPI_Cart_rank(comm_column, coords_root, &rank_root);
        
        //Determine the local size for each processor in the first column
        int local_size = ele_num[coords[0]];
        temp_mat = new double[local_size];
        
        //Distribute the values
        MPI_Scatterv(input_matrix, ele_num, displace, MPI_DOUBLE, temp_mat, local_size, MPI_DOUBLE, rank_root, comm_column);
        
        //Free up memories after the function finish the computation
        delete [] ele_num;
        delete [] displace;
    }
    
    // Initialize new communicator
    MPI_Comm comm_row;
    int row_dims[2] = {false, true};
    MPI_Cart_sub(comm, row_dims, &comm_row);
    
    int n_row = block_decompose(n, dims[0], coords[0]);
    int n_col = block_decompose(n, dims[1], coords[1]);
    
    //Determine the number of elements each processor will receive and the displacements
    int* ele_num = new int[dims[1]];
    int* displace = new int[dims[1]];
    
    for (int i = 0; i < dims[1]; i++) {
        ele_num[i] = block_decompose(n, dims[1], i);
    }
    
    displace[0] = 0;
    for (int i = 1; i < dims[1]; i++) {
        displace[i] = displace[i-1] + ele_num[i-1];
    }
    
    //Determine which processor is the root in the first column
    int rank_root;
    int coords_root[] = {0};
    MPI_Cart_rank(comm_row, coords_root, &rank_root);
    
    //Memory allocation for local matrix
    *local_matrix = new double[n_row * n_col];
    
    for (int i = 0; i < n_row; i++) {
        MPI_Scatterv((temp_mat + i * n), ele_num, displace, MPI_DOUBLE, (*local_matrix + i * n_col), n_col, MPI_DOUBLE, rank_root, comm_row);
    }
    //Free up memories after the function finish the computation
    delete [] ele_num;
    delete [] displace;
    delete [] temp_mat;

    return;
}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
    //Determine the rank and the number of processors in each column and row
    int rank, p, q;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);
    q = (int)sqrt(p);
    
    MPI_Comm comm_column, comm_row;
    int col_dims[2] = {true, false};
    int row_dims[2] = {false, true};
    MPI_Cart_sub(comm, col_dims, &comm_column);
    MPI_Cart_sub(comm, row_dims, &comm_row);
    
    //Determine the column and row rank
    int coords[2], rank_col, rank_row;
    MPI_Cart_coords(comm, rank, 2, coords);
    rank_row = coords[0];
    rank_col = coords[1];
    
    //Send and Receive the rows or column (transposing) to the corresponding processors
    if (rank == 0) {
        int cnt = block_decompose(n, q, rank);
        std::memcpy(row_vector, col_vector, cnt * sizeof(double));
    }
    else if (rank_col == 0) {
        int cnt_snt = block_decompose(n, q, rank_row);
        MPI_Send(col_vector, cnt_snt, MPI_DOUBLE, rank_row, 0, comm_row);
    }
    else if (rank_row == rank_col) {
        int cnt_rcv = block_decompose(n, q, rank_row);
        MPI_Recv(row_vector, cnt_rcv, MPI_DOUBLE, 0, 0, comm_row, MPI_STATUS_IGNORE);
    }
    
    int cnt_bcast = block_decompose(n, q, rank_col);
    MPI_Bcast(row_vector, cnt_bcast, MPI_DOUBLE, rank_col, comm_column);
}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    // TODO
    // Initialize the variables for the Cartesian topology
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    
    // Initialize new communicator
    MPI_Comm comm_row;
    int row_dims[2] = {false, true};
    MPI_Cart_sub(comm, row_dims, &comm_row);
    
    int n_row = block_decompose(n, dims[0], coords[0]);
    int n_col = block_decompose(n, dims[1], coords[1]);
    
    //Determine the matrix multiplication for the local vector
    double* local_X = new double[n_col];
    transpose_bcast_vector(n, local_x, local_X, comm);
    
    //Determine the result after multiplication
    double* local_Y = new double[n_row];
    
    if (n_col == n_row) {
        matrix_vector_mult(n_row, local_A, local_X, local_Y);
    }
    else {
        matrix_vector_mult(n_row, n_col, local_A, local_X, local_Y);
    }
    
    //Determine which processor is the root in the first column
    int rank_root;
    int coords_root[2] = {0, 0};
    MPI_Cart_rank(comm_row, coords_root, &rank_root);
    MPI_Reduce(local_Y, local_y, n_row, MPI_DOUBLE, MPI_SUM, rank_root, comm_row);
    
    //Free up memories after the function finish the computation
    delete [] local_X;
    delete [] local_Y;
 
    return;
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{
    // TODO
    double l2_norm_tot;
    // Initialize the variables for the Cartesian topology
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    
    // Initialize new communicators
    MPI_Comm comm_column;
    int col_dims[2] = {true, false};
    MPI_Cart_sub(comm, col_dims, &comm_column);
    
    MPI_Comm comm_row;
    int row_dims[2] = {false, true};
    MPI_Cart_sub(comm, row_dims, &comm_row);
    
    int n_row = block_decompose(n, dims[0], coords[0]);
    int n_col = block_decompose(n, dims[1], coords[1]);
    
    //Determine the diagonal elements, D and the Matrix R, where R = A - D
    double* local_D = new double[n_row];
    double* local_R = new double[n_row * n_col];
    
    if (coords[0] == coords[1]) {
        std::copy(local_A, local_A + n_row * n_col, local_R);
        
        for (int i = 0; i < n_row; i++) {
            local_D[i] = local_A[i * (n_col + 1)];
            local_R[i * (n_col + 1)] = 0;
        }
        
        if (coords[0] != 0) {
            int rank_dest, coords_dest[] = {0};
            MPI_Cart_rank(comm_row, coords_dest, &rank_dest);
            MPI_Send(local_D, n_row, MPI_DOUBLE, rank_dest, 0, comm_row);
        }
    }
    else {
        std::copy(local_A, local_A + n_row * n_col, local_R);
    }
    
    if (coords[1] == 0 && coords[0] != 0) {
        int rank_source, coords_source[] = {coords[0]};
        MPI_Cart_rank(comm_row, coords_source, &rank_source);
        MPI_Recv(local_D, n_row, MPI_DOUBLE, rank_source, 0, comm_row, MPI_STATUS_IGNORE);
    }
    
    if (coords[1] == 0) {
        for (int i = 0; i < n_row; i++) {
            local_x[i] = 0;
        }
    }
    
    int i = 1;
    while (i <= max_iter) {
        
        double* local_y = new double[n_row];
        
        distributed_matrix_vector_mult(n, local_R, local_x, local_y, comm);
        if (coords[1] == 0) {
            for (int j = 0; j < n_row; j++) {
                local_x[j] = (local_b[j] - local_y[j]) / local_D[j];
            }
        }
        
        distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);
        if (coords[1] == 0) {
            double l2_norm = 0.0;
            for (int j = 0; j < n_row; j++) {
                l2_norm = l2_norm + (local_b[j] - local_y[j]) * (local_b[j] - local_y[j]);
            }
            MPI_Allreduce(&l2_norm, &l2_norm_tot, 1, MPI_DOUBLE, MPI_SUM, comm_column);
        }
        
        l2_norm_tot = sqrt(l2_norm_tot);
        int rank_root, coords_root[] = {0};
        MPI_Cart_rank(comm_row, coords_root, &rank_root);
        MPI_Bcast(&l2_norm_tot, 1, MPI_DOUBLE, rank_root, comm_row);
        
        delete [] local_y;
        if (l2_norm_tot <= l2_termination) {
            return;
        }
        i++;
    }
    
    //Free up memories after the function finish the computation
    delete [] local_D;
    delete [] local_R;
}


// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double* A,
                            double* x, double* y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_x = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &x[0], &local_x, comm);

    // allocate local result space
    double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
                int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);

    // allocate local result space
    double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}
