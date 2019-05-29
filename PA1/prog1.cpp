#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <math.h>

int dboard(int N);      //Initializing the dartboard function

int main(int argc, char* argv[]) {
    //Initializing the MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    double initT = MPI_Wtime();    //Initialing time for the random seeds
    srand(time(NULL) + rank * 10000);   //Randomizing the seed for each rank
    int N, R;
    float Pi = 0;

    //Checking for valid arguements entered from command line. Master procesor
    //will store the arguements.
    if (rank == 0) {
        if (argc < 3) {
            printf("\nNot enough arguements. Program will be TERMINATED\n");
            exit (EXIT_FAILURE);
        }
        N = atoi(argv[1]);      //Converting string to int type
        R = atoi(argv[2]);      //Converting string to int type
        if (N == 0 || R == 0) {
            printf("\nInvalid arguements. Program will be TERMINATED\n");
            exit (EXIT_FAILURE);
        }
    }
    //Master processor broadcast the value of N and R to all processors.
    MPI_Bcast(&N, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&R, 1, MPI_DOUBLE, 0, comm);
    double M, MSum, t0, t1, difTime, maxTime;
    double totTime = 0;
    double piSum[R];

    
    //Number of rounds the dartboard algorithm is called for each processor
    for (int i = 0; i < R; i++) {
        t0 = MPI_Wtime();        //The beginning runtime
        M = dboard(N);          //Each processor call the dartboard algorithm
        //MPI_Reduce is used to sum up the M values from each processor to
        //the Master processor
        MPI_Reduce(&M, &MSum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        t1 = MPI_Wtime();       //The end runtime
        //Max runtime among all processors is passed to the Master processor
        difTime = t1 - t0;
        MPI_Reduce(&difTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        //Master processor using the formula to estimate the value of Pi
        //and totalling the runtime
        if (rank == 0) {
            totTime = totTime + maxTime;
            piSum[i] = 2.0f * (N/MSum);
        }
    }
    //Master processor finding the average of the Pi values from the R rounds.
    if (rank == 0) {
        for (int i = 0; i < R; i++) {
            Pi = Pi + piSum[i];
        }
        //Printing out the output
        printf("N = %d, R = %d, P = %d, Pi = %f\n", N, R, size, Pi/R);
        printf("Time taken is %0.3fs\n", totTime);      //The final runtime
    }
    return MPI_Finalize();
}

//The dartboard function
int dboard(int N) {
    //Initializing variables
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n = N/size; //Each processor will throw N/size darts
    int M = 0;
    int r = 1;      //The radius of the dartboard is set to 1
    float num_r[n], theta[n], x[n], y[n];
    //Generate random numbers for X and Y coordinates
    for (int i = 0; i < n; i++) {
        num_r[i] = 1.0f * rand()/RAND_MAX;
        theta[i] = 2.0f * rand()/RAND_MAX * M_PI;
        x[i] = sqrt(num_r[i]) * cos(theta[i]);
        y[i] = sqrt(num_r[i]) * sin(theta[i]);
    }
    // Checking whether dart lands in square
    for (int i = 0; i < n; i++) {
        if (sqrt(x[i]*x[i]) <= r / sqrt(2) && sqrt(y[i]*y[i]) <= r / sqrt(2)) {
            M += 1;
        }
    }
    return M;   //Returning the number of darts that land in the square
}
