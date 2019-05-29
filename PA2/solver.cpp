#include "solver.h"


/*************************** DECLARE YOUR HELPER FUNCTIONS HERE ************************/
// A recursive function to solve the N Queen problem based on sequential www.geeksforgeeks.org/n-queen-problem-backtracking-3/
// which was greatly modified to work with parralel and MPI
bool solveNQ(std::vector<std::vector<unsigned int> >& board, int col, int n, int k, int& solCounter, int& irecvCounter, int& flag, int& solSize, std::vector<std::vector<unsigned int> >& part_solns, std::vector<unsigned int>& str_solns, MPI_Request& req, MPI_Status& status);
// Function to check if the Queen is safe from other attacking Queen
bool safeCheck(std::vector<std::vector<unsigned int> >& board, int row, int col, int n);
// When all Queens have been safely placed, the solution will be printed
void printSol(std::vector<std::vector<unsigned int> >& board, int n, std::vector<std::vector<unsigned int> >& part_solns);




/*************************** solver.h functions ************************/

void seq_solver(unsigned int n, std::vector<std::vector<unsigned int> >& all_solns)
{
    int numSol;     // Number of total Valid solutions
    int solCounter, irecvCounter, flag, solSize;    // To count number of partial solutions. Not used in seq_solver
    
    MPI_Request req;
    MPI_Status status;
    // Initialize the vectors
    std::vector<std::vector<unsigned int> > part_solns;        // Not use in seq_solver
    std::vector<unsigned int> str_solns;        // The solutions are saved in one long vector
    std::vector<std::vector<unsigned int> > board;      // 2d vectors of the board
    
    // Initilize the board to have NxN zeroes
    board.resize(n);
    for (int i = 0; i < n; i++)
        board[i].resize(n, 0);
    
    // Calling the solver function
    solveNQ(board, 0, n, n, solCounter, irecvCounter, flag, solSize, part_solns, str_solns, req, status);
    
    // Cutting the long vector solution and pushing the value into all_solns which is a 2d vector
    numSol = str_solns.size()/n;
    
    all_solns.resize(numSol);       // Resizing the all_solns according to the number of solutions
    int m = 0;
    for (int i = 0; i < numSol; i++)
    {
        all_solns[i].resize(n);
        for (int j = 0; j < n; j++)
        {
            all_solns[i][j] = str_solns[m];
            m++;
        }
    }
}


void nqueen_master(    unsigned int n,
                   unsigned int k,
                   std::vector<std::vector<unsigned int> >& all_solns)
{
    
    // Initialize some variables
    int size;
    int solCounter, irecvCounter, ctr;  // To count number of partial solutions
    int solSize;                        // The size of valid solutions from each worker
    int numSol;                         // Number of total valid solutions
    int flag = -1;                      // Flag is used when using MPI_Irecv and MPI_Test
    int partSolSize;                    // Number of k partial solutions
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    MPI_Request req;
    
    // Initialize the vectors
    std::vector<std::vector<unsigned int> > part_solns;     // Vector containing the partial solutions
    std::vector<std::vector<unsigned int> > board;          // 2d vectors of the board
    std::vector<unsigned int> str_solns;                    // Vector containing complete solution from workers
    
    // Initilize the board to have NxN zeroes
    board.resize(n);
    for (int i = 0; i < n; i++)
        board[i].resize(n, 0);
    
    solCounter = 0;         // This counter is to count the number of partial solutions sent to worker
    irecvCounter = 0;       // This counter is to count the number of MPI_Irecv is called in Master, it should
    // match with the number of send from workers
    
    // Calling the solver function
    // Note, all the distribution of partial solutions are done in the solveNQ function
    // since the function is recursive and partial solutions have to be sent to worker if
    // it is found without waiting for Master to find all the partial solutions. The Master will continuosly
    // find the partial solution without waiting for worker to send completed solution.
    solveNQ(board, 0, n, k, solCounter, irecvCounter, flag, solSize, part_solns, str_solns, req, status);
    // NOTE: MASTER can finish finding all the k partial solutions before workers completed their work, thus
    // the MPI_Irecv &status and &request is passed to the nqueen_master function to continue waiting and
    // receieving the valid solutions from workers
    
    // Since Master can potentionally finish finding all the k partial solutions before workers, the number
    // of total k partial solutions have to be checked.
    partSolSize = part_solns.size();
    // Master will receive the size of the completed solution from the worker first. If no valid solution
    // exist from worker, the size will be 0
    
    // This is where Master will loop if it finishes finding all the k partial solutions before sending them to
    // all workers. Master keeps track of the number of k partial solutions that has been sent and also the number
    // of MPI_Irecv that it calls.
    while (irecvCounter < (partSolSize-(size-1)))
    {
        // Ensuring that Master will only create a new memory block for MPI_Irecv when there is a worker
        // that will sent it's work
        if (flag != 0)
        {
            // Master has to receive the size of the valid solution from worker first
            MPI_Irecv(&solSize, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &req);
            flag = 0;
        }
        
        // Master will keep looping until a worker sends their completed solution and Master will receive it
        // and push it into a string vector solutions. Immediately after that, Master will send the next
        // k partial solution to the worker that finishes its work
        MPI_Test(&req, &flag, &status);
        if (flag != 0)
        {
            
            // MPI_Irecv counter to counter the total succesfully Receive and Send from Master and Worker
            irecvCounter++;
            
            // If workers cannot find a valid solution, it will send a zero integer for the solSize
            if (solSize != 0)
            {
                // When valid size obtained from workers, the str_solns size will be increased to accomodate
                // new completed solutions from workers
                str_solns.resize(str_solns.size() + solSize);
                
                // Master receive the completed solution from the worker that sent the size of the solution
                // and will store it in one long vector str_solns
                MPI_Recv(&str_solns[str_solns.size() - solSize], solSize, MPI_UNSIGNED, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
            }
            
            // Master will send the next k partial solutions to the worker that just completed its previous work
            MPI_Send(&part_solns[solCounter][0], n, MPI_UNSIGNED, status.MPI_SOURCE, 111, MPI_COMM_WORLD);
            solCounter++;   // Counting a successful send of k partial solution to workers
            flag = -1;
        }
        
    }
    
    // This is important for the special case where number of partial solutions is less
    // than the number of processor used
    if (solCounter < size)
        ctr = solCounter + 1;
    else
        ctr = size;
    
    // This will be the final receive from the Master and no more sending of k partial solutions, instead
    // Master will send Kill Signal to workers
    int finalRecv = 1;
    while (finalRecv < ctr)
    {
        // Master will receive the size of the completed solution from the worker first. If no valid solution
        // exist from worker, the size will be 0
        if (flag != 0)
        {
            MPI_Irecv(&solSize, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &req);
            flag = 0;
        }
        MPI_Test(&req, &flag, &status);
        if (flag != 0)
        {
            finalRecv++;
            if (solSize != 0)
            {
                // When valid size obtained from workers, the str_solns size will be increased to accomodate
                // new completed solutions from workers
                str_solns.resize(str_solns.size() + solSize);
                //                printf("MASTER : str_solns size in MASTER %lu\n", str_solns.size());
                
                // Master receive the completed solution from the worker that sent the size of the solution
                // and will store it in one long vector str_solns
                MPI_Recv(&str_solns[str_solns.size() - solSize], solSize, MPI_UNSIGNED, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
                flag = -1;
            }
        }
    }
    // Master send the kill signal to ALL workers in the form of MPI_TAG 999
    for (int i = 1; i < size; i++)
        MPI_Send(0, 0, MPI_INT, i, 999, MPI_COMM_WORLD);
    
    // Cutting the long vector solution and pushing the value into all_solns which is a 2d vector
    numSol = str_solns.size()/n;
    all_solns.resize(numSol);       // Resizing the all_solns according to the number of solutions
    int m = 0;
    for (int i = 0; i < numSol; i++)
    {
        all_solns[i].resize(n);
        for (int j = 0; j < n; j++)
        {
            all_solns[i][j] = str_solns[m];
            m++;
        }
    }
}



void nqueen_worker(    unsigned int n,
                   unsigned int k)
{
    // Initialize some variables
    int rank;
    int solSize;                         // The size of valid solutions from each worker
    int solCounter, irecvCounter;        // To count number of partial solutions. Not used in seq_solver
    int flag = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    MPI_Request req;
    
    
    std::vector<std::vector<unsigned int> > part_solns;     // 2d Vector to store each valid solutions
    std::vector<unsigned int> str_solns;                    // 1d Vector to store complete solution
    std::vector<unsigned int> single_solns;                 // Vector containing the partial solutions from Master
    // Initilize the board to have NxN zeroes
    std::vector<std::vector<unsigned int> > board;
    board.resize(n);
    for (int i = 0; i < n; i++)
        board[i].resize(n, 0);
    
    single_solns.resize(n);
    
    // All workers will constantly loop and received a partial solution from Master and find
    // a valid solution and send the solution back to Master. Workers will repeat until a kill
    // signal from Master is sent to each workers to stop
    while (1)
    {
        // All the vectors are clear and re-initialize at the begining of each loop
        str_solns.clear();
        single_solns.clear();
        single_solns.resize(n);
        board.clear();
        board.resize(n);
        for (int i = 0; i < n; i++)
            board[i].resize(n, 0);
        
        // Workers received the partial solutions from Master, this can also be a kill signal
        // depending on the MPI_TAG.
        MPI_Recv(&single_solns[0], n, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        // Checking if kill signal is sent from Master in the form of MPI_TAG 999
        if (status.MPI_TAG == 999)
            return;     // worker exit function on receiving kill signal
        
        // Workers create the board based on the partial solutions received from Master
        for (int i = 0; i < k; i++)
            board[single_solns[i]][i] = 1;
        
        // Initialize the bool variable
        bool solFound = true;
        
        // Workers run the solveNQ function based on the partially completed solution from
        // Master and continue to find valid solutions
        if (solveNQ(board, k, n, k, solCounter, irecvCounter, flag, solSize, part_solns, str_solns, req, status) == false)
        {
            // If no valid solution found, worker will send a 0 size solution
            solSize = 0;
            MPI_Send(&solSize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            // Solution not found
            solFound = false ;
        }
        // If solution is found, send the solution to Master
        if (solFound)
        {
            
            // Workers will convert the 2d vector of complete solution into 1d contiguos vector
            // MPI cannot send non-contiguous 2d vectors without using special datatype
            solSize = str_solns.size();

            // Workers will sent the size of the solution to Master first. This is due to the
            // fact that the partial solution received from Master can have multiple valid
            // solutions and some workers will have longer valid solutions than others
            MPI_Send(&solSize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&str_solns[0], solSize, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
        }
        // Workers will loop and continue working and receiving new partial solutions from Master
    }
    return;
}



/*************************** DEFINE YOUR HELPER FUNCTIONS HERE ************************/



// A recursive function to solve the N Queen problem based on sequential www.geeksforgeeks.org/n-queen-problem-backtracking-3/
// which was greatly modified to work with parralel and MPI
bool solveNQ(std::vector<std::vector<unsigned int> >& board, int col, int n, int k, int& solCounter, int& irecvCounter, int& flag, int& solSize, std::vector<std::vector<unsigned int> >& part_solns, std::vector<unsigned int>& str_solns, MPI_Request& req, MPI_Status& status)
{
    // Initialize some variables
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    static int r = 1;   // This variable corresponds to the worker rank for the Master to sent partial solution to
    
    // If Master and not seq_solver will do this part. Master will continuosly find k partial solutions for the N x N problem
    if (rank == 0 && size > 1)
    {
        // Each k partial solutions will first be sent to each worker. Then, the Master will receive the valid solutions
        // from worker and the next partial solutions will be sent to each worker again
        // Master will continuos generate k partial solutions and sending to available worker
        if (col == k)
        {
            // The k partial solution is stored in part_solns vector
            printSol(board, n, part_solns);

            // Master will initiate the first send of the k partial solutions to all workers
            if (r < size)
            {
                MPI_Send(&part_solns[solCounter][0], n, MPI_UNSIGNED, r, 111, MPI_COMM_WORLD);
                solCounter++;       // Counting the k partial solutions that have been sent to workers
            }
            
            // Then, Master will continue and find the next k partial solutions without waiting for the workers to
            // send their completed solutions
            else
            {
                // Ensuring that Master will only create a new memory block for MPI_Irecv when there is a worker
                // that will sent it's work
                if (flag != 0)
                {
                    MPI_Irecv(&solSize, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &req);
                    flag = 0;
                }
                
                // Master will keep looping until a worker sends their completed solution and Master will receive it
                // and push it into a string vector solutions. Immediately after that, Master will send the next
                // k partial solution to the worker that finishes its work
                MPI_Test(&req, &flag, &status);
                if (flag != 0)
                {
                    
                    // MPI_Irecv counter to counter the total succesfully Receive and Send from Master and Worker
                    irecvCounter++;
                    
                    // If workers cannot find a valid solution, it will send a zero integer for the solSize
                    if (solSize != 0)
                    {
                        // When valid size obtained from workers, the str_solns size will be increased to accomodate
                        // new completed solutions from workers
                        str_solns.resize(str_solns.size() + solSize);
                        
                        // Master receive the completed solution from the worker that sent the size of the solution
                        // and will store it in one long vector str_solns
                        MPI_Recv(&str_solns[str_solns.size() - solSize], solSize, MPI_UNSIGNED, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
                    }
                    
                    // Master will send the next k partial solutions to the worker that just completed its previous work
                    MPI_Send(&part_solns[solCounter][0], n, MPI_UNSIGNED, status.MPI_SOURCE, 111, MPI_COMM_WORLD);
                    solCounter++;   // Counting a successful send of k partial solution to workers
                    flag = -1;
                }
            }
            r++;        // The increment rank of workers. Master will send the first partial work incrementally
            return true;
        }
    }
    // For workers and seq_solver will do this part.
    else
    {
        if (col == n)
        {
            // The completed solution is stored in str_solns vector. If worker has more complete solutions, the result
            // will be appended to str_solns. This will save communication time as workers only have to send one long
            // vector of completed solutions instead of sending multiple completed solutions to Master
            part_solns.resize(0);
            printSol(board, n, part_solns);
            for (int j = 0; j < n; j++)
                str_solns.push_back(part_solns[0][j]);
            return true;
        }
    }
    
    // Consider this column and try placing this queen in all rows one by one
    bool res = false;
    // Iterate over the rows
    for (int i = 0; i < n; i++)
    {
        // Check if queen can be placed on this position
        if (safeCheck(board, i, col, n) )
        {
            // Place the queen in this position
            board[i][col] = 1;
            
            // Make result true if placement is possible, then the solveNQ is called recursively for the next column
            res = solveNQ(board, col + 1, n, k, solCounter, irecvCounter, flag, solSize, part_solns, str_solns, req, status) || res;
            
            // If no valid solutions, then remove the queen in the position and go to the next row
            board[i][col] = 0; // Backtracking
        }
    }
    
    // If queen can not be place in any row in this column col then return false
    return res;
}

/* A utility function to check if a queen can
 be placed on board[row][col]. Note that this
 function is called when "col" queens are
 already placed in columns from 0 to col -1.
 So we need to check only left side for
 attacking queens */
bool safeCheck(std::vector<std::vector<unsigned int> >& board, int row, int col, int n)
{
    int i, j;
    
    /* Check this row on left side */
    for (i = 0; i < col; i++)
        if (board[row][i])
            return false;
    
    /* Check upper diagonal on left side */
    for (i=row, j=col; i>=0 && j>=0; i--, j--)
        if (board[i][j])
            return false;
    
    /* Check lower diagonal on left side */
    for (i=row, j=col; j>=0 && i<n; i++, j--)
        if (board[i][j])
            return false;
    
    return true;
}

// Function to store the valid position of queens
void printSol(std::vector<std::vector<unsigned int> >& board, int n, std::vector<std::vector<unsigned int> >& part_solns)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // For Master, it will store the k partial solutions in a 2d vectors
    if (rank == 0 && size != 1)
    {
        static int m = 0;
        part_solns.push_back(std::vector<unsigned int>(n));
        for (int i = 0; i < n; i++)
        {
            part_solns[m].resize(n);
            for (int j = 0; j < n; j++)
                if (board[j][i] == 1)
                    part_solns[m][i] = j;
        }
        m++;
    }
    // For workers and seq_solver will store the completed solutions in a 1d vector
    else
    {
        part_solns.push_back(std::vector<unsigned int>(n));
        for (int i = 0; i < n; i++)
        {
            part_solns[0].resize(n);
            for (int j = 0; j < n; j++)
                if (board[j][i] == 1)
                    part_solns[0][i] = j;
        }
    }
    
}
