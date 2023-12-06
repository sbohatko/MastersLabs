#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <iostream>

int ProcNum;
int ProcRank;
int* pParallelPivotPos;
int* pProcPivotIter;
int* pProcInd;
int* pProcNum;

void RandomDataInitialization(double* pMatrix, double* pVector, int Size)
{
    srand(unsigned(clock()));
    for (int i = 0; i < Size; i++)
    {
        pVector[i] = rand() / double(1000);
        for (int j = 0; j < Size; j++)
        {
            if (j <= i)
            {
                pMatrix[i * Size + j] = rand() / double(1000);
            }
            else
            {
                pMatrix[i * Size + j] = 0;
            }
        }
    }
}

void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult, double*& pProcRows, double*& pProcVector, double*& pProcResult, int& Size, int& RowNum)
{
    int RestRows;
    if (ProcRank == 0)
    {
        do
        {
            printf("\nEnter the size of the matrix and the vector: ");
            std::cin >> Size;
            if (Size < ProcNum)
            {
                printf("Size must be greater than the number of processes! \n");
            }
        } while (Size < ProcNum);
    }
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Corrected the calculation of RowNum
    int localRows = Size / ProcNum;
    int extraRows = Size % ProcNum;

    RowNum = (ProcRank < extraRows) ? localRows + 1 : localRows;

    RestRows = Size - (ProcRank * localRows + std::min(ProcRank, extraRows) * (localRows + 1));

    pProcRows = new double[RowNum * Size];
    pProcVector = new double[RowNum];
    pProcResult = new double[RowNum];

    pParallelPivotPos = new int[Size];
    pProcPivotIter = new int[RowNum];

    pProcInd = new int[ProcNum];
    pProcNum = new int[ProcNum];

    for (int i = 0; i < RowNum; i++)
    {
        pProcPivotIter[i] = -1;
    }
    if (ProcRank == 0)
    {
        pMatrix = new double[Size * Size];
        pVector = new double[Size];
        pResult = new double[Size];
        RandomDataInitialization(pMatrix, pVector, Size);
    }
}


void DataDistribution(double* pMatrix, double* pProcRows, double* pVector, double* pProcVector, int Size, int RowNum)
{
    int* pSendNum; // Number of the elements sent to the process
    int* pSendInd; // Index of the first data element sent to the process
    int RestRows = Size; // Number of rows that have not been distributed yet
    int i; // Loop variable
    // Alloc memory for temporary objects
    pSendInd = new int[ProcNum];
    pSendNum = new int[ProcNum];
    // Define the disposition of the matrix rows for the current process
    RowNum = (Size / ProcNum);
    pSendNum[0] = RowNum * Size;
    pSendInd[0] = 0;
    for (i = 1; i < ProcNum; i++)
    {
        RestRows -= RowNum;
        RowNum = RestRows / (ProcNum - i);
        pSendNum[i] = RowNum * Size;
        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
    }
    // Scatter the rows
    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows,
        pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Define the disposition of the matrix rows for the current process
    RestRows = Size;
    pProcInd[0] = 0;
    pProcNum[0] = Size / ProcNum;
    for (i = 1; i < ProcNum; i++)
    {
        RestRows -= pProcNum[i - 1];
        pProcNum[i] = RestRows / (ProcNum - i);
        pProcInd[i] = pProcInd[i - 1] + pProcNum[i - 1];
    }
    MPI_Scatterv(pVector, pProcNum, pProcInd, MPI_DOUBLE, pProcVector, pProcNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    delete[] pSendNum;
    delete[] pSendInd;
}

void ResultCollection(double* pProcResult, double* pResult)
{
    MPI_Gatherv(pProcResult, pProcNum[ProcRank], MPI_DOUBLE, pResult, pProcNum, pProcInd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void ParallelEliminateColumns(double* pProcRows, double* pProcVector, double* pPivotRow, int Size, int RowNum, int Iter)
{
    double multiplier;
    for (int i = 0; i < RowNum; i++)
    {
        if (pProcPivotIter[i] == -1)
        {
            multiplier = pProcRows[i * Size + Iter] / pPivotRow[Iter];
            for (int j = Iter; j < Size; j++)
            {
                pProcRows[i * Size + j] -= pPivotRow[j] * multiplier;
            }
            pProcVector[i] -= pPivotRow[Size] * multiplier;
        }
    }
}

void ParallelGaussianElimination(double* pProcRows, double* pProcVector, int Size, int RowNum)
{
    double MaxValue; // Value of the pivot element of th� process
    int PivotPos;
    struct { double MaxValue; int ProcRank; } ProcPivot, Pivot;
    // pPivotRow is used for storing the pivot row and the corresponding
    // element of the vector b
    double* pPivotRow = new double[Size + 1];
    // The iterations of the Gaussian elimination stage
    for (int i = 0; i < Size; i++)
    {
        // Calculating the local pivot row
        double MaxValue = 0;
        for (int j = 0; j < RowNum; j++)
        {
            if ((pProcPivotIter[j] == -1) && (MaxValue < fabs(pProcRows[j * Size + i])))
            {
                MaxValue = fabs(pProcRows[j * Size + i]);
                PivotPos = j;
            }
        }
        ProcPivot.MaxValue = MaxValue;
        ProcPivot.ProcRank = ProcRank;
        // Finding the pivot process
        // (process with the maximum value of MaxValue)
        MPI_Allreduce(&ProcPivot, &Pivot, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        // Broadcasting the pivot row
        if (ProcRank == Pivot.ProcRank)
        {
            pProcPivotIter[PivotPos] = i; //iteration number
            pParallelPivotPos[i] = pProcInd[ProcRank] + PivotPos;
        }
        MPI_Bcast(&pParallelPivotPos[i], 1, MPI_INT, Pivot.ProcRank,
            MPI_COMM_WORLD);
        if (ProcRank == Pivot.ProcRank)
        {
            // Fill the pivot row
            for (int j = 0; j < Size; j++)
            {
                pPivotRow[j] = pProcRows[PivotPos * Size + j];
            }
            pPivotRow[Size] = pProcVector[PivotPos];
        }
        MPI_Bcast(pPivotRow, Size + 1, MPI_DOUBLE, Pivot.ProcRank, MPI_COMM_WORLD);
        ParallelEliminateColumns(pProcRows, pProcVector, pPivotRow, Size, RowNum, i);
    }
}

void FindBackPivotRow(int RowIndex, int& IterProcRank, int& IterPivotPos)
{
    for (int i = 0; i < ProcNum - 1; i++)
    {
        if ((pProcInd[i] <= RowIndex) && (RowIndex < pProcInd[i + 1]))
            IterProcRank = i;
    }
    if (RowIndex >= pProcInd[ProcNum - 1])
        IterProcRank = ProcNum - 1;
    IterPivotPos = RowIndex - pProcInd[IterProcRank];
}

void ParallelBackSubstitution(double* pProcRows, double* pProcVector, double* pProcResult, int Size, int RowNum)
{
    int IterProcRank; // Rank of the process with the current pivot row
    int IterPivotPos; // Position of the pivot row of the process
    double IterResult; // Calculated value of the current unknown
    double val;
    // Iterations of the back substitution stage
    for (int i = Size - 1; i >= 0; i--)
    {
        // Calculating the rank of the process, which holds the pivot row
        FindBackPivotRow(pParallelPivotPos[i], IterProcRank, IterPivotPos);
        // Calculating the unknown
        if (ProcRank == IterProcRank)
        {
            IterResult = pProcVector[IterPivotPos] / pProcRows[IterPivotPos * Size + i];
            pProcResult[IterPivotPos] = IterResult;
        }
        // Broadcasting the value of the current unknown
        MPI_Bcast(&IterResult, 1, MPI_DOUBLE, IterProcRank, MPI_COMM_WORLD);
        // Updating the values of the vector b
        for (int j = 0; j < RowNum; j++)
            if (pProcPivotIter[j] < i)
            {
                val = pProcRows[j * Size + i] * IterResult;
                pProcVector[j] = pProcVector[j] - val;
            }
    }
}

void ParallelResultCalculation(double* pProcRows, double* pProcVector, double* pProcResult, int Size, int RowNum)
{
    ParallelGaussianElimination(pProcRows, pProcVector, Size, RowNum);
    ParallelBackSubstitution(pProcRows, pProcVector, pProcResult, Size, RowNum);
}

void ProcessTermination(double* pMatrix, double* pVector, double* pResult, double* pProcRows, double* pProcVector, double* pProcResult)
{
    if (ProcRank == 0)
    {
        delete[] pMatrix;
        delete[] pVector;
        delete[] pResult;
    }
    delete[] pProcRows;
    delete[] pProcVector;
    delete[] pProcResult;

    delete[] pParallelPivotPos;
    delete[] pProcPivotIter;

    delete[] pProcInd;
    delete[] pProcNum;
}

int main(int argc, char* argv[])
{
    double* pMatrix; // Matrix of the linear system
    double* pVector; // Right parts of the linear system
    double* pResult; // Result vector
    double* pProcRows; // Rows of the matrix A
    double* pProcVector; // Elements of the vector b
    double* pProcResult; // Elements of the vector x
    int Size; // Sizes of the matrix and the vectors
    int RowNum; // Number of the matrix rows
    double start, finish, duration;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    if (ProcRank == 0)
        printf("Parallel Gauss algorithm for solving linear systems\n");
    // Memory allocation and data initialization
    ProcessInitialization(pMatrix, pVector, pResult, pProcRows, pProcVector, pProcResult, Size, RowNum);
    // The execution of the parallel Gauss algorithm
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    DataDistribution(pMatrix, pProcRows, pVector, pProcVector, Size, RowNum);
    ParallelResultCalculation(pProcRows, pProcVector, pProcResult, Size, RowNum);
    ResultCollection(pProcResult, pResult);
    finish = MPI_Wtime();
    duration = finish - start;
    // Printing the time spent by the Gauss algorithm
    if (ProcRank == 0)
    {
        printf("\n Time of execution: %f\n", duration);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Computational process termination
    ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcVector, pProcResult);
    MPI_Finalize();
    return 0;
}
