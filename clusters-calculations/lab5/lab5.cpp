//#include <cstdlib>
//#include <cstdio>
//#include <ctime>
//#include <algorithm>
//#include <mpi.h>
//#include <iostream>
//
//using namespace std;
//int ProcRank; // Rank of current process
//int ProcNum; // Number of processes
//const double InfinitiesPercent = 50.0;
//const double RandomDataMultiplier = 10;
//
//
//int Min(int A, int B) {
//	int Result = (A < B) ? A : B;
//	if ((A < 0) && (B >= 0)) Result = B;
//	if ((B < 0) && (A >= 0)) Result = A;
//	if ((A < 0) && (B < 0)) Result = -1;
//	return Result;
//}
//
//void CopyMatrix(int* pMatrix, int Size, int* pMatrixCopy) {
//	copy(pMatrix, pMatrix + Size * Size, pMatrixCopy);
//}
//
//// Function for comparing the matrices
//bool CompareMatrices(int* pMatrix1, int* pMatrix2, int Size) {
//	return equal(pMatrix1, pMatrix1 + Size * Size, pMatrix2);
//}
//
//// Function for the serial Floyd algorithm
//void SerialFloyd(int* pMatrix, int Size) {
//	int t1, t2;
//	for (int k = 0; k < Size; k++)
//		for (int i = 0; i < Size; i++)
//			for (int j = 0; j < Size; j++)
//				if ((pMatrix[i * Size + k] != -1) &&
//					(pMatrix[k * Size + j] != -1)) {
//					t1 = pMatrix[i * Size + j];
//					t2 = pMatrix[i * Size + k] + pMatrix[k * Size + j];
//					pMatrix[i * Size + j] = Min(t1, t2);
//				}
//}
//
//// Function for formatted matrix output
//void PrintMatrix(int* pMatrix, int RowCount, int ColCount) {
//	for (int i = 0; i < RowCount; i++) {
//		for (int j = 0; j < ColCount; j++) {
//			printf("%7d", pMatrix[i * ColCount + j]);
//			fflush(stdout);
//		}
//		printf("\n");
//		fflush(stdout);
//	}
//}
//
//
//// Function for computational process termination
//void ProcessTermination(int* pMatrix, int* pProcRows) {
//	if (ProcRank == 0)
//		delete[]pMatrix;
//	delete[]pProcRows;
//}
//
//// Function for simple setting the initial data
//void DummyDataInitialization(int* pMatrix, int Size) {
//	for (int i = 0; i < Size; i++)
//		for (int j = i; j < Size; j++) {
//			if (i == j) pMatrix[i * Size + j] = 0;
//			else
//				if (i == 0) pMatrix[i * Size + j] = j;
//				else pMatrix[i * Size + j] = -1;
//			pMatrix[j * Size + i] = pMatrix[i * Size + j];
//		}
//}
//
//// Function for setting the data by the random generator
//void RandomDataInitialization(int* pMatrix, int Size) {
//	srand((unsigned)time(0));
//	for (int i = 0; i < Size; i++)
//		for (int j = 0; j < Size; j++)
//			if (i != j) {
//				if ((rand() % 100) < InfinitiesPercent)
//					pMatrix[i * Size + j] = -1;
//				else
//					pMatrix[i * Size + j] = rand() + 1;
//			}
//			else
//				pMatrix[i * Size + j] = 0;
//}
//
//// Function for the data distribution among the processes
//void DataDistribution(int* pMatrix, int* pProcRows, int Size, int RowNum) {
//	int* pSendNum; // Number of elements sent to the process
//	int* pSendInd; // Index of the first data element sent to the process
//	int RestRows = Size; // Number of rows, that havenÒt been distributed yet
//	// Allocate memory for temporary objects
//	pSendInd = new int[ProcNum];
//	pSendNum = new int[ProcNum];
//	// Define the disposition of the matrix rows for current process
//	RowNum = Size / ProcNum;
//	pSendNum[0] = RowNum * Size;
//	pSendInd[0] = 0;
//	for (int i = 1; i < ProcNum; i++) {
//		RestRows -= RowNum;
//		RowNum = RestRows / (ProcNum - i);
//		pSendNum[i] = RowNum * Size;
//		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
//	}
//	// Scatter the rows
//	MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_INT,
//		pProcRows, pSendNum[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);
//	// Free allocated memory
//	delete[]pSendNum;
//	delete[]pSendInd;
//}
//
//// Function for allocating the memory and setting the initial values
//void ProcessInitialization(int*& pMatrix, int*& pProcRows, int& Size,
//	int& RowNum) {
//	setvbuf(stdout, 0, _IONBF, 0);
//	if (ProcRank == 0) {
//		do {
//			printf("Enter the number of vertices: ");
//			std::cin >> Size;
//			if (Size < ProcNum)
//				printf("The number of vertices should be greater then"
//					"the number of processes\n");
//		} while (Size < ProcNum);
//		printf("Using the graph with %d vertices\n", Size);
//	}
//	// Broadcast the number of vertices
//	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	// Number of rows for each process
//	int RestRows = Size;
//	for (int i = 0; i < ProcRank; i++)
//		RestRows = RestRows - RestRows / (ProcNum - i);
//	RowNum = RestRows / (ProcNum - ProcRank);
//	// Allocate memory for the current process rows
//	pProcRows = new int[Size * RowNum];
//	if (ProcRank == 0) {
//		// Allocate memory for the adjacency matrix
//		pMatrix = new int[Size * Size];
//		// Data initalization
//		DummyDataInitialization(pMatrix, Size);
//		//RandomDataInitialization(pMatrix, Size);
//	}
//}
//
//// Function for process result collection
//void ResultCollection(int* pMatrix, int* pProcRows, int Size, int RowNum) {
//	int* pReceiveNum; // Number of elements, that current process sends
//	int* pReceiveInd; // Offset for storing the data from current process
//	int RestRows = Size; // Number of rows, that haven't been gathered yet
//	// Allocate memory for temporary objects
//	pReceiveNum = new int[ProcNum];
//	pReceiveInd = new int[ProcNum];
//	// Determine the disposition of the result data block of current process
//	RowNum = Size / ProcNum;
//	pReceiveInd[0] = 0;
//	pReceiveNum[0] = RowNum * Size;
//	for (int i = 1; i < ProcNum; i++) {
//		RestRows -= RowNum;
//		RowNum = RestRows / (ProcNum - i);
//		pReceiveNum[i] = RowNum * Size;
//		pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
//	}
//	// Gather the whole matrix on the process 0
//	MPI_Gatherv(pProcRows, pReceiveNum[ProcRank], MPI_INT,
//		pMatrix, pReceiveNum, pReceiveInd, MPI_INT, 0, MPI_COMM_WORLD);
//	// Free allocated memory
//	delete[]pReceiveNum;
//	delete[]pReceiveInd;
//}
//
//// Function for row broadcasting among all processes
//void RowDistribution(int* pProcRows, int Size, int RowNum, int k, int
//	* pRow) {
//	int ProcRowRank; // Process rank with the row k
//	int ProcRowNum; // Process row number
//	// Finding the process rank with the row k
//	int RestRows = Size;
//	int Ind = 0;
//	int Num = Size / ProcNum;
//	for (ProcRowRank = 1; ProcRowRank < ProcNum + 1; ProcRowRank++) {
//		if (k < Ind + Num) break;
//		RestRows -= Num;
//		Ind += Num;
//		Num = RestRows / (ProcNum - ProcRowRank);
//	}
//	ProcRowRank = ProcRowRank - 1;
//	ProcRowNum = k - Ind;
//	if (ProcRowRank == ProcRank)
//		// Copy the row to pRow array
//		copy(&pProcRows[ProcRowNum * Size], &pProcRows[(ProcRowNum + 1) * Size], pRow);
//	// Broadcast row to all processes
//	MPI_Bcast(pRow, Size, MPI_INT, ProcRowRank, MPI_COMM_WORLD);
//}
//
//// Function for the parallel Floyd algorithm
//void ParallelFloyd(int* pProcRows, int Size, int RowNum) {
//	int* pRow = new int[Size];
//	int t1, t2;
//	for (int k = 0; k < Size; k++) {
//		// Distribute row among all processes
//		RowDistribution(pProcRows, Size, RowNum, k, pRow);
//		// Update adjacency matrix elements
//		for (int i = 0; i < RowNum; i++)
//			for (int j = 0; j < Size; j++)
//				if ((pProcRows[i * Size + k] != -1) &&
//					(pRow[j] != -1)) {
//					t1 = pProcRows[i * Size + j];
//					t2 = pProcRows[i * Size + k] + pRow[j];
//					pProcRows[i * Size + j] = Min(t1, t2);
//				}
//	}
//	delete[]pRow;
//}
//
//// Function for formatted output of all stripes
//void ParallelPrintMatrix(int* pProcRows, int Size, int RowNum) {
//	for (int i = 0; i < ProcNum; i++) {
//		MPI_Barrier(MPI_COMM_WORLD);
//		if (ProcRank == i) {
//			printf("ProcRank = %d\n", ProcRank);
//			fflush(stdout);
//			printf("Proc rows:\n");
//			fflush(stdout);
//			PrintMatrix(pProcRows, RowNum, Size);
//			fflush(stdout);
//		}
//		MPI_Barrier(MPI_COMM_WORLD);
//	}
//}
//
//// Function for testing the data distribution
//void TestDistribution(int* pMatrix, int* pProcRows, int Size, int RowNum) {
//	MPI_Barrier(MPI_COMM_WORLD);
//	if (ProcRank == 0) {
//		printf("Initial adjacency matrix:\n");
//		PrintMatrix(pMatrix, Size, Size);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	ParallelPrintMatrix(pProcRows, Size, RowNum);
//}
//
//// Testing the result of parallel Floyd algorithm
//void TestResult(int* pMatrix, int* pSerialMatrix, int Size) {
//	MPI_Barrier(MPI_COMM_WORLD);
//	if (ProcRank == 0) {
//		SerialFloyd(pSerialMatrix, Size);
//		if (!CompareMatrices(pMatrix, pSerialMatrix, Size)) {
//			printf("Results of serial and parallel algorithms are "
//				"NOT identical. Check your code\n");
//		}
//		else {
//			printf("Results of serial and parallel algorithms are "
//				"identical\n");
//		}
//	}
//}
//
//int main(int argc, char* argv[]) {
//	int* pMatrix; // Adjacency matrix
//	int Size; // Size of adjacency matrix
//	int* pProcRows; // Process rows
//	int RowNum; // Number of process rows
//	double start, finish;
//	double duration = 0.0;
//	int* pSerialMatrix = 0;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//	if (ProcRank == 0)
//		printf("Parallel Floyd algorithm\n");
//	// Process initialization
//	ProcessInitialization(pMatrix, pProcRows, Size, RowNum);
//	if (ProcRank == 0) {
//		// Matrix copying
//		pSerialMatrix = new int[Size * Size];
//		CopyMatrix(pMatrix, Size, pSerialMatrix);
//	}
//	start = MPI_Wtime();
//	// Distributing the initial data between processes
//	DataDistribution(pMatrix, pProcRows, Size, RowNum);
//	// Testing the distribution
//	//TestDistribution(pMatrix, pProcRows, Size, RowNum);
//	// Parallel Floyd algorithm
//	ParallelFloyd(pProcRows, Size, RowNum);
//	//ParallelPrintMatrix(pProcRows, Size, RowNum);
//	// Process data collection
//	ResultCollection(pMatrix, pProcRows, Size, RowNum);
//	//if(ProcRank == 0)
//	// PrintMatrix(pMatrix, Size, Size);
//	finish = MPI_Wtime();
//	//TestResult(pMatrix, pSerialMatrix, Size);
//	duration = finish - start;
//	if (ProcRank == 0)
//		printf("Time of execution: %f\n", duration);
//	if (ProcRank == 0)
//		delete[]pSerialMatrix;
//	// Process termination
//	ProcessTermination(pMatrix, pProcRows);
//	MPI_Finalize();
//	return 0;
//}
//


//#include <stdio.h>
//#include <stdlib.h>
//#include "mpi.h"
//#include <math.h>
//#define SEED 35791246
//
//int main(int argc, char* argv[])
//{
//    long niter = 100000;
//    int myid;                       //holds process's rank id
//    double x, y;                     //x,y value for the random coordinate
//    int i, count = 0;                 //Count holds all the number of how many good coordinates
//    double z;                       //Used to check if x^2+y^2<=1
//    double pi;                      //holds approx value of pi
//    int nodenum;
//
//    MPI_Init(&argc, &argv);                 //Start MPI
//    MPI_Comm_rank(MPI_COMM_WORLD, &myid);           //get rank of node's process
//    MPI_Comm_size(MPI_COMM_WORLD, &nodenum);
//    int* recieved = new int[nodenum-1];
//    long* recvniter = new long[nodenum-1];
//    srand(SEED + myid);                       //Give rand() a seed value. Needs to be different on each node
//
//    if (myid != 0)
//    {
//        for (i = 0; i < niter; ++i)                  //main loop
//        {
//            x = ((double)rand()) / RAND_MAX;           //gets a random x coordinate
//            y = ((double)rand()) / RAND_MAX;           //gets a random y coordinate
//            z = sqrt(x * x + y * y);                  //Checks to see if number in inside unit circle
//            if (z <= 1)
//            {
//                count++;                //if it is, consider it a valid random point
//            }
//        }
//        MPI_Send(&count,
//            1,
//            MPI_INT,
//            0,
//            1,
//            MPI_COMM_WORLD);
//        MPI_Send(&niter,
//            1,
//            MPI_LONG,
//            0,
//            2,
//            MPI_COMM_WORLD);
//    }
//
//    if (myid == 0)
//    {
//        for (i = 0; i < nodenum-1; ++i)
//        {
//            MPI_Recv(&recieved[i],
//                1,
//                MPI_INT,
//                MPI_ANY_SOURCE,
//                1,
//                MPI_COMM_WORLD,
//                MPI_STATUS_IGNORE);
//            MPI_Recv(&recvniter[i],
//                1,
//                MPI_LONG,
//                MPI_ANY_SOURCE,
//                2,
//                MPI_COMM_WORLD,
//                MPI_STATUS_IGNORE);
//        }
//
//        int finalcount = 0;
//        long finalniter = 0;
//        for (i = 0; i < nodenum-1; ++i)
//        {
//            finalcount += recieved[i];
//            finalniter += recvniter[i];
//        }
//
//        pi = ((double)finalcount / (double)finalniter) * 4.0;               //p = 4(m/n)
//        printf("Pi: %f\nfinalcount: %i\nfinalniter: %u\n", pi, finalcount, finalniter);             //Print the calculated value of pi
//
//    }
//
//    MPI_Finalize();                     //Close the MPI instance
//    return 0;
//}

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <time.h>

int main(int argc, char* argv[])
{
    int niter = 1000000;
    int myid;                       //holds process's rank id
    double x, y;                     //x,y value for the random coordinate
    int i;
    int count = 0;                //Count holds all the number of how many good coordinates
    double z;                       //Used to check if x^2+y^2<=1
    double pi;                      //holds approx value of pi
    int reducedcount;                   //total number of "good" points from all nodes
    int reducedniter;                   //total number of ALL points from all nodes

    MPI_Init(&argc, &argv);                 //Start MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);           //get rank of node's process

    /* Everyone can do the work */
    srand(time(NULL) + myid);                  //Give rand() a seed value unique on each node (times are synced)

    for (i = 0; i < niter; ++i)                  //main loop
    {
        x = (double)rand() / RAND_MAX;          //gets a random x coordinate
        y = (double)rand() / RAND_MAX;          //gets a random y coordinate
        z = sqrt((x * x) + (y * y));              //Checks to see if number in inside unit circle
        if (z <= 1)
        {
            ++count;                //if it is, consider it a valid random point
        }
    }


    MPI_Reduce(&count,
        &reducedcount,
        1,
        MPI_INT,
        MPI_SUM,
        0,
        MPI_COMM_WORLD);
    MPI_Reduce(&niter,
        &reducedniter,
        1,
        MPI_INT,
        MPI_SUM,
        0,
        MPI_COMM_WORLD);

    if (myid == 0)                      //if root process
    {
        pi = ((double)reducedcount / (double)reducedniter) * 4.0;               //p = 4(m/n)
        printf("Pi: %f\n%i\n%d\n", pi, reducedcount, reducedniter);
        //Print the calculated value of pi

    }

    MPI_Finalize();                     //Close the MPI instance
    return 0;
}