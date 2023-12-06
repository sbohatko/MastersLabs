// MatrixMultiplication.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "iostream"
#include <stdlib.h>
#include <conio.h>
#include <time.h>

void RandomDataInitialization(double* pAMatrix, double* pBMatrix,
    int Size) {
    int i, j; // Loop variables
    srand(unsigned(clock()));
    for (i = 0; i < Size; i++)
        for (j = 0; j < Size; j++) {
            pAMatrix[i * Size + j] = rand() / double(1000);
            pBMatrix[i * Size + j] = rand() / double(1000);
        }
}

void ProcessInitialization(double*& pAMatrix, double*& pBMatrix,
    double*& pCMatrix, int& Size) {
    // Setting the size of matrices
    do {
        printf("\nEnter the size of matrices: ");
        std::cin >> Size;
        printf("\nChosen matrices' size = %d\n", Size);
        if (Size <= 0)
            printf("\nSize of objects must be greater than 0!\n");
    } while (Size <= 0);
    // Memory allocation
    pAMatrix = new double[Size * Size];
    pBMatrix = new double[Size * Size];
    pCMatrix = new double[Size * Size];
    // Initialization of matrix elements
    RandomDataInitialization(pAMatrix, pBMatrix, Size);
    for (int i = 0; i < Size * Size; i++) {
        pCMatrix[i] = 0;
    }
}

void SerialResultCalculation(double* pAMatrix, double* pBMatrix,
    double* pCMatrix, int Size) {
    int i, j, k; // Loop variables
    for (i = 0; i < Size; i++) {
        for (j = 0; j < Size; j++)
            for (k = 0; k < Size; k++)
                pCMatrix[i * Size + j] += pAMatrix[i * Size + k] * pBMatrix[k * Size + j];
    }
}

void ProcessTermination(double* pAMatrix, double* pBMatrix,
    double* pCMatrix) {
    delete[] pAMatrix;
    delete[] pBMatrix;
    delete[] pCMatrix;
}

int main()
{
    double* pAMatrix; // First argument of matrix multiplication
    double* pBMatrix; // Second argument of matrix multiplication
    double* pCMatrix; // Result matrix
    int Size; // Size of matrices
    time_t start, finish;
    double duration;
    printf("Serial matrix multiplication program\n");
    // Memory allocation and initialization of matrix elements
    ProcessInitialization(pAMatrix, pBMatrix, pCMatrix, Size);
    // Matrix multiplication
    start = clock();
    SerialResultCalculation(pAMatrix, pBMatrix, pCMatrix, Size);
    finish = clock();
    duration = (finish - start) / double(CLOCKS_PER_SEC);
    // Printing the time spent by matrix multiplication
    printf("\n Time of execution: %f\n", duration);
    // Computational process termination
    ProcessTermination(pAMatrix, pBMatrix, pCMatrix);
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
