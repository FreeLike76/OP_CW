#include <iostream>
#include <cmath>

using namespace std;

double* Jacobi(double** A, double* B, double* X, int size, double diff);
int** getNullM(int rows, int cols);
void deleteM(int** Matrx, int rows);

int main()
{
	int size;
	double diff;
}

int** getNullM(int rows, int cols)
{
	int** res = new int* [rows];
	for (int i = 0; i < rows; i++)
	{
		res[i] = new int[cols];
		for (int j = 0; j < cols; j++)
		{
			res[i][j] = 0;
		}
	}
	return res;
}
void deleteM(int** Matrx, int rows)
{
	for (int i = 0; i < rows; i++)
	{
		delete[] Matrx[i];
	}
	delete[] Matrx;
}

double* Jacobi( double** A,  double* B,  double* X,int size,double diff) //Need to test with double X as B
{
	double* tempX = new double[size];
	double norm;

	do {
		for (int i = 0; i < size; i++) {
			tempX[i] = B[i];
			for (int j = 0; j < size; j++) {
				if (i != j)
					tempX[i] -= A[i][j] * X[j];
			}
			tempX[i] /= A[i][i];
		}
		norm = fabs(X[0] - tempX[0]);
		for (int i = 0; i < size; i++) {
			if (fabs(X[i] - tempX[i]) > norm)
				norm = fabs(X[i] - tempX[i]);
			X[i] = tempX[i];
		}
	} while (norm > diff);
	delete[] tempX;
}
double* G_Z(double** A, double* B, double* X, int size, double diff) //Need to test with double X as B
{
	double* tempX = new double[size];
	double norm;

	do {
		for (int i = 0; i < size; i++) {
			tempX[i] = X[i];						//tempx(newX)=Old X
			for (int j = 0; j < size; j++) {
				if (i != j)
					tempX[i] -= A[i][j] * tempX[j];	// As we go we use the new parameters
			}
			tempX[i] /= A[i][i];
		}
		norm = fabs(X[0] - tempX[0]);
		for (int i = 0; i < size; i++) {
			if (fabs(X[i] - tempX[i]) > norm)
				norm = fabs(X[i] - tempX[i]);
			X[i] = tempX[i];
		}
	} while (norm > diff);
	delete[] tempX;
}