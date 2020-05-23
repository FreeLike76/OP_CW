#include <iostream>
#include <cmath>

using namespace std;

double* Jacobi(double** A, double* B, int size, double diff);
double* G_Z(double** A, double* B, int size, double diff);
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
bool zeroIJ(int** Matrx,int size)
{
	for (int i = 0; i < size; i++)
	{
		if (Matrx[i][i] == 0)
			return true;
	}
	return false;
}
void zeroFix(int** A,int* B,int size)
{
	for (int iter = 0; iter < size; iter++)			//Run through matrix
	{
		if(A[iter][iter]==0)						//If there is 0  on [i][i]
		{
			for (int i = 0; i < size; i++)			//Look for another row
			{
				if(A[i][iter]!=0)					//That does not have zero on required position
				{
					for (int j = 0; j < size; j++)	//Add it's coef to initial row
					{
						A[iter][j] += A[i][j];
					}
					B[iter] += B[i];				//And the ...=B(i)
					break;							//End
				}
			}
		}
	}
}


double* Jacobi( double** A,  double* B,int size,double diff) //Need to test with double X as B
{
	double* X = new double[size];
	for (int i = 0; i < size; i++)
	{
		X[i] = B[i];
	}
	double* tempX = new double[size];
	double norm;

	do 
	{
		for (int i = 0; i < size; i++) 
		{
			tempX[i] = X[i];
			for (int j = 0; j < size; j++) 
			{
				if (i != j)
					tempX[i] -= A[i][j] * X[j];
			}
			tempX[i] /= A[i][i];
		}
		norm = fabs(X[0] - tempX[0]);
		for (int i = 0; i < size; i++) 
		{
			if (fabs(X[i] - tempX[i]) > norm)
				norm = fabs(X[i] - tempX[i]);
			X[i] = tempX[i];
		}
	} while (norm > diff);
	delete[] tempX;
	return X;
}
double* G_Z(double** A, double* B,  int size, double diff) //Need to test with double X as B
{
	double* X = new double[size];
	for (int i = 0; i < size; i++)
	{
		X[i] = B[i];
	}
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
	return X;
}