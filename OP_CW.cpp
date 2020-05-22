#include <iostream>
#include <cmath>

using namespace std;

int main()
{
    std::cout << "Hello World!\n";
}

void Jacobi( double** A,  double* B,  double* X,int size,double diff)
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