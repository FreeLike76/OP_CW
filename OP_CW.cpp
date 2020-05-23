#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class LinEquasSys
{
	double** Acoef;
	double* Bcoef;
	double* Result;
	int size;
	void clearABR()
	{
		if (Acoef != nullptr)
		{
			for (int i = 0; i < size; i++)
			{
				delete[] Acoef[i];
			}
			delete[] Acoef;
		}
		if (Bcoef != nullptr)
		{
			delete[] Bcoef;
		}
		if (Result != nullptr)
		{
			delete[] Result;
		}
		size = 0;
	}
	void Jacobi(double diff)
	{
		if (Result == nullptr)
			Result = new double[size];

		for (int i = 0; i < size; i++)
		{
			Result[i] = Bcoef[i];
		}
		double* tempX = new double[size];
		double norm;

		do
		{
			for (int i = 0; i < size; i++)
			{
				tempX[i] = Result[i];
				for (int j = 0; j < size; j++)
				{
					if (i != j)
						tempX[i] -= Acoef[i][j] * Result[j];
				}
				tempX[i] /= Acoef[i][i];
			}
			norm = fabs(Result[0] - tempX[0]);
			for (int i = 0; i < size; i++)
			{
				if (fabs(Result[i] - tempX[i]) > norm)
					norm = fabs(Result[i] - tempX[i]);
				Result[i] = tempX[i];
			}
		} while (norm > diff);
		delete[] tempX;
	}
	void G_Z(double diff)
	{
		if(Result==nullptr)
			Result = new double[size];

		for (int i = 0; i < size; i++)
		{
			Result[i] = Bcoef[i];
		}
		double* tempX = new double[size];
		double norm;

		do {
			for (int i = 0; i < size; i++) {
				tempX[i] = Result[i];						//tempx(newX)=Old X
				for (int j = 0; j < size; j++) {
					if (i != j)
						tempX[i] -= Acoef[i][j] * tempX[j];	// As we go we use the new parameters
				}
				tempX[i] /= Acoef[i][i];
			}
			norm = fabs(Result[0] - tempX[0]);
			for (int i = 0; i < size; i++) {
				if (fabs(Result[i] - tempX[i]) > norm)
					norm = fabs(Result[i] - tempX[i]);
				Result[i] = tempX[i];
			}
		} while (norm > diff);
		delete[] tempX;
	}
public:
	LinEquasSys()
	{
		Acoef = nullptr;
		Bcoef = nullptr;
		Result = nullptr;
		size = 0;
	}
	~LinEquasSys()
	{
		clearABR();
	}
	void setEquasion(double** A,double* B, int cur_size)
	{
		clearABR();
		size = cur_size;
		//									New A and B coef
		Acoef = new double* [size];
		Bcoef = new double[size];
		//									Applying new data
		for (int i = 0; i < size; i++)				
		{
			Acoef[i] = new double[size];
			Bcoef[i] = B[i];
			for (int j = 0; j < size; j++)
			{
				Acoef[i][j] = A[i][j];
			}
		}	  
	}
	void doMath(int method,double accuracy=0.001)
	{
		switch (method)
		{
		case 1:
			Jacobi(accuracy);
			break;
		case 2:
			G_Z(accuracy);
			break;
		case 3:
			//Gradient(accuracy...
			break;
		default:
			break;
		}
	}
	vector<double> getResult()
	{
		vector<double> get;
		if(Result!=nullptr)
		{
			for (int i = 0; i < size; i++)
			{
				get.push_back(Result[i]);
			}
			return get;
		}
		else
		{
			cout << "Error! Results are unavalible, because no calculations has been done!" << endl;
			get.push_back(0);
			return get;
		}
	}
};

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