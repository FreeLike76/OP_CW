#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

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
	void Gradient(){}
	bool zeroIJ()
	{
		for (int i = 0; i < size; i++)
		{
			if (Acoef[i][i] == 0)
				return true;
		}
		return false;
	}
	void zeroFix()
	{
		for (int iter = 0; iter < size; iter++)						//Run through matrix
		{
			if (Acoef[iter][iter] == 0)								//If there is 0  on [i][i]
			{
				for (int i = 0; i < size; i++)						//Look for another row
				{
					if (Acoef[i][iter] != 0)						//That does not have zero on required position
					{
						for (int j = 0; j < size; j++)				//Add it's coef to initial row
						{
							Acoef[iter][j] += Acoef[i][j];
						}
						Bcoef[iter] += Bcoef[i];					//And the ...=B(i)
						break;										//End
					}
				}
			}
		}
	}
	void getMatrixWithoutRowAndCol(double** matrix, int _size, int row, int col, double** newMatrix) {
		int offsetRow = 0; //Смещение индекса строки в матрице
		int offsetCol = 0; //Смещение индекса столбца в матрице
		for (int i = 0; i < _size - 1; i++) {
			//Пропустить row-ую строку
			if (i == row) {
				offsetRow = 1; //Как только встретили строку, которую надо пропустить, делаем смещение для исходной матрицы
			}

			offsetCol = 0; //Обнулить смещение столбца
			for (int j = 0; j < _size - 1; j++) {
				//Пропустить col-ый столбец
				if (j == col) {
					offsetCol = 1; //Встретили нужный столбец, проускаем его смещением
				}

				newMatrix[i][j] = matrix[i + offsetRow][j + offsetCol];
			}
		}
	}
	double matrixDet(double** matrix, int _size) {
		double det = 0;
		int degree = 1; // (-1)^(1+j) из формулы определителя

		//Условие выхода из рекурсии
		if (_size == 1) {
			return matrix[0][0];
		}
		//Условие выхода из рекурсии
		else if (_size == 2) {
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
		}
		else {
			//Матрица без строки и столбца
			double** newMatrix = new double* [_size - 1];
			for (int i = 0; i < _size - 1; i++) {
				newMatrix[i] = new double[_size - 1];
			}

			//Раскладываем по 0-ой строке, цикл бежит по столбцам
			for (int j = 0; j < _size; j++) {
				//Удалить из матрицы i-ю строку и j-ый столбец
				//Результат в newMatrix
				getMatrixWithoutRowAndCol(matrix, _size, 0, j, newMatrix);

				//Рекурсивный вызов
				//По формуле: сумма по j, (-1)^(1+j) * matrix[0][j] * minor_j (это и есть сумма из формулы)
				//где minor_j - дополнительный минор элемента matrix[0][j]
				// (напомню, что минор это определитель матрицы без 0-ой строки и j-го столбца)
				det = det + (degree * matrix[0][j] * matrixDet(newMatrix, _size - 1));
				//"Накручиваем" степень множителя
				degree = -degree;
			}

			//Чистим память на каждом шаге рекурсии(важно!)
			for (int i = 0; i < _size - 1; i++) {
				delete[] newMatrix[i];
			}
			delete[] newMatrix;
		}
		return det;
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
			//WRONG
			break;
		}
	}
	vector<double> getvectorResult()
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
	void CheckFixDiagonals()
	{
		if (zeroIJ())
			zeroFix();
	}
	double getDet()
	{
		return matrixDet(Acoef, size);
	}
	bool is_semetrical()
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (Acoef[i][j] != Acoef[j][i])
					return false;
			}
		}
		return true;
	}
	bool is_positive()
	{
		for (int i = 1; i <= size; i++)
		{
			if (matrixDet(Acoef, i) <= 0)
				return false;
		}
		return true;
	}
};

int main()
{
	int size;
	double diff;
	LinEquasSys test;
}
bool readLESff(LinEquasSys cur,string path)
{
	ifstream input(path);
	if (!input.is_open())
		return false;
	else
	{
		//Reading size
		int _size;
		input >> _size;
		//Creating A and B
		double ** A= new double* [_size];
		for (int i = 0; i < _size; i++)
		{
			A[i] = new double[_size];
		}
		double* B = new double[_size];
		//Reading A and B
		for (int i = 0; i < _size; i++)
		{
			for (int j = 0; j < _size; j++)
			{
				input >> A[i][j];
			}
			input >> B[i];
		}
		cur.setEquasion(A, B, _size);
		input.close();
	}


}