#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#define _PATH "input.txt"

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
				tempX[i] = Bcoef[i];
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
				tempX[i] = Bcoef[i];
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
	void conj_grad(double diff)
	{
		if (Result == nullptr)
			Result = new double[size];	//Result== Xk	(original)

		double* Zk = new double[size];
		double* Rk = new double[size];
		double* Sz = new double[size];
		double Spr, Spr1, Spz, alpha, beta, mf;
		int i, j, kl = 1;
		double max_iter = 100000;

		/* Вычисляем сумму квадратов элементов вектора F*/
		for (mf = 0, i = 0; i < size; i++) {
			mf += Bcoef[i] * Bcoef[i];
		}

		/* Задаем начальное приближение корней. В Result хранятся значения корней
		 * к-й итерации. */
		for (i = 0; i < size; i++) {
			Result[i] = 0.2;
		}

		/* Задаем начальное значение r0 и z0. */
		for (i = 0; i < size; i++) {
			for (Sz[i] = 0, j = 0; j < size; j++)
				Sz[i] += Acoef[i][j] * Result[j];
			Rk[i] = Bcoef[i] - Sz[i];
			Zk[i] = Rk[i];
		}

		int Iteration = 0;

		do {
			Iteration++;
			/* Вычисляем числитель и знаменатель для коэффициента
			 * alpha = (rk-1,rk-1)/(Azk-1,zk-1) */
			Spz = 0;
			Spr = 0;
			for (i = 0; i < size; i++) {
				for (Sz[i] = 0, j = 0; j < size; j++) {
					Sz[i] += Acoef[i][j] * Zk[j];
				}
				Spz += Sz[i] * Zk[i];
				Spr += Rk[i] * Rk[i];
			}
			alpha = Spr / Spz;             /*  alpha    */


			/* Вычисляем вектор решения: xk = xk-1+ alpha * zk-1,
			вектор невязки: rk = rk-1 - alpha * A * zk-1 и числитель для beta равный (rk,rk) */
			Spr1 = 0;
			for (i = 0; i < size; i++) {
				Result[i] += alpha * Zk[i];
				Rk[i] -= alpha * Sz[i];
				Spr1 += Rk[i] * Rk[i];
			}
			kl++;

			/* Вычисляем  beta  */
			beta = Spr1 / Spr;

			/* Вычисляем вектор спуска: zk = rk+ beta * zk-1 */
			for (i = 0; i < size; i++)
				Zk[i] = Rk[i] + beta * Zk[i];
		}
		/* Проверяем условие выхода из итерационного цикла  */
		while (Spr1 / mf > diff* diff&& Iteration < max_iter);
		delete[] Zk;
		delete[] Rk;
		delete[] Sz;
	}
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
			conj_grad(accuracy);
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
	void coutAB()
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				cout << Acoef[i][j] << '\t';
			}
			cout << "=   " << Bcoef[i]<<endl;
		}
	}
	void coutRes()
	{
		for (int i = 0; i < size; i++)
		{
			cout << Result[i] << endl;
		}
	}
	bool readLESff(string path)
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
			double** A = new double* [_size];
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
			this->setEquasion(A, B, _size);
			input.close();
			for (int i = 0; i < _size; i++)
			{
				delete[] A[i];
			}
			delete[] A;
			delete[] B;
			return true;
		}
	}
};

bool readLESff(LinEquasSys &cur, string path);

int main()
{
	int size;
	double diff;
	LinEquasSys test;
	readLESff(test, _PATH);
	test.coutAB();
	test.doMath(3);
	test.coutRes();
}
