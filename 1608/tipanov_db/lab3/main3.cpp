#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <fstream>

using namespace std;

int* CreateVectorMatrix(int N) {
	int *matrix = new int[N*N];
	return matrix;
}

void PrintMatrix(int* matrix, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			cout << matrix[i*N + j] << " ";
		cout << endl;
	}
	cout << endl;
}

void fPrintMatrix(int* matrix, int N, char name[20]) {
	ofstream fout(name, ios_base::trunc);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			fout << matrix[i*N + j] << " ";
		fout << endl;
	}
	fout.close();
}

void RandMatrix(int* matrix1, int* matrix2, int N) {
	srand(time(0));
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			matrix1[i*N + j] = rand() % 10;
			matrix2[i*N + j] = rand() % 10;
		}
}

int* Trivial_alghorithm(int* matrix1, int* matrix2, int N) {
	int* Rez = CreateVectorMatrix(N);
	int sum;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			sum = 0;
			for (int k = 0; k < N; k++)
				sum += matrix1[i*N + k] * matrix2[k*N + j];
			Rez[i*N + j] = sum;
		}
	return Rez;
}

int* Add(int* matrix1, int* matrix2, int N) {
	int* Rez = CreateVectorMatrix(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Rez[i*N + j] = matrix1[i*N + j] + matrix2[i*N + j];

	return Rez;
}

int* Add(int* matrix1, int* matrix2, int* matrix3, int* matrix4, int N) {
	int* Rez = CreateVectorMatrix(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Rez[i*N + j] = matrix1[i*N + j] + matrix2[i*N + j] + matrix3[i*N + j] + matrix4[i*N + j];

	return Rez;
}

int* Sub(int* matrix1, int* matrix2, int N) {
	int* Rez = CreateVectorMatrix(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Rez[i*N + j] = matrix1[i*N + j] - matrix2[i*N + j];

	return Rez;
}

int* Sub(int* matrix1, int* matrix2, int* matrix3, int* matrix4, int N) {//сумма 3-х матриц вычесть четвертую
	int* Rez = CreateVectorMatrix(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Rez[i*N + j] = matrix1[i*N + j] + matrix2[i*N + j] + matrix3[i*N + j] - matrix4[i*N + j];

	return Rez;
}

int* Shtr_alg(int* matrix1, int* matrix2, int N, int threshold) {
	int* Rez;

	if (N <= threshold) {
		Rez = Trivial_alghorithm(matrix1, matrix2, N);
	}
	else {
		Rez = CreateVectorMatrix(N);
		N = N / 2;
		int* A[4];
		int* B[4];
		int* C[4];
		int* P[7];

		//Выделяем память под вспомогательные матрицы
		for (int i = 0; i < 4; i++) {
			A[i] = CreateVectorMatrix(N);
			B[i] = CreateVectorMatrix(N);
		}

		//Разбиваем матрицы на 4 блока
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++) {
				int index_new = i*N + j, index_old = 2 * i*N + j, N_N = 2 * N*N;
				A[0][index_new] = matrix1[index_old];
				A[1][index_new] = matrix1[index_old + N];
				A[2][index_new] = matrix1[index_old + N_N];
				A[3][index_new] = matrix1[index_old + N_N + N];

				B[0][index_new] = matrix2[index_old];
				B[1][index_new] = matrix2[index_old + N];
				B[2][index_new] = matrix2[index_old + N_N];
				B[3][index_new] = matrix2[index_old + N_N + N];
			}

		//Выполняем умножения 7 шт(рекурсивно)
		int *TMP = Add(A[0], A[3], N);
		int *_TMP = Add(B[0], B[3], N);
		P[0] = Shtr_alg(TMP, _TMP, N, threshold);
		delete[] TMP;
		delete[] _TMP;

		TMP = Add(A[2], A[3], N);
		P[1] = Shtr_alg(TMP, B[0], N, threshold);
		delete[] TMP;

		TMP = Sub(B[1], B[3], N);
		P[2] = Shtr_alg(A[0], TMP, N, threshold);
		delete[] TMP;

		TMP = Sub(B[2], B[0], N);
		P[3] = Shtr_alg(A[3], TMP, N, threshold);
		delete[] TMP;

		TMP = Add(A[0], A[1], N);
		P[4] = Shtr_alg(TMP, B[3], N, threshold);
		delete[] TMP;

		TMP = Sub(A[2], A[0], N);
		_TMP = Add(B[0], B[1], N);
		P[5] = Shtr_alg(TMP, _TMP, N, threshold);
		delete[] TMP;
		delete[] _TMP;

		TMP = Sub(A[1], A[3], N);
		_TMP = Add(B[2], B[3], N);
		P[6] = Shtr_alg(TMP, _TMP, N, threshold);
		delete[] TMP;
		delete[] _TMP;

		//Находим результирующие значения(блоки)
		C[0] = Sub(P[0], P[3], P[6], P[4], N);
		C[1] = Add(P[2], P[4], N);
		C[2] = Add(P[1], P[3], N);
		C[3] = Sub(P[0], P[2], P[5], P[1], N);

		//Формируем результирующую матрицу
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++) {
				Rez[i * 2 * N + j] = C[0][i*N + j];
				Rez[i * 2 * N + j + N] = C[1][i*N + j];
				Rez[i * 2 * N + j + 2 * N*N] = C[2][i*N + j];
				Rez[i * 2 * N + j + 2 * N*N + N] = C[3][i*N + j];
			}

		//Освобождаем выделенную память
		for (int i = 0; i < 4; i++) {
			delete[] A[i];
			delete[] B[i];
			delete[] C[i];
		}
		for (int i = 0; i < 7; i++)
			delete[] P[i];

	}

	return Rez;
}

int main(int argc, char** argv) {
	int *matr_A = nullptr, *matr_B = nullptr, *matr_Rez_Seq = nullptr, *matr_Rez_Par = nullptr, //исходные матрицы 
		**A = nullptr, **B = nullptr, **TMP_Rez = nullptr;
	int  N, thr = 64, sqr, new_N;

	int ProcNum, ProcRank;
	MPI_Status Status;

	double EndSeqAlg = 0; // Время завершения последовательного алгоритма
	double StartSeqAlg = 0;	// Время старта последовательного алгоритма
	double EndParAlg = 0; //Время завершения параллельного алгоритма
	double StartParAlg = 0;	 //Время старта паралельного алгоритма
	double TimeSeqAlg = 0;// Время, затраченное на работу последовательной версии алгоритма
	double TimeParAlg = 0;// -//- параллельной версии алгоритма

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		//Считывание данных из командной строки
		int k = 0;
		if (argc > 1)
			k = atoi(argv[1]);
		N = (int)pow(2.0, k);
		if (argc > 2)
			thr = atoi(argv[2]);

		//Создание и заполнение матриц
		matr_A = CreateVectorMatrix(N);
		matr_B = CreateVectorMatrix(N);
		matr_Rez_Par = CreateVectorMatrix(N);
		RandMatrix(matr_A, matr_B, N);

		//Вывод сформированных данных
		/*		char fA[20] = "C:\MatrA.txt";
		char fB[20] = "D:\MatrB.txt";

		if (N < 20){
		cout << "Matr A:" << endl;
		PrintMatrix(matr_A, N);
		cout << "Matr B:" << endl;
		PrintMatrix(matr_B, N);
		}
		else{
		fPrintMatrix(matr_A, N, fA);
		fPrintMatrix(matr_B, N, fB);
		}
		*/
		//Параллельный алгоритм
		StartParAlg = MPI_Wtime();

		//ВЫделение памяти под вспомогательные матрицы и разбиение матриц на блоки
		sqr = (int)sqrt((double)ProcNum), new_N = N / sqr;
		A = new int*[ProcNum], B = new int*[ProcNum];


		for (int i = 0; i < ProcNum; i++) {
			A[i] = CreateVectorMatrix(new_N);
			B[i] = CreateVectorMatrix(new_N);
		}

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++) {
				A[sqr*(i / new_N) + j / new_N][(i%new_N)*new_N + (j%new_N)] = matr_A[i*N + j];
				B[sqr*(i / new_N) + j / new_N][(i%new_N)*new_N + (j%new_N)] = matr_B[i*N + j];
			}

		//Рассылка данных другим процессам
		MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&thr, 1, MPI_INT, 0, MPI_COMM_WORLD);

		for (int i = 1; i < ProcNum; i++) {
			int coef_A = sqr*(i / sqr), coef_B = i % sqr;
			for (int j = 0; j < sqr; j++) {
				MPI_Send(A[coef_A], new_N*new_N, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(B[coef_B], new_N*new_N, MPI_INT, i, 0, MPI_COMM_WORLD);
				coef_A++;
				coef_B += sqr;
			}
		}

		//Swap указателей для однообразных вычислений
		for (int i = 0; i < sqr; i++) {
			int* TMP = B[i];
			B[i] = B[i*sqr];
			B[i*sqr] = TMP;
		}
	}

	//Прием данных от процесса-root и формирование нужных данных
	if (ProcRank != 0) {
		MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&thr, 1, MPI_INT, 0, MPI_COMM_WORLD);

		sqr = (int)sqrt((double)ProcNum), new_N = N / sqr;

		A = new int*[sqr], B = new int*[sqr];
		for (int i = 0; i < sqr; i++) {
			A[i] = CreateVectorMatrix(new_N);
			B[i] = CreateVectorMatrix(new_N);
		}

		for (int i = 0; i < sqr; i++) {
			MPI_Recv(A[i], new_N*new_N, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
			MPI_Recv(B[i], new_N*new_N, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
		}
	}

	//вычисление каждым процессом своего куска матрицы
	TMP_Rez = new int*[sqr + 1];
	for (int i = 0; i < sqr; i++) {
		TMP_Rez[i + 1] = Shtr_alg(A[i], B[i], new_N, thr);
	}
	if (ProcNum == 4)
		TMP_Rez[0] = Add(TMP_Rez[1], TMP_Rez[2], new_N);
	if (ProcNum == 16)
		TMP_Rez[0] = Add(TMP_Rez[1], TMP_Rez[2], TMP_Rez[3], TMP_Rez[4], new_N);


	//Освобождение вспомогательной памяти
	if (ProcRank == 0) {
		for (int i = 0; i < ProcNum; i++) {
			delete[] A[i];
			delete[] B[i];
		}
	}
	else {
		for (int i = 0; i < sqr; i++) {
			delete[] A[i];
			delete[] B[i];
		}
	}
	for (int i = 1; i < sqr + 1; i++) {
		delete[] TMP_Rez[i];
	}
	delete[] A;
	delete[] B;


	//Отправка результата на 0 процесс
	if (ProcRank != 0) {
		MPI_Send(TMP_Rez[0], new_N*new_N, MPI_INT, 0, 0, MPI_COMM_WORLD);
		delete[] TMP_Rez[0];
	}


	if (ProcRank == 0) {
		int coef = (int)sqrt((double)ProcNum);
		//ззаписываем результат совей работы
		for (int i = 0; i < new_N; i++)
			for (int j = 0; j < new_N; j++)
				matr_Rez_Par[coef*i*new_N + j] = TMP_Rez[0][i*new_N + j];
		for (int k = 1; k < ProcNum; k++) {
			//принимаем и записываем результаты работы других процессов
			MPI_Recv(TMP_Rez[0], new_N*new_N, MPI_INT, k, 0, MPI_COMM_WORLD, &Status);
			for (int i = 0; i < new_N; i++)
				for (int j = 0; j < new_N; j++)
					matr_Rez_Par[(k / coef)*new_N*N + (k%coef)*new_N + coef*i*new_N + j] = TMP_Rez[0][i*new_N + j];
		}

		EndParAlg = MPI_Wtime();
		TimeParAlg = EndParAlg - StartParAlg;

		//Вывод результата параллельного алгоритма
		cout << "ParAlg time = " << TimeParAlg << endl;
		/*		char fParAlg[20] = "C:\ParAlg.txt";
		if (N < 20){
		cout << "Matr C:" << endl;
		PrintMatrix(matr_Rez_Par, N);
		}
		else{
		fPrintMatrix(matr_Rez_Par, N, fParAlg);
		}
		*/
		delete[] TMP_Rez[0];

		//Последовательный алгоритм
		StartSeqAlg = MPI_Wtime();
		matr_Rez_Seq = Shtr_alg(matr_A, matr_B, N, thr);
		EndSeqAlg = MPI_Wtime();
		TimeSeqAlg = EndSeqAlg - StartSeqAlg;

		//Вывод результата последовательного алгоритма
		cout << "SeqAlg time = " << TimeSeqAlg << endl;
		/*		char fSeqAlg[20] = "C:\SeqAlg.txt";
		if (N < 20){
		cout << "Matr C:" << endl;
		PrintMatrix(matr_Rez_Seq, N);
		}
		else{
		fPrintMatrix(matr_Rez_Seq, N, fSeqAlg);
		}
		*/
		//Сравнение последовательной и параллельной версий
		if (TimeParAlg <= TimeSeqAlg)
			cout << "Parallel version(" << TimeParAlg << ") faster, then sequence(" << TimeSeqAlg << ")" << endl << "Speed-up = " << TimeSeqAlg / TimeParAlg << endl << "Efficiency = " << (TimeSeqAlg / TimeParAlg) / ProcNum << endl;
		else
			cout << "Sequence version(" << TimeSeqAlg << ") faster, then parallel(" << TimeParAlg << ")" << endl;

		bool flag = false;
		for (int k = 0; k < N; k++)
			for (int l = 0; l < N; l++)
				if (matr_Rez_Seq[k*N + l] != matr_Rez_Par[k*N + l])
					flag = true;

		//Сравнение результатов работы алгоритмов
		if (flag)
			cout << "matr_Rez_Seq != matr_Rez_Par" << endl;
		if (!flag)
			cout << "matr_Rez_Seq == matr_Rez_Par" << endl;
		delete[] matr_Rez_Par;
		delete[] matr_Rez_Seq;
		delete[] matr_A;
		delete[] matr_B;
	}

	delete[] TMP_Rez;

	MPI_Finalize();

	return 0;
}