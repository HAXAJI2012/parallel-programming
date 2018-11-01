#define _CRT_NO__WARNINGS
#include "mpi.h" 
#include <Windows.h>
#include <time.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	int sum = 0;
	int paralSum = 0;
	int tmpSum = 0;
	int totalSum = 0;
	int ProcNum, ProcRank;
	double parTimeStart, parTimeEnd, timeStart, timeEnd;

	// ���� ���� MPI
	MPI_Init(&argc, &argv); //�������������
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum); //���-�� ���������
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);//���� �������� 
	MPI_Status Status; //���������� ��� ����� ��� ������ � ���������
	int n = atoi(argv[1]);
	//int n = 666666;
	int countelem = 0;
	countelem = n / ProcNum;
	//��� 0-�� ��������
	if (ProcRank == 0)
	{
		srand(time(0));
		int *mas = new int[n];
		for (int i = 0; i < n; i++)
		{
			mas[i] = rand() % 10; //��������� ������ ��������� �������
		}

		//���������������� ������
		timeStart = MPI_Wtime();
		for (int i = 0; i < n; i++)
			sum += mas[i];
		timeEnd = MPI_Wtime();
		printf("sum = %d\ntime = %.9f\n\n", sum, timeEnd - timeStart);

		parTimeStart = MPI_Wtime(); // ��������� �������� ������� ���������� ���������
		//���������� ������� ���������
		for (int i = 1; i < ProcNum; i++)
		{
			int *tmp = mas + countelem * (i - 1);
			MPI_Send(tmp, countelem, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		for (int i = (ProcNum - 1)*countelem; i < n; i++)
		{
			paralSum += mas[i];
		}
		delete[] mas;
	}

	//��� ���� ��������� ���������
	if (ProcRank != 0)
	{
		int *tmp = new int[countelem];
		MPI_Recv(tmp, countelem, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
		for (int i = 0; i<countelem; i++)
			paralSum += tmp[i];
		delete[] tmp;
	}
	MPI_Reduce(&paralSum, &totalSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		parTimeEnd = MPI_Wtime();
		printf("Parallel sum = %d\ntime = %.9f\n", totalSum, parTimeEnd - parTimeStart);
	}
	MPI_Finalize(); //��������� ����� mpi
	return 0;
}