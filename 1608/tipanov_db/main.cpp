#include "mpi.h" 
#include <Windows.h>
#include <time.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	int n = 666666;
	int sum = 0;
	int paralSum = 0;
	int tmpSum = 0;
	int ProcNum, ProcRank;
	double parTimeStart, parTimeEnd, timeStart, timeEnd;

	// ���� ���� MPI
	MPI_Init(&argc, &argv); //�������������
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum); //���-�� ���������
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);//���� �������� 
	MPI_Status Status; //���������� ��� ����� ��� ������ � ���������
	int countelem = 0;
	srand(time(0));
	rand();
	int *mas = new int[n];
	for (int i = 0; i<n; i++)
	{
		mas[i] = rand() % 10; //��������� ������ ��������� �������
	}
	countelem = n / ProcNum;
	int *tmp = new int[countelem];
	//��� 0-�� ��������
	if (ProcRank == 0)
	{
		parTimeStart = MPI_Wtime(); // ��������� �������� ������� ���������� ���������
		//���������� ������� ���������
		for (int i = 1; i<ProcNum; i++)
		{
			tmp = mas + countelem*(i - 1);
			MPI_Send(tmp, countelem, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		for (int i = (ProcNum - 1)*countelem; i<n; i++)
		{
			paralSum += mas[i];
		}
		//��������� ����������
		for (int i = 1; i<ProcNum; i++)
		{
			MPI_Recv(&tmpSum, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &Status);
			paralSum += tmpSum;
		}
		parTimeEnd = MPI_Wtime(); // ����� ��������

		//���������������� ������
		timeStart = MPI_Wtime();
		for (int i = 0; i<n; i++)
			sum += mas[i];
		timeEnd = MPI_Wtime();

		//����� �����������
		printf("sum = %d\ntime = %.9f\n\n", sum, timeEnd - timeStart);
		printf("Parallel sum = %d\ntime = %.9f\n", paralSum, parTimeEnd - parTimeStart);

	}

	//��� ���� ��������� ���������
	if (ProcRank != 0)
	{
		MPI_Recv(tmp, countelem, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
		for (int i = 0; i<countelem; i++)
			tmpSum += tmp[i];
		MPI_Send(&tmpSum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize(); //��������� ����� mpi
	delete mas;
	return 0;
}