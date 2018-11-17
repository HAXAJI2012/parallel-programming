#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define ROOT_NUM 0
#define type int

using namespace std;

void mysum(const void *arr1, void *arr2, int n, MPI_Datatype t);
int mylog2(int x);
void myallreduce(void* send, void* recv, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int main(int argc, char* argv[])
{
	int n;
	type *arr1 = nullptr, *arr2 = nullptr;
	int ProcRank, ProcNum;
	double timeStart, timeFinish;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	srand(time(0));
	if (argc == 2)
		n = atoi(argv[1]);
	else
		n = 10;
	arr1 = new type[n];
	arr2 = new type[n];
	for (int i = 0; i < n; i++)
	{
		arr1[i] = (type)rand();
		arr2[i] = 0;
		//arr1[i] = i; // used for tested!! [changed]
		//arr2[i] = 0;
	}

	//MPI_Allreduce(arr1, arr2, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // used for tests!! [deleted]

	timeStart = MPI_Wtime();
	myallreduce(arr1, arr2, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	timeFinish = MPI_Wtime() - timeStart;

	if (ProcRank == ROOT_NUM)
	{
		for (int i = 0; i < n; i++)
			cout << (type)arr1[i] << " ";
		cout << "\nTime: " << timeFinish << endl;
		// used one more FOR-cycle for tests (allreduce comparison) [deleted]
	}
	MPI_Finalize();
	return 0;
}

void myallreduce(void* send, void* recv, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
	int ProcNum, ProcRank;
	MPI_Comm_rank(comm, &ProcRank);
	MPI_Comm_size(comm, &ProcNum);
	int h = (int)(log((float)ProcNum) / log(2.0));
	if (mylog2(ProcNum) == 0)
		h++;
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < ProcNum; j += pow(2, (float)(i + 1)))
		{
			int k = j + pow(2, (float)i);
			if (k < ProcNum)
			{
				if (ProcRank == j)
				{
					MPI_Status status;
					MPI_Recv(recv, count, datatype, k, 0, MPI_COMM_WORLD, &status);
					mysum(recv, send, count, datatype);
				}
				if (ProcRank == k)
					MPI_Send(send, count, datatype, j, 0, MPI_COMM_WORLD);
			}
		}
	}
	
	for (int i = 1; i < ProcNum; i++)
	{
		for (int j = 0; j < i; j++)
		{
			int t = 1;
			for (int k = 1; k < i; k++)
				t *= 2;
			if (ProcRank == j)
				if ((t + j != 0) && (t + j < ProcNum))
					MPI_Send(send, count, datatype, t + j, 0, comm);
			if ((t + j != 0) && (t + j < ProcNum) && (t + j == ProcRank))
			{
				MPI_Status status;
				MPI_Recv(send, count, datatype, j, 0, comm, &status);
			}
		}
		MPI_Barrier(comm);
	}
}

int mylog2(int x)
{
	int k = x;
	int i = 0;
	while (k % 2 == 0)
	{
		i++;
		k = k / 2;
	}
	if ((int)(pow(2, (float)(i))) == x)
		return i;
	else
		return 0;
}

void mysum(const void *arr1, void *arr2, int n, MPI_Datatype t)
{
	if (t == MPI_INT)
	{
		for (int i = 0; i < n; i++)
			((int*)arr2)[i] += ((int*)arr1)[i];
	}
	if (t == MPI_FLOAT)
	{
		for (int i = 0; i < n; i++)
			((float*)arr2)[i] += ((float*)arr1)[i];
	}
	if (t == MPI_DOUBLE)
	{
		for (int i = 0; i < n; i++)
			((double*)arr2)[i] += ((double*)arr1)[i];
	}
}