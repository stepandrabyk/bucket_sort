// MPI Sort.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <vector>

using namespace std;

int ProcNum;
int ProcRank;

void ProcessInitialization(int*& input, int& size, int& n_of_blocks) {
	if (ProcRank == 0) {
		cout << "Enter number of input elements\n";
		cin >> size;
		cout << "Enter number of blocks\n";
		cin >> n_of_blocks;

		/*while (true)
		{
			cout << "Enter number of blocks\n";
			cin >> n_of_blocks;
			if (n_of_blocks % ProcNum == 0)
				break;
			else
			{
				printf("Number of blocks must be divisible by %d \n", ProcNum);
				continue;
			}
		}*/
	}

	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_of_blocks, 1, MPI_INT, 0, MPI_COMM_WORLD);
	input = new int[size];
	//MPI_Bcast(&input, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void DataInitialization(int* input, int size, int*& blocks_per_proc, int num_of_blocks) {

	if (ProcRank == 0) {
		// input initialization
		for (int i = 0; i < size; i++) {
			input[i] = rand() % size;
		}

		// numbers of blocks per process
		for (int i = 0; i < ProcNum; i++) {
			if (i == ProcNum - 1)
				blocks_per_proc[i] = int(num_of_blocks / ProcNum) + num_of_blocks % ProcNum;
			else
				blocks_per_proc[i] = int(num_of_blocks / ProcNum);
		}
	}
	MPI_Bcast(input, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(blocks_per_proc, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);

}

void ParallelDataInitialization(int* input, int size, int*& blocks_per_proc, int num_of_blocks) {
	// numbers of blocks per process
	for (int i = 0; i < ProcNum; i++) {
		if (i == ProcNum - 1)
			blocks_per_proc[i] = int(num_of_blocks / ProcNum) + num_of_blocks % ProcNum;
		else
			blocks_per_proc[i] = int(num_of_blocks / ProcNum);
	}

	int to_generate = 0;
	if (ProcRank == ProcNum - 1)
		to_generate = int(size / ProcNum) + size % ProcNum;
	else
		to_generate = int(size / ProcNum);
	int* generated = new int[to_generate];
	for (int i = 0; i < to_generate; i++) {
		generated[i] = rand() % size;
	}

	MPI_Gather(generated, to_generate, MPI_INT, input, to_generate, MPI_INT, 0, MPI_COMM_WORLD);
}

void DataDistribution(int* input, int size, int* blocks_per_proc, int num_of_blocks, int* length_of_blocks, int*& p_blocks, int*& p_blocks_length) {
	if (ProcRank == 0) {
		// first step - find maximum and minimum values
		int max = input[0];
		int	min = input[0];
		for (int i = 0; i < size; i++) {
			if (max < input[i])
				max = input[i];
			if (min > input[i])
				min = input[i];
		}

		// second step - creating blocks and inserting data in them
		double range_of_block = (double(max) - double(min)) / double(num_of_blocks);
		vector<int> block;
		vector<int> blocks;
		double lower_bound;
		double upper_bound;
		for (int i = 0; i < num_of_blocks; i++) {
			block.clear();
			lower_bound = min + i * range_of_block;
			upper_bound = lower_bound + range_of_block;
			for (int j = 0; j < size; j++) {
				if (input[j] == -1)
					continue;

				if (lower_bound <= double(input[j]) && double(input[j]) <= upper_bound) {
					block.push_back(input[j]);
					blocks.push_back(input[j]);
					input[j] = -1;
				}
			}
			length_of_blocks[i] = block.size();
		}
		for (int i = 0; i < size; i++) {
			input[i] = blocks[i];
		}
		blocks.clear();
		block.clear();
	}

	MPI_Bcast(input, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(length_of_blocks, num_of_blocks, MPI_INT, 0, MPI_COMM_WORLD);

	// send length of blocks
	int* displ = new int[ProcNum];
	displ[0] = 0;
	for (int i = 1; i < ProcNum; i++) {
		displ[i] = displ[i - 1] + blocks_per_proc[i - 1];
	}
	p_blocks_length = new int[blocks_per_proc[ProcRank]];
	MPI_Scatterv(length_of_blocks, blocks_per_proc, displ, MPI_INT, p_blocks_length, blocks_per_proc[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);

	// send blocks
	int* send_num = new int[ProcNum];
	for (int i = 0; i < ProcNum; i++) {
		send_num[i] = 0;
		for (int j = 0; j < blocks_per_proc[i]; j++) {
			send_num[i] += length_of_blocks[j + displ[i]];
		}
	}

	displ[0] = 0;
	for (int i = 1; i < ProcNum; i++) {
		displ[i] = displ[i - 1] + send_num[i - 1];
	}

	p_blocks = new int[send_num[ProcRank]];
	MPI_Scatterv(input, send_num, displ, MPI_INT, p_blocks, send_num[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);

}

void SortBlock(int* blocks, int length, int start_position) {
	int tmp = 0;

	for (int i = 0; i < length; i++) {
		for (int j = length - 1; j >= i + 1; j--) {
			if (blocks[start_position + j] < blocks[start_position + j - 1]) {
				tmp = blocks[start_position + j];
				blocks[start_position + j] = blocks[start_position + j - 1];
				blocks[start_position + j - 1] = tmp;
			}
		}
	}
}

void Sort(/*int*& p_res, */int* p_blocks, int* p_blocks_length, int n_of_blocks) {
	int start_position = 0;
	/*
	vector<int>sorted_block;
	//vector<int>p_result;
	//
	int all_blocks_length = 0;
	for(int i = 0; i < n_of_blocks; i++) {
		all_blocks_length += p_blocks_length[i];
	}
	p_res = new int[all_blocks_length];
	//

	for (int i = 0; i < n_of_blocks; i++) {
		sorted_block = SortBlock(p_blocks, p_blocks_length[i], start_position);
		for (int j = 0; j < p_blocks_length[i]; i++) {
			p_res[start_position + j] = sorted_block[j];
		}
		start_position += p_blocks_length[i];
	}*/

	for (int i = 0; i < n_of_blocks; i++) {
		SortBlock(p_blocks, p_blocks_length[i], start_position);
		start_position += p_blocks_length[i];
	}
}

void ResultReplication(int* sorted_input, int* p_blocks, int* p_blocks_length, int* blocks_per_proc, int* blocks_length) {
	int send_count = 0;
	for (int i = 0; i < blocks_per_proc[ProcRank]; i++)
		send_count += p_blocks_length[i];

	int* recv_count = new int[ProcNum];
	int start_pos = 0;
	for (int i = 0; i < ProcNum; i++) {
		recv_count[i] = 0;
		for (int j = 0; j < blocks_per_proc[i]; j++) {
			recv_count[i] += blocks_length[start_pos + j];
		}
		start_pos += blocks_per_proc[i];
	}

	int* displ = new int[ProcNum];
	displ[0] = 0;
	for (int i = 1; i < ProcNum; i++) {
		displ[i] = displ[i - 1] + recv_count[i - 1];
	}


	/*
	if (ProcRank == 0) {
		cout << endl << "Send count " << send_count << endl;
		cout << "recv count\n";
		for (int i = 0; i < ProcNum; i++) {
			cout << recv_count[i] << " ";
		}

		cout << "\ndispl\n";
		for (int i = 0; i < ProcNum; i++) {
			cout << displ[i] << " ";
		}
	}
	*/

	MPI_Allgatherv(p_blocks, send_count, MPI_INT, sorted_input, recv_count, displ, MPI_INT, MPI_COMM_WORLD);

}

// test functions
void TestPrint_input(int* input, int size) {
	cout << "\nProc rank >>" << ProcRank << endl;

	for (int i = 0; i < size; i++) {
		cout << input[i] << " ";
	}
}

void TestPrint_b_p_p(int* blocks_per_proc) {
	cout << "Proc rank >>" << ProcRank << endl;

	for (int i = 0; i < ProcNum; i++) {
		cout << blocks_per_proc[i] << " ";
	}
}

void TestPrint_p_blocks(int* p_blocks, int* p_blocks_length, int* blocks_per_proc) {
	int n = 0;
	for (int i = 0; i < blocks_per_proc[ProcRank]; i++) {
		n += p_blocks_length[i];
	}
	cout << "\nProc rank " << ProcRank << endl;
	for (int i = 0; i < n; i++) {
		cout << p_blocks[i] << " ";
	}
}
//


int main(int argc, char* argv[]) {
 
	// place for variables
	int* input;
	int size;
	int num_of_blocks;
	int* length_of_blocks;
	int* blocks_per_proc = new int[ProcNum];

	int* p_blocks;
	int* p_blocks_length;
	// time variables
	double Start, Finish, Duration;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	//

	ProcessInitialization(input, size, num_of_blocks);

	Start = MPI_Wtime();

	length_of_blocks = new int[num_of_blocks];

	//DataInitialization(input, size, blocks_per_proc, num_of_blocks);
	ParallelDataInitialization(input, size, blocks_per_proc, num_of_blocks);

	/*
	if (ProcRank == 0) {
		TestPrint_input(input, size);
	}
	*/
	//TestPrint_b_p_p(blocks_per_proc);

	DataDistribution(input, size, blocks_per_proc, num_of_blocks, length_of_blocks, p_blocks, p_blocks_length);

	//TestPrint_input(input, size);
	//TestPrint_p_blocks(p_blocks, p_blocks_length, blocks_per_proc);

	Sort(p_blocks, p_blocks_length, blocks_per_proc[ProcRank]);

	//TestPrint_p_blocks(p_blocks, p_blocks_length, blocks_per_proc);

	int* res = new int[size];
	ResultReplication(res, p_blocks, p_blocks_length, blocks_per_proc, length_of_blocks);

	//TestPrint_input(res, size);

	/*
	if (ProcRank == 0) {
		TestPrint_input(res, size);
	}
	*/
	//

	Finish = MPI_Wtime();
	Duration = Finish - Start;
	if (ProcRank == 0)
		cout << "\nTime of computing >> " << Duration << endl;
	MPI_Finalize();
}