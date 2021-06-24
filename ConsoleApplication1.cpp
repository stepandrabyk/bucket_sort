// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <vector>

using namespace std;

void Initialize(int& size, int& number_of_blocks, int*& input, int*& length_of_blocks) {
	cout << "Enter number of input elements\n";
	cin >> size;
	cout << "Enter number of blocks\n";
	cin >> number_of_blocks;
	input = new int[size];
	length_of_blocks = new int[number_of_blocks];
}

void DataInitialization(int size, int* input) {
	for (int i = 0; i < size; i++) {
		input[i] = rand() % size;
	}
}

void SplitByBlocks(int size, int number_of_blocks, int* input, int* length_of_blocks) {
	int max_val = input[0];
	int min_val = input[0];
	for (int i = 1; i < size; i++) {
		if (max_val < input[i])
			max_val = input[i];
		if (min_val > input[i])
			min_val = input[i];
	}
	double range_of_block = double(max_val - min_val) / double(number_of_blocks);
	vector<int> block;
	vector<int> blocks;
	double lower_bound;
	double upper_bound;
	for (int i = 0; i < number_of_blocks; i++) {
		block.clear();
		lower_bound = min_val + i * range_of_block;
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

void BlockSort(int* input, int length, int start_position) {
	int tmp = 0;

	for (int i = 0; i < length; i++) {
		for (int j = length - 1; j >= i; j--) {
			if (input[start_position + j] < input[start_position + j - 1]) {
				tmp = input[start_position + j];
				input[start_position + j] = input[start_position + j - 1];
				input[start_position + j - 1] = tmp;
			}
		}
	}
}

void Sort(int size, int number_of_blocks, int* input, int* length_of_blocks) {
	int start_position = 0;
	for (int i = 0; i < number_of_blocks; i++) {
		BlockSort(input, length_of_blocks[i], start_position);
		start_position += length_of_blocks[i];
	}
}

void testPrintInput(int size, int* input) {
	cout << endl;
	for (int i = 0; i < size; i++) {
		cout << input[i] << " ";
	}
}

int main() {
	int size;
	int number_of_blocks;
	int* input;
	int* length_of_blocks;

	Initialize(size, number_of_blocks, input, length_of_blocks);

	double start = clock();

	DataInitialization(size, input);

	//testPrintInput(size, input);

	SplitByBlocks(size, number_of_blocks, input, length_of_blocks);

	Sort(size, number_of_blocks, input, length_of_blocks);

	double end = clock();

	//testPrintInput(size, input);

	cout << "\nTime of computing >> " << (end - start) / CLOCKS_PER_SEC;

	return 0;
}