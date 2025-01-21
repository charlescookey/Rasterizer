#pragma once
#include <chrono>
#include "matrix.h"
#include "zbuffer.h"

void redundantIdentityCall() {
	matrix temp;
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 2000; i++) {
		temp = matrix::makeTranslation_old(0.0f, 0.0f, 4.0f);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "identity call: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 2000; i++) {
		temp = matrix::makeTranslation(0.0f, 0.0f, 4.0f);
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << "no identity call: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";
}

void memsetOverLoops() {
	matrix temp;
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 2000; i++) {
		temp = matrix::makeIdentity_old();
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "loop: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 2000; i++) {
		temp = matrix::makeIdentity();
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << "Using Memset: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";
}

