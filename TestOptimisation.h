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

void clearSIMD() {
	Zbuffer<float> zbuffer;
	Zbuffer<float> zbuffer2;
	zbuffer.create(1024, 768);
	zbuffer2.create(1024, 768);
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 2000; i++) {
		zbuffer.clear();
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "SIMD clear: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 2000; i++) {
		zbuffer2.clear_old();
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << "Loop Clear: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";
}

void reserveVector() {

	std::vector<Mesh*> scene;
			
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 200; i++) {
		Mesh* m = new Mesh();
		*m = Mesh::makeRectangle(1.f, 1.f, 5.f,5.f);
		scene.push_back(m);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Using reserve: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";
	scene.clear();

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 200; i++) {
		Mesh* m = new Mesh();
		*m = Mesh::makeRectangle_old(1.f, 1.f, 5.f, 5.f);
		scene.push_back(m);
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << "No reserve: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";
}

void unrollRotateXYZ() {
	RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
	struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
	std::vector<rRot> rotations;

	matrix world = matrix::makeTranslation(4, 0.f, -6.f);
	matrix temp;

	// Create a grid of cubes with random rotations
	for (unsigned int y = 0; y < 200; y++) {
		rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
		rotations.emplace_back(r);
	}

	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 200; i++) {
		temp = world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Using unrolled: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";


	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 200; i++) {
		temp = world * matrix::makeRotateXYZ_old(rotations[i].x, rotations[i].y, rotations[i].z);
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << "multiplying x,y,z: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";
}

void render_old(Renderer&, Mesh*, matrix&, Light&);
void render(Renderer&, Mesh*, matrix&, Light&);

void precomputeTranformations() {
	Renderer renderer;
	matrix camera;
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

	//Mesh Cube = Mesh::makeCube(1.f);
	Mesh Cube = Mesh::makeSphere(1.0f, 10, 20);
	Cube.world = matrix::makeTranslation(-7.0f + (static_cast<float>(1) * 2.f), 5.0f - (static_cast<float>(1) * 2.f), -8.f);

	renderer.clear();


	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 300; i++) {
		render(renderer, &Cube, camera, L);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Precomputed Transformation: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";

	renderer.clear();

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 300; i++) {
		render_old(renderer, &Cube, camera, L);
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << "Computing transformation on the go: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";

}

void unrollMultiplication() {
	matrix test = matrix::makeRotateXYZ(0.0139100105, 0.0927862525, 0.09);
	matrix test2 = matrix::makeRotateXYZ(0.0139100105, 0.0927862525, 0.0777);


	matrix mull1;
	matrix mull2;;

	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 200; i++) {
		mull1 = test * test2;
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Using unrolled: "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";


	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 200; i++) {
		mull2 = test % test2;
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << "using loops "
		<< std::chrono::duration<double, std::milli>(end - start).count()
		<< "ms...\n";
}

}