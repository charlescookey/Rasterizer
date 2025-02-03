#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>
#include <thread>

#include <cmath>
#include "matrix.h"
#include "colour.h"
#include "mesh.h"
#include "zbuffer.h"
#include "renderer.h"
#include "RNG.h"
#include "light.h"
#include "triangle.h"


#include "TestOptimisation.h"


#define numberThreads 10

// Main rendering function that processes a mesh, transforms its vertices, applies lighting, and draws triangles on the canvas.
// Input Variables:
// - renderer: The Renderer object used for drawing.
// - mesh: Pointer to the Mesh object containing vertices and triangles to render.
// - camera: Matrix representing the camera's transformation.
// - L: Light object representing the lighting parameters.
void render_old(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
    // Combine perspective, camera, and world transformations for the mesh
    matrix p = renderer.perspective * camera * mesh->world;

    // Iterate through all triangles in the mesh
    for (triIndices& ind : mesh->triangles) {
        Vertex t[3]; // Temporary array to store transformed triangle vertices

        // Transform each vertex of the triangle
        for (unsigned int i = 0; i < 3; i++) {
            t[i].p = p * mesh->vertices[ind.v[i]].p; // Apply transformations
            t[i].p.divideW(); // Perspective division to normalize coordinates

            // Transform normals into world space for accurate lighting
            // no need for perspective correction as no shearing or non-uniform scaling
            t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal; 
            t[i].normal.normalise();

            // Map normalized device coordinates to screen space
            t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
            t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
            t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1]; // Invert y-axis

            // Copy vertex colours
            t[i].rgb = mesh->vertices[ind.v[i]].rgb;
        }

        // Clip triangles with Z-values outside [-1, 1]
        if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

        // Create a triangle object and render it
        triangle tri(t[0], t[1], t[2]);
        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

std::vector<Vertex> transformations(400); 

void render(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
    // Compute full transformation matrix
    matrix p = renderer.perspective * camera * mesh->world;

    size_t numVertices = mesh->vertices.size();
    float width = static_cast<float>(renderer.canvas.getWidth());
    float height = static_cast<float>(renderer.canvas.getHeight());

    // Step 1: Transform each unique vertex once
    for (size_t i = 0; i < numVertices; i++) {
        transformations[i].p = p * mesh->vertices[i].p;
        transformations[i].p.divideW();

        transformations[i].normal = mesh->world * mesh->vertices[i].normal;
        transformations[i].normal.normalise();

        transformations[i].p[0] = (transformations[i].p[0] + 1.f) * 0.5f * width;
        transformations[i].p[1] = (transformations[i].p[1] + 1.f) * 0.5f * height;
        transformations[i].p[1] = height - transformations[i].p[1]; // Invert y-axis

        transformations[i].rgb = mesh->vertices[i].rgb;
    }

    // Step 2: Process triangles 
    for (triIndices& ind : mesh->triangles) {
        // Directly pass transformed vertices instead of copying them into t[3]
        if (fabs(transformations[ind.v[0]].p[2]) > 1.0f ||
            fabs(transformations[ind.v[1]].p[2]) > 1.0f ||
            fabs(transformations[ind.v[2]].p[2]) > 1.0f) continue;

        // Create the triangle directly with transformed vertices
        triangle tri(transformations[ind.v[0]], transformations[ind.v[1]], transformations[ind.v[2]]);

        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

#include "test.h"

void transformVertex(Vertex& transform, Vertex& other, matrix &p, matrix &world, float height, float width) {
    transform.p = p * other.p;
    transform.p.divideW();

    transform.normal = world * other.normal;
    transform.normal.normalise();

    transform.p[0] = (transform.p[0] + 1.f) * 0.5f * width;
    transform.p[1] = (transform.p[1] + 1.f) * 0.5f * height;
    transform.p[1] = height - transform.p[1]; // Invert y-axis

    transform.rgb = other.rgb;
}

void transformVerticesRange(std::vector<Vertex>& transformations, std::vector<Vertex>& vertices, matrix& p, matrix& world, float height, float width, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        transformVertex(transformations[i], vertices[i], p, world, height, width);
    }
}

void render_mt_pool(Renderer& renderer, Mesh* mesh,  matrix& camera,  Light& L, ThreadPool& pool) {
    // Compute full transformation matrix
    matrix p = renderer.perspective * camera * mesh->world;

    size_t numVertices = mesh->vertices.size();
    std::vector<Vertex> transformations(numVertices);
    float width = static_cast<float>(renderer.canvas.getWidth());
    float height = static_cast<float>(renderer.canvas.getHeight());

    int middle = numVertices / 2;

    //std::thread receivingThread = std::thread(transformVerticesRange, std::ref(transformations), &mesh->vertices, &p, &mesh->world, height, width, 0, numVertices);

    auto future = pool.enqueue([&transformations, &mesh, &p, height, width, numVertices] {
        transformVerticesRange(transformations, mesh->vertices, p, mesh->world, height, width, 0, numVertices);
        });
    future.get();
    //transformVerticesRange(transformations, mesh->vertices, p, mesh->world, height, width, 0, numVertices);
    //receivingThread.join();

    for (triIndices& ind : mesh->triangles) {
        // Directly pass transformed vertices instead of copying them into t[3]
        if (fabs(transformations[ind.v[0]].p[2]) > 1.0f ||
            fabs(transformations[ind.v[1]].p[2]) > 1.0f ||
            fabs(transformations[ind.v[2]].p[2]) > 1.0f) continue;

        // Create the triangle directly with transformed vertices
        triangle tri(transformations[ind.v[0]], transformations[ind.v[1]], transformations[ind.v[2]]);

        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

void render_mt_pool_chunk(Renderer& renderer, Mesh* mesh,  matrix& camera,  Light& L, ThreadPool& pool) {
    // Compute full transformation matrix
    matrix p = renderer.perspective * camera * mesh->world;

    size_t numVertices = mesh->vertices.size();
    std::vector<Vertex> transformations(numVertices);
    float width = static_cast<float>(renderer.canvas.getWidth());
    float height = static_cast<float>(renderer.canvas.getHeight());

    int middle = numVertices / 2;
    std::vector<std::future<void>> voidResults;
    //std::thread receivingThread = std::thread(transformVerticesRange, std::ref(transformations), &mesh->vertices, &p, &mesh->world, height, width, 0, numVertices);
    size_t chunkSize = numVertices / numberThreads;
    for (int i = 0; i < numberThreads; i++) {
        size_t start = i * chunkSize;
        size_t end = (i == numberThreads - 1) ? numVertices : start + chunkSize;

        voidResults.emplace_back(pool.enqueue([&transformations, &mesh, &p, height, width, numVertices] {
            transformVerticesRange(transformations, mesh->vertices, p, mesh->world, height, width, 0, numVertices);
            }));
    }

    for (auto& t : voidResults) t.get();

    for (triIndices& ind : mesh->triangles) {
        // Directly pass transformed vertices instead of copying them into t[3]
        if (fabs(transformations[ind.v[0]].p[2]) > 1.0f ||
            fabs(transformations[ind.v[1]].p[2]) > 1.0f ||
            fabs(transformations[ind.v[2]].p[2]) > 1.0f) continue;

        // Create the triangle directly with transformed vertices
        triangle tri(transformations[ind.v[0]], transformations[ind.v[1]], transformations[ind.v[2]]);

        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

void render_mt_chunk(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
    // Compute full transformation matrix
    matrix p = renderer.perspective * camera * mesh->world;

    size_t numVertices = mesh->vertices.size();
    std::vector<Vertex> transformations(numVertices);
    float width = static_cast<float>(renderer.canvas.getWidth());
    float height = static_cast<float>(renderer.canvas.getHeight());

    std::vector<std::thread> threads;
    size_t chunkSize = numVertices / numberThreads;
    for (int i = 0; i < numberThreads; i++) {
        size_t start = i * chunkSize;
        size_t end = (i == numberThreads - 1) ? numVertices : start + chunkSize;
        threads.emplace_back(transformVerticesRange, std::ref(transformations), std::ref(mesh->vertices), std::ref(p), std::ref(mesh->world), height, width, start, end);
    }
    for (auto& t : threads) t.join();

    for (triIndices& ind : mesh->triangles) {
        // Directly pass transformed vertices instead of copying them into t[3]
        if (fabs(transformations[ind.v[0]].p[2]) > 1.0f ||
            fabs(transformations[ind.v[1]].p[2]) > 1.0f ||
            fabs(transformations[ind.v[2]].p[2]) > 1.0f) continue;

        // Create the triangle directly with transformed vertices
        triangle tri(transformations[ind.v[0]], transformations[ind.v[1]], transformations[ind.v[2]]);

        tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

// Test scene function to demonstrate rendering with user-controlled transformations
// No input variables
void sceneTest() {
    Renderer renderer;
    // create light source {direction, diffuse intensity, ambient intensity}
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    L.omega_i.normalise();
    // camera is just a matrix
    matrix camera;// = matrix::makeIdentity(); // Initialize the camera with identity matrix

    bool running = true; // Main loop control variable

    std::vector<Mesh*> scene; // Vector to store scene objects

    // Create a sphere and a rectangle mesh
    Mesh mesh = Mesh::makeSphere(1.0f, 10, 20);
    //Mesh mesh2 = Mesh::makeRectangle(-2, -1, 2, 1);

    // add meshes to scene
    scene.push_back(&mesh);
   // scene.push_back(&mesh2); 

    float x = 0.0f, y = 0.0f, z = -4.0f; // Initial translation parameters
    mesh.world = matrix::makeTranslation(x, y, z);
    //mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput(); // Handle user input
        renderer.clear(); // Clear the canvas for the next frame

        // Apply transformations to the meshes
     //   mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);
        mesh.world = matrix::makeTranslation(x, y, z);

        // Handle user inputs for transformations
        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;
        if (renderer.canvas.keyPressed('A')) x += -0.1f;
        if (renderer.canvas.keyPressed('D')) x += 0.1f;
        if (renderer.canvas.keyPressed('W')) y += 0.1f;
        if (renderer.canvas.keyPressed('S')) y += -0.1f;
        if (renderer.canvas.keyPressed('Q')) z += 0.1f;
        if (renderer.canvas.keyPressed('E')) z += -0.1f;

        // Render each object in the scene
        for (auto& m : scene)
            render(renderer, m, camera, L);

        renderer.present(); // Display the rendered frame
    }
}

// Utility function to generate a random rotation matrix
// No input variables
matrix makeRandomRotation() {
    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
    unsigned int r = rng.getRandomInt(0, 3);

    switch (r) {
    case 0: return matrix::makeRotateX(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 1: return matrix::makeRotateY(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 2: return matrix::makeRotateZ(rng.getRandomFloat(0.f, 2.0f * M_PI));
    default: return matrix();//::makeIdentity();
    }
}

// Function to render a scene with multiple objects and dynamic transformations
// No input variables
void scene1() {
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    L.omega_i.normalise();

    bool running = true;

    std::vector<Mesh*> scene;
    scene.reserve(40);

    // Create a scene of 40 cubes with random rotations
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.emplace_back(m);
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.emplace_back(m);
    }

    float zoffset = 8.0f; // Initial camera Z-offset
    float step = -0.1f;  // Step size for camera movement

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

        // Rotate the first two cubes in the scene
        scene[0]->world *=  matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
        scene[1]->world *=  matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        zoffset += step;
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        for (auto& m : scene)
            render(renderer, m, camera, L);
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

void scene1_mt() {
    //ThreadPool pool(numberThreads);
    std::vector<std::thread> threads;
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    L.omega_i.normalise();

    bool running = true;

    std::vector<Mesh*> scene;
    scene.reserve(40);

    // Create a scene of 40 cubes with random rotations
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.emplace_back(m);
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.emplace_back(m);
    }

    float zoffset = 8.0f; // Initial camera Z-offset
    float step = -0.1f;  // Step size for camera movement

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    std::vector<std::future<void>> voidResults;
    size_t chunkSize = 0;
    size_t Rotsize = 0;

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

        // Rotate the first two cubes in the scene
        scene[0]->world *= matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
        scene[1]->world *= matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        zoffset += step;
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        std::atomic<size_t> RenderCounter(0);
        chunkSize = scene.size() / numberThreads;
        Rotsize = scene.size();
        for (int j = 0; j < numberThreads; j++) {
            threads.emplace_back([&scene, &renderer, Rotsize, &RenderCounter, &camera, &L] {
                        size_t i;
                        while ((i = RenderCounter.fetch_add(1, std::memory_order_relaxed)) < Rotsize) {
                            render(renderer, scene[i], camera, L);
                        };
                        });
                }
        for (auto& t : threads) t.join();
        threads.clear();

        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

void scene1_mt_pool() {
    ThreadPool pool(numberThreads);
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    L.omega_i.normalise();

    bool running = true;

    std::vector<Mesh*> scene;
    scene.reserve(40);

    // Create a scene of 40 cubes with random rotations
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.emplace_back(m);
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.emplace_back(m);
    }

    float zoffset = 8.0f; // Initial camera Z-offset
    float step = -0.1f;  // Step size for camera movement

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    std::vector<std::future<void>> voidResults;
    size_t chunkSize = 0;
    size_t Rotsize = 0;

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

        // Rotate the first two cubes in the scene
        scene[0]->world *= matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
        scene[1]->world *= matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        zoffset += step;
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        std::atomic<size_t> RenderCounter(0);
        chunkSize = scene.size() / numberThreads;
        Rotsize = scene.size();

        for (int i = 0; i < numberThreads; i++) {
            voidResults.emplace_back(pool.enqueue([&scene, &renderer, Rotsize, &RenderCounter, &camera, &L] {
                size_t i;
                while ((i = RenderCounter.fetch_add(1, std::memory_order_relaxed)) < Rotsize) {
                    render(renderer, scene[i], camera, L);
                };
                }));
        }

        for (auto& t : voidResults) t.get();
        voidResults.clear();

        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
    Renderer renderer;
    matrix camera;// = matrix::makeIdentity();//UNROLL??//ALSO USE RESEVER SO IT DOESNT CREATE AND CAPACITY INSCREWASE
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    L.omega_i.normalise();

    std::vector<Mesh*> scene;
    scene.reserve(49);

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;
    rotations.reserve(48);

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    // Create a grid of cubes with random rotations
    for (unsigned int y = 0; y < 6; y++) {
        for (unsigned int x = 0; x < 8; x++) {
            Mesh* m = new Mesh();
            *m = Mesh::makeCube(1.f);
            scene.emplace_back(m);
            m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
            rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
            rotations.emplace_back(r);
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.emplace_back(sphere);
    float sphereOffset = -6.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    bool running = true;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        // Rotate each cube in the grid
        for (unsigned int i = 0; i < rotations.size(); i++)
            scene[i]->world *=  matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        for (auto& m : scene)
            render(renderer, m, camera, L);
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

void scene2_mt() {
    ThreadPool pool(numberThreads);

    Renderer renderer;
    matrix camera;// = matrix::makeIdentity();//UNROLL??//ALSO USE RESEVER SO IT DOESNT CREATE AND CAPACITY INSCREWASE
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    L.omega_i.normalise();

    std::vector<Mesh*> scene;
    scene.reserve(49);

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;
    rotations.reserve(48);

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    // Create a grid of cubes with random rotations
    for (unsigned int y = 0; y < 6; y++) {
        for (unsigned int x = 0; x < 8; x++) {
            Mesh* m = new Mesh();
            *m = Mesh::makeCube(1.f);
            scene.emplace_back(m);
            m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
            rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
            rotations.emplace_back(r);
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.emplace_back(sphere);
    float sphereOffset = -6.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    bool running = true;
    std::vector<std::future<void>> voidResults;
    size_t chunkSize = 0;
    size_t Rotsize = 0;

    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        std::atomic<size_t> counter(0);
        chunkSize = rotations.size() / numberThreads;
        Rotsize = rotations.size();
        for (int j = 0; j < numberThreads; j++) {
            voidResults.emplace_back(pool.enqueue([&scene, &rotations, Rotsize, &counter] {
                size_t i;
                while ((i = counter.fetch_add(1, std::memory_order_relaxed)) < Rotsize) {
                    scene[i]->world *= matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);
                };
                }));
        }

        for (auto& t : voidResults) t.get();
        voidResults.clear();

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        std::atomic<size_t> RenderCounter(0);
        chunkSize = scene.size() / numberThreads;
        Rotsize = scene.size();

        for (int j = 0; j < numberThreads; j++) {
            voidResults.emplace_back(pool.enqueue([&scene, &renderer, Rotsize, &RenderCounter, &camera, &L] {
                size_t i;
                while ((i = RenderCounter.fetch_add(1, std::memory_order_relaxed)) < Rotsize) {
                    render(renderer, scene[i], camera, L);
                };
                }));
        }

        for (auto& t : voidResults) t.get();
        voidResults.clear();

        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

void customScene() {
    Renderer renderer;
    matrix camera =  matrix::makeTranslation(0.f, 0.f, -10.f);;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    L.omega_i.normalise();

    std::vector<Mesh*> scene;
    scene.reserve(10001);

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;
    rotations.reserve(1000);

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    for (unsigned int z = 0; z < 10; z++) {          // 10 layers (z-axis)
        for (unsigned int y = 0; y < 10; y++) {      // 10 rows (y-axis)
            for (unsigned int x = 0; x < 10; x++) {  // 10 columns (x-axis)
                Mesh* m = new Mesh();
                *m = Mesh::makeCube(0.5f);            
                scene.emplace_back(m);               

                // Position the cube in a 10x10x10 grid
                m->world = matrix::makeTranslation(
                    -9.0f + (static_cast<float>(x) * 2.f),
                    9.0f - (static_cast<float>(y) * 2.f),   
                    -10.0f + (static_cast<float>(z) * 2.f)  
                );
                rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
                rotations.emplace_back(r);
            }
        }
    }
    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.emplace_back(sphere);
    float sphereOffset = -15.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    bool running = true;
    std::vector<std::future<void>> voidResults;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

         //Rotate each cube in the grid
        for (unsigned int i = 0; i < rotations.size(); i++)
            scene[i]->world *= matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 15.0f || sphereOffset < -15.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        for (auto& m : scene)
            render(renderer, m, camera, L);

        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

void customScene_mt() {
    ThreadPool pool(numberThreads);

    Renderer renderer;
    matrix camera = matrix::makeTranslation(0.f, 0.f, -10.f);
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    L.omega_i.normalise();

    std::vector<Mesh*> scene;
    scene.reserve(10001);

    struct rRot { float x; float y; float z; }; 
    std::vector<rRot> rotations;
    rotations.reserve(1000);

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    for (unsigned int z = 0; z < 10; z++) {          // 10 layers (z-axis)
        for (unsigned int y = 0; y < 10; y++) {      // 10 rows (y-axis)
            for (unsigned int x = 0; x < 10; x++) {  // 10 columns (x-axis)
                Mesh* m = new Mesh();
                *m = Mesh::makeCube(0.5f);            // Create a cube
                scene.emplace_back(m);               // Add the cube to the scene

                // Position the cube in a 10x10x10 grid
                m->world = matrix::makeTranslation(
                    -9.0f + (static_cast<float>(x) * 2.f),
                    9.0f - (static_cast<float>(y) * 2.f),   
                    -10.0f + (static_cast<float>(z) * 2.f)  
                );

                // Add random rotations
                rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
                rotations.emplace_back(r);
            }
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.emplace_back(sphere);
    float sphereOffset = -15.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    bool running = true;
    std::vector<std::future<void>> voidResults;
    size_t chunkSize = 0;
    size_t Rotsize = 0;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

         //Rotate each cube in the grid
        std::atomic<size_t> counter(0);
        chunkSize = rotations.size() / numberThreads;
        Rotsize = rotations.size();
        for (int j = 0; j < numberThreads; j++) {
            voidResults.emplace_back(pool.enqueue([&scene, &rotations, Rotsize, &counter] {
                size_t i;
                while ((i = counter.fetch_add(1, std::memory_order_relaxed)) < Rotsize) {
                    scene[i]->world *= matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);
                };
                }));
        }

        for (auto& t : voidResults) t.get();
        voidResults.clear();

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 15.0f || sphereOffset < -15.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        for (auto& m : scene)
            render(renderer, m, camera, L);

        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

void customScene_mt2() {
    ThreadPool pool(numberThreads);

    Renderer renderer;
    matrix camera = matrix::makeTranslation(0.f, 0.f, -10.f);
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
    L.omega_i.normalise();

    std::vector<Mesh*> scene;
    scene.reserve(10001);

    struct rRot { float x; float y; float z; };
    std::vector<rRot> rotations;
    rotations.reserve(1000);

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    for (unsigned int z = 0; z < 10; z++) {          // 10 layers (z-axis)
        for (unsigned int y = 0; y < 10; y++) {      // 10 rows (y-axis)
            for (unsigned int x = 0; x < 10; x++) {  // 10 columns (x-axis)
                Mesh* m = new Mesh();
                *m = Mesh::makeCube(0.5f);            // Create a cube
                scene.emplace_back(m);               // Add the cube to the scene

                // Position the cube in a 10x10x10 grid
                m->world = matrix::makeTranslation(
                    -9.0f + (static_cast<float>(x) * 2.f),
                    9.0f - (static_cast<float>(y) * 2.f),
                    -10.0f + (static_cast<float>(z) * 2.f)
                );

                // Add random rotations
                rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
                rotations.emplace_back(r);
            }
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.emplace_back(sphere);
    float sphereOffset = -15.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;

    bool running = true;
    std::vector<std::future<void>> voidResults;
    size_t chunkSize = 0;
    size_t Rotsize = 0;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        //Rotate each cube in the grid
        std::atomic<size_t> counter(0);
        chunkSize = rotations.size() / numberThreads;
        Rotsize = rotations.size();
        for (int j = 0; j < numberThreads; j++) {
            voidResults.emplace_back(pool.enqueue([&scene, &rotations, Rotsize, &counter] {
                size_t i;
                while ((i = counter.fetch_add(1, std::memory_order_relaxed)) < Rotsize) {
                    scene[i]->world *= matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);
                };
                }));
        }

        for (auto& t : voidResults) t.get();
        voidResults.clear();

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 15.0f || sphereOffset < -15.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << cycle / 2 << " :" << std::chrono::duration<double, std::milli>(end - start).count() << "ms\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        std::atomic<size_t> RenderCounter(0);
        chunkSize = scene.size() / numberThreads;
        Rotsize = scene.size();

        for (int j = 0; j < numberThreads; j++) {
            voidResults.emplace_back(pool.enqueue([&scene, &renderer, Rotsize, &RenderCounter, &camera, &L] {
                size_t i;
                while ((i = RenderCounter.fetch_add(1, std::memory_order_relaxed)) < Rotsize) {
                    render(renderer, scene[i], camera, L);
                };
                }));
        }

        for (auto& t : voidResults) t.get();
        voidResults.clear();

        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Entry point of the application
// No input variables
int main() {
    // Uncomment the desired scene function to run
    //scene1();
    //scene1_mt();
    //scene1_mt_pool();
    //scene2();
    //scene2_mt();
    //sceneTest();
    //customScene();
    customScene_mt2();
    //customScene_mt();

    //redundantIdentityCall();
    //memsetOverLoops();
    //clearSIMD();
    //reserveVector();
    //unrollRotateXYZ();
    //precomputeTranformations();
    //unrollMultiplication();
    //usingMove();
    //avoidInterpolations();

    return 0;
}