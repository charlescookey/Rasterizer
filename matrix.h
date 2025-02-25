#pragma once

#include <iostream>
#include <vector>
#include "vec4.h"

// Matrix class for 4x4 transformation matrices
class matrix {
    union {
        float m[4][4]; // 2D array representation of the matrix
        float a[16];   // 1D array representation of the matrix for linear access
    };

public:
    // Default constructor initializes the matrix as an identity matrix
    matrix() {
        identity();
    }

    // Access matrix elements by row and column
    float& operator()(unsigned int row, unsigned int col) { return m[row][col]; }

    // Display the matrix elements in a readable format
    void display() {
        for (unsigned int i = 0; i < 4; i++) {
            for (unsigned int j = 0; j < 4; j++)
                std::cout << m[i][j] << '\t';
            std::cout << std::endl;
        }
    }

    // Multiply the matrix by a 4D vector
    // Input Variables:
    // - v: vec4 object to multiply with the matrix
    // Returns the resulting transformed vec4
    vec4 operator * (const vec4& v) const {
        vec4 result;
        result[0] = a[0] * v[0] + a[1] * v[1] + a[2] * v[2] + a[3] * v[3];
        result[1] = a[4] * v[0] + a[5] * v[1] + a[6] * v[2] + a[7] * v[3];
        result[2] = a[8] * v[0] + a[9] * v[1] + a[10] * v[2] + a[11] * v[3];
        result[3] = a[12] * v[0] + a[13] * v[1] + a[14] * v[2] + a[15] * v[3];
        return result;
    }

    // Multiply the matrix by another matrix
    // Input Variables:
    // - mx: Another matrix to multiply with
    // Returns the resulting matrix
    matrix operator % (const matrix& mx) const {
        matrix ret;
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                ret.a[row * 4 + col] =
                    a[row * 4 + 0] * mx.a[0 * 4 + col] +
                    a[row * 4 + 1] * mx.a[1 * 4 + col] +
                    a[row * 4 + 2] * mx.a[2 * 4 + col] +
                    a[row * 4 + 3] * mx.a[3 * 4 + col];
            }
        }
        return ret;
    }

    matrix operator *(const matrix& mx) const
    {
        matrix ret;
        ret.a[0] = mx.a[0] * a[0] + mx.a[4] * a[1] + mx.a[8] * a[2] + mx.a[12] * a[3];
        ret.a[1] = mx.a[1] * a[0] + mx.a[5] * a[1] + mx.a[9] * a[2] + mx.a[13] * a[3];
        ret.a[2] = mx.a[2] * a[0] + mx.a[6] * a[1] + mx.a[10] * a[2] + mx.a[14] * a[3];
        ret.a[3] = mx.a[3] * a[0] + mx.a[7] * a[1] + mx.a[11] * a[2] + mx.a[15] * a[3];
        ret.a[4] = mx.a[0] * a[4] + mx.a[4] * a[5] + mx.a[8] * a[6] + mx.a[12] * a[7];
        ret.a[5] = mx.a[1] * a[4] + mx.a[5] * a[5] + mx.a[9] * a[6] + mx.a[13] * a[7];
        ret.a[6] = mx.a[2] * a[4] + mx.a[6] * a[5] + mx.a[10] * a[6] + mx.a[14] * a[7];
        ret.a[7] = mx.a[3] * a[4] + mx.a[7] * a[5] + mx.a[11] * a[6] + mx.a[15] * a[7];
        ret.a[8] = mx.a[0] * a[8] + mx.a[4] * a[9] + mx.a[8] * a[10] + mx.a[12] * a[11];
        ret.a[9] = mx.a[1] * a[8] + mx.a[5] * a[9] + mx.a[9] * a[10] + mx.a[13] * a[11];
        ret.a[10] = mx.a[2] * a[8] + mx.a[6] * a[9] + mx.a[10] * a[10] + mx.a[14] * a[11];
        ret.a[11] = mx.a[3] * a[8] + mx.a[7] * a[9] + mx.a[11] * a[10] + mx.a[15] * a[11];
        ret.a[12] = mx.a[0] * a[12] + mx.a[4] * a[13] + mx.a[8] * a[14] + mx.a[12] * a[15];
        ret.a[13] = mx.a[1] * a[12] + mx.a[5] * a[13] + mx.a[9] * a[14] + mx.a[13] * a[15];
        ret.a[14] = mx.a[2] * a[12] + mx.a[6] * a[13] + mx.a[10] * a[14] + mx.a[14] * a[15];
        ret.a[15] = mx.a[3] * a[12] + mx.a[7] * a[13] + mx.a[11] * a[14] + mx.a[15] * a[15];
        return ret;
    }

    matrix& operator *=(const matrix& mx)
    {
        matrix ret;
        ret.a[0] = mx.a[0] * a[0] + mx.a[4] * a[1] + mx.a[8] * a[2] + mx.a[12] * a[3];
        ret.a[1] = mx.a[1] * a[0] + mx.a[5] * a[1] + mx.a[9] * a[2] + mx.a[13] * a[3];
        ret.a[2] = mx.a[2] * a[0] + mx.a[6] * a[1] + mx.a[10] * a[2] + mx.a[14] * a[3];
        ret.a[3] = mx.a[3] * a[0] + mx.a[7] * a[1] + mx.a[11] * a[2] + mx.a[15] * a[3];
        ret.a[4] = mx.a[0] * a[4] + mx.a[4] * a[5] + mx.a[8] * a[6] + mx.a[12] * a[7];
        ret.a[5] = mx.a[1] * a[4] + mx.a[5] * a[5] + mx.a[9] * a[6] + mx.a[13] * a[7];
        ret.a[6] = mx.a[2] * a[4] + mx.a[6] * a[5] + mx.a[10] * a[6] + mx.a[14] * a[7];
        ret.a[7] = mx.a[3] * a[4] + mx.a[7] * a[5] + mx.a[11] * a[6] + mx.a[15] * a[7];
        ret.a[8] = mx.a[0] * a[8] + mx.a[4] * a[9] + mx.a[8] * a[10] + mx.a[12] * a[11];
        ret.a[9] = mx.a[1] * a[8] + mx.a[5] * a[9] + mx.a[9] * a[10] + mx.a[13] * a[11];
        ret.a[10] = mx.a[2] * a[8] + mx.a[6] * a[9] + mx.a[10] * a[10] + mx.a[14] * a[11];
        ret.a[11] = mx.a[3] * a[8] + mx.a[7] * a[9] + mx.a[11] * a[10] + mx.a[15] * a[11];
        ret.a[12] = mx.a[0] * a[12] + mx.a[4] * a[13] + mx.a[8] * a[14] + mx.a[12] * a[15];
        ret.a[13] = mx.a[1] * a[12] + mx.a[5] * a[13] + mx.a[9] * a[14] + mx.a[13] * a[15];
        ret.a[14] = mx.a[2] * a[12] + mx.a[6] * a[13] + mx.a[10] * a[14] + mx.a[14] * a[15];
        ret.a[15] = mx.a[3] * a[12] + mx.a[7] * a[13] + mx.a[11] * a[14] + mx.a[15] * a[15];

        *this = std::move(ret);
        return *this;
    }

    // Create a perspective projection matrix
    // Input Variables:
    // - fov: Field of view in radians
    // - aspect: Aspect ratio of the viewport
    // - n: Near clipping plane
    // - f: Far clipping plane
    // Returns the perspective matrix
    static matrix makePerspective(float fov, float aspect, float n, float f) {
        matrix m;
        m.zero();
        float tanHalfFov = std::tan(fov / 2.0f);

        m.a[0] = 1.0f / (aspect * tanHalfFov);
        m.a[5] = 1.0f / tanHalfFov;
        m.a[10] = -f / (f - n);
        m.a[11] = -(f * n) / (f - n);
        m.a[14] = -1.0f;
        return m;
    }

    // Create a translation matrix
    // Input Variables:
    // - tx, ty, tz: Translation amounts along the X, Y, and Z axes
    // Returns the translation matrix
    static matrix makeTranslation(float tx, float ty, float tz) {
        matrix m;
        m.a[3] = tx;
        m.a[7] = ty;
        m.a[11] = tz;
        return m;
    }
    static matrix makeTranslation_old(float tx, float ty, float tz) {
        matrix m;
        m.identity();
        m.a[3] = tx;
        m.a[7] = ty;
        m.a[11] = tz;
        return m;
    }

    // Create a rotation matrix around the Z-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateZ(float aRad) {
        matrix m;
        m.a[0] = std::cos(aRad);
        m.a[1] = -std::sin(aRad);
        m.a[4] = std::sin(aRad);
        m.a[5] = std::cos(aRad);
        return m;
    }

    // Create a rotation matrix around the X-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateX(float aRad) {
        matrix m;
        m.a[5] = std::cos(aRad);
        m.a[6] = -std::sin(aRad);
        m.a[9] = std::sin(aRad);
        m.a[10] = std::cos(aRad);
        return m;
    }

    // Create a rotation matrix around the Y-axis
    // Input Variables:
    // - aRad: Rotation angle in radians
    // Returns the rotation matrix
    static matrix makeRotateY(float aRad) {
        matrix m;
        m.a[0] = std::cos(aRad);
        m.a[2] = std::sin(aRad);
        m.a[8] = -std::sin(aRad);
        m.a[10] = std::cos(aRad);
        return m;
    }

    // Create a composite rotation matrix from X, Y, and Z rotations
    // Input Variables:
    // - x, y, z: Rotation angles in radians around each axis
    // Returns the composite rotation matrix
    static matrix makeRotateXYZ_old(float x, float y, float z) {
        return matrix::makeRotateX(x) * matrix::makeRotateY(y) * matrix::makeRotateZ(z);
    }

    static matrix makeRotateXYZ(float x, float y, float z) {
        matrix m;
        float cosx = std::cos(x);
        float cosy = std::cos(y);
        float cosz = std::cos(z);

        float sinx = std::sin(x);
        float siny = std::sin(y);
        float sinz = std::sin(z);

        m.m[0][0] = cosy * cosz;
        m.m[0][1] = -cosy * sinz;
        m.m[0][2] = siny;

        m.m[1][0] = (siny * sinx * cosz) + (cosx * sinz);
        m.m[1][1] = -(siny * sinx * sinz) + (cosx * cosz);
        m.m[1][2] = -sinx * cosy;

        m.m[2][0] = -(cosx * siny * cosz) + (sinx * sinz);
        m.m[2][1] = (cosx * siny * sinz) + (sinx * cosz);
        m.m[2][2] = cosx * cosy;

        return m;
    }

    // Create a scaling matrix
    // Input Variables:
    // - s: Scaling factor
    // Returns the scaling matrix
    static matrix makeScale(float s) {
        matrix m;
        s = max(s, 0.01f); // Ensure scaling factor is not too small
        m.a[0] = s;
        m.a[5] = s;
        m.a[10] = s;
        return m;
    }

    // Create an identity matrix
    // Returns an identity matrix
    static matrix makeIdentity_old() {
        matrix m;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                m.m[i][j] = (i == j) ? 1.0f : 0.0f;
            }
        }
        return m;
    }
    static matrix makeIdentity() {
        matrix m;
        memset(m.a, 0, 16 * sizeof(float));
        m.m[0][0] = 1.0f;
        m.m[1][1] = 1.0f;
        m.m[2][2] = 1.0f;
        m.m[3][3] = 1.0f;
        return m;
    }

private:
    // Set all elements of the matrix to 0
    void zero() {
        memset(a, 0, 16 * sizeof(float));
    }

    // Set the matrix as an identity matrix
    void identity() {
        memset(a, 0, 16 * sizeof(float));
        m[0][0] = 1.0f;
        m[1][1] = 1.0f;
        m[2][2] = 1.0f;
        m[3][3] = 1.0f;
    }
};


