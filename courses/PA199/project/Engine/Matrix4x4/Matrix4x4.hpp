#pragma once

#include <array>
#include "../Vector4D/Vector4D.hpp"

class Matrix4x4{
public:
    // Constructors
    explicit Matrix4x4(double valArray[16]);
    explicit Matrix4x4(double diagonalValue = 0.0);
    Matrix4x4(const Vector4D& col1, const Vector4D& col2, const Vector4D& col3, const Vector4D& col4);

    // Basic Matrix Operations
    Matrix4x4 Multiply(Matrix4x4 m2);
    Matrix4x4 Transpose() const;
    
    // Utility Functions
    float* ToArray() const;
    std::array<double, 16> ToStdArray() const;
    
    // Static utility functions for transformations
    static Matrix4x4 Identity();
    static Matrix4x4 Scale(double sx, double sy, double sz);
    static Matrix4x4 Rotate(double angleRad, const Vector4D& axis);
    static Matrix4x4 Translate(double tx, double ty, double tz);
    
    // Conversion between Rotation Matrix and Orthonormal Basis
    Matrix4x4 RMatrix2OB();
    Matrix4x4 OB2RMatrix(Vector4D v1, Vector4D v2, Vector4D v3, Vector4D v4);

    // Vector Transformation
    Vector4D MultiplyByVector(Vector4D v);

    // Operator overloads
    Matrix4x4& operator=(const Matrix4x4& other);
    Matrix4x4 operator*(const Matrix4x4& other) const;
    Vector4D operator*(const Vector4D& v) const;
    bool operator==(const Matrix4x4& other) const;
    
private:
    double m[16]; // Represents a 4x4 matrix
};