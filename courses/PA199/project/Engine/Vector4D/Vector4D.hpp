#pragma once

#include <array>

class Vector4D {
public:
    double x, y, z, w;

    Vector4D(double _x, double _y, double _z, double _w);
    Vector4D(double _x, double _y, double _z);
    Vector4D();

    Vector4D UnitVector() const;
    double MagnitudeOfAVector() const;
    Vector4D OppositeVector() const;
    Vector4D AddVectors(const Vector4D& addend) const;
    Vector4D SubtractVectors(const Vector4D& subtrahend) const;
    double DotProduct(const Vector4D& vec) const;
    Vector4D CrossProduct(const Vector4D& vec) const;
    std::array<double, 4> ToStdArray() const;
    double AngleBetween(const Vector4D& v1, const Vector4D& v2);
    Vector4D Homogenize() const;
    bool IsNormalized() const;
    bool IsParallel(const Vector4D& v1) const;
    static bool Equals(const Vector4D& v1, const Vector4D& v2, int precision);

	/* Operator overloads */
    Vector4D operator+(const Vector4D& addend) const; // Vector addition
    Vector4D operator-(const Vector4D& subtrahend) const; // Vector subtraction
    Vector4D operator*(double scalar) const;       // Scalar multiplication
    Vector4D operator/(double scalar) const;       // Scalar division
    bool operator<(const Vector4D& rhs) const;
        
    float* ToArray() const;
};