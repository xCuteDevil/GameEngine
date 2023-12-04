#pragma once

#include <array>
#include "../Matrix4x4/Matrix4x4.hpp"

# define M_PI           3.14159265358979323846  /* pi */

// Forward declaration
class AxisAngle;

class Quaternion {
public:
    // Constructor
    Quaternion(double _real, double _i, double _j, double _k);

    // Basic Quaternion operations
    Quaternion Conjugate() const;
    double Magnitude() const;
    Quaternion Normalise() const;
    Quaternion Inverse() const;
    double DotProduct(const Quaternion& q2) const;
    Quaternion Logarithm() const;

    // Conversion methods
    Matrix4x4 ToRotationMatrix() const;
    AxisAngle ToAxisAngle();
    std::array<double, 3> ToEulerAngles() const;

    // Quaternion Interpolation
    static Quaternion Slerp(const Quaternion& q1, const Quaternion& q2, double t);

    /* Operator overloads */
    Quaternion operator*(const Quaternion& q2) const;
    Quaternion operator+(const Quaternion& q2) const;
    Quaternion operator-(const Quaternion& q2) const;

    // Static factory methods
    static Quaternion Identity();
    static Quaternion FromEulerAngles(double roll, double pitch, double yaw);

    // Utility methods
    Quaternion Scale(double scalar) const;

    // Getters for member variables
    double GetRealPart() const { return real; }
    double GetI() const { return i; }
    double GetJ() const { return j; }
    double GetK() const { return k; }

private:
    // Member variables
    double real;
    double i;
    double j;
    double k;
};
