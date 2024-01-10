#include "Vector4D.hpp"
#include "Vector4D.hpp"
#include <cmath>
#include <stdexcept>

using namespace std;

Vector4D::Vector4D(double _x, double _y, double _z, double _w) : x(_x), y(_y), z(_z), w(_w) {}

Vector4D::Vector4D(double _x, double _y, double _z) : x(_x), y(_y), z(_z), w(0.0) {}

Vector4D::Vector4D() : x(0), y(0), z(0), w(0) {}

Vector4D Vector4D::UnitVector() const {
    double magnitude = MagnitudeOfAVector();
    if (magnitude < 1e-8) {
        throw std::runtime_error("Cannot calculate the unit vector of a zero vector.");
    }
    return { x / magnitude, y / magnitude, z / magnitude, 0 }; // Normalize all components
}

double Vector4D::MagnitudeOfAVector() const {
    return sqrt(x * x + y * y + z * z);
}

Vector4D Vector4D::OppositeVector() const {
    return { -x, -y, -z, -w }; // Assuming we want to negate all components including w
}

Vector4D Vector4D::AddVectors(const Vector4D& addend) const {
    return { x + addend.x, y + addend.y, z + addend.z, w + addend.w };
}

Vector4D Vector4D::SubtractVectors(const Vector4D& subtrahend) const {
    return { x - subtrahend.x, y - subtrahend.y, z - subtrahend.z, w - subtrahend.w };
}

double Vector4D::DotProduct(const Vector4D& vec) const {
    return x * vec.x + y * vec.y + z * vec.z + w * vec.w;
}

Vector4D Vector4D::CrossProduct(const Vector4D& vec) const {
    if (w != 0.0 || vec.w != 0.0) {
        throw std::invalid_argument("Cross product expects both w components to be 0.");
    }

    double newX = y * vec.z - z * vec.y;
    double newY = z * vec.x - x * vec.z;
    double newZ = x * vec.y - y * vec.x;

    return Vector4D(newX, newY, newZ, 0.0);
}


array<double, 4> Vector4D::ToStdArray() const {
    return { x, y, z, w };
}

double Vector4D::AngleBetween(const Vector4D& v1, const Vector4D& v2) {
    double dot = v1.DotProduct(v2);
    double magV1 = v1.MagnitudeOfAVector();
    double magV2 = v2.MagnitudeOfAVector();
    return acos(dot / (magV1 * magV2));
}

bool Vector4D::IsNormalized() const {
    const double epsilon = 1e-5; // Tolerance for floating-point comparison
    return std::abs(1.0 - MagnitudeOfAVector()) < epsilon;
}

bool Vector4D::IsParallel(const Vector4D& v1) const {
    // Checks whether the vectors are parallel
	// If the angle between them is 0 or PI, they are parallel
	Vector4D v1Normalised = v1.UnitVector();
	Vector4D v2Normalised = Vector4D(x,y,z,w).UnitVector();
	double dot = v1Normalised.DotProduct(v2Normalised);
    
    return std::abs(dot) > 0.999999f;
}

/* Operator overloads */

Vector4D Vector4D::operator+(const Vector4D& addend) const {
    return Vector4D(x + addend.x, y + addend.y, z + addend.z, w + addend.w);
}

Vector4D Vector4D::operator-(const Vector4D& subtrahend) const {
    return Vector4D(x - subtrahend.x, y - subtrahend.y, z - subtrahend.z, w - subtrahend.w);
}

// Scalar multiplication
Vector4D Vector4D::operator*(double scalar) const {
    return Vector4D(x * scalar, y * scalar, z * scalar, w * scalar);
}

// Scalar division
Vector4D Vector4D::operator/(double scalar) const {
    if (scalar == 0) {
        throw std::invalid_argument("Division by zero exception.");
    }
    return Vector4D(x / scalar, y / scalar, z / scalar, w / scalar);
}

float* Vector4D::ToArray() const {
    float* array = new float[4];
    array[0] = static_cast<float>(x);
    array[1] = static_cast<float>(y);
    array[2] = static_cast<float>(z);
    array[3] = static_cast<float>(w);
    return array;
}

bool Vector4D::Equals(const Vector4D& v1, const Vector4D& v2, int precision) {
    double tolerance = std::pow(10.0, -precision);
    return std::abs(v1.x - v2.x) < tolerance &&
        std::abs(v1.y - v2.y) < tolerance &&
        std::abs(v1.z - v2.z) < tolerance &&
        std::abs(v1.w - v2.w) < tolerance;
}