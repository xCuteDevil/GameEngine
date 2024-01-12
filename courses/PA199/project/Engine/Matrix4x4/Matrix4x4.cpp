#include "Matrix4x4.hpp"
#include <cmath>
#include <array>
#include <stdexcept>

#include "../AxisAngle/AxisAngle.hpp"

using namespace std;

// Constructors
Matrix4x4::Matrix4x4(double valArray[16]) { 
    if (!valArray) {
        throw invalid_argument("Constructor argument cannot be null");
    }
    copy(valArray, valArray + 16, m);
}

Matrix4x4::Matrix4x4(double diagonalValue)
{
    for (int i = 0; i < 16; ++i) {
        m[i] = (i % 5 == 0) ? diagonalValue : 0.0;
    }
}

Matrix4x4::Matrix4x4(const Vector4D& col1, const Vector4D& col2, const Vector4D& col3, const Vector4D& col4) {
    m[0] = col1.x; m[1] = col2.x; m[2] = col3.x; m[3] = col4.x;
    m[4] = col1.y; m[5] = col2.y; m[6] = col3.y; m[7] = col4.y;
    m[8] = col1.z; m[9] = col2.z; m[10] = col3.z; m[11] = col4.z;
    m[12] = col1.w; m[13] = col2.w; m[14] = col3.w; m[15] = col4.w;
}

Matrix4x4 Matrix4x4::Identity() {
    return Matrix4x4(1.0);
}

Matrix4x4 Matrix4x4::Scale(double sx, double sy, double sz) {
    Matrix4x4 scale;
    scale.m[0] = sx;
    scale.m[5] = sy;
    scale.m[10] = sz;
    scale.m[15] = 1.0;
    return scale;
}

Matrix4x4 Matrix4x4::Rotate(double angleRad, const Vector4D& axis) {
    Vector4D nAxis = axis.UnitVector();
    double cosTheta = cos(angleRad);
    double sinTheta = sin(angleRad);

    Matrix4x4 rotation;
    rotation.m[0] = cosTheta + (1 - cosTheta) * nAxis.x * nAxis.x;
    rotation.m[1] = (1 - cosTheta) * nAxis.x * nAxis.y - sinTheta * nAxis.z;
    rotation.m[2] = (1 - cosTheta) * nAxis.x * nAxis.z + sinTheta * nAxis.y;

    rotation.m[4] = (1 - cosTheta) * nAxis.y * nAxis.x + sinTheta * nAxis.z;
    rotation.m[5] = cosTheta + (1 - cosTheta) * nAxis.y * nAxis.y;
    rotation.m[6] = (1 - cosTheta) * nAxis.y * nAxis.z - sinTheta * nAxis.x;

    rotation.m[8] = (1 - cosTheta) * nAxis.z * nAxis.x - sinTheta * nAxis.y;
    rotation.m[9] = (1 - cosTheta) * nAxis.z * nAxis.y + sinTheta * nAxis.x;
    rotation.m[10] = cosTheta + (1 - cosTheta) * nAxis.z * nAxis.z;

    rotation.m[15] = 1.0;

    return rotation;
}

Matrix4x4 Matrix4x4::Translate(double tx, double ty, double tz) {
    Matrix4x4 translation = Identity();
    translation.m[3] = tx;
    translation.m[7] = ty;
    translation.m[11] = tz;
    return translation;
}

Matrix4x4 Matrix4x4::Multiply(Matrix4x4 m2) {
    double result[16];
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            result[row * 4 + col] = 0;
            for (int i = 0; i < 4; ++i) {
                result[row * 4 + col] += m[row * 4 + i] * m2.m[i * 4 + col];
            }
        }
    }
    return Matrix4x4(result);
}

float* Matrix4x4::ToArray() const {
    float* a = new float[16];
    a[0] = static_cast<float>(m[0]); a[1] = static_cast<float>(m[4]); a[2] = static_cast<float>(m[8]); a[3] = static_cast<float>(m[12]);
    a[4] = static_cast<float>(m[1]); a[5] = static_cast<float>(m[5]); a[6] = static_cast<float>(m[9]); a[7] = static_cast<float>(m[13]);
    a[8] = static_cast<float>(m[2]); a[9] = static_cast<float>(m[6]); a[10] = static_cast<float>(m[10]); a[11] = static_cast<float>(m[14]);
    a[12] = static_cast<float>(m[3]); a[13] = static_cast<float>(m[7]); a[14] = static_cast<float>(m[11]); a[15] = static_cast<float>(m[15]);
    return a;
}

Matrix4x4 Matrix4x4::Transpose() const {
    double result[16];
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            result[col * 4 + row] = m[row * 4 + col];
        }
    }
    return Matrix4x4(result);
}

array<double, 16> Matrix4x4::ToStdArray() const{
    
    return { m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15] };
}

Vector4D Matrix4x4::MultiplyByVector(Vector4D v) {
    return Vector4D(
        m[0] * v.x + m[1] * v.y + m[2] * v.z + m[3] * v.w,
        m[4] * v.x + m[5] * v.y + m[6] * v.z + m[7] * v.w,
        m[8] * v.x + m[9] * v.y + m[10] * v.z + m[11] * v.w,
        m[12] * v.x + m[13] * v.y + m[14] * v.z + m[15] * v.w
    );
}

Matrix4x4 Matrix4x4::RMatrix2OB() {
    Vector4D v1 = Vector4D(m[0], m[4], m[8], m[12]).UnitVector();
    Vector4D v2 = Vector4D(m[1], m[5], m[9], m[13]).UnitVector();
    Vector4D v3 = Vector4D(m[2], m[6], m[10], m[14]).UnitVector();
    Vector4D v4 = Vector4D(m[3], m[7], m[11], m[15]).UnitVector();

    double vals[] = { v1.x, v2.x, v3.x, v4.x,
                 v1.y, v2.y, v3.y, v4.y,
                 v1.z, v2.z, v3.z, v4.z,
                 v1.w, v2.w, v3.w, v4.w };
    return Matrix4x4(vals);
}

Matrix4x4 Matrix4x4::OB2RMatrix(Vector4D v1, Vector4D v2, Vector4D v3, Vector4D v4) {
    if (v1.DotProduct(v2) != 0 || v1.DotProduct(v3) != 0 || v2.DotProduct(v3) != 0) {
        throw logic_error("Vectors are not orthonormal.");
    }
    
    double vals[16] = {
        v1.x, v1.y, v1.z, v1.w,
        v2.x, v2.y, v2.z, v2.w,
        v3.x, v3.y, v3.z, v3.w,
        v4.x, v4.y, v4.z, v4.w
    };
    return Matrix4x4(vals);
}

// Assignment operator
Matrix4x4& Matrix4x4::operator=(const Matrix4x4& other) {
    if (this != &other) {
        std::copy(std::begin(other.m), std::end(other.m), std::begin(m));
    }
    return *this;
}

// Matrix multiplication operator
Matrix4x4 Matrix4x4::operator*(const Matrix4x4& other) const {
    Matrix4x4 result;

    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            result.m[row * 4 + col] = 0.0;
            for (int i = 0; i < 4; ++i) {
                result.m[row * 4 + col] += m[row * 4 + i] * other.m[i * 4 + col];
            }
        }
    }
    return result;
}

// Matrix-vector multiplication operator
Vector4D Matrix4x4::operator*(const Vector4D& v) const {
    double x = m[0] * v.x + m[1] * v.y + m[2] * v.z + m[3] * v.w;
    double y = m[4] * v.x + m[5] * v.y + m[6] * v.z + m[7] * v.w;
    double z = m[8] * v.x + m[9] * v.y + m[10] * v.z + m[11] * v.w;
    double w = m[12] * v.x + m[13] * v.y + m[14] * v.z + m[15] * v.w;

    return Vector4D(x, y, z, w);
}

// Equality comparison operator
bool Matrix4x4::operator==(const Matrix4x4& other) const {
    for (int i = 0; i < 16; ++i) {
        if (m[i] != other.m[i]) return false;
    }
    return true;
}