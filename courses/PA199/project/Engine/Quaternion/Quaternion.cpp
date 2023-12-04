#include <cmath>
#include "Quaternion.hpp"
#include "../AxisAngle/AxisAngle.hpp"
#include <stdexcept>

# define M_PI           3.14159265358979323846  /* pi */

using namespace std;

Quaternion::Quaternion(double _real, double _i, double _j, double _k)
{
	// a + bi + cj + dk
	real = _real;
	i = _i;
	j = _j;
	k = _k;
}

Quaternion Quaternion::Conjugate() const
{
	return Quaternion(real, -i, -j, -k);
}

double Quaternion::Magnitude() const
{
	return sqrt(pow(GetRealPart(), 2) + pow(GetI(), 2) + pow(GetJ(), 2) + pow(GetK(), 2));
}

Quaternion Quaternion::Inverse() const
{
	// q_inverse = conjugate / magnitude^2
	Quaternion q_conjugate = this->Conjugate();

	// Calculate magnitude squared
	double q_magnitudeSquared = this->real * this->real +
		this->i * this->i +
		this->j * this->j +
		this->k * this->k;

	return Quaternion(q_conjugate.real / q_magnitudeSquared,
		q_conjugate.i / q_magnitudeSquared,
		q_conjugate.j / q_magnitudeSquared,
		q_conjugate.k / q_magnitudeSquared);
}

double Quaternion::DotProduct(const Quaternion& q2) const
{
    return real * q2.real + i * q2.i + j * q2.j + k * q2.k;
}

Quaternion Quaternion::Normalise() const
{
	/*
	*	q_normalised = q/magnitude
	*/
	double q_magnitude = this->Magnitude();
	if (q_magnitude < 1e-8) {
		throw std::runtime_error("Cannot normalize a zero-magnitude quaternion.");
	}
	return Quaternion(this->real/q_magnitude, this->i /q_magnitude, this->j /q_magnitude, this->k / q_magnitude);
}

Matrix4x4 Quaternion::ToRotationMatrix() const {
	double a11 = 1 - 2 * (j * j) - 2 * (k * k);
	double a12 = 2 * i * j - 2 * real * k;
	double a13 = 2 * i * k + 2 * real * j;
	double a14 = 0;

	double a21 = 2 * i * j + 2 * real * k;
	double a22 = 1 - 2 * (i * i) - 2 * (k * k);
	double a23 = 2 * j * k - 2 * real * i;
	double a24 = 0;

	double a31 = 2 * i * k - 2 * real * j;
	double a32 = 2 * j * k + 2 * real * i;
	double a33 = 1 - 2 * (i * i) - 2 * (j * j);
	double a34 = 0;

	double a41 = 0;
	double a42 = 0;
	double a43 = 0;
	double a44 = 1.0;

	double mResArray[16] = { a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44 };
	return Matrix4x4(mResArray);
}

Quaternion Quaternion::Identity() {
	// Returns the identity quaternion
	return Quaternion(1.0, 0.0, 0.0, 0.0);
}

Quaternion Quaternion::operator*(const Quaternion& q2) const {
	// Hamilton product of two quaternions
	return Quaternion(
		real * q2.real - i * q2.i - j * q2.j - k * q2.k,
		real * q2.i + i * q2.real + j * q2.k - k * q2.j,
		real * q2.j - i * q2.k + j * q2.real + k * q2.i,
		real * q2.k + i * q2.j - j * q2.i + k * q2.real
	);
}

Quaternion Quaternion::operator+(const Quaternion& q2) const {
	return Quaternion(real + q2.real, i + q2.i, j + q2.j, k + q2.k);
}

Quaternion Quaternion::operator-(const Quaternion& q2) const {
	return Quaternion(real - q2.real, i - q2.i, j - q2.j, k - q2.k);
}

Quaternion Quaternion::FromEulerAngles(double roll, double pitch, double yaw) {
	// Convert Euler angles to quaternion
	double cy = cos(yaw * 0.5);
	double sy = sin(yaw * 0.5);
	double cp = cos(pitch * 0.5);
	double sp = sin(pitch * 0.5);
	double cr = cos(roll * 0.5);
	double sr = sin(roll * 0.5);

	return Quaternion(
		cr * cp * cy + sr * sp * sy,
		sr * cp * cy - cr * sp * sy,
		cr * sp * cy + sr * cp * sy,
		cr * cp * sy - sr * sp * cy
	);
}

std::array<double, 3> Quaternion::ToEulerAngles() const {
    std::array<double, 3> angles;

    // Roll (x-axis rotation)
    double sinr_cosp = 2 * (real * i + j * k);
    double cosr_cosp = 1 - 2 * (i * i + j * j);
    angles[0] = std::atan2(sinr_cosp, cosr_cosp);

    // Pitch (y-axis rotation)
    double sinp = 2 * (real * j - k * i);
    if (std::abs(sinp) >= 1)
        angles[1] = std::copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    else
        angles[1] = std::asin(sinp);

    // Yaw (z-axis rotation)
    double siny_cosp = 2 * (real * k + i * j);
    double cosy_cosp = 1 - 2 * (j * j + k * k);
    angles[2] = std::atan2(siny_cosp, cosy_cosp);

    return angles;
}


Quaternion Quaternion::Scale(double scalar) const {
	// Scale the quaternion by a scalar
	return Quaternion(real * scalar, i * scalar, j * scalar, k * scalar);
}

Quaternion Quaternion::Logarithm() const {

	double magnitude = this->Magnitude();
	double vectorMagnitude = sqrt(i * i + j * j + k * k);
	double logReal = log(magnitude);
	double logI = 0.0, logJ = 0.0, logK = 0.0;
	
	if (vectorMagnitude > 1e-8) {
		double factor = acos(real / magnitude) / vectorMagnitude;
		
		logI = factor * i;
		logJ = factor * j;
		logK = factor * k;
	}

	return Quaternion(logReal, logI, logJ, logK);
}

Quaternion Quaternion::Slerp(const Quaternion& q1, const Quaternion& q2, double t) {
	double dot = q1.DotProduct(q2);
	if (dot > 0.9995) {
		Quaternion result = q1 + (q2 - q1).Scale(t);
		return result.Normalise();
	}

	if (dot < 0.0) {
		dot = -dot;
	}

	double theta_0 = acos(dot);
	double theta = theta_0 * t;
	double sin_theta = sin(theta);
	double sin_theta_0 = sin(theta_0);

	double s0 = cos(theta) - dot * sin_theta / sin_theta_0;
	double s1 = sin_theta / sin_theta_0;

	return (q1.Scale(s0) + q2.Scale(s1)).Normalise();
}



AxisAngle Quaternion::ToAxisAngle() {
	double norm = sqrt(real * real + i * i + j * j + k * k);
	double normalized_real = real / norm;
	double normalized_i = i / norm;
	double normalized_j = j / norm;
	double normalized_k = k / norm;

	double angle = 2 * acos(normalized_real);
	double s = sqrt(1 - normalized_real * normalized_real);
	
	const double threshold = 1e-8;

	double x, y, z;

	if (s < threshold) {
		x = normalized_i;
		y = normalized_j;
		z = normalized_k;
	}
	else {
		x = normalized_i / s;
		y = normalized_j / s;
		z = normalized_k / s;
	}
	
	return AxisAngle(angle, x, y, z, 0.0);
}
