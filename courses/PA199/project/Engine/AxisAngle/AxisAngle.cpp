#include <cmath>
#include "AxisAngle.hpp"
#include "../Quaternion/Quaternion.hpp"
#include "../Matrix4x4/Matrix4x4.hpp"

AxisAngle::AxisAngle(double _angle, double _x, double _y, double _z, double _w)
{
	angle = _angle;
	x = _x;
	y = _y;
	z = _z;
	w = _w;
}

AxisAngle::AxisAngle(const Vector4D& v, double angle)
{
	this->angle = angle;
	x = v.x;
	y = v.y;
	z = v.z;
	w = v.w;
}

Quaternion AxisAngle::ToQuaternion() {
	double magnitude = sqrt(x * x + y * y + z * z);
	double normX = x / magnitude;
	double normY = y / magnitude;
	double normZ = z / magnitude;
	
	double halfAngle = angle * 0.5;
	double s = sin(halfAngle);
	double qw = cos(halfAngle);
	double qx = normX * s;
	double qy = normY * s;
	double qz = normZ * s;

	return Quaternion(qw, qx, qy, qz);
}

Matrix4x4 AxisAngle::GetRotation() const {
    // Using Rodrigues' rotation formula

    double c = cos(angle);
    double s = sin(angle);
    double t = 1.0 - c;

    // Components of the normalized axis
    double nx = x;
    double ny = y;
    double nz = z;

    // Components of the matrix
    double r[16] = {
        t * nx * nx + c,        t * nx * ny - s * nz,  t * nx * nz + s * ny, 0.0,
        t * nx * ny + s * nz,  t * ny * ny + c,        t * ny * nz - s * nx, 0.0,
        t * nx * nz - s * ny,  t * ny * nz + s * nx,  t * nz * nz + c,      0.0,
        0.0,                   0.0,                   0.0,                  1.0
    };

    return Matrix4x4(r);
}


