#pragma once
#include <array>
#include "../Vector4D/Vector4D.hpp"
#include "../Matrix4x4/Matrix4x4.hpp"



class Quaternion;

class AxisAngle {

public:
    double angle;
    double x;
    double y;
    double z;
    double w;


    AxisAngle(double _angle, double _x, double _y, double _z, double _w);
    AxisAngle(const Vector4D& v, double angle);
	Matrix4x4 GetRotation() const;
    Quaternion ToQuaternion();
};