#include "PolarCoords.hpp"
#include <math.h>

# define M_PI           3.14159265358979323846  /* pi */

PolarCoords::PolarCoords(double _radius, double _angle) { // Constructor
    radius = _radius;
    angle = _angle; // pi radians
}

// Getters
double PolarCoords::GetRadius() const {
    return radius;
}

double PolarCoords::GetAngle() const {
    return angle;
}

PolarCoords PolarCoords::NormaliseAngle()
{
    double normalized_angle = fmod(angle, 2 * M_PI);
    if (normalized_angle < 0) {
        normalized_angle += 2 * M_PI;
    }
    return PolarCoords(radius, normalized_angle);
}

void PolarCoords::PC2Cartesian(double &x, double &y)
{
    x = radius * cos(angle);
    y = radius * sin(angle);
    return; 
}

 PolarCoords PolarCoords::Cartesian2PC(double x, double y)
{
   double r = sqrt(x * x + y * y);
   double theta = atan2(y, x);
   return PolarCoords(r, theta);
}