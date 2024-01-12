#pragma once
#include <array>

class PolarCoords {

public:
    PolarCoords(double _radius, double _angle);
    PolarCoords NormaliseAngle();
    void PC2Cartesian(double &x, double &y);
    static PolarCoords Cartesian2PC(double x, double y);
    
    // Accessor methods
    double GetRadius() const;
    double GetAngle() const;  
    
private:
    double radius;
    double angle;
};

