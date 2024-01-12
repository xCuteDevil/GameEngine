#pragma once

#include "../../Engine/Vector4D/Vector4D.hpp"
#include "../../Engine/Matrix4x4/Matrix4x4.hpp"

class Camera{

public:
    double angle;
    double near;
    double far;
    double aspect;
    Vector4D position;
	Vector4D cameraDirection;
    
    
    bool orthographic;

    // Constructors
    Camera(double _angle, double _aspect, double _near, double _far);
    Camera();
    
    // Setters for projection matrices
    Matrix4x4 setprojectionMatrixToOrthographic(double b, double t, double l, double r, double n, double f);
    Matrix4x4 setprojectionMatrixToPerspective(double angle, double aspect, double n, double f);
    Matrix4x4 SetProjection(Matrix4x4 m);
    
    // Getters
    Matrix4x4 GetProjection() const;
    Matrix4x4 GetView() const;
    
    void UpdateViewMatrix(const Vector4D& cameraPos, const Vector4D& cameraTarget, const Vector4D& upVector = Vector4D(0.0, 1.0, 0.0, 0.0));
    

public:
    Matrix4x4 view;
    Matrix4x4 projection;
};

