#include "camera.hpp"
#include "Camera.hpp"
#include <array>

# define M_PI           3.14159265358979323846  /* pi */

	
// CONSTRUCTOR
Camera::Camera(double _angle, double _aspect, double _near, double _far) {
	angle = _angle;
	near = _near;
	far = _far;
	aspect = _aspect;
	projection = setprojectionMatrixToPerspective(angle, aspect, near, far);
	//projection = setprojectionMatrixToOrthographic(-4.0, 1.0, -4.0, 1.0, near, far, projection);
}

Camera::Camera() {
	angle = 45.0 * M_PI / 180.0;
	near = 0.1f;
	far = 100.0f;
	aspect = 1.0f;
	projection = setprojectionMatrixToPerspective(angle, aspect, near, far);
	//projection = setprojectionMatrixToOrthographic(-4.0, 1.0, -4.0, 1.0, near, far, projection);
}

// METHODS

/*Matrix4x4 Camera::setprojectionMatrixToPerspective(double fov, double aspect, double near, double far)
{
	double M[16] = {};

	double top = near * tan((M_PI / 180.0) * (fov * 0.5));
	double bottom = -top;
	double right = top * aspect;
	double left = -right;

	M[0] = 2 * near / (right - left);
	M[5] = 2 * near / (top - bottom);
	M[8] = (right + left) / (right - left);
	M[9] = (top + bottom) / (top - bottom);
	M[10] = -(far + near) / (far - near);
	M[11] = -1;
	M[14] = -2 * far * near / (far - near);
	M[15] = 0;

	projection = Matrix4x4(M);
	orthographic = false;

	return projection;
}*/

Matrix4x4 Camera::SetProjection(Matrix4x4 m)
{
	projection = m;
	return projection;
}


Matrix4x4 Camera::setprojectionMatrixToPerspective(double fov, double aspect, double near, double far)
{
	double M[16] = {};
	M[0] = 1 / (aspect * tan(fov * M_PI / 180.0 * 0.5));
	M[1] = 0;
	M[2] = 0;
	M[3] = 0;

	M[4] = 0;
	M[5] = 1 / (tan(fov * M_PI / 180.0 * 0.5));
	M[6] = 0;
	M[7] = 0;

	M[8] = 0;
	M[9] = 0;
	M[10] = -(far + near) / (far - near);
	M[11] = -(2 * far * near) / (far - near);

	M[12] = 0;
	M[13] = 0;
	M[14] = -1;
	M[15] = 0;

	projection = Matrix4x4(M);

	return projection;
}

Matrix4x4 Camera::setprojectionMatrixToOrthographic(double b, double t, double l, double r, double n, double f)
{
	// bottom, top, left, right, near, far
	double M[16] = {};
	M[0] = 2 / (r - l);
	M[1] = 0;
	M[2] = 0;
	M[3] = -((r + l) / (r - l));

	M[4] = 0;
	M[5] = 2 / (t - b);
	M[6] = 0;
	M[7] = -((t + b) / (t - b));

	M[8] = 0;
	M[9] = 0;
	M[10] = -2 / (f - n);
	M[11] = -((f + n) / (f - n));

	M[12] = 0; 
	M[13] = 0; 
	M[14] = 0; 
	M[15] = 1;

	projection = Matrix4x4(M);
	orthographic = false;

	return projection;
}


/* This function takes the camera position and camera target and returns corresponding view matrix */
/*Matrix4x4 CreateView(Vector4D _cameraPos, Vector4D _cameraTarget) {
	
	// Camera Direction
	Vector4D cameraDirection = _cameraTarget.SubtractVectors(_cameraPos);
	cameraDirection = cameraDirection.UnitVector();

	// Right axis
	Vector4D up = Vector4D(0.0, 1.0, 0.0, 0.0);
	Vector4D cameraRight = up.CrossProduct(cameraDirection);
	cameraRight = cameraRight.UnitVector();
	 
	// Up axis
	Vector4D cameraUp = cameraDirection.CrossProduct(cameraRight);
	cameraUp = cameraUp.UnitVector();

	// View matrix
	double viewArray[16] = { cameraRight.x, cameraRight.y, cameraRight.z, 0, cameraUp.x, cameraUp.y, cameraUp.z, 0, cameraDirection.x, cameraDirection.y, cameraDirection.z, 0, _cameraPos.x, _cameraPos.y, _cameraPos.z, 1.0 };
	Matrix4x4 view = Matrix4x4(viewArray);

	return view;
}*/

/*void Camera::UpdateViewMatrix(const Vector4D& cameraPos, const Vector4D& cameraTarget, const Vector4D& upVector) {
	cameraDirection = (cameraPos - cameraTarget).UnitVector();
	position = cameraPos;
	Vector4D cameraRight = upVector.CrossProduct(cameraDirection).UnitVector();
	Vector4D cameraUp = cameraRight.CrossProduct(cameraDirection).UnitVector();
	double viewArray[16] = {
    cameraRight.x, cameraRight.y, cameraRight.z, 0.0,
    cameraUp.x, cameraUp.y, cameraUp.z, 0.0,
    cameraDirection.x, cameraDirection.y, cameraDirection.z, 0.0,
    -cameraRight.DotProduct(cameraPos), -cameraUp.DotProduct(cameraPos), -cameraDirection.DotProduct(cameraPos), 1.0
};


	view = Matrix4x4(viewArray);
}*/

void Camera::UpdateViewMatrix(const Vector4D& cameraPos, const Vector4D& cameraTarget, const Vector4D& worldUp) {
	
	// Must be pointing in the reverse direction of what it is targeting
	cameraDirection = (cameraPos - cameraTarget).UnitVector();

	// In case camera direction is parallel to world up vector (orthographic camera)
	/*Vector4D upVector = Vector4D(0.0, 1.0, 0.0001, 0.0);
	if (!worldUp.IsParallel(cameraDirection)) {
		upVector = worldUp;
	}*/
		

	
	Vector4D cameraRight = worldUp.CrossProduct(cameraDirection).UnitVector();
	Vector4D cameraUp = cameraDirection.CrossProduct(cameraRight).UnitVector();
	double viewArray[16] = {
	cameraRight.x, cameraRight.y, cameraRight.z, 0,
	cameraUp.x, cameraUp.y, cameraUp.z, 0,
	cameraDirection.x, cameraDirection.y, cameraDirection.z, 0,
	0.0, 0.0, 0.0, 1.0
	};


	view = Matrix4x4(viewArray);
}

Matrix4x4 Camera::GetView() const {
	return view;
}

Matrix4x4 Camera::GetProjection() const
{
	return projection;
}
