// This is a tutorial file. Feel free to remove it.

//#include "../some_my_library_dir/some_library_file.hpp"
#include <gtest/gtest.h>
#include "../Engine/Vector4D/Vector4D.hpp"
#include "../Engine/Matrix4x4/Matrix4x4.hpp"
#include "../Engine/PolarCoords/PolarCoords.hpp"
#include "../Engine/Quaternion/Quaternion.hpp"
#include "../Engine/AxisAngle/AxisAngle.hpp"
#include <array>
#include <cmath>

#define M_PI 3.14159265358979323846

// Helper function to test matrix equality
void ExpectMatrixEqual(const Matrix4x4& matrix, const double expected[16], const std::string& message) {
    for (int i = 0; i < 16; ++i) {
        EXPECT_EQ(matrix.ToStdArray()[i], expected[i]) << message << " Matrix element at index " << i << " should be " << expected[i];
    }
}


TEST(some_test_suite, unitVector_test)
{
    Vector4D vec(3, 0, 0, 0);
    Vector4D vecUnit = vec.::Vector4D::UnitVector();
    Vector4D expectedRes(1, 0, 0, 0);

    EXPECT_EQ(vecUnit.x, expectedRes.x) << "1";
    EXPECT_EQ(vecUnit.y, expectedRes.y) << "0";
    EXPECT_EQ(vecUnit.z, expectedRes.z) << "0";
    EXPECT_EQ(vecUnit.w, expectedRes.w) << "0";
}

TEST(some_test_suite, magnitudeOfAVector_test) {
    Vector4D vec(2, 3, 4, 0);
    double expectedMagnitude = 5.385164; // sqrt(2^2 + 3^2 + 4^2)
    EXPECT_NEAR(vec.MagnitudeOfAVector(), expectedMagnitude, 1e-5) << "Magnitude should be close to " << expectedMagnitude;
}

TEST(some_test_suite, oppositeVector_test)
{
    Vector4D vec(2, 3, 4, 0);
    Vector4D vecOpposite = vec.::Vector4D::OppositeVector();
    Vector4D expectedRes(-2, -3, -4, 0);

    EXPECT_EQ(vecOpposite.x, expectedRes.x) << "-2";
    EXPECT_EQ(vecOpposite.y, expectedRes.y) << "-3";
    EXPECT_EQ(vecOpposite.z, expectedRes.z) << "-4";
    EXPECT_EQ(vecOpposite.w, expectedRes.w) << "0";
}

TEST(some_test_suite, addVectors_test)
{
    Vector4D vec1(2, 3, 4, 0);
    Vector4D vec2(1, 2, 3, 0);
    Vector4D res = vec1.::Vector4D::AddVectors(vec2);
    Vector4D expectedRes(3, 5, 7, 0);

    EXPECT_EQ(res.x, expectedRes.x) << "3";
    EXPECT_EQ(res.y, expectedRes.y) << "5";
    EXPECT_EQ(res.z, expectedRes.z) << "7";
    EXPECT_EQ(res.w, expectedRes.w) << "0";
}

TEST(some_test_suite, subtractVectors_test)
{
    Vector4D vec1(2, 3, 4, 0);
    Vector4D vec2(1, 2, 3, 0);
    Vector4D res = vec1.::Vector4D::SubtractVectors(vec2);
    Vector4D expectedRes(1, 1, 1, 0);

    EXPECT_EQ(res.x, expectedRes.x) << "1";
    EXPECT_EQ(res.y, expectedRes.y) << "1";
    EXPECT_EQ(res.z, expectedRes.z) << "1";
    EXPECT_EQ(res.w, expectedRes.w) << "0";
}

TEST(some_test_suite, dotProduct_test)
{
    Vector4D vec1(1, 3, -5, 0);
    Vector4D vec2(4, -2, -1, 0);
    double dotProduct = vec1.::Vector4D::DotProduct(vec2);
    double expectedRes = 3;

    EXPECT_EQ(dotProduct, expectedRes) << "3";
}

TEST(some_test_suite, crossProduct_test)
{
    Vector4D vec1(0, 1, 2, 0);
    Vector4D vec2(-5, 0, 4, 0);
    Vector4D res = vec1.::Vector4D::CrossProduct(vec2);
    Vector4D expectedRes(4, -10, 5, 0);

    EXPECT_EQ(res.x, expectedRes.x) << "4";
    EXPECT_EQ(res.y, expectedRes.y) << "-10";
    EXPECT_EQ(res.z, expectedRes.z) << "5";
    EXPECT_EQ(res.w, expectedRes.w) << "1";
}

TEST(matrix_suite, ArrayConstructor) {
    double vals[16] = {
        1.0, 2.0, 3.0, 4.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0
    };
    Matrix4x4 matrix(vals);

    ExpectMatrixEqual(matrix, vals, "ArrayConstructor");
}

TEST(matrix_suite, ConstructorWithDiagonalValue) {
    Matrix4x4 matrix(5.0);
    double expected[16] = {
        5, 0, 0, 0,
        0, 5, 0, 0,
        0, 0, 5, 0,
        0, 0, 0, 5
    };

    ExpectMatrixEqual(matrix, expected, "ConstructorWithDiagonalValue");
}

TEST(matrix_suite, ConstructorWithVectors) {
    Vector4D col1(1, 0, 0, 0);
    Vector4D col2(0, 1, 0, 0);
    Vector4D col3(0, 0, 1, 0);
    Vector4D col4(0, 0, 0, 1);
    Matrix4x4 matrix(col1, col2, col3, col4);

    double expected[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };

    ExpectMatrixEqual(matrix, expected, "ConstructorWithVectors");
}

TEST(matrix_suite, IdentityMatrix) {
    Matrix4x4 identity = Matrix4x4::Identity();
    double expected[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };

    ExpectMatrixEqual(identity, expected, "IdentityMatrix");
}

TEST(matrix_suite, ScaleMatrix) {
    Matrix4x4 scale = Matrix4x4::Scale(2, 3, 4);
    double expected[16] = {
        2, 0, 0, 0,
        0, 3, 0, 0,
        0, 0, 4, 0,
        0, 0, 0, 1
    };

    ExpectMatrixEqual(scale, expected, "ScaleMatrix");
}

TEST(matrix_suite, RotationMatrix) {
    double angleRad = M_PI / 4; // 45 degrees
    Vector4D axis(0, 0, 1, 0);
    Matrix4x4 rotation = Matrix4x4::Rotate(angleRad, axis);
    double cosTheta = std::cos(angleRad);
    double sinTheta = std::sin(angleRad);
    double expected[16] = {
        cosTheta, -sinTheta, 0, 0,
        sinTheta, cosTheta, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };

    ExpectMatrixEqual(rotation, expected, "RotationMatrix");
}

TEST(matrix_suite, TranslationMatrix) {
    Matrix4x4 translation = Matrix4x4::Translate(1, 2, 3);
    double expected[16] = {
        1, 0, 0, 1,
        0, 1, 0, 2,
        0, 0, 1, 3,
        0, 0, 0, 1
    };

    ExpectMatrixEqual(translation, expected, "TranslationMatrix");
}

TEST(matrix_suite, MultiplyByMatrix) {
    double m1Array[16] = { 5, 2, 6, 1, 0, 6, 2, 0, 3, 8, 1, 4, 1, 8, 5, 6 };
    double m2Array[16] = { 7, 5, 8, 0, 1, 8, 2, 6, 9, 4, 3, 8, 5, 3, 7, 9 };
    std::array<double, 16> expectedArray{ 96, 68, 69, 69, 24, 56, 18, 52, 58, 95, 71, 92, 90, 107, 81, 142 };
    Matrix4x4 m1(m1Array);
    Matrix4x4 m2(m2Array);
    Matrix4x4 res = m1.Multiply(m2);
    auto resArray = res.ToStdArray();
    for (int i = 0; i < 16; ++i) {
        EXPECT_NEAR(resArray[i], expectedArray[i], 1e-5) << "Matrix multiplication result at index " << i << " should be close to " << expectedArray[i];
    }
}

TEST(matrix_suite, ToArray) {
    double m1Array[16] = { 5, 2, 6, 1, 0, 6, 2, 0, 3, 8, 1, 4, 1, 8, 5, 6 };
    float expectedArray[16] = { 5, 2, 6, 1, 0, 6, 2, 0, 3, 8, 1, 4, 1, 8, 5, 6 };
    Matrix4x4 m1(m1Array);
    const float* array = m1.ToArray();
    for (int i = 0; i < 16; i++) {
        EXPECT_EQ(m1.ToStdArray()[i], expectedArray[i]) << "Matrix element at index " << i << " should be " << expectedArray[i];
    }
}

TEST(matrix_suite, ToStdArray)
{
    double m1Array[16] = { 5,2,6,1,0,6,2,0,3,8,1,4,1,8,5,6 };
    float expectedArray[16] = { 5,2,6,1,0,6,2,0,3,8,1,4,1,8,5,6 };
    Matrix4x4 m1(m1Array);
    for (int i = 0; i < 16; i++) {
        EXPECT_EQ(m1.ToStdArray()[i], expectedArray[i]);
    }
}

TEST(matrix_suite, Transpose) {
    double mArray[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    double expectedArray[16] = {1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16};
    std::array<double, 16> expectedMArray{ 1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16 };
    Matrix4x4 m1(mArray);
    Matrix4x4 res = m1.::Matrix4x4::Transpose();

    EXPECT_EQ(res.ToStdArray(), expectedMArray);
}

TEST(matrix_suite, MultiplyByVector) {
    double mArray[16] = { 1,0,2,0,0,3,0,4,0,0,5,0,6,0,0,7 };
    Matrix4x4 m1(mArray);

    Vector4D v(2, 5, 1, 8);
    std::array<double, 4> expectedVArray{ 4,47,5,68 };
    Vector4D res = m1.MultiplyByVector(v);
    EXPECT_EQ(res.ToStdArray(), expectedVArray);
}

// QUATERNIONS

TEST(quaternion_suite, Conjugate) {
    Quaternion q(1.0, 2.0, 3.0, 4.0);
    Quaternion conjugate = q.Conjugate();

    EXPECT_DOUBLE_EQ(conjugate.GetRealPart(), 1.0) << "Real part of conjugate should be 1.0";
    EXPECT_DOUBLE_EQ(conjugate.GetI(), -2.0) << "I component of conjugate should be -2.0";
    EXPECT_DOUBLE_EQ(conjugate.GetJ(), -3.0) << "J component of conjugate should be -3.0";
    EXPECT_DOUBLE_EQ(conjugate.GetK(), -4.0) << "K component of conjugate should be -4.0";
}


TEST(quaternion_suite, Magnitude) {
    Quaternion q(1.0, 2.0, 3.0, 4.0);
    
    double magnitude = q.Magnitude();
    
    double expectedMagnitude = std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0 + 4.0 * 4.0);
    EXPECT_NEAR(magnitude, expectedMagnitude, 1e-6) << "Magnitude should be " << expectedMagnitude;
}

TEST(quaternion_suite, Normalise) {
    Quaternion q(1.0, 2.0, 3.0, 4.0);

    Quaternion normalised = q.Normalise();
    double mag = std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0 + 4.0 * 4.0);

    EXPECT_NEAR(normalised.GetRealPart(), 1.0 / mag, 1e-6) << "Real part of normalised should be " << 1.0 / mag;
    EXPECT_NEAR(normalised.GetI(), 2.0 / mag, 1e-6) << "I component of normalised should be " << 2.0 / mag;
    EXPECT_NEAR(normalised.GetJ(), 3.0 / mag, 1e-6) << "J component of normalised should be " << 3.0 / mag;
    EXPECT_NEAR(normalised.GetK(), 4.0 / mag, 1e-6) << "K component of normalised should be " << 4.0 / mag;
}

TEST(quaternion_suite, Inverse) {
    Quaternion q(1.0, 2.0, 3.0, 4.0);

    Quaternion inverse = q.Inverse();
    double magSquared = 1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0 + 4.0 * 4.0;

    EXPECT_NEAR(inverse.GetRealPart(), 1.0 / magSquared, 1e-6) << "Real part of inverse should be " << 1.0 / magSquared;
    EXPECT_NEAR(inverse.GetI(), -2.0 / magSquared, 1e-6) << "I component of inverse should be " << -2.0 / magSquared;
    EXPECT_NEAR(inverse.GetJ(), -3.0 / magSquared, 1e-6) << "J component of inverse should be " << -3.0 / magSquared;
    EXPECT_NEAR(inverse.GetK(), -4.0 / magSquared, 1e-6) << "K component of inverse should be " << -4.0 / magSquared;
}

TEST(quaternion_suite, DotProduct) {
    Quaternion q1(1.0, 2.0, 3.0, 4.0);
    Quaternion q2(2.0, 3.0, 4.0, 1.0);

    double dotProduct = q1.DotProduct(q2);
    double expectedDotProduct = 1.0 * 2.0 + 2.0 * 3.0 + 3.0 * 4.0 + 4.0 * 1.0;

    EXPECT_DOUBLE_EQ(dotProduct, expectedDotProduct) << "Dot product should be " << expectedDotProduct;
}

TEST(quaternion_suite, Logarithm) {
    Quaternion q(1.0, 0.0, 0.0, 0.0);

    Quaternion logQ = q.Logarithm();
    
    double expectedReal = 0.0;
    std::array<double, 3> expectedImag = { 0.0, 0.0, 0.0 };

    EXPECT_NEAR(logQ.GetRealPart(), expectedReal, 1e-6) << "Real part of logarithm should be " << expectedReal;
    EXPECT_NEAR(logQ.GetI(), expectedImag[0], 1e-6) << "I component of logarithm should be " << expectedImag[0];
    EXPECT_NEAR(logQ.GetJ(), expectedImag[1], 1e-6) << "J component of logarithm should be " << expectedImag[1];
    EXPECT_NEAR(logQ.GetK(), expectedImag[2], 1e-6) << "K component of logarithm should be " << expectedImag[2];
}

TEST(quaternion_suite, ToRotationMatrix) {
    Quaternion q(1.0, 0.0, 0.0, 0.0);
    Matrix4x4 matrix = q.ToRotationMatrix();
    
    Matrix4x4 expectedMatrix = Matrix4x4::Identity();

    for (int i = 0; i < 16; ++i) {
        EXPECT_NEAR(matrix.ToStdArray()[i], expectedMatrix.ToStdArray()[i], 1e-6)
            << "Matrix element at index " << i << " should be " << expectedMatrix.ToStdArray()[i];
    }
}

TEST(quaternion_suite, ToEulerAngles) {
    Quaternion q = Quaternion::Identity();
    std::array<double, 3> eulerAngles = q.ToEulerAngles();

    std::array<double, 3> expectedAngles = { 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(eulerAngles[i], expectedAngles[i], 1e-6)
            << "Euler angle at index " << i << " should be " << expectedAngles[i];
    }
}

TEST(quaternion_suite, FromEulerAngles) {
    double roll = M_PI / 6;  // 30 degrees
    double pitch = M_PI / 4; // 45 degrees
    double yaw = M_PI / 3;   // 60 degrees

    Quaternion quaternion = Quaternion::FromEulerAngles(roll, pitch, yaw);

    auto angles = quaternion.ToEulerAngles();

    EXPECT_NEAR(angles[0], roll, 1e-6) << "Roll angle should be close to " << roll;
    EXPECT_NEAR(angles[1], pitch, 1e-6) << "Pitch angle should be close to " << pitch;
    EXPECT_NEAR(angles[2], yaw, 1e-6) << "Yaw angle should be close to " << yaw;
}

TEST(quaternion_suite, Slerp) {
    Quaternion q1 = Quaternion::Quaternion(0.966,0,0.259,0);
    Quaternion q2 = Quaternion::Quaternion(0.866, 0, 0, 0.5);
    double t = 0.33;
    
    Quaternion slerped = Quaternion::Slerp(q1, q2, t);
    
    Quaternion expected = Quaternion::Quaternion(0.968463539816881, 0, 0.1790642406252351, 0.17355792033802755);

    EXPECT_NEAR(slerped.GetRealPart(), expected.GetRealPart(), 1e-4);
    EXPECT_NEAR(slerped.GetI(), expected.GetI(), 1e-4);
    EXPECT_NEAR(slerped.GetJ(), expected.GetJ(), 1e-4);
    EXPECT_NEAR(slerped.GetK(), expected.GetK(), 1e-4);
}


TEST(quaternion_suite, Multiplication) {
    Quaternion q1(1.0, 2.0, 3.0, 4.0);
    Quaternion q2(5.0, 6.0, 7.0, 8.0);

    Quaternion product = q1 * q2;
    Quaternion expected(-60, 12, 30, 24);

    EXPECT_NEAR(product.GetRealPart(), expected.GetRealPart(), 1e-6);
    EXPECT_NEAR(product.GetI(), expected.GetI(), 1e-6);
    EXPECT_NEAR(product.GetJ(), expected.GetJ(), 1e-6);
    EXPECT_NEAR(product.GetK(), expected.GetK(), 1e-6);
}

TEST(quaternion_suite, Addition) {
    Quaternion q1(1.0, 2.0, 3.0, 4.0);
    Quaternion q2(5.0, 6.0, 7.0, 8.0);

    Quaternion sum = q1 + q2;
    Quaternion expected(6.0, 8.0, 10.0, 12.0);

    EXPECT_DOUBLE_EQ(sum.GetRealPart(), expected.GetRealPart());
    EXPECT_DOUBLE_EQ(sum.GetI(), expected.GetI());
    EXPECT_DOUBLE_EQ(sum.GetJ(), expected.GetJ());
    EXPECT_DOUBLE_EQ(sum.GetK(), expected.GetK());
}

TEST(quaternion_suite, Subtraction) {
    Quaternion q1(5.0, 6.0, 7.0, 8.0);
    Quaternion q2(1.0, 2.0, 3.0, 4.0);

    Quaternion difference = q1 - q2;
    Quaternion expected(4.0, 4.0, 4.0, 4.0);

    EXPECT_DOUBLE_EQ(difference.GetRealPart(), expected.GetRealPart());
    EXPECT_DOUBLE_EQ(difference.GetI(), expected.GetI());
    EXPECT_DOUBLE_EQ(difference.GetJ(), expected.GetJ());
    EXPECT_DOUBLE_EQ(difference.GetK(), expected.GetK());
}

// AXIS ANGLE

TEST(axisangle_suite, GetRotation) {
    AxisAngle aa(M_PI / 4.0, 0.0, 0.0, 1.0, 0.0); // 45 degrees around z-axis
    
    Matrix4x4 rotationMatrix = aa.GetRotation();
    
    Vector4D vector(1, 0, 0);
    Vector4D rotatedVector = rotationMatrix.MultiplyByVector(vector);

    // Check the result
    EXPECT_NEAR(rotatedVector.x, sqrt(2) / 2, 1e-6);
    EXPECT_NEAR(rotatedVector.y, sqrt(2) / 2, 1e-6);
    EXPECT_NEAR(rotatedVector.z, 0.0, 1e-6);
}

TEST(axisangle_suite, ToQuaternion) {
    // Define axis-angle for 90 degrees around the y-axis
    AxisAngle aa(M_PI / 2, 0.0, 1.0, 0.0, 0.0);

    // Convert to quaternion
    Quaternion q = aa.ToQuaternion();

    // For a 90 degree rotation around the y-axis, the quaternion should be (0, sqrt(2)/2, 0, sqrt(2)/2)
    EXPECT_NEAR(q.GetRealPart(), sqrt(2) / 2, 1e-6);
    EXPECT_NEAR(q.GetI(), 0.0, 1e-6);
    EXPECT_NEAR(q.GetJ(), sqrt(2) / 2, 1e-6);
    EXPECT_NEAR(q.GetK(), 0.0, 1e-6);
}

// POLAR COORDS

bool PolarCoordsEqual(const PolarCoords& a, const PolarCoords& b, double tolerance = 1e-5) {
    return std::fabs(a.GetRadius() - b.GetRadius()) < tolerance &&
        std::fabs(a.GetAngle() - b.GetAngle()) < tolerance;
}

TEST(polarcoords_suite, NormaliseAngle) {
    PolarCoords withinRange(5, M_PI);
    EXPECT_TRUE(PolarCoordsEqual(withinRange.NormaliseAngle(), PolarCoords(5, M_PI)))
        << "NormaliseAngle() failed to normalize an angle within the range.";
    
    PolarCoords overRange(5, 3 * M_PI); // angle greater than 2*PI
    EXPECT_TRUE(PolarCoordsEqual(overRange.NormaliseAngle(), PolarCoords(5, M_PI)))
        << "NormaliseAngle() failed to normalize an angle over the 2*PI range.";

    PolarCoords negativeAngle(5, -M_PI / 2); // negative angle
    EXPECT_TRUE(PolarCoordsEqual(negativeAngle.NormaliseAngle(), PolarCoords(5, 3 * M_PI / 2)))
        << "NormaliseAngle() failed to normalize a negative angle.";
}

TEST(polarcoords_suite, PC2Cartesian) {
    PolarCoords polar(5, M_PI / 4);
    double x, y;
    polar.PC2Cartesian(x, y);

    double expectedX = 5 * std::cos(M_PI / 4);
    double expectedY = 5 * std::sin(M_PI / 4);

    EXPECT_NEAR(x, expectedX, 1e-5) << "PC2Cartesian() X coordinate mismatch: Expected " << expectedX << ", got " << x;
    EXPECT_NEAR(y, expectedY, 1e-5) << "PC2Cartesian() Y coordinate mismatch: Expected " << expectedY << ", got " << y;
}

TEST(polarcoords_suite, PC2Cartesian2) {
    PolarCoords pc(5.0, 1.04719755);
    double x, y;
    double xRes = 2.5;
    double yRes = 4.3301;
    pc.PC2Cartesian(x, y);

    EXPECT_EQ(round(x * 1000.0) / 1000.0, round(2.5 * 1000.0) / 1000.0) << "2.5";
    EXPECT_EQ(round(y * 1000.0) / 1000.0, round(4.3301 * 1000.0) / 1000.0) << "4.3301";
}

TEST(polarcoords_suite, Cartesian2PC) {
    double x = 5.0;
    double y = 5.0;
    PolarCoords polar = PolarCoords::Cartesian2PC(x, y);
    
    double expectedRadius = std::sqrt(x * x + y * y);
    double expectedTheta = std::atan2(y, x);

    EXPECT_NEAR(polar.GetRadius(), expectedRadius, 1e-5)
        << "Cartesian2PC() Radius mismatch: Expected " << expectedRadius << ", got " << polar.GetRadius();
    EXPECT_NEAR(polar.GetAngle(), expectedTheta, 1e-5)
        << "Cartesian2PC() Angle mismatch: Expected " << expectedTheta << ", got " << polar.GetAngle();
}

TEST(quaternion_suite, ToAxisAngle) {
    Quaternion q(0.70710678118, 0.0, 0.0, 0.70710678118);

    AxisAngle axisAngle = q.ToAxisAngle();
    
    double expected_angle = M_PI / 2;
    double expected_x = 0.0;
    double expected_y = 0.0;
    double expected_z = 1.0;

    EXPECT_NEAR(axisAngle.angle, expected_angle, 1e-6);
    EXPECT_NEAR(axisAngle.x, expected_x, 1e-6);
    EXPECT_NEAR(axisAngle.y, expected_y, 1e-6);
    EXPECT_NEAR(axisAngle.z, expected_z, 1e-6);
    EXPECT_NEAR(axisAngle.w, 0.0, 1e-6);
}
