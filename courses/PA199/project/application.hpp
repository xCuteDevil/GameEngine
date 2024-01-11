// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) 2021-2022 Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#pragma once
#include "iapplication.h"
#include "glad/glad.h"
#include "Opengl/Camera/Camera.hpp"
#include "Opengl/Shape/Shape.hpp"
#include "Opengl/Mesh/Mesh.hpp"
#include "Engine/Vector4D/Vector4D.hpp"
#include "Engine/Matrix4x4/Matrix4x4.hpp"
#include "Engine/PolarCoords/PolarCoords.hpp"
#include <tuple>


class Application : public IApplication {
    // ----------------------------------------------------------------------------
    // Variables
    // ----------------------------------------------------------------------------
private:
    bool orthographicProj = false;
    bool left = false, right = false;

    // Brick configuration
    const int bricksPerStory = 12;
    const int numberOfStories = 4;
    const float brickHeight = 1.0f;
    const float brickWidth = 1.5f;
    const int brickDetail = 5;
    const float brickInnerRadius = 2.5f;
    const float brickOuterRadius = brickInnerRadius + brickWidth;


    // Paddle configuration
    const int paddleCount = 3;
    const float paddleWidth = 1.0f;
    const float paddleHeight = 1.0f;
    const float paddleInnerRadius = 11.0f;
    float paddleOuterRadius = paddleInnerRadius + paddleWidth;
    const int paddleDetail = 36;

    // Ground configuration
    const float groundDiameter = 28.0f;

    bool collisionDetected = false;
    double prevRad = 0.0f;

    GLuint vertex_shader;
    GLuint fragment_shader;
    GLuint shader_program;
    
	GLuint texture_vertex_shader;
	GLuint texture_fragment_shader;
    GLuint texture_program;
    
    GLuint vertex_array;
    GLuint vertex_buffer;
    GLuint index_buffer;

    // Texture program
    GLuint texture_vertex_array;
    GLuint texture_vertex_buffer;
    GLuint texture_index_buffer;
    

    
    Camera cam;

    // Control
	bool isPaused = false;
	bool isGameOver = false;
    bool isGameWon = false;
    int playerLives = 3;
	int bricksLeft = bricksPerStory * numberOfStories;

    // Textures
    GLuint textureWin;
    GLuint textureLoss;
    GLuint texturePause;

    bool isBallInGame = false;
    const float ballRadius = 0.5f;
    const float ballStartPosCoefOffset = 1.0f;
    const int ballShapesVectorIndex = numberOfStories * bricksPerStory + 1;
    float ballSpeed = 0.01f; // Should be between 0.01f and 0.02f
	float paddleSpeed = 0.015f;

    float brickCooldownDuration = 150.0f; // time to wait
    float paddleCooldownDuration = 1000.0f; // time to wait


    // ----------------------------------------------------------------------------
    // Variables (Geometry)
    // ----------------------------------------------------------------------------
    std::vector<Shape> shapes;
    std::vector<Shape*> paddles;
    Shape* ball = nullptr;
    std::vector<Shape*> groundLevelBricks;

    // ----------------------------------------------------------------------------
    // Constructors & Destructors
    // ----------------------------------------------------------------------------
public:
    Application(int initial_width, int initial_height, std::vector<std::string> arguments = {});

    /** Destroys the {@link Application} and releases the allocated resources. */
    virtual ~Application();

    // ----------------------------------------------------------------------------
    // Methods
    // ----------------------------------------------------------------------------

    /** @copydoc IApplication::update */
    void update(float delta) override;

    /** @copydoc IApplication::render */
    void render() override;

    /** @copydoc IApplication::render_ui */
    void render_ui() override;

    /** @copydoc IApplication::on_resize */
    void on_resize(int width, int height) override;

    /** @copydoc IApplication::on_mouse_move */
    void on_mouse_move(double x, double y) override;

    /** @copydoc IApplication::on_mouse_button */
    void on_mouse_button(int button, int action, int mods) override;

    /** @copydoc IApplication::on_key_pressed */
    void on_key_pressed(int key, int scancode, int action, int mods) override;

    void startGame();

    void BallPhysicsUpdate(float delta, Shape& shape);

    void BroadPhaseDetection(Shape& ball);

    void CollisionWithBricks(Shape& ball);
    bool ProcessBrickCollision(Shape& ball, Shape* brick, float ballAngle, float distanceFromCenter, int colId);
    Vector4D CalculateCollisionNormal(const Vector4D& ballPosition, const Vector4D& paddlePosition);
    Vector4D Reflect(const Vector4D& direction, const Vector4D& normal);
    void ReflectBall(Shape& ball, const Vector4D& normal, float speed);

    void CollisionWithPaddles(Shape& ball);
    bool ProcessPaddleCollision(Shape& ball, Shape* paddle, float ballAngle, float distanceFromCenter);

    Vector4D GetTransformedVertex(Shape* shape, const Vertex& vertex);
    bool IsBallWithinObstacleRange(float ballAngle, float obstacleStartAngle, float obstacleEndAngle);
    void SetDirection(Shape& shape, Vector4D newDirection, float speed);

    std::tuple<float, float, float> NormalizeCollisionAngles(float ballAngle, float obstacleStartAngle, float obstacleEndAngle);
	Vector4D ClosestPointOnTheLine(Vector4D lineStart, Vector4D lineEnd, Vector4D point);
    void CheckWinLossConditions();
	void ResetCollisionCooldowns(float delta);
    
    std::vector<Vector4D> colors = {
        Vector4D(1, 1, 0, 1), // Yellow
        Vector4D(0, 1, 0, 1), // Green
        Vector4D(0, 0, 1, 1), // Blue
        Vector4D(1, 0, 0, 1), // Red
        Vector4D(0, 0, 0, 1)  // Black
    };
    void renderTexture(GLuint texture);
};

