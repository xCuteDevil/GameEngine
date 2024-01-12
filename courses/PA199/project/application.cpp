// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
// 
// Copyright (c) 2021-2022 Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#include "application.hpp"
#include "glad/glad.h"
#include "lodepng.h"
#include <filesystem>
#include <iostream>

static GLuint load_shader(std::filesystem::path const& path, GLenum const shader_type)
{
    std::filesystem::path const current = std::filesystem::current_path();
    GLuint const shader = glCreateShader(shader_type);
    assert(glGetError() == 0U && shader != 0);
    std::ifstream ifs(path);
    assert(ifs.is_open());
    std::string const str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    char const* code = str.c_str();
    glShaderSource(shader, 1, &code, nullptr);
    assert(glGetError() == 0U);
    glCompileShader(shader);
    assert(glGetError() == 0U);
    return shader;
}

static GLuint load_texture(std::filesystem::path const& path)
{
    std::vector<unsigned char> texels;
    unsigned int width, height;
    unsigned int error_code = lodepng::decode(texels, width, height, path.string(), LCT_RGBA);
    assert(error_code == 0);

    //flip the image vertically
    for (unsigned int lo = 0, hi = height - 1; lo < hi; ++lo, --hi)
        for (unsigned int*  lo_ptr = (unsigned int*)texels.data() + lo * width,
                         *  lo_end = lo_ptr + width,
                         *  hi_ptr = (unsigned int*)texels.data() + hi * width;
                lo_ptr != lo_end; ++lo_ptr, ++hi_ptr)
            std::swap(*lo_ptr, *hi_ptr);

    GLuint texture;
	glGenTextures(1, &texture);
    assert(glGetError() == 0U);
	glBindTexture(GL_TEXTURE_2D, texture);
    assert(glGetError() == 0U);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, (void*)texels.data());
    assert(glGetError() == 0U);
    glGenerateMipmap(GL_TEXTURE_2D);
    assert(glGetError() == 0U);
    return texture;
}

// ----------------------------------------------------------------------------
// Constructors & Destructors
// ----------------------------------------------------------------------------

Application::Application(int initial_width, int initial_height, std::vector<std::string> arguments)
    : IApplication(initial_width, initial_height, arguments)
    , vertex_shader(load_shader(lecture_folder_path / "data" / "shaders" / "sample.vert", GL_VERTEX_SHADER))
    , fragment_shader(load_shader(lecture_folder_path / "data" / "shaders" / "sample.frag", GL_FRAGMENT_SHADER))
    , texture_vertex_shader(load_shader(lecture_folder_path / "data" / "shaders" / "texture_vertex_shader.vert", GL_VERTEX_SHADER))
    , texture_fragment_shader(load_shader(lecture_folder_path / "data" / "shaders" / "texture_fragment_shader.frag", GL_FRAGMENT_SHADER))
    , shader_program([](GLuint const  vertex_shader, GLuint const  fragment_shader) -> GLuint {
            GLuint const  shader_program = glCreateProgram();
            assert(glGetError() == 0U && shader_program != 0);
            glAttachShader(shader_program,vertex_shader);
            assert(glGetError() == 0U);
            glAttachShader(shader_program,fragment_shader);
            assert(glGetError() == 0U);
            glLinkProgram(shader_program);
            assert(glGetError() == 0U);
            glDetachShader(shader_program, vertex_shader);
            assert(glGetError() == 0U);
            glDetachShader(shader_program, fragment_shader);
            assert(glGetError() == 0U);
            return shader_program;
    }(vertex_shader, fragment_shader))
    , texture_program([](GLuint const  texture_vertex_shader, GLuint const  texture_fragment_shader) -> GLuint {
            GLuint const  texture_program = glCreateProgram();
            assert(glGetError() == 0U && texture_program != 0);
            glAttachShader(texture_program, texture_vertex_shader);
            assert(glGetError() == 0U);
            glAttachShader(texture_program, texture_fragment_shader);
            assert(glGetError() == 0U);
            glLinkProgram(texture_program);
            assert(glGetError() == 0U);
            glDetachShader(texture_program, texture_vertex_shader);
            assert(glGetError() == 0U);
            glDetachShader(texture_program, texture_fragment_shader);
            assert(glGetError() == 0U);
            return texture_program;
    }(texture_vertex_shader, texture_fragment_shader))
    , textureWin(load_texture(lecture_folder_path / "data" / "textures" / "you_win.png"))
    , textureLoss(load_texture(lecture_folder_path / "data" / "textures" / "game_over.png"))
    , texturePause(load_texture(lecture_folder_path / "data" / "textures" / "pause.png"))
    , texture_vertex_array([]() -> GLuint {
            GLuint vertex_array;
            glGenVertexArrays(1, &vertex_array);
            glBindVertexArray(vertex_array);
            return vertex_array;
    }())
    , texture_vertex_buffer([]() -> GLuint {
            struct Vertex { float x, y, z, u, v; };
            std::vector<Vertex> const vertices{
                { -0.5f, -0.125f, 0.0f,       0.0f, 0.0f },
                {  0.5f, -0.125f, 0.0f,       1.0f, 0.0f },
                {  0.5f,  0.125f, 0.0f,       1.0f, 1.0f },
                { -0.5f,  0.125f, 0.0f,       0.0f, 1.0f },
            };
            GLuint vertex_buffer;
            glGenBuffers(1, &vertex_buffer);
            assert(glGetError() == 0U);
            glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
            assert(glGetError() == 0U);
            glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), &vertices.front(), GL_STATIC_DRAW);
            assert(glGetError() == 0U);
            glEnableVertexAttribArray(0);
            assert(glGetError() == 0U);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), nullptr);
            assert(glGetError() == 0U);
            glEnableVertexAttribArray(1);
            assert(glGetError() == 0U);
            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3 * sizeof(float)));
            assert(glGetError() == 0U);
            return vertex_buffer;
    }())
    , texture_index_buffer([]() -> GLuint {
            std::vector<GLuint> const indices{ 0U, 1U, 2U, 0U, 2U, 3U };
            GLuint index_buffer;
            glGenBuffers(1, &index_buffer);
            assert(glGetError() == 0U);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
            assert(glGetError() == 0U);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), &indices.front(), GL_STATIC_DRAW);
            assert(glGetError() == 0U);
            return index_buffer;
    }())
{
    glViewport(0, 0, width, height);
    glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_CULL_FACE);
    
    startGame();
}

Application::~Application()
{
    glDeleteVertexArrays(1U, &vertex_array);
    glDeleteTextures(1, &textureWin);
    glDeleteTextures(1, &textureLoss);
    glDeleteTextures(1, &texturePause);

    glDeleteProgram(shader_program);
    glDeleteShader(fragment_shader);
    glDeleteShader(vertex_shader);
    
    glDeleteProgram(texture_program);
    glDeleteShader(texture_fragment_shader);
    glDeleteShader(texture_vertex_shader);
    
    
    for (int i = 0; i < shapes.size(); i++) {
        GLuint array = shapes[i].GetVertexArray();
        glDeleteVertexArrays(1U, &array);
    }
}

// ----------------------------------------------------------------------------
// Methods
// ----------------------------------------------------------------------------

bool Application::CheckGameWon() {
    for (int i = 0; i < bricksPerStory; i++) {
        if (!groundLevelBricks[i]->IsColumnDestroyed()) {
            return false; // Game continues if any column is not fully destroyed
        }
    }
    return true; // All columns are destroyed, player wins
}

void Application::CheckWinLossConditions() {
    // Win Condition
	if (CheckGameWon()) {
		isGameWon = true;
	}
    // Loss Condition
    else if (playerLives <= 0)
    {
        isGameOver = true;
    }
    else
    {
        return;
    }
    // Stop Rendering ball
    shapes[ballShapesVectorIndex].SetDestroyed(true);
    isBallInGame = false;
}

void Application::MovePaddle(Shape* paddle, bool moveLeft, float delta) {
    auto toRadians = [](float degrees) { return degrees * 3.14159265358979323846f / 180.0f; };
    Matrix4x4 model = paddle->GetModelMatrix();
    float rotationDirection = moveLeft ? 1.0f : -1.0f; // Positive for left, negative for right

    // Assume paddleSpeed is in degrees per second
    float rotationAngle = paddleSpeed * delta; // Rotation angle per frame
    model = model * Matrix4x4(1.0).Rotate(toRadians(rotationAngle) * rotationDirection, Vector4D(0, 1, 0, 0));
    paddle->SetModelMatrix(model);
}

void Application::update(float delta) {
    
	// Win/Lose Condition
    CheckWinLossConditions();
	if (isGameWon || isGameOver) return;

    // Movement left and right
    if (left && !isPaused) {
        for (Shape* p : paddles) {
            MovePaddle(p, true, delta);
        }
    }
    else if (right && !isPaused) {
        for (Shape* p : paddles) {
            MovePaddle(p, false, delta);
        }
    }

    // Waiting to cast the ball
    if (!isBallInGame)
    {
        Shape* p = paddles[0];
        Vector4D paddlePos = p->CalculatePosition();
		paddlePos = paddlePos - paddlePos.UnitVector() * ballRadius * ballStartPosCoefOffset;
        shapes[ballShapesVectorIndex].SetModelMatrix(Matrix4x4(1.0).Translate(paddlePos.x, ballRadius, paddlePos.z));
        shapes[ballShapesVectorIndex].isDirectionSet = false;
    }
    else
    {
        // Space press, set initial position and direction
        if (!shapes[ballShapesVectorIndex].isDirectionSet)
        {
            Vector4D currentPos = shapes[ballShapesVectorIndex].CalculatePosition();
            shapes[ballShapesVectorIndex].SetPosition(currentPos);
            Vector4D dir = Vector4D(-currentPos.x+1.0f, 0, -currentPos.z-1.0f);
            SetDirection(shapes[ballShapesVectorIndex], dir.UnitVector(), ballSpeed);
            shapes[ballShapesVectorIndex].isDirectionSet = true;
        }
        
        BallPhysicsUpdate(delta, shapes[ballShapesVectorIndex]);

        PowerUpsPhysicsUpdate(delta);

        BroadPhaseDetection(shapes[ballShapesVectorIndex]);

        for (int i = 0; i < powerUps.size(); i++) {
            if (!powerUps[i].IsDestroyed())
            {
                bool didCatch = BroadPhasePowerUpDetection(powerUps[i]);
                if (didCatch) {
                    
                    
                    // GREEN - ENLARGE PADDLES
                    if (Vector4D::Equals(powerUps[i].GetColour(), Vector4D(0, 1, 0, 1), 3))
                    {
                        if (paddleCount * (paddleLength * 1.05f) < 360.0f)
                        {
                            paddleLength = paddleLength * 1.05f;
                            RegeneratePaddles();
                        }
                    }
                    // RED - LESSEN PADDLES
                    else if (Vector4D::Equals(powerUps[i].GetColour(), Vector4D(1, 0, 0, 1), 3))
                    {
                        if (paddleCount * (paddleLength * 0.95f) > 100.0f)
                        {
                            paddleLength = paddleLength * 0.95f;
                            RegeneratePaddles();
                        }
                    }
                    // YELLOW - MAKE BALL FASTER
                    else if (Vector4D::Equals(powerUps[i].GetColour(), Vector4D(1, 1, 0, 1), 3))
                    {
                        if (speedModifier * 1.0f < 1.5f)
                        {
                            speedModifier = speedModifier * 1.05f;
                        }
                    }
                    // BLUE - MAKE BALL SLOWER
                    else if (Vector4D::Equals(powerUps[i].GetColour(), Vector4D(0, 0, 1, 1), 3))
                    {
                        if (speedModifier * 0.95f > 0.5f)
                        {
                            speedModifier = speedModifier * 0.95f;
                        }
                    }
                    powerUps.erase(powerUps.begin() + i);
                }
            }
        }
        
        // Handles the colission cooldown period for each brick/paddle after a collision.
        ResetCollisionCooldowns(delta);
    }
}

void Application::ResetCollisionCooldowns(float delta)
{
    int i = 0;
    while (i < bricksPerStory) {
        groundLevelBricks[i]->Update(delta);
        i++;
    }
    i = 0;
    while (i < paddles.size())
    {
        paddles[i]->Update(delta);
        i++;
    }
}
bool Application::BroadPhasePowerUpDetection(Shape& powerUp) {
    PolarCoords pUpCoords = PolarCoords::Cartesian2PC(powerUp.position.x, powerUp.position.z);
    float pUpAngle = pUpCoords.GetAngle();
    float distanceFromCenter = pUpCoords.GetRadius();
    float pUpRadiusPlus = distanceFromCenter + powerUpRadius;
    float pUpRadiusMinus = distanceFromCenter - powerUpRadius;
        
    if (pUpRadiusMinus < paddleOuterRadius && pUpRadiusPlus > paddleInnerRadius)
    {
        Vector4D ballPos = powerUp.position;
        PolarCoords ballPolarCoords = PolarCoords::Cartesian2PC(ballPos.x, ballPos.z);

        for (int p = 0; p < paddles.size(); p++)
        {
            if (!paddles[p]->isOnCooldown)
            {
                bool collided = ProcessPaddleCollision(powerUp, paddles[p], ballPolarCoords.GetAngle(), ballPolarCoords.GetRadius());
                // Colission found, stop 
                if (collided) {
                    powerUp.SetDestroyed(true);
                    return true;
                }
            }
        }
    }
    return false;
}

void Application::BroadPhaseDetection(Shape& ball) {

    PolarCoords ballPolarCoords = PolarCoords::Cartesian2PC(ball.position.x, ball.position.z);
    float ballAngle = ballPolarCoords.GetAngle();
    float distanceFromCenter = ballPolarCoords.GetRadius();
    float ballRadiusPlus = distanceFromCenter + ballRadius;
    float ballRadiusMinus = distanceFromCenter - ballRadius;

    if (distanceFromCenter > (groundDiameter * 0.5f))
    {
        isBallInGame = false;
        playerLives--;
        // remove power-ups
        paddleLength = initialPaddleLength;
        speedModifier = initialSpeedModifierValue;
        RegeneratePaddles();
        
        powerUps.clear();
    }
    // Potential Collision with bricks
    else if (ballRadiusMinus <= (brickOuterRadius) && ballRadiusPlus >= (brickInnerRadius))
    {
        CollisionWithBricks(ball);

    }
    // Potential Collision with paddles
    else if (ballRadiusMinus < paddleOuterRadius && ballRadiusPlus > paddleInnerRadius)
    {
        CollisionWithPaddles(ball);
    }
}

void Application::GeneratePowerUp(Vector4D colour, Vector4D position) {
    // Create the power-up
    Matrix4x4 model = Matrix4x4(1.0f);
    Sphere b(position, powerUpRadius, powerUpStacks, powerUpSectors, colour);
    b.SetModelMatrix(model);
    b.SetPosition(position);

    // Calculate and set initial velocity
    Vector4D direction = position;
    direction = direction.UnitVector();
    b.velocity = direction * powerUpSpeed;
    b.SetDestroyed(true); // stands for do not render
    powerUps.push_back(b);
}

Shape* Application::GetBrickToBeDestroyed(int columnID) {
    Shape* b = groundLevelBricks[columnID];
    while (b != nullptr) {
        if (!b->IsDestroyed())
        {
            return b;
        }
        b = b->next;
    }
}

bool Application::ProcessBrickCollision(Shape& ball, Shape* brick, float ballAngle, float distanceFromCenter, int colId)
{
    int nextBrickId = (colId + 1) % bricksPerStory;
    int prevBrickId = colId - 1;
    if (prevBrickId < 0) {
        prevBrickId = bricksPerStory - 1;
    }

    Vector4D brickPos = brick->CalculatePosition();
    PolarCoords brickPolarCoords = PolarCoords::Cartesian2PC(brickPos.x, brickPos.z);

    float brickAngle = brickPolarCoords.GetAngle();
    std::vector<Vertex> verticesVector = brick->GetVertices();

    Vector4D startVertex = GetTransformedVertex(brick, verticesVector.front()); // outer radius
    Vector4D endVertex = GetTransformedVertex(brick, verticesVector.back()); // inner radius

    PolarCoords startPC = PolarCoords::Cartesian2PC(startVertex.x, startVertex.z);
    PolarCoords endPC = PolarCoords::Cartesian2PC(endVertex.x, endVertex.z);

    float brickStartAngle = startPC.GetAngle();
    float brickEndAngle = endPC.GetAngle();
    float ballAngleNormalised, brickStartNormalised, brickEndNormalised;

    Vector4D ballPosition = ball.position;
    Vector4D ballDirection = ball.velocity.UnitVector();

    Vector4D newDirection;
    Vector4D ballPosXY(ball.position.x, 0, ball.position.z, 0);
    Vector4D center(0, 0, 0, 1);

    float BallRadiusAngle = atan(ballRadius / distanceFromCenter);

    // Brick front/back collision
    std::tie(ballAngleNormalised, brickStartNormalised, brickEndNormalised) = NormalizeCollisionAngles(ballAngle, brickStartAngle, brickEndAngle);
    if (IsBallWithinObstacleRange(ballAngleNormalised, brickStartNormalised, brickEndNormalised) && !brick->IsColumnDestroyed())
    {
        Vector4D collisionNormal = CalculateCollisionNormal(ballPosition, Vector4D(brickPos.x, ballPosition.y, brickPos.z));
        ReflectBall(ball, collisionNormal, ballSpeed);
        
        Shape *b = GetBrickToBeDestroyed(colId);
        brick->DestroyBrick();
        HandlePowerUpGeneration(b, brick->CalculatePosition());
        return true;
    }

    // Brick start-side collision
    std::tie(ballAngleNormalised, brickStartNormalised, brickEndNormalised) = NormalizeCollisionAngles(ballAngle, brickStartAngle - BallRadiusAngle, brickStartAngle);
    if (IsBallWithinObstacleRange(ballAngleNormalised, brickStartNormalised, brickEndNormalised) && !brick->IsColumnDestroyed() && !(!groundLevelBricks[prevBrickId]->IsColumnDestroyed()))
    {

        Vector4D outerBrickRadiusStartVertex = Vector4D(startVertex.x, 0, startVertex.z);
        Vector4D directionToCenter = center - outerBrickRadiusStartVertex;
        directionToCenter = directionToCenter.UnitVector();
        Vector4D innerBrickRadiusStartVertex = outerBrickRadiusStartVertex + directionToCenter * brickWidth;
        Vector4D closestPointOnTheBrickSide = ClosestPointOnTheLine(innerBrickRadiusStartVertex, outerBrickRadiusStartVertex, ballPosXY);
        newDirection = (ballPosXY - closestPointOnTheBrickSide).UnitVector();
        if (sqrt(pow(ballPosXY.x - closestPointOnTheBrickSide.x, 2) + pow(ballPosXY.z - closestPointOnTheBrickSide.z, 2)) > ballRadius) {
            ;
            //return false;
        }
        /*if (Vector4D::Equals(innerBrickRadiusStartVertex, closestPointOnTheBrickSide, 3)) {
            // Corner collision
            return false;
        }*/

        SetDirection(ball, newDirection, ballSpeed);

        Shape* b = GetBrickToBeDestroyed(colId);
        brick->DestroyBrick();
        HandlePowerUpGeneration(b, brick->CalculatePosition());
        
        return true;
    }

    // Brick end-side collision
    std::tie(ballAngleNormalised, brickStartNormalised, brickEndNormalised) = NormalizeCollisionAngles(ballAngle, brickEndAngle, brickEndAngle + BallRadiusAngle);
    if (IsBallWithinObstacleRange(ballAngleNormalised, brickStartNormalised, brickEndNormalised) && !brick->IsColumnDestroyed() && !(!groundLevelBricks[nextBrickId]->IsColumnDestroyed()))
    {
        Vector4D outerBrickRadiusStartVertex = Vector4D(startVertex.x, 0, startVertex.z);
        Vector4D directionToCenter = center - outerBrickRadiusStartVertex;
        directionToCenter = directionToCenter.UnitVector();
        Vector4D innerBrickRadiusStartVertex = outerBrickRadiusStartVertex + directionToCenter * brickWidth;
        Vector4D closestPointOnTheBrickSide = ClosestPointOnTheLine(innerBrickRadiusStartVertex, outerBrickRadiusStartVertex, ballPosXY);
        newDirection = (ballPosXY - closestPointOnTheBrickSide).UnitVector();
        if (sqrt(pow(ballPosXY.x - closestPointOnTheBrickSide.x, 2) + pow(ballPosXY.z - closestPointOnTheBrickSide.z, 2)) > ballRadius) {
            ;
            //return false;
        }
        /*if (Vector4D::Equals(innerBrickRadiusStartVertex, closestPointOnTheBrickSide, 3)) {
            // Corner collision
            return false;
        }*/

        SetDirection(ball, newDirection, ballSpeed);

        Shape* b = GetBrickToBeDestroyed(colId);
        brick->DestroyBrick();
        HandlePowerUpGeneration(b, brick->CalculatePosition());
        
        return true;
    }

    // Brick front/back collision
    /*std::tie(ballAngleNormalised, brickStartNormalised, brickEndNormalised) = NormalizeCollisionAngles(ballAngle, brickStartAngle, brickEndAngle);
    if (IsBallWithinObstacleRange(ballAngleNormalised, brickStartNormalised, brickEndNormalised) && !brick->IsColumnDestroyed())
    {
        Vector4D collisionNormal = CalculateCollisionNormal(ballPosition, Vector4D(brickPos.x, ballPosition.y, brickPos.z));
        ReflectBall(ball, collisionNormal, ballSpeed);
        Shape* b = GetBrickToBeDestroyed(colId);
        brick->DestroyBrick();
        HandlePowerUpGeneration(b, brick->CalculatePosition());
    }*/
    return false;
}

void Application::HandlePowerUpGeneration(Shape *brick, Vector4D position) {
    if (Vector4D::Equals(brick->GetColour(), colors[3], 3)) {
        colorDestructionCount[colors[3]]++;
        if (colorDestructionCount[colors[3]] == 3) {
            GeneratePowerUp(colors[3], position);
            colorDestructionCount[colors[3]] = 0;
        }
    }
    // GREEN
    else if (Vector4D::Equals(brick->GetColour(), colors[1], 3)) {
        colorDestructionCount[colors[1]]++;
        if (colorDestructionCount[colors[1]] == 3) {
            GeneratePowerUp(colors[1], position);
            colorDestructionCount[colors[1]] = 0;
        }
    }
    // BLUE
    else if (Vector4D::Equals(brick->GetColour(), colors[2], 3)) {
        colorDestructionCount[colors[2]]++;
        if (colorDestructionCount[colors[2]] == 3) {
            GeneratePowerUp(colors[2], position);
            colorDestructionCount[colors[2]] = 0;
        }
    }
    // YELLOW
    else if (Vector4D::Equals(brick->GetColour(), colors[0], 3)) {
        colorDestructionCount[colors[0]]++;
        if (colorDestructionCount[colors[0]] == 3) {
            GeneratePowerUp(colors[0], position);
            colorDestructionCount[colors[0]] = 0;
        }
    }
}

Vector4D Application::CalculateCollisionNormal(const Vector4D& ballPosition, const Vector4D& obstaclePosition)
{
    Vector4D normal = ballPosition.UnitVector();
    PolarCoords ballPolarCoords = PolarCoords::Cartesian2PC(ballPosition.x, ballPosition.z);
    PolarCoords brickPolarCoords = PolarCoords::Cartesian2PC(obstaclePosition.x, obstaclePosition.z);

    // Front or back collision
    if ((ballPolarCoords.GetRadius()) < brickPolarCoords.GetRadius())
    {
        normal = normal.OppositeVector();
    }
    return normal;
}

void Application::ReflectBall(Shape& ball, const Vector4D& normal, float speed)
{
    Vector4D direction = ball.velocity.UnitVector();
    Vector4D newDirection = Reflect(direction, normal);
    newDirection.y = 0;
    SetDirection(ball, newDirection, speed);
}

Vector4D Application::Reflect(const Vector4D& direction, const Vector4D& normal)
{
    return direction - normal * 2 * (direction.DotProduct(normal));
}

std::tuple<float, float, float> Application::NormalizeCollisionAngles(float ballAngle, float paddleStartAngle, float paddleEndAngle) {
    if (paddleStartAngle > paddleEndAngle) {
        paddleEndAngle += 2 * M_PI;
        if (ballAngle < paddleStartAngle) {
            ballAngle += 2 * M_PI;
        }
    }
    return std::make_tuple(ballAngle, paddleStartAngle, paddleEndAngle);
}

void Application::CollisionWithPaddles(Shape& ball)
{
	Vector4D ballPos = ball.position;
    PolarCoords ballPolarCoords = PolarCoords::Cartesian2PC(ballPos.x, ballPos.z);
    
	for (int p = 0; p < paddles.size(); p++)
	{
        if (!paddles[p]->isOnCooldown)
        {
            bool collided = ProcessPaddleCollision(ball, paddles[p], ballPolarCoords.GetAngle(), ballPolarCoords.GetRadius());
            // Colission found, stop 
            if (collided) {
                paddles[p]->StartCooldown(paddleCooldownDuration);
                return;
            }
        }
	}
}

bool Application::ProcessPaddleCollision(Shape& ball, Shape* paddle, float ballAngle, float ballDistanceFromCenter)
{
    Vector4D paddlePos = paddle->CalculatePosition();
    PolarCoords paddlePolarCoords = PolarCoords::Cartesian2PC(paddlePos.x, paddlePos.z);
    
    // Paddle's angular velocity if the paddle is moving
    const float angularVelocity = (right || left) ? 0.1f : 0.0f;

    // Determine the direction of the paddle's movement (clockwise or counterclockwise)
    Vector4D tangentialDirection = Vector4D(sin(paddlePolarCoords.GetAngle()), 0, cos(paddlePolarCoords.GetAngle()), 0) * (right ? -1.0f : 1.0f);
    Vector4D paddleVelocity = tangentialDirection * angularVelocity;

    float paddleAngle = paddlePolarCoords.GetAngle();
    std::vector<Vertex> verticesVector = paddle->GetVertices();
    
    Vector4D startVertex = GetTransformedVertex(paddle, verticesVector.front());
    Vector4D endVertex = GetTransformedVertex(paddle, verticesVector.back());

    PolarCoords startPC = PolarCoords::Cartesian2PC(startVertex.x, startVertex.z);
    PolarCoords endPC = PolarCoords::Cartesian2PC(endVertex.x, endVertex.z);

    Vector4D newDirection;
    float BallRadiusAngle = atan(ballRadius / ballDistanceFromCenter);
    Vector4D ballPosXY(ball.position.x, 0, ball.position.z, 0);
    Vector4D center(0, 0, 0, 1);

    float paddleStartAngle = startPC.GetAngle(); // on the outer radius
    float paddleEndAngle = endPC.GetAngle(); // on the inner radius
    float ballAngleNormalised, paddleStartNormalised, paddleEndNormalised;

    // Front-side paddle collision
    std::tie(ballAngleNormalised, paddleStartNormalised, paddleEndNormalised) = NormalizeCollisionAngles(ballAngle, paddleStartAngle, paddleEndAngle);
    if (IsBallWithinObstacleRange(ballAngleNormalised, paddleStartNormalised, paddleEndNormalised))
    {
        // Reflect the ball's direction based on the collision
        Vector4D collisionNormal = CalculateCollisionNormal(Vector4D(ball.position.x, 0, ball.position.z), Vector4D(paddlePos.x, 0, paddlePos.z));
        newDirection = Reflect(ball.velocity, collisionNormal);

        // Apply the slice effect based on the paddle's velocity
        float frictionCoefficient = 0.06f;
        Vector4D sliceForce = paddleVelocity * frictionCoefficient;
        newDirection = newDirection + sliceForce;
        newDirection = newDirection.UnitVector();

		// If the newDirection is too close to the paddle's direction, adjust it
        newDirection = AdjustDeflectionDirection(ball, newDirection, 0.30f, 0.1f);
        
        SetDirection(ball, newDirection.UnitVector(), ballSpeed);
        return true;
    }
    // Start-side paddle collision
    std::tie(ballAngleNormalised, paddleStartNormalised, paddleEndNormalised) = NormalizeCollisionAngles(ballAngle, paddleStartAngle - BallRadiusAngle, paddleEndAngle);
    if (IsBallWithinObstacleRange(ballAngleNormalised, paddleStartNormalised, paddleEndNormalised))
    {
        Vector4D outerPaddleRadiusStartVertex = Vector4D(startVertex.x, 0, startVertex.z);
        Vector4D directionToCenter = center - outerPaddleRadiusStartVertex;
        directionToCenter = directionToCenter.UnitVector(); 
        Vector4D innerPaddleRadiusStartVertex = outerPaddleRadiusStartVertex + directionToCenter * paddleWidth;
        Vector4D closestPointOnThePaddleSide = ClosestPointOnTheLine(innerPaddleRadiusStartVertex, outerPaddleRadiusStartVertex, ballPosXY);
        newDirection = (ballPosXY - closestPointOnThePaddleSide).UnitVector();
        /*if (Vector4D::Equals(innerPaddleRadiusStartVertex, closestPointOnThePaddleSide, 1)) {
            // Corner collision
            SetDirection(ball, newDirection, ballSpeed);
        }*/
        // So that ball is faster than paddle if it goes along with it 
        SetDirection(ball, newDirection, ballSpeed*1.2f);
        return true;
    }

    // End-side paddle collision
    std::tie(ballAngleNormalised, paddleStartNormalised, paddleEndNormalised) = NormalizeCollisionAngles(ballAngle, paddleEndAngle, paddleEndAngle + BallRadiusAngle);
    if (IsBallWithinObstacleRange(ballAngleNormalised, paddleStartNormalised, paddleEndNormalised))
    {
        Vector4D innerPaddleRadiusStartVertex = Vector4D(endVertex.x, 0, endVertex.z);
        Vector4D directionToCenter = center - innerPaddleRadiusStartVertex;
        directionToCenter = directionToCenter.UnitVector();
        Vector4D outerPaddleRadiusStartVertex = innerPaddleRadiusStartVertex - directionToCenter * paddleWidth;
        Vector4D closestPointOnThePaddleSide = ClosestPointOnTheLine(outerPaddleRadiusStartVertex, innerPaddleRadiusStartVertex, ballPosXY);
        newDirection = (ballPosXY - closestPointOnThePaddleSide).UnitVector();
        /*if (Vector4D::Equals(innerPaddleRadiusStartVertex, closestPointOnThePaddleSide, 1)) {
            // Corner collision
            SetDirection(ball, newDirection, ballSpeed);
        }*/
        
        // So that ball is faster than paddle if it goes along with it 
        SetDirection(ball, newDirection, ballSpeed*1.2f);
        return true;
    }
    return false;
}

Vector4D Application::AdjustDeflectionDirection(Shape& ball, Vector4D realDirection, float similarityThreshold, float blendFactor) {

    // Calculate the dot product to measure similarity
    float similarity = realDirection.UnitVector().DotProduct(ball.velocity.UnitVector());

    if (abs(similarity) < similarityThreshold) {
        Vector4D ballPositionNormalized = Vector4D(ball.position.x, 0, ball.position.z);
        ballPositionNormalized = Vector4D(-ball.position.x, 0, -ball.position.z).UnitVector();
        
        realDirection = realDirection * (1.0f - blendFactor) + ballPositionNormalized * blendFactor;
    }
    return realDirection.UnitVector();
}

Vector4D Application::ClosestPointOnTheLine(Vector4D lineStart, Vector4D lineEnd, Vector4D point)
{    
    Vector4D lineVec = lineEnd - lineStart;
    Vector4D pointVec = point - lineStart;
    double lineLengthSquared = lineVec.DotProduct(lineVec);
    // Project pointVec onto lineVec
    double t = pointVec.DotProduct(lineVec) / lineLengthSquared;
    t = std::max(0.0, std::min(1.0, t));
    Vector4D closestPoint = lineStart + lineVec * t;

    return closestPoint;
}

bool Application::IsBallWithinObstacleRange(float ballAngle, float obstacleStartAngle, float obstacleEndAngle)
{
	return ballAngle >= obstacleStartAngle && ballAngle <= obstacleEndAngle;
}

Vector4D Application::GetTransformedVertex(Shape* shape, const Vertex& vertex) {
    return shape->GetModelMatrix() * Vector4D(vertex.x, vertex.y, vertex.z, 1);
}

void Application::CollisionWithBricks(Shape& ball)
{
    Vector4D ballPos = ball.position;
    float ballDistanceFromCenter = sqrt(ballPos.x * ballPos.x + ballPos.z * ballPos.z);
	PolarCoords ballPolarCoords = PolarCoords::Cartesian2PC(ballPos.x, ballPos.z);

    for (int b = 0; b < groundLevelBricks.size(); b++)
    {
        if (!groundLevelBricks[b]->IsColumnDestroyed() && !groundLevelBricks[b]->isOnCooldown)
        {
            bool collided = ProcessBrickCollision(ball, groundLevelBricks[b], ballPolarCoords.GetAngle(), ballDistanceFromCenter, b);
            // Colission found, stop 
            if (collided) {
                groundLevelBricks[b]->StartCooldown(brickCooldownDuration);
                return;
            }
        }
    }
}

void Application::PowerUpsPhysicsUpdate(float delta) {
    if (isPaused) {
        return;
    }

    for (int i = 0; i < powerUps.size(); i++) {
        float radius = PolarCoords::Cartesian2PC(powerUps[i].position.x, powerUps[i].position.z).GetRadius();
        if (radius < (groundDiameter*0.5f) && radius > brickInnerRadius) {
            powerUps[i].SetDestroyed(false); // stands for do render
        }
        else {
            powerUps[i].SetDestroyed(true); // stands for do not render
        }
        // Update position based on current velocity
        Vector4D newPosition = powerUps[i].position + (powerUps[i].velocity * delta);
        powerUps[i].SetPosition(newPosition);

        // Update model matrix
        Matrix4x4 model = powerUps[i].GetModelMatrix();
        model = Matrix4x4(1.0f).Translate(newPosition.x, ballRadius, newPosition.z);
        powerUps[i].SetModelMatrix(model);
    }
}

void Application::BallPhysicsUpdate(float delta, Shape& shape)
{   
    if (isPaused) {
        return;
    }
    
    // Update position
    Vector4D newPosition = shape.position + (shape.velocity * speedModifier * delta);
	shape.SetPosition(newPosition);

    // Update model matrix
    Matrix4x4 model = shape.GetModelMatrix();
    model = Matrix4x4(1.0f).Translate(newPosition.x, ballRadius, newPosition.z);
    shape.SetModelMatrix(model);
}

void Application::SetDirection(Shape& shape, Vector4D newDirection, float speed) {
    Vector4D unitDirection = newDirection.UnitVector();
    shape.velocity = unitDirection * speed;
}

void Application::startGame() {
    // Reset variables
    ball = nullptr;
	paddles = std::vector<Shape*>();
	groundLevelBricks = std::vector<Shape*>();
    powerUps = std::vector<Shape>();
    shapes.clear();
    shapes.reserve(bricksPerStory*numberOfStories + paddleCount + 2);
    paddleLength = initialPaddleLength;
    speedModifier = initialSpeedModifierValue;
        
    orthographicProj = false;
    float aspectRatio = static_cast<float>(width) / static_cast<float>(height);
    cam.setprojectionMatrixToPerspective(45.0f, aspectRatio, 0.1f, 100.0f);
    
    Matrix4x4 model = Matrix4x4(1.0f);

    /*** SHAPES ***/

    // Bricks
    
    // Conversion helper
    auto toRadians = [](float degrees) { return degrees * 3.14159265358979323846f / 180.0f; };
    
    powerUps.clear();
    powerUps.reserve(bricksPerStory*numberOfStories);
    
    // Initialize bricks to hold all bricks in all stories
    groundLevelBricks.clear();
    groundLevelBricks.reserve(bricksPerStory);

    // Create bricks with alternating colors
    for (int story = 0; story < numberOfStories; ++story) {
        float yOffset = story * brickHeight;  // Vertical offset for each story along the Y-axis

        for (int i = 0; i < bricksPerStory; ++i) {
            float currentAngleDegrees = i * (360.0f / bricksPerStory);
            float currentAngleRadians = toRadians(currentAngleDegrees);
            
            Vector4D position(
                (brickWidth + brickInnerRadius) * sin(currentAngleRadians),
                yOffset,  // Adjust Y for story offset
                (brickWidth + brickInnerRadius) * cos(currentAngleRadians),
                1
            );

            // Loop through the color array
            Vector4D brickColor = colors[(i+story+1) % colors.size()];

            Brick brick(Vector4D(0,0,0,1), brickInnerRadius, brickWidth, brickHeight, bricksPerStory, brickDetail, brickColor);
            if (Vector4D::Equals(brickColor, Vector4D(0, 0, 0, 1), 3)) {
                // if black
                brick.SetDurability(2);
            }
            // Initial transformation: place the brick at the calculated position and rotate to face outward
            Matrix4x4 brickModel = Matrix4x4::Translate(0, position.y, 0)
                * Matrix4x4(1.0).Rotate(-currentAngleRadians, Vector4D(0, 1, 0, 0));
            brick.SetModelMatrix(brickModel);

            shapes.push_back(brick);

            shapes.back().next = nullptr;

            // Link each brick to the corresponding brick in the story above
            if (story > 0)
            {
                shapes[(story - 1) * bricksPerStory + i].next = &shapes.back();
            }
            else
            {
                groundLevelBricks.push_back(&shapes.back());
            }
        }
    }

    // Ground
    model = Matrix4x4(1.0);
    model = model.Rotate(-M_PI / 2, Vector4D(1, 0, 0, 0));

    Circle ground = *new Circle(Vector4D(0, 0, 0, 1), groundDiameter*0.5f, groundDetail, Vector4D(0, 0, 1, 1));
    ground.SetTexture(load_texture(lecture_folder_path / "data" / "textures" / "ground.png"));

    ground.SetModelMatrix(model);
    shapes.push_back(ground);

    // Ball
    model = Matrix4x4(1.0).Translate(1, 0, 1);
    
    Sphere b = *new Sphere(Vector4D(0, 0, 0, 1), ballRadius, ballStacks, ballSectors, Vector4D(0.5, 0.5, 0.5, 1));

    b.SetModelMatrix(model);

    shapes.push_back(b);
    ball = &shapes.back();

    // Paddle creation
    for (int i = 0; i < paddleCount; ++i) {
        float angleDegrees = i * (360.0f / paddleCount);
        float angleRadians = toRadians(angleDegrees);
        
        float x = paddleInnerRadius * cos(angleRadians);
        float z = paddleInnerRadius * sin(angleRadians);
        
        Paddle paddle(Vector4D(x, 0.0f, z, 1), paddleInnerRadius, paddleWidth, paddleHeight, toRadians(paddleLength), paddleDetail, Vector4D(1.0f, 0.0f, 0.0f, 0.0f));

        Matrix4x4 model = Matrix4x4::Rotate(angleRadians, Vector4D(0, 1, 0, 0));
        paddle.SetModelMatrix(model);
        shapes.push_back(paddle);
        paddles.push_back(&shapes.back());
    }
}

void Application::RegeneratePaddles() {
    if (shapes.size() < paddleCount) {
        return;
    }
    
    auto toRadians = [](float degrees) { return degrees * 3.14159265358979323846f / 180.0f; };
    paddles.clear();
    std::vector<Matrix4x4> paddleModelMatrices;
    paddleModelMatrices.reserve(paddleCount);
    paddleModelMatrices.push_back(shapes[shapes.size() - 1].GetModelMatrix());
    paddleModelMatrices.push_back(shapes[shapes.size() - 2].GetModelMatrix());
    paddleModelMatrices.push_back(shapes[shapes.size() - 3].GetModelMatrix());
    shapes.erase(shapes.end() - paddleCount, shapes.end());
    
    for (int i = 0; i < paddleCount; ++i) {
        float angleDegrees = i * (360.0f / paddleCount);
        float angleRadians = toRadians(angleDegrees);
        
        float x = paddleInnerRadius * cos(angleRadians);
        float z = paddleInnerRadius * sin(angleRadians);
        
        Paddle paddle(Vector4D(x, 0.0f, z, 1), paddleInnerRadius, paddleWidth, paddleHeight, toRadians(paddleLength), paddleDetail, Vector4D(1.0f, 0.0f, 0.0f, 0.0f));
        
        Matrix4x4 model = Matrix4x4::Rotate(angleRadians, Vector4D(0, 1, 0, 0));
        
        paddle.SetModelMatrix(paddleModelMatrices[i]);
        
        shapes.push_back(paddle);
        paddles.push_back(&shapes.back());
    }
}

void Application::render() {
    // Sets the clear color.
    glClearColor(0, 0, 0, 1.0f);
    // Clears the window using the above color.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(shader_program);
    
    GLint viewUniformLocation = glGetUniformLocation(shader_program, "view");
    GLint projectionUniformLocation = glGetUniformLocation(shader_program, "projection");
    GLint modelUniformLocation = glGetUniformLocation(shader_program, "model");

    // Check if the uniform locations were found
    if (viewUniformLocation == -1 || projectionUniformLocation == -1 || modelUniformLocation == -1) {
        std::cerr << "Error retrieving uniform location." << std::endl;
    }

    // VIEW
    float distance = 30.0f * sqrt(2.0f);
    float pitchAngleRadians = 45.0f * M_PI / 180.0f;

    Vector4D cameraPosition = (orthographicProj) ?
        Vector4D(0.0f, 5.0f, 0.0001f, 1.0f) :
        Vector4D(
            0.0f, // x
            distance * sin(pitchAngleRadians), // y
            distance * cos(pitchAngleRadians), // z
            1.0f  // w
        );

    Vector4D upVector = Vector4D(0.0f, 1.0f, 0.0f, 0.0f);

    Vector4D cameraTarget = Vector4D(0.0f, 0.0f, 0.0f, 1.0f);
    
    cam.UpdateViewMatrix(cameraPosition, cameraTarget, upVector);

	Matrix4x4 vw = cam.GetView();
    
    if (orthographicProj) {
		vw = vw * vw.Translate(-cameraPosition.x, -cameraPosition.y, -cameraPosition.z);
    }
    else {
        vw = vw * vw.Translate(-cameraPosition.x, -cameraPosition.y, -cameraPosition.z);
    }
    
    float* viewMatrix = vw.ToArray();
    assert(glGetError() == 0U);
    glUniformMatrix4fv(viewUniformLocation, 1, GL_FALSE, viewMatrix);

    // PROJECTION
    /* Projection is by default set to perspective. Projection can be chosen vie keyboard.
     * Key '2' stands for orthographic
     * Kez '1' stands for perspective
     */
    glUseProgram(shader_program);
    float* projMatrix = cam.GetProjection().ToArray();
    assert(glGetError() == 0U);
    glUniformMatrix4fv(projectionUniformLocation, 1, GL_FALSE, projMatrix);

    // RENDERING OBJECTS
    RenderVectorOfObjects(shapes);
    // RENDERING POWERUPS
    RenderVectorOfObjects(powerUps);

}

void Application::RenderVectorOfObjects(std::vector<Shape> vectorOfObjects) {
    GLint useTextureLocation = glGetUniformLocation(shader_program, "useTexture");
    GLint solidColorLocation = glGetUniformLocation(shader_program, "solidColor");
    GLint modelUniformLocation = glGetUniformLocation(shader_program, "model");
    
    for (int i = 0; i < vectorOfObjects.size(); i++) {
        if (vectorOfObjects[i].IsDestroyed())
        {
            continue;
        }
        glDisable(GL_DEPTH_TEST);
        if (i >= 0) {
            glEnable(GL_DEPTH_TEST);
        }
        glBindVertexArray(vectorOfObjects[i].GetVertexArray());
        assert(glGetError() == 0U);
        glBindBuffer(GL_ARRAY_BUFFER, vectorOfObjects[i].GetVertexBuffer());
        assert(glGetError() == 0U);

        if (vectorOfObjects[i].GetTexture() != 0) {
            glUniform1i(useTextureLocation, 1);
            glBindTexture(GL_TEXTURE_2D, vectorOfObjects[i].GetTexture());
        }
        else {
            glUniform1i(useTextureLocation, 0);
            Vector4D colour = vectorOfObjects[i].GetColour();
            glUniform4f(solidColorLocation, colour.x, colour.y, colour.z, colour.w);
        }

        assert(glGetError() == 0U);
        Matrix4x4 model = vectorOfObjects[i].GetModelMatrix();

        glUniformMatrix4fv(modelUniformLocation, 1, GL_FALSE, model.ToArray());
        assert(glGetError() == 0U);

        glUniform3f(glGetUniformLocation(shader_program, "lightColor"), 1, 1, 1);
        glUniform3f(glGetUniformLocation(shader_program, "lightPos"), 0, 12, 10);
        glUniform3fv(glGetUniformLocation(shader_program, "viewPos"), 1, cam.position.ToArray());

        glDrawElements(GL_TRIANGLES, vectorOfObjects[i].GetIndices().size(), GL_UNSIGNED_INT, nullptr);
        assert(glGetError() == 0U);
    }
}

void Application::render_ui() {
	if (isPaused){
        renderTexture(texturePause);
    }
    if (isGameWon) {
        renderTexture(textureWin);
    }
	if (isGameOver){
        renderTexture(textureLoss);
    }
}

void Application::renderTexture(GLuint texture) {
    glUseProgram(texture_program);
    assert(glGetError() == 0U);
    glBindVertexArray(texture_vertex_array);
    assert(glGetError() == 0U);
    glActiveTexture(GL_TEXTURE0);
    assert(glGetError() == 0U);
    glBindTexture(GL_TEXTURE_2D, texture);
    assert(glGetError() == 0U);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
    assert(glGetError() == 0U);
}

void Application::on_resize(int width, int height) {
    float aspectRatio = static_cast<float>(width) / static_cast<float>(height);
    float scale = height / 48.0f;
    // Calls the default implementation to set the class variables.
    IApplication::on_resize(width, height);
    if (orthographicProj) {
        cam.setprojectionMatrixToOrthographic(-scale, scale, -scale * aspectRatio, scale * aspectRatio, -10.0f, 10.0f);
    }
    else {
        cam.setprojectionMatrixToPerspective(45.0f, aspectRatio, 0.1f, 100.0f);
    }
    // Changes the viewport.
    glViewport(0, 0, width, height);
}

void Application::on_mouse_move(double x, double y) {}

void Application::on_mouse_button(int button, int action, int mods) {}

void Application::resetGame() {
    isGameOver = false;
    isGameWon = false;
    isPaused = false;
    playerLives = 3;
    paddleSpeed = initialPaddleSpeed;
    startGame();
}

void Application::on_key_pressed(int key, int scancode, int action, int mods) {
    float aspectRatio = static_cast<float>(width) / static_cast<float>(height);
    float scale = height / 48.0f;
        
    // KEY REALEASE
    if (action == GLFW_RELEASE) {
        switch (key) {
        case GLFW_KEY_LEFT:
            left = false;
            break;
        case GLFW_KEY_RIGHT:
            right = false;
            break;
        }
    }
    
    // KEY PRESS
    if (action == GLFW_PRESS) {
		if (isGameOver || isGameWon) {
            resetGame();
			return;
		}
        switch (key) {
        case GLFW_KEY_2:
            orthographicProj = true;
            cam.setprojectionMatrixToOrthographic(-scale, scale, -scale * aspectRatio, scale * aspectRatio, -10.0f, 10.0f);
            break;

        case GLFW_KEY_1:
            orthographicProj = false;
            cam.setprojectionMatrixToPerspective(45.0f, aspectRatio, 0.1f, 100.0f);
            break;
        case GLFW_KEY_LEFT:
            if (!isGameOver && !isGameWon) {
                left = true;
            }
            break;
        case GLFW_KEY_RIGHT:
            if (!isGameOver && !isGameWon) {
                right = true;
            }
            break;
        case GLFW_KEY_SPACE:
            isBallInGame = true;
            break;
        case GLFW_KEY_P:
            if (!isGameOver && !isGameWon) {
                isPaused = !isPaused;
            }
            break;
        }
    }
}