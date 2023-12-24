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
    , texture(load_texture(lecture_folder_path / "data" / "textures" / "you_win.png"))
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
    glDeleteVertexArrays(1U, &vertex_arrays);
    glDeleteTextures(1, &texture);
    glDeleteProgram(shader_program);
    glDeleteShader(fragment_shader);
    glDeleteShader(vertex_shader);
    for (int i = 0; i < shapes.size(); i++) {
        GLuint array = shapes[i].GetVertexArray();
        glDeleteVertexArrays(1U, &array);
    }
}

// ----------------------------------------------------------------------------
// Methods
// ----------------------------------------------------------------------------

void Application::update(float delta) {}

void Application::startGame() {
    orthographicProj = true;
    cam.setprojectionMatrixToOrthographic(-static_cast<double>(width) / height, static_cast<double>(width) / height, -1.0f, 1.0f, 0.1f, 100.0f);

    shapes.clear();
    Matrix4x4 model = Matrix4x4(1.0f);

    /*** SHAPES ***/

    // Ground
    model = Matrix4x4(1.0);
    model = model.Rotate(-M_PI/2, Vector4D(1, 0, 0, 0));
    
	Circle ground = *new Circle(Vector4D(0, 0, 0, 1), 3.0f, 360, Vector4D(0, 0, 1, 1));
	ground.SetTexture(load_texture(lecture_folder_path / "data" / "textures" / "ground.png"));
    
    ground.SetModelMatrix(model);
    shapes.push_back(ground);

    // Ball
    model = Matrix4x4(1.0).Translate(1, 0, 1) * Matrix4x4(1.0).Scale(0.75, 0.75, 0.75);
    
	Sphere ball = *new Sphere(Vector4D(0, 0, 0, 1), 0.2f, 16, 16, Vector4D(0, 0.5, 1, 1));
    
    ball.SetModelMatrix(model);
    shapes.push_back(ball);
    
    // Conversion helper
    auto toRadians = [](float degrees) { return degrees * 3.14159265358979323846f / 180.0f; };

    // Brick configuration
    const int bricksPerStory = 7;
    const int numberOfStories = 3;
    const float brickHeight = 0.2f;
    const float brickWidth = 0.2f;
    const float radius = 1.5f;
    const int brickDetail = 5;
    const float brickInnerRadius = 0.4f;
    const Vector4D brickColor(0, 0.5, 1, 1);

    // Paddle configuration
    const int paddleCount = 3;
    const float paddleWidth = 0.2f;
    const float paddleHeight = 0.2f;
    const float paddleInnerRadius = 2.5f;
    const int paddleDetail = 36;
    const Vector4D paddleColor(1, 0, 0, 0);

    // Create bricks with alternating colors
    for (int story = 0; story < numberOfStories; ++story) {
        float yOffset = story * brickHeight;  // Vertical offset for each story along the Y-axis

        for (int i = 0; i < bricksPerStory; ++i) {
            float currentAngleDegrees = i * (360.0f / bricksPerStory);
            float currentAngleRadians = toRadians(currentAngleDegrees);

            // Calculate the position in polar coordinates and then convert to Cartesian
            Vector4D position(
                (radius + brickInnerRadius) * sin(currentAngleRadians), // X coordinate
                yOffset,  // Adjust Y for story offset
                (radius + brickInnerRadius) * cos(currentAngleRadians), // Z coordinate
                1
            );

            // Cycle through the color array
            Vector4D brickColor = colors[(i+story+1) % colors.size()];

            Brick brick(Vector4D(0,0,0,1), brickInnerRadius, brickWidth, brickHeight, bricksPerStory, brickDetail, brickColor);

            // Initial transformation: place the brick at the calculated position and rotate to face outward
            Matrix4x4 brickModel = Matrix4x4::Translate(0,position.y,0)
                * Matrix4x4(1.0).Rotate(-currentAngleRadians, Vector4D(0, 1, 0, 0))
                * Matrix4x4::Scale(1.5, 1, 1.5);
            brick.SetModelMatrix(brickModel);

            shapes.push_back(brick);
        }
    }

    // Paddle creation
    for (int i = 0; i < paddleCount; ++i) {
        float angleDegrees = i * (360.0f / paddleCount);
        float angleRadians = toRadians(angleDegrees);

        // Calculate the position on the XZ plane using polar coordinates
        float x = paddleInnerRadius * cos(angleRadians);
        float z = paddleInnerRadius * sin(angleRadians);

        // Create paddle and position it on the XZ plane at the calculated position
        Paddle paddle(Vector4D(x, 0.0f, z, 1), paddleInnerRadius, paddleWidth, paddleHeight, toRadians(60.0f), paddleDetail, paddleColor);

        // Rotate paddle around the global Y-axis to face outward from the center
        Matrix4x4 model = Matrix4x4::Rotate(angleRadians, Vector4D(0, 1, 0, 0));

        // Apply the model matrix to the paddle
        paddle.SetModelMatrix(model);

        // Add the paddle to the shapes vector
        shapes.push_back(paddle);
    }
}

void Application::render() {
    // Sets the clear color.
    glClearColor(red, green, blue, 1.0f);
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
    float distance = 5.0f * sqrt(2.0f);
    float pitchAngleRadians = 45.0f * M_PI / 180.0f;

    Vector4D cameraPosition = (orthographicProj) ?
        Vector4D(0.0f, 2.0f, 0.0000001f, 1.0f) :
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
        //vw = vw.Rotate(M_PI/2, Vector4D(0, 1, 0, 0));
        //vw = vw.Rotate(-M_PI, Vector4D(0, 1, 0, 0));
    }
    else {
        vw = vw * vw.Translate(-cameraPosition.x, -cameraPosition.y, -cameraPosition.z);
        //vw = vw.Rotate(-M_PI, Vector4D(1, 0, 0, 0));
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

    GLint useTextureLocation = glGetUniformLocation(shader_program, "useTexture");
    GLint solidColorLocation = glGetUniformLocation(shader_program, "solidColor");
    
    // RENDERING OBJECTS
    for (int i = 0; i < shapes.size(); i++) {
        glDisable(GL_DEPTH_TEST);
        if (i > 0) {
            glEnable(GL_DEPTH_TEST);
        }
        glBindVertexArray(shapes[i].GetVertexArray());
        assert(glGetError() == 0U);
		glBindBuffer(GL_ARRAY_BUFFER, shapes[i].GetVertexBuffer());
        assert(glGetError() == 0U);
        
		if (shapes[i].GetTexture() != 0){
            glUniform1i(useTextureLocation, 1);
            glBindTexture(GL_TEXTURE_2D, shapes[i].GetTexture());
        }
		else{
            glUniform1i(useTextureLocation, 0);
            Vector4D colour = shapes[i].GetColour();
            glUniform4f(solidColorLocation, colour.x, colour.y, colour.z, colour.w);
        }

        assert(glGetError() == 0U);
        Matrix4x4 model = shapes[i].GetModelMatrix();
        
        glUniformMatrix4fv(modelUniformLocation, 1, GL_FALSE, model.ToArray());
        assert(glGetError() == 0U);

        glUniform3f(glGetUniformLocation(shader_program, "lightColor"), 1, 1, 1);
        glUniform3f(glGetUniformLocation(shader_program, "lightPos"), 0, 12, 10);
        glUniform3fv(glGetUniformLocation(shader_program, "viewPos"), 1, cam.position.ToArray());

        glDrawElements(GL_TRIANGLES, shapes[i].GetIndices().size(), GL_UNSIGNED_INT, nullptr);
        assert(glGetError() == 0U);
    }
}

void Application::render_ui() {}

void Application::on_resize(int width, int height) {
    // Calls the default implementation to set the class variables.
    IApplication::on_resize(width, height);
    // Changes the viewport.
    glViewport(0, 0, width, height);
}

void Application::on_mouse_move(double x, double y) {}

void Application::on_mouse_button(int button, int action, int mods) {}

void Application::on_key_pressed(int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS) {
        red = 0;
        green = 0;
        blue = 0;

        float aspectRatio = static_cast<float>(width) / static_cast<float>(height);
        float scale = height / 256.0f;
        
        //Camera cam(45.0f, 1.0f, 0.1f, 10.0f);
        switch (key) {
            case GLFW_KEY_R: {
                red = 0;
                green = 0;
                blue = 0;
                break;
            }
            case GLFW_KEY_G: {
                green = 1;
                break;
            }
            case GLFW_KEY_B: {
                blue = 1;
                break;
            }
            case GLFW_KEY_2: {
                orthographicProj = true;
                cam.setprojectionMatrixToOrthographic(-scale, scale, -scale * aspectRatio, scale * aspectRatio, -10.0f, 10.0f);
                break;
            }
            case GLFW_KEY_1:
                orthographicProj = false;
                cam.setprojectionMatrixToPerspective(45.0f, aspectRatio, 0.1f, 100.0f);
                break;
            }
    }
}

/* The function write given string into a .txt file for debugging purposes.
 * Destination 'out\build\x64 - Debug\courses\PA199\project'
 */
void Application::debug(std::string message) {
    std::ofstream myfile;
    myfile.open("LAUS.txt", std::ios::out | std::ios::app | std::ios::binary);

    myfile << message + "\n";
    myfile.close();
}

/* The function returns 'vector<GLuint>', which contains circle vertex indexing.
 * The sequence is (0,x,x+1), where x goes from 1 to 'numberOfVertices'
 */
std::vector<GLuint> CircleIndexing(int numberOfVertices)
{
    std::vector<GLuint> indices;

    for (GLuint i = 1; i < numberOfVertices + 1; i++) {
        
        indices.push_back(0U);
        indices.push_back(i);
        indices.push_back(i + 1);
    }
    return indices;
}

std::vector<GLuint> SphereIndexing(Sphere s)
{
    std::vector<GLuint> indices;
    std::vector<int> lineIndices;
    int k1, k2;
    for (int i = 0; i < s.stackCount; ++i)
    {
        k1 = i * (s.sectorCount + 1);     // beginning of current stack
        k2 = k1 + s.sectorCount + 1;      // beginning of next stack

        for (int j = 0; j < s.sectorCount; ++j, ++k1, ++k2)
        {
            // 2 triangles per sector excluding first and last stacks
            // k1 => k2 => k1+1
            if (i != 0)
            {
                indices.push_back(k1);
                indices.push_back(k2);
                indices.push_back(k1 + 1);
            }

            // k1+1 => k2 => k2+1
            if (i != (s.stackCount - 1))
            {
                indices.push_back(k1 + 1);
                indices.push_back(k2);
                indices.push_back(k2 + 1);
            }

            // store indices for lines
            // vertical lines for all stacks, k1 => k2
            lineIndices.push_back(k1);
            lineIndices.push_back(k2);
            if (i != 0)  // horizontal lines except 1st stack, k1 => k+1
            {
                lineIndices.push_back(k1);
                lineIndices.push_back(k1 + 1);
            }
        }
    }
    return indices;
}
