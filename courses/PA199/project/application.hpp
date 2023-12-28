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
//#include "Engine/Vertex/Vertex.hpp"


class Application : public IApplication {
    // ----------------------------------------------------------------------------
    // Variables
    // ----------------------------------------------------------------------------
private:
    bool orthographicProj = false;
	
    
    GLuint vertex_shader;
    GLuint fragment_shader;
    GLuint shader_program;
    GLuint vertex_arrays;
    GLuint vertex_buffer;
    GLuint index_buffer;
    GLuint texture;
    Camera cam;
    bool left = false, right = false;
    bool isBallMoving = false;
    float ballRadius = 0.2f;
    
    // ----------------------------------------------------------------------------
    // Variables (Geometry)
    // ----------------------------------------------------------------------------
    std::vector<Shape> shapes;
    std::vector<Shape*> paddles;
    Shape* ball = nullptr;
    
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

    std::vector<Vector4D> colors = {
        Vector4D(1, 1, 0, 1), // Yellow
        Vector4D(0, 1, 0, 1), // Green
        Vector4D(0, 0, 1, 1), // Blue
        Vector4D(1, 0, 0, 1), // Red
        Vector4D(0, 0, 0, 1)  // Black
    };
    
};

