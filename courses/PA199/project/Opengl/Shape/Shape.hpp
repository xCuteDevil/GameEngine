#pragma once

#include <iostream> 
#include "../../Engine/Vector4D/Vector4D.hpp"
#include "../../Engine/Vertex/Vertex.hpp"
#include "../../Engine/Matrix4x4/Matrix4x4.hpp"
#include <vector>
#include <cassert>
#include "glad/glad.h"


//using namespace std;
 
//struct Vertex { float x, y, z, u, v; };

class Shape {
   protected:
      Vector4D origin;
      Vector4D colour;
	  unsigned int texture = 0;
      Matrix4x4 shapeMatrix;

      unsigned int vertex_array;
      unsigned int vertex_buffer;
      unsigned int index_buffer;
      
      std::vector<Vertex> vertices;
      std::vector<unsigned int> indices;
      void CalculateNormals(std::vector<Vertex>& vertices, const std::vector<unsigned int>& indices);
      bool destroyed = false;
      
   public:
	  Vector4D velocity = Vector4D(0, 0, 0);
	  Vector4D force = Vector4D(0.0f, 0.0f, 0.0f);
	  Vector4D position = Vector4D(0.0f, 0.0f, 0.0f);
      bool isDirectionSet = false;
      Shape* next;
       
      Shape(Vector4D _origin, Vector4D _colour){
         origin = _origin;
		 colour = _colour;
      }
      virtual ~Shape() {}  // Virtual destructor
      std::vector<Vertex> GenerateMesh();

      // Getters
      Matrix4x4 GetModelMatrix();
      const std::vector<unsigned int>& GetIndices() const;
      std::vector<Vertex> GetVertices();
      unsigned int GetVertexArray() const;
      unsigned int GetVertexBuffer() const;
      unsigned int GetIndexBuffer() const;
      unsigned int GetTexture() const;
	  Vector4D GetColour() const;
      Vector4D GetOrigin() const;
      Vector4D CalculatePosition();
	  bool IsDestroyed() const;
      bool IsColumnDestroyed();
      
      // Setters
      void SetTexture(unsigned int texture);
      void SetModelMatrix(Matrix4x4 modelMatrix);
      void SetIndexBuffer(unsigned int _index_buffer);
      void GenerateAndBindBuffers();
      void ConfigureVertexAttributes();
      void SetArrays();
	  void SetPosition(Vector4D position);
      
      void DestroyBrick();
	  void RecursiveBrickFall(Matrix4x4 prevBrickModelMatrix, bool alreadyDestroyedOne);
      void Update(float delta);
      void StartCooldown();
      bool isOnCooldown = false;
      float cooldownTimer = 0.0f; // time elapsed
	  float cooldownDuration = 100.0f; // time to wait
};

class Circle : public Shape {
public:
    double radius;
    int numberOfSides;
    Circle(Vector4D _origin, double _radius, int _numberOfSides, Vector4D _colour) :Shape(_origin, _colour)
    {
        radius = _radius;
        numberOfSides = _numberOfSides;
        colour = _colour;
        std::pair<std::vector<Vertex>, std::vector<unsigned int>> meshData = GenerateMesh();
        vertices = std::move(meshData.first);
        indices = std::move(meshData.second);
        SetArrays();
    }
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> GenerateMesh();
};

class Square : public Shape {
public:
    double side;
    Square(Vector4D _origin, double _side, Vector4D _colour) :Shape(_origin, _colour)
    {
        side = _side;
        colour = _colour;
        std::pair<std::vector<Vertex>, std::vector<unsigned int>> meshData = GenerateMesh();
        vertices = std::move(meshData.first);
        indices = std::move(meshData.second);
        SetArrays();
    }
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> GenerateMesh();
};

class Sphere : public Shape {
public:
    double radius;
    int stackCount; // latitude
    int sectorCount; // longitude
    Sphere(Vector4D _origin, double _radius, int _stacks, int _sectors, Vector4D _colour) :Shape(_origin, _colour)
    {
        radius = _radius;
        stackCount = _stacks;
        sectorCount = _sectors;
        colour = _colour;
        std::pair<std::vector<Vertex>, std::vector<unsigned int>> meshData = GenerateMesh();
        vertices = std::move(meshData.first);
        indices = std::move(meshData.second);
        SetArrays();
    }
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> GenerateMesh();
};

class Brick : public Shape {
public:
    float innerRadius;
    float width; // Renamed from thickness for consistency
    float height;
    int detail; // Assuming detail represents the number of subdivisions or detail level
    int count; 

    // Aligning the parameter names and ordering with Paddle
    Brick(Vector4D _origin, float _innerRadius, float _width, float _height, int _count, int _detail, Vector4D _colour)
        : Shape(_origin, _colour), innerRadius(_innerRadius), width(_width), height(_height), count(_count), detail(_detail) {
        std::pair<std::vector<Vertex>, std::vector<unsigned int>> meshData = GenerateMesh(innerRadius, width, height, count, detail);
        vertices.reserve( (detail+1) * 8 + 8);
        vertices = std::move(meshData.first);
        indices = std::move(meshData.second);
        SetArrays();
    }

    std::pair<std::vector<Vertex>, std::vector<unsigned int>> GenerateMesh(float _innerRadius, float _width, float _height, int _count, int _detail);
};

class Paddle : public Shape {
public:
    float innerRadius;
    float width;
    float height;
    int detail;
    float paddleAngle; // Angle in radians for a single paddle

    // Constructor with parameters for paddle geometry
    Paddle(Vector4D _origin, float _innerRadius, float _width, float _height, float _paddleAngle, int _detail, Vector4D _colour)
        : Shape(_origin, _colour), innerRadius(_innerRadius), width(_width), height(_height), paddleAngle(_paddleAngle), detail(_detail) {
        // Generates the mesh data for a single paddle
        std::pair<std::vector<Vertex>, std::vector<unsigned int>> meshData = GenerateMesh(innerRadius, width, height, detail, paddleAngle);
        vertices = std::move(meshData.first);
        indices = std::move(meshData.second);
        SetArrays();
    }

    // Method to generate the mesh for a single paddle
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> GenerateMesh(
        float innerRadius, float width, float height, int detail, float paddleAngle);
};



