#include "Shape.hpp"

# define M_PI           3.14159265358979323846  /* pi */

/*std::pair<std::vector<Vertex>, std::vector<unsigned int>> Square::GenerateMesh()
{

        {   -(float)side/2, -(float)side / 2, 0.0f,       0.0f, 0.0f},
        {  (float)side / 2, -(float)side / 2, 0.0f,       1.0f, 0.0f },
        {  (float)side / 2,  (float)side / 2, 0.0f,       1.0f, 1.0f },
        { -(float)side / 2,  (float)side / 2, 0.0f,       0.0f, 1.0f },
    };
    
    std::vector<unsigned int> indices{
        0, 1, 2, // First triangle
        2, 3, 0  // Second triangle
    };

    return { vertices, indices };
}*/

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Circle::GenerateMesh()
{
    std::vector<Vertex> vertices;
    int numberOfVertices = numberOfSides + 1;

    // Normal for a flat circle facing upwards on the XY plane is (0, 0, 1)
    float nz = 1.0f; // Use -1.0f if the circle should face downwards

    // Center vertex with normal pointing up
    Vertex centerVertex = { 0.0f, 0.0f, 0.0f, 0.5f, 0.5f, 0.0f, 0.0f, nz };
    vertices.push_back(centerVertex);

    // Outer vertices
    for (int i = 0; i < numberOfVertices; ++i) {
        float theta = i * 2.0f * M_PI / numberOfSides;
        float x = radius * cosf(theta);
        float y = radius * sinf(theta);

        // Texture coordinates could be more accurately calculated if needed
        float u = (cosf(theta) + 1.0f) / 2.0f;
        float v = (sinf(theta) + 1.0f) / 2.0f;

        Vertex outerVertex = { x, y, 0.0f, u, v, 0.0f, 0.0f, nz };
        vertices.push_back(outerVertex);
    }

    // Triangle fan indices
    std::vector<unsigned int> indices;
    for (int i = 1; i < numberOfVertices; ++i) {
        indices.push_back(0);
        indices.push_back(i);
        indices.push_back(i % numberOfSides + 1);
    }

    // Correct the last triangle
    indices[indices.size() - 1] = 1;

    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Sphere::GenerateMesh()
{
    std::vector<Vertex> vertices;
    float x, y, z, xy;                              // vertex position
    float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
    float s, t;                                     // vertex texCoord

    float sectorStep = 2 * M_PI / sectorCount;
    float stackStep = M_PI / stackCount;
    float sectorAngle, stackAngle;

    for (int i = 0; i <= stackCount; ++i) {
        stackAngle = M_PI / 2.0f - i * stackStep;  // starting from pi/2 to -pi/2
        xy = radius * cosf(stackAngle);            // r * cos(u)
        z = radius * sinf(stackAngle);             // r * sin(u)

        for (int j = 0; j <= sectorCount; ++j) {
            sectorAngle = j * sectorStep;          // starting from 0 to 2pi

            x = xy * cosf(sectorAngle);            // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle);            // r * cos(u) * sin(v)

            nx = x * lengthInv;
            ny = y * lengthInv;
            nz = z * lengthInv;

            s = (float)j / sectorCount;
            t = (float)i / stackCount;

            // Constructing vertex with position, tex coord, and normal
            Vertex vertex = { x, y, z, s, t, nx, ny, -nz };
            vertices.push_back(vertex);
        }
    }

    std::vector<unsigned int> indices;
    for (int i = 0; i < stackCount; ++i) {
        for (int j = 0; j < sectorCount; ++j) {
            int first = (i * (sectorCount + 1)) + j;
            int second = first + sectorCount + 1;

            // first triangle
            indices.push_back(first);
            indices.push_back(second);
            indices.push_back(first + 1);

            // second triangle
            indices.push_back(first + 1);
            indices.push_back(second);
            indices.push_back(second + 1);
        }
    }
    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Brick::GenerateMesh(
    float innerRadius, float width, float height, int count, int detail) {

    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    float outerRadius = innerRadius + width;
    float brickAngle = (2.0f * M_PI) / count;
    float angleIncrement = brickAngle / detail;

    // Normals for top and bottom faces
    Vector4D normalTop(0.0f, 0.0f, 1.0f, 0.0f);
    Vector4D normalBottom(0.0f, 0.0f, -1.0f, 0.0f);

    // Top and bottom vertices
    for (int i = 0; i <= detail; ++i) {
        float angle = i * angleIncrement;
        float cosA = cos(angle);
        float sinA = sin(angle);

        // Calculate the normals for the curved surfaces
        Vector4D normalOuter(cosA, sinA, 0.0f, 0.0f);  // Radially outward for the outer surface
        Vector4D normalInner(-cosA, -sinA, 0.0f, 0.0f);  // Radially inward for the inner surface

        // Normalize the vectors
        normalOuter = normalOuter.UnitVector();
        normalInner = normalInner.UnitVector();

        // Outer vertices (top)
        vertices.push_back({ outerRadius * cosA, height, outerRadius * sinA, cosA, sinA, (float)normalTop.x, (float)normalTop.y, (float)normalTop.z });

        // Inner vertices (top)
        vertices.push_back({ innerRadius * cosA, height, innerRadius * sinA, -cosA, sinA, (float)normalTop.x, (float)normalTop.y, (float)normalTop.z });

        // Outer vertices (bottom)
        vertices.push_back({ outerRadius * cosA, 0.0f, outerRadius * sinA, cosA, -sinA, (float)normalBottom.x, (float)normalBottom.y, (float)normalBottom.z });

        // Inner vertices (bottom)
        vertices.push_back({ innerRadius * cosA, 0.0f, innerRadius * sinA, -cosA, -sinA, (float)normalBottom.x, (float)normalBottom.y, (float)normalBottom.z });
        
        // End cap normals (for the first and last segments)
        if (i == 0 || i == detail) {
            // The end cap normals are tangential to the end segments, so they will
            // be perpendicular to the radial direction defined by cosA and sinA
            Vector4D normalEndCap = (i == 0) ? Vector4D(-sinA, cosA, 0.0f, 0.0f) : Vector4D(sinA, -cosA, 0.0f, 0.0f);
            normalEndCap = normalEndCap.UnitVector();
            
        }
    
    }

    // Create indices
    for (int i = 0; i < detail; ++i) {
        int outerTop1 = i * 4;
        int outerBottom1 = outerTop1 + 1;
        int innerTop1 = outerTop1 + 2;
        int innerBottom1 = outerTop1 + 3;

        int outerTop2 = outerTop1 + 4;
        int outerBottom2 = outerTop2 + 1;
        int innerTop2 = outerTop2 + 2;
        int innerBottom2 = outerTop2 + 3;

        // Outer face
        indices.push_back(outerTop1);
        indices.push_back(outerBottom2);
        indices.push_back(outerBottom1);

        indices.push_back(outerTop1);
        indices.push_back(outerTop2);
        indices.push_back(outerBottom2);

        // Inner face
        indices.push_back(innerTop1);
        indices.push_back(innerBottom1);
        indices.push_back(innerBottom2);

        indices.push_back(innerTop1);
        indices.push_back(innerBottom2);
        indices.push_back(innerTop2);

        // Top face
        indices.push_back(innerTop1);
        indices.push_back(outerTop2);
        indices.push_back(outerTop1);

        indices.push_back(innerTop1);
        indices.push_back(innerTop2);
        indices.push_back(outerTop2);

        // Bottom face
        indices.push_back(innerBottom1);
        indices.push_back(outerBottom1);
        indices.push_back(outerBottom2);

        indices.push_back(innerBottom1);
        indices.push_back(outerBottom2);
        indices.push_back(innerBottom2);

        // End caps
        if (i == 0) { // First end cap
            indices.push_back(innerTop1);
            indices.push_back(innerBottom1);
            indices.push_back(outerBottom1);

            indices.push_back(innerTop1);
            indices.push_back(outerBottom1);
            indices.push_back(outerTop1);
        }

        if (i == detail - 1) { // Second end cap
            indices.push_back(outerTop2);
            indices.push_back(outerBottom2);
            indices.push_back(innerBottom2);

            indices.push_back(outerTop2);
            indices.push_back(innerBottom2);
            indices.push_back(innerTop2);
        }
    }

    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Paddle::GenerateMesh(
    float innerRadius, float width, float height, int detail, float paddleAngle) {

    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    float outerRadius = innerRadius + width;
    float angleIncrement = paddleAngle / detail; // Increment based on detail

    // Normals for top and bottom faces
    Vector4D normalTop(0.0f, 1.0f, 0.0f);
    Vector4D normalBottom(0.0f, -1.0f, 0.0f);

    // Create vertices for the paddle
    for (int i = 0; i <= detail; ++i) {
        float angle = i * angleIncrement;
        float cosA = cos(angle);
        float sinA = sin(angle);

        // Normals for outer and inner edge
        Vector4D normalOuter(cosA, 0.0f, sinA); // Radially outward
        normalOuter = normalOuter.UnitVector();
        Vector4D normalInner(-cosA, 0.0f, -sinA); // Radially inward
        normalInner = normalInner.UnitVector();
        
        // Outer edge vertices (top and bottom)
        vertices.push_back({ outerRadius * cosA, height, outerRadius * sinA, (float)normalTop.x, (float)normalTop.y, (float)normalTop.z });
        vertices.push_back({ outerRadius * cosA, 0.0f, outerRadius * sinA, (float)normalBottom.x, (float)normalBottom.y, (float)normalBottom.z });

        // Inner edge vertices (top and bottom)
        vertices.push_back({ innerRadius * cosA, height, innerRadius * sinA, (float)normalTop.x, (float)normalTop.y, (float)normalTop.z });
        vertices.push_back({ innerRadius * cosA, 0.0f, innerRadius * sinA, (float)normalBottom.x, (float)normalBottom.y, (float)normalBottom.z });
    }

    // Generate indices for the paddle faces, excluding end caps
    for (int i = 0; i < detail; ++i) {
        int baseIndex = i * 4;

        // Outer face
        indices.push_back(baseIndex);
        indices.push_back(baseIndex + 4);
        indices.push_back(baseIndex + 1);

        indices.push_back(baseIndex + 1);
        indices.push_back(baseIndex + 4);
        indices.push_back(baseIndex + 5);

        // Inner face
        indices.push_back(baseIndex + 2);
        indices.push_back(baseIndex + 3);
        indices.push_back(baseIndex + 6);

        indices.push_back(baseIndex + 3);
        indices.push_back(baseIndex + 7);
        indices.push_back(baseIndex + 6);

        // Top face
        indices.push_back(baseIndex);
        indices.push_back(baseIndex + 2);
        indices.push_back(baseIndex + 6);

        indices.push_back(baseIndex);
        indices.push_back(baseIndex + 6);
        indices.push_back(baseIndex + 4);

        // Bottom face
        indices.push_back(baseIndex + 3);
        indices.push_back(baseIndex + 1);
        indices.push_back(baseIndex + 5);

        indices.push_back(baseIndex + 3);
        indices.push_back(baseIndex + 5);
        indices.push_back(baseIndex + 7);
    }
    
    return { vertices, indices };
}

void Shape::SetTexture(unsigned int texture) {
    this->texture = texture;
}

void Shape::SetModelMatrix(Matrix4x4 modelMatrix)
{
	this->shapeMatrix = modelMatrix;
}

void Shape::SetIndexBuffer(unsigned int _index_buffer)
{
    this->index_buffer = _index_buffer;
}

Matrix4x4 Shape::GetModelMatrix() {
	return shapeMatrix;
}

const std::vector<unsigned int>& Shape::GetIndices() const
{
    return indices;
}

unsigned int Shape::GetVertexArray() const
{
    return vertex_array;
}

unsigned int Shape::GetVertexBuffer() const
{
    return vertex_buffer;
}

unsigned int Shape::GetIndexBuffer() const
{
    return index_buffer;
}

unsigned int Shape::GetTexture() const
{
	return texture;
}

Vector4D Shape::GetColour() const 
{
	return colour;
}

void Shape::GenerateAndBindBuffers() {
    // Generate and bind vertex arrays
    glGenVertexArrays(1, &vertex_array);
    assert(glGetError() == 0U);
    glBindVertexArray(vertex_array);
    assert(glGetError() == 0U);

    // Generate and bind vertex buffers
    glGenBuffers(1, &vertex_buffer);
    assert(glGetError() == 0U);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
    assert(glGetError() == 0U);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), &vertices.front(), GL_STATIC_DRAW);
    assert(glGetError() == 0U);

    // Generate and bind index buffers
    glGenBuffers(1, &index_buffer);
    assert(glGetError() == 0U);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
    assert(glGetError() == 0U);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), &indices.front(), GL_STATIC_DRAW);
    assert(glGetError() == 0U);
}

void Shape::ConfigureVertexAttributes() {
    // Assume that the VAO is already bound
    glEnableVertexAttribArray(0);
    assert(glGetError() == 0U);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), nullptr);
    assert(glGetError() == 0U);

    glEnableVertexAttribArray(1);
    assert(glGetError() == 0U);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3 * sizeof(float)));
    assert(glGetError() == 0U);

    glEnableVertexAttribArray(2);
    assert(glGetError() == 0U);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(5 * sizeof(float)));
    assert(glGetError() == 0U);
}

void Shape::SetArrays() {
    GenerateAndBindBuffers();
    ConfigureVertexAttributes();
}