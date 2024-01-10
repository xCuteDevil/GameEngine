#include "Shape.hpp"

# define M_PI           3.14159265358979323846  /* pi */

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

    /*std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    float outerRadius = innerRadius + width;
    float brickAngle = (2.0f * M_PI) / count;
    float angleIncrement = brickAngle / detail;*/
    
	Paddle p(GetOrigin(), innerRadius, width, height, 2*M_PI/count, detail, GetColour());
    return { p.GetVertices(), p.GetIndices()};
}

// GenerateMesh function creates a 3D mesh for a paddle-shaped object.
// Parameters:
// innerRadius: The inner radius of the paddle.
// width: Width of the paddle from its inner to outer edge.
// height: Height of the paddle.
// detail: Number of segments used to create the paddle (higher detail means smoother paddle).
// paddleAngle: The angular extent of the paddle in radians.
// Returns a pair of vectors, one for vertices and one for indices, defining the paddle mesh.
std::pair<std::vector<Vertex>, std::vector<unsigned int>> Paddle::GenerateMesh(
    float innerRadius, float width, float height, int detail, float paddleAngle) {

    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    float outerRadius = innerRadius + width;
    float angleIncrement = paddleAngle / detail;
    
	// Create side vertices for the paddle (radially inward and outward)
    for (int i = 0; i <= detail; ++i) {
        float angle = i * angleIncrement;
        float cosA = cos(angle);
        float sinA = sin(angle);

        // Add top and bottom vertices for both outer and inner radius
        vertices.push_back({ outerRadius * cosA, height, outerRadius * sinA, 0, 0, 0 });
        vertices.push_back({ innerRadius * cosA, height, innerRadius * sinA, 0, 0, 0 });
        vertices.push_back({ outerRadius * cosA, 0.0f, outerRadius * sinA, 0, 0, 0 });
        vertices.push_back({ innerRadius * cosA, 0.0f, innerRadius * sinA, 0, 0, 0 });
    }

    // Generate indices for the side faces
    for (int i = 0; i < detail; ++i) {
        int topIndex = i * 4;
        int bottomIndex = topIndex + 2;

        // Outer face
        indices.push_back(topIndex);
        indices.push_back(topIndex + 4);
        indices.push_back(bottomIndex);

        indices.push_back(bottomIndex);
        indices.push_back(topIndex + 4);
        indices.push_back(bottomIndex + 4);

        // Inner face
        indices.push_back(topIndex + 1);
        indices.push_back(bottomIndex + 1);
        indices.push_back(topIndex + 5);

        indices.push_back(bottomIndex + 1);
        indices.push_back(bottomIndex + 5);
        indices.push_back(topIndex + 5);
    }

    // Create vertices for the top and bottom faces of the paddle
    for (int i = 0; i <= detail; ++i) {
        float angle = i * angleIncrement;
        float cosA = cos(angle);
        float sinA = sin(angle);

        // Top face vertices
        vertices.push_back({ outerRadius * cosA, height, outerRadius * sinA, 0, 0, 0 });
        vertices.push_back({ innerRadius * cosA, height, innerRadius * sinA, 0, 0, 0 });

        // Bottom face vertices
        vertices.push_back({ outerRadius * cosA, 0.0f, outerRadius * sinA, 0, 0, 0 });
        vertices.push_back({ innerRadius * cosA, 0.0f, innerRadius * sinA, 0, 0, 0 });
    }

    // Generate indices for the top and bottom faces
    int faceStartIndex = (detail + 1) * 4;
    for (int i = 0; i < (detail); ++i) {
        int topIndex = faceStartIndex + i * 4;
        int bottomIndex = topIndex + 2;
        
        // Top face
        indices.push_back(topIndex);
        indices.push_back(topIndex + 1);
        indices.push_back(topIndex + 5);

        indices.push_back(topIndex);
        indices.push_back(topIndex + 5);
        indices.push_back(topIndex + 4);

        // Bottom face
        indices.push_back(bottomIndex);
        indices.push_back(bottomIndex + 1);
        indices.push_back(bottomIndex + 5);

        indices.push_back(bottomIndex);
        indices.push_back(bottomIndex + 5);
        indices.push_back(bottomIndex + 4);
    }

    // Side caps
    float angle = 0;
    float cosA = cos(angle);
    float sinA = sin(angle);
	int startIndex = vertices.size();
    
    // Top outer and inner vertices
    vertices.push_back({ outerRadius * cosA, height, outerRadius * sinA, 0, 0, 0 });
    vertices.push_back({ innerRadius * cosA, height, innerRadius * sinA, 0, 0, 0 });

    // Bottom outer and inner vertices
    vertices.push_back({ outerRadius * cosA, 0.0f, outerRadius * sinA, 0, 0, 0 });
    vertices.push_back({ innerRadius * cosA, 0.0f, innerRadius * sinA, 0, 0, 0 });
    
    indices.push_back(startIndex);
    indices.push_back(startIndex + 2);
    indices.push_back(startIndex + 1);

    indices.push_back(startIndex + 1);
    indices.push_back(startIndex + 2);
    indices.push_back(startIndex + 3);

    angle = detail * angleIncrement;
    cosA = cos(angle);
    sinA = sin(angle);
    startIndex = vertices.size();

    // Top outer and inner vertices
    vertices.push_back({ outerRadius * cosA, height, outerRadius * sinA, 0, 0, 0 });
    vertices.push_back({ innerRadius * cosA, height, innerRadius * sinA, 0, 0, 0 });

    // Bottom outer and inner vertices
    vertices.push_back({ outerRadius * cosA, 0.0f, outerRadius * sinA, 0, 0, 0 });
    vertices.push_back({ innerRadius * cosA, 0.0f, innerRadius * sinA, 0, 0, 0 });

    indices.push_back(startIndex);
    indices.push_back(startIndex + 1);
    indices.push_back(startIndex + 2);
    

    indices.push_back(startIndex + 1);
    indices.push_back(startIndex + 3);
    indices.push_back(startIndex + 2);
    

    CalculateNormals(vertices, indices);

    return { vertices, indices };
}

Vector4D calculateLocalCenter(const std::vector<Vertex>& vertices) {
    Vector4D center(0, 0, 0, 1);
    for (const auto& vertex : vertices) {
        center.x += vertex.x;
        center.y += vertex.y;
        center.z += vertex.z;
    }
    center.x /= vertices.size();
    center.y /= vertices.size();
    center.z /= vertices.size();
    return center;
}

void Shape::CalculateNormals(std::vector<Vertex>& vertices, const std::vector<unsigned int>& indices) {
    
	// Compute the normal for each face
    for (size_t i = 0; i < indices.size(); i += 3) {
        Vertex& v0 = vertices[indices[i]];
        Vertex& v1 = vertices[indices[i + 1]];
        Vertex& v2 = vertices[indices[i + 2]];

        Vector4D vec0 = { v1.x - v0.x, v1.y - v0.y, v1.z - v0.z };
        Vector4D vec1 = { v2.x - v0.x, v2.y - v0.y, v2.z - v0.z };

        Vector4D normal = vec0.CrossProduct(vec1).UnitVector();
        
        v0.nx = normal.x; v0.ny = normal.y; v0.nz = normal.z;
        v1.nx = normal.x; v1.ny = normal.y; v1.nz = normal.z;
        v2.nx = normal.x; v2.ny = normal.y; v2.nz = normal.z;
    }
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
	return this->shapeMatrix;
}

const std::vector<unsigned int>& Shape::GetIndices() const
{
    return indices;
}

std::vector<Vertex> Shape::GetVertices()
{
    return vertices;
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

Vector4D Shape::GetOrigin() const
{
	return origin;
}

Vector4D Shape::CalculatePosition()
{
    Vector4D localCenter = calculateLocalCenter(vertices);
    auto m = this->shapeMatrix * localCenter;
    auto n = m.ToArray();
	SetPosition(Vector4D(n[0], n[1], n[2], 1.0f));
	return Vector4D(n[0], n[1], n[2], 1.0f);
}

bool Shape::IsDestroyed() const
{
	return destroyed;
}

bool Shape::IsColumnDestroyed()
{
    while (next != nullptr)
    {
        if (!next->destroyed) {
            return false;
        }
        next = next->next;
    }
    return true;
}

void Shape::DestroyBrick() { 
    RecursiveBrickFall(GetModelMatrix(), true);
}

void Shape::RecursiveBrickFall(Matrix4x4 prevBrickModelMatrix, bool destroyThisBrick)
{
    // If we are to destroy this brick and it's not already destroyed, do so.
    if (destroyThisBrick && !destroyed)
    {
        destroyed = true;
        // The current brick has been destroyed, so we don't want to destroy the next one.
        destroyThisBrick = false;
    }

    if (next != nullptr) {
        next->RecursiveBrickFall(GetModelMatrix(), destroyThisBrick);
    }

    if (!destroyed) {
        SetModelMatrix(prevBrickModelMatrix);
    }
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

void Shape::SetPosition(Vector4D p)
{
    position = p;
}

void Shape::StartCooldown(float cooldownDuration)
{
    isOnCooldown = true;
    cooldownTimer = cooldownDuration;
}

void Shape::Update(float delta)
{
    // If the brick is on cooldown, decrement the timer
    if (isOnCooldown)
    {
        cooldownTimer -= delta; // delta is the time passed since the last frame
        if (cooldownTimer <= 0)
        {
            isOnCooldown = false;
        }
    }
}