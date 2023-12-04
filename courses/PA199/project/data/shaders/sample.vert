#version 430 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec2 texcoord;
layout(location = 2) in vec3 normal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out VertexData {
    vec3 position_ws;    // The vertex position in world space.
    vec3 normal_ws;      // The vertex normal in world space.
    vec2 tex_coord;      // The vertex texture coordinates.
} out_data;

void main() {
    // Transform the vertex position into world space.
    vec4 position_world = model * vec4(position, 1.0f);
    out_data.position_ws = position_world.xyz;
    
    out_data.normal_ws = normalize(mat3(model) * normal);
    
    out_data.tex_coord = texcoord;

    // Calculate the position of the vertex in clip space.
    //gl_Position = view * projection * position_world;
    gl_Position = projection * view * model * vec4(position, 1.0f);
    //gl_Position = vec4(position, 1.0f) * model * view * projection;
    //gl_Position = position_world * view * projection;
    //gl_Position = position_world * projection * view;
    //gl_Position = position_world * view * projection;
}
