#version 430 core

layout(location = 0) out vec4 final_color;

in vec2 tex_coord;
uniform sampler2D sample_texture;

void main()
{
    final_color = texture(sample_texture, tex_coord);
}
