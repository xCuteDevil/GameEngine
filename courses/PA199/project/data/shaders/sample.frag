#version 430 core

in VertexData {
    vec3 position_ws;    // The vertex position in world space.
    vec3 normal_ws;      // The vertex normal in world space.
    vec2 tex_coord;      // The vertex texture coordinates.
} in_data;

layout(location = 0) out vec4 final_color;
uniform sampler2D sample_texture;
uniform int useTexture;
uniform vec4 solidColor;

// Light properties
uniform vec3 lightPos; 
uniform vec3 viewPos; 
uniform vec3 lightColor;

void main()
{
    // Apply texture or solid color based on useTexture flag.
    vec3 color;
    if (useTexture == 1) {
        color = texture(sample_texture, in_data.tex_coord).rgb;
    } else {
        color = solidColor.rgb;
    }

    // Ambient light component
    float ambientStrength = 0.1;
    vec3 ambient = ambientStrength * lightColor;
    
    // Diffuse light component
    vec3 norm = normalize(in_data.normal_ws);
    vec3 lightDir = normalize(lightPos - in_data.position_ws);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;
    
    // Specular light component
    float specularStrength = 0.5;
    vec3 viewDir = normalize(viewPos - in_data.position_ws);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = specularStrength * spec * lightColor;  
    
    // Combine
    vec3 result = (ambient + diffuse + specular) * color;
    final_color = vec4(result, 1.0);
}
