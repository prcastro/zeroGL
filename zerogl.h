#ifndef ZGL_H
#define ZGL_H

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

/* Utils */
// Define ZGL_DEBUG to enable debug print
#ifdef ZGL_DEBUG
#define ZGL_DEBUG_PRINT(...) printf(__VA_ARGS__)
#else
#define ZGL_DEBUG_PRINT(...) do {} while (0)
#endif

/* Vectors and Matrices */
typedef struct {
  float x, y, z;
} zgl_vec3_t;

typedef struct {
  float x, y, z, w;
} zgl_vec4_t;

typedef struct {
  float data[4][4];
} zgl_mat4x4_t;

static const zgl_mat4x4_t ZGL_IDENTITY_M4 = {{
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
}};

static inline zgl_vec3_t zgl_cross(zgl_vec3_t a, zgl_vec3_t b);
static inline float zgl_dot(zgl_vec3_t a, zgl_vec3_t b);
static inline float zgl_dot_v4(zgl_vec4_t a, zgl_vec4_t b);
static inline float zgl_magnitude(zgl_vec3_t v);
static inline zgl_vec3_t zgl_sub(zgl_vec3_t a, zgl_vec3_t b);
static inline zgl_vec4_t zgl_sub_v4(zgl_vec4_t a, zgl_vec4_t b);
static inline zgl_vec3_t zgl_add(zgl_vec3_t a, zgl_vec3_t b);
static inline zgl_vec3_t zgl_add_three_vec3(zgl_vec3_t a, zgl_vec3_t b, zgl_vec3_t c);
static inline zgl_vec4_t zgl_add_v4(zgl_vec4_t a, zgl_vec4_t b);
static inline zgl_vec3_t zgl_normalize(zgl_vec3_t v);
static inline zgl_vec3_t zgl_mul_scalar(float k, zgl_vec3_t v);
static inline zgl_vec4_t zgl_mul_mat_v4(zgl_mat4x4_t mat4x4, zgl_vec4_t vec4);
static inline zgl_vec3_t zgl_mul_mat_v3(zgl_mat4x4_t mat4x4, zgl_vec3_t v);
static inline zgl_mat4x4_t zgl_mul_mat(zgl_mat4x4_t m1, zgl_mat4x4_t m2);
static inline zgl_mat4x4_t zgl_transpose(zgl_mat4x4_t m);
static inline zgl_mat4x4_t zgl_inverse(zgl_mat4x4_t matrix);
static inline zgl_mat4x4_t zgl_translation_mat(zgl_vec3_t vector);
static inline zgl_mat4x4_t zgl_scale_mat(float scale);
static inline zgl_mat4x4_t zgl_rotx_mat(float degrees);
static inline zgl_mat4x4_t zgl_roty_mat(float degrees);
static inline zgl_mat4x4_t zgl_rotz_mat(float degrees);

/* Quaternions */
typedef struct {
    float w, x, y, z;
} zgl_quaternion_t;

static inline zgl_quaternion_t zgl_quaternion(float degrees, zgl_vec3_t axis);
static inline zgl_quaternion_t zgl_mul_quat(zgl_quaternion_t q1, zgl_quaternion_t q2);
static inline zgl_vec3_t zgl_rotate(zgl_vec3_t v, zgl_quaternion_t q);

/* Colors */

// Pixel formats
#define ZGL_PIXELFORMAT_RGBA8888 0
#define ZGL_PIXELFORMAT_ARGB8888 1
#define ZGL_PIXELFORMAT_BGRA8888 2
#define ZGL_PIXELFORMAT_ABGR8888 3

// Configure pixel format
#ifndef ZGL_PIXELFORMAT
#define ZGL_PIXELFORMAT ZGL_PIXELFORMAT_RGBA8888
#endif // ZGL_PIXELFORMAT

// Default colors
#if  ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_RGBA8888
static const uint32_t ZGL_COLOR_WHITE  = 0xFFFFFFFF;
static const uint32_t ZGL_COLOR_BLACK  = 0x000000FF;
static const uint32_t ZGL_COLOR_RED    = 0xFF0000FF;
static const uint32_t ZGL_COLOR_GREEN  = 0x00FF00FF;
static const uint32_t ZGL_COLOR_BLUE   = 0x0000FFFF;
#elif ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_ARGB8888
static const uint32_t ZGL_COLOR_WHITE  = 0x00FFFFFF;
static const uint32_t ZGL_COLOR_BLACK  = 0x00000000;
static const uint32_t ZGL_COLOR_RED    = 0x00FF0000;
static const uint32_t ZGL_COLOR_GREEN  = 0x0000FF00;
static const uint32_t ZGL_COLOR_BLUE   = 0x000000FF;
#elif ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_BGRA8888
static const uint32_t ZGL_COLOR_WHITE  = 0xFFFFFF00;
static const uint32_t ZGL_COLOR_BLACK  = 0x00000000;
static const uint32_t ZGL_COLOR_RED    = 0x0000FF00;
static const uint32_t ZGL_COLOR_GREEN  = 0x00FF0000;
static const uint32_t ZGL_COLOR_BLUE   = 0xFF000000;
#elif ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_ABGR8888
static const uint32_t ZGL_COLOR_WHITE  = 0x00FFFFFF;
static const uint32_t ZGL_COLOR_BLACK  = 0x00000000;
static const uint32_t ZGL_COLOR_RED    = 0x000000FF;
static const uint32_t ZGL_COLOR_GREEN  = 0x0000FF00;
static const uint32_t ZGL_COLOR_BLUE   = 0x00FF0000;
#endif

static inline uint32_t zgl_color(uint8_t r, uint8_t g, uint8_t b);
static inline void zgl_color_components(uint32_t c, uint8_t* r, uint8_t* g, uint8_t* b);
static inline uint32_t zgl_mul_scalar_color(double x, uint32_t color);
static inline uint32_t zgl_mul_vec3_color(zgl_vec3_t v, uint32_t color);
static inline uint32_t zgl_add_colors(uint32_t c0, uint32_t c1);
static inline uint32_t zgl_color_from_floats(float r, float g, float b);
static inline void zgl_color_to_floats(uint32_t color, float* r, float* g, float* b);

/* 3D Objects */
typedef struct {
    uint32_t* pixels;
    int       width;
    int       height;
} zgl_texture_t;

typedef struct {
    zgl_texture_t color;   // color attachment; also usable as a zgl_texture_t
    float*        depth;   // optional; NULL means no depth buffer
} zgl_framebuffer_t;

typedef struct {
  int     v0, v1, v2;
  int     t0, t1, t2;
  int     n0, n1, n2;
  int     materialIndex;
} zgl_triangle_t;

typedef struct {
    char*        name;
    uint32_t     ambientColor;
    uint32_t     diffuseColor;
    uint32_t     specularColor;
    float        specularExponent;
    zgl_texture_t ambientTexture;
    zgl_texture_t diffuseTexture;
    zgl_texture_t specularTexture;
} zgl_material_t;

typedef struct {
    char*           name;
    int             numVertices;
    zgl_vec3_t*     vertices;
    int             numTextureCoords;
    zgl_vec3_t*     textureCoords;
    int             numNormals;
    zgl_vec3_t*     normals;
    float*          invMagnitudeNormals;
    int             numTriangles;
    zgl_triangle_t* triangles;
    int             numMaterials;
    zgl_material_t* materials;
    zgl_vec3_t      center;
    float           boundsRadius;
} zgl_mesh_t;

typedef struct {
    zgl_mesh_t*  mesh;
    zgl_vec3_t   translation;
    float        scale;
    zgl_mat4x4_t rotation;
    zgl_mat4x4_t transform;
} zgl_object3D_t;

static inline zgl_object3D_t zgl_object(zgl_mesh_t *mesh, zgl_vec3_t translation, float scale, zgl_mat4x4_t rotation);
static inline zgl_vec3_t zgl_mesh_center(zgl_vec3_t* vertices, int numVertices);
static inline float zgl_mesh_bound_radius(zgl_vec3_t* vertices, int numVertices, zgl_vec3_t center);
uint32_t zgl_sample_texture(float u, float v, zgl_texture_t texture, int bilinearFiltering);

/* Lighting */
typedef struct {
    zgl_vec3_t intensity;
} zgl_ambient_light_t;

typedef struct {
    zgl_vec3_t intensity;
    zgl_vec3_t direction;
} zgl_dir_light_t;

typedef struct {
    zgl_vec3_t intensity;
    zgl_vec3_t position;
} zgl_point_light_t;

typedef struct {
    int                  numAmbientLights;
    zgl_ambient_light_t* ambientLights;
    int                  numDirectionalLights;
    zgl_dir_light_t*     directionalLights;
    int                  numPointLights;
    zgl_point_light_t*   pointLights;
} zgl_light_sources_t;

typedef struct {
    zgl_vec3_t diffuse;
    zgl_vec3_t specular;
    zgl_vec3_t ambient;
} zgl_lighting_result_t;

static inline zgl_lighting_result_t zgl_lighting(zgl_vec3_t position, zgl_vec3_t normal, float invMagnitudeNormal, float specularExponent,
                                      zgl_light_sources_t lightSources, uint8_t renderOptions);

/* Camera */
// Camera for left handed coordinate system
typedef struct {
    zgl_vec3_t   position;         // x, y, z coordinates of the camera
    zgl_vec3_t   direction;        // Direction in which the camera is looking
    zgl_vec3_t   up;               // Up direction for the camera (usually [0, 1, 0])
    float        fov;              // Field of view (in degrees)
    float        aspectRatio;      // Aspect ratio of the viewport
    float        nearPlane;        // Distance to the near clipping plane
    float        farPlane;         // Distance to the far clipping plane
    zgl_mat4x4_t viewMatrix;       // View matrix
    zgl_mat4x4_t projectionMatrix; // Projection matrix
    zgl_mat4x4_t viewProjMatrix;   // View * Projection matrix
    zgl_vec4_t   frustumPlanes[6]; // View frustum planes
    float        movementSpeed;
    float        turningSpeed;
} zgl_camera_t;

static inline zgl_camera_t zgl_camera(zgl_vec3_t position, zgl_vec3_t direction, zgl_vec3_t up,
                                      float fov, float aspectRatio, float near, float far,
                                      float movementSpeed, float turningSpeed);

/* Shaders */

// Configuration
#ifndef ZGL_MAX_VERTEX_SHADER_ATTRIBUTES
#define ZGL_MAX_VERTEX_SHADER_ATTRIBUTES 50
#endif // ZGL_MAX_VERTEX_SHADER_ATTRIBUTES

typedef struct {
    zgl_vec3_t position;
    zgl_vec3_t normal;
    zgl_vec3_t textureCoord;
    int        materialIndex;
} zgl_vertex_input_t;

typedef struct {
    zgl_vec4_t position;
    int        numAttributes;
    float      attributes[ZGL_MAX_VERTEX_SHADER_ATTRIBUTES];
    int        numFlatAttributes;
    float      flatAttributes[ZGL_MAX_VERTEX_SHADER_ATTRIBUTES];
} zgl_shader_context_t;

typedef zgl_shader_context_t zgl_vertex_shader_t(void* inputVertex, void* uniformData);
typedef uint32_t zgl_fragment_shader_t(const zgl_shader_context_t* input, void* uniformData);

// Named typed views over zgl_shader_context_t::attributes/flatAttributes.
// Each built-in shader defines a struct of plain floats and casts the float[]
// storage to it, so vertex/fragment pairs can no longer drift on a numeric
// index. The cast is routed through void* (every supported compiler aliases
// a struct-of-floats over a float[]).
#define ZGL_VARYING_FLOATS(T)              (sizeof(T) / sizeof(float))
#define ZGL_VARYING_AS(T, ctx_attrs)       ((T*)       (void*)       (ctx_attrs))
#define ZGL_VARYING_CONST_AS(T, ctx_attrs) ((const T*) (const void*) (ctx_attrs))

typedef struct {
    zgl_mat4x4_t modelViewProjection;
} zgl_basic_uniform_t;

static inline zgl_shader_context_t zgl_basic_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_basic_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

typedef struct {
    zgl_mat4x4_t    modelViewProjection;
    zgl_material_t* materials;
} zgl_colored_uniform_t;

typedef struct {
    float ambient_r,  ambient_g,  ambient_b;
    float diffuse_r,  diffuse_g,  diffuse_b;
    float specular_r, specular_g, specular_b;
} zgl_colored_varying_t;

// Used by zgl_render_triangle, which calls zgl_colored_fragment_shader with
// only the ambient colour slots populated.
typedef struct {
    float r, g, b;
} zgl_render_triangle_varying_t;

static inline zgl_shader_context_t zgl_colored_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_colored_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

typedef struct {
    zgl_mat4x4_t        modelMatrix;
    zgl_mat4x4_t        modelInvRotationMatrixTransposed;
    zgl_mat4x4_t        viewProjectionMatrix;
    int                 bilinearFiltering;
    zgl_material_t*     materials;
} zgl_unlit_uniform_t;

typedef struct {
    float u, v;
} zgl_unlit_varying_t;

typedef struct {
    float nx, ny, nz;
    float ambient_r,  ambient_g,  ambient_b;
    float diffuse_r,  diffuse_g,  diffuse_b;
    float specular_r, specular_g, specular_b;
    float specular_exponent;
    float material_index;
} zgl_unlit_per_tri_t;

static inline zgl_shader_context_t zgl_unlit_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_unlit_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

typedef struct {
    zgl_mat4x4_t        modelMatrix;
    zgl_mat4x4_t        modelInvRotationMatrixTransposed;
    zgl_mat4x4_t        viewProjectionMatrix;
    zgl_light_sources_t lightSources;
    int                 bilinearFiltering;
    zgl_material_t*     materials;
} zgl_flat_uniform_t;

typedef struct {
    float u, v;
} zgl_flat_varying_t;

typedef struct {
    float nx, ny, nz;
    float ambient_r,    ambient_g,    ambient_b;
    float diffuse_r,    diffuse_g,    diffuse_b;
    float specular_r,   specular_g,   specular_b;
    float specular_exponent;
    float light_amb_r,  light_amb_g,  light_amb_b;
    float light_diff_r, light_diff_g, light_diff_b;
    float light_spec_r, light_spec_g, light_spec_b;
    float material_index;
} zgl_flat_per_tri_t;

static inline zgl_shader_context_t zgl_flat_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_flat_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

typedef struct {
    zgl_mat4x4_t        modelMatrix;
    zgl_mat4x4_t        modelInvRotationMatrixTransposed;
    zgl_mat4x4_t        viewProjectionMatrix;
    zgl_light_sources_t lightSources;
    int                 bilinearFiltering;
    zgl_material_t*     materials;
} zgl_gouraud_uniform_t;

typedef struct {
    float nx, ny, nz;
    float u, v;
    float ambient_r,    ambient_g,    ambient_b;
    float diffuse_r,    diffuse_g,    diffuse_b;
    float specular_r,   specular_g,   specular_b;
    float specular_exponent;
    float light_amb_r,  light_amb_g,  light_amb_b;
    float light_diff_r, light_diff_g, light_diff_b;
    float light_spec_r, light_spec_g, light_spec_b;
} zgl_gouraud_varying_t;

typedef struct {
    float material_index;
} zgl_gouraud_per_tri_t;

static inline zgl_shader_context_t zgl_gouraud_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_gouraud_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

typedef struct {
    zgl_mat4x4_t        modelMatrix;
    zgl_mat4x4_t        modelInvRotationMatrixTransposed;
    zgl_mat4x4_t        viewProjectionMatrix;
    zgl_light_sources_t lightSources;
    int                 bilinearFiltering;
    zgl_material_t*     materials;
} zgl_phong_uniform_t;

typedef struct {
    float nx, ny, nz;
    float wx, wy, wz;
    float u, v;
    float ambient_r,  ambient_g,  ambient_b;
    float diffuse_r,  diffuse_g,  diffuse_b;
    float specular_r, specular_g, specular_b;
    float specular_exponent;
} zgl_phong_varying_t;

typedef struct {
    float material_index;
} zgl_phong_per_tri_t;

static inline zgl_shader_context_t zgl_phong_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_phong_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

_Static_assert(ZGL_VARYING_FLOATS(zgl_colored_varying_t)         <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_colored_varying_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");
_Static_assert(ZGL_VARYING_FLOATS(zgl_render_triangle_varying_t) <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_render_triangle_varying_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");
_Static_assert(ZGL_VARYING_FLOATS(zgl_unlit_varying_t)           <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_unlit_varying_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");
_Static_assert(ZGL_VARYING_FLOATS(zgl_unlit_per_tri_t)           <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_unlit_per_tri_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");
_Static_assert(ZGL_VARYING_FLOATS(zgl_flat_varying_t)            <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_flat_varying_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");
_Static_assert(ZGL_VARYING_FLOATS(zgl_flat_per_tri_t)            <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_flat_per_tri_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");
_Static_assert(ZGL_VARYING_FLOATS(zgl_gouraud_varying_t)         <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_gouraud_varying_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");
_Static_assert(ZGL_VARYING_FLOATS(zgl_gouraud_per_tri_t)         <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_gouraud_per_tri_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");
_Static_assert(ZGL_VARYING_FLOATS(zgl_phong_varying_t)           <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_phong_varying_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");
_Static_assert(ZGL_VARYING_FLOATS(zgl_phong_per_tri_t)           <= ZGL_MAX_VERTEX_SHADER_ATTRIBUTES, "zgl_phong_per_tri_t exceeds ZGL_MAX_VERTEX_SHADER_ATTRIBUTES");

/* Rendering */

//  Rendering options
#define ZGL_DIFFUSE_LIGHTING (1 << 0)
#define ZGL_SPECULAR_LIGHTING (1 << 1)
#define ZGL_BACKFACE_CULLING (1 << 2)
#define ZGL_FRUSTUM_CULLING (1 << 3)
#define ZGL_BILINEAR_FILTERING (1 << 4) // TODO: Implement bilinear filtering

static inline void zgl_clear_depth_buffer(zgl_framebuffer_t fb);
static inline void zgl_render_pixel(int i, int j, float z, uint32_t color, zgl_framebuffer_t fb);
static inline void zgl_render_fill(uint32_t color, zgl_framebuffer_t fb);
static inline void zgl_render_line(int x0, int x1, int y0, int y1, uint32_t color, zgl_framebuffer_t fb);
static inline void zgl_render_circle(int x, int y, int r, uint32_t color, zgl_framebuffer_t fb);
static inline void zgl_render_triangle(int x0, int y0, uint32_t color0,
                                       int x1, int y1, uint32_t color1,
                                       int x2, int y2, uint32_t color2,
                                       zgl_framebuffer_t fb);
static inline void zgl_render_object3D(zgl_object3D_t* object, void *uniformData, zgl_camera_t camera, zgl_framebuffer_t fb,
                                       zgl_vertex_shader_t vertexShader, zgl_fragment_shader_t fragmentShader, uint16_t renderOptions);

#endif // ZGL_H


// Include implementation if the client sets ZGL_IMPLEMENTATION, but not twice
#ifdef ZGL_IMPLEMENTATION
#ifndef ZGL_IMPLEMENTATION_INCLUDED
#define ZGL_IMPLEMENTATION_INCLUDED

/* Utils */
#define ZGL__PI 3.14159265358979323846264338327950288
#define ZGL__MIN(a,b) (((a)<(b))?(a):(b))
#define ZGL__MAX(a,b) (((a)>(b))?(a):(b))

/* Vectors and Matrices */

static inline zgl_vec3_t zgl_cross(zgl_vec3_t a, zgl_vec3_t b) {
    zgl_vec3_t result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

static inline float zgl_dot(zgl_vec3_t a, zgl_vec3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline float zgl_dot_v4(zgl_vec4_t a, zgl_vec4_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

static inline float zgl_magnitude(zgl_vec3_t v) {
    float result = sqrt(zgl_dot(v, v));
    assert(result >= 0);
    return result;
}

static inline zgl_vec3_t zgl_sub(zgl_vec3_t a, zgl_vec3_t b) {
    return (zgl_vec3_t) {a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline zgl_vec4_t zgl_sub_v4(zgl_vec4_t a, zgl_vec4_t b) {
    return (zgl_vec4_t) {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}

static inline zgl_vec3_t zgl_add(zgl_vec3_t a, zgl_vec3_t b) {
    return (zgl_vec3_t) {a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline zgl_vec3_t zgl_add_three_vec3(zgl_vec3_t a, zgl_vec3_t b, zgl_vec3_t c) {
    return (zgl_vec3_t) {a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z};
}

static inline zgl_vec4_t zgl_add_v4(zgl_vec4_t a, zgl_vec4_t b) {
    return (zgl_vec4_t) {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

static inline zgl_vec3_t zgl_normalize(zgl_vec3_t v) {
    float mag = zgl_magnitude(v);
    return (zgl_vec3_t) {v.x / mag, v.y / mag, v.z / mag};
}

static inline zgl_vec3_t zgl_mul_scalar(float k, zgl_vec3_t v) {
    return (zgl_vec3_t) {k*v.x, k*v.y, k*v.z};
}

static inline zgl_vec4_t zgl_mul_mat_v4(zgl_mat4x4_t mat4x4, zgl_vec4_t vec4) {
  float result[4] = {0};
  for (int i = 0; i < 4; i++) {
    result[i] += mat4x4.data[i][0]*vec4.x;
    result[i] += mat4x4.data[i][1]*vec4.y;
    result[i] += mat4x4.data[i][2]*vec4.z;
    result[i] += mat4x4.data[i][3]*vec4.w;
  }
  return (zgl_vec4_t) {result[0], result[1], result[2], result[3]};
}

static inline zgl_vec3_t zgl_mul_mat_v3(zgl_mat4x4_t mat4x4, zgl_vec3_t v) {
    zgl_vec4_t vec4 = {v.x, v.y, v.z, 1};
    zgl_vec4_t result = zgl_mul_mat_v4(mat4x4, vec4);
    return (zgl_vec3_t) {result.x, result.y, result.z};
}

static inline zgl_mat4x4_t zgl_mul_mat(zgl_mat4x4_t m1, zgl_mat4x4_t m2) {
    zgl_mat4x4_t result = {0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                result.data[i][j] += m1.data[i][k]*m2.data[k][j];
            }
        }
    }
    return result;
}

static inline zgl_mat4x4_t zgl_transpose(zgl_mat4x4_t m) {
    zgl_mat4x4_t result = {0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result.data[i][j] = m.data[j][i];
        }
    }
    return result;
}

static inline zgl_mat4x4_t zgl_inverse(zgl_mat4x4_t matrix) {
    float m[16];
    for (int i = 0; i < 16; i++) {
        m[i] = matrix.data[i/4][i%4];
    }

    float invOut[16];
    float inv[16], det;

    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return ZGL_IDENTITY_M4;

    det = 1.0 / det;

    for (int i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    zgl_mat4x4_t invMatrix;
    for (int i = 0; i < 4; i++) {
        invMatrix.data[i][0] = invOut[i*4];
        invMatrix.data[i][1] = invOut[i*4+1];
        invMatrix.data[i][2] = invOut[i*4+2];
        invMatrix.data[i][3] = invOut[i*4+3];
    }

    return invMatrix;
}

static inline zgl_mat4x4_t zgl_translation_mat(zgl_vec3_t vector) {
    return (zgl_mat4x4_t) {{
        {1, 0, 0, vector.x},
        {0, 1, 0, vector.y},
        {0, 0, 1, vector.z},
        {0, 0, 0,        1}
    }};
}

static inline zgl_mat4x4_t zgl_scale_mat(float scale) {
    return (zgl_mat4x4_t) {{
        {scale, 0,     0,     0},
        {0,     scale, 0,     0},
        {0,     0,     scale, 0},
        {0,     0,     0,     1}
    }};
}

// TODO: Only use quaternion for rotation
static inline zgl_mat4x4_t zgl_rotx_mat(float degrees) {
    float radians = degrees * ZGL__PI / 180.0f;
    float cos = cosf(radians);
    float sin = sinf(radians);
    return (zgl_mat4x4_t) {{
        { 1, 0,    0,   0 },
        { 0, cos,  sin, 0 },
        { 0, -sin, cos, 0 },
        { 0, 0,    0,   1 }
    }};
}

static inline zgl_mat4x4_t zgl_roty_mat(float degrees) {
    float radians = degrees * ZGL__PI / 180.0f;
    float cos = cosf(radians);
    float sin = sinf(radians);
    return (zgl_mat4x4_t) {{
        { cos, 0, -sin, 0 },
        { 0,   1, 0,    0 },
        { sin, 0, cos,  0 },
        { 0,   0, 0,    1 }
    }};
}

static inline zgl_mat4x4_t zgl_rotz_mat(float degrees) {
    float radians = degrees * ZGL__PI / 180.0f;
    float cos = cosf(radians);
    float sin = sinf(radians);
    return (zgl_mat4x4_t) {{
        { cos,  sin, 0, 0 },
        { -sin, cos, 0, 0 },
        { 0,    0,   1, 0 },
        { 0,    0,   0, 1 }
    }};
}

/* Quaternions */

static inline zgl_quaternion_t zgl_quaternion(float degrees, zgl_vec3_t axis) {
    zgl_vec3_t normalizedAxis = zgl_normalize(axis);
    float angle = degrees * ZGL__PI / 180.0f;
    float halfAngle = angle * 0.5f;
    float sinHalfAngle = sinf(halfAngle);
    return (zgl_quaternion_t) {cosf(halfAngle), normalizedAxis.x * sinHalfAngle, normalizedAxis.y * sinHalfAngle, normalizedAxis.z * sinHalfAngle};
}

static inline zgl_quaternion_t zgl_mul_quat(zgl_quaternion_t q1, zgl_quaternion_t q2) {
    zgl_quaternion_t result;
    result.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    result.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    result.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    result.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return result;
}

static inline zgl_vec3_t zgl_rotate(zgl_vec3_t v, zgl_quaternion_t q) {
    zgl_quaternion_t q_v = {0, v.x, v.y, v.z};
    zgl_quaternion_t q_conjugate = {q.w, -q.x, -q.y, -q.z};
    zgl_quaternion_t rotated = zgl_mul_quat(zgl_mul_quat(q, q_v), q_conjugate);
    return (zgl_vec3_t){rotated.x, rotated.y, rotated.z};
}

/* Colors */

static inline uint32_t zgl_color(uint8_t r, uint8_t g, uint8_t b) {
    #if  ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_RGBA8888
    return 0x000000FF | ((uint32_t)r << 24) | ((uint32_t)g << 16) | ((uint32_t)b << 8);
    #elif ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_ARGB8888
    return 0xFF000000 | ((uint32_t)r << 16) | ((uint32_t)g << 8) |  (uint32_t)b;
    #elif ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_BGRA8888
    return 0x000000FF | ((uint32_t)b << 24) | ((uint32_t)g << 16) | ((uint32_t)r << 8);
    #elif ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_ABGR8888
    return 0xFF000000 | ((uint32_t)b << 16) | ((uint32_t)g << 8) |  (uint32_t)r;
    #endif
}

static inline void zgl_color_components(uint32_t c, uint8_t* r, uint8_t* g, uint8_t* b) {
    #if  ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_RGBA8888
    *r = (c & 0xFF000000) >> 24;
    *g = (c & 0x00FF0000) >> 16;
    *b = (c & 0x0000FF00) >> 8;
    #elif ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_ARGB8888
    *r = (c & 0x00FF0000) >> 16;
    *g = (c & 0x0000FF00) >> 8;
    *b = c & 0x000000FF;
    #elif ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_BGRA8888
    *r = (c & 0x0000FF00) >> 8;
    *g = (c & 0x00FF0000) >> 16;
    *b = (c & 0xFF000000) >> 24;
    #elif ZGL_PIXELFORMAT == ZGL_PIXELFORMAT_ABGR8888
    *r = c & 0x000000FF;
    *g = (c & 0x0000FF00) >> 8;
    *b = (c & 0x00FF0000) >> 16;
    #endif
}

static inline float zgl__clamp(float v, float max) {
    return v > max ? max : v;
}

static inline uint32_t zgl_mul_scalar_color(double x, uint32_t color) {
    uint8_t r, g, b;
    zgl_color_components(color, &r, &g, &b);
    return zgl_color(zgl__clamp(x * r, 255.0), zgl__clamp(x * g, 255.0), zgl__clamp(x * b, 255.0));
}

static inline uint32_t zgl_mul_vec3_color(zgl_vec3_t v, uint32_t color) {
    uint8_t r, g, b;
    zgl_color_components(color, &r, &g, &b);
    return zgl_color(zgl__clamp(v.x * r, 255.0), zgl__clamp(v.y * g, 255.0), zgl__clamp(v.z * b, 255.0));
}

static inline uint32_t zgl_add_colors(uint32_t c0, uint32_t c1) {
    uint8_t r0, g0, b0;
    uint8_t r1, g1, b1;
    zgl_color_components(c0, &r0, &g0, &b0);
    zgl_color_components(c1, &r1, &g1, &b1);
    return zgl_color(zgl__clamp(r0 + r1, 255.0), zgl__clamp(g0 + g1, 255.0), zgl__clamp(b0 + b1, 255.0));
}

static inline uint32_t zgl_mul_colors(uint32_t c0, uint32_t c1) {
    uint8_t r0, g0, b0;
    uint8_t r1, g1, b1;
    zgl_color_components(c0, &r0, &g0, &b0);
    zgl_color_components(c1, &r1, &g1, &b1);
    return zgl_color(zgl__clamp(r0 * r1 / 255.0, 255.0), zgl__clamp(g0 * g1 / 255.0, 255.0), zgl__clamp(b0 * b1 / 255.0, 255.0));
}

static inline uint32_t zgl_color_from_floats(float r, float g, float b) {
    return zgl_color(zgl__clamp(r * 255.0, 255.0), zgl__clamp(g * 255.0, 255.0), zgl__clamp(b * 255.0, 255.0));
}

static inline void zgl_color_to_floats(uint32_t color, float* r, float* g, float* b) {
    uint8_t r8, g8, b8;
    zgl_color_components(color, &r8, &g8, &b8);
    *r = r8 / 255.0f;
    *g = g8 / 255.0f;
    *b = b8 / 255.0f;
}

/* 3D Objects */

static inline zgl_object3D_t zgl_object(zgl_mesh_t *mesh, zgl_vec3_t translation, float scale, zgl_mat4x4_t rotation) {
    zgl_mat4x4_t translationMatrix = zgl_translation_mat(translation);
    zgl_mat4x4_t scaleMatrix = zgl_scale_mat(scale);
    zgl_mat4x4_t transform = zgl_mul_mat(translationMatrix, zgl_mul_mat(rotation, scaleMatrix));
    return (zgl_object3D_t) {mesh, translation, scale, rotation, transform};
}

static inline zgl_vec3_t zgl_mesh_center(zgl_vec3_t* vertices, int numVertices) {
    zgl_vec3_t result = {0, 0, 0};
    for (int i = 0; i < numVertices; i++) {
        result = zgl_add(result, vertices[i]);
    }
    return zgl_mul_scalar(1.0f / numVertices, result);
}

static inline float zgl_mesh_bound_radius(zgl_vec3_t* vertices, int numVertices, zgl_vec3_t center) {
    float result = 0.0f;
    for (int i = 0; i < numVertices; i++) {
        float distance = zgl_magnitude(zgl_sub(vertices[i], center));
        if (distance > result) {
            result = distance;
        }
    }
    return result;
}

uint32_t zgl_sample_texture(float u, float v, zgl_texture_t texture, int bilinearFiltering) {
    int textureWidth = texture.width;
    int textureHeight = texture.height;
    float tex_u = ZGL__MIN(fabs(u * textureWidth), textureWidth - 1);
    float tex_v = ZGL__MIN(fabs(v * textureHeight), textureHeight - 1);
    int floor_u = floor(tex_u);
    int floor_v = floor(tex_v);

    if (bilinearFiltering) {
        float ratio_u = tex_u - floor_u;
        float ratio_v = tex_v - floor_v;
        int next_u = ZGL__MIN(floor_u + 1, textureWidth - 1);
        int next_v = ZGL__MIN(floor_v + 1, textureHeight - 1);
        uint32_t color00 = texture.pixels[floor_v * textureWidth + floor_u];
        uint32_t color10 = texture.pixels[floor_v * textureWidth + next_u];
        uint32_t color01 = texture.pixels[next_v * textureWidth + floor_u];
        uint32_t color11 = texture.pixels[next_v * textureWidth + next_u];
        uint32_t color0 = zgl_add_colors(zgl_mul_scalar_color(1.0f - ratio_u, color00), zgl_mul_scalar_color(ratio_u, color10));
        uint32_t color1 = zgl_add_colors(zgl_mul_scalar_color(1.0f - ratio_u, color01), zgl_mul_scalar_color(ratio_u, color11));
        return zgl_add_colors(zgl_mul_scalar_color(1.0f - ratio_v, color0), zgl_mul_scalar_color(ratio_v, color1));
    } else {
        return texture.pixels[floor_v * textureWidth + floor_u];
    }
}

/* Lighting */

// TODO: Code here is a bit repeated between directional and point lights. Maybe refactor?
static inline zgl_lighting_result_t zgl_lighting(zgl_vec3_t position, zgl_vec3_t normal, float invMagnitudeNormal, float specularExponent,
                   zgl_light_sources_t lightSources, uint8_t renderOptions) {
    int numDirectionalLights = lightSources.numDirectionalLights;
    zgl_dir_light_t* directionalLights = lightSources.directionalLights;
    int numPointLights = lightSources.numPointLights;
    zgl_point_light_t* pointLights = lightSources.pointLights;
    int numAmbientLights = lightSources.numAmbientLights;
    zgl_ambient_light_t* ambientLights = lightSources.ambientLights;

    zgl_vec3_t diffuseIntensity  = {0.0, 0.0, 0.0};
    zgl_vec3_t specularIntensity = {0.0, 0.0, 0.0};
    zgl_vec3_t ambientIntensity  = {0.0, 0.0, 0.0};

    // Directional lights
    for (int i = 0; i < numDirectionalLights; i++) {
        zgl_vec3_t lightDirection = directionalLights[i].direction;
        float magnitudeLightDirection = zgl_magnitude(lightDirection);
        float invMagnitudeLightDirection = 1.0f / magnitudeLightDirection;
        if (renderOptions & ZGL_DIFFUSE_LIGHTING) {
            float cos_alpha = -zgl_dot(lightDirection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            diffuseIntensity.x += ZGL__MAX(cos_alpha, 0.0f) * directionalLights[i].intensity.x;
            diffuseIntensity.y += ZGL__MAX(cos_alpha, 0.0f) * directionalLights[i].intensity.y;
            diffuseIntensity.z += ZGL__MAX(cos_alpha, 0.0f) * directionalLights[i].intensity.z;
        }

        if (renderOptions & ZGL_SPECULAR_LIGHTING) {
            zgl_vec3_t reflection = zgl_sub(zgl_mul_scalar(2 * -zgl_dot(lightDirection, normal), normal), lightDirection);
            float cos_beta = -zgl_dot(reflection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            specularIntensity.x += pow(ZGL__MAX(cos_beta, 0.0f), specularExponent) * directionalLights[i].intensity.x;
            specularIntensity.y += pow(ZGL__MAX(cos_beta, 0.0f), specularExponent) * directionalLights[i].intensity.y;
            specularIntensity.z += pow(ZGL__MAX(cos_beta, 0.0f), specularExponent) * directionalLights[i].intensity.z;
        }
    }

    // Point lights
    for (int i = 0; i < numPointLights; i++) {
        zgl_vec3_t lightDirection = zgl_sub(pointLights[i].position, position);
        float invMagnitudeLightDirection = 1.0f / zgl_magnitude(lightDirection);
        if (renderOptions & ZGL_DIFFUSE_LIGHTING) {
            float cos_alpha = zgl_dot(lightDirection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            diffuseIntensity.x += ZGL__MAX(cos_alpha, 0) * pointLights[i].intensity.x;
            diffuseIntensity.y += ZGL__MAX(cos_alpha, 0) * pointLights[i].intensity.y;
            diffuseIntensity.z += ZGL__MAX(cos_alpha, 0) * pointLights[i].intensity.z;
        }

        if (renderOptions & ZGL_SPECULAR_LIGHTING) {
            zgl_vec3_t reflection = zgl_sub(zgl_mul_scalar(2 * zgl_dot(lightDirection, normal), normal), lightDirection);
            float cos_beta = zgl_dot(reflection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            specularIntensity.x += pow(ZGL__MAX(cos_beta, 0), specularExponent) * pointLights[i].intensity.x;
            specularIntensity.y += pow(ZGL__MAX(cos_beta, 0), specularExponent) * pointLights[i].intensity.y;
            specularIntensity.z += pow(ZGL__MAX(cos_beta, 0), specularExponent) * pointLights[i].intensity.z;
        }
    }

    // Ambient light
    for (int i = 0; i < numAmbientLights; i++) {
        ambientIntensity.x += ambientLights[i].intensity.x;
        ambientIntensity.y += ambientLights[i].intensity.y;
        ambientIntensity.z += ambientLights[i].intensity.z;
    }

    zgl_lighting_result_t result = {diffuseIntensity, specularIntensity, ambientIntensity};
    return result;
}


/* Camera */

static inline zgl_vec4_t zgl__normalize_plane(zgl_vec4_t plane) {
    float mag = sqrt(plane.x*plane.x + plane.y*plane.y + plane.z*plane.z);
    return (zgl_vec4_t) {plane.x / mag, plane.y / mag, plane.z / mag, plane.w / mag};
}

static inline zgl_camera_t zgl_camera(zgl_vec3_t position, zgl_vec3_t direction, zgl_vec3_t up,
                                      float fov, float aspectRatio, float near, float far,
                                      float movementSpeed, float turningSpeed) {
    direction = zgl_normalize(direction);
    up = zgl_normalize(up);
    zgl_vec3_t right = zgl_normalize(zgl_cross(up, direction));
    zgl_vec3_t correctedUp = zgl_cross(direction, right);

    zgl_mat4x4_t rotationMatrix = (zgl_mat4x4_t) {{
        {right.x,       right.y,       right.z,       0},
        {correctedUp.x, correctedUp.y, correctedUp.z, 0},
        {direction.x,   direction.y,   direction.z,   0},
        {0,             0,             0,             1}
    }};

    // Create the translation matrix
    zgl_mat4x4_t translationMatrix = (zgl_mat4x4_t) {{
        {1, 0, 0, -position.x},
        {0, 1, 0, -position.y},
        {0, 0, 1, -position.z},
        {0, 0, 0, 1          }
    }};

    // Create the view matrix
    zgl_mat4x4_t viewMatrix = zgl_mul_mat(rotationMatrix, translationMatrix);

    float fovRadians = fov * ZGL__PI / 180.0;
    float yScale = 1.0 / tan(fovRadians / 2.0);
    float xScale = yScale / aspectRatio;
    float zScale = far / (far - near);
    float zTranslate = -near * zScale;

    zgl_mat4x4_t projectionMatrix = (zgl_mat4x4_t) {{
        {xScale, 0,      0,          0         },
        {0,      yScale, 0,          0         },
        {0,      0,      zScale,     zTranslate},
        {0,      0,      1,          0         }
    }};

    zgl_mat4x4_t viewProjMatrix = zgl_mul_mat(projectionMatrix, viewMatrix);

    // Compute the view frustum planes
    zgl_vec4_t viewProjMatrixCol0 = {viewProjMatrix.data[0][0], viewProjMatrix.data[0][1], viewProjMatrix.data[0][2], viewProjMatrix.data[0][3]};
    zgl_vec4_t viewProjMatrixCol1 = {viewProjMatrix.data[1][0], viewProjMatrix.data[1][1], viewProjMatrix.data[1][2], viewProjMatrix.data[1][3]};
    zgl_vec4_t viewProjMatrixCol2 = {viewProjMatrix.data[2][0], viewProjMatrix.data[2][1], viewProjMatrix.data[2][2], viewProjMatrix.data[2][3]};
    zgl_vec4_t viewProjMatrixCol3 = {viewProjMatrix.data[3][0], viewProjMatrix.data[3][1], viewProjMatrix.data[3][2], viewProjMatrix.data[3][3]};
    zgl_vec4_t leftPlane = zgl__normalize_plane(zgl_add_v4(viewProjMatrixCol3, viewProjMatrixCol0));
    zgl_vec4_t rightPlane = zgl__normalize_plane(zgl_sub_v4(viewProjMatrixCol3, viewProjMatrixCol0));
    zgl_vec4_t bottomPlane = zgl__normalize_plane(zgl_add_v4(viewProjMatrixCol3, viewProjMatrixCol1));
    zgl_vec4_t topPlane = zgl__normalize_plane(zgl_sub_v4(viewProjMatrixCol3, viewProjMatrixCol1));
    zgl_vec4_t nearPlane = zgl__normalize_plane(viewProjMatrixCol2);
    zgl_vec4_t farPlane = zgl__normalize_plane(zgl_add_v4(viewProjMatrixCol3, viewProjMatrixCol2));

    return (zgl_camera_t) {
        .position = position,
        .direction = direction,
        .up = up,
        .fov = fov,
        .aspectRatio = aspectRatio,
        .nearPlane = near,
        .farPlane = far,
        .viewMatrix = viewMatrix,
        .projectionMatrix = projectionMatrix,
        .viewProjMatrix = viewProjMatrix,
        .frustumPlanes = {leftPlane, rightPlane, bottomPlane, topPlane, nearPlane, farPlane},
        .movementSpeed = movementSpeed,
        .turningSpeed = turningSpeed
    };
}

static inline int zgl__tri_in_frustum(zgl_vec4_t v1, zgl_vec4_t v2, zgl_vec4_t v3) {
    // Using NDC coordinates
    if (v1.x < -1 && v2.x < -1 && v3.x < -1) return 0;
    if (v1.x >  1 && v2.x >  1 && v3.x >  1) return 0;
    if (v1.y < -1 && v2.y < -1 && v3.y < -1) return 0;
    if (v1.y >  1 && v2.y >  1 && v3.y >  1) return 0;
    if (v1.z <  0 && v2.z <  0 && v3.z <  0) return 0;
    if (v1.z >  1 && v2.z >  1 && v3.z >  1) return 0;
    return 1;
}

/* Rendering */

static inline int zgl__edge_cross(int ax, int ay, int bx, int by, int px, int py) {
  int abx = bx - ax;
  int aby = by - ay;
  int apx = px - ax;
  int apy = py - ay;
  return abx * apy - aby * apx;
}

static inline void zgl_clear_depth_buffer(zgl_framebuffer_t fb) {
    if (fb.depth == NULL) return;
    for (int i = 0; i < fb.color.width * fb.color.height; i++) {
        fb.depth[i] = FLT_MAX;
    }
}

static inline void zgl_render_pixel(int i, int j, float z, uint32_t color, zgl_framebuffer_t fb) {
    if ((i >= 0) && (i < fb.color.width) && (j >= 0) && (j < fb.color.height)) {
        int position = j * fb.color.width + i;
        fb.color.pixels[position] = color;
        if (fb.depth) fb.depth[position] = z;
    }
}

static inline void zgl_render_fill(uint32_t color, zgl_framebuffer_t fb) {
    for (int i = 0; i < fb.color.width * fb.color.height; i++) {
        fb.color.pixels[i] = color;
    }
}

static inline void zgl_render_line(int x0, int x1, int y0, int y1, uint32_t color, zgl_framebuffer_t fb) {
    int delta_x = (x1 - x0);
    int delta_y = (y1 - y0);
    int longest_side_length = (abs(delta_x) >= abs(delta_y)) ? abs(delta_x) : abs(delta_y);
    float x_inc = delta_x / (float)longest_side_length;
    float y_inc = delta_y / (float)longest_side_length;
    float current_x = x0;
    float current_y = y0;
    for (int i = 0; i <= longest_side_length; i++) {
        zgl_render_pixel(round(current_x), round(current_y), 0.0, color, fb);
        current_x += x_inc;
        current_y += y_inc;
    }
}

static inline void zgl_render_circle(int x, int y, int r, uint32_t color, zgl_framebuffer_t fb) {
    int x1 = x - r;
    int x2 = x + r;
    int y1 = y - r;
    int y2 = y + r;
    for (int j = y1; j < y2; j++) {
        for (int i = x1; i < x2; i++) {
            int dx = i - x;
            int dy = j - y;
            if (dx * dx + dy * dy <= r * r) {
                zgl_render_pixel(i, j, 0.0, color, fb);
            }
        }
    }
}

typedef zgl_shader_context_t zgl_vertex_shader_t(void* inputVertex, void* uniformData);
typedef uint32_t zgl_fragment_shader_t(const zgl_shader_context_t* input, void* uniformData);

static inline void zgl__rasterize_triangle(int x0, int x1, int x2,
                                           int y0, int y1, int y2,
                                           float z0, float z1, float z2,
                                           float invw0, float invw1, float invw2,
                                           int area,
                                           zgl_shader_context_t vertexShaderOutput[3],
                                           zgl_fragment_shader_t fragmentShader, void* uniformData,
                                           zgl_framebuffer_t fb) {
    int x_min = ZGL__MAX(ZGL__MIN(ZGL__MIN(x0, x1), x2), 0);
    int x_max = ZGL__MIN(ZGL__MAX(ZGL__MAX(x0, x1), x2), fb.color.width - 1);
    int y_min = ZGL__MAX(ZGL__MIN(ZGL__MIN(y0, y1), y2), 0);
    int y_max = ZGL__MIN(ZGL__MAX(ZGL__MAX(y0, y1), y2), fb.color.height - 1);

    int delta_w0_col = (y1 - y2);
    int delta_w1_col = (y2 - y0);
    int delta_w2_col = (y0 - y1);
    int delta_w0_row = (x2 - x1);
    int delta_w1_row = (x0 - x2);
    int delta_w2_row = (x1 - x0);

    int w0_row = zgl__edge_cross(x1, y1, x2, y2, x_min, y_min);
    int w1_row = zgl__edge_cross(x2, y2, x0, y0, x_min, y_min);
    int w2_row = zgl__edge_cross(x0, y0, x1, y1, x_min, y_min);

    float invArea = 1.0f / area;

    for (int y = y_min; y <= y_max; y++) {
        int was_inside = 0;
        int w0 = w0_row;
        int w1 = w1_row;
        int w2 = w2_row;
        for (int x = x_min; x <= x_max; x++) {
            // Perspective correct baricentric coordinates
            // TODO: Maybe I can avoid dividing by sum here?
            float alpha = w0 * invArea * invw0;
            float beta  = w1 * invArea * invw1;
            float gamma = w2 * invArea * invw2;
            float sum = alpha + beta + gamma;
            alpha /= sum;
            beta  /= sum;
            gamma /= sum;

            // Check if the fragment is inside the triangle
            int is_inside = alpha >= 0 && beta >= 0 && gamma >= 0;
            if (is_inside) {
                was_inside = 1;

                // Interpolate z
                float z = alpha * z0 + beta * z1 + gamma * z2;

                // Depth test
                if (!fb.depth || z < fb.depth[y * fb.color.width + x]) {
                    // Compute fragment input attributes from the outputs of the vertex shader
                    zgl_shader_context_t fragmentShaderInput = {0};

                    // Interpolate clip position
                    fragmentShaderInput.position.x = alpha * vertexShaderOutput[0].position.x +
                                                        beta  * vertexShaderOutput[1].position.x +
                                                        gamma * vertexShaderOutput[2].position.x;
                    fragmentShaderInput.position.y = alpha * vertexShaderOutput[0].position.y +
                                                        beta  * vertexShaderOutput[1].position.y +
                                                        gamma * vertexShaderOutput[2].position.y;
                    fragmentShaderInput.position.z = alpha * vertexShaderOutput[0].position.z +
                                                        beta  * vertexShaderOutput[1].position.z +
                                                        gamma * vertexShaderOutput[2].position.z;
                    fragmentShaderInput.position.w = 1.0f; // w is always one, because we already did the perspective divide

                    // Interpolate other attributes
                    fragmentShaderInput.numAttributes = vertexShaderOutput[0].numAttributes;
                    for (int i = 0; i < fragmentShaderInput.numAttributes; i++) {
                        fragmentShaderInput.attributes[i] = alpha * vertexShaderOutput[0].attributes[i] +
                                                            beta  * vertexShaderOutput[1].attributes[i] +
                                                            gamma * vertexShaderOutput[2].attributes[i];
                    }

                    // Set flat attributes using the values of the first vertex
                    fragmentShaderInput.numFlatAttributes = vertexShaderOutput[0].numFlatAttributes;
                    for (int i = 0; i < fragmentShaderInput.numFlatAttributes; i++) {
                        fragmentShaderInput.flatAttributes[i] = vertexShaderOutput[0].flatAttributes[i];
                    }

                    uint32_t color = fragmentShader(&fragmentShaderInput, uniformData);
                    zgl_render_pixel(x, y, z, color, fb); // TODO: Avoid scissor test in zgl_render_pixel
                }
            }

            if (!is_inside && was_inside) {
                break;
            }

            w0 += delta_w0_col;
            w1 += delta_w1_col;
            w2 += delta_w2_col;
        }
        w0_row += delta_w0_row;
        w1_row += delta_w1_row;
        w2_row += delta_w2_row;
    }
}

static inline void zgl_render_object3D(zgl_object3D_t* object, void *uniformData, zgl_camera_t camera, zgl_framebuffer_t fb,
                                       zgl_vertex_shader_t vertexShader, zgl_fragment_shader_t fragmentShader, uint16_t renderOptions) {
    zgl_mesh_t* mesh = object->mesh;

    // Don't draw if the object is fully outside the camera frustum
    if (renderOptions & ZGL_FRUSTUM_CULLING) {
        for (int p = 0; p < 6; p++) {
            zgl_vec4_t plane = camera.frustumPlanes[p];
            int isInside = 1;
            zgl_vec4_t center = zgl_mul_mat_v4(object->transform, (zgl_vec4_t) {mesh->center.x, mesh->center.y, mesh->center.z, 1});
            if (zgl_dot_v4(plane, center) < -(object->scale * mesh->boundsRadius)) {
                ZGL_DEBUG_PRINT("DEBUG: Culled object using frustum culling {plane %d}\n", p);
                return;
            }
        }
    }

    for (int tri = 0; tri < mesh->numTriangles; tri++) {
        zgl_triangle_t triangle = mesh->triangles[tri];
        zgl_shader_context_t vertexShaderOutput[3];
        int xs[3];
        int ys[3];
        float zs[3];
        float invws[3];

        // Get vertex data
        zgl_vec3_t vertices[3] = {mesh->vertices[triangle.v0], mesh->vertices[triangle.v1], mesh->vertices[triangle.v2]};
        zgl_vec3_t normals[3] = {0};
        if (mesh->numNormals != 0) {
            normals[0] = mesh->normals[triangle.n0];
            normals[1] = mesh->normals[triangle.n1];
            normals[2] = mesh->normals[triangle.n2];
        }
        zgl_vec3_t textureCoords[3] = {0};
        if (mesh->numTextureCoords != 0) {
            textureCoords[0] = mesh->textureCoords[triangle.t0];
            textureCoords[1] = mesh->textureCoords[triangle.t1];
            textureCoords[2] = mesh->textureCoords[triangle.t2];
        }

        // Get material data
        float diffR = 1.0f;
        float diffG = 1.0f;
        float diffB = 1.0f;
        float specR = 1.0f;
        float specG = 1.0f;
        float specB = 1.0f;
        float specularExponent = 1.0f;
        if (mesh->numMaterials != 0) {
            zgl_material_t material = mesh->materials[triangle.materialIndex];
            zgl_color_to_floats(material.diffuseColor, &diffR, &diffG, &diffB);
            zgl_color_to_floats(material.specularColor, &specR, &specG, &specB);
            specularExponent = material.specularExponent;
        }

        for (int v = 0; v < 3; v++) {
            zgl_vertex_input_t inputVertex = {
                .position = vertices[v],
                .normal = normals[v],
                .textureCoord = textureCoords[v],
                .materialIndex = triangle.materialIndex
            };

            // Vertex shader (local space -> clip space and compute attributes)
            vertexShaderOutput[v] = vertexShader(&inputVertex, uniformData);

            float invw = 1.0f / vertexShaderOutput[v].position.w;

            // Perspective divide (clip space -> NDC)
            vertexShaderOutput[v].position.x *= invw;
            vertexShaderOutput[v].position.y *= invw;
            vertexShaderOutput[v].position.z *= invw;
            vertexShaderOutput[v].position.w = 1.0f;

            // Viewport transform (NDC -> screen space)
            xs[v] = (vertexShaderOutput[v].position.x + 1.0f) * fb.color.width / 2.0f;
            ys[v] = (1.0f - vertexShaderOutput[v].position.y) * fb.color.height / 2.0f;
            zs[v] = vertexShaderOutput[v].position.z; // For z-buffer
            invws[v] = invw; // Store 1/w to avoid divisions later when performing perspective correct interpolation
        }

        int area = zgl__edge_cross(xs[0], ys[0], xs[1], ys[1], xs[2], ys[2]);

        // Backface culling
        if ((renderOptions & ZGL_BACKFACE_CULLING) && area <= 0) {
            ZGL_DEBUG_PRINT("DEBUG: Culled triangle using backface culling\n");
            continue;
        }

        if ((renderOptions & ZGL_FRUSTUM_CULLING) && !zgl__tri_in_frustum(vertexShaderOutput[0].position,
                                                                          vertexShaderOutput[1].position,
                                                                          vertexShaderOutput[2].position)) {
            ZGL_DEBUG_PRINT("DEBUG: Culled triangle using frustum culling\n");
            continue;
        }

        // Rasterization
        zgl__rasterize_triangle(xs[0], xs[1], xs[2],
                                ys[0], ys[1], ys[2],
                                zs[0], zs[1], zs[2],
                                invws[0], invws[1], invws[2],
                                area, vertexShaderOutput, fragmentShader, uniformData, fb);
    }
}

static inline void zgl_render_triangle(int x0, int y0, uint32_t color0,
                                       int x1, int y1, uint32_t color1,
                                       int x2, int y2, uint32_t color2,
                                       zgl_framebuffer_t fb) {
    int area = zgl__edge_cross(x0, y0, x1, y1, x2, y2);

    zgl_shader_context_t shaderContext0 = {0};
    shaderContext0.position = (zgl_vec4_t) {x0, y0, 0, 1};
    shaderContext0.numAttributes = ZGL_VARYING_FLOATS(zgl_render_triangle_varying_t);
    zgl_render_triangle_varying_t* v0 = ZGL_VARYING_AS(zgl_render_triangle_varying_t, shaderContext0.attributes);
    zgl_color_to_floats(color0, &v0->r, &v0->g, &v0->b);

    zgl_shader_context_t shaderContext1 = {0};
    shaderContext1.position = (zgl_vec4_t) {x1, y1, 0, 1};
    shaderContext1.numAttributes = ZGL_VARYING_FLOATS(zgl_render_triangle_varying_t);
    zgl_render_triangle_varying_t* v1 = ZGL_VARYING_AS(zgl_render_triangle_varying_t, shaderContext1.attributes);
    zgl_color_to_floats(color1, &v1->r, &v1->g, &v1->b);

    zgl_shader_context_t shaderContext2 = {0};
    shaderContext2.position = (zgl_vec4_t) {x2, y2, 0, 1};
    shaderContext2.numAttributes = ZGL_VARYING_FLOATS(zgl_render_triangle_varying_t);
    zgl_render_triangle_varying_t* v2 = ZGL_VARYING_AS(zgl_render_triangle_varying_t, shaderContext2.attributes);
    zgl_color_to_floats(color2, &v2->r, &v2->g, &v2->b);

    zgl_shader_context_t shaderContexts[3] = {shaderContext0, shaderContext1, shaderContext2};
    zgl__rasterize_triangle(x0, x1, x2, y0, y1, y2, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                            area, shaderContexts, zgl_colored_fragment_shader, NULL, fb);
}

/* Shaders */

/* Basic shading */
// Draw with a single color, no lighting or textures

static inline zgl_shader_context_t zgl_basic_vertex_shader(void* inputVertex, void* uniformData) {
    zgl_vertex_input_t* inputVertexData = (zgl_vertex_input_t*) inputVertex;
    zgl_shader_context_t result = {0};
    zgl_basic_uniform_t* basicUniformData = (zgl_basic_uniform_t*) uniformData;
    zgl_vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};
    result.position = zgl_mul_mat_v4(basicUniformData->modelViewProjection, inputVertex4);
    return result;
}

static inline uint32_t zgl_basic_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    return ZGL_COLOR_WHITE;
}

/* Colored Shading */
// Draw with the color of the vertex, no lighting or textures

static inline zgl_shader_context_t zgl_colored_vertex_shader(void* inputVertex, void* uniformData) {
    zgl_vertex_input_t* inputVertexData = (zgl_vertex_input_t*) inputVertex;
    zgl_shader_context_t result = {0};
    zgl_colored_uniform_t* basicUniformData = (zgl_colored_uniform_t*) uniformData;
    zgl_vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};
    result.position = zgl_mul_mat_v4(basicUniformData->modelViewProjection, inputVertex4);
    result.numAttributes = ZGL_VARYING_FLOATS(zgl_colored_varying_t);

    zgl_colored_varying_t* v = ZGL_VARYING_AS(zgl_colored_varying_t, result.attributes);
    zgl_material_t mat = basicUniformData->materials[inputVertexData->materialIndex];
    zgl_color_to_floats(mat.ambientColor,  &v->ambient_r,  &v->ambient_g,  &v->ambient_b);
    zgl_color_to_floats(mat.diffuseColor,  &v->diffuse_r,  &v->diffuse_g,  &v->diffuse_b);
    zgl_color_to_floats(mat.specularColor, &v->specular_r, &v->specular_g, &v->specular_b);

    return result;
}

static inline uint32_t zgl_colored_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    const zgl_colored_varying_t* v = ZGL_VARYING_CONST_AS(zgl_colored_varying_t, input->attributes);
    uint32_t ambientColor  = zgl_color_from_floats(v->ambient_r,  v->ambient_g,  v->ambient_b);
    uint32_t diffuseColor  = zgl_color_from_floats(v->diffuse_r,  v->diffuse_g,  v->diffuse_b);
    uint32_t specularColor = zgl_color_from_floats(v->specular_r, v->specular_g, v->specular_b);

    return zgl_mul_colors(ambientColor, zgl_mul_colors(diffuseColor, specularColor));
}

/* Unlit shading */
// Draw with textures, no lighting

static inline zgl_shader_context_t zgl_unlit_vertex_shader(void* inputVertex, void* uniformData) {
    zgl_vertex_input_t* inputVertexData = (zgl_vertex_input_t*) inputVertex;
    zgl_shader_context_t result = {0};
    zgl_unlit_uniform_t* defaultUniformData = (zgl_unlit_uniform_t*) uniformData;
    zgl_vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};
    zgl_vec4_t worldSpaceVertex = zgl_mul_mat_v4(defaultUniformData->modelMatrix, inputVertex4); // Local to world space
    result.position = zgl_mul_mat_v4(defaultUniformData->viewProjectionMatrix, worldSpaceVertex); // World to clip space

    result.numAttributes     = ZGL_VARYING_FLOATS(zgl_unlit_varying_t);
    result.numFlatAttributes = ZGL_VARYING_FLOATS(zgl_unlit_per_tri_t);

    zgl_unlit_varying_t* v = ZGL_VARYING_AS(zgl_unlit_varying_t, result.attributes);
    zgl_unlit_per_tri_t* f = ZGL_VARYING_AS(zgl_unlit_per_tri_t, result.flatAttributes);

    v->u = inputVertexData->textureCoord.x;
    v->v = inputVertexData->textureCoord.y;

    zgl_vec3_t worldSpaceNormal = zgl_mul_mat_v3(defaultUniformData->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    f->nx = worldSpaceNormal.x;
    f->ny = worldSpaceNormal.y;
    f->nz = worldSpaceNormal.z;
    zgl_material_t mat = defaultUniformData->materials[inputVertexData->materialIndex];
    zgl_color_to_floats(mat.ambientColor,  &f->ambient_r,  &f->ambient_g,  &f->ambient_b);
    zgl_color_to_floats(mat.diffuseColor,  &f->diffuse_r,  &f->diffuse_g,  &f->diffuse_b);
    zgl_color_to_floats(mat.specularColor, &f->specular_r, &f->specular_g, &f->specular_b);
    f->specular_exponent = mat.specularExponent;
    f->material_index    = (float) inputVertexData->materialIndex;
    return result;
}

static inline uint32_t zgl_unlit_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    zgl_unlit_uniform_t* uniform = (zgl_unlit_uniform_t*) uniformData;
    const zgl_unlit_varying_t* v = ZGL_VARYING_CONST_AS(zgl_unlit_varying_t, input->attributes);
    const zgl_unlit_per_tri_t* f = ZGL_VARYING_CONST_AS(zgl_unlit_per_tri_t, input->flatAttributes);
    int materialIndex = (int) f->material_index;

    zgl_texture_t ambientTexture = uniform->materials[materialIndex].ambientTexture;
    uint32_t ambientColor = zgl_color_from_floats(f->ambient_r, f->ambient_g, f->ambient_b);
    if (ambientTexture.width != 0 && ambientTexture.height != 0) {
        ambientColor = zgl_mul_colors(ambientColor, zgl_sample_texture(v->u, v->v, ambientTexture, uniform->bilinearFiltering));
    }

    zgl_texture_t diffuseTexture = uniform->materials[materialIndex].diffuseTexture;
    uint32_t diffuseColor = zgl_color_from_floats(f->diffuse_r, f->diffuse_g, f->diffuse_b);
    if (diffuseTexture.width != 0 && diffuseTexture.height != 0) {
        diffuseColor = zgl_mul_colors(diffuseColor, zgl_sample_texture(v->u, v->v, diffuseTexture, uniform->bilinearFiltering));
    }

    zgl_texture_t specularTexture = uniform->materials[materialIndex].specularTexture;
    uint32_t specularColor = zgl_color_from_floats(f->specular_r, f->specular_g, f->specular_b);
    if (specularTexture.width != 0 && specularTexture.height != 0) {
        specularColor = zgl_mul_colors(specularColor, zgl_sample_texture(v->u, v->v, specularTexture, uniform->bilinearFiltering));
    }

    return zgl_mul_colors(ambientColor, zgl_mul_colors(diffuseColor, specularColor));
}


/* Flat shading */
// Compute the lighting at one vertex and use it for the whole triangle

static inline zgl_shader_context_t zgl_flat_vertex_shader(void* inputVertex, void* uniformData) {
    zgl_vertex_input_t* inputVertexData = (zgl_vertex_input_t*) inputVertex;
    zgl_shader_context_t result = {0};
    zgl_flat_uniform_t* defaultUniformData = (zgl_flat_uniform_t*) uniformData;
    zgl_vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};
    zgl_vec4_t worldSpaceVertex = zgl_mul_mat_v4(defaultUniformData->modelMatrix, inputVertex4); // Local to world space
    result.position = zgl_mul_mat_v4(defaultUniformData->viewProjectionMatrix, worldSpaceVertex); // World to clip space

    result.numAttributes     = ZGL_VARYING_FLOATS(zgl_flat_varying_t);
    result.numFlatAttributes = ZGL_VARYING_FLOATS(zgl_flat_per_tri_t);

    zgl_flat_varying_t* v = ZGL_VARYING_AS(zgl_flat_varying_t, result.attributes);
    zgl_flat_per_tri_t* f = ZGL_VARYING_AS(zgl_flat_per_tri_t, result.flatAttributes);

    v->u = inputVertexData->textureCoord.x;
    v->v = inputVertexData->textureCoord.y;

    zgl_vec3_t worldSpaceNormal = zgl_mul_mat_v3(defaultUniformData->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    float invMagNormal = 1.0f / zgl_magnitude(worldSpaceNormal);
    f->nx = worldSpaceNormal.x;
    f->ny = worldSpaceNormal.y;
    f->nz = worldSpaceNormal.z;

    zgl_material_t mat = defaultUniformData->materials[inputVertexData->materialIndex];
    zgl_color_to_floats(mat.ambientColor,  &f->ambient_r,  &f->ambient_g,  &f->ambient_b);
    zgl_color_to_floats(mat.diffuseColor,  &f->diffuse_r,  &f->diffuse_g,  &f->diffuse_b);
    zgl_color_to_floats(mat.specularColor, &f->specular_r, &f->specular_g, &f->specular_b);
    f->specular_exponent = mat.specularExponent;

    zgl_lighting_result_t lightResult = zgl_lighting((zgl_vec3_t) {worldSpaceVertex.x, worldSpaceVertex.y, worldSpaceVertex.z},
                                                     worldSpaceNormal, invMagNormal,
                                                     mat.specularExponent,
                                                     defaultUniformData->lightSources, ZGL_DIFFUSE_LIGHTING | ZGL_SPECULAR_LIGHTING);
    f->light_amb_r  = lightResult.ambient.x;  f->light_amb_g  = lightResult.ambient.y;  f->light_amb_b  = lightResult.ambient.z;
    f->light_diff_r = lightResult.diffuse.x;  f->light_diff_g = lightResult.diffuse.y;  f->light_diff_b = lightResult.diffuse.z;
    f->light_spec_r = lightResult.specular.x; f->light_spec_g = lightResult.specular.y; f->light_spec_b = lightResult.specular.z;

    f->material_index = (float) inputVertexData->materialIndex;
    return result;
}

static inline uint32_t zgl_flat_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    zgl_flat_uniform_t* uniform = (zgl_flat_uniform_t*) uniformData;
    const zgl_flat_varying_t* v = ZGL_VARYING_CONST_AS(zgl_flat_varying_t, input->attributes);
    const zgl_flat_per_tri_t* f = ZGL_VARYING_CONST_AS(zgl_flat_per_tri_t, input->flatAttributes);
    int materialIndex = (int) f->material_index;

    zgl_texture_t ambientTexture = uniform->materials[materialIndex].ambientTexture;
    uint32_t ambientColor = zgl_color_from_floats(f->ambient_r, f->ambient_g, f->ambient_b);
    if (ambientTexture.width != 0 && ambientTexture.height != 0) {
        ambientColor = zgl_mul_colors(ambientColor, zgl_sample_texture(v->u, v->v, ambientTexture, uniform->bilinearFiltering));
    }
    ambientColor = zgl_mul_vec3_color((zgl_vec3_t){f->light_amb_r, f->light_amb_g, f->light_amb_b}, ambientColor);

    zgl_texture_t diffuseTexture = uniform->materials[materialIndex].diffuseTexture;
    uint32_t diffuseColor = zgl_color_from_floats(f->diffuse_r, f->diffuse_g, f->diffuse_b);
    if (diffuseTexture.width != 0 && diffuseTexture.height != 0) {
        diffuseColor = zgl_mul_colors(diffuseColor, zgl_sample_texture(v->u, v->v, diffuseTexture, uniform->bilinearFiltering));
    }
    diffuseColor = zgl_mul_vec3_color((zgl_vec3_t){f->light_diff_r, f->light_diff_g, f->light_diff_b}, diffuseColor);

    zgl_texture_t specularTexture = uniform->materials[materialIndex].specularTexture;
    uint32_t specularColor = zgl_color_from_floats(f->specular_r, f->specular_g, f->specular_b);
    if (specularTexture.width != 0 && specularTexture.height != 0) {
        specularColor = zgl_mul_colors(specularColor, zgl_sample_texture(v->u, v->v, specularTexture, uniform->bilinearFiltering));
    }
    specularColor = zgl_mul_vec3_color((zgl_vec3_t){f->light_spec_r, f->light_spec_g, f->light_spec_b}, specularColor);

    return zgl_add_colors(ambientColor, zgl_add_colors(diffuseColor, specularColor));
}


/* Gouraud shading */
// Compute the lighting at each vertex and interpolate the values at each fragment

static inline zgl_shader_context_t zgl_gouraud_vertex_shader(void* inputVertex, void* uniformData) {
    zgl_vertex_input_t* inputVertexData = (zgl_vertex_input_t*) inputVertex;
    zgl_shader_context_t result = {0};
    zgl_gouraud_uniform_t* defaultUniformData = (zgl_gouraud_uniform_t*) uniformData;
    zgl_vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};
    zgl_vec4_t worldSpaceVertex = zgl_mul_mat_v4(defaultUniformData->modelMatrix, inputVertex4); // Local to world space
    result.position = zgl_mul_mat_v4(defaultUniformData->viewProjectionMatrix, worldSpaceVertex); // World to clip space

    result.numAttributes     = ZGL_VARYING_FLOATS(zgl_gouraud_varying_t);
    result.numFlatAttributes = ZGL_VARYING_FLOATS(zgl_gouraud_per_tri_t);

    zgl_gouraud_varying_t* v = ZGL_VARYING_AS(zgl_gouraud_varying_t, result.attributes);
    zgl_gouraud_per_tri_t* f = ZGL_VARYING_AS(zgl_gouraud_per_tri_t, result.flatAttributes);

    zgl_vec3_t worldSpaceNormal = zgl_mul_mat_v3(defaultUniformData->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    float invMagNormal = 1.0f / zgl_magnitude(worldSpaceNormal);
    v->nx = worldSpaceNormal.x;
    v->ny = worldSpaceNormal.y;
    v->nz = worldSpaceNormal.z;
    v->u  = inputVertexData->textureCoord.x;
    v->v  = inputVertexData->textureCoord.y;

    zgl_material_t mat = defaultUniformData->materials[inputVertexData->materialIndex];
    zgl_color_to_floats(mat.ambientColor,  &v->ambient_r,  &v->ambient_g,  &v->ambient_b);
    zgl_color_to_floats(mat.diffuseColor,  &v->diffuse_r,  &v->diffuse_g,  &v->diffuse_b);
    zgl_color_to_floats(mat.specularColor, &v->specular_r, &v->specular_g, &v->specular_b);
    v->specular_exponent = mat.specularExponent;

    zgl_lighting_result_t lightResult = zgl_lighting((zgl_vec3_t) {worldSpaceVertex.x, worldSpaceVertex.y, worldSpaceVertex.z},
                                                     worldSpaceNormal, invMagNormal,
                                                     mat.specularExponent,
                                                     defaultUniformData->lightSources, ZGL_DIFFUSE_LIGHTING | ZGL_SPECULAR_LIGHTING);
    v->light_amb_r  = lightResult.ambient.x;  v->light_amb_g  = lightResult.ambient.y;  v->light_amb_b  = lightResult.ambient.z;
    v->light_diff_r = lightResult.diffuse.x;  v->light_diff_g = lightResult.diffuse.y;  v->light_diff_b = lightResult.diffuse.z;
    v->light_spec_r = lightResult.specular.x; v->light_spec_g = lightResult.specular.y; v->light_spec_b = lightResult.specular.z;

    f->material_index = (float) inputVertexData->materialIndex;
    return result;
}


static inline uint32_t zgl_gouraud_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    zgl_gouraud_uniform_t* uniform = (zgl_gouraud_uniform_t*) uniformData;
    const zgl_gouraud_varying_t* v = ZGL_VARYING_CONST_AS(zgl_gouraud_varying_t, input->attributes);
    const zgl_gouraud_per_tri_t* f = ZGL_VARYING_CONST_AS(zgl_gouraud_per_tri_t, input->flatAttributes);
    int materialIndex = (int) f->material_index;

    zgl_texture_t ambientTexture = uniform->materials[materialIndex].ambientTexture;
    uint32_t ambientColor = zgl_color_from_floats(v->ambient_r, v->ambient_g, v->ambient_b);
    if (ambientTexture.width != 0 && ambientTexture.height != 0) {
        ambientColor = zgl_mul_colors(ambientColor, zgl_sample_texture(v->u, v->v, ambientTexture, uniform->bilinearFiltering));
    }
    ambientColor = zgl_mul_vec3_color((zgl_vec3_t){v->light_amb_r, v->light_amb_g, v->light_amb_b}, ambientColor);

    zgl_texture_t diffuseTexture = uniform->materials[materialIndex].diffuseTexture;
    uint32_t diffuseColor = zgl_color_from_floats(v->diffuse_r, v->diffuse_g, v->diffuse_b);
    if (diffuseTexture.width != 0 && diffuseTexture.height != 0) {
        diffuseColor = zgl_mul_colors(diffuseColor, zgl_sample_texture(v->u, v->v, diffuseTexture, uniform->bilinearFiltering));
    }
    diffuseColor = zgl_mul_vec3_color((zgl_vec3_t){v->light_diff_r, v->light_diff_g, v->light_diff_b}, diffuseColor);

    zgl_texture_t specularTexture = uniform->materials[materialIndex].specularTexture;
    uint32_t specularColor = zgl_color_from_floats(v->specular_r, v->specular_g, v->specular_b);
    if (specularTexture.width != 0 && specularTexture.height != 0) {
        specularColor = zgl_mul_colors(specularColor, zgl_sample_texture(v->u, v->v, specularTexture, uniform->bilinearFiltering));
    }
    specularColor = zgl_mul_vec3_color((zgl_vec3_t){v->light_spec_r, v->light_spec_g, v->light_spec_b}, specularColor);

    return zgl_add_colors(specularColor, zgl_add_colors(ambientColor, diffuseColor));
}


/* Phong shading */
// Compute the lighting at each fragment

static inline zgl_shader_context_t zgl_phong_vertex_shader(void* inputVertex, void* uniformData) {
    zgl_vertex_input_t* inputVertexData = (zgl_vertex_input_t*) inputVertex;
    zgl_shader_context_t result = {0};
    zgl_phong_uniform_t* uniform = (zgl_phong_uniform_t*) uniformData;
    zgl_vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};
    zgl_vec4_t worldSpaceVertex = zgl_mul_mat_v4(uniform->modelMatrix, inputVertex4); // Local to world space
    result.position = zgl_mul_mat_v4(uniform->viewProjectionMatrix, worldSpaceVertex); // World to clip space

    result.numAttributes     = ZGL_VARYING_FLOATS(zgl_phong_varying_t);
    result.numFlatAttributes = ZGL_VARYING_FLOATS(zgl_phong_per_tri_t);

    zgl_phong_varying_t* v = ZGL_VARYING_AS(zgl_phong_varying_t, result.attributes);
    zgl_phong_per_tri_t* f = ZGL_VARYING_AS(zgl_phong_per_tri_t, result.flatAttributes);

    zgl_vec3_t worldSpaceNormal = zgl_mul_mat_v3(uniform->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    v->nx = worldSpaceNormal.x; v->ny = worldSpaceNormal.y; v->nz = worldSpaceNormal.z;
    v->wx = worldSpaceVertex.x; v->wy = worldSpaceVertex.y; v->wz = worldSpaceVertex.z;
    v->u  = inputVertexData->textureCoord.x;
    v->v  = inputVertexData->textureCoord.y;

    zgl_material_t mat = uniform->materials[inputVertexData->materialIndex];
    zgl_color_to_floats(mat.ambientColor,  &v->ambient_r,  &v->ambient_g,  &v->ambient_b);
    zgl_color_to_floats(mat.diffuseColor,  &v->diffuse_r,  &v->diffuse_g,  &v->diffuse_b);
    zgl_color_to_floats(mat.specularColor, &v->specular_r, &v->specular_g, &v->specular_b);
    v->specular_exponent = mat.specularExponent;

    f->material_index = (float) inputVertexData->materialIndex;
    return result;
}

static inline uint32_t zgl_phong_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    zgl_phong_uniform_t* uniform = (zgl_phong_uniform_t*) uniformData;
    const zgl_phong_varying_t* v = ZGL_VARYING_CONST_AS(zgl_phong_varying_t, input->attributes);
    const zgl_phong_per_tri_t* f = ZGL_VARYING_CONST_AS(zgl_phong_per_tri_t, input->flatAttributes);

    zgl_vec3_t normal   = {v->nx, v->ny, v->nz};
    zgl_vec3_t position = {v->wx, v->wy, v->wz};
    float invMagNormal = 1.0f / zgl_magnitude(normal);
    zgl_lighting_result_t lightingResult = zgl_lighting(position, normal, invMagNormal, v->specular_exponent, uniform->lightSources, ZGL_DIFFUSE_LIGHTING | ZGL_SPECULAR_LIGHTING);
    int materialIndex = (int) f->material_index;

    zgl_texture_t ambientTexture = uniform->materials[materialIndex].ambientTexture;
    uint32_t ambientColor = zgl_color_from_floats(v->ambient_r, v->ambient_g, v->ambient_b);
    if (ambientTexture.width != 0 && ambientTexture.height != 0) {
        ambientColor = zgl_mul_colors(ambientColor, zgl_sample_texture(v->u, v->v, ambientTexture, uniform->bilinearFiltering));
    }
    ambientColor = zgl_mul_vec3_color(lightingResult.ambient, ambientColor);

    zgl_texture_t diffuseTexture = uniform->materials[materialIndex].diffuseTexture;
    uint32_t diffuseColor = zgl_color_from_floats(v->diffuse_r, v->diffuse_g, v->diffuse_b);
    if (diffuseTexture.width != 0 && diffuseTexture.height != 0) {
        diffuseColor = zgl_mul_colors(diffuseColor, zgl_sample_texture(v->u, v->v, diffuseTexture, uniform->bilinearFiltering));
    }
    diffuseColor = zgl_mul_vec3_color(lightingResult.diffuse, diffuseColor);

    zgl_texture_t specularTexture = uniform->materials[materialIndex].specularTexture;
    uint32_t specularColor = zgl_color_from_floats(v->specular_r, v->specular_g, v->specular_b);
    if (specularTexture.width != 0 && specularTexture.height != 0) {
        specularColor = zgl_mul_colors(specularColor, zgl_sample_texture(v->u, v->v, specularTexture, uniform->bilinearFiltering));
    }
    specularColor = zgl_mul_vec3_color(lightingResult.specular, specularColor);

    return zgl_add_colors(specularColor, zgl_add_colors(ambientColor, diffuseColor));
}

#endif // ZGL_IMPLEMENTATION_INCLUDED
#endif // ZGL_IMPLEMENTATION
