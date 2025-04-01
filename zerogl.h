#ifndef ZEROGL_H
#define ZEROGL_H

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

static const zgl_mat4x4_t IDENTITY_M4x4 = {{
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
#define ZEROGL_PIXELFORMAT_RGBA8888 0
#define ZEROGL_PIXELFORMAT_ARGB8888 1
#define ZEROGL_PIXELFORMAT_BGRA8888 2
#define ZEROGL_PIXELFORMAT_ABGR8888 3

// Configure pixel format
#ifndef ZEROGL_PIXELFORMAT
#define ZEROGL_PIXELFORMAT ZEROGL_PIXELFORMAT_RGBA8888
#endif // ZEROGL_PIXELFORMAT

// Default colors
#if  ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_RGBA8888
static const uint32_t ZGL_COLOR_WHITE  = 0xFFFFFFFF;
static const uint32_t ZGL_COLOR_BLACK  = 0x000000FF;
static const uint32_t ZGL_COLOR_RED    = 0xFF0000FF;
static const uint32_t ZGL_COLOR_GREEN  = 0x00FF00FF;
static const uint32_t ZGL_COLOR_BLUE   = 0x0000FFFF;
#elif ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_ARGB8888
static const uint32_t ZGL_COLOR_WHITE  = 0x00FFFFFF;
static const uint32_t ZGL_COLOR_BLACK  = 0x00000000;
static const uint32_t ZGL_COLOR_RED    = 0x00FF0000;
static const uint32_t ZGL_COLOR_GREEN  = 0x0000FF00;
static const uint32_t ZGL_COLOR_BLUE   = 0x000000FF;
#elif ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_BGRA8888
static const uint32_t ZGL_COLOR_WHITE  = 0xFFFFFF00;
static const uint32_t ZGL_COLOR_BLACK  = 0x00000000;
static const uint32_t ZGL_COLOR_RED    = 0x0000FF00;
static const uint32_t ZGL_COLOR_GREEN  = 0x00FF0000;
static const uint32_t ZGL_COLOR_BLUE   = 0xFF000000;
#elif ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_ABGR8888
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
    uint32_t* frameBuffer;
    int       width;
    int       height;
    int       hasDepthBuffer;
    float*    depthBuffer;
} zgl_canvas_t;

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
    zgl_canvas_t ambientTexture;
    zgl_canvas_t diffuseTexture;
    zgl_canvas_t specularTexture;
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
uint32_t zgl_sample_texture(float u, float v, zgl_canvas_t texture, int bilinearFiltering);

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
    zgl_vec4_t   fustrumPlanes[6]; // View frustum planes
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

typedef struct {
    zgl_mat4x4_t modelviewprojection;
} zgl_basic_uniform_t;

static inline zgl_shader_context_t zgl_basic_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_basic_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

typedef struct {
    zgl_mat4x4_t    modelviewprojection;
    zgl_material_t* materials;
} zgl_colored_uniform_t;

static inline zgl_shader_context_t zgl_colored_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_colored_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

typedef struct {
    zgl_mat4x4_t        modelMatrix;
    zgl_mat4x4_t        modelInvRotationMatrixTransposed;
    zgl_mat4x4_t        viewProjectionMatrix;
    int                 bilinearFiltering;
    zgl_material_t*     materials;
} zgl_unlit_uniform_t;

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

static inline zgl_shader_context_t zgl_flat_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_flat_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

typedef struct {
    zgl_mat4x4_t        modelMatrix;
    zgl_mat4x4_t        modelInvRotationMatrixTransposed;
    zgl_mat4x4_t        viewProjectionMatrix;
    zgl_light_sources_t lightSources;
    int                 bilinearFiltering;
    zgl_material_t*     materials;
} zgl_gourard_uniform_t;

static inline zgl_shader_context_t zgl_gourard_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_gourard_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

typedef struct {
    zgl_mat4x4_t        modelMatrix;
    zgl_mat4x4_t        modelInvRotationMatrixTransposed;
    zgl_mat4x4_t        viewProjectionMatrix;
    zgl_light_sources_t lightSources;
    int                 bilinearFiltering;
    zgl_material_t*     materials;
} zgl_phong_uniform_t;

static inline zgl_shader_context_t zgl_phong_vertex_shader(void* inputVertex, void* uniformData);
static inline uint32_t zgl_phong_fragment_shader(const zgl_shader_context_t* input, void* uniformData);

/* Rendering */

//  Rendering options
#define ZGL_DIFFUSE_LIGHTING (1 << 0)
#define ZGL_SPECULAR_LIGHTING (1 << 1)
#define ZGL_BACKFACE_CULLING (1 << 2)
#define ZGL_FUSTRUM_CULLING (1 << 3)
#define ZGL_BILINEAR_FILTERING (1 << 4) // TODO: Implement bilinear filtering

static inline void zgl_clear_depth_buffer(zgl_canvas_t canvas);
static inline void zgl_render_pixel(int i, int j, float z, uint32_t color, zgl_canvas_t canvas);
static inline void zgl_render_fill(uint32_t color, zgl_canvas_t canvas);
static inline void zgl_render_line(int x0, int x1, int y0, int y1, uint32_t color, zgl_canvas_t canvas);
static inline void zgl_render_circle(int x, int y, int r, uint32_t color, zgl_canvas_t canvas);
static inline void zgl_render_triangle(int x0, int y0, uint32_t color0,
                                       int x1, int y1, uint32_t color1,
                                       int x2, int y2, uint32_t color2,
                                       zgl_canvas_t canvas);
static inline void zgl_render_object3D(zgl_object3D_t* object, void *uniformData, zgl_camera_t camera, zgl_canvas_t canvas,
                                       zgl_vertex_shader_t vertexShader, zgl_fragment_shader_t fragmentShader, uint16_t renderOptions);

#endif // ZEROGL_H


// Include implementation if the client sets ZEROGL_IMPLEMENTATION, but not twice
#ifdef ZEROGL_IMPLEMENTATION
#ifndef ZEROGL_IMPLEMENTATION_INCLUDED
#define ZEROGL_IMPLEMENTATION_INCLUDED

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
        return IDENTITY_M4x4;

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
    #if  ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_RGBA8888
    return 0x000000FF | ((uint32_t)r << 24) | ((uint32_t)g << 16) | ((uint32_t)b << 8);
    #elif ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_ARGB8888
    return 0xFF000000 | ((uint32_t)r << 16) | ((uint32_t)g << 8) |  (uint32_t)b;
    #elif ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_BGRA8888
    return 0x000000FF | ((uint32_t)b << 24) | ((uint32_t)g << 16) | ((uint32_t)r << 8);
    #elif ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_ABGR8888
    return 0xFF000000 | ((uint32_t)b << 16) | ((uint32_t)g << 8) |  (uint32_t)r;
    #endif
}

static inline void zgl_color_components(uint32_t c, uint8_t* r, uint8_t* g, uint8_t* b) {
    #if  ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_RGBA8888
    *r = (c & 0xFF000000) >> 24;
    *g = (c & 0x00FF0000) >> 16;
    *b = (c & 0x0000FF00) >> 8;
    #elif ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_ARGB8888
    *r = (c & 0x00FF0000) >> 16;
    *g = (c & 0x0000FF00) >> 8;
    *b = c & 0x000000FF;
    #elif ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_BGRA8888
    *r = (c & 0x0000FF00) >> 8;
    *g = (c & 0x00FF0000) >> 16;
    *b = (c & 0xFF000000) >> 24;
    #elif ZEROGL_PIXELFORMAT == ZEROGL_PIXELFORMAT_ABGR8888
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

uint32_t zgl_sample_texture(float u, float v, zgl_canvas_t texture, int bilinearFiltering) {
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
        uint32_t color00 = texture.frameBuffer[floor_v * textureWidth + floor_u];
        uint32_t color10 = texture.frameBuffer[floor_v * textureWidth + next_u];
        uint32_t color01 = texture.frameBuffer[next_v * textureWidth + floor_u];
        uint32_t color11 = texture.frameBuffer[next_v * textureWidth + next_u];
        uint32_t color0 = zgl_add_colors(zgl_mul_scalar_color(1.0f - ratio_u, color00), zgl_mul_scalar_color(ratio_u, color10));
        uint32_t color1 = zgl_add_colors(zgl_mul_scalar_color(1.0f - ratio_u, color01), zgl_mul_scalar_color(ratio_u, color11));
        return zgl_add_colors(zgl_mul_scalar_color(1.0f - ratio_v, color0), zgl_mul_scalar_color(ratio_v, color1));
    } else {
        return texture.frameBuffer[floor_v * textureWidth + floor_u];
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
        .fustrumPlanes = {leftPlane, rightPlane, bottomPlane, topPlane, nearPlane, farPlane},
        .movementSpeed = movementSpeed,
        .turningSpeed = turningSpeed
    };
}

static inline int zgl__tri_in_fustrum(zgl_vec4_t v1, zgl_vec4_t v2, zgl_vec4_t v3) {
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

static inline void zgl_clear_depth_buffer(zgl_canvas_t canvas) {
    for (int i = 0; i < canvas.width * canvas.height; i++) {
        canvas.depthBuffer[i] = FLT_MAX;
    }
}

static inline void zgl_render_pixel(int i, int j, float z, uint32_t color, zgl_canvas_t canvas) {
    if ((i >= 0) && (i < canvas.width) && (j >= 0) && (j < canvas.height)) {
        int position = j * canvas.width + i;
        canvas.frameBuffer[position] = color;
        canvas.depthBuffer[position] = z;
    }
}

static inline void zgl_render_fill(uint32_t color, zgl_canvas_t canvas) {
    for (int i = 0; i < canvas.width * canvas.height; i++) {
        canvas.frameBuffer[i] = color;
    }
}

static inline void zgl_render_line(int x0, int x1, int y0, int y1, uint32_t color, zgl_canvas_t canvas) {
    int delta_x = (x1 - x0);
    int delta_y = (y1 - y0);
    int longest_side_length = (abs(delta_x) >= abs(delta_y)) ? abs(delta_x) : abs(delta_y);
    float x_inc = delta_x / (float)longest_side_length;
    float y_inc = delta_y / (float)longest_side_length;
    float current_x = x0;
    float current_y = y0;
    for (int i = 0; i <= longest_side_length; i++) {
        zgl_render_pixel(round(current_x), round(current_y), 0.0, color, canvas);
        current_x += x_inc;
        current_y += y_inc;
    }
}

static inline void zgl_render_circle(int x, int y, int r, uint32_t color, zgl_canvas_t canvas) {
    int x1 = x - r;
    int x2 = x + r;
    int y1 = y - r;
    int y2 = y + r;
    for (int j = y1; j < y2; j++) {
        for (int i = x1; i < x2; i++) {
            int dx = i - x;
            int dy = j - y;
            if (dx * dx + dy * dy <= r * r) {
                zgl_render_pixel(i, j, 0.0, color, canvas);
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
                                           zgl_canvas_t canvas) {
    int x_min = ZGL__MAX(ZGL__MIN(ZGL__MIN(x0, x1), x2), 0);
    int x_max = ZGL__MIN(ZGL__MAX(ZGL__MAX(x0, x1), x2), canvas.width - 1);
    int y_min = ZGL__MAX(ZGL__MIN(ZGL__MIN(y0, y1), y2), 0);
    int y_max = ZGL__MIN(ZGL__MAX(ZGL__MAX(y0, y1), y2), canvas.height - 1);

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
                if (z < canvas.depthBuffer[y * canvas.width + x]) {
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
                    zgl_render_pixel(x, y, z, color, canvas); // TODO: Avoid scissor test in zgl_render_pixel
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

static inline void zgl_render_object3D(zgl_object3D_t* object, void *uniformData, zgl_camera_t camera, zgl_canvas_t canvas,
                                       zgl_vertex_shader_t vertexShader, zgl_fragment_shader_t fragmentShader, uint16_t renderOptions) {
    zgl_mesh_t* mesh = object->mesh;

    // Don't draw if the object is fully outside the camera fustrum
    if (renderOptions & ZGL_FUSTRUM_CULLING) {
        for (int p = 0; p < 6; p++) {
            zgl_vec4_t plane = camera.fustrumPlanes[p];
            int isInside = 1;
            zgl_vec4_t center = zgl_mul_mat_v4(object->transform, (zgl_vec4_t) {mesh->center.x, mesh->center.y, mesh->center.z, 1});
            if (zgl_dot_v4(plane, center) < -(object->scale * mesh->boundsRadius)) {
                ZGL_DEBUG_PRINT("DEBUG: Culled object using fustrum culling {plane %d}\n", p);
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
            xs[v] = (vertexShaderOutput[v].position.x + 1.0f) * canvas.width / 2.0f;
            ys[v] = (1.0f - vertexShaderOutput[v].position.y) * canvas.height / 2.0f;
            zs[v] = vertexShaderOutput[v].position.z; // For z-buffer
            invws[v] = invw; // Store 1/w to avoid divisions later when performing perspective correct interpolation
        }

        int area = zgl__edge_cross(xs[0], ys[0], xs[1], ys[1], xs[2], ys[2]);

        // Backface culling
        if ((renderOptions & ZGL_BACKFACE_CULLING) && area <= 0) {
            ZGL_DEBUG_PRINT("DEBUG: Culled triangle using backface culling\n");
            continue;
        }

        if ((renderOptions & ZGL_FUSTRUM_CULLING) && !zgl__tri_in_fustrum(vertexShaderOutput[0].position,
                                                                          vertexShaderOutput[1].position,
                                                                          vertexShaderOutput[2].position)) {
            ZGL_DEBUG_PRINT("DEBUG: Culled triangle using fustrum culling\n");
            continue;
        }

        // Rasterization
        zgl__rasterize_triangle(xs[0], xs[1], xs[2],
                                ys[0], ys[1], ys[2],
                                zs[0], zs[1], zs[2],
                                invws[0], invws[1], invws[2],
                                area, vertexShaderOutput, fragmentShader, uniformData, canvas);
    }
}

static inline void zgl_render_triangle(int x0, int y0, uint32_t color0,
                                       int x1, int y1, uint32_t color1,
                                       int x2, int y2, uint32_t color2,
                                       zgl_canvas_t canvas) {
    int area = zgl__edge_cross(x0, y0, x1, y1, x2, y2);
    float r, g, b;

    zgl_shader_context_t shaderContext0 = {0};
    shaderContext0.position = (zgl_vec4_t) {x0, y0, 0, 1};
    shaderContext0.numAttributes = 3;
    zgl_color_to_floats(color0, &r, &g, &b);
    shaderContext0.attributes[0] = r;
    shaderContext0.attributes[1] = g;
    shaderContext0.attributes[2] = b;

    zgl_shader_context_t shaderContext1 = {0};
    shaderContext1.position = (zgl_vec4_t) {x1, y1, 0, 1};
    shaderContext1.numAttributes = 3;
    zgl_color_to_floats(color1, &r, &g, &b);
    shaderContext1.attributes[0] = r;
    shaderContext1.attributes[1] = g;
    shaderContext1.attributes[2] = b;

    zgl_shader_context_t shaderContext2 = {0};
    shaderContext2.position = (zgl_vec4_t) {x2, y2, 0, 1};
    shaderContext2.numAttributes = 3;
    zgl_color_to_floats(color2, &r, &g, &b);
    shaderContext2.attributes[0] = r;
    shaderContext2.attributes[1] = g;
    shaderContext2.attributes[2] = b;

    zgl_shader_context_t shaderContexts[3] = {shaderContext0, shaderContext1, shaderContext2};
    zgl__rasterize_triangle(x0, x1, x2, y0, y1, y2, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                            area, shaderContexts, zgl_colored_fragment_shader, NULL, canvas);
}

/* Shaders */

/* Basic shading */
// Draw with a single color, no lighting or textures

static inline zgl_shader_context_t zgl_basic_vertex_shader(void* inputVertex, void* uniformData) {
    zgl_vertex_input_t* inputVertexData = (zgl_vertex_input_t*) inputVertex;
    zgl_shader_context_t result = {0};
    zgl_basic_uniform_t* basicUniformData = (zgl_basic_uniform_t*) uniformData;
    zgl_vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};
    result.position = zgl_mul_mat_v4(basicUniformData->modelviewprojection, inputVertex4);
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
    result.position = zgl_mul_mat_v4(basicUniformData->modelviewprojection, inputVertex4);
    result.numAttributes = 9;
    zgl_color_to_floats(basicUniformData->materials[inputVertexData->materialIndex].ambientColor, &result.attributes[0], &result.attributes[1], &result.attributes[2]);
    zgl_color_to_floats(basicUniformData->materials[inputVertexData->materialIndex].diffuseColor, &result.attributes[3], &result.attributes[4], &result.attributes[5]);
    zgl_color_to_floats(basicUniformData->materials[inputVertexData->materialIndex].specularColor, &result.attributes[6], &result.attributes[7], &result.attributes[8]);

    return result;
}

static inline uint32_t zgl_colored_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    uint32_t ambientColor = zgl_color_from_floats(input->attributes[0], input->attributes[1], input->attributes[2]);
    uint32_t diffuseColor = zgl_color_from_floats(input->attributes[3], input->attributes[4], input->attributes[5]);
    uint32_t specularColor = zgl_color_from_floats(input->attributes[6], input->attributes[7], input->attributes[8]);

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


    result.numAttributes = 2;
    result.attributes[0] = inputVertexData->textureCoord.x;   // u
    result.attributes[1] = inputVertexData->textureCoord.y;   // v

    // Set other vertex attributes
    result.numFlatAttributes = 14;

    zgl_vec3_t worldSpaceNormal = zgl_mul_mat_v3(defaultUniformData->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    float invMagNormal = 1.0f / zgl_magnitude(worldSpaceNormal);
    result.flatAttributes[0] = worldSpaceNormal.x;
    result.flatAttributes[1] = worldSpaceNormal.y;
    result.flatAttributes[2] = worldSpaceNormal.z;
    zgl_color_to_floats(defaultUniformData->materials[inputVertexData->materialIndex].ambientColor, &result.flatAttributes[3], &result.flatAttributes[4], &result.flatAttributes[5]);
    zgl_color_to_floats(defaultUniformData->materials[inputVertexData->materialIndex].diffuseColor, &result.flatAttributes[6], &result.flatAttributes[7], &result.flatAttributes[8]);
    zgl_color_to_floats(defaultUniformData->materials[inputVertexData->materialIndex].specularColor, &result.flatAttributes[9], &result.flatAttributes[10], &result.flatAttributes[11]);
    result.flatAttributes[12] = defaultUniformData->materials[inputVertexData->materialIndex].specularExponent;

    result.flatAttributes[13] = inputVertexData->materialIndex;
    return result;
}

static inline uint32_t zgl_unlit_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    zgl_unlit_uniform_t* uniform = (zgl_unlit_uniform_t*) uniformData;
    int materialIndex = input->flatAttributes[13];

    zgl_canvas_t ambientTexture = uniform->materials[materialIndex].ambientTexture;
    uint32_t ambientColor = zgl_color_from_floats(input->flatAttributes[3], input->flatAttributes[4], input->flatAttributes[5]);
    if (ambientTexture.width != 0 && ambientTexture.height != 0) {
        float u = input->attributes[0];
        float v = input->attributes[1];
        ambientColor = zgl_mul_colors(ambientColor, zgl_sample_texture(u, v, ambientTexture, uniform->bilinearFiltering));
    }

    zgl_canvas_t diffuseTexture = uniform->materials[materialIndex].diffuseTexture;
    uint32_t diffuseColor = zgl_color_from_floats(input->flatAttributes[6], input->flatAttributes[7], input->flatAttributes[8]);
    if (diffuseTexture.width != 0 && diffuseTexture.height != 0) {
        float u = input->attributes[0];
        float v = input->attributes[1];
        diffuseColor = zgl_mul_colors(diffuseColor, zgl_sample_texture(u, v, diffuseTexture, uniform->bilinearFiltering));
    }

    zgl_canvas_t specularTexture = uniform->materials[materialIndex].specularTexture;
    uint32_t specularColor = zgl_color_from_floats(input->flatAttributes[9], input->flatAttributes[10], input->flatAttributes[11]);
    if (specularTexture.width != 0 && specularTexture.height != 0) {
        float u = input->attributes[0];
        float v = input->attributes[1];
        specularColor = zgl_mul_colors(specularColor, zgl_sample_texture(u, v, specularTexture, uniform->bilinearFiltering));
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

    result.numAttributes = 2;
    result.attributes[0] = inputVertexData->textureCoord.x;   // u
    result.attributes[1] = inputVertexData->textureCoord.y;   // v

    // Set other vertex attributes
    result.numFlatAttributes = 23;

    zgl_vec3_t worldSpaceNormal = zgl_mul_mat_v3(defaultUniformData->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    float invMagNormal = 1.0f / zgl_magnitude(worldSpaceNormal);
    result.flatAttributes[0] = worldSpaceNormal.x;
    result.flatAttributes[1] = worldSpaceNormal.y;
    result.flatAttributes[2] = worldSpaceNormal.z;
    zgl_color_to_floats(defaultUniformData->materials[inputVertexData->materialIndex].ambientColor, &result.flatAttributes[3], &result.flatAttributes[4], &result.flatAttributes[5]);
    zgl_color_to_floats(defaultUniformData->materials[inputVertexData->materialIndex].diffuseColor, &result.flatAttributes[6], &result.flatAttributes[7], &result.flatAttributes[8]);
    zgl_color_to_floats(defaultUniformData->materials[inputVertexData->materialIndex].specularColor, &result.flatAttributes[9], &result.flatAttributes[10], &result.flatAttributes[11]);
    result.flatAttributes[12] = defaultUniformData->materials[inputVertexData->materialIndex].specularExponent;

    zgl_lighting_result_t lightResult = zgl_lighting((zgl_vec3_t) {worldSpaceVertex.x, worldSpaceVertex.y, worldSpaceVertex.z},
                                                     worldSpaceNormal, invMagNormal,
                                                     defaultUniformData->materials[inputVertexData->materialIndex].specularExponent,
                                                     defaultUniformData->lightSources, ZGL_DIFFUSE_LIGHTING | ZGL_SPECULAR_LIGHTING);
    result.flatAttributes[13] = lightResult.ambient.x;
    result.flatAttributes[14] = lightResult.ambient.y;
    result.flatAttributes[15] = lightResult.ambient.z;
    result.flatAttributes[16] = lightResult.diffuse.x;
    result.flatAttributes[17] = lightResult.diffuse.y;
    result.flatAttributes[18] = lightResult.diffuse.z;
    result.flatAttributes[19] = lightResult.specular.x;
    result.flatAttributes[20] = lightResult.specular.y;
    result.flatAttributes[21] = lightResult.specular.z;

    result.flatAttributes[22] = inputVertexData->materialIndex;
    return result;
}

static inline uint32_t zgl_flat_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    zgl_flat_uniform_t* uniform = (zgl_flat_uniform_t*) uniformData;
    int materialIndex = input->flatAttributes[19];

    zgl_canvas_t ambientTexture = uniform->materials[materialIndex].ambientTexture;
    uint32_t ambientColor = zgl_color_from_floats(input->flatAttributes[3], input->flatAttributes[4], input->flatAttributes[5]);
    if (ambientTexture.width != 0 && ambientTexture.height != 0) {
        float u = input->attributes[0];
        float v = input->attributes[1];
        ambientColor = zgl_mul_colors(ambientColor, zgl_sample_texture(u, v, ambientTexture, uniform->bilinearFiltering));
    }
    zgl_vec3_t ambientLight = {input->flatAttributes[13], input->flatAttributes[14], input->flatAttributes[15]};
    ambientColor = zgl_mul_vec3_color(ambientLight, ambientColor);

    zgl_canvas_t diffuseTexture = uniform->materials[materialIndex].diffuseTexture;
    uint32_t diffuseColor = zgl_color_from_floats(input->flatAttributes[6], input->flatAttributes[7], input->flatAttributes[8]);
    if (diffuseTexture.width != 0 && diffuseTexture.height != 0) {
        float u = input->attributes[0];
        float v = input->attributes[1];
        diffuseColor = zgl_mul_colors(diffuseColor, zgl_sample_texture(u, v, diffuseTexture, uniform->bilinearFiltering));
    }
    zgl_vec3_t diffuseLight = {input->flatAttributes[16], input->flatAttributes[17], input->flatAttributes[18]};
    diffuseColor = zgl_mul_vec3_color(diffuseLight, diffuseColor);

    zgl_canvas_t specularTexture = uniform->materials[materialIndex].specularTexture;
    uint32_t specularColor = zgl_color_from_floats(input->flatAttributes[9], input->flatAttributes[10], input->flatAttributes[11]);
    if (specularTexture.width != 0 && specularTexture.height != 0) {
        float u = input->attributes[0];
        float v = input->attributes[1];
        specularColor = zgl_mul_colors(specularColor, zgl_sample_texture(u, v, specularTexture, uniform->bilinearFiltering));
    }
    zgl_vec3_t specularLight = {input->flatAttributes[19], input->flatAttributes[20], input->flatAttributes[21]};
    specularColor = zgl_mul_vec3_color(specularLight, specularColor);

    return zgl_add_colors(ambientColor, zgl_add_colors(diffuseColor, specularColor));
}


/* Gouraud shading */
// Compute the lighting at each vertex and interpolate the values at each fragment

static inline zgl_shader_context_t zgl_gourard_vertex_shader(void* inputVertex, void* uniformData) {
    zgl_vertex_input_t* inputVertexData = (zgl_vertex_input_t*) inputVertex;
    zgl_shader_context_t result = {0};
    zgl_gourard_uniform_t* defaultUniformData = (zgl_gourard_uniform_t*) uniformData;
    zgl_vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};
    zgl_vec4_t worldSpaceVertex = zgl_mul_mat_v4(defaultUniformData->modelMatrix, inputVertex4); // Local to world space
    result.position = zgl_mul_mat_v4(defaultUniformData->viewProjectionMatrix, worldSpaceVertex); // World to clip space

    // Set other vertex attributes
    result.numAttributes = 24;

    zgl_vec3_t worldSpaceNormal = zgl_mul_mat_v3(defaultUniformData->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    float invMagNormal = 1.0f / zgl_magnitude(worldSpaceNormal);
    result.attributes[0] = worldSpaceNormal.x;
    result.attributes[1] = worldSpaceNormal.y;
    result.attributes[2] = worldSpaceNormal.z;
    result.attributes[3] = inputVertexData->textureCoord.x;   // u
    result.attributes[4] = inputVertexData->textureCoord.y;   // v
    zgl_color_to_floats(defaultUniformData->materials[inputVertexData->materialIndex].ambientColor, &result.attributes[5], &result.attributes[6], &result.attributes[7]);
    zgl_color_to_floats(defaultUniformData->materials[inputVertexData->materialIndex].diffuseColor, &result.attributes[8], &result.attributes[9], &result.attributes[10]);
    zgl_color_to_floats(defaultUniformData->materials[inputVertexData->materialIndex].specularColor, &result.attributes[11], &result.attributes[12], &result.attributes[13]);
    result.attributes[14] = defaultUniformData->materials[inputVertexData->materialIndex].specularExponent;

    zgl_lighting_result_t lightResult = zgl_lighting((zgl_vec3_t) {worldSpaceVertex.x, worldSpaceVertex.y, worldSpaceVertex.z},
                                                     worldSpaceNormal, invMagNormal,
                                                     defaultUniformData->materials[inputVertexData->materialIndex].specularExponent,
                                                     defaultUniformData->lightSources, ZGL_DIFFUSE_LIGHTING | ZGL_SPECULAR_LIGHTING);
    result.attributes[15] = lightResult.ambient.x;
    result.attributes[16] = lightResult.ambient.y;
    result.attributes[17] = lightResult.ambient.z;
    result.attributes[18] = lightResult.diffuse.x;
    result.attributes[19] = lightResult.diffuse.y;
    result.attributes[20] = lightResult.diffuse.z;
    result.attributes[21] = lightResult.specular.x;
    result.attributes[22] = lightResult.specular.y;
    result.attributes[23] = lightResult.specular.z;

    result.numFlatAttributes = 1;
    result.flatAttributes[0] = inputVertexData->materialIndex;
    return result;
}


static inline uint32_t zgl_gourard_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    zgl_gourard_uniform_t* uniform = (zgl_gourard_uniform_t*) uniformData;
    int materialIndex = input->flatAttributes[0];

    zgl_canvas_t ambientTexture = uniform->materials[materialIndex].ambientTexture;
    uint32_t ambientColor = zgl_color_from_floats(input->attributes[5], input->attributes[6], input->attributes[7]);
    if (ambientTexture.width != 0 && ambientTexture.height != 0) {
        float u = input->attributes[3];
        float v = input->attributes[4];
        ambientColor = zgl_mul_colors(ambientColor, zgl_sample_texture(u, v, ambientTexture, uniform->bilinearFiltering));
    }
    zgl_vec3_t ambientLight = {input->attributes[15], input->attributes[16], input->attributes[17]};
    ambientColor = zgl_mul_vec3_color(ambientLight, ambientColor);

    zgl_canvas_t diffuseTexture = uniform->materials[materialIndex].diffuseTexture;
    uint32_t diffuseColor = zgl_color_from_floats(input->attributes[8], input->attributes[9], input->attributes[10]);
    if (diffuseTexture.width != 0 && diffuseTexture.height != 0) {
        float u = input->attributes[3];
        float v = input->attributes[4];
        diffuseColor = zgl_mul_colors(diffuseColor, zgl_sample_texture(u, v, diffuseTexture, uniform->bilinearFiltering));
    }
    zgl_vec3_t diffuseLight = {input->attributes[18], input->attributes[19], input->attributes[20]};
    diffuseColor = zgl_mul_vec3_color(diffuseLight, diffuseColor);

    zgl_canvas_t specularTexture = uniform->materials[materialIndex].specularTexture;
    uint32_t specularColor = zgl_color_from_floats(input->attributes[11], input->attributes[12], input->attributes[13]);
    if (specularTexture.width != 0 && specularTexture.height != 0) {
        float u = input->attributes[3];
        float v = input->attributes[4];
        specularColor = zgl_mul_colors(specularColor, zgl_sample_texture(u, v, specularTexture, uniform->bilinearFiltering));
    }
    zgl_vec3_t specularLight = {input->attributes[21], input->attributes[22], input->attributes[23]};
    specularColor = zgl_mul_vec3_color(specularLight, specularColor);

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

    // Set other vertex attributes
    result.numAttributes = 18;
    zgl_vec3_t worldSpaceNormal = zgl_mul_mat_v3(uniform->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    result.attributes[0] = worldSpaceNormal.x;
    result.attributes[1] = worldSpaceNormal.y;
    result.attributes[2] = worldSpaceNormal.z;
    result.attributes[3] = worldSpaceVertex.x;
    result.attributes[4] = worldSpaceVertex.y;
    result.attributes[5] = worldSpaceVertex.z;
    result.attributes[6] = inputVertexData->textureCoord.x;    // u
    result.attributes[7] = inputVertexData->textureCoord.y;    // v

    zgl_color_to_floats(uniform->materials[inputVertexData->materialIndex].ambientColor, &result.attributes[8], &result.attributes[9], &result.attributes[10]);
    zgl_color_to_floats(uniform->materials[inputVertexData->materialIndex].diffuseColor, &result.attributes[11], &result.attributes[12], &result.attributes[13]);
    zgl_color_to_floats(uniform->materials[inputVertexData->materialIndex].specularColor, &result.attributes[14], &result.attributes[15], &result.attributes[16]);
    result.attributes[17] = uniform->materials[inputVertexData->materialIndex].specularExponent;

    result.numFlatAttributes = 1;
    result.flatAttributes[0] = inputVertexData->materialIndex;
    return result;
}

static inline uint32_t zgl_phong_fragment_shader(const zgl_shader_context_t* input, void* uniformData) {
    zgl_phong_uniform_t* uniform = (zgl_phong_uniform_t*) uniformData;

    zgl_vec3_t normal = {input->attributes[0], input->attributes[1], input->attributes[2]};
    zgl_vec3_t position = {input->attributes[3], input->attributes[4], input->attributes[5]};
    float specularExponent = input->attributes[17];
    float invMagNormal = 1.0f / zgl_magnitude(normal);
    zgl_lighting_result_t lightingResult = zgl_lighting(position, normal, invMagNormal, specularExponent, uniform->lightSources, ZGL_DIFFUSE_LIGHTING | ZGL_SPECULAR_LIGHTING);
    int materialIndex = input->flatAttributes[0];
    
    float u = input->attributes[6];
    float v = input->attributes[7];

    zgl_canvas_t ambientTexture = uniform->materials[materialIndex].ambientTexture;
    uint32_t ambientColor = zgl_color_from_floats(input->attributes[8], input->attributes[9], input->attributes[10]);
    if (ambientTexture.width != 0 && ambientTexture.height != 0) {    
        ambientColor = zgl_mul_colors(ambientColor, zgl_sample_texture(u, v, ambientTexture, uniform->bilinearFiltering));
    }
    ambientColor = zgl_mul_vec3_color(lightingResult.ambient, ambientColor);
    
    zgl_canvas_t diffuseTexture = uniform->materials[materialIndex].diffuseTexture;
    uint32_t diffuseColor = zgl_color_from_floats(input->attributes[11], input->attributes[12], input->attributes[13]);
    if (diffuseTexture.width != 0 && diffuseTexture.height != 0) {    
        diffuseColor = zgl_mul_colors(diffuseColor, zgl_sample_texture(u, v, diffuseTexture, uniform->bilinearFiltering));
    }
    diffuseColor = zgl_mul_vec3_color(lightingResult.diffuse, diffuseColor);

    zgl_canvas_t specularTexture = uniform->materials[materialIndex].specularTexture;
    uint32_t specularColor = zgl_color_from_floats(input->attributes[14], input->attributes[15], input->attributes[16]);
    if (specularTexture.width != 0 && specularTexture.height != 0) {
        specularColor = zgl_mul_colors(specularColor, zgl_sample_texture(u, v, specularTexture, uniform->bilinearFiltering));
    }
    specularColor = zgl_mul_vec3_color(lightingResult.specular, specularColor);

    return zgl_add_colors(specularColor, zgl_add_colors(ambientColor, diffuseColor));
}

#endif // ZEROGL_IMPLEMENTATION_INCLUDED
#endif // ZEROGL_IMPLEMENTATION
