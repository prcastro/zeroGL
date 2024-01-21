#ifndef ZEROGL_H
#define ZEROGL_H

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#ifdef ZGL_DEBUG
#define DEBUG_PRINT(...) printf(__VA_ARGS__)
#else
#define DEBUG_PRINT(...) do {} while (0)
#endif

/* CONFIGURATION */

#ifndef MAX_VERTEX_ATTRIBUTES
#define MAX_VERTEX_ATTRIBUTES 50
#endif // MAX_VERTEX_ATTRIBUTES

// Rendering options
#define DIFFUSE_LIGHTING (1 << 0)
#define SPECULAR_LIGHTING (1 << 1)
#define BACKFACE_CULLING (1 << 2)
#define FUSTRUM_CULLING (1 << 3)
#define BILINEAR_FILTERING (1 << 4) // TODO: Implement bilinear filtering
#define FLAT_SHADING (1 << 5)

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/* VECTORS AND MATRICES */

#define M_PI 3.14159265358979323846264338327950288

typedef struct {
  float x, y, z;
} vec3_t;

typedef struct {
  float x, y, z, w;
} vec4_t;

typedef struct {
  float data[4][4];
} mat4x4_t;

static const mat4x4_t IDENTITY_M4x4 = {{
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
}};

static inline vec3_t crossProduct(vec3_t a, vec3_t b) {
    vec3_t result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

static inline float dot(vec3_t a, vec3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline float dotV4(vec4_t a, vec4_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

static inline float magnitude(vec3_t v) {
    float result = sqrt(dot(v, v));
    assert(result >= 0);
    return result;
}

static inline vec3_t sub(vec3_t a, vec3_t b) {
    return (vec3_t) {a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline vec4_t subV4(vec4_t a, vec4_t b) {
    return (vec4_t) {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}

static inline vec3_t add(vec3_t a, vec3_t b) {
    return (vec3_t) {a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline vec4_t addV4(vec4_t a, vec4_t b) {
    return (vec4_t) {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

static inline vec3_t normalize(vec3_t v) {
    float mag = magnitude(v);
    return (vec3_t) {v.x / mag, v.y / mag, v.z / mag};
}

static inline vec3_t mulScalarV3(float k, vec3_t v) {
    return (vec3_t) {k*v.x, k*v.y, k*v.z};
}

static inline vec4_t mulMV4(mat4x4_t mat4x4, vec4_t vec4) {
  float result[4] = {0};
  for (int i = 0; i < 4; i++) {
    result[i] += mat4x4.data[i][0]*vec4.x;
    result[i] += mat4x4.data[i][1]*vec4.y;
    result[i] += mat4x4.data[i][2]*vec4.z;
    result[i] += mat4x4.data[i][3]*vec4.w;
  }
  return (vec4_t) {result[0], result[1], result[2], result[3]};
}

static inline vec3_t mulMV3(mat4x4_t mat4x4, vec3_t v) {
    vec4_t vec4 = {v.x, v.y, v.z, 1};
    vec4_t result = mulMV4(mat4x4, vec4);
    return (vec3_t) {result.x, result.y, result.z};
}

static inline mat4x4_t mulMM4(mat4x4_t m1, mat4x4_t m2) {
    mat4x4_t result = {0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                result.data[i][j] += m1.data[i][k]*m2.data[k][j];
            }
        }
    }
    return result;
}

static inline mat4x4_t transposeM4(mat4x4_t m) {
    mat4x4_t result = {0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result.data[i][j] = m.data[j][i];
        }
    }
    return result;
}

static inline float determinant(float a, float b, float c, float d, float e, float f, float g, float h, float i) {
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}

static inline mat4x4_t inverseM4(mat4x4_t matrix) {
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
    
    mat4x4_t invMatrix;
    for (int i = 0; i < 4; i++) {
        invMatrix.data[i][0] = invOut[i*4];
        invMatrix.data[i][1] = invOut[i*4+1];
        invMatrix.data[i][2] = invOut[i*4+2];
        invMatrix.data[i][3] = invOut[i*4+3];
    }

    return invMatrix;
}

static inline mat4x4_t translationToMatrix(vec3_t vector) {
    return (mat4x4_t) {{
        {1, 0, 0, vector.x},
        {0, 1, 0, vector.y},
        {0, 0, 1, vector.z},
        {0, 0, 0,        1}
    }};
}

static inline mat4x4_t scaleToMatrix(float scale) {
    return (mat4x4_t) {{
        {scale, 0,     0,     0},
        {0,     scale, 0,     0},
        {0,     0,     scale, 0},
        {0,     0,     0,     1}
    }};
}

// TODO: Only use quaternion for rotation
static inline mat4x4_t rotationX(float degrees) {
    float radians = degrees * M_PI / 180.0f;
    float cos = cosf(radians);
    float sin = sinf(radians);
    return (mat4x4_t) {{
        { 1, 0,    0,   0 },
        { 0, cos,  sin, 0 },
        { 0, -sin, cos, 0 },
        { 0, 0,    0,   1 }
    }};
}

static inline mat4x4_t rotationY(float degrees) {
    float radians = degrees * M_PI / 180.0f;
    float cos = cosf(radians);
    float sin = sinf(radians);
    return (mat4x4_t) {{
        { cos, 0, -sin, 0 },
        { 0,   1, 0,    0 },
        { sin, 0, cos,  0 },
        { 0,   0, 0,    1 }
    }};
}

static inline mat4x4_t rotationZ(float degrees) {
    float radians = degrees * M_PI / 180.0f;
    float cos = cosf(radians);
    float sin = sinf(radians);
    return (mat4x4_t) {{
        { cos,  sin, 0, 0 },
        { -sin, cos, 0, 0 },
        { 0,    0,   1, 0 },
        { 0,    0,   0, 1 }
    }};
}

/* QUATERNIONS */

typedef struct {
    float w, x, y, z;
} quaternion_t;

quaternion_t quaternionMul(quaternion_t q1, quaternion_t q2) {
    quaternion_t result;
    result.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    result.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    result.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    result.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return result;
}

quaternion_t quaternionFromAngleAxis(float degrees, vec3_t axis) {
    vec3_t normalizedAxis = normalize(axis);
    float angle = degrees * M_PI / 180.0f;
    float halfAngle = angle * 0.5f;
    float sinHalfAngle = sinf(halfAngle);
    return (quaternion_t) {cosf(halfAngle), normalizedAxis.x * sinHalfAngle, normalizedAxis.y * sinHalfAngle, normalizedAxis.z * sinHalfAngle};
}

vec3_t rotateVectorByQuaternion(vec3_t v, quaternion_t q) {
    quaternion_t q_v = {0, v.x, v.y, v.z};
    quaternion_t q_conjugate = {q.w, -q.x, -q.y, -q.z};
    quaternion_t rotated = quaternionMul(quaternionMul(q, q_v), q_conjugate);
    return (vec3_t){rotated.x, rotated.y, rotated.z};
}

/* COLORS */

static const uint32_t COLOR_WHITE  = 0x00FFFFFF;
static const uint32_t COLOR_BLACK  = 0x00000000;
static const uint32_t COLOR_RED    = 0x00FF0000;
static const uint32_t COLOR_GREEN  = 0x0000FF00;
static const uint32_t COLOR_BLUE   = 0x000000FF;
static const uint32_t COLOR_YELLOW = 0x00FFFF00;
static const uint32_t COLOR_PURPLE = 0x00FF00FF;
static const uint32_t COLOR_CYAN   = 0x0000FFFF;

static inline uint32_t colorToUint32(uint8_t r, uint8_t g, uint8_t b) {
    return 0x00000000 | (r << 16) | (g << 8) |  b;    
}

static inline void colorFromUint32(uint32_t c, uint8_t* r, uint8_t* g, uint8_t* b) {
    *r = (c & 0x00FF0000) >> 16;
    *g = (c & 0x0000FF00) >> 8;
    *b = c & 0x000000FF;
}

static inline float clamp(float v, float max) {
    return v > max ? max : v;
}

static inline uint32_t mulScalarColor(double x, uint32_t color) {
    uint8_t r, g, b;
    colorFromUint32(color, &r, &g, &b);
    return colorToUint32(clamp(x * r, 255.0), clamp(x * g, 255.0), clamp(x * b, 255.0));
}

static inline uint32_t sumColors(uint32_t c0, uint32_t c1) {
    uint8_t r0, g0, b0;
    uint8_t r1, g1, b1;
    colorFromUint32(c0, &r0, &g0, &b0);
    colorFromUint32(c1, &r1, &g1, &b1);
    return colorToUint32(clamp(r0 + r1, 255.0), clamp(g0 + g1, 255.0), clamp(b0 + b1, 255.0));
}

static inline uint32_t colorFromFloats(float r, float g, float b) {
    return colorToUint32(clamp(r * 255.0, 255.0), clamp(g * 255.0, 255.0), clamp(b * 255.0, 255.0));
}

static inline void colorToFloats(uint32_t color, float* r, float* g, float* b) {
    uint8_t r8, g8, b8;
    colorFromUint32(color, &r8, &g8, &b8);
    *r = r8 / 255.0f;
    *g = g8 / 255.0f;
    *b = b8 / 255.0f;
}

/* 3D OBJECTS */

typedef struct {
  int     v0, v1, v2;
  int     t0, t1, t2;
  int     n0, n1, n2;
  int     materialIndex;
} triangle_t;

typedef struct {
    char*   name;
    uint32_t diffuseColor;
    uint32_t specularColor;
    float   specularExponent;
    int     textureWidth;
    int     textureHeight;
    uint32_t* texture;
} material_t;

typedef struct {
    char*       name;
    int         numVertices;
    vec3_t*     vertices;
    int         numTextureCoords;
    vec3_t*     textureCoords;
    int         numNormals;
    vec3_t*     normals;
    float*      invMagnitudeNormals;
    int         numTriangles;
    triangle_t* triangles;
    int         numMaterials;
    material_t* materials;
    vec3_t      center;
    float       boundsRadius;
} mesh_t;

typedef struct {
    mesh_t*  mesh;
    vec3_t   translation;
    float    scale;
    mat4x4_t rotation;
    mat4x4_t transform;
} object3D_t;

static inline object3D_t makeObject(mesh_t *mesh, vec3_t translation, float scale, mat4x4_t rotation) {
    mat4x4_t translationMatrix = translationToMatrix(translation);
    mat4x4_t scaleMatrix = scaleToMatrix(scale);
    mat4x4_t transform = mulMM4(translationMatrix, mulMM4(rotation, scaleMatrix));
    return (object3D_t) {mesh, translation, scale, rotation, transform};
}

static inline vec3_t meshCenter(vec3_t* vertices, int numVertices) {
    vec3_t result = {0, 0, 0};
    for (int i = 0; i < numVertices; i++) {
        result = add(result, vertices[i]);
    }
    return mulScalarV3(1.0f / numVertices, result);
}

static inline float meshBoundsRadius(vec3_t* vertices, int numVertices, vec3_t center) {
    float result = 0.0f;
    for (int i = 0; i < numVertices; i++) {
        float distance = magnitude(sub(vertices[i], center));
        if (distance > result) {
            result = distance;
        }
    }
    return result;
}

/* LIGHTING */

typedef struct {
    float intensity;
} ambient_light_t;

typedef struct {
    float intensity;
    vec3_t direction;
} dir_light_t;

typedef struct {
    float intensity;
    vec3_t position;
} point_light_t;

typedef struct {
    int              numAmbientLights;
    ambient_light_t* ambientLights;
    int              numDirectionalLights;
    dir_light_t*     directionalLights;
    int              numPointLights;
    point_light_t*   pointLights;
} light_sources_t;

// TODO: Code here is a bit repeated between directional and point lights. Maybe refactor?
float computeLighting(vec3_t position, vec3_t normal, float invMagnitudeNormal, float specularExponent,
                      light_sources_t lightSources, uint8_t renderOptions) {
    int numDirectionalLights = lightSources.numDirectionalLights;
    dir_light_t* directionalLights = lightSources.directionalLights;
    int numPointLights = lightSources.numPointLights;
    point_light_t* pointLights = lightSources.pointLights;
    int numAmbientLights = lightSources.numAmbientLights;
    ambient_light_t* ambientLights = lightSources.ambientLights;

    float diffuseIntensity  = 0.0;
    float specularIntensity = 0.0;
    float ambientIntensity  = 0.0;
    
    // Directional lights
    for (int i = 0; i < numDirectionalLights; i++) {
        vec3_t lightDirection = directionalLights[i].direction;
        float magnitudeLightDirection = magnitude(lightDirection);
        float invMagnitudeLightDirection = 1.0f / magnitudeLightDirection;
        if (renderOptions & DIFFUSE_LIGHTING) {
            float cos_alpha = -dot(lightDirection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            diffuseIntensity += MAX(cos_alpha, 0.0f) * directionalLights[i].intensity;
        }

        if (renderOptions & SPECULAR_LIGHTING) {
            vec3_t reflection = sub(mulScalarV3(2 * -dot(lightDirection, normal), normal), lightDirection);
            float cos_beta = -dot(reflection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            specularIntensity += pow(MAX(cos_beta, 0.0f), specularExponent) * directionalLights[i].intensity;
        }
    }

    // Point lights
    for (int i = 0; i < numPointLights; i++) {
        vec3_t lightDirection = sub(pointLights[i].position, position);
        float invMagnitudeLightDirection = 1.0f / magnitude(lightDirection);
        if (renderOptions & DIFFUSE_LIGHTING) {
            float cos_alpha = dot(lightDirection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            diffuseIntensity += MAX(cos_alpha, 0) * pointLights[i].intensity;
        }

        if (renderOptions & SPECULAR_LIGHTING) {
            vec3_t reflection = sub(mulScalarV3(2 * dot(lightDirection, normal), normal), lightDirection);
            float cos_beta = dot(reflection, normal) * invMagnitudeLightDirection * invMagnitudeNormal;
            specularIntensity += pow(MAX(cos_beta, 0), specularExponent) * pointLights[i].intensity;
        }
    }

    // Ambient light
    for (int i = 0; i < numAmbientLights; i++) {
        ambientIntensity += ambientLights[i].intensity;
    }

    return (diffuseIntensity + specularIntensity + ambientIntensity);
}

/* 3D CAMERA */

// Camera for left handed coordinate system
typedef struct {
    vec3_t position;            // x, y, z coordinates of the camera
    vec3_t direction;           // Direction in which the camera is looking
    vec3_t up;                  // Up direction for the camera (usually [0, 1, 0])
    float fov;                  // Field of view (in degrees)
    float aspectRatio;          // Aspect ratio of the viewport
    float nearPlane;            // Distance to the near clipping plane
    float farPlane;             // Distance to the far clipping plane
    mat4x4_t viewMatrix;        // View matrix
    mat4x4_t projectionMatrix;  // Projection matrix
    mat4x4_t viewProjMatrix;    // View * Projection matrix
    vec4_t fustrumPlanes[6];    // View frustum planes
    float movementSpeed;
    float turningSpeed;
} camera_t;

vec4_t normalizePlane(vec4_t plane) {
    float mag = sqrt(plane.x*plane.x + plane.y*plane.y + plane.z*plane.z);
    return (vec4_t) {plane.x / mag, plane.y / mag, plane.z / mag, plane.w / mag};
}

static inline camera_t makeCamera(vec3_t position, vec3_t direction, vec3_t up,
                                  float fov, float aspectRatio, float near, float far,
                                  float movementSpeed, float turningSpeed) {
    direction = normalize(direction);
    up = normalize(up);
    vec3_t right = normalize(crossProduct(up, direction));
    vec3_t correctedUp = crossProduct(direction, right);
    
    mat4x4_t rotationMatrix = (mat4x4_t) {{
        {right.x,       right.y,       right.z,       0},
        {correctedUp.x, correctedUp.y, correctedUp.z, 0},
        {direction.x,   direction.y,   direction.z,   0},
        {0,             0,             0,             1}
    }};

    // Create the translation matrix
    mat4x4_t translationMatrix = (mat4x4_t) {{
        {1, 0, 0, -position.x},
        {0, 1, 0, -position.y},
        {0, 0, 1, -position.z},
        {0, 0, 0, 1          }
    }};

    // Create the view matrix
    mat4x4_t viewMatrix = mulMM4(rotationMatrix, translationMatrix);

    float fovRadians = fov * M_PI / 180.0;
    float yScale = 1.0 / tan(fovRadians / 2.0);
    float xScale = yScale / aspectRatio;
    float zScale = far / (far - near);
    float zTranslate = -near * zScale;

    mat4x4_t projectionMatrix = (mat4x4_t) {{
        {xScale, 0,      0,          0         },
        {0,      yScale, 0,          0         },
        {0,      0,      zScale,     zTranslate},
        {0,      0,      1,          0         }
    }};

    mat4x4_t viewProjMatrix = mulMM4(projectionMatrix, viewMatrix);

    // Compute the view frustum planes
    vec4_t viewProjMatrixCol0 = {viewProjMatrix.data[0][0], viewProjMatrix.data[0][1], viewProjMatrix.data[0][2], viewProjMatrix.data[0][3]};
    vec4_t viewProjMatrixCol1 = {viewProjMatrix.data[1][0], viewProjMatrix.data[1][1], viewProjMatrix.data[1][2], viewProjMatrix.data[1][3]};
    vec4_t viewProjMatrixCol2 = {viewProjMatrix.data[2][0], viewProjMatrix.data[2][1], viewProjMatrix.data[2][2], viewProjMatrix.data[2][3]};
    vec4_t viewProjMatrixCol3 = {viewProjMatrix.data[3][0], viewProjMatrix.data[3][1], viewProjMatrix.data[3][2], viewProjMatrix.data[3][3]};
    vec4_t leftPlane = normalizePlane(addV4(viewProjMatrixCol3, viewProjMatrixCol0));
    vec4_t rightPlane = normalizePlane(subV4(viewProjMatrixCol3, viewProjMatrixCol0));
    vec4_t bottomPlane = normalizePlane(addV4(viewProjMatrixCol3, viewProjMatrixCol1));
    vec4_t topPlane = normalizePlane(subV4(viewProjMatrixCol3, viewProjMatrixCol1));
    vec4_t nearPlane = normalizePlane(viewProjMatrixCol2);
    vec4_t farPlane = normalizePlane(addV4(viewProjMatrixCol3, viewProjMatrixCol2));

    return (camera_t) {
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

static inline int triangleInFustrum(vec4_t v1, vec4_t v2, vec4_t v3) {
    // Using NDC coordinates
    if (v1.x < -1 && v2.x < -1 && v3.x < -1) return 0;
    if (v1.x >  1 && v2.x >  1 && v3.x >  1) return 0;
    if (v1.y < -1 && v2.y < -1 && v3.y < -1) return 0;
    if (v1.y >  1 && v2.y >  1 && v3.y >  1) return 0;
    if (v1.z <  0 && v2.z <  0 && v3.z <  0) return 0;
    if (v1.z >  1 && v2.z >  1 && v3.z >  1) return 0;
    return 1;
}

/* DRAWING */

typedef struct {
    uint32_t* frameBuffer;
    int       width;
    int       height;
    int       hasDepthBuffer;
    float*    depthBuffer;
} canvas_t;

static inline int edgeCross(int ax, int ay, int bx, int by, int px, int py) {
  int abx = bx - ax;
  int aby = by - ay;
  int apx = px - ax;
  int apy = py - ay;
  return abx * apy - aby * apx;
}

static inline void drawPixel(int i, int j, float z, uint32_t color, canvas_t canvas) {
    if ((i >= 0) && (i < canvas.width) && (j >= 0) && (j < canvas.height)) {
        int position = j * canvas.width + i;
        canvas.frameBuffer[position] = color;
        canvas.depthBuffer[position] = z;
    }
}

static inline void drawLine(int x0, int x1, int y0, int y1, uint32_t color, canvas_t canvas) {
    int delta_x = (x1 - x0);
    int delta_y = (y1 - y0);
    int longest_side_length = (abs(delta_x) >= abs(delta_y)) ? abs(delta_x) : abs(delta_y);
    float x_inc = delta_x / (float)longest_side_length; 
    float y_inc = delta_y / (float)longest_side_length;
    float current_x = x0;
    float current_y = y0;
    for (int i = 0; i <= longest_side_length; i++) {
        drawPixel(round(current_x), round(current_y), 0.0, color, canvas);
        current_x += x_inc;
        current_y += y_inc;
    }
}

static inline void drawTriangleWireframe(int x0, int x1, int x2,
                           int y0, int y1, int y2,
                           uint32_t color, canvas_t canvas) {
    drawLine(x0, x1, y0, y1, color, canvas);
    drawLine(x1, x2, y1, y2, color, canvas);
    drawLine(x2, x0, y2, y0, color, canvas);
}

typedef struct {
    vec3_t position;
    vec3_t normal;
    vec3_t textureCoord;
    vec3_t diffuseColor;
    vec3_t specularColor;
    float specularExponent;
} vertex_input_t;

typedef struct {
    vec4_t position;
    int numAttributes;
    float attributes[MAX_VERTEX_ATTRIBUTES];
} shader_context_t;

typedef shader_context_t vertexShader_t(void* inputVertex, void* uniformData);
// TODO: Should we pass the texture as a parameter as we're doing now? What happens when we have multiple textures?
typedef uint32_t fragmentShader_t(shader_context_t* input, void* uniformData, int textureWidth, int textureHeight, uint32_t* texture);

// TODO: Pass texture as a canvas_t
void drawObject(object3D_t* object, void *uniformData, camera_t camera, canvas_t canvas, vertexShader_t vertexShader, fragmentShader_t fragmentShader, uint16_t renderOptions) {
    mesh_t* mesh = object->mesh;

    // Don't draw if the object is fully outside the camera fustrum
    if (renderOptions & FUSTRUM_CULLING) {
        for (int p = 0; p < 6; p++) {
            vec4_t plane = camera.fustrumPlanes[p];
            int isInside = 1;
            vec4_t center = mulMV4(object->transform, (vec4_t) {mesh->center.x, mesh->center.y, mesh->center.z, 1});
            if (dotV4(plane, center) < -(object->scale * mesh->boundsRadius)) {
                DEBUG_PRINT("DEBUG: Culled object using fustrum culling {plane %d}\n", p);
                return;
            }
        }
    }

    for (int tri = 0; tri < mesh->numTriangles; tri++) {
        triangle_t triangle = mesh->triangles[tri];
        shader_context_t vertexShaderOutput[3];
        int xs[3];
        int ys[3];
        float zs[3];
        float invws[3];

        // Get vertex data
        vec3_t vertices[3] = {mesh->vertices[triangle.v0], mesh->vertices[triangle.v1], mesh->vertices[triangle.v2]};
        vec3_t normals[3] = {0};
        if (mesh->numNormals != 0) {
            normals[0] = mesh->normals[triangle.n0];
            normals[1] = mesh->normals[triangle.n1];
            normals[2] = mesh->normals[triangle.n2];
        }
        vec3_t textureCoords[3] = {0};
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
        int textureWidth, textureHeight = 0;
        uint32_t* texture = NULL;
        if (mesh->numMaterials != 0) {
            material_t material = mesh->materials[triangle.materialIndex];
            colorToFloats(material.diffuseColor, &diffR, &diffG, &diffB);
            colorToFloats(material.specularColor, &specR, &specG, &specB);
            specularExponent = material.specularExponent;
            textureWidth = material.textureWidth;
            textureHeight = material.textureHeight;
            texture = material.texture;
        }
        

        for (int v = 0; v < 3; v++) {
            vertex_input_t inputVertex = {
                .position = vertices[v],
                .normal = normals[v],
                .textureCoord = textureCoords[v],
                .diffuseColor = {diffR, diffG, diffB},
                .specularColor = {specR, specG, specB},
                .specularExponent = specularExponent
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

        int area = edgeCross(xs[0], ys[0], xs[1], ys[1], xs[2], ys[2]);

        // Backface culling
        if ((renderOptions & BACKFACE_CULLING) && area <= 0) {
            DEBUG_PRINT("DEBUG: Culled triangle using backface culling\n");
            continue;
        }

        if ((renderOptions & FUSTRUM_CULLING) && !triangleInFustrum(vertexShaderOutput[0].position,
                                                                       vertexShaderOutput[1].position,
                                                                       vertexShaderOutput[2].position)) {
            DEBUG_PRINT("DEBUG: Culled triangle using fustrum culling\n");
            continue;
        }

        // Rasterization
        int x_min = MAX(MIN(MIN(xs[0], xs[1]), xs[2]), 0);
        int x_max = MIN(MAX(MAX(xs[0], xs[1]), xs[2]), canvas.width - 1);
        int y_min = MAX(MIN(MIN(ys[0], ys[1]), ys[2]), 0);
        int y_max = MIN(MAX(MAX(ys[0], ys[1]), ys[2]), canvas.height - 1);

        int delta_w0_col = (ys[1] - ys[2]);
        int delta_w1_col = (ys[2] - ys[0]);
        int delta_w2_col = (ys[0] - ys[1]);
        int delta_w0_row = (xs[2] - xs[1]);
        int delta_w1_row = (xs[0] - xs[2]);
        int delta_w2_row = (xs[1] - xs[0]);

        int w0_row = edgeCross(xs[1], ys[1], xs[2], ys[2], x_min, y_min);
        int w1_row = edgeCross(xs[2], ys[2], xs[0], ys[0], x_min, y_min);
        int w2_row = edgeCross(xs[0], ys[0], xs[1], ys[1], x_min, y_min);

        float invArea = 1.0f / area;

        for (int y = y_min; y <= y_max; y++) {
            int was_inside = 0;
            int w0 = w0_row;
            int w1 = w1_row;
            int w2 = w2_row;
            for (int x = x_min; x <= x_max; x++) {
                // Perspective correct baricentric coordinates
                // TODO: Maybe I can avoid dividing by sum here?
                float alpha = w0 * invArea * invws[0];
                float beta  = w1 * invArea * invws[1];
                float gamma = w2 * invArea * invws[2];
                float sum = alpha + beta + gamma;
                alpha /= sum;
                beta  /= sum;
                gamma /= sum;

                // Check if the fragment is inside the triangle
                int is_inside = alpha >= 0 && beta >= 0 && gamma >= 0;
                if (is_inside) {
                    was_inside = 1;

                    // Interpolate z
                    float z = alpha * zs[0] + beta * zs[1] + gamma * zs[2];

                    // Depth test
                    if (z < canvas.depthBuffer[y * canvas.width + x]) {
                        // Compute fragment input attributes from the outputs of the vertex shader
                        shader_context_t fragmentShaderInput = {0};
                        if (renderOptions & FLAT_SHADING) {
                            fragmentShaderInput.position = vertexShaderOutput[0].position;
                            fragmentShaderInput.numAttributes = vertexShaderOutput[0].numAttributes;
                            for (int i = 0; i < fragmentShaderInput.numAttributes; i++) {
                                fragmentShaderInput.attributes[i] = vertexShaderOutput[0].attributes[i];
                            }
                        } else {
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
                        }

                        uint32_t color = fragmentShader(&fragmentShaderInput, uniformData, textureWidth, textureHeight, texture);
                        drawPixel(x, y, z, color, canvas); // TODO: Avoid scissor test in drawPixel
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
}

/* SHADER IMPLEMENTATIONS */

/* Basic shading */
// Draw with a single color, no lighting or textures
typedef struct {
    mat4x4_t modelviewprojection;
} basic_uniform_t;

static inline shader_context_t basicVertexShader(void* inputVertex, void* uniformData) {
    vertex_input_t* inputVertexData = (vertex_input_t*) inputVertex;
    shader_context_t result = {0};
    basic_uniform_t* basicUniformData = (basic_uniform_t*) uniformData;
    vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};
    result.position = mulMV4(basicUniformData->modelviewprojection, inputVertex4);
    return result;
}

static inline uint32_t basicFragmentShader(const shader_context_t* input, void* uniformData, int textureWidth, int textureHeight, uint32_t* texture) {
    return COLOR_WHITE;
}

/* Gouraud shading */
// Compute the lighting at each vertex
typedef struct {
    mat4x4_t modelMatrix;
    mat4x4_t modelInvRotationMatrixTransposed;
    mat4x4_t viewProjectionMatrix;
    light_sources_t lightSources;
    int bilinearFiltering;
} gourard_uniform_t;

static inline shader_context_t gourardVertexShader(void* inputVertex, void* uniformData) {
    vertex_input_t* inputVertexData = (vertex_input_t*) inputVertex;
    shader_context_t result = {0};
    gourard_uniform_t* defaultUniformData = (gourard_uniform_t*) uniformData;
    vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};    
    vec4_t worldSpaceVertex = mulMV4(defaultUniformData->modelMatrix, inputVertex4); // Local to world space
    result.position = mulMV4(defaultUniformData->viewProjectionMatrix, worldSpaceVertex); // World to clip space
    
    // Set other vertex attributes
    result.numAttributes = 13;

    vec3_t worldSpaceNormal = mulMV3(defaultUniformData->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    float invMagNormal = 1.0f / magnitude(worldSpaceNormal);
    result.attributes[0] = worldSpaceNormal.x;
    result.attributes[1] = worldSpaceNormal.y;
    result.attributes[2] = worldSpaceNormal.z;
    result.attributes[3] = inputVertexData->textureCoord.x;   // u
    result.attributes[4] = inputVertexData->textureCoord.y;   // v
    result.attributes[5] = inputVertexData->diffuseColor.x;   // R
    result.attributes[6] = inputVertexData->diffuseColor.y;   // G
    result.attributes[7] = inputVertexData->diffuseColor.z;   // B
    result.attributes[8] = inputVertexData->specularColor.x;  // R
    result.attributes[9] = inputVertexData->specularColor.y;  // G
    result.attributes[10] = inputVertexData->specularColor.z; // B
    result.attributes[11] = inputVertexData->specularExponent;
    result.attributes[12] = computeLighting((vec3_t) {worldSpaceVertex.x, worldSpaceVertex.y, worldSpaceVertex.z}, worldSpaceNormal, invMagNormal, inputVertexData->specularExponent, defaultUniformData->lightSources, DIFFUSE_LIGHTING | SPECULAR_LIGHTING);
    return result;
}

static inline uint32_t gourardFragmentShader(const shader_context_t* input, void* uniformData, int textureWidth, int textureHeight, uint32_t* texture) {
    gourard_uniform_t* uniform = (gourard_uniform_t*) uniformData;

    uint32_t unshadedColor;
    if (textureHeight == 0 || textureWidth == 0) {
        unshadedColor = colorFromFloats(input->attributes[5], input->attributes[6], input->attributes[7]);
    } else {
        float u = input->attributes[3];
        float v = input->attributes[4];
        float tex_u = MIN(fabs(u * textureWidth), textureWidth - 1);
        float tex_v = MIN(fabs(v * textureHeight), textureHeight - 1);
        int floor_u = floor(tex_u);
        int floor_v = floor(tex_v);
        
        if (uniform->bilinearFiltering) {
            // Bilinear filtering
            float ratio_u = tex_u - floor_u;
            float ratio_v = tex_v - floor_v;
            int next_u = MIN(floor_u + 1, textureWidth - 1);
            int next_v = MIN(floor_v + 1, textureHeight - 1);
            uint32_t color00 = texture[floor_v * textureWidth + floor_u];
            uint32_t color10 = texture[floor_v * textureWidth + next_u];
            uint32_t color01 = texture[next_v * textureWidth + floor_u];
            uint32_t color11 = texture[next_v * textureWidth + next_u];
            uint32_t color0 = sumColors(mulScalarColor(1.0f - ratio_u, color00), mulScalarColor(ratio_u, color10));
            uint32_t color1 = sumColors(mulScalarColor(1.0f - ratio_u, color01), mulScalarColor(ratio_u, color11));
            unshadedColor = sumColors(mulScalarColor(1.0f - ratio_v, color0), mulScalarColor(ratio_v, color1));
        } else {
            unshadedColor = texture[floor_v * textureWidth + floor_u];
        }
    }
    return mulScalarColor(input->attributes[12], unshadedColor);
}


/* Phong shading */
// Compute the lighting at each fragment
typedef struct {
    mat4x4_t modelMatrix;
    mat4x4_t modelInvRotationMatrixTransposed;
    mat4x4_t viewProjectionMatrix;
    light_sources_t lightSources;
    int bilinearFiltering;
} phong_uniform_t;

static inline shader_context_t phongVertexShader(void* inputVertex, void* uniformData) {
    vertex_input_t* inputVertexData = (vertex_input_t*) inputVertex;
    shader_context_t result = {0};
    phong_uniform_t* uniform = (phong_uniform_t*) uniformData;
    vec4_t inputVertex4 = {inputVertexData->position.x, inputVertexData->position.y, inputVertexData->position.z, 1.0f};    
    vec4_t worldSpaceVertex = mulMV4(uniform->modelMatrix, inputVertex4); // Local to world space
    result.position = mulMV4(uniform->viewProjectionMatrix, worldSpaceVertex); // World to clip space
    
    // Set other vertex attributes
    result.numAttributes = 15;
    vec3_t worldSpaceNormal = mulMV3(uniform->modelInvRotationMatrixTransposed, inputVertexData->normal); // Local to world space
    result.attributes[0] = worldSpaceNormal.x;
    result.attributes[1] = worldSpaceNormal.y;
    result.attributes[2] = worldSpaceNormal.z;
    result.attributes[3] = worldSpaceVertex.x;
    result.attributes[4] = worldSpaceVertex.y;
    result.attributes[5] = worldSpaceVertex.z;
    result.attributes[6] = inputVertexData->textureCoord.x;   // u
    result.attributes[7] = inputVertexData->textureCoord.y;   // v
    result.attributes[8] = inputVertexData->diffuseColor.x;   // R
    result.attributes[9] = inputVertexData->diffuseColor.y;   // G
    result.attributes[10] = inputVertexData->diffuseColor.z;   // B
    result.attributes[11] = inputVertexData->specularColor.x;  // R
    result.attributes[12] = inputVertexData->specularColor.y;  // G
    result.attributes[13] = inputVertexData->specularColor.z; // B
    result.attributes[14] = inputVertexData->specularExponent;
    return result;
}

static inline uint32_t phongFragmentShader(const shader_context_t* input, void* uniformData, int textureWidth, int textureHeight, uint32_t* texture) {
    phong_uniform_t* uniform = (phong_uniform_t*) uniformData;

    vec3_t normal = {input->attributes[0], input->attributes[1], input->attributes[2]};
    vec3_t position = {input->attributes[3], input->attributes[4], input->attributes[5]};
    float specularExponent = input->attributes[14];
    float invMagNormal = 1.0f / magnitude(normal);
    float lighting = computeLighting(position, normal, invMagNormal, specularExponent, uniform->lightSources, DIFFUSE_LIGHTING | SPECULAR_LIGHTING);
    
    uint32_t unshadedColor;
    if (textureHeight == 0 || textureWidth == 0) {
        unshadedColor = colorFromFloats(input->attributes[8], input->attributes[9], input->attributes[10]);
    } else {
        float u = input->attributes[6];
        float v = input->attributes[7];
        float tex_u = MIN(fabs(u * textureWidth), textureWidth - 1);
        float tex_v = MIN(fabs(v * textureHeight), textureHeight - 1);
        int floor_u = floor(tex_u);
        int floor_v = floor(tex_v);
        
        if (uniform->bilinearFiltering) {
            float ratio_u = tex_u - floor_u;
            float ratio_v = tex_v - floor_v;
            int next_u = MIN(floor_u + 1, textureWidth - 1);
            int next_v = MIN(floor_v + 1, textureHeight - 1);
            uint32_t color00 = texture[floor_v * textureWidth + floor_u];
            uint32_t color10 = texture[floor_v * textureWidth + next_u];
            uint32_t color01 = texture[next_v * textureWidth + floor_u];
            uint32_t color11 = texture[next_v * textureWidth + next_u];
            uint32_t color0 = sumColors(mulScalarColor(1.0f - ratio_u, color00), mulScalarColor(ratio_u, color10));
            uint32_t color1 = sumColors(mulScalarColor(1.0f - ratio_u, color01), mulScalarColor(ratio_u, color11));
            unshadedColor = sumColors(mulScalarColor(1.0f - ratio_v, color0), mulScalarColor(ratio_v, color1));
        } else {
            unshadedColor = texture[floor_v * textureWidth + floor_u];
        }
    }
    return mulScalarColor(lighting, unshadedColor);
}

#endif // ZEROGL_H
