#ifndef SIMPLERENDERER_H
#define SIMPLERENDERER_H

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#ifdef SR_DEBUG
#define DEBUG_PRINT(...) printf(__VA_ARGS__)
#else
#define DEBUG_PRINT(...) do {} while (0)
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/* VECTORS AND MATRICES */

typedef struct vec3_t {
  float x, y, z;
} vec3_t;

typedef struct vec4_t {
  float x, y, z, w;
} vec4_t;

typedef struct mat4x4_t {
  float data[4][4];
} mat4x4_t;

typedef struct plane_t {
    vec3_t   normal;
    float    distance;
} plane_t;

static const struct mat4x4_t IDENTITY_M4x4 = {{
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
}};

static inline void printMatrix(mat4x4_t matrix) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            DEBUG_PRINT("%f ", matrix.data[i][j]);
        }
        DEBUG_PRINT("\n");
    }
}

static inline void printVertex(vec3_t vertex) {
    DEBUG_PRINT("%f %f %f\n", vertex.x, vertex.y, vertex.z);
}

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

static inline float magnitude(vec3_t v) {
    float result = sqrt(dot(v, v));
    assert(result >= 0);
    return result;
}

static inline vec3_t sub(vec3_t a, vec3_t b) {
    return (vec3_t) {a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline vec3_t add(vec3_t a, vec3_t b) {
    return (vec3_t) {a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline vec3_t normalize(vec3_t v) {
    float mag = magnitude(v);
    return (vec3_t) {v.x / mag, v.y / mag, v.z / mag};
}

static inline vec3_t mulScalarV3(float k, vec3_t v) {
    return (vec3_t) {k*v.x, k*v.y, k*v.z};
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
    int i;

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

    for (i = 0; i < 16; i++)
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

static inline float distancePlaneV3(plane_t plane, vec3_t v) {
    return dot(plane.normal, v) + plane.distance;
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

/* 3D OBJECTS */

#define M_PI 3.14159265358979323846264338327950288

typedef struct triangle_t {
  int     v0, v1, v2;
  int     t0, t1, t2;
  int     n0, n1, n2;
  int     materialIndex;
} triangle_t;

typedef struct material_t {
    char*   name;
    uint32_t diffuseColor;
    uint32_t specularColor;
    float   specularExponent;
    int     textureWidth;
    int     textureHeight;
    uint32_t* texture;
} material_t;

typedef struct mesh_t {
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

typedef struct object3D_t {
    mesh_t*  mesh;
    vec3_t   translation;
    float    scale;
    mat4x4_t rotation;
    mat4x4_t transform;
} object3D_t;

typedef struct camera_t {
    vec3_t   translation;
    mat4x4_t rotation;
    mat4x4_t transform;
    int      numPlanes;
    plane_t* planes;
    float    viewportWidth;
    float    viewportHeight;
    float    viewportDistance;
    float    movementSpeed;
    float    turningSpeed;
} camera_t;

typedef struct ambient_light_t {
    float intensity;
} ambient_light_t;

typedef struct dir_light_t {
    float intensity;
    vec3_t direction;
} dir_light_t;

typedef struct point_light_t {
    float intensity;
    vec3_t position;
} point_light_t;

typedef struct light_sources_t {
    int              numAmbientLights;
    ambient_light_t* ambientLights;
    int              numDirectionalLights;
    dir_light_t*     directionalLights;
    int              numPointLights;
    point_light_t*   pointLights;
} light_sources_t;

static inline vec3_t triangleNormal(vec3_t v0, vec3_t v1, vec3_t v2) {
    return crossProduct(sub(v1, v0), sub(v2, v0));
}

static inline vec3_t triangleCenter(vec3_t v0, vec3_t v1, vec3_t v2) {
    struct vec3_t result;

    result.x = (v0.x + v1.x + v2.x) / 3.0f;
    result.y = (v0.y + v1.y + v2.y) / 3.0f;
    result.z = (v0.z + v1.z + v2.z) / 3.0f;

    return result;
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

static inline object3D_t makeObject(mesh_t *mesh, vec3_t translation, float scale, mat4x4_t rotation) {
    mat4x4_t translationMatrix = translationToMatrix(translation);
    mat4x4_t scaleMatrix = scaleToMatrix(scale);
    mat4x4_t transform = mulMM4(translationMatrix, mulMM4(rotation, scaleMatrix));
    return (object3D_t) {mesh, translation, scale, rotation, transform};
}

static inline camera_t makeCamera(vec3_t translation, mat4x4_t rotation,
                    float viewportDist, float viewportWidth, float viewportHeight,
                    float movSpeed, float turnSpeed) {
    mat4x4_t rotationMatrix = transposeM4(rotation);
    mat4x4_t translationMatrix = translationToMatrix(mulScalarV3(-1.0, translation));
    mat4x4_t transform = mulMM4(rotationMatrix, translationMatrix);
    int numPlanes = 5;
    const float sqrt2 = 1/sqrt(2);

    plane_t* planes = (plane_t*) malloc(numPlanes * sizeof(plane_t)); // allocate memory dynamically
    if (!planes) {
        fprintf(stderr, "ERROR: Failed to allocate memory for planes.\n");
        exit(1);
    }

    // FIXME: This is wrong. The FOV depends on the viewport size and distance
    //        and with the current value it's 53 not 90 degrees. But the plane
    //        normals (l, r, t and b) are set to 90 degree FOV. We should compute
    //        the right values based on the camera parameters.
    planes[0] = (plane_t) {{0,      0,      1    }, viewportDist}; // Near
    planes[1] = (plane_t) {{sqrt2,  0,      sqrt2}, 0           }; // Left
    planes[2] = (plane_t) {{-sqrt2, 0,      sqrt2}, 0           }; // Right
    planes[3] = (plane_t) {{0,      sqrt2,  sqrt2}, 0           }; // Top
    planes[4] = (plane_t) {{0,      -sqrt2, sqrt2}, 0           }; // Bottom

    return (camera_t) {
        translation, rotation, transform, numPlanes,  planes,
        viewportDist, viewportWidth, viewportHeight,
        movSpeed, turnSpeed
    };
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

/* DRAWING */

// Rendering options
#define DIFFUSE_LIGHTING (1 << 0)
#define SPECULAR_LIGHTING (1 << 1)
#define BACKFACE_CULLING (1 << 2)
#define BILINEAR_FILTERING (1 << 3)
#define SHADED (1 << 4)
#define SHADED_FLAT (1 << 5)
#define SHADED_GOURAUD (1 << 6)
#define SHADED_PHONG (1 << 7)
#define DRAW_WIREFRAME (1 << 8)
#define DRAW_FILLED (1 << 9)

typedef struct canvas_t {
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

static inline vec3_t unprojectPoint(int x, int y, float invz, canvas_t canvas, camera_t cam) {
  return (vec3_t) {
    (x - canvas.width/2) * (cam.viewportWidth / canvas.width) / (invz * cam.viewportDistance),
    (canvas.height/2 - y - 1) / (cam.viewportDistance * invz * canvas.height/cam.viewportHeight),
    1.0f / invz
  };
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

// TODO: Code here is a bit repeated between directional and point lights. Maybe refactor?
float shadeVertex(vec3_t v, vec3_t normal, float invMagnitudeNormal, float specularExponent,
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
        vec3_t lightDirection = sub(pointLights[i].position, v);
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

// TODO: It's very weird to have the camera as a parameter here. This is because the phong shading
//       needs the camera position. Maybe I should get rid of phong shading until I implement
//       programmable shaders.
// TODO: Pass texture as a canvas_t
void drawTriangleFilled(int x0, int x1, int x2,
                        int y0, int y1, int y2,
                        float invz0, float invz1, float invz2,
                        float i0, float i1, float i2,
                        vec3_t n0, vec3_t n1, vec3_t n2,
                        vec3_t t0, vec3_t t1, vec3_t t2,
                        uint32_t c0, uint32_t c1, uint32_t c2,
                        float specularExponent,
                        uint32_t* texture, int textureWidth, int textureHeight,
                        mat4x4_t invCameraTransform,
                        int area,
                        light_sources_t lightSources, camera_t camera,
                        canvas_t canvas, uint16_t renderOptions) {
    int x_min = MAX(MIN(MIN(x0, x1), x2), 0);
    int x_max = MIN(MAX(MAX(x0, x1), x2), canvas.width - 1); 
    int y_min = MAX(MIN(MIN(y0, y1), y2), 0);
    int y_max = MIN(MAX(MAX(y0, y1), y2), canvas.height - 1);

    float invArea = 1.0f / area;

    // Compute the constant delta_s that will be used for the horizontal and vertical steps
    int delta_w0_col = (y1 - y2);
    int delta_w1_col = (y2 - y0);
    int delta_w2_col = (y0 - y1);
    int delta_w0_row = (x2 - x1);
    int delta_w1_row = (x0 - x2);
    int delta_w2_row = (x1 - x0);

    // Compute the edge functions for the fist (top-left) point
    int w0_row = edgeCross(x1, y1, x2, y2, x_min, y_min);
    int w1_row = edgeCross(x2, y2, x0, y0, x_min, y_min);
    int w2_row = edgeCross(x0, y0, x1, y1, x_min, y_min);

    // Illuminate each vertex
    if ((renderOptions & SHADED) && !(renderOptions & SHADED_PHONG)) {
        c0 = mulScalarColor(i0, c0);
        c1 = mulScalarColor(i1, c1);
        c2 = mulScalarColor(i2, c2);
    }
    
    for (int y = y_min; y <= y_max; y++) {
        bool was_inside = false;
        int w0 = w0_row;
        int w1 = w1_row;
        int w2 = w2_row;
        for (int x = x_min; x <= x_max; x++) {
            bool is_inside = (w0 | w1 | w2) >= 0;
            if (is_inside) {
                was_inside = true;
            
                float alpha = w0 * invArea;
                float beta  = w1 * invArea;
                float gamma = w2 * invArea;
                float invz = alpha * invz0 + beta * invz1 + gamma * invz2;
                if (invz > canvas.depthBuffer[y * canvas.width + x]) {
                    uint32_t color = c0; // Fallback in case of no texture and no shading
                    float light = 1;

                    if (renderOptions & SHADED_PHONG) {
                        vec3_t v = mulMV3(invCameraTransform, unprojectPoint(x, y, invz, canvas, camera));
                        vec3_t normal = add(add(mulScalarV3(alpha, n0), mulScalarV3(beta, n1)), mulScalarV3(gamma, n2));
                        light = shadeVertex(v, normal , 1/magnitude(normal), specularExponent, lightSources, renderOptions);
                    } else if (renderOptions & SHADED) {
                        light = alpha * i0 + beta * i1 + gamma * i2;
                    }
                    
                    if (textureWidth != 0 && textureHeight != 0) {
                        // Interpolate u/z and v/z to get perspective correct texture coordinates
                        float u_over_z = alpha * (t0.x * invz0) + beta * (t1.x * invz1) + gamma * (t2.x * invz2);
                        float v_over_z = alpha * (t0.y * invz0) + beta * (t1.y * invz1) + gamma * (t2.y * invz2);
                        uint32_t unshaded_color;
                        // TODO: Fix crash when we have overflow here
                        if (renderOptions & BILINEAR_FILTERING) {
                            float tex_u = u_over_z/invz;
                            if (tex_u < 0) {
                                tex_u = 1 + tex_u;
                            }
                            tex_u = MIN(tex_u * textureWidth, textureWidth - 1);

                            float tex_v = v_over_z/invz;
                            if (tex_v < 0) {
                                tex_v = 1 + tex_v;
                            }
                            tex_v = MIN(tex_v * textureHeight, textureHeight - 1);

                            int floor_u = floor(tex_u);
                            int floor_v = floor(tex_v);
                            int next_u = MIN(floor_u + 1, textureWidth - 1);
                            int next_v = MIN(floor_v + 1, textureHeight - 1);
                            float frac_u = tex_u - floor_u;
                            float frac_v = tex_v - floor_v;
                            uint32_t color_tl = texture[floor_v * textureWidth + floor_u];
                            uint32_t color_tr = texture[floor_v * textureWidth + next_u];
                            uint32_t color_bl = texture[next_v * textureWidth + floor_u];
                            uint32_t color_br = texture[next_v * textureWidth + next_u];
                            uint32_t color_b = sumColors(mulScalarColor(1 - frac_u, color_bl), mulScalarColor(frac_u, color_br));
                            uint32_t color_tp = sumColors(mulScalarColor(1 - frac_u, color_tl), mulScalarColor(frac_u, color_tr));
                            unshaded_color = sumColors(mulScalarColor(1 - frac_v, color_b), mulScalarColor(frac_v, color_tp));
                        } else {
                            int tex_x = MIN(abs((int)((u_over_z/invz) * textureWidth)), textureWidth - 1);
                            int tex_y = MIN(abs((int)((v_over_z/invz) * textureHeight)), textureHeight - 1);
                            unshaded_color = texture[tex_y * textureWidth + tex_x];
                        }
                        
                        color = mulScalarColor(light, unshaded_color);
                    } else if (renderOptions & SHADED) {
                        // TODO: There should be an easier route to interpolate colors
                        uint8_t c0r, c0g, c0b;
                        uint8_t c1r, c1g, c1b;
                        uint8_t c2r, c2g, c2b;
                        colorFromUint32(c0, &c0r, &c0g, &c0b);
                        colorFromUint32(c1, &c1r, &c1g, &c1b);
                        colorFromUint32(c2, &c2r, &c2g, &c2b);
                        uint32_t unshaded_color = colorToUint32((c0r * alpha + c1r * beta + c2r * gamma),
                                                                (c0g * alpha + c1g * beta + c2g * gamma),
                                                                (c0b * alpha + c1b * beta + c2b * gamma));
                        color = mulScalarColor(light, unshaded_color);
                    }
                    drawPixel(x, y, invz, color, canvas);
                }
            }

            // Go to next row if we jumped outside the triangle
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

// TODO: This function is very dumb, as it is allocating transformed
//       vertices just to allocate the projections right after.
void drawObject(object3D_t* object, light_sources_t lightSources, camera_t camera, canvas_t canvas, uint16_t renderOptions) {
    mesh_t* mesh = object->mesh;
    material_t* materials = mesh->materials;
    mat4x4_t invCameraTransform = inverseM4(camera.transform);
    
    // Transform the bounding sphere and don't draw when object is fully out of the camera volume
    mat4x4_t transform = mulMM4(camera.transform, object->transform);
    vec3_t transformedCenter = mulMV3(transform, mesh->center);
    float transformedBoundsRadius = mesh->boundsRadius * object->scale;
    for (int p = 0; p < camera.numPlanes; p++) {
        float distance = distancePlaneV3(camera.planes[p], transformedCenter);
        if (distance < -transformedBoundsRadius) {
            DEBUG_PRINT("DEBUG: Clipped an object fully outside of the clipping volume\n");
            return;
        }
    }

    // If the object is not discarded, transform and project all vertices
    vec3_t *transformed =    (vec3_t*) malloc(mesh->numVertices * sizeof(vec3_t));
    vec3_t *camTransformed = (vec3_t*) malloc(mesh->numVertices * sizeof(vec3_t));
    int   *projected_x =     (int*) malloc(mesh->numVertices * sizeof(int));
    int   *projected_y =     (int*) malloc(mesh->numVertices * sizeof(int));
    float *projected_invz =  (float*) malloc(mesh->numVertices * sizeof(float));

    vec3_t *transformedNormals = (vec3_t*) malloc(mesh->numNormals * sizeof(vec3_t));
    if (projected_x == NULL || projected_y == NULL || projected_invz == NULL || transformed == NULL || camTransformed == NULL || transformedNormals == NULL) {
        fprintf(stderr, "ERROR: Transformed vertices/normals memory couldn't be allocated.\n");
        exit(-1);
    }

    for (int i = 0; i < mesh->numVertices; i++) {
        transformed[i] = mulMV3(object->transform, mesh->vertices[i]);
        vec3_t v = mulMV3(camera.transform, transformed[i]);
        camTransformed[i] = v;
        // Project the vertex into screen coordinates
        projected_x[i] = (int) (v.x * camera.viewportDistance / v.z  * canvas.width/camera.viewportWidth + canvas.width/2),
        projected_y[i] = (int) (canvas.height/2 - (v.y * camera.viewportDistance / v.z * canvas.height/camera.viewportHeight) - 1),
        projected_invz[i] = 1.0f / v.z;
    }

    for (int i = 0; i < mesh->numNormals; i++) {
        transformedNormals[i] = mulMV3(object->transform, mesh->normals[i]);
    }

    // Cull, shade and draw each triangle
    for (int i = 0; i < mesh->numTriangles; i++) {
        triangle_t triangle = mesh->triangles[i];

        bool discarded = false;

        // Backface culling
        int p0x      = projected_x[triangle.v0];
        int p0y      = projected_y[triangle.v0];
        float p0invz = projected_invz[triangle.v0];
        int p1x      = projected_x[triangle.v1];
        int p1y      = projected_y[triangle.v1];
        float p1invz = projected_invz[triangle.v1];
        int p2x      = projected_x[triangle.v2];
        int p2y      = projected_y[triangle.v2];
        float p2invz = projected_invz[triangle.v2];
        int area     = edgeCross(p0x, p0y, p1x, p1y, p2x, p2y);
        if (area < 0 && (renderOptions & BACKFACE_CULLING)) {
            discarded = true;
        }

        // Fustrum culling
        // TODO: Add config for this
        for (int p = 0; !discarded && p < camera.numPlanes; p++) {
            plane_t plane = camera.planes[p];
            if (distancePlaneV3(plane, camTransformed[triangle.v0]) < 0 &&
                distancePlaneV3(plane, camTransformed[triangle.v1]) < 0 &&
                distancePlaneV3(plane, camTransformed[triangle.v2]) < 0) {
                    DEBUG_PRINT("DEBUG: Clipped triangle fully outside of the camera clipping volume\n");
                    discarded = true;
                    break;
            }

            // TODO: Deal with the case where only one or two vertexes of the triangle
            //       are out of the volume. In this case, we should split the triangle
            //       and create new ones.
        }

        // TODO: Check if this is right, but it fixes overflows when we
        //       interpolate zs and we have z close to 0
        // Don't draw if the triangle has any vertice behind the camera
        if (camTransformed[triangle.v0].z < 0 ||
            camTransformed[triangle.v1].z < 0 ||
            camTransformed[triangle.v2].z < 0) {
            DEBUG_PRINT("DEBUG: Clipped triangle with a vertice behind the camera\n");
            discarded = true;
        }

        if (!discarded) {
            float i0 = 1.0;
            float i1 = 1.0;
            float i2 = 1.0;

            // Shading
            if ((renderOptions & SHADED) && (renderOptions & SHADED_FLAT)) {
                vec3_t v0 = transformed[triangle.v0];
                vec3_t v1 = transformed[triangle.v1];
                vec3_t v2 = transformed[triangle.v2];
                vec3_t normal = triangleNormal(v0, v1, v2);
                float invMag = 1.0f / magnitude(normal);
                vec3_t center = {(v0.x + v1.x + v2.x)/3.0f, (v0.y + v1.y + v2.y)/3.0f, (v0.z + v1.z + v2.z)/3.0f};
                float intensity = shadeVertex(center, normal, invMag, materials[triangle.materialIndex].specularExponent, lightSources, renderOptions);
                i0 = intensity;
                i1 = intensity;
                i2 = intensity;
            } else if ((renderOptions & SHADED)  && (renderOptions & SHADED_GOURAUD)) {
                float specularExponent = 900.0f;
                if (mesh->numMaterials != 0) {
                    specularExponent = materials[triangle.materialIndex].specularExponent;
                }

                i0 = shadeVertex(transformed[triangle.v0], transformedNormals[triangle.n0], mesh->invMagnitudeNormals[triangle.n0], specularExponent, lightSources, renderOptions);
                i1 = shadeVertex(transformed[triangle.v1], transformedNormals[triangle.n1], mesh->invMagnitudeNormals[triangle.n1], specularExponent, lightSources, renderOptions);
                i2 = shadeVertex(transformed[triangle.v2], transformedNormals[triangle.n2], mesh->invMagnitudeNormals[triangle.n2], specularExponent, lightSources, renderOptions);
            }

            // Drawing
            if (renderOptions & DRAW_WIREFRAME) {
                drawTriangleWireframe(p0x, p1x, p2x,
                                      p0y, p1y, p2y,
                                      materials[triangle.materialIndex].diffuseColor,
                                      canvas);
            }
            
            if (renderOptions & DRAW_FILLED) {
                vec3_t t0 = {0};
                vec3_t t1 = {0};
                vec3_t t2 = {0};
                int textureWidth = 0;
                int textureHeight = 0;
                uint32_t* texture = NULL;
                if (mesh->numTextureCoords > 0) {
                    t0 = mesh->textureCoords[triangle.t0];
                    t1 = mesh->textureCoords[triangle.t1];
                    t2 = mesh->textureCoords[triangle.t2];
                    textureWidth = materials[triangle.materialIndex].textureWidth;
                    textureHeight = materials[triangle.materialIndex].textureHeight;
                    texture = materials[triangle.materialIndex].texture;
                }
                
                uint32_t color = COLOR_WHITE;
                float specularExponent = 900.0f;

                if (mesh->numMaterials != 0) {
                    color = materials[triangle.materialIndex].diffuseColor;
                    specularExponent = materials[triangle.materialIndex].specularExponent;
                }
                
                drawTriangleFilled(p0x, p1x, p2x,
                                   p0y, p1y, p2y,
                                   p0invz, p1invz, p2invz,
                                   i0, i1, i2,
                                   transformedNormals[triangle.n0], transformedNormals[triangle.n1], transformedNormals[triangle.n2],
                                   t0, t1, t2,
                                   color, color, color,
                                   specularExponent,
                                   texture, textureWidth, textureHeight,
                                   invCameraTransform,
                                   area,
                                   lightSources, camera, canvas, renderOptions);
            }
        }
    }

    free(transformed);
    free(camTransformed);
    free(projected_x);
    free(projected_y);
    free(projected_invz);
    free(transformedNormals);
}

#endif // SIMPLERENDERER_H
