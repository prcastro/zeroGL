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

typedef struct color_t {
  uint8_t r;
  uint8_t g;
  uint8_t b;
} color_t;

static const color_t COLOR_WHITE      = {255, 255, 255};
static const color_t COLOR_BLACK      = {0,   0,   0  };
static const color_t COLOR_GREEN      = {0,   255, 0  };
static const color_t COLOR_BLUE       = {0,   0,   255};
static const color_t COLOR_RED        = {255, 0,   0  };
static const color_t COLOR_YELLOW     = {255, 255, 0  };
static const color_t COLOR_PURPLE     = {255, 0,   255};
static const color_t COLOR_CYAN       = {0,   255, 255};

static inline uint32_t colorToUint32(color_t c) {
    return 0x00000000 | (c.r << 16) | (c.g << 8) |  c.b;    
}

static inline color_t colorFromUint32(uint32_t c) {
    return (color_t) {
        (uint8_t) ((c & 0x00FF0000) >> 16),
        (uint8_t) ((c & 0x0000FF00) >> 8),
        (uint8_t) (c & 0x000000FF)
    };
}

static inline float clamp(float v, float max) {
    return v > max ? max : v;
}

static inline color_t mulScalarColor(double x, color_t color) {
    color_t result = {
        (uint8_t) clamp(x * color.r, 255.0),
        (uint8_t) clamp(x * color.g, 255.0),
        (uint8_t) clamp(x * color.b, 255.0)
    };
    return result;
}

static inline color_t sumColors3(color_t c0, color_t c1, color_t c2) {
    return (color_t) {
        (uint8_t) clamp(c0.r + c1.r + c2.r, 255.0),
        (uint8_t) clamp(c0.g + c1.g + c2.g, 255.0),
        (uint8_t) clamp(c0.b + c1.b + c2.b, 255.0)
    };
}

static inline color_t sumColors(color_t c0, color_t c1) {
    return (color_t) {
        (uint8_t) clamp(c0.r + c1.r, 255.0),
        (uint8_t) clamp(c0.g + c1.g, 255.0),
        (uint8_t) clamp(c0.b + c1.b, 255.0)
    };
}

static inline color_t colorFromFloats(float r, float g, float b) {
    return (color_t) {
        (uint8_t) clamp(r * 255.0, 255.0),
        (uint8_t) clamp(g * 255.0, 255.0),
        (uint8_t) clamp(b * 255.0, 255.0)
    };
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
    color_t diffuseColor;
    color_t specularColor;
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
    vec3_t translation;
    mat4x4_t rotation;
    mat4x4_t transform;
    int      numPlanes;
    plane_t* planes;
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
                    float viewportDist, float movSpeed, float turnSpeed) {
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

    return (camera_t) {translation, rotation, transform, numPlanes,  planes,
            viewportDist, movSpeed, turnSpeed};
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

typedef struct canvas_t {
    uint32_t* frameBuffer;
    int       width;
    int       height;
    int       hasDepthBuffer;
    float*    depthBuffer;
} canvas_t;

typedef struct point_t {
  int   x, y;
  float invz;
} point_t;

static inline void drawPixel(int i, int j, uint32_t color, canvas_t canvas) {
    if ((i >= 0) && (i < canvas.width) && (j >= 0) && (j < canvas.height)) {
        canvas.frameBuffer[j * canvas.width + i] = color;
    }
}

static inline void drawPixelDepthBuffer(int i, int j, float z, uint32_t color, canvas_t canvas) {
    if ((i >= 0) && (i < canvas.width) && (j >= 0) && (j < canvas.height)) {
        int position = j * canvas.width + i;
        canvas.frameBuffer[position] = color;
        canvas.depthBuffer[position] = z;
    }
}

static inline void drawLine(int x0, int x1, int y0, int y1, color_t color, canvas_t canvas) {
    int delta_x = (x1 - x0);
    int delta_y = (y1 - y0);

    int longest_side_length = (abs(delta_x) >= abs(delta_y)) ? abs(delta_x) : abs(delta_y);

    float x_inc = delta_x / (float)longest_side_length; 
    float y_inc = delta_y / (float)longest_side_length;

    float current_x = x0;
    float current_y = y0;
    for (int i = 0; i <= longest_side_length; i++) {
        drawPixel(round(current_x), round(current_y), colorToUint32(color), canvas);
        current_x += x_inc;
        current_y += y_inc;
    }
}

// TODO: Receive viewport info as parameters
static inline point_t projectVertex(vec3_t v, canvas_t canvas) {
  float viewportWidth = canvas.width/(float) canvas.height;
  float viewportHeight = 1.0f;
  float viewportDistance = 1.0f;
  return (point_t) {
    (int) (v.x * viewportDistance / v.z  * canvas.width/viewportWidth + canvas.width/2),
    (int) (canvas.height/2 - (v.y * viewportDistance / v.z * canvas.height/viewportHeight) - 1),
    1.0f / v.z
  };
}

// TODO: Receive viewport info as parameters
static inline vec3_t unprojectPoint(point_t p, canvas_t canvas) {
  float viewportWidth = canvas.width/(float) canvas.height;
  float viewportHeight = 1.0f;
  float viewportDistance = 1.0f;
  return (vec3_t) {
    (p.x - canvas.width/2) * (viewportWidth / canvas.width) / (p.invz * viewportDistance),
    (canvas.height/2 - p.y - 1) / (viewportDistance * p.invz * canvas.height/viewportHeight),
    1.0f / p.invz
  };
}

#endif // SIMPLERENDERER_H
