#ifndef VECTORS_H
#define VECTORS_H

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

// Print matrix for debugging
void printVertex(vec3_t vertex);
void printMatrix(mat4x4_t matrix);

// Vector ops
vec3_t crossProduct(vec3_t a, vec3_t b);
float dot(vec3_t a, vec3_t b);
float magnitude(vec3_t v);
vec3_t sub(vec3_t a, vec3_t b);
vec3_t add(vec3_t a, vec3_t b);
vec3_t normalize(vec3_t v);
vec3_t mulScalarV3(float k, vec3_t v);

// Matrix ops
mat4x4_t transposeM4(mat4x4_t m);
mat4x4_t inverseM4(mat4x4_t matrix);
vec4_t mulMV4(mat4x4_t mat4x4, vec4_t vec4);
vec3_t mulMV3(mat4x4_t mat4x4, vec3_t v);
mat4x4_t mulMM4(mat4x4_t m1, mat4x4_t m2);

// Plane ops
float distancePlaneV3(plane_t plane, vec3_t v);

#endif