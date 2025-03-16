#ifndef OBJLOADER_H
#define OBJLOADER_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define STB_IMAGE_IMPLEMENTATION
#include "external/stb/stb_image.h"

#include "zerogl.h"

static inline zgl_canvas_t loadTexture(char* filename) {
    ZGL_DEBUG_PRINT("DEBUG: Loading texture %s\n", filename);
    zgl_canvas_t texture = {0};
    texture.frameBuffer = NULL;

    int channels;
    unsigned char* data = stbi_load(filename, &texture.width, &texture.height, &channels, 3);
    if (data == NULL) {
        fprintf(stderr, "ERROR: Couldn't load texture file %s\n", filename);
        exit(-1);
    }

    texture.frameBuffer = (uint32_t*) malloc(texture.width * texture.height * sizeof(uint32_t));
    if (texture.frameBuffer == NULL) {
        fprintf(stderr, "ERROR: Texture memory allocation failed.\n");
        exit(-1);
    }

    for (int i = 0; i < texture.width * texture.height; i++) {
        int bufferIndex = i * 3;
        uint8_t r = data[bufferIndex];
        uint8_t g = data[bufferIndex + 1];
        uint8_t b = data[bufferIndex + 2];
        texture.frameBuffer[i] = zgl_color(r, g, b);
    }

    stbi_image_free(data);
    ZGL_DEBUG_PRINT("DEBUG: Loaded texture %s\n", filename);
    return texture;
}


static inline char* getPath(const char* filename) {
    char* path = (char*) malloc(128 * sizeof(char));
    if (!path) {
        fprintf(stderr, "ERROR: Failed to allocate memory for path.\n");
        exit(1);
    }

    strcpy(path, filename);
    char* lastSlash = strrchr(path, '/');
    if (lastSlash) {
        *(lastSlash + 1) = '\0';
    } else {
        strcpy(path, "");
    }

    return path;
}

static inline zgl_material_t* loadMtlFile(const char* filename, int* numMaterials) {
    *numMaterials = 0;
    zgl_material_t* materials = NULL;

    char line[128];
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "WARN: Failed to open file %s.\n", filename);
        return NULL;
    }

    int hasKa = 0;
    int hasMapKa = 0;
    float r, g, b;
    while (fgets(line, 128, fp) != NULL) {
        if (line[0] == 'n' && line[1] == 'e' && line[2] == 'w') {
            // Cleanup last material
            if (*numMaterials > 0) {
                if (!hasKa) {
                    materials[*numMaterials - 1].ambientColor = materials[*numMaterials - 1].diffuseColor;
                }

                if (!hasMapKa) {
                    // Deep copy diffuse texture to ambient texture if no ambient texture is present
                    materials[*numMaterials - 1].ambientTexture.height = materials[*numMaterials - 1].diffuseTexture.height;
                    materials[*numMaterials - 1].ambientTexture.width = materials[*numMaterials - 1].diffuseTexture.width;
                    materials[*numMaterials - 1].ambientTexture.frameBuffer = (uint32_t*) malloc(materials[*numMaterials - 1].diffuseTexture.width * materials[*numMaterials - 1].diffuseTexture.height * sizeof(uint32_t));
                    memcpy(materials[*numMaterials - 1].ambientTexture.frameBuffer, materials[*numMaterials - 1].diffuseTexture.frameBuffer, materials[*numMaterials - 1].diffuseTexture.width * materials[*numMaterials - 1].diffuseTexture.height * sizeof(uint32_t));
                }
            }

            *numMaterials = *numMaterials + 1;
            materials = (zgl_material_t*) realloc(materials, *numMaterials * sizeof(zgl_material_t));
            if (materials == NULL) {
                fprintf(stderr, "ERROR: Material memory allocation failed.\n");
                exit(1);
            }

            char* name = (char*) malloc(128 * sizeof(char));
            sscanf(line, "newmtl %s\n", name);
            materials[*numMaterials - 1].name = name;
            materials[*numMaterials - 1].diffuseColor = zgl_color(255, 255, 255);
            materials[*numMaterials - 1].specularColor = zgl_color(0, 0, 0);
            materials[*numMaterials - 1].specularExponent = 10.0f;
            materials[*numMaterials - 1].diffuseTexture = (zgl_canvas_t) {NULL, 0, 0,};
            materials[*numMaterials - 1].specularTexture = (zgl_canvas_t) {NULL, 0, 0,};
        }

        if (line[0] == 'K' && line[1] == 'a') {
            sscanf(line, "Ka %f %f %f\n", &r, &g, &b);
            materials[*numMaterials - 1].ambientColor = zgl_color_from_floats(r, g, b);
            hasKa = 1;
        }

        if (line[0] == 'K' && line[1] == 'd') {
            sscanf(line, "Kd %f %f %f\n", &r, &g, &b);
            materials[*numMaterials - 1].diffuseColor = zgl_color_from_floats(r, g, b);
        }

        if (line[0] == 'K' && line[1] == 's') {
            sscanf(line, "Ks %f %f %f\n", &r, &g, &b);
            materials[*numMaterials - 1].specularColor = zgl_color_from_floats(r, g, b);
        }

        if (line[0] == 'N' && line[1] == 's') {
            sscanf(line, "Ns %f\n", &materials[*numMaterials - 1].specularExponent);
        }

        if (line[0] == 'm' && line[1] == 'a' && line[2] == 'p' && line[3] == '_' && line[4] == 'K' && line[5] == 'a') {
            char* textureFilename = (char*) malloc(128 * sizeof(char));
            sscanf(line, "map_Ka %s\n", textureFilename);
            char* path = getPath(filename);
            strcat(path, textureFilename);
            ZGL_DEBUG_PRINT("DEBUG: Loading ambient texture %s\n", textureFilename);
            materials[*numMaterials - 1].ambientTexture = loadTexture(path);
            hasMapKa = 1;
        }

        if (line[0] == 'm' && line[1] == 'a' && line[2] == 'p' && line[3] == '_' && line[4] == 'K' && line[5] == 'd') {
            char* textureFilename = (char*) malloc(128 * sizeof(char));
            sscanf(line, "map_Kd %s\n", textureFilename);
            char* path = getPath(filename);
            strcat(path, textureFilename);
            ZGL_DEBUG_PRINT("DEBUG: Loading diffuse texture %s\n", textureFilename);
            materials[*numMaterials - 1].diffuseTexture = loadTexture(path);
        }

        if (line[0] == 'm' && line[1] == 'a' && line[2] == 'p' && line[3] == '_' && line[4] == 'K' && line[5] == 's') {
            char* textureFilename = (char*) malloc(128 * sizeof(char));
            sscanf(line, "map_Ks %s\n", textureFilename);
            char* path = getPath(filename);
            strcat(path, textureFilename);
            ZGL_DEBUG_PRINT("DEBUG: Loading specular texture %s\n", textureFilename);
            materials[*numMaterials - 1].specularTexture = loadTexture(path);
        }
    }

    // Cleanup last material
    if (!hasKa) {
        materials[*numMaterials - 1].ambientColor = materials[*numMaterials - 1].diffuseColor;
    }

    if (!hasMapKa) {
        // Deep copy diffuse texture to ambient texture if no ambient texture is present
        materials[*numMaterials - 1].ambientTexture.height = materials[*numMaterials - 1].diffuseTexture.height;
        materials[*numMaterials - 1].ambientTexture.width = materials[*numMaterials - 1].diffuseTexture.width;
        materials[*numMaterials - 1].ambientTexture.frameBuffer = (uint32_t*) malloc(materials[*numMaterials - 1].diffuseTexture.width * materials[*numMaterials - 1].diffuseTexture.height * sizeof(uint32_t));
        memcpy(materials[*numMaterials - 1].ambientTexture.frameBuffer, materials[*numMaterials - 1].diffuseTexture.frameBuffer, materials[*numMaterials - 1].diffuseTexture.width * materials[*numMaterials - 1].diffuseTexture.height * sizeof(uint32_t));
    }

    fclose(fp);

    return materials;
}

static inline zgl_mesh_t* loadObjFile(const char* filename, bool flipTexturesVertically) {
    char name[128];

    zgl_vec3_t* vertices = NULL;
    int num_vertices = 0;

    zgl_vec3_t *textureCoords = NULL;
    int num_textureCoords = 0;

    zgl_vec3_t* normals = NULL;
    int num_normals = 0;

    zgl_triangle_t* triangles = NULL;
    int num_triangles = 0;

    zgl_material_t* materials = NULL;
    int num_materials = 0;

    int currentMaterial;

    char line[128];

    FILE* fp = fopen(filename, "r");

    while (fgets(line, 128, fp) != NULL) {
        if (line[0] == 'o') {
            sscanf(line, "o %s\n", name);
            ZGL_DEBUG_PRINT("DEBUG: Loading object %s.\n", name);
        }

        if (line[0] == 'v' && line[1] == ' ') {
            num_vertices++;
            vertices = (zgl_vec3_t*) realloc(vertices, num_vertices * sizeof(zgl_vec3_t));
            if (vertices == NULL) {
                fprintf(stderr, "ERROR: Vertex memory couldn't be allocated.\n");
                exit(-1);
            }
            sscanf(line, "v %f %f %f\n", &vertices[num_vertices - 1].x, &vertices[num_vertices - 1].y, &vertices[num_vertices - 1].z);
        }

        if (line[0] == 'v' && line[1] == 't') {
            num_textureCoords++;
            textureCoords = (zgl_vec3_t*) realloc(textureCoords, num_textureCoords * sizeof(zgl_vec3_t));
            if (textureCoords == NULL) {
                fprintf(stderr, "ERROR: Texture coordinate memory couldn't be allocated.\n");
                exit(-1);
            }
            float u, v;
            sscanf(line, "vt %f %f\n", &u, &v);
            textureCoords[num_textureCoords - 1].x = u;
            textureCoords[num_textureCoords - 1].y = flipTexturesVertically ? 1 - v : v;
            textureCoords[num_textureCoords - 1].z = 0.0f;
        }

        if (line[0] == 'v' && line[1] == 'n') {
            num_normals++;
            normals = (zgl_vec3_t*) realloc(normals, num_normals * sizeof(zgl_vec3_t));
            if (normals == NULL) {
                fprintf(stderr, "ERROR: Normal memory couldn't be allocated.\n");
                exit(-1);
            }
            sscanf(line, "vn %f %f %f\n", &normals[num_normals - 1].x, &normals[num_normals - 1].y, &normals[num_normals - 1].z);
        }

        if (line[0] == 'm' && line[1] == 't' && line[2] == 'l' && line[3] == 'l') {
            char mtl_filename[128];
            sscanf(line, "mtllib %s\n", mtl_filename);
            char* path = getPath(filename);
            strcat(path, mtl_filename);
            ZGL_DEBUG_PRINT("DEBUG: Loading MTL %s\n", path);
            materials = loadMtlFile(path, &num_materials);
            ZGL_DEBUG_PRINT("DEBUG: Loaded %d materials\n", num_materials);
        }

        if (line[0] == 'u' && line[1] == 's' && line[2] == 'e') {
            char material_name[128];
            sscanf(line, "usemtl %s\n", material_name);

            for (int i = 0; i < num_materials; i++) {
                if (strcmp(material_name, materials[i].name) == 0) {
                    ZGL_DEBUG_PRINT("DEBUG: Using material %s\n", materials[i].name);
                    uint8_t r, g, b;
                    zgl_color_components(materials[i].diffuseColor, &r, &g, &b);
                    ZGL_DEBUG_PRINT("DEBUG: Color %d %d %d\n", r, g, b);
                    currentMaterial = i;
                }
            }
        }

        if (line[0] == 'f' && line[1] == ' ') {
            num_triangles++;
            triangles = (zgl_triangle_t*) realloc(triangles, num_triangles * sizeof(zgl_triangle_t));
            if (triangles == NULL) {
                fprintf(stderr, "ERROR: Triangle memory couldn't be allocated.\n");
                exit(-1);
            }

            // Check how many slashes there are to determine the format of the face
            int num_slashes = 0;
            for (int i = 0; i < strlen(line); i++) {
                if (line[i] == '/') {
                    num_slashes++;
                }
            }

            int v0, v1, v2;
            int n0, n1, n2;
            int t0, t1, t2;
            if (num_slashes == 6) {
                sscanf(line, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", &v0, &t0, &n0, &v1, &t1, &n1, &v2, &t2, &n2);
                triangles[num_triangles - 1].v0 = v0 - 1;
                triangles[num_triangles - 1].v1 = v1 - 1;
                triangles[num_triangles - 1].v2 = v2 - 1;
                triangles[num_triangles - 1].t0 = t0 - 1;
                triangles[num_triangles - 1].t1 = t1 - 1;
                triangles[num_triangles - 1].t2 = t2 - 1;
                triangles[num_triangles - 1].n0 = n0 - 1;
                triangles[num_triangles - 1].n1 = n1 - 1;
                triangles[num_triangles - 1].n2 = n2 - 1;
            } else if (num_slashes == 3) {
                sscanf(line, "f %d/%d %d/%d %d/%d\n", &v0, &t0, &v1, &t1, &v2, &t2);
                triangles[num_triangles - 1].v0 = v0 - 1;
                triangles[num_triangles - 1].v1 = v1 - 1;
                triangles[num_triangles - 1].v2 = v2 - 1;
                triangles[num_triangles - 1].t0 = t0 - 1;
                triangles[num_triangles - 1].t1 = t1 - 1;
                triangles[num_triangles - 1].t2 = t2 - 1;

                // Compute normals
                zgl_vec3_t v0v1 = zgl_sub(vertices[v1 - 1], vertices[v0 - 1]);
                zgl_vec3_t v0v2 = zgl_sub(vertices[v2 - 1], vertices[v0 - 1]);
                zgl_vec3_t normal = zgl_cross(v0v1, v0v2);
                normal = zgl_normalize(normal);
                num_normals++;
                normals = (zgl_vec3_t*) realloc(normals, num_normals * sizeof(zgl_vec3_t));
                if (normals == NULL) {
                    fprintf(stderr, "ERROR: Normal memory couldn't be allocated.\n");
                    exit(-1);
                }
                normals[num_normals - 1] = normal;
                triangles[num_triangles - 1].n0 = num_normals - 1;
                triangles[num_triangles - 1].n1 = num_normals - 1;
                triangles[num_triangles - 1].n2 = num_normals - 1;
            } else if (num_slashes == 0) {
                sscanf(line, "f %d %d %d\n", &v0, &v1, &v2);
                triangles[num_triangles - 1].v0 = v0 - 1;
                triangles[num_triangles - 1].v1 = v1 - 1;
                triangles[num_triangles - 1].v2 = v2 - 1;

                // Compute normals
                zgl_vec3_t v0v1 = zgl_sub(vertices[v1 - 1], vertices[v0 - 1]);
                zgl_vec3_t v0v2 = zgl_sub(vertices[v2 - 1], vertices[v0 - 1]);
                zgl_vec3_t normal = zgl_cross(v0v1, v0v2);
                normal = zgl_normalize(normal);
                num_normals++;
                normals = (zgl_vec3_t*) realloc(normals, num_normals * sizeof(zgl_vec3_t));
                if (normals == NULL) {
                    fprintf(stderr, "ERROR: Normal memory couldn't be allocated.\n");
                    exit(-1);
                }
                ZGL_DEBUG_PRINT("DEBUG: Normal %d: ", num_normals - 1);
                normals[num_normals - 1] = normal;
                triangles[num_triangles - 1].n0 = num_normals - 1;
                triangles[num_triangles - 1].n1 = num_normals - 1;
                triangles[num_triangles - 1].n2 = num_normals - 1;
            }

            triangles[num_triangles - 1].materialIndex = currentMaterial;
        }
    }

    fclose(fp);

    zgl_mesh_t* mesh = (zgl_mesh_t*) malloc(sizeof(zgl_mesh_t));
    if (mesh == NULL) {
        fprintf(stderr, "ERROR: Mesh memory couldn't be allocated.\n");
        exit(-1);
    }

    float* invMagnitudeNormals = (float*) malloc(num_normals * sizeof(float));
    if (invMagnitudeNormals == NULL) {
        fprintf(stderr, "ERROR: Normal magnitudes couldn't be allocated.\n");
        exit(-1);
    }

    for (int i = 0; i < num_normals; i++) {
        invMagnitudeNormals[i] = 1.0f / zgl_magnitude(normals[i]);
    }

    mesh->numVertices = num_vertices;
    mesh->vertices = vertices;
    mesh->numTextureCoords = num_textureCoords;
    mesh->textureCoords = textureCoords;
    mesh->numNormals = num_normals;
    mesh->normals = normals;
    mesh->invMagnitudeNormals = invMagnitudeNormals;
    mesh->numTriangles = num_triangles;
    mesh->triangles = triangles;
    mesh->numMaterials = num_materials;
    mesh->materials = materials;

    mesh->center = zgl_mesh_center(vertices, num_vertices);
    mesh->boundsRadius = zgl_mesh_bound_radius(vertices, num_vertices, mesh->center);

    // Mesh name
    mesh->name = (char*) malloc(strlen(name) * sizeof(char) + 1);
    if (mesh->name == NULL) {
        fprintf(stderr, "ERROR: Mesh name memory couldn't be allocated.\n");
        exit(-1);
    }
    strcpy(mesh->name, name);

    ZGL_DEBUG_PRINT("DEBUG: Loaded mesh %s\n", mesh->name);
    return mesh;
}

#endif // OBJLOADER_H
