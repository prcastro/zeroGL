#ifndef OBJLOADER_H
#define OBJLOADER_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "simplerenderer.h"
#include "external/upng/upng.h"

static inline uint32_t* loadTexture(char* filename, int* textureWidth, int* textureHeight) {
    DEBUG_PRINT("DEBUG: Loading texture %s\n", filename);
    uint32_t* texture = NULL;

    upng_t* upng = upng_new_from_file(filename);
    if (upng  == NULL) {
        fprintf(stderr, "ERROR: Couldn't load texture file %s\n", filename);
        exit(-1);
    }

    upng_decode(upng);
    uint32_t  a, r, g, b;
    if (upng_get_error(upng) == UPNG_EOK) {
        *textureWidth = upng_get_width(upng);
        *textureHeight = upng_get_height(upng);

        DEBUG_PRINT("DEBUG: Texture size %d x %d\n", *textureWidth, *textureHeight);
        uint32_t* buffer = (uint32_t*) upng_get_buffer(upng);
        texture = (uint32_t*) malloc((*textureWidth) * (*textureHeight) * sizeof(uint32_t));
        if (texture == NULL) {
            fprintf(stderr, "ERROR: Texture memory allocation failed.\n");
            exit(-1);
        }
        for (int i = 0; i < (*textureWidth) * (*textureHeight); i++) {
            uint32_t pixel = buffer[i];
            a = (pixel & 0xFF000000) >> 24;
            b = (pixel & 0x00FF0000) >> 16;
            g = (pixel & 0x0000FF00) >> 8;
            r = (pixel & 0x000000FF);
            texture[i]  = (a << 24) | (r << 16) | (g << 8) |  b;
        }   
    } else {
        fprintf(stderr, "ERROR: Couldn't decode texture file\n");
        exit(-1);
    }

    upng_free(upng);
    DEBUG_PRINT("DEBUG: Loaded texture %s\n", filename);
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

static inline material_t* loadMtlFile(const char* filename, int* numMaterials) {
    *numMaterials = 0;
    material_t* materials = NULL;

    char line[128];
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "WARN: Failed to open file %s.\n", filename);
        return NULL;
    }

    float r, g, b;
    while (fgets(line, 128, fp) != NULL) {
        if (line[0] == 'n' && line[1] == 'e' && line[2] == 'w') {
            *numMaterials = *numMaterials + 1;
            materials = (material_t*) realloc(materials, *numMaterials * sizeof(material_t));
            if (materials == NULL) {
                fprintf(stderr, "ERROR: Material memory allocation failed.\n");
                exit(1);
            }

            char* name = (char*) malloc(128 * sizeof(char));
            sscanf(line, "newmtl %s\n", name);
            materials[*numMaterials - 1].name = name;
            
            materials[*numMaterials - 1].texture = NULL;
            materials[*numMaterials - 1].textureWidth = 0;
            materials[*numMaterials - 1].textureHeight = 0;
        }

        if (line[0] == 'K' && line[1] == 'd') {
            sscanf(line, "Kd %f %f %f\n", &r, &g, &b);
            materials[*numMaterials - 1].diffuseColor = colorFromFloats(r, g, b);
        }

        if (line[0] == 'K' && line[1] == 's') {
            sscanf(line, "Ks %f %f %f\n", &r, &g, &b);
            materials[*numMaterials - 1].specularColor = colorFromFloats(r, g, b);
        }

        if (line[0] == 'N' && line[1] == 's') {
            sscanf(line, "Ns %f\n", &materials[*numMaterials - 1].specularExponent);
        }

        
        if (line[0] == 'm' && line[1] == 'a' && line[2] == 'p' && line[3] == '_' && line[4] == 'K' && line[5] == 'd') {
            char* textureFilename = (char*) malloc(128 * sizeof(char));
            sscanf(line, "map_Kd %s\n", textureFilename);
            char* path = getPath(filename);
            strcat(path, textureFilename);
            DEBUG_PRINT("DEBUG: Loading texture %s\n", textureFilename);
            materials[*numMaterials - 1].texture = loadTexture(path, &materials[*numMaterials - 1].textureWidth, &materials[*numMaterials - 1].textureHeight);
        }
    }

    fclose(fp);

    return materials;
}

static inline mesh_t* loadObjFile(const char* filename, bool flipTexturesVertically) {
    char name[128];

    vec3_t* vertices = NULL;
    int num_vertices = 0;

    vec3_t *textureCoords = NULL;
    int num_textureCoords = 0;

    vec3_t* normals = NULL;
    int num_normals = 0;

    triangle_t* triangles = NULL;
    int num_triangles = 0;

    material_t* materials = NULL;
    int num_materials = 0;

    int currentMaterial;

    char line[128];
    
    FILE* fp = fopen(filename, "r");

    while (fgets(line, 128, fp) != NULL) {
        if (line[0] == 'o') {
            sscanf(line, "o %s\n", name);
            DEBUG_PRINT("DEBUG: Loading object %s.\n", name);
        }

        if (line[0] == 'v' && line[1] == ' ') {
            num_vertices++;
            vertices = (vec3_t*) realloc(vertices, num_vertices * sizeof(vec3_t));
            if (vertices == NULL) {
                fprintf(stderr, "ERROR: Vertex memory couldn't be allocated.\n");
                exit(-1);
            }
            sscanf(line, "v %f %f %f\n", &vertices[num_vertices - 1].x, &vertices[num_vertices - 1].y, &vertices[num_vertices - 1].z);
        }

        if (line[0] == 'v' && line[1] == 't') {
            num_textureCoords++;
            textureCoords = (vec3_t*) realloc(textureCoords, num_textureCoords * sizeof(vec3_t));
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
            normals = (vec3_t*) realloc(normals, num_normals * sizeof(vec3_t));
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
            DEBUG_PRINT("DEBUG: Loading MTL %s\n", path);
            materials = loadMtlFile(path, &num_materials);
            DEBUG_PRINT("DEBUG: Loaded %d materials\n", num_materials);
        }

        if (line[0] == 'u' && line[1] == 's' && line[2] == 'e') {
            char material_name[128];
            sscanf(line, "usemtl %s\n", material_name);

            for (int i = 0; i < num_materials; i++) {
                if (strcmp(material_name, materials[i].name) == 0) {
                    DEBUG_PRINT("DEBUG: Using material %s\n", materials[i].name);
                    uint8_t r, g, b;
                    colorFromUint32(materials[i].diffuseColor, &r, &g, &b);
                    DEBUG_PRINT("DEBUG: Color %d %d %d\n", r, g, b);
                    currentMaterial = i;
                }
            }
        }

        if (line[0] == 'f' && line[1] == ' ') {
            num_triangles++;
            triangles = (triangle_t*) realloc(triangles, num_triangles * sizeof(triangle_t));
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
                vec3_t v0v1 = sub(vertices[v1 - 1], vertices[v0 - 1]);
                vec3_t v0v2 = sub(vertices[v2 - 1], vertices[v0 - 1]);
                vec3_t normal = crossProduct(v0v1, v0v2);
                normal = normalize(normal);
                num_normals++;
                normals = (vec3_t*) realloc(normals, num_normals * sizeof(vec3_t));
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
                vec3_t v0v1 = sub(vertices[v1 - 1], vertices[v0 - 1]);
                vec3_t v0v2 = sub(vertices[v2 - 1], vertices[v0 - 1]);
                vec3_t normal = crossProduct(v0v1, v0v2);
                normal = normalize(normal);
                num_normals++;
                normals = (vec3_t*) realloc(normals, num_normals * sizeof(vec3_t));
                if (normals == NULL) {
                    fprintf(stderr, "ERROR: Normal memory couldn't be allocated.\n");
                    exit(-1);
                }
                DEBUG_PRINT("DEBUG: Normal %d: ", num_normals - 1);
                normals[num_normals - 1] = normal;
                triangles[num_triangles - 1].n0 = num_normals - 1;
                triangles[num_triangles - 1].n1 = num_normals - 1;
                triangles[num_triangles - 1].n2 = num_normals - 1;
            }
            
            triangles[num_triangles - 1].materialIndex = currentMaterial;
        }
    }

    fclose(fp);

    mesh_t* mesh = (mesh_t*) malloc(sizeof(mesh_t));
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
        invMagnitudeNormals[i] = 1.0f / magnitude(normals[i]);
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

    mesh->center = meshCenter(vertices, num_vertices);
    mesh->boundsRadius = meshBoundsRadius(vertices, num_vertices, mesh->center);

    // Mesh name
    mesh->name = (char*) malloc(strlen(name) * sizeof(char) + 1);
    if (mesh->name == NULL) {
        fprintf(stderr, "ERROR: Mesh name memory couldn't be allocated.\n");
        exit(-1);
    }
    strcpy(mesh->name, name);
    
    DEBUG_PRINT("DEBUG: Loaded mesh %s\n", mesh->name);
    return mesh;
}

#endif // OBJLOADER_H
