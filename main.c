// #define ZGL_DEBUG
#define DEBUGUI
#define SDL_MAIN_HANDLED

#ifdef DEBUGUI
#define NK_INCLUDE_FIXED_TYPES
#define NK_INCLUDE_STANDARD_IO
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_STANDARD_BOOL
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_INCLUDE_VERTEX_BUFFER_OUTPUT
#define NK_INCLUDE_FONT_BAKING
#define NK_INCLUDE_DEFAULT_FONT
#define NK_IMPLEMENTATION
#define NK_SDL_RENDERER_IMPLEMENTATION
#include "external/nuklear/nuklear.h"
#include "external/nuklear/nuklear_sdl_renderer.h"
#endif // DEBUGUI

#include <SDL.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define ZEROGL_IMPLEMENTATION
#include "zerogl.h"
#include "objloader.h"


#define WIDTH 1066
#define HEIGHT 600
#define ROTATION_SPEED 0.0f //15.0f // degrees per second
#define VIEWPORT_WIDTH (WIDTH /(float) HEIGHT)
#define VIEWPORT_HEIGHT 1.0f
#define VIEWPORT_DISTANCE 1.0f
#define PIXEL_DEPTH 4
#define PITCH (PIXEL_DEPTH * WIDTH)

typedef enum {
    BASIC_SHADER,
    COLORED_SHADER,
    UNLIT_SHADER,
    FLAT_SHADER,
    GOURAUD_SHADER,
    PHONG_SHADER,
} shader_type_t;

typedef struct game_state_t {
    // Loop control
    int          running;
    double       elapsedTime;
    uint64_t     lastTime;

    // Rendering
    SDL_Event     event;
    SDL_Window*   window;
    SDL_Renderer* renderer;
    SDL_Texture*  texture;
    zgl_canvas_t  canvas;
    uint32_t      backgroundColor;
    uint16_t      renderOptions;
    int           drawLights;
    int           draw3DObjects;
    int           bilinearFiltering;
    shader_type_t shaderType;

    // Game objects
    int                  numMeshes;
    zgl_mesh_t*          meshes;
    int                  numObjects;
    zgl_object3D_t*      objects;
    zgl_light_sources_t  lightSources;
    zgl_object3D_t*      pointLightObjects;
    zgl_camera_t         camera;
    float                rotationSpeed;
    const uint8_t*       keys;

    // GUI
    #ifdef DEBUGUI
    struct nk_context* nuklearContext;
    int                showGUI;
    int                toggleGUIKeyPressed;
    #endif // DEBUG
} game_state_t;

game_state_t* init() {
    ZGL_DEBUG_PRINT("INFO: Initializing game objects\n");

    ZGL_DEBUG_PRINT("INFO: Initializing SDL\n");
    if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
        fprintf(stderr, "ERROR: error initializing SDL: %s\n", SDL_GetError());
        exit(-1);
    }

    ZGL_DEBUG_PRINT("INFO: Loading meshes and objects\n");
    int numObjects = 1;
    int numMeshes = 7;
    zgl_mesh_t* meshes = (zgl_mesh_t*) malloc(numMeshes * sizeof(zgl_mesh_t));
    zgl_object3D_t *objects = (zgl_object3D_t*) malloc(numObjects * sizeof(zgl_object3D_t));
    if (meshes == NULL || objects == NULL) {
        fprintf(stderr, "ERROR: 3D objects memory couldn't be allocated.\n");
        exit(-1);
    }

    // Define a debug mesh
    zgl_vec3_t* vertices = (zgl_vec3_t*) malloc(3 * sizeof(zgl_vec3_t));
    zgl_triangle_t* triangles = (zgl_triangle_t*) malloc(1 * sizeof(zgl_triangle_t));
    zgl_material_t* materials = (zgl_material_t*) malloc(1 * sizeof(zgl_material_t));
    if (vertices == NULL || triangles == NULL || materials == NULL) {
        fprintf(stderr, "ERROR: Debug mesh memory couldn't be allocated.\n");
        exit(-1);
    }

    vertices[0] = (zgl_vec3_t) {0, 0, 0};
    vertices[1] = (zgl_vec3_t) {1, 0, 0};
    vertices[2] = (zgl_vec3_t) {0, 1, 0};
    triangles[0] = (zgl_triangle_t) {0, 2, 1, 0, 0, 0, 0, 0, 0, 0};
    materials[0] = (zgl_material_t) {"RedMaterial", ZGL_COLOR_RED, ZGL_COLOR_RED, ZGL_COLOR_RED, 0.0f, (zgl_canvas_t) {NULL, 0, 0, 0, NULL}};
    meshes[0] = (zgl_mesh_t) {
        .name = "Debug",
        .numVertices = 3,
        .numTriangles = 1,
        .numTextureCoords = 0,
        .vertices = vertices,
        .triangles = triangles,
        .materials = materials
    };

    meshes[1] = *loadObjFile("assets/light.obj", false);
    meshes[2] = *loadObjFile("assets/snake/snake.obj", true);
    meshes[3] = *loadObjFile("assets/engineer/engineer.obj", false);
    meshes[4] = *loadObjFile("assets/cube.obj", false);
    meshes[5] = *loadObjFile("assets/woodcube/woodcube.obj", false);
    meshes[6] = *loadObjFile("assets/sphere.obj", false);

    objects[0] = zgl_object(&meshes[3], (zgl_vec3_t) {0, 0, 0}, 1.0 , IDENTITY_M4x4);

    ZGL_DEBUG_PRINT("INFO: Loading lights\n");
    int numAmbientLights = 1;
    int numDirLights = 1;
    int numPointLights = 1;
    zgl_ambient_light_t* ambientLights = (zgl_ambient_light_t*) malloc(numAmbientLights * sizeof(zgl_ambient_light_t));
    zgl_dir_light_t* directionalLights = (zgl_dir_light_t*) malloc(numDirLights * sizeof(zgl_dir_light_t));
    zgl_point_light_t* pointLights = (zgl_point_light_t*) malloc(numPointLights * sizeof(zgl_point_light_t));
    zgl_object3D_t *pointLightObjects = (zgl_object3D_t*) malloc(numPointLights * sizeof(zgl_object3D_t));
    if (pointLights == NULL || directionalLights == NULL || ambientLights == NULL || pointLightObjects == NULL) {
        fprintf(stderr, "ERROR: Lights memory couldn't be allocated.\n");
        exit(-1);
    }

    ambientLights[0] = (zgl_ambient_light_t) {
        .intensity = {0.4, 0.4, 0.4}
    };
    directionalLights[0] = (zgl_dir_light_t) {{0.0, 0.0, 0.0}, {0.0, -1.0, 1.0}};
    pointLights[0] = (zgl_point_light_t) {{0.9, 0.9, 0.9}, {-0.5, 1.5, -2.0}};

    pointLightObjects[0] = zgl_object(&meshes[1], pointLights[0].position, 0.05, IDENTITY_M4x4);

    ZGL_DEBUG_PRINT("INFO: Initializing game state\n");
    uint32_t *frameBuffer = (uint32_t*) malloc(WIDTH * HEIGHT * sizeof(uint32_t));
    float *depthBuffer = (float*) malloc(WIDTH * HEIGHT * sizeof(float));
    game_state_t* game = (game_state_t*) malloc(sizeof(game_state_t));
    if (game == NULL || frameBuffer == NULL || depthBuffer == NULL) {
        fprintf(stderr, "ERROR: Game state memory and buffers couldn't be allocated.\n");
        exit(-1);
    }

    zgl_canvas_t canvas = {
        .width = WIDTH,
        .height = HEIGHT,
        .hasDepthBuffer = 1,
        .frameBuffer = frameBuffer,
        .depthBuffer = depthBuffer
    };

    game->running = 1;
    game->elapsedTime = 0;
    game->lastTime = SDL_GetPerformanceCounter();
    game->window = SDL_CreateWindow("Rasterizer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
    game->renderer = SDL_CreateRenderer(game->window, -1, 0);
    game->texture = SDL_CreateTexture(game->renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);
    game->canvas = canvas;
    game->backgroundColor = ZGL_COLOR_BLACK;
    game->drawLights    = 1;
    game->draw3DObjects = 1;
    game->bilinearFiltering = 0;
    game->shaderType = GOURAUD_SHADER;
    game->renderOptions = ZGL_DIFFUSE_LIGHTING | ZGL_SPECULAR_LIGHTING | ZGL_BACKFACE_CULLING | ZGL_FUSTRUM_CULLING;
    game->numMeshes = numMeshes;
    game->meshes = meshes;
    game->numObjects = numObjects;
    game->objects = objects;

    // Lights
    game->lightSources = (zgl_light_sources_t) {
        .ambientLights = ambientLights,
        .numAmbientLights = numAmbientLights,
        .directionalLights = directionalLights,
        .numDirectionalLights = numDirLights,
        .pointLights = pointLights,
        .numPointLights = numPointLights
    };

    game->pointLightObjects = pointLightObjects;

    game->camera = zgl_camera(
        (zgl_vec3_t) {0, 0, -5}, // position
        (zgl_vec3_t) {0, 0, 1}, // direction
        (zgl_vec3_t) {0, 1, 0}, // up
        53.0f, // fov
        VIEWPORT_WIDTH / VIEWPORT_HEIGHT, // aspect ratio
        VIEWPORT_DISTANCE, // near plane
        100.0f, // far plane
        5.0f, // movement speed
        90.0f // turning speed
    );

    game->rotationSpeed = ROTATION_SPEED;
    game->keys = SDL_GetKeyboardState(NULL);

    #ifdef DEBUGUI
    ZGL_DEBUG_PRINT("INFO:  Initializing Dear ImGui\n");
    struct nk_context *ctx = nk_sdl_init(game->window, game->renderer);
    struct nk_font_atlas *atlas;
    struct nk_font_config config = nk_font_config(0);
    struct nk_font *font;
    nk_sdl_font_stash_begin(&atlas);
    font = nk_font_atlas_add_default(atlas, 10.0, &config);
    nk_sdl_font_stash_end();
    nk_style_set_font(ctx, &font->handle);
    game->nuklearContext = ctx;
    game->showGUI = 1;
    game->toggleGUIKeyPressed = 0;
    #endif // DEBUGUI

    return game;
}

void handleEvents(game_state_t* game) {
    ZGL_DEBUG_PRINT("INFO: Handle events\n");

    #ifdef DEBUGUI
    nk_input_begin(game->nuklearContext);
    #endif // DEBUGUI

    while (SDL_PollEvent(&game->event)) {
        #ifdef DEBUGUI
        nk_sdl_handle_event(&game->event);
        #endif // DEBUGUI
        switch (game->event.type) {
        case SDL_QUIT:
            // handling of close button
            ZGL_DEBUG_PRINT("INFO: Quitting application\n");
            game->running = 0;
            break;
        }
    }
}

#ifdef DEBUGUI
void updateDebugUI(game_state_t *game) {
    // Check for toggle key
    if (game->keys[SDL_SCANCODE_SPACE]) {
        if (!game->toggleGUIKeyPressed) {
            game->toggleGUIKeyPressed = 1;
            game->showGUI = !game->showGUI; // Toggle state of ImGui display
        }
    } else {
        game->toggleGUIKeyPressed = 0;
    }


    if (game->showGUI) {
        nk_window_show(game->nuklearContext, "Settings", NK_SHOWN);
    }


    if (game->showGUI) {
        ZGL_DEBUG_PRINT("INFO: Updating GUI\n");
        int row_size = 12;
        struct nk_context *ctx = game->nuklearContext;
        if (nk_begin(ctx, "Settings", nk_rect(0, 0, 220, HEIGHT),
                     NK_WINDOW_BORDER|NK_WINDOW_MOVABLE|NK_WINDOW_CLOSABLE|NK_WINDOW_SCALABLE)) {

            nk_layout_row_static(ctx, row_size, 150, 1);
            nk_labelf(ctx, NK_TEXT_LEFT, "FPS: %.2f (%.3lf ms)", floor(1000.0f / game->elapsedTime), game->elapsedTime);

            if (nk_tree_push(ctx, NK_TREE_NODE, "Meshes", NK_MINIMIZED)) {
                for (int i = 0; i < game->numMeshes; i++) {
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_labelf(ctx, NK_TEXT_LEFT, "%s: (vertices %d, triangles: %d, uvs: %d)", game->meshes[i].name, game->meshes[i].numVertices, game->meshes[i].numTriangles, game->meshes[i].numTextureCoords);
                }
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Objects", NK_MAXIMIZED)) {
                for (int i = 0; i < game->numObjects; i++) {
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_labelf(ctx, NK_TEXT_LEFT, "%d: %s at (%.0f, %.0f, %.0f)", i, game->objects[i].mesh->name, game->objects[i].translation.x, game->objects[i].translation.y, game->objects[i].translation.z);
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_property_float(ctx, "x", -10.0f, &game->objects[i].translation.x, 10.0f, 0.1f, 0.1f);
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_property_float(ctx, "y", -10.0f, &game->objects[i].translation.y, 10.0f, 0.1f, 0.1f);
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_property_float(ctx, "z", -10.0f, &game->objects[i].translation.z, 10.0f, 0.1f, 0.1f);
                    nk_layout_row_dynamic(ctx, row_size, 1);
                    nk_property_float(ctx, "scale", -0.01f, &game->objects[i].scale, 10.0f, 0.01f, 0.01f);
                }
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Lights", NK_MINIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 2);

                nk_bool isDiffuse = game->renderOptions & ZGL_DIFFUSE_LIGHTING;
                nk_checkbox_label(ctx, "Difuse", &isDiffuse);
                game->renderOptions = isDiffuse ? game->renderOptions | ZGL_DIFFUSE_LIGHTING : game->renderOptions & ~ZGL_DIFFUSE_LIGHTING;

                nk_bool isSpecular = game->renderOptions & ZGL_SPECULAR_LIGHTING;
                nk_checkbox_label(ctx, "Specular", &isSpecular);
                game->renderOptions = isSpecular ? game->renderOptions | ZGL_SPECULAR_LIGHTING : game->renderOptions & ~ZGL_SPECULAR_LIGHTING;

                if (nk_tree_push(ctx, NK_TREE_NODE, "Ambient", NK_MAXIMIZED)) {
                    for (int i = 0; i < game->lightSources.numAmbientLights; i++) {
                        nk_layout_row_dynamic(ctx, row_size * 7, 1);
                        char label[100];
                        sprintf(label, "Ambient Light %d", i);
                        if (nk_group_begin(ctx, label, NK_WINDOW_TITLE|NK_WINDOW_BORDER|NK_WINDOW_NO_SCROLLBAR)) {
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity R", 0.0f, &game->lightSources.ambientLights[i].intensity.x, 3.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity G", 0.0f, &game->lightSources.ambientLights[i].intensity.y, 3.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity B", 0.0f, &game->lightSources.ambientLights[i].intensity.z, 3.0f, 0.1f, 0.1f);
                            nk_group_end(ctx);
                        }
                    }
                    nk_tree_pop(ctx);
                }

                if (nk_tree_push(ctx, NK_TREE_NODE, "Directional", NK_MAXIMIZED)) {
                    for (int i = 0; i < game->lightSources.numDirectionalLights; i++) {
                        nk_layout_row_dynamic(ctx, row_size * 11, 1);
                        char label[100];
                        sprintf(label, "Directional Light %d", i);
                        if (nk_group_begin(ctx, label, NK_WINDOW_TITLE|NK_WINDOW_BORDER|NK_WINDOW_NO_SCROLLBAR)) {
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity R", 0.0f, &game->lightSources.directionalLights[i].intensity.x, 3.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity G", 0.0f, &game->lightSources.directionalLights[i].intensity.y, 3.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity B", 0.0f, &game->lightSources.directionalLights[i].intensity.z, 3.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "x", -1.0f, &game->lightSources.directionalLights[i].direction.x, 1.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "y", -1.0f, &game->lightSources.directionalLights[i].direction.y, 1.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "z", -1.0f, &game->lightSources.directionalLights[i].direction.z, 1.0f, 0.1f, 0.1f);
                            nk_group_end(ctx);
                        }
                    }
                    nk_tree_pop(ctx);
                }

                if (nk_tree_push(ctx, NK_TREE_NODE, "Point", NK_MAXIMIZED)) {
                    for (int i = 0; i < game->lightSources.numPointLights; i++) {
                        nk_layout_row_dynamic(ctx, row_size * 11, 1);
                        char label[100];
                        sprintf(label, "Point Light %d", i);
                        if (nk_group_begin(ctx, label, NK_WINDOW_TITLE|NK_WINDOW_BORDER|NK_WINDOW_NO_SCROLLBAR)) {
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity R", 0.0f, &game->lightSources.pointLights[i].intensity.x, 3.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity G", 0.0f, &game->lightSources.pointLights[i].intensity.y, 3.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity B", 0.0f, &game->lightSources.pointLights[i].intensity.z, 3.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "x", -10.0f, &game->lightSources.pointLights[i].position.x, 10.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "y", -10.0f, &game->lightSources.pointLights[i].position.y, 10.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "z", -10.0f, &game->lightSources.pointLights[i].position.z, 10.0f, 0.1f, 0.1f);
                            game->pointLightObjects[i] = zgl_object(&game->meshes[0], game->lightSources.pointLights[i].position, 0.05, IDENTITY_M4x4);
                            nk_group_end(ctx);
                        }
                    }
                    nk_tree_pop(ctx);
                }
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Camera", NK_MINIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "Camera: (%.1f, %.1f, %.1f)", game->camera.position.x, game->camera.position.y, game->camera.position.z);
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "Direction: (%.1f, %.1f, %.1f)", game->camera.direction.x, game->camera.direction.y, game->camera.direction.z);
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "Up: (%.1f, %.1f, %.1f)", game->camera.up.x, game->camera.up.y, game->camera.up.z);
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "FOV: %.1f", game->camera.fov);
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "Aspect ratio: %.1f", game->camera.aspectRatio);
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "Near plane: %.1f", game->camera.nearPlane);
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "Far plane: %.1f", game->camera.farPlane);
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "Movement speed: %.1f", game->camera.movementSpeed);
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_labelf(ctx, NK_TEXT_LEFT, "Turning speed: %.1f", game->camera.turningSpeed);
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Scene", NK_MAXIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_label(ctx, "What to draw", NK_TEXT_LEFT);
                nk_layout_row_dynamic(ctx, row_size, 3);

                nk_bool draw3DObjects = game->draw3DObjects;
                nk_checkbox_label(ctx, "3D Obj", &draw3DObjects);
                game->draw3DObjects = draw3DObjects;

                nk_bool drawLights = game->drawLights;
                nk_checkbox_label(ctx, "Lights", &drawLights);
                game->drawLights = drawLights;

                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Render Options", NK_MAXIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 4);
                nk_bool isBasic = game->shaderType == BASIC_SHADER;
                if (nk_radio_label(ctx, "Basic", &isBasic)) {
                    game->shaderType = BASIC_SHADER;
                };

                nk_bool isColored = game->shaderType == COLORED_SHADER;
                if (nk_radio_label(ctx, "Colored", &isColored)) {
                    game->shaderType = COLORED_SHADER;
                };

                nk_bool isUnlit = game->shaderType == UNLIT_SHADER;
                if (nk_radio_label(ctx, "Unlit", &isUnlit)) {
                    game->shaderType = UNLIT_SHADER;
                };

                nk_bool isFlat = game->shaderType == FLAT_SHADER;
                if (nk_radio_label(ctx, "Flat", &isFlat)) {
                    game->shaderType = FLAT_SHADER;
                };

                nk_bool isGouraud = game->shaderType == GOURAUD_SHADER;
                if (nk_radio_label(ctx, "Gouraud", &isGouraud)) {
                    game->shaderType = GOURAUD_SHADER;
                };

                nk_bool isPhong = game->shaderType == PHONG_SHADER;
                if (nk_radio_label(ctx, "Phong", &isPhong)) {
                    game->shaderType = PHONG_SHADER;
                };

                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_bool isBackfaceCulling = game->renderOptions & ZGL_BACKFACE_CULLING;
                nk_checkbox_label(ctx, "Backface culling", &isBackfaceCulling);
                game->renderOptions = isBackfaceCulling ? game->renderOptions | ZGL_BACKFACE_CULLING : game->renderOptions & ~ZGL_BACKFACE_CULLING;

                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_bool isFustrumCulling = game->renderOptions & ZGL_FUSTRUM_CULLING;
                nk_checkbox_label(ctx, "Fustrum culling", &isFustrumCulling);
                game->renderOptions = isFustrumCulling ? game->renderOptions | ZGL_FUSTRUM_CULLING : game->renderOptions & ~ZGL_FUSTRUM_CULLING;

                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_bool isBilinearFiltering = game->bilinearFiltering;
                nk_checkbox_label(ctx, "Bilinear filtering", &isBilinearFiltering);
                game->bilinearFiltering = isBilinearFiltering;

                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Background", NK_MAXIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_label(ctx, "Color", NK_TEXT_LEFT);
                nk_layout_row_dynamic(ctx, row_size * 10, 1);
                uint8_t r, g, b;
                zgl_color_components(game->backgroundColor, &r, &g, &b);
                struct nk_colorf nkBackgroundColor = nk_color_cf(nk_rgb(r, g, b));
                nk_color_pick(ctx, &nkBackgroundColor, NK_RGBA);
                game->backgroundColor = zgl_color_from_floats(nkBackgroundColor.r, nkBackgroundColor.g, nkBackgroundColor.b);
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Animation", NK_MAXIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_property_float(ctx, "Rotation speed", -360.0f, &game->rotationSpeed, 360.0f, 0.1f, 0.1f);
                nk_tree_pop(ctx);
            }
        }

        nk_end(game->nuklearContext);
    }

    if (game->showGUI && nk_window_is_hidden(game->nuklearContext, "Settings")) {
        game->showGUI = 0;
    }
}
#endif // DEBUGUI

void updateCameraPosition(game_state_t* game) {
    const uint8_t* keys = game->keys;
    zgl_camera_t *camera = &game->camera;

    // TODO: Store this in the camera struct
    zgl_vec3_t cameraRight = zgl_cross(camera->up, camera->direction);

    float elapsedTime = game->elapsedTime / 1000.0f;
    float movementSpeed = camera->movementSpeed * elapsedTime;
    float turningSpeed = camera->turningSpeed * elapsedTime;
    zgl_vec3_t newCameraPosition = camera->position;
    zgl_vec3_t newCameraDirection = camera->direction;
    zgl_vec3_t newCameraUp = camera->up;


    if (keys[SDL_SCANCODE_A]) {
        // Translate left
        newCameraPosition = zgl_add(newCameraPosition, zgl_mul_scalar(-movementSpeed, cameraRight));
    }

    if (keys[SDL_SCANCODE_D]) {
        // Translate right
        newCameraPosition = zgl_add(newCameraPosition, zgl_mul_scalar(movementSpeed, cameraRight));
    }

    if (keys[SDL_SCANCODE_S]) {
        // Translate back
        newCameraPosition = zgl_add(newCameraPosition, zgl_mul_scalar(-movementSpeed, newCameraDirection));
    }

    if (keys[SDL_SCANCODE_W]) {
        // Translate forward
        newCameraPosition = zgl_add(newCameraPosition, zgl_mul_scalar(movementSpeed, newCameraDirection));
    }

    if (keys[SDL_SCANCODE_PAGEDOWN]) {
        // Translate down
        newCameraPosition = zgl_add(newCameraPosition, zgl_mul_scalar(-movementSpeed, newCameraUp));
    }

    if (keys[SDL_SCANCODE_PAGEUP]) {
        // Translate up
        newCameraPosition = zgl_add(newCameraPosition, zgl_mul_scalar(movementSpeed, newCameraUp));
    }

    if (keys[SDL_SCANCODE_RIGHT]) {
        // Rotate right around up axis using quaternion
        zgl_quaternion_t rotation = zgl_quaternion(turningSpeed, newCameraUp);
        newCameraDirection = zgl_rotate(newCameraDirection, rotation);
        cameraRight = zgl_cross(newCameraUp, newCameraDirection);
    }

    if (keys[SDL_SCANCODE_LEFT]) {
        // Rotate left around up axis
        zgl_quaternion_t rotation = zgl_quaternion(-turningSpeed, newCameraUp);
        newCameraDirection = zgl_rotate(newCameraDirection, rotation);
        cameraRight = zgl_cross(newCameraUp, newCameraDirection);
    }

    if (keys[SDL_SCANCODE_UP]) {
        // Rotate up around right axis
        zgl_quaternion_t rotation = zgl_quaternion(-turningSpeed, cameraRight);
        newCameraDirection = zgl_rotate(newCameraDirection, rotation);
        newCameraUp = zgl_rotate(newCameraUp, rotation);
    }

    if (keys[SDL_SCANCODE_DOWN]) {
        // Rotate down around right axis
        zgl_quaternion_t rotation = zgl_quaternion(turningSpeed, cameraRight);
        newCameraDirection = zgl_rotate(newCameraDirection, rotation);
        newCameraUp = zgl_rotate(newCameraUp, rotation);
    }

    if (keys[SDL_SCANCODE_Q]) {
        // Rotate left around direction axis
        zgl_quaternion_t rotation = zgl_quaternion(turningSpeed, newCameraDirection);
        newCameraUp = zgl_rotate(newCameraUp, rotation);
        cameraRight = zgl_cross(newCameraUp, newCameraDirection);
    }

    if (keys[SDL_SCANCODE_E]) {
        // Rotate right around direction axis
        zgl_quaternion_t rotation = zgl_quaternion(-turningSpeed, newCameraDirection);
        newCameraUp = zgl_rotate(newCameraUp, rotation);
        cameraRight = zgl_cross(newCameraUp, newCameraDirection);
    }

    game->camera = zgl_camera(
        newCameraPosition,
        newCameraDirection,
        newCameraUp,
        camera->fov,
        camera->aspectRatio,
        camera->nearPlane,
        camera->farPlane,
        camera->movementSpeed,
        camera->turningSpeed
    );
}

zgl_object3D_t rotateObjectY(zgl_object3D_t object, float degrees) {
    return zgl_object(object.mesh, object.translation, object.scale, zgl_mul_mat(zgl_roty_mat(degrees), object.rotation));
}

// TODO: Use quaternions for rotation
void animateObjects(game_state_t* game) {
    float degrees = game->rotationSpeed * (game->elapsedTime / 1000.0f);
    game->objects[0] = rotateObjectY(game->objects[0], fmod(degrees, 360.0f));
}

void update(game_state_t* game) {
    ZGL_DEBUG_PRINT("INFO: Update game state\n");
    updateCameraPosition(game);
    animateObjects(game);
}

void drawObjects(game_state_t* game) {
    for (int i = 0; i < game->numObjects; i++) {
        zgl_object3D_t object = game->objects[i];
        switch (game->shaderType) {
            case BASIC_SHADER:
                zgl_basic_uniform_t basicUniformData = {
                    .modelviewprojection = zgl_mul_mat(game->camera.viewProjMatrix, object.transform),
                };
                zgl_render_object3D(&object, &basicUniformData, game->camera, game->canvas, zgl_basic_vertex_shader, zgl_basic_fragment_shader, game->renderOptions);
            case COLORED_SHADER:
                zgl_colored_uniform_t uniformData = {
                    .modelviewprojection = zgl_mul_mat(game->camera.viewProjMatrix, object.transform),
                    .materials = object.mesh->materials
                };
                zgl_render_object3D(&object, &uniformData, game->camera, game->canvas, zgl_colored_vertex_shader, zgl_colored_fragment_shader, game->renderOptions);
                break;
            case UNLIT_SHADER:
                zgl_unlit_uniform_t unlitUniformData = {
                    .modelMatrix = object.transform,
                    .modelInvRotationMatrixTransposed = zgl_transpose(zgl_inverse(object.rotation)),
                    .viewProjectionMatrix = game->camera.viewProjMatrix,
                    .bilinearFiltering = game->bilinearFiltering,
                    .materials = object.mesh->materials
                };
                zgl_render_object3D(&object, &unlitUniformData, game->camera, game->canvas, zgl_unlit_vertex_shader, zgl_unlit_fragment_shader, game->renderOptions);
                break;
            case FLAT_SHADER:
                zgl_flat_uniform_t flatUniformData = {
                    .modelMatrix = object.transform,
                    .modelInvRotationMatrixTransposed = zgl_transpose(zgl_inverse(object.rotation)),
                    .viewProjectionMatrix = game->camera.viewProjMatrix,
                    .lightSources = game->lightSources,
                    .bilinearFiltering = game->bilinearFiltering,
                    .materials = object.mesh->materials
                };
                zgl_render_object3D(&object, &flatUniformData, game->camera, game->canvas, zgl_flat_vertex_shader, zgl_flat_fragment_shader, game->renderOptions);
                break;
            case GOURAUD_SHADER:
                zgl_gourard_uniform_t gourardUniformData = {
                    .modelMatrix = object.transform,
                    .modelInvRotationMatrixTransposed = zgl_transpose(zgl_inverse(object.rotation)),
                    .viewProjectionMatrix = game->camera.viewProjMatrix,
                    .lightSources = game->lightSources,
                    .bilinearFiltering = game->bilinearFiltering,
                    .materials = object.mesh->materials
                };
                zgl_render_object3D(&object, &gourardUniformData, game->camera, game->canvas, zgl_gourard_vertex_shader, zgl_gourard_fragment_shader, game->renderOptions);
                break;
            case PHONG_SHADER:
                zgl_phong_uniform_t phongUniformData = {
                    .modelMatrix = object.transform,
                    .modelInvRotationMatrixTransposed = zgl_transpose(zgl_inverse(object.rotation)),
                    .viewProjectionMatrix = game->camera.viewProjMatrix,
                    .lightSources = game->lightSources,
                    .bilinearFiltering = game->bilinearFiltering,
                    .materials = object.mesh->materials
                };
                zgl_render_object3D(&object, &phongUniformData, game->camera, game->canvas, zgl_phong_vertex_shader, zgl_phong_fragment_shader, game->renderOptions);
                break;
        }
    }
}

void drawLights(game_state_t* game) {
    // Use basic shader to draw lights
    for (int i = 0; i < game->lightSources.numPointLights; i++) {
        zgl_basic_uniform_t uniformData = {
            .modelviewprojection = zgl_mul_mat(game->camera.viewProjMatrix, game->pointLightObjects[i].transform),
        };
        zgl_render_object3D(&game->pointLightObjects[i], &uniformData, game->camera, game->canvas, zgl_basic_vertex_shader, zgl_basic_fragment_shader, game->renderOptions);
    }
}

void render(game_state_t* game) {
    ZGL_DEBUG_PRINT("INFO: Rendering scene\n");
    zgl_canvas_t canvas = game->canvas;

    zgl_clear_depth_buffer(canvas);

    ZGL_DEBUG_PRINT("INFO: Drawing background\n");
    zgl_render_fill(game->backgroundColor, canvas);

    // Draw 3D Objects
    if (game->draw3DObjects) {
        ZGL_DEBUG_PRINT("INFO: Drawing 3D Objects\n");
        drawObjects(game);
    }

    // Draw lights
    if (game->drawLights) {
        ZGL_DEBUG_PRINT("INFO: Drawing lights\n");
        drawLights(game);
    }

    ZGL_DEBUG_PRINT("INFO: Update backbuffer\n");
    SDL_UpdateTexture(game->texture, NULL, game->canvas.frameBuffer, PITCH);
    SDL_RenderCopy(game->renderer, game->texture, NULL, NULL);
}

#ifdef DEBUGUI
void renderDebugUI(game_state_t* game) {
    if (game->showGUI) {
        ZGL_DEBUG_PRINT("INFO: Rendering GUI\n");
        nk_sdl_render(NK_ANTI_ALIASING_ON);
    }
}
#endif // DEBUGUI

void destroy(game_state_t* game) {
    free(game->canvas.frameBuffer);
    free(game->canvas.depthBuffer);

    // Free all triangles, vertices and materials from meshes
    for (int i = 0; i < game->numMeshes; i++) {
        free(game->meshes[i].materials);
        free(game->meshes[i].triangles);
        free(game->meshes[i].vertices);
    }

    free(game->meshes);
    free(game->objects);
    SDL_DestroyTexture(game->texture);
    SDL_DestroyWindow(game->window);

    #ifdef DEBUGUI
    nk_sdl_shutdown();
    #endif // DEBUGUI

    SDL_Quit();
}

int main(int argc, char* argv[])
{
    game_state_t* game = init();

    while (game->running) {
        handleEvents(game);

        #ifdef DEBUGUI
        updateDebugUI(game);
        #endif // DEBUGUI

        update(game);
        render(game);

        #ifdef DEBUGUI
        renderDebugUI(game);
        #endif // DEBUGUI

        // Present frame and compute FPS
        ZGL_DEBUG_PRINT("INFO: Present frame\n");
        SDL_RenderPresent(game->renderer);
        uint64_t currentTime = SDL_GetPerformanceCounter();
        game->elapsedTime = 1000.0 * (currentTime - game->lastTime) / SDL_GetPerformanceFrequency();
        ZGL_DEBUG_PRINT("INFO: Frame rendered in %.00lf ms (%.0f FPS)\n", game->elapsedTime, floor(1000.0f / game->elapsedTime));
        game->lastTime = currentTime;
    }

    ZGL_DEBUG_PRINT("INFO: Closing\n");
    destroy(game);
    return 0;
}
