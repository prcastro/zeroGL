// #define SR_DEBUG
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
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "simplerenderer.h"
#include "objloader.h"

#define WIDTH 1066
#define HEIGHT 600
#define ROTATION_SPEED 15.0f // degrees per second
#define VIEWPORT_WIDTH (WIDTH /(float) HEIGHT)
#define VIEWPORT_HEIGHT 1.0f
#define VIEWPORT_DISTANCE 1.0f
#define PIXEL_DEPTH 4
#define PITCH (PIXEL_DEPTH * WIDTH)

typedef enum {
    BASIC_SHADER,
    GOURAUD_SHADER,
    PHONG_SHADER // TODO: Uninplemented
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
    canvas_t      canvas;
    uint32_t      backgroundColor;
    uint16_t      renderOptions;
    int           drawLights;
    int           draw3DObjects;
    shader_type_t shaderType;
    
    // Game objects
    int              numMeshes;
    mesh_t*          meshes;
    int              numObjects;
    object3D_t*      objects;
    light_sources_t  lightSources;
    object3D_t*      pointLightObjects;
    camera_t        camera;
    float            rotationSpeed;
    const uint8_t*   keys;

    // GUI
    #ifdef DEBUGUI
    struct nk_context* nuklearContext;
    int                showGUI;
    int                toggleGUIKeyPressed;
    #endif // DEBUG
} game_state_t;

game_state_t* init() {
    DEBUG_PRINT("INFO: Initializing game objects\n");

    DEBUG_PRINT("INFO: Initializing SDL\n");
    if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
        fprintf(stderr, "ERROR: error initializing SDL: %s\n", SDL_GetError());
        exit(-1);
    }

    DEBUG_PRINT("INFO: Loading meshes and objects\n");
    int numObjects = 1;
    int numMeshes = 5;
    mesh_t* meshes = (mesh_t*) malloc(numMeshes * sizeof(mesh_t));
    object3D_t *objects = (object3D_t*) malloc(numObjects * sizeof(object3D_t));
    if (meshes == NULL || objects == NULL) {
        fprintf(stderr, "ERROR: 3D objects memory couldn't be allocated.\n");
        exit(-1);
    }
        
    meshes[0] = *loadObjFile("assets/light.obj", false);
    meshes[1] = *loadObjFile("assets/snake/snake.obj", true);
    meshes[2] = *loadObjFile("assets/engineer/engineer.obj", false);
    meshes[3] = *loadObjFile("assets/cube.obj", false);

    // Define a debug mesh
    vec3_t* vertices = (vec3_t*) malloc(3 * sizeof(vec3_t));
    triangle_t* triangles = (triangle_t*) malloc(1 * sizeof(triangle_t));
    material_t* materials = (material_t*) malloc(1 * sizeof(material_t));
    if (vertices == NULL || triangles == NULL || materials == NULL) {
        fprintf(stderr, "ERROR: Debug mesh memory couldn't be allocated.\n");
        exit(-1);
    }

    vertices[0] = (vec3_t) {0, 0, 0};
    vertices[1] = (vec3_t) {1, 0, 0};
    vertices[2] = (vec3_t) {0, 1, 0};
    triangles[0] = (triangle_t) {0, 2, 1, 0, 0, 0, 0, 0, 0, 0};
    materials[0] = (material_t) {"RedMaterial", COLOR_RED, COLOR_RED, 0.0f, 0, 0, NULL};
    meshes[4] = (mesh_t) {
        .name = "Debug",
        .numVertices = 3,
        .numTriangles = 1,
        .numTextureCoords = 0,
        .vertices = vertices,
        .triangles = triangles,
        .materials = materials
    };


    objects[0] = makeObject(&meshes[2], (vec3_t) {0, 0, 0}, 1.0 , IDENTITY_M4x4);
    // objects[0] = makeObject(&meshes[3], (vec3_t) {0, 0, 0}, 1.0 , IDENTITY_M4x4);
    // objects[0] = makeObject(&meshes[4], (vec3_t) {0, 0, 0}, 1.0 , IDENTITY_M4x4);

    DEBUG_PRINT("INFO: Loading lights\n");
    int numAmbientLights = 1;
    int numDirLights = 1;
    int numPointLights = 1;
    ambient_light_t* ambientLights = (ambient_light_t*) malloc(numAmbientLights * sizeof(ambient_light_t));
    dir_light_t* directionalLights = (dir_light_t*) malloc(numDirLights * sizeof(dir_light_t));
    point_light_t* pointLights = (point_light_t*) malloc(numPointLights * sizeof(point_light_t));
    object3D_t *pointLightObjects = (object3D_t*) malloc(numPointLights * sizeof(object3D_t));
    if (pointLights == NULL || directionalLights == NULL || ambientLights == NULL || pointLightObjects == NULL) {
        fprintf(stderr, "ERROR: Lights memory couldn't be allocated.\n");
        exit(-1);
    }

    ambientLights[0] = (ambient_light_t) {0.4};
    directionalLights[0] = (dir_light_t) {0.0, {0.0, -1.0, 1.0}};
    pointLights[0] = (point_light_t) {0.9, {-0.5, 1.5, -2.0}};

    pointLightObjects[0] = makeObject(&meshes[0], pointLights[0].position, 0.05, IDENTITY_M4x4);

    DEBUG_PRINT("INFO: Initializing game state\n");
    uint32_t *frameBuffer = (uint32_t*) malloc(WIDTH * HEIGHT * sizeof(uint32_t));
    float *depthBuffer = (float*) malloc(WIDTH * HEIGHT * sizeof(float));
    game_state_t* game = (game_state_t*) malloc(sizeof(game_state_t));
    if (game == NULL || frameBuffer == NULL || depthBuffer == NULL) {
        fprintf(stderr, "ERROR: Game state memory and buffers couldn't be allocated.\n");
        exit(-1);
    }

    struct canvas_t canvas = {
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
    game->texture = SDL_CreateTexture(game->renderer, SDL_PIXELFORMAT_BGRA32, SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);
    game->canvas = canvas;
    game->backgroundColor = COLOR_BLACK;
    game->drawLights    = 1;
    game->draw3DObjects = 1;
    game->shaderType = GOURAUD_SHADER;
    game->renderOptions = DIFFUSE_LIGHTING | SPECULAR_LIGHTING | BACKFACE_CULLING | FUSTRUM_CULLING | SHADED_FLAT;
    game->numMeshes = numMeshes;
    game->meshes = meshes;
    game->numObjects = numObjects;
    game->objects = objects;

    // Lights
    game->lightSources = (light_sources_t) {
        .ambientLights = ambientLights,
        .numAmbientLights = numAmbientLights,
        .directionalLights = directionalLights,
        .numDirectionalLights = numDirLights,
        .pointLights = pointLights,
        .numPointLights = numPointLights
    };
    
    game->pointLightObjects = pointLightObjects;

    game->camera = makeCamera(
        (vec3_t) {0, 0, -5}, // position
        (vec3_t) {0, 0, 1}, // direction
        (vec3_t) {0, 1, 0}, // up
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
    DEBUG_PRINT("INFO:  Initializing Dear ImGui\n");
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
    DEBUG_PRINT("INFO: Handle events\n");

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
            DEBUG_PRINT("INFO: Quitting application\n");
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
        DEBUG_PRINT("INFO: Updating GUI\n");
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
                nk_layout_row_dynamic(ctx, row_size, 1);

                nk_layout_row_dynamic(ctx, row_size, 3);
                nk_bool isBasic = game->shaderType == BASIC_SHADER;
                if (nk_radio_label(ctx, "Basic", &isBasic)) {
                    game->shaderType = BASIC_SHADER; 
                };

                nk_bool isGouraud = game->shaderType == GOURAUD_SHADER;
                if (nk_radio_label(ctx, "Gouraud", &isGouraud)) {
                    game->shaderType = GOURAUD_SHADER;
                };

                nk_bool isPhong = game->shaderType == PHONG_SHADER;
                if (nk_radio_label(ctx, "Phong", &isPhong)) {
                    game->shaderType = PHONG_SHADER;
                };

                nk_layout_row_dynamic(ctx, row_size, 2);

                nk_bool isDiffuse = game->renderOptions & DIFFUSE_LIGHTING;
                nk_checkbox_label(ctx, "Difuse", &isDiffuse);
                game->renderOptions = isDiffuse ? game->renderOptions | DIFFUSE_LIGHTING : game->renderOptions & ~DIFFUSE_LIGHTING;

                nk_bool isSpecular = game->renderOptions & SPECULAR_LIGHTING;
                nk_checkbox_label(ctx, "Specular", &isSpecular);
                game->renderOptions = isSpecular ? game->renderOptions | SPECULAR_LIGHTING : game->renderOptions & ~SPECULAR_LIGHTING;

                if (nk_tree_push(ctx, NK_TREE_NODE, "Ambient", NK_MAXIMIZED)) {
                    for (int i = 0; i < game->lightSources.numAmbientLights; i++) {
                        nk_layout_row_dynamic(ctx, row_size * 4, 1);
                        char label[100];
                        sprintf(label, "Ambient Light %d", i);
                        if (nk_group_begin(ctx, label, NK_WINDOW_TITLE|NK_WINDOW_BORDER|NK_WINDOW_NO_SCROLLBAR)) {
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity", 0.0f, &game->lightSources.ambientLights[i].intensity, 3.0f, 0.1f, 0.1f);
                            nk_group_end(ctx);
                        }
                    }
                    nk_tree_pop(ctx);
                }

                if (nk_tree_push(ctx, NK_TREE_NODE, "Directional", NK_MAXIMIZED)) {
                    for (int i = 0; i < game->lightSources.numDirectionalLights; i++) {
                        nk_layout_row_dynamic(ctx, row_size * 8, 1);
                        char label[100];
                        sprintf(label, "Directional Light %d", i);
                        if (nk_group_begin(ctx, label, NK_WINDOW_TITLE|NK_WINDOW_BORDER|NK_WINDOW_NO_SCROLLBAR)) {
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity", 0.0f, &game->lightSources.directionalLights[i].intensity, 3.0f, 0.1f, 0.1f);
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
                        nk_layout_row_dynamic(ctx, row_size * 8, 1);
                        char label[100];
                        sprintf(label, "Point Light %d", i);
                        if (nk_group_begin(ctx, label, NK_WINDOW_TITLE|NK_WINDOW_BORDER|NK_WINDOW_NO_SCROLLBAR)) {
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "Intensity", 0.0f, &game->lightSources.pointLights[i].intensity, 3.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "x", -10.0f, &game->lightSources.pointLights[i].position.x, 10.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "y", -10.0f, &game->lightSources.pointLights[i].position.y, 10.0f, 0.1f, 0.1f);
                            nk_layout_row_dynamic(ctx, row_size, 1);
                            nk_property_float(ctx, "z", -10.0f, &game->lightSources.pointLights[i].position.z, 10.0f, 0.1f, 0.1f);
                            game->pointLightObjects[i] = makeObject(&game->meshes[0], game->lightSources.pointLights[i].position, 0.05, IDENTITY_M4x4);
                            nk_group_end(ctx);
                        }
                    }
                    nk_tree_pop(ctx);
                }
                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Camera", NK_MAXIMIZED)) {
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
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_bool isBackfaceCulling = game->renderOptions & BACKFACE_CULLING;
                nk_checkbox_label(ctx, "Backface culling", &isBackfaceCulling);
                game->renderOptions = isBackfaceCulling ? game->renderOptions | BACKFACE_CULLING : game->renderOptions & ~BACKFACE_CULLING;

                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_bool isFustrumCulling = game->renderOptions & FUSTRUM_CULLING;
                nk_checkbox_label(ctx, "Fustrum culling", &isFustrumCulling);
                game->renderOptions = isFustrumCulling ? game->renderOptions | FUSTRUM_CULLING : game->renderOptions & ~FUSTRUM_CULLING;

                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_bool isBilinearFiltering = game->renderOptions & BILINEAR_FILTERING;
                nk_checkbox_label(ctx, "Bilinear filtering", &isBilinearFiltering);
                game->renderOptions = isBilinearFiltering ? game->renderOptions | BILINEAR_FILTERING : game->renderOptions & ~BILINEAR_FILTERING;

                nk_tree_pop(ctx);
            }

            if (nk_tree_push(ctx, NK_TREE_NODE, "Background", NK_MAXIMIZED)) {
                nk_layout_row_dynamic(ctx, row_size, 1);
                nk_label(ctx, "Color", NK_TEXT_LEFT);
                nk_layout_row_dynamic(ctx, row_size * 10, 1);
                uint8_t r, g, b;
                colorFromUint32(game->backgroundColor, &r, &g, &b);
                struct nk_colorf nkBackgroundColor = nk_color_cf(nk_rgb(r, g, b));
                nk_color_pick(ctx, &nkBackgroundColor, NK_RGBA);
                game->backgroundColor = colorFromFloats(nkBackgroundColor.r, nkBackgroundColor.g, nkBackgroundColor.b);
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
    camera_t *camera = &game->camera;

    // TODO: Store this in the camera struct
    vec3_t cameraRight = crossProduct(camera->up, camera->direction);
    
    float elapsedTime = game->elapsedTime / 1000.0f;
    float movementSpeed = camera->movementSpeed * elapsedTime;
    float turningSpeed = camera->turningSpeed * elapsedTime;
    vec3_t newCameraPosition = camera->position;
    vec3_t newCameraDirection = camera->direction;
    vec3_t newCameraUp = camera->up;
    

    if (keys[SDL_SCANCODE_A]) {
        // Translate left       
        newCameraPosition = add(newCameraPosition, mulScalarV3(-movementSpeed, cameraRight));
    }

    if (keys[SDL_SCANCODE_D]) {
        // Translate right
        newCameraPosition = add(newCameraPosition, mulScalarV3(movementSpeed, cameraRight));
    }

    if (keys[SDL_SCANCODE_S]) {
        // Translate back
        newCameraPosition = add(newCameraPosition, mulScalarV3(-movementSpeed, newCameraDirection));
    }

    if (keys[SDL_SCANCODE_W]) {
        // Translate forward
        newCameraPosition = add(newCameraPosition, mulScalarV3(movementSpeed, newCameraDirection));
    }

    if (keys[SDL_SCANCODE_PAGEDOWN]) {
        // Translate down
        newCameraPosition = add(newCameraPosition, mulScalarV3(-movementSpeed, newCameraUp));
    }

    if (keys[SDL_SCANCODE_PAGEUP]) {
        // Translate up
        newCameraPosition = add(newCameraPosition, mulScalarV3(movementSpeed, newCameraUp));
    }

    if (keys[SDL_SCANCODE_RIGHT]) {
        // Rotate right around up axis using quaternion
        quaternion_t rotation = quaternionFromAngleAxis(turningSpeed, newCameraUp);
        newCameraDirection = rotateVectorByQuaternion(newCameraDirection, rotation);
        cameraRight = crossProduct(newCameraUp, newCameraDirection);
    }

    if (keys[SDL_SCANCODE_LEFT]) {
        // Rotate left around up axis
        quaternion_t rotation = quaternionFromAngleAxis(-turningSpeed, newCameraUp);
        newCameraDirection = rotateVectorByQuaternion(newCameraDirection, rotation);
        cameraRight = crossProduct(newCameraUp, newCameraDirection);
    }

    if (keys[SDL_SCANCODE_UP]) {
        // Rotate up around right axis
        quaternion_t rotation = quaternionFromAngleAxis(-turningSpeed, cameraRight);
        newCameraDirection = rotateVectorByQuaternion(newCameraDirection, rotation);
        newCameraUp = rotateVectorByQuaternion(newCameraUp, rotation);

    }

    if (keys[SDL_SCANCODE_DOWN]) {
        // Rotate down around right axis
        quaternion_t rotation = quaternionFromAngleAxis(turningSpeed, cameraRight);
        newCameraDirection = rotateVectorByQuaternion(newCameraDirection, rotation);
        newCameraUp = rotateVectorByQuaternion(newCameraUp, rotation);
    }

    if (keys[SDL_SCANCODE_Q]) {
        // Rotate left around direction axis
        quaternion_t rotation = quaternionFromAngleAxis(turningSpeed, newCameraDirection);
        newCameraUp = rotateVectorByQuaternion(newCameraUp, rotation);
        cameraRight = crossProduct(newCameraUp, newCameraDirection);
    }

    if (keys[SDL_SCANCODE_E]) {
        // Rotate right around direction axis
        quaternion_t rotation = quaternionFromAngleAxis(-turningSpeed, newCameraDirection);
        newCameraUp = rotateVectorByQuaternion(newCameraUp, rotation);
        cameraRight = crossProduct(newCameraUp, newCameraDirection);
    }

    game->camera = makeCamera(
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

object3D_t rotateObjectY(object3D_t object, float degrees) {
    return makeObject(object.mesh, object.translation, object.scale, mulMM4(rotationY(degrees), object.rotation));
}

// TODO: Use quaternions for rotation
void animateObjects(game_state_t* game) {
    float degrees = game->rotationSpeed * (game->elapsedTime / 1000.0f);
    game->objects[0] = rotateObjectY(game->objects[0], fmod(degrees, 360.0f));
}

void update(game_state_t* game) {
    DEBUG_PRINT("INFO: Update game state\n");
    updateCameraPosition(game);
    animateObjects(game);
}

void drawObjects(game_state_t* game) {
    for (int i = 0; i < game->numObjects; i++) {
        if (game->shaderType == GOURAUD_SHADER) {
            gourardUniformData_t uniformData = {
                .modelMatrix = game->objects[i].transform,
                .viewProjectionMatrix = game->camera.viewProjMatrix,
                .lightSources = game->lightSources,
            };
            drawObject(&game->objects[i], &uniformData, game->camera, game->canvas, gourardVertexShader, gourardFragmentShader, game->renderOptions);
        } else if (game->shaderType == BASIC_SHADER) {
            basicUniformData_t uniformData = {
                .modelviewprojection = mulMM4(game->camera.viewProjMatrix, game->objects[i].transform),
            };
            drawObject(&game->objects[i], &uniformData, game->camera, game->canvas, basicVertexShader, basicFragmentShader, game->renderOptions);
        } else if (game->shaderType == PHONG_SHADER) {
            printf("ERROR: Phong shader not implemented\n");
        }
    }
}

void drawLights(game_state_t* game) {
    // Use basic shader to draw lights
    for (int i = 0; i < game->lightSources.numPointLights; i++) {
        basicUniformData_t uniformData = {
            .modelviewprojection = mulMM4(game->camera.viewProjMatrix, game->pointLightObjects[i].transform),
        };
        drawObject(&game->pointLightObjects[i], &uniformData, game->camera, game->canvas, basicVertexShader, basicFragmentShader, game->renderOptions);
    }
}

void render(game_state_t* game) {
    DEBUG_PRINT("INFO: Rendering scene\n");
    canvas_t canvas = game->canvas;

    // Init depthBuffer
    for (int i = 0; i < canvas.width * canvas.height; i++) {
        canvas.depthBuffer[i] = 0.0;
    }

    DEBUG_PRINT("INFO: Drawing background\n");
    for (int i = 0; i < canvas.width * canvas.height; i++) {
        game->canvas.frameBuffer[i] = game->backgroundColor;
    }

    // Draw 3D Objects
    if (game->draw3DObjects) {
        DEBUG_PRINT("INFO: Drawing 3D Objects\n");
        drawObjects(game);
    }

    // Draw lights
    if (game->drawLights) {
        DEBUG_PRINT("INFO: Drawing lights\n");
        drawLights(game);
    }
    
    DEBUG_PRINT("INFO: Update backbuffer\n");
    SDL_UpdateTexture(game->texture, NULL, game->canvas.frameBuffer, PITCH);
    SDL_RenderCopy(game->renderer, game->texture, NULL, NULL);
}

#ifdef DEBUGUI
void renderDebugUI(game_state_t* game) {
    if (game->showGUI) {
        DEBUG_PRINT("INFO: Rendering GUI\n");
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
        DEBUG_PRINT("INFO: Present frame\n");
        SDL_RenderPresent(game->renderer);
        uint64_t currentTime = SDL_GetPerformanceCounter();
        game->elapsedTime = 1000.0 * (currentTime - game->lastTime) / SDL_GetPerformanceFrequency();
        DEBUG_PRINT("INFO: Frame rendered in %.00lf ms (%.0f FPS)\n", game->elapsedTime, floor(1000.0f / game->elapsedTime));
        game->lastTime = currentTime;
    }

    DEBUG_PRINT("INFO: Closing\n");
    destroy(game);
    return 0;
}
