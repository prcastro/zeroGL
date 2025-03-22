#define SDL_MAIN_HANDLED

#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#define WIDTH 800
#define HEIGHT 600
#define PIXEL_DEPTH 4
#define PITCH (PIXEL_DEPTH * WIDTH)

int main(int argc, char* argv[]) {
    // Create window and renderer
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_CreateWindowAndRenderer("Hello World", WIDTH, HEIGHT, SDL_WINDOW_RESIZABLE, &window, &renderer);

    SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);
    int32_t* pixels = (int32_t*) malloc(800 * 600 * sizeof(int32_t));

    // Main loop
    int running = 1;
    while (running) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) {
                running = 0;
            }
        }

        // Clear texture
        for (int i = 0; i < WIDTH * HEIGHT; i++) {
            pixels[i] = 0xFF000000;
        }

        // Draw a red rectangle in texture
        for (int y = 100; y < 300; y++) {
            for (int x = 100; x < 300; x++) {
                pixels[y * 800 + x] = 0xFF0000FF; // RGBA
            }
        }

        SDL_UpdateTexture(texture, NULL, pixels, PITCH);
        
        // Render texture
        SDL_RenderTexture(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    return 0;
}
