#ifndef DRAW_H
#define DRAW_H

#include "color.h"
#include "vectors.h"

#define WIDTH 1066
#define HEIGHT 600
#define VIEWPORT_WIDTH (WIDTH / (float) HEIGHT)
#define VIEWPORT_HEIGHT 1
#define VIEWPORT_DISTANCE 1.0f
#define PIXEL_DEPTH 4
#define PITCH (PIXEL_DEPTH * WIDTH)

typedef struct point_t {
  int   x, y;
  float invz;
} point_t;


void drawPixel(int i, int j, uint32_t color, uint32_t* frameBuffer);
void drawPixelDepthBuffer(int i, int j, float z, uint32_t color, float* depthBuffer, uint32_t* frameBuffer);
void drawLine(int x0, int y0, int x1, int y1, color_t color, uint32_t* frameBuffer);
point_t projectVertex(vec3_t v);
vec3_t unprojectPoint(point_t p);

#endif