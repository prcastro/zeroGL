#include "draw.h"
#include "color.h"
#include "vectors.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void drawPixel(int i, int j, uint32_t color, uint32_t* frameBuffer) {
    if ((i >= 0) && (i < WIDTH) && (j >= 0) && (j < HEIGHT)) {
        frameBuffer[j * WIDTH + i] = color;
    }
}

void drawPixelDepthBuffer(int i, int j, float z, uint32_t color, float* depthBuffer, uint32_t* frameBuffer) {
    if ((i >= 0) && (i < WIDTH) && (j >= 0) && (j < HEIGHT)) {
        int position = j * WIDTH + i;
        frameBuffer[position] = color;
        depthBuffer[position] = z;
    }
}

void drawLine(int x0, int x1, int y0, int y1, color_t color, uint32_t* frameBuffer) {
    int delta_x = (x1 - x0);
    int delta_y = (y1 - y0);

    int longest_side_length = (abs(delta_x) >= abs(delta_y)) ? abs(delta_x) : abs(delta_y);

    float x_inc = delta_x / (float)longest_side_length; 
    float y_inc = delta_y / (float)longest_side_length;

    float current_x = x0;
    float current_y = y0;
    for (int i = 0; i <= longest_side_length; i++) {
        drawPixel(round(current_x), round(current_y), colorToUint32(color), frameBuffer);
        current_x += x_inc;
        current_y += y_inc;
    }
}

point_t projectVertex(vec3_t v) {
  return (point_t) {
    (int) (v.x * VIEWPORT_DISTANCE / v.z  * WIDTH/VIEWPORT_WIDTH + WIDTH/2),
    (int) (HEIGHT/2 - (v.y * VIEWPORT_DISTANCE / v.z * HEIGHT/VIEWPORT_HEIGHT) - 1),
    1.0f / v.z
  };
}

vec3_t unprojectPoint(point_t p) {
  return (vec3_t) {
    (p.x - WIDTH/2) * (VIEWPORT_WIDTH / WIDTH) / (p.invz * VIEWPORT_DISTANCE),
    (HEIGHT/2 - p.y - 1) / (VIEWPORT_DISTANCE * p.invz * HEIGHT/VIEWPORT_HEIGHT),
    1.0f / p.invz
  };
}