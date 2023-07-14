#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

typedef struct color_t {
  uint8_t r;
  uint8_t g;
  uint8_t b;
} color_t;

const color_t COLOR_WHITE      = {255, 255, 255};
const color_t COLOR_BLACK      = {0,   0,   0  };
const color_t COLOR_GREEN      = {0,   255, 0  };
const color_t COLOR_BLUE       = {0,   0,   255};
const color_t COLOR_RED        = {255, 0,   0  };
const color_t COLOR_YELLOW     = {255, 255, 0  };
const color_t COLOR_PURPLE     = {255, 0,   255};
const color_t COLOR_CYAN       = {0,   255, 255};

uint32_t colorToUint32(color_t c);
color_t colorFromUint32(uint32_t c);
float clamp(float v, float max);
color_t mulScalarColor(double x, color_t color);
color_t sumColors3(color_t c0, color_t c1, color_t c2);
color_t sumColors(color_t c0, color_t c1);
color_t colorFromFloats(float r, float g, float b);

#endif