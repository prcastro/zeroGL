#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

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

#endif // COLOR_H